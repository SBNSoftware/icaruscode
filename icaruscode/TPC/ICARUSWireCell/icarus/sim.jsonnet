local wc = import "wirecell.jsonnet";
local g = import "pgraph.jsonnet";
local f = import "pgrapher/common/funcs.jsonnet";
local sim_maker = import "pgrapher/common/sim/nodes.jsonnet";


// return some nodes, includes base sim nodes.
function(params, tools) {
    local sim = sim_maker(params, tools),

    local nanodes = std.length(tools.anodes),

    // I rue the day that we must have an (anode) X (field) cross product!
    local ductors = sim.make_detector_ductors("nominal", tools.anodes, tools.pirs[0]),


    local zippers = [sim.make_depozipper("depozipper-"+tools.anodes[n].name, tools.anodes[n], tools.pirs[0])
                     for n in std.range(0, nanodes-1)],
    local transforms = [sim.make_depotransform("depotransform-"+tools.anodes[n].name, tools.anodes[n], tools.pirs[0])
                        for n in std.range(0, nanodes-1)],


    local transformsyz = [sim.make_depotransform_withplane("depotransform-%d-"%n+tools.anodes[std.floor(n/45)].name+"-plane%d"%std.mod(std.floor(n/15),3),tools.anodes[std.floor(n/45)], [std.mod(std.floor(n/15),3)],tools.pirs[std.mod(n,15)])
	                  for n in std.range(0, 359)],
//    local transformsyz = [sim.make_depotransform_withplane("depotransform-%d-"%n+tools.anodes[std.floor(n/45)].name+"-plane%d"%std.mod(std.floor(n/15),3),tools.anodes[std.floor(n/45)], [std.mod(std.floor(n/15),3)],tools.pirs[0])



    local depos2traces = transforms,
    local depos2tracesyz = transformsyz,
    //local depos2traces = zippers,

    local digitizers = [
        sim.digitizer(tools.anodes[n], name="digitizer-" + tools.anodes[n].name, tag="orig%d"%n)
        for n in std.range(0,nanodes-1)],

    local reframers = [
        g.pnode({
            type: 'Reframer',
            name: 'reframer-'+tools.anodes[n].name,
            data: {
                anode: wc.tn(tools.anodes[n]),
                tags: [],           // ?? what do?
                fill: 0.0,
                tbin: params.sim.reframer.tbin,
                toffset: 0,
                nticks: params.sim.reframer.nticks,
            },
        }, nin=1, nout=1) for n in std.range(0, nanodes-1)],

    local reframersyz = [
        g.pnode({
            type: 'Reframer',
            name: 'reframer-%d-'%n+tools.anodes[std.floor(n/45)].name,
            data: {
                anode: wc.tn(tools.anodes[std.floor(n/45)]),
                tags: [],           // ?? what do?
                fill: 0.0,
                tbin: params.sim.reframer.tbin,
                toffset: 0,
                nticks: params.sim.reframer.nticks,
            },
	    }, nin=1, nout=1) for n in std.range(0, 359)],
    

    // fixme: see https://github.com/WireCell/wire-cell-gen/issues/29
    local make_noise_model = function(anode, csdb=null) {
        type: "EmpiricalNoiseModel",
        name: "empericalnoise-" + anode.name,
        data: {
            anode: wc.tn(anode),
            dft: wc.tn(tools.dft),
            chanstat: if std.type(csdb) == "null" then "" else wc.tn(csdb),
            spectra_file: params.files.noise,
            nsamples: params.daq.nticks,
            period: params.daq.tick,
            wire_length_scale: 1.0*wc.cm, // optimization binning
        },
        uses: [anode, tools.dft] + if std.type(csdb) == "null" then [] else [csdb],
    },
    local noise_models = [make_noise_model(anode) for anode in tools.anodes],


    local add_noise = function(model) g.pnode({
        type: "AddNoise",
        name: "addnoise-" + model.name,
        data: {
            rng: wc.tn(tools.random),
            dft: wc.tn(tools.dft),
            model: wc.tn(model),
	    nsamples: params.daq.nticks,
            replacement_percentage: 0.02, // random optimization
        }}, nin=1, nout=1, uses=[tools.random, tools.dft, model]),

    local noises = [add_noise(model) for model in noise_models],
    
    local outtags = ["orig%d"%n for n in std.range(0, nanodes-1)],

    ret : {

        analog_pipelines: [g.pipeline([depos2traces[n], reframers[n]],
                                      name="simanalogpipe-" + tools.anodes[n].name) for n in std.range(0, nanodes-1)],
        signal_pipelines: [g.pipeline([depos2traces[n], reframers[n],  digitizers[n]],
                                      name="simsigpipe-" + tools.anodes[n].name) for n in std.range(0, nanodes-1)],
        splusn_pipelines:  [g.pipeline([depos2traces[n], reframers[n], noises[n], digitizers[n]],
                                       name="simsignoipipe-" + tools.anodes[n].name) for n in std.range(0, nanodes-1)],

        analog: f.fanpipe('DepoSetFanout', self.analog_pipelines, 'FrameFanin', "simanaloggraph", outtags),
        signal: f.fanpipe('DepoSetFanout', self.signal_pipelines, 'FrameFanin', "simsignalgraph", outtags),
        splusn: f.fanpipe('DepoSetFanout', self.splusn_pipelines, 'FrameFanin', "simsplusngraph", outtags),

        analog_pipelinesyz: [g.pipeline([depos2tracesyz[n]],
                                        name="simanalogpipe-%d-"%n + tools.anodes[std.floor(n/45)].name) for n in std.range(0, 359)],
        signal_pipelinesyz: [g.pipeline([depos2tracesyz[n], reframersyz[n],  digitizers[n]],
                                      name="simsigpipe-" + tools.anodes[n].name) for n in std.range(0, nanodes-1)],
        splusn_pipelinesyz:  [g.pipeline([depos2tracesyz[n], reframersyz[n], noises[n], digitizers[n]],
                                       name="simsignoipipe-" + tools.anodes[n].name) for n in std.range(0, nanodes-1)],

        analogyz: f.fanpipe('DepoSetFanout', self.analog_pipelinesyz, 'FrameFanin', "simanaloggraph", outtags),
        signalyz: f.fanpipe('DepoSetFanout', self.signal_pipelinesyz, 'FrameFanin', "simsignalgraph", outtags),
        splusnyz: f.fanpipe('DepoSetFanout', self.splusn_pipelinesyz, 'FrameFanin', "simsplusngraph", outtags),

    } + sim,                    // tack on base for user sugar.
}.ret
