// Same configuration as in wcls-multitpc-sim-drift-simchannel.jsonnet
// except the sed_label variable name in wclsDepoFluxWriter set to 'ionization'

local g = import 'pgraph.jsonnet';
local f = import 'pgrapher/common/funcs.jsonnet';
local wc = import 'wirecell.jsonnet';

local io = import 'pgrapher/common/fileio.jsonnet';
local tools_maker = import 'pgrapher/experiment/icarus/icarus_tools.jsonnet';
local base = import 'pgrapher/experiment/icarus/simparams.jsonnet';

// load the electronics response parameters
local er_params = [
  {
    gain: std.extVar('gain0')*wc.mV/wc.fC,
    shaping: std.extVar('shaping0')*wc.us,
  },

  {
    gain: std.extVar('gain1')*wc.mV/wc.fC,
    shaping: std.extVar('shaping1')*wc.us,
  },

  {
    gain: std.extVar('gain2')*wc.mV/wc.fC,
    shaping: std.extVar('shaping2')*wc.us,
  },
];

local params = base {
  lar: super.lar {
    // Longitudinal diffusion constant
    DL: std.extVar('DL') * wc.cm2 / wc.ns,
    // Transverse diffusion constant
    DT: std.extVar('DT') * wc.cm2 / wc.ns,
    // Electron lifetime
    lifetime: std.extVar('lifetime') * wc.us,
    // Electron drift speed, assumes a certain applied E-field
    // drift_speed: std.extVar('driftSpeed') * wc.mm / wc.us,
  },
  files: super.files {
    fields: [ std.extVar('files_fields'), ],
  },

  rc_resp: if std.extVar('file_rcresp') != "" then
  {
    // "icarus_fnal_rc_tail.json"
    filename: std.extVar('file_rcresp'),
    postgain: 1.0,
    start: 0.0,
    tick: 0.4*wc.us,
    nticks: 4255,
    type: "JsonElecResponse",
    rc_layers: 1
  }
  else super.rc_resp,

  elec: std.mapWithIndex(function (n, eparam)
    super.elec[n] + {
      gain: eparam.gain,
      shaping: eparam.shaping,
    }, er_params),
};

local tools = tools_maker(params);

local sim_maker = import 'pgrapher/experiment/icarus/sim.jsonnet';
local sim = sim_maker(params, tools);

local nanodes = std.length(tools.anodes);
local anode_iota = std.range(0, nanodes - 1);


local output = 'wct-sim-ideal-sig.npz';


//local depos = g.join_sources(g.pnode({type:"DepoMerger", name:"BlipTrackJoiner"}, nin=2, nout=1),
//                             [sim.ar39(), sim.tracks(tracklist)]);
// local depos = sim.tracks(tracklist, step=1.0 * wc.mm);

local wcls_maker = import 'pgrapher/ui/wcls/nodes.jsonnet';
local wcls = wcls_maker(params, tools);
//local wcls_input = {
//  // depos: wcls.input.depos(name="", art_tag="ionization"),
//  depos: wcls.input.depos(name='electron', art_tag='ionization'),  // default art_tag="blopper"
//};

//Haiwang DepoSetSource Implementation:
local wcls_input = {
	depos: wcls.input.depos(name="", art_tag="IonAndScint"),
	deposet: g.pnode({
        	type: 'wclsSimDepoSetSource',
        	name: "electron",
        	data: {
            	model: "",
            	scale: -1, //scale is -1 to correct a sign error in the SimDepoSource converter.
		art_tag: "ionization", //name of upstream art producer of depos "label:instance:processName"
            	assn_art_tag: "",
              id_is_track: false,    // Use this for "id-is-index" in the output
        	},
    	}, nin=0, nout=1),
};

// Collect all the wc/ls output converters for use below.  Note the
// "name" MUST match what is used in theh "outputers" parameter in the
// FHiCL that loads this file.
local mega_anode = {
  type: 'MegaAnodePlane',
  name: 'meganodes',
  data: {
    anodes_tn: [wc.tn(anode) for anode in tools.anodes],
  },
};

// A ``duo'' anode consists of two ``splits''
local duoanodes = [
{
  type: 'MegaAnodePlane',
  name: 'duoanode%d' %n,
  data: {
    // anodes_tn: ["AnodePlane:anode110", "AnodePlane:anode120"],
    anodes_tn: [wc.tn(a) for a in tools.anodes[2*n:2*(n+1)]],
    // anodes_tn: [wc.tn(tools.anodes[2*n]), wc.tn(tools.anodes[2*n+1])],
  }, 
}
for n in std.range(0,3)];
local volname = ["EE", "EW", "WE", "WW"];
local wcls_output = {
  // ADC output from simulation
  // sim_digits: wcls.output.digits(name="simdigits", tags=["orig"]),
  sim_digits: [ 
  g.pnode({
    type: 'wclsFrameSaver',
    name: 'simdigits%d' %n,
    data: {
      // anode: wc.tn(tools.anode),
      anode: wc.tn(duoanodes[n]),
      digitize: true,  // true means save as RawDigit, else recob::Wire
      //frame_tags: ['daq%d' %n],
      frame_tags: ['TPC%s' %volname[n]],
      // Three options for nticks:
      // - If nonzero, force number of ticks in output waveforms.
      // - If zero, use whatever input data has. (default)
      // - If -1, use value as per LS's detector properties service.
      // nticks: params.daq.nticks,
      // nticks: -1,
      // chanmaskmaps: ['bad'],
    },
  }, nin=1, nout=1, uses=[duoanodes[n]])
  for n in std.range(0,3)],

  // The noise filtered "ADC" values.  These are truncated for
  // art::Event but left as floats for the WCT SP.  Note, the tag
  // "raw" is somewhat historical as the output is not equivalent to
  // "raw data".
  nf_digits: wcls.output.digits(name='nfdigits', tags=['raw']),

  // The output of signal processing.  Note, there are two signal
  // sets each created with its own filter.  The "gauss" one is best
  // for charge reconstruction, the "wiener" is best for S/N
  // separation.  Both are used in downstream WC code.
  sp_signals: wcls.output.signals(name='spsignals', tags=['gauss', 'wiener']),

  // save "threshold" from normal decon for each channel noise
  // used in imaging
  sp_thresholds: wcls.output.thresholds(name='spthresholds', tags=['threshold']),
};

//local deposio = io.numpy.depos(output);
local drifter = sim.drifter;
local setdrifter = g.pnode({
            type: 'DepoSetDrifter',
            data: {
                drifter: "Drifter"
            }
        }, nin=1, nout=1,
        uses=[drifter]);
local bagger = sim.make_bagger();

// signal plus noise pipelines
//local sn_pipes = sim.signal_pipelines;
// local sn_pipes = sim.splusn_pipelines;
local analog_pipes = sim.analog_pipelines;

local perfect = import 'pgrapher/experiment/icarus/chndb-base.jsonnet';
local chndb = [{
  type: 'OmniChannelNoiseDB',
  name: 'ocndbperfect%d' % n,
  data: perfect(params, tools.anodes[n], tools.field, n){dft:wc.tn(tools.dft)},
  uses: [tools.anodes[n], tools.field, tools.dft],
} for n in anode_iota];


// local nf_maker = import 'pgrapher/experiment/icarus/nf.jsonnet';
// local nf_pipes = [nf_maker(params, tools.anodes[n], chndb_pipes[n]) for n in std.range(0, std.length(tools.anodes)-1)];
// local nf_pipes = [nf_maker(params, tools.anodes[n], chndb[n], n, name='nf%d' % n) for n in anode_iota];

local sp_maker = import 'pgrapher/experiment/icarus/sp.jsonnet';
local sp = sp_maker(params, tools);
local sp_pipes = [sp.make_sigproc(a) for a in tools.anodes];

local rng = tools.random;
local wcls_simchannel_sink_old = 
  g.pnode({
    type: 'wclsDepoSetSimChannelSink',
    name: 'postdriftold',
    data: {
      artlabel: 'simpleSCOld',  // where to save in art::Event
      anodes_tn: [wc.tn(anode) for anode in tools.anodes],
      rng: wc.tn(rng),
      tick: params.daq.tick,
      start_time: -0.34 * wc.ms, // TriggerOffsetTPC from detectorclocks_icarus.fcl
      readout_time: params.daq.readout_time,
      nsigma: 3.0,
      drift_speed: params.lar.drift_speed,
      u_to_rp: 100 * wc.mm,
      v_to_rp: 100 * wc.mm,
      y_to_rp: 100 * wc.mm,

      // GP: The shaping time of the electronics response (1.3us) shifts the peak
      //     of the field response time. Eyeballing simulation times, it does this
      //     by a bit less than the 1.3us (1us).
      //
      //     N.B. for future: there is likely an additional offset on the two induction
      //     planes due to where the deconvolution precisely defines where the "peak"
      //     of the pulse is. One may want to refine these parameters to account for that.
      //     This perturbation shouldn't be more than a tick or two.
      u_time_offset: std.extVar('time_offset_u') * wc.us,
      v_time_offset: std.extVar('time_offset_v') * wc.us,
      y_time_offset: std.extVar('time_offset_y') * wc.us,

      g4_ref_time: -1500 * wc.us, // G4RefTime from detectorclocks_icarus.fcl
      use_energy: true,
    },
  },nin=1, nout=1, uses=tools.anodes);

local wcls_simchannel_sink =
  g.pnode({
    type: 'wclsDepoFluxWriter',
    name: 'postdrift',
    data: {
      anodes: [wc.tn(anode) for anode in tools.anodes],
      field_response: wc.tn(tools.field),

      // time binning
      tick: params.daq.tick,
      window_start: -340 * wc.us, // TriggerOffsetTPC from detectorclocks_icarus.fcl
      window_duration: params.daq.readout_time,

      nsigma: 3.0,

      reference_time: -1500 * wc.us - self.window_start, // G4RefTime from detectorclocks_icarus.fcl less window start as per Brett Viren

      smear_long: 0.0,
      smear_tran: 0.0,

      time_offsets: [std.extVar('time_offset_u') * wc.us, std.extVar('time_offset_v') * wc.us, std.extVar('time_offset_y') * wc.us],

      // input from art::Event
      sed_label: 'ionization',

      // output to art::Event
      simchan_label: 'simpleSC',
    },
  },   nin=1, nout=1, uses=tools.anodes+[tools.field]);

local nicks = ["incoTPCEE","incoTPCEW","incoTPCWE","incoTPCWW", "coheTPCEE","coheTPCEW","coheTPCWE","coheTPCWW"];
local scale_int = std.extVar('int_noise_scale');
local scale_coh = std.extVar('coh_noise_scale');
local models = [
    {
        type: "GroupNoiseModel",
        name: nicks[n],
        data: {
            // This can also be given as a JSON/Jsonnet file
            spectra: params.files.noisegroups[n],
            groups: params.files.wiregroups,
            scale: if n<4 then scale_int else scale_coh,
            nsamples: params.daq.nticks,
            tick: params.daq.tick,
        }
    } for n in std.range(0,7)];

local add_noise = function(model, n,t) g.pnode({
    type: t,
    name: "addnoise%d-" %n + model.name,
    data: {
        rng: wc.tn(tools.random),
        dft: wc.tn(tools.dft),
        model: wc.tn(model),
        nsamples: params.daq.nticks,
    }}, nin=1, nout=1, uses=[tools.random, tools.dft, model]);
local noises = [add_noise(models[n], n,"IncoherentAddNoise") for n in std.range(0,3)];
local coh_noises = [add_noise(models[n],n,"CoherentAddNoise") for n in std.range(4,7)];

// local digitizer = sim.digitizer(mega_anode, name="digitizer", tag="orig");
local digitizers = [
    sim.digitizer(mega_anode, name="digitizer%d-" %n + mega_anode.name, tag="TPC%s"%volname[n])
    for n in std.range(0,3)];

local retaggers = [
g.pnode({
  type: 'Retagger',
  name: 'retagger%d' %n,
  data: {
    // Note: retagger keeps tag_rules an array to be like frame fanin/fanout.
    tag_rules: [{
      // Retagger also handles "frame" and "trace" like fanin/fanout
      // merge separately all traces like origN to orig.
      frame: {
        '.*': 'orig',
      },
      merge: {
        'orig\\d': 'daq%d' %n,
      },
    }],
  },
}, nin=1, nout=1)
for n in std.range(0, 3)];

local frame_summers = [
    g.pnode({
        type: 'FrameSummer',
        name: 'framesummer%d' %n,
        data: {
            align: true,
            offset: 0.0*wc.s,
        },
    }, nin=2, nout=1) for n in std.range(0, 3)];

local actpipes = [g.pipeline([noises[n], coh_noises[n], digitizers[n], /*retaggers[n],*/ wcls_output.sim_digits[n]], name="noise-digitizer%d" %n) for n in std.range(0,3)];
local util = import 'pgrapher/experiment/icarus/funcs.jsonnet';
local outtags = ['orig%d' % n for n in std.range(0, 3)];
local pipe_reducer = util.fansummer('DepoSetFanout', analog_pipes, frame_summers, actpipes, 'FrameFanin', 'fansummer', outtags);

// local retagger = g.pnode({
//   type: 'Retagger',
//   data: {
//     // Note: retagger keeps tag_rules an array to be like frame fanin/fanout.
//     tag_rules: [{
//       // Retagger also handles "frame" and "trace" like fanin/fanout
//       // merge separately all traces like origN to orig.
//       frame: {
//         '.*': 'orig',
//       },
//       merge: {
//         'orig\\d': 'daq',
//       },
//     }],
//   },
// }, nin=1, nout=1);

//local frameio = io.numpy.frames(output);
local sink = sim.frame_sink;

// local graph = g.pipeline([wcls_input.depos, drifter,  wcls_simchannel_sink.simchannels, bagger, pipe_reducer, retagger, wcls_output.sim_digits, sink]);
//local graph = g.pipeline([wcls_input.depos, drifter,  wcls_simchannel_sink, bagger, pipe_reducer, sink]);
local graph = g.pipeline([wcls_input.deposet, setdrifter, wcls_simchannel_sink_old, wcls_simchannel_sink, pipe_reducer, sink]);

local app = {
  type: 'Pgrapher',
  data: {
    edges: g.edges(graph),
  },
};


// Finally, the configuration sequence which is emitted.

g.uses(graph) + [app]
