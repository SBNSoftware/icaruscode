local wc = import "wirecell.jsonnet";
local g = import 'pgraph.jsonnet';

// The approximated sim+sigproc
function(params, tools, anode, name=null) {
    local sufix = name,
    local sp = g.pnode({
        type: 'DepoFluxSplat',
        name: sufix,
        data: {
            anode: wc.tn(anode),
            field_response: wc.tn(tools.field), // for speed and origin
            sparse: true,
            tick: params.daq.tick,
            window_start: -340 * wc.us, // TriggerOffsetTPC from detectorclocks_icarus.fcl
            window_duration: params.daq.readout_time,
            reference_time: -1500 * wc.us - self.window_start, // G4RefTime from detectorclocks_icarus.fcl less window start as per Brett Viren
            // Run wirecell-gen morse-* to find these numbers that match the extra
            // spread the sigproc induces.
            "smear_long": [
                4.55,
                4.55,
                4.55,
            ],
            "smear_tran": [
                1.55,
                1.55,
                0.175,
            ]
        },
    }, nin=1, nout=1, uses=[anode, tools.field]),
    local hio = g.pnode({
      type: 'HDF5FrameTap',
      name: 'hio_true%s' % name,
      data: {
        anode: wc.tn(anode),
        trace_tags: ['deposplat%s' % name],
        filename: "wc-true-%s.h5" % name,
        chunk: [0, 0], // ncol, nrow
        gzip: 2,
        high_throughput: true,
      },
    }, nin=1, nout=1),
    local rt = g.pnode({
        type: 'Retagger',
        name: sufix,
        data: {
            // Note: retagger keeps tag_rules an array to be like frame fanin/fanout.
            tag_rules: [{
                // Retagger also handles "frame" and "trace" like fanin/fanout
                // merge separately all traces like gaussN to gauss.
                frame: {
                ".*": "deposplat%s" % sufix
                },
                merge: {
                ".*": "deposplat%s" % sufix
                },
            }],
        },
    }, nin=1, nout=1),
    ret: g.pipeline([sp, rt, hio], name=sp.name),
}.ret
