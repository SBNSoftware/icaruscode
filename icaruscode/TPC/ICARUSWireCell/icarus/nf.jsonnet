// This provides some noise filtering related pnodes,

local g = import 'pgraph.jsonnet';
local wc = import 'wirecell.jsonnet';

function(params, anode, chndbobj, n, name='')
  {

    local single = {
      type: 'pdOneChannelNoise',
      name: name,
      data: {
        noisedb: wc.tn(chndbobj),
        anode: wc.tn(anode),
      },
    },
    local grouped = {
      type: 'mbCoherentNoiseSub',
      name: name,
      data: {
        noisedb: wc.tn(chndbobj),
        anode: wc.tn(anode),
        rms_threshold: 0.0,
      },
    },

    local obnf = g.pnode({
      type: 'OmnibusNoiseFilter',
      name: name,
      data: {

        // Nonzero forces the number of ticks in the waveform
        nticks: 0,

        // channel bin ranges are ignored
        // only when the channelmask is merged to `bad`
        // maskmap: {sticky: "bad", ledge: "bad", noisy: "bad"},
        channel_filters: [
          wc.tn(single),
        ],
        grouped_filters: [
          // wc.tn(grouped),
        ],
        channel_status_filters: [
        ],
        noisedb: wc.tn(chndbobj),
        intraces: 'orig%d' % anode.data.ident,  // frame tag get all traces
        outtraces: 'raw%d' % anode.data.ident,
      },
    }, uses=[chndbobj, anode, single, grouped], nin=1, nout=1),


    pipe: g.pipeline([obnf], name=name),
  }.pipe
