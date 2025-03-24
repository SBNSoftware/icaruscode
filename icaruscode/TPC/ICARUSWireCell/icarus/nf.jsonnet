// This provides some noise filtering related pnodes,

local g = import 'pgraph.jsonnet';
local wc = import 'wirecell.jsonnet';

function(anode, chndbobj, tools, name='')
  {
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
        ],
        grouped_filters: [
        ],
        channel_status_filters: [
        ],
        noisedb: wc.tn(chndbobj),
        intraces: 'orig%d' % anode.data.ident,
        outtraces: 'orig%d' % anode.data.ident,
      },
    }, uses=[chndbobj, anode], nin=1, nout=1),


    pipe: g.pipeline([obnf], name=name),
  }.pipe
