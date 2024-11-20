// This is a main entry point to configure a WC/LS job that applies
// noise filtering and signal processing to existing RawDigits.  The
// FHiCL is expected to provide the following parameters as attributes
// in the "params" structure.
//
// epoch: the hardware noise fix expoch: "before", "after", "dynamic" or "perfect"
// reality: whether we are running on "data" or "sim"ulation.
// raw_input_label: the art::Event inputTag for the input RawDigit
//
// see the .fcl of the same name for an example Version 2
//
// Manual testing, eg:
//
// jsonnet -V reality=data -V epoch=dynamic -V raw_input_label=daq \\
//         -V signal_output_form=sparse \\
//         -J cfg cfg/pgrapher/experiment/uboone/wcls-nf-sp.jsonnet
//
// jsonnet -V reality=sim -V epoch=perfect -V raw_input_label=daq \\
//         -V signal_output_form=sparse \\
//         -J cfg cfg/pgrapher/experiment/uboone/wcls-nf-sp.jsonnet


local epoch = std.extVar('epoch');  // eg "dynamic", "after", "before", "perfect"
local reality = std.extVar('reality');
local sigoutform = std.extVar('signal_output_form');  // eg "sparse" or "dense"

local wc = import 'wirecell.jsonnet';
local g = import 'pgraph.jsonnet';

local raw_input_label = std.extVar('raw_input_label');  // eg "daq"
local volume_label = std.extVar('tpc_volume_label');  // eg "",0,1,2,3
local volume = if volume_label == '' then -1 else std.parseInt(volume_label);

// local data_params = import 'params.jsonnet';
// local simu_params = import 'simparams.jsonnet';
// local params_init = if reality == 'data' then data_params else simu_params;
local params_twofaced = import 'pgrapher/experiment/icarus/params_twofaced.jsonnet';

# Load the sim-params, overwrite the volume config for the two-faced version
local base_params = import 'pgrapher/experiment/icarus/simparams.jsonnet';
local base = base_params + params_twofaced; 

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
  files: super.files {
    fields: [ std.extVar('files_fields'), ],
    chresp: null,
  },

  rc_resp: if std.extVar('file_rcresp') != "" then
  {
    // "icarus_fnal_rc_tail.json"
    filename: std.extVar('file_rcresp'),
    postgain: 1.0,
    start: 0.0,
    tick: 0.4*wc.us,
    nticks: params.daq.nticks,// 4255,
    type: "JsonElecResponse",
  }
  else super.rc_resp,

  elec: std.mapWithIndex(function (n, eparam)
    super.elec[n] + {
      gain: eparam.gain,
      shaping: eparam.shaping,
    }, er_params),

};

// local tools_maker = import 'pgrapher/common/tools.jsonnet';
local tools_maker = import 'pgrapher/experiment/icarus/icarus_tools.jsonnet';
local tools = tools_maker(params);


local wcls_maker = import 'pgrapher/ui/wcls/nodes.jsonnet';
local wcls = wcls_maker(params, tools);

local sp_maker = import 'pgrapher/experiment/icarus/sp.jsonnet';

//local chndbm = chndb_maker(params, tools);
//local chndb = if epoch == "dynamic" then chndbm.wcls_multi(name="") else chndbm.wct(epoch);


// Collect the WC/LS input converters for use below.  Make sure the
// "name" argument matches what is used in the FHiCL that loads this
// file.  In particular if there is no ":" in the inputer then name
// must be the emtpy string.
local wcls_input = {
  adc_digits: g.pnode({
    type: 'wclsRawFrameSource',
    name: 'rfsrc%d' %volume, // to use multiple wirecell instances in a fhicl job
    data: {
      art_tag: raw_input_label,
      frame_tags: ['orig'],  // this is a WCT designator
      tick: params.daq.tick,
      // nticks: params.daq.nticks,
    },
  }, nin=0, nout=1),

};

// Collect all the wc/ls output converters for use below.  Note the
// "name" MUST match what is used in theh "outputers" parameter in the
// FHiCL that loads this file.

local this_anode = tools.anodes[volume];

local wcls_output = {
  // The noise filtered "ADC" values.  These are truncated for
  // art::Event but left as floats for the WCT SP.  Note, the tag
  // "raw" is somewhat historical as the output is not equivalent to
  // "raw data".
  nf_digits: g.pnode({
    type: 'wclsFrameSaver',
    name: 'nfsaver',
    data: {
      // anode: wc.tn(tools.anode),
      anode: wc.tn(this_anode),
      digitize: true,  // true means save as RawDigit, else recob::Wire
      frame_tags: ['raw'],
      // nticks: params.daq.nticks,
      chanmaskmaps: ['bad'],
    },
  }, nin=1, nout=1, uses=[this_anode]),


  // The output of signal processing.  Note, there are two signal
  // sets each created with its own filter.  The "gauss" one is best
  // for charge reconstruction, the "wiener" is best for S/N
  // separation.  Both are used in downstream WC code.
  sp_signals: g.pnode({
    type: 'wclsFrameSaver',
    name: 'spsaver%d' %volume, // to use multiple wirecell instances in a fhicl job
    data: {
      // anode: wc.tn(tools.anode),
      anode: wc.tn(this_anode),
      digitize: false,  // true means save as RawDigit, else recob::Wire
      // frame_tags: ['gauss', 'wiener', 'looseLf','shrinkROI','extendROI'],
      // frame_scale: [0.1, 0.1, 0.1],
      // frame_tags: ['gauss','wiener','looseLf','shrinkROI','extendROI','mp3ROI','mp2ROI', 'cleanupROI'],
      // frame_scale: [0.009,0.009,0.009,0.009,0.009,0.009,0.009,0.009],

      frame_tags: ['dnnsp'],
      frame_scale: [std.extVar('gain_ADC_per_e')],

      // nticks: params.daq.nticks,
      chanmaskmaps: [],
      nticks: -1,
    },
  }, nin=1, nout=1, uses=[this_anode]),

  h5io: g.pnode({
      type: 'HDF5FrameTap',
      name: 'hio_sp%d' % volume,
      data: {
        anode: wc.tn(this_anode),
        trace_tags: ['gauss'
	, 'wiener'
	, 'tightLf'
	, 'looseLf'
	, 'decon'
        , 'cleanupROI'
        , 'breakROI1'
        , 'breakROI2'
        , 'shrinkROI'
        , 'extendROI'
	, 'mp3ROI' 
	, 'mp2ROI' 
        , 'dnnsp'
        ],
        filename: "wc-sp-%d.h5" % volume,
        chunk: [0, 0], // ncol, nrow
        gzip: 2,
        high_throughput: true,
      },
    }, nin=1, nout=1, uses=[this_anode]),

};

// local perfect = import 'chndb-perfect.jsonnet';
local base = import 'pgrapher/experiment/icarus/chndb-base.jsonnet';
local chndb = [{
  type: 'OmniChannelNoiseDB',
  name: 'ocndbperfect%d' % n,
  // data: perfect(params, tools.anodes[n], tools.field, n),
  data: base(params, tools.anodes[n], tools.field, n){dft:wc.tn(tools.dft)},
  uses: [tools.anodes[n], tools.field, tools.dft],  // pnode extension
} for n in std.range(0, std.length(tools.anodes) - 1)];

local nf_maker = import 'pgrapher/experiment/icarus/nf.jsonnet';
local nf_pipes = [nf_maker(params, tools.anodes[n], chndb[n], tools, name='nf%d' % n) for n in std.range(0, std.length(tools.anodes) - 1)];

local sp_override = { // assume all tages sets in base sp.jsonnet
    sparse: sigoutform == 'sparse',
    use_roi_refinement: false,
    use_roi_debug_mode: false,
    wiener_tag: "",
    // gauss_tag: "",
    // tight_lf_tag: "",
    loose_lf_tag: "",
    break_roi_loop1_tag: "",
    break_roi_loop2_tag: "",
    shrink_roi_tag: "",
    extend_roi_tag: "",
    m_decon_charge_tag: "",
    cleanup_roi_tag: "",
    // mp2_roi_tag: "",
    // mp3_roi_tag: "",
    use_multi_plane_protection: true,
    mp_tick_resolution: 8,
    process_planes: [0, 1, 2],
    isWrapped: true,
    nwires_separate_planes: [
      [1056, 1056], [5600], [5600]
    ]
};
local sp = sp_maker(params, tools, sp_override);
local sp_pipes = [sp.make_sigproc(a) for a in tools.anodes];

local util = import 'pgrapher/experiment/icarus/funcs.jsonnet';
local chsel_pipes = [
  g.pnode({
    type: 'ChannelSelector',
    name: 'chsel%d' % n,
    data: {
      channels: util.anode_channels_twofaced(n),
      //tags: ['orig%d' % n], // traces tag
    },
  }, nin=1, nout=1)
  for n in std.range(0, std.length(tools.anodes) - 1)
];

local magoutput = 'icarus-data-check.root';
local magnify = import 'pgrapher/experiment/icarus/magnify-sinks.jsonnet';
local magnifyio = magnify(tools, magoutput);

local dnnroi = import 'pgrapher/experiment/icarus/dnnroi.jsonnet';
local ts_u = {
    type: "TorchService",
    name: "dnnroi_u",
    data: {
        model: "NNs/UNet-optfilter-opaqueMC_U_Plane.ts",
        device: "cpu",
        concurrency: 1,
    },
};

local ts_v = {
    type: "TorchService",
    name: "dnnroi_v",
    data: {
        model: "NNs/UNet-optfilter-opaqueMC_V_Plane.ts",
        device: "cpu",
        concurrency: 1,
    },
};

local nfsp_pipes = [
  g.pipeline([
               chsel_pipes[n],
               // magnifyio.orig_pipe[n],

               // nf_pipes[n],
               // magnifyio.raw_pipe[n],
               sp_pipes[n],
               dnnroi(tools.anodes[n], ts_u, ts_v, output_scale=1),
               // magnifyio.decon_pipe[n],
               // magnifyio.threshold_pipe[n],
               // magnifyio.debug_pipe[n], // use_roi_debug_mode: true in sp.jsonnet
             ],
             'nfsp_pipe_%d' % n)
  for n in [volume]
];

local fanout_tag_rules = [ 
          {
            frame: {
              '.*': 'orig%d' % tools.anodes[n].data.ident,
            },
            trace: {
              // fake doing Nmult SP pipelines
              //orig: ['wiener', 'gauss'],
              //'.*': 'orig',
            },
          }
          for n in [volume]
        ];

local fanin_tag_rules = [
          {
            frame: {
              //['number%d' % n]: ['output%d' % n, 'output'],
              '.*': 'framefanin',
            },
            trace: {
              ['extend_roi%d'%ind]:'extend_roi%d'%ind,
              ['shrink_roi%d'%ind]:'shrink_roi%d'%ind,
          //    ['break_roi_2nd%d'%ind]:'break_roi_2nd%d'%ind,
          //    ['break_roi_1st%d'%ind]:'break_roi_1st%d'%ind,
              ['cleanup_roi%d'%ind]:'cleanup_roi%d'%ind,
              ['mp2_roi%d'%ind]:'mp2_roi%d'%ind,
              ['mp3_roi%d'%ind]:'mp3_roi%d'%ind,
              ['gauss%d'%ind]:'gauss%d'%ind,
              ['dnnsp%d'%ind]:'dnnsp%d'%ind,
              ['wiener%d'%ind]:'wiener%d'%ind,
            //  ['threshold%d'%ind]:'threshold%d'%ind,
            //  ['tight_lf%d'%ind]:'tight_lf%d'%ind,
              ['loose_lf%d'%ind]:'loose_lf%d'%ind,
             // ['decon%d'%ind]:'decon%d'%ind,
            },

          }
          for ind in [this_anode.data.ident]
        ];
// local fanpipe = util.fanpipe('FrameFanout', nfsp_pipes, 'FrameFanin', 'nfsp', [], fanout_tag_rules, fanin_tag_rules);

local retagger = g.pnode({
  type: 'Retagger',
  data: {
    // Note: retagger keeps tag_rules an array to be like frame fanin/fanout.
    tag_rules: [{
      // Retagger also handles "frame" and "trace" like fanin/fanout
      // merge separately all traces like gaussXYZ to gauss.
      frame: {
        '.*': 'retagger',
      },
      merge: {
        'dnnsp\\d': 'dnnsp',
        'gauss\\d': 'gauss',
        'wiener\\d': 'wiener',
       // 'tight_lf\\d': 'tightLf',
        'loose_lf\\d': 'looseLf',
       // 'decon\\d': 'decon',
        'cleanup_roi\\d': 'cleanupROI',
       // 'break_roi_1st\\d': 'breakROI1',
       // 'break_roi_2nd\\d': 'breakROI2',
        'shrink_roi\\d': 'shrinkROI',
        'extend_roi\\d': 'extendROI',
	'mp3_roi\\d': 'mp3ROI',
	'mp2_roi\\d': 'mp2ROI',
      },
    }],
  },
}, nin=1, nout=1);

local sink = g.pnode({ type: 'DumpFrames' }, nin=1, nout=0);

local graph = g.pipeline([wcls_input.adc_digits, chsel_pipes[volume], sp_pipes[volume], dnnroi(this_anode, ts_u, ts_v, output_scale=1), retagger, wcls_output.h5io, wcls_output.sp_signals, sink]);

local app = {
  type: 'Pgrapher',
  data: {
    edges: g.edges(graph),
  },
};

// Finally, the configuration sequence
g.uses(graph) + [app]
