local g = import 'pgraph.jsonnet';
local f = import 'pgrapher/common/funcs.jsonnet';
local wc = import 'wirecell.jsonnet';

local io = import 'pgrapher/common/fileio.jsonnet';
local tools_maker = import 'pgrapher/experiment/icarus/icarus_tools.jsonnet';
local base = import 'pgrapher/experiment/icarus/simparams.jsonnet';
local splat = import 'pgrapher/experiment/icarus/splat.jsonnet';

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


local wcls_maker = import 'pgrapher/ui/wcls/nodes.jsonnet';
local wcls = wcls_maker(params, tools);

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

local anode_names = ["EES", "EEN", "EWS", "EWN", "WES", "WEN", "WWS", "WWN"];

local drifter = sim.drifter;
local setdrifter = g.pnode({
            type: 'DepoSetDrifter',
            data: {
                drifter: "Drifter"
            }
        }, nin=1, nout=1,
        uses=[drifter]);

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
      sed_label: 'largeant:TPCActive',

      // output to art::Event
      simchan_label: 'simpleSC',
    },
  },   nin=1, nout=1, uses=tools.anodes+[tools.field]);

local util = import 'pgrapher/experiment/icarus/funcs.jsonnet';

local deposplats = [splat(params, tools, tools.anodes[n], name="A%d" % n) for n in anode_iota] ;
local makesplat = util.fanpipe("DepoSetFanout", deposplats, "FrameFanin");

local sink = sim.frame_sink;

local graph = g.pipeline([wcls_input.deposet, setdrifter, wcls_simchannel_sink, makesplat, sink]);

local app = {
  type: 'Pgrapher',
  data: {
    edges: g.edges(graph),
  },
};


// Finally, the configuration sequence which is emitted.

g.uses(graph) + [app]
