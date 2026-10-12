// Per-TPC variant of wcls-multitpc-sim-drift-simchannel-yzsim-refactored.jsonnet.
//
// Reads `tpc_idx` (0..3) from an external variable.  Builds only the slice of
// the graph that belongs to that one TPC by narrowing `tools.anodes` to the
// per-TPC two-element subset and letting sim_maker scope its derived per-
// anode pipelines (transformsyz, analog_pipelinesyz, etc.) accordingly.
//
// Designed to be driven by 4 separate `WireCellToolkit` art modules, one per
// TPC, so each module's `WCLS_tool::process()` flushes its TPC's frame to the
// art Event before the next module starts and the float SimpleTraces can be
// reclaimed.
//
// Per-instance component names that need to be unique across the 4 modules
// (Pgrapher, fanout, summer, actpipe, FrameSaver, etc.) carry `%d %tpc_idx`.
// The shared heavy components (PIR, FieldResponse, AnodePlane, DFT) stay
// singletons across modules - WCT's NamedFactory de-duplicates them by name
// and their `configure()` is now idempotent (PIR's `m_bywire`/`m_ir` are
// cleared at the top of `build_responses()`).
//
// Instance names that show up in the art Event are kept identical to the
// single-module setup (`simdigits<tpc_idx>` and `simpleSC<n>` with global n in
// [tpc_idx*90, tpc_idx*90+90)).

local g = import 'pgraph.jsonnet';
local f = import 'pgrapher/common/funcs.jsonnet';
local wc = import 'wirecell.jsonnet';

local io = import 'pgrapher/common/fileio.jsonnet';
local tools_maker = import 'pgrapher/experiment/icarus/icarus_tools.jsonnet';
local base = import 'pgrapher/experiment/icarus/simparams.jsonnet';

local er_params = [
  { gain: std.extVar('gain0')*wc.mV/wc.fC, shaping: std.extVar('shaping0')*wc.us },
  { gain: std.extVar('gain1')*wc.mV/wc.fC, shaping: std.extVar('shaping1')*wc.us },
  { gain: std.extVar('gain2')*wc.mV/wc.fC, shaping: std.extVar('shaping2')*wc.us },
];

local params = base {
  lar: super.lar {
    DL: std.extVar('DL') * wc.cm2 / wc.ns,
    DT: std.extVar('DT') * wc.cm2 / wc.ns,
    lifetime: std.extVar('lifetime') * wc.us,
  },
  files: super.files {
    fields: [
      "icarus_fnal_fit_ks_P0nom_P1bin0.json.bz2",
      "icarus_fnal_fit_ks_P0nom_P1bin1.json.bz2",
      "icarus_fnal_fit_ks_P0nom_P1bin2.json.bz2",
      "icarus_fnal_fit_ks_P0nom_P1bin3.json.bz2",
      "icarus_fnal_fit_ks_P0nom_P1bin4.json.bz2",
      "icarus_fnal_fit_ks_P0nom_P1bin5.json.bz2",
      "icarus_fnal_fit_ks_P0nom_P1bin6.json.bz2",
      "icarus_fnal_fit_ks_P0nom_P1bin7.json.bz2",
      "icarus_fnal_fit_ks_P0nom_P1bin8.json.bz2",
      "icarus_fnal_fit_ks_P0nom_P1bin9.json.bz2",
      "icarus_fnal_fit_ks_P0nom_P1bin10.json.bz2",
      "icarus_fnal_fit_ks_P0nom_P1bin11.json.bz2",
      "icarus_fnal_fit_ks_P0nom_P1bin12.json.bz2",
      "icarus_fnal_fit_ks_P0nom_P1bin13.json.bz2",
      "icarus_fnal_fit_ks_P0nom_P1bin14.json.bz2",
      "icarus_fnal_fit_ks_P0nom_P1bin15.json.bz2",
    ],
  },
  rc_resp: if std.extVar('file_rcresp') != "" then {
    filename: std.extVar('file_rcresp'),
    postgain: 1.0,
    start: 0.0,
    tick: 0.4*wc.us,
    nticks: 4255,
    type: "JsonElecResponse",
    rc_layers: 1,
  } else super.rc_resp,
  elec: std.mapWithIndex(function (n, eparam)
    super.elec[n] + { gain: eparam.gain, shaping: eparam.shaping }, er_params),
};

// ---------------------------------------------------------------------------
// Per-TPC selection.
// tools_all is the full 8-anode tools object; `tools` is the same object with
// anodes narrowed to this TPC's pair.  sim_maker then builds depotransforms
// and analog pipelines only for those anodes (×45 per anode = 90 entries).
// ---------------------------------------------------------------------------
local tpc_idx = std.parseInt(std.extVar('tpc_idx'));
local apa_lo = tpc_idx * 2;     // 2 anodes per TPC
local local_iota = std.range(0, 89);  // 90 per-anode-plane-response branches
local volname = ["EE", "EW", "WE", "WW"];

local tools_all = tools_maker(params);
local tools = tools_all { anodes: tools_all.anodes[apa_lo : apa_lo + 2] };

local sim_maker = import 'pgrapher/experiment/icarus/sim.jsonnet';
local sim = sim_maker(params, tools);   // analog_pipelinesyz has 90 entries

// ---------------------------------------------------------------------------
// Input: depo source from art Event (identical across all 4 per-TPC modules)
// ---------------------------------------------------------------------------
local wcls_maker = import 'pgrapher/ui/wcls/nodes.jsonnet';
local wcls = wcls_maker(params, tools);
local wcls_input = {
  deposet: g.pnode({
    type: 'wclsSimDepoSetSource',
    name: "electron",
    data: {
      model: "",
      scale: -1,
      art_tag: std.extVar('SimEnergyDepositLabel'),
      assn_art_tag: "",
      id_is_track: false,
    },
  }, nin=0, nout=1),
};

// ---------------------------------------------------------------------------
// MegaAnode / DuoAnode for this TPC.  `mega_anode` here is a *per-TPC*
// MegaAnodePlane (not the shared 8-anode one) - it covers exactly the two
// anodes this TPC processes.
// ---------------------------------------------------------------------------
local mega_anode = {
  type: 'MegaAnodePlane',
  name: 'meganodes%d' % tpc_idx,
  data: { anodes_tn: [wc.tn(a) for a in tools.anodes] },
};

// duoanode is here for compatibility with wclsFrameSaver, which also takes
// this TPC's pair of anodes.
local duoanode = mega_anode;

// ---------------------------------------------------------------------------
// Output: this TPC's RawDigit saver (instance name stays simdigits<tpc_idx>)
// ---------------------------------------------------------------------------
local sim_digits = g.pnode({
  type: 'wclsFrameSaver',
  name: 'simdigits%d' % tpc_idx,
  data: {
    anode: wc.tn(duoanode),
    digitize: true,
    frame_tags: ['TPC%s' % volname[tpc_idx]],
  },
}, nin=1, nout=1, uses=[duoanode]);

// ---------------------------------------------------------------------------
// Helpers: global index for the n-th branch of this TPC.
// ---------------------------------------------------------------------------
local n_lo = tpc_idx * 90;
local glob(i) = n_lo + i;

// ---------------------------------------------------------------------------
// 90 drifters / scalers / simchannel sinks for this TPC.
// Instance names keep the global numbering (n_lo..n_hi-1) so SimChannel
// product instances match the single-module reference run.
// ---------------------------------------------------------------------------
local overlay_drifter = std.extVar("overlay_drifter");
local localeLiftime = [std.extVar('lifetime') * wc.us,std.extVar('lifetime') * wc.us,std.extVar('lifetime') * wc.us,std.extVar('lifetime') * wc.us,std.extVar('lifetime') * wc.us,std.extVar('lifetime') * wc.us,std.extVar('lifetime') * wc.us,std.extVar('lifetime') * wc.us];
local drifter_data = if overlay_drifter then sim.overlay_drifter_data else sim.drifter_data;

local drifters = [{
  type: if overlay_drifter then "wclsICARUSDrifter" else "Drifter",
  name: "drifter%d" % glob(i),
  data: params.lar + drifter_data {
    lifetime: localeLiftime[std.floor(glob(i)/45)],
    TPC: tpc_idx,
    charge_scale: 1,
  },
} for i in local_iota];

local setdrifters = [g.pnode({
  type: 'DepoSetDrifter',
  name: 'setdrifters%d' % glob(i),
  data: { drifter: wc.tn(drifters[i]) },
}, nin=1, nout=1, uses=[drifters[i]]) for i in local_iota];

local yzmap_filter = {
  type: "YZMap", name: "yzmap_filter",
  data: {
    filename: 'yzmap_icarus_v3_run1.json',
    bin_width: 10*wc.cm, bin_height: 10*wc.cm,
    yoffset: 180*wc.cm, zoffset: 900*wc.cm,
    nbinsy: 31, nbinsz: 180,
  },
};

local yzmap_gain = {
  type: "YZMap", name: "yzmap_gain",
  data: {
    filename: 'yzmap_gain_icarus_v3_run1.json',
    bin_width: 10*wc.cm, bin_height: 10*wc.cm,
    yoffset: 180*wc.cm, zoffset: 900*wc.cm,
    nbinsy: 31, nbinsz: 180,
  },
};

local scalers = [{
  type: "Scaler",
  name: "scaler%d" % glob(i),
  data: {
    yzmap: wc.tn(yzmap_gain),
    anode: wc.tn(tools.anodes[std.floor(i/45)]),
    plane: std.mod(std.floor(i/15), 3),
  },
} for i in local_iota];

local setscaler = [g.pnode({
  type: 'DepoSetScaler',
  name: 'setscaler%d' % glob(i),
  data: { scaler: wc.tn(scalers[i]) },
}, nin=1, nout=1, uses=[scalers[i], yzmap_gain]) for i in local_iota];

local wcls_simchannel_sink = [g.pnode({
  type: 'wclsDepoFluxWriter',
  name: 'postdrift%d' % glob(i),
  data: {
    anodes: [wc.tn(a) for a in tools.anodes],
    field_response: wc.tn(tools.field),
    tick: params.daq.tick,
    window_start: -340 * wc.us,
    window_duration: params.daq.readout_time,
    nsigma: 3.0,
    reference_time: -1500 * wc.us - self.window_start,
    smear_long: 0.0,
    smear_tran: 0.0,
    time_offsets: [
      std.extVar('time_offset_u') * wc.us,
      std.extVar('time_offset_v') * wc.us,
      std.extVar('time_offset_y') * wc.us,
    ],
    process_planes: [std.mod(std.floor(i/15), 3)],
    sed_label: std.extVar('SimEnergyDepositLabel'),
    simchan_label: 'simpleSC%d' % glob(i),
  },
}, nin=1, nout=1, uses=tools.anodes + [tools.field]) for i in local_iota];

local deposetfilteryz = [g.pnode({
  type: 'DepoSetFilterYZ',
  name: 'deposetfilteryz_resp%d-plane%d-%s' % [
    std.mod(i, 15),
    std.mod(std.floor(i/15), 3),
    tools.anodes[std.floor(i/45)].name,
  ],
  data: {
    yzmap: wc.tn(yzmap_filter),
    resp: std.mod(i, 15),
    anode: wc.tn(tools.anodes[std.floor(i/45)]),
    plane: std.mod(std.floor(i/15), 3),
  },
}, nin=1, nout=1, uses=tools.anodes + [yzmap_filter]) for i in local_iota];

// ---------------------------------------------------------------------------
// 90 analog pipelines for this TPC, sourced from the narrowed sim object.
// ---------------------------------------------------------------------------
local sigpipes = sim.analog_pipelinesyz;
assert std.length(sigpipes) == 90 :
       "expected 90 analog pipelines for one TPC, got %d" % std.length(sigpipes);

// ---------------------------------------------------------------------------
// 1 FrameSummerYZ for this TPC (sums 90 → 1)
// ---------------------------------------------------------------------------
local frame_summer = g.pnode({
  type: 'FrameSummerYZ',
  name: 'framesummer%d' % tpc_idx,
  data: { multiplicity: 90 },
}, nin=90, nout=1);

// ---------------------------------------------------------------------------
// Noise + digitizer for this TPC (one of each)
// ---------------------------------------------------------------------------
local nicks = ["incoTPCEE","incoTPCEW","incoTPCWE","incoTPCWW",
               "coheTPCEE","coheTPCEW","coheTPCWE","coheTPCWW"];
local scale_int = std.extVar('int_noise_scale');
local scale_coh = std.extVar('coh_noise_scale');

local inco_model = {
  type: "GroupNoiseModel",
  name: nicks[tpc_idx],
  data: {
    spectra: params.files.noisegroups[tpc_idx],
    groups: params.files.wiregroups,
    scale: scale_int,
    nsamples: params.daq.nticks,
    tick: params.daq.tick,
  },
};

local cohe_model = {
  type: "GroupNoiseModel",
  name: nicks[tpc_idx + 4],
  data: {
    spectra: params.files.noisegroups[tpc_idx + 4],
    groups: params.files.wiregroups,
    scale: scale_coh,
    nsamples: params.daq.nticks,
    tick: params.daq.tick,
  },
};

local add_noise_node = function(model, t) g.pnode({
  type: t,
  name: "addnoise%d-" % tpc_idx + model.name,
  data: {
    rng: wc.tn(tools.random),
    dft: wc.tn(tools.dft),
    model: wc.tn(model),
    nsamples: params.daq.nticks,
  },
}, nin=1, nout=1, uses=[tools.random, tools.dft, model]);

local inco_noise = add_noise_node(inco_model, "IncoherentAddNoise");
local cohe_noise = add_noise_node(cohe_model, "CoherentAddNoise");

local digitizer = sim.digitizer(mega_anode,
                                name="digitizer%d-" % tpc_idx + mega_anode.name,
                                tag="TPC%s" % volname[tpc_idx]);

local reframer = g.pnode({
  type: 'Reframer',
  name: 'reframer-%d-' % tpc_idx + mega_anode.name,
  data: {
    anode: wc.tn(mega_anode),
    tags: "TPC%s" % volname[tpc_idx],
    fill: 0.0,
    tbin: params.sim.reframer.tbin,
    toffset: 0,
    nticks: params.sim.reframer.nticks,
  },
}, nin=1, nout=1);

local actpipe = g.pipeline(
  [reframer, inco_noise, cohe_noise, digitizer, sim_digits],
  name="noise-digitizer%d" % tpc_idx);

// ---------------------------------------------------------------------------
// Drift-side pipelines: one per anode/plane (filter → drifter → scaler → simchan-sink)
// ---------------------------------------------------------------------------
local driftpipes = [g.pipeline(
  [deposetfilteryz[i], setdrifters[i], setscaler[i], wcls_simchannel_sink[i]],
  name="depo-set-drifter%d" % glob(i)) for i in local_iota];

// ---------------------------------------------------------------------------
// Top-level per-TPC topology:
//   wcls_input.deposet
//      → fanout (1 → 90)
//      → driftpipes[i]  (90×; each ends in simchan-sink which writes to art)
//      → sigpipes[i]    (90× DepoTransform → analog frame)
//      → frame_summer   (90 → 1)
//      → actpipe        (reframer → noise → cohe-noise → digitizer → wclsFrameSaver)
//      → final DumpFrames sink
// ---------------------------------------------------------------------------
local fanout = g.pnode({
  type: 'DepoSetFanout',
  name: 'fandrifter%d' % tpc_idx,
  data: { multiplicity: 90, tag_rules: [] },
}, nin=1, nout=90);

local actpipe_sink = g.pnode(
  { type: 'DumpFrames', name: 'actpipe%d-dump' % tpc_idx },
  nin=1, nout=0);

local drift = g.intern(
  innodes=driftpipes,
  outnodes=[actpipe],
  centernodes=sigpipes + [frame_summer],
  edges=
    [g.edge(driftpipes[i], sigpipes[i]) for i in local_iota] +
    [g.edge(sigpipes[i], frame_summer, 0, i) for i in local_iota] +
    [g.edge(frame_summer, actpipe)],
  name='drift-tpc%d' % tpc_idx);

local fandrifter = g.intern(
  innodes=[fanout],
  outnodes=[actpipe_sink],
  centernodes=[drift],
  edges=
    [g.edge(fanout, driftpipes[i], i, 0) for i in local_iota] +
    [g.edge(drift, actpipe_sink, 0, 0)],
  name='fandrifter%d-top' % tpc_idx);

local graph = g.pipeline([wcls_input.deposet, fandrifter]);

// ---------------------------------------------------------------------------
// Per-TPC Pgrapher name.  WCT's NamedFactory is a singleton across art
// modules; naming per-TPC keeps each module's graph isolated.  The fcl's
// `apps:` for daq<k> must reference "Pgrapher:pgrapher<k>" to match.
// ---------------------------------------------------------------------------
local app = {
  type: 'Pgrapher',
  name: 'pgrapher%d' % tpc_idx,
  data: { edges: g.edges(graph) },
};

// ---------------------------------------------------------------------------
// Memoized DFS post-order for graph.uses (same as the upstream jsonnet).
// ---------------------------------------------------------------------------
local key_of(x) = if std.objectHas(x, 'name') && x.name != ''
                  then x.type + ':' + x.name
                  else x.type;
local strip_uses(x) = { [k]: x[k] for k in std.objectFields(x) if k != 'uses' };

local visit(node, state) =
  if node.type == 'Pnode' then
    if std.objectHas(node, 'uses')
    then std.foldl(function(s, c) visit(c, s), node.uses, state)
    else state
  else
    local k = key_of(node);
    if std.objectHas(state.seen, k) then state
    else
      local after_children =
        if std.objectHas(node, 'uses')
        then std.foldl(function(s, c) visit(c, s), node.uses, state)
        else state;
      {
        seen: after_children.seen { [k]: true },
        result: after_children.result + [strip_uses(node)],
      };

local resolve_uses_unique(roots) =
  std.foldl(function(s, n) visit(n, s), roots, { seen: {}, result: [] }).result;

resolve_uses_unique(graph.uses) + [app]
