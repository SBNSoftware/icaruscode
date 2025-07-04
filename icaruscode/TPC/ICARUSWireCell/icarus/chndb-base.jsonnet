// Base channel noise DB object configuration for microboone
// This does not include any run dependent RMS cuts.
// See chndb.jsonnet

local handmade = import 'chndb-resp.jsonnet';
local wc = import 'wirecell.jsonnet';
local util = import 'pgrapher/experiment/icarus/funcs.jsonnet';

function(params, anode, field, n, rms_cuts=[])
  {
    anode: wc.tn(anode),
    field_response: wc.tn(field),

    tick: params.daq.tick,

    // This sets the number of frequency-domain bins used in the noise
    // filtering.  It is not necessarily true that the time-domain
    // waveforms have the same number of ticks.  This must be non-zero.
    nsamples: params.nf.nsamples,

    // For MicroBooNE, channel groups is a 2D list.  Each element is
    // one group of channels which should be considered together for
    // coherent noise filtering.
    groups: [std.range(g*32, (g+1)*32-1) for g in std.range(0,171)],
    // groups: [std.range(n * 2560 + u * 40, n * 2560 + (u + 1) * 40 - 1) for u in std.range(0, 19)]
    //         + [std.range(n * 2560 + 800 + v * 40, n * 2560 + 800 + (v + 1) * 40 - 1) for v in std.range(0, 19)]
    //         + [std.range(n * 2560 + 1600 + w * 48, n * 2560 + 1600 + (w + 1) * 48 - 1) for w in std.range(0, 19)],


    // Externally determined "bad" channels.
    bad: [
        167, 1795, 1796, 1824, 1826, 1828, 1830, 1860, 1888, 1889, 1890, 1891, 1922, 2013, 2806, 3725, 3743, 4262, 5572, 8523, 8524, 8525, 8526, 8527, 8945, 9180, 9387, 9391, 9664, 10138, 12051, 14028, 14403, 15307, 15737, 16017, 16018, 16019, 16020, 16021, 16022, 16023, 16024, 16025, 16026, 16027, 16028, 16029, 16030, 16031, 16937, 17201, 19580, 23225, 23457, 23984, 23985, 25662, 27219, 27231, 27335, 27967, 28640, 28804, 28805, 28806, 28807, 28808, 28809, 28810, 30496, 30497, 30498, 30499, 30500, 30501, 30502, 30503, 30504, 30505, 30506, 30507, 30508, 30509, 30510, 30511, 30512, 30513, 30514, 30515, 30516, 30517, 30518, 30519, 30520, 30521, 30522, 30523, 30524, 30525, 30526, 30527, 30705, 31316, 31317, 31726, 34220, 38294, 41501, 41975, 42851, 43231, 43518, 44051, 44410, 49216, 49792, 49793, 49794, 49795, 49796, 49797, 49798, 49799, 49800, 49801, 49802, 49803, 49804, 49805, 49806, 49807, 49808, 49809, 49810, 49811, 49812, 49813, 49814, 49815, 49816, 49817, 49818, 49819, 49820, 49821, 49822, 49823, 50080, 50081, 50082, 50083, 50084, 50085, 50086, 50087, 50088, 50089, 50090, 50091, 50092, 50093, 50094, 50095, 50096, 50097, 50098, 50099, 50100, 50101, 50102, 50103, 50104, 50105, 50106, 50107, 50108, 50109, 50110, 50111, 50673, 51900, // v3r3
       1858, 1887, 1892, 3355, 2525, 6289, 7402, 47914, 18132, 18609, // added
    ],

    // Overide defaults for specific channels.  If an info is
    // mentioned for a particular channel in multiple objects in this
    // list then last mention wins.
    channel_info: [

      // First entry provides default channel info across ALL
      // channels.  Subsequent entries override a subset of channels
      // with a subset of these entries.  There's no reason to
      // repeat values found here in subsequent entries unless you
      // wish to change them.
      {
        channels: util.anode_channels_twofaced(n),
        nominal_baseline: 2048.0,  // adc count
        gain_correction: 1.0,  // unitless
        response_offset: 0.0,  // ticks?
        pad_window_front: 10,  // ticks?
        pad_window_back: 10,  // ticks?
        decon_limit: 0.02,
        decon_limit1: 0.09,
        adc_limit: 15,
        roi_min_max_ratio: 0.8, // default 0.8
        min_rms_cut: 1.0,  // units???
        max_rms_cut: 30.0,  // units???

        // parameter used to make "rcrc" spectrum
        rcrc: 1.1 * wc.millisecond, // 1.1 for collection, 3.3 for induction
        rc_layers: 1, // default 2

        // parameters used to make "config" spectrum
        reconfig: {},

        // list to make "noise" spectrum mask
        freqmasks: [],

        // field response waveform to make "response" spectrum.
        response: {},

      },

    ] + rms_cuts,
  }
