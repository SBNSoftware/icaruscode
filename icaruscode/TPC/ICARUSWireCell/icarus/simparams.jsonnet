// Here we override params.jsonnet to provide simulation-specific params.

local base = import 'pgrapher/experiment/icarus/params.jsonnet';
local wc = import 'wirecell.jsonnet';

base {
  lar: super.lar {
    // Longitudinal diffusion constant
    DL :  4.0 * wc.cm2/wc.s,
    // Transverse diffusion constant
    DT :  8.8 * wc.cm2/wc.s,
    // Electron lifetime
    lifetime : 3.5*wc.ms,
    // Electron drift speed, assumes a certain applied E-field
    drift_speed : 1.5756*wc.mm/wc.us, // at 500 V/cm
  },
  daq: super.daq {

    // Number of readout ticks.  See also sim.response.nticks.
    // In MB LArSoft simulation, they expect a different number of
    // ticks than acutal data.
    // nticks: 4096,
  },

  sim: super.sim {

    // For running in LArSoft, the simulation must be in fixed time mode.
    fixed: true,
    continuous: false,
    fluctuate: true,

    //ductor : super.ductor {
    //    start_time: $.daq.start_time - $.elec.fields.drift_dt + $.trigger.time,
    //},

  },

  //files: super.files{
  //    chresp: null,
  //},

  sys_status: false,
  sys_resp: {
    // overall_short_padding should take into account this offset "start".
    start: -10 * wc.us,
    magnitude: 1.0,
    time_smear: 1.0 * wc.us,
  },

   rc_resp: {
     width: 1.1*wc.ms,
     rc_layers: 1,
   }
}
