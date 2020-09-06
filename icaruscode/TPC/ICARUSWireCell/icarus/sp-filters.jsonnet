// WARNING: the SP C++ code has a lot of hard coded names for various
// filter components.  Until this is cleaned up, one MUST configure
// the filter competents with matching type names and not change their
// instance names.


local wc = import 'wirecell.jsonnet';

local lf(name, data={}) = {
  type: 'LfFilter',
  name: name,
  data: {
    max_freq: 1 * wc.megahertz,
    tau: 0.0 * wc.megahertz,
  } + data,
};
local hf(name, data={}) = {
  type: 'HfFilter',
  name: name,
  data: {
    max_freq: 1 * wc.megahertz,
    sigma: 0.0 * wc.megahertz,
    power: 2,
    flag: true,
  } + data,
};
// All "wire" filters are Hf with different base values.
local wf(name, data={}) = {
  type: 'HfFilter',
  name: name,
  data: {
    max_freq: 1,  // warning: units
    power: 2,
    flag: false,
    sigma: 0.0,  // caller should provide
  } + data,
};

// Zeus take my eyes! Magic numbers are everywhere!
[
  lf('ROI_tight_lf', { tau: 0.014 * wc.megahertz }),  // 0.02 
  lf('ROI_tighter_lf', { tau: 0.06 * wc.megahertz }),  // 0.1 
  lf('ROI_loose_lf', { tau: 0.0025 * wc.megahertz }),  // 0.0025 

  hf('Gaus_tight'),
  hf('Gaus_wide', { sigma: 0.12 * wc.megahertz }), 


  hf('Wiener_tight_U', {
    sigma: 0.148788  * wc.megahertz,
    power: 3.76194,
  }),
  hf("Wiener_tight_V", {
    sigma: 0.1596568 * wc.megahertz,
    power: 4.36125 }),
  hf('Wiener_tight_W', {
    sigma: 0.13623 * wc.megahertz,
    power: 3.35324,
  }),

  hf('Wiener_wide_U', {
    sigma: 0.186765  * wc.megahertz,
    power: 5.05429,
  }),
  hf("Wiener_wide_V", {
    sigma: 0.1936 * wc.megahertz,
    power: 5.77422,
  }),
  hf('Wiener_wide_W', {
    sigma: 0.175722  * wc.megahertz,
    power: 4.37928,
  }),

  wf('Wire_ind', { sigma: 1.0 / wc.sqrtpi * 0.75 }), 
  wf('Wire_col', { sigma: 1.0 / wc.sqrtpi * 3.0 }),
]
