// This file specifies the paramter configuration for the ICARUS detector. It
// inherit from the base params.jsonnet and override the relevant parameters

local wc = import "wirecell.jsonnet";
local base = import "pgrapher/common/params.jsonnet";

base {
    det : {

        // define the 4 APAs.  This must use the coordinate system
        // defined by the wire geometry file.
        //
        // The "faces" is consumed by, at least, the Drifter and
        // AnodePlane.  The "wires" number is used to set
        // AnodePlane.ident used to lookup the anode in WireSchema.
        // It corresponds to the anode number.
        //
        // Also see:
        //   lar -c dump_icarus_geometry.fcl
        //   wirecell-util wire-volumes icarus-wires-dualanode.json.bz2
        // to help with defining these parameters.


        // The "response" plane is where the field response functions
        // start.  Garfield calcualtions start somewhere relative to
        // something, here's where that is made concrete.  This MUST
        // match what field response functions also used.
        response_plane: 10*wc.cm, // relative to collection wires


        // Each wire gets identified with the id number in the geometry
        // file. The horizontal induction is split into two, both
        // sides are enumerated by "s", while "n" counts the physical anodes.
        // NOTE: the actual physical volumes in ICARUS are only 4

        local xanode = [-369.33*wc.cm, -71.1*wc.cm, 71.1*wc.cm, 369.33*wc.cm],
        local offset_response = [if a%2==0 then +10*wc.cm else -10*wc.cm for a in std.range(0,3)],
        local xresponse = [xanode[a] + offset_response[a] for a in std.range(0,3)],
        local xcathode = [-220.29*wc.cm, -220.29*wc.cm, 220.29*wc.cm, 220.29*wc.cm],
        volumes : [
            {
                local world = 100,  // identify this geometry
                local split = s*10, // identify anode side (1 left, 2 right)
                local anode = a,    // physical anode number

                wires: (world+split+anode),
                name: "anode%d"%(world+split+anode),

                faces: [
                        {
                            anode: xanode[a],
                            response: xresponse[a],
                            cathode: xcathode[a],
                        },

                        null
                ],
            } for a in std.range(0,3) for s in std.range(1,2)
        ],

        // This describes some rough, overall bounding box.  It's not
        // directly needed but can be useful on the Jsonnet side, for
        // example when defining some simple kinematics.  It is
        // represented by a ray going from extreme corners of a
        // rectangular solid.  Again "wirecell-util wires-info" helps
        // to choose something. //FIXME -- ARE CORRECT?
        bounds : {
            tail: wc.point(-3.65, -1.7, -9.1, wc.m),
            head: wc.point(+3.65, +1.4, +8.8, wc.m),
        }
    },

    daq: super.daq {
        tick: 0.4*wc.us, // 2.5 MHz
        nticks: 4096,
    },

    adc: super.adc {
        // fix baseline at 2048 (induction), 400 (collection)
        baselines: [1650.0*wc.millivolt, 1650.0*wc.millivolt, 322.3*wc.millivolt],

        // From ICARUS paper: https://iopscience.iop.org/article/10.1088/1748-0221/13/12/P12007/pdf
        //check (values taken from the FE calibration shown in pg. 7 of the paper)
        fullscale: [0.8*wc.millivolt, 3.3*wc.volt],
    },

    elec: super.elec {
        type: "WarmElecResponse",
        gain: 17.8075*wc.mV/wc.fC, // 0.027 fC/(ADC*us)
        shaping: 1.3*wc.us,
        postgain: 1.0,
        start: 0,
    },


    sim: super.sim {

        // For running in LArSoft, the simulation must be in fixed time mode. 
        fixed: true,

        // The "absolute" time (ie, in G4 time) that the lower edge of
        // of final readout tick #0 should correspond to.  This is a
        // "fixed" notion.
        local tick0_time = -340*wc.us, // TriggerOffsetTPC from detectorclocks_icarus.fcl

        // Open the ductor's gate a bit early.
        local response_time_offset = $.det.response_plane / $.lar.drift_speed,
        local response_nticks = wc.roundToInt(response_time_offset / $.daq.tick),

        ductor : {
            nticks: $.daq.nticks + response_nticks,
            readout_time: self.nticks * $.daq.tick,
            start_time: tick0_time - response_time_offset,
        },

        // To counter the enlarged duration of the ductor, a Reframer
        // chops off the little early, extra time.  Note, tags depend on how 
        reframer: {
            tbin: response_nticks,
            nticks: $.daq.nticks,
        }
        
    },

    files: {
        wires: "icarus-wires-dualanode-v2.json.bz2",

        fields: ["garfield-icarus-fnal-commissioning.json.bz2"],

        noise: "sbn_fd_incoherent_noise.json.bz2",
        coherent_noise: "sbn_fd_coherent_noise.json.bz2",


        chresp: null,
    },

}

