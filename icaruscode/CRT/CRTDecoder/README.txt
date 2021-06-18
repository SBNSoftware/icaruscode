/////////////////////////////////////////////////////////////
//  Brief overview of the CRTDecoder code
//  Chris.Hilgenberg@colostate.edu
//  last updated 2 July 2020
//
/////////////////////////////////////////////////////////////

    //                 DAQ Overview                           //
    ////////////////////////////////////////////////////////////
    All code contained in /icaruscode/CRT/CRTDecoder is designed for
    analysis and processing of raw CRT data files acquired using a
    pseudo-production version of the CRT DAQ system based on artDAQ/sbndaq.

    Currently, DAQInterface is not functional. This software provides a central
    DAQ configuration/control interface for managing multiple front-end board
    network daisychains spanning multiple DAQ servers. The data streams from the
    individual daisychains are forwarded to the event builder server for
    pre-processing and merging with other detector subsystem data streams (e.g.
    PMTs).

    In the current setup, we must deal with primitve DAQ processes known as
    artdaqDriver. This driver handles communication with the front-end board
    daisychains with one instance per chain including transmission of
    configuration bitstreams and receipt of data packets. The operator must
    instantiate each instance manually.

    The current CRT system has 2 walls (North and West Rolling) with 2 layers,
    each with its own daisychain. Each wall is controlled by one DAQ server with
    two network ports each.

    artdaqDriver each data fragment, corresponding to one front-end board
    trigger/readout is saved in a single art::event, thus one event contains
    that single fragment only. The events are ordered such that fragments for
    a particular front-end board are time ordered. However, the data packets are
    received via a single poll of the DAQ containing all fragment series for each
    of the front-end boards on the daisychain. Thus, accross a poll, the data is
    not time-ordered.

    DAQInterface (which is how we run DAQ in Control Room) in turn collects
    multiple fragments from given FEB in fragment containers. An event will
    therefore contain multiple fragment containers, and each fragment container
    will contain zero or multiple fragments.

    The data collected by artdaqDriver and DAQInterface therefore has different
    structure and must be read differently. In the programs described below
    you can find function analyze() which distinguishes between fragments and
    fragment containers in each event and properly loops over them.

     -- poll structure with one time ordered sequence per FEB --

     ^      -     -     -     -
     |     -     -     -     -
     |    -     -     -     -
  t  |   -     -     -     -
  i  |  -     -     -     -
  m  | -     -     -     -
  e  |________________________>
       event number

    A fragment header file provides the accessors we need to extract the
    fields and metadata from the individual fragments.



    //              CRTDecoder workflows                     //
    ///////////////////////////////////////////////////////////
    There are different possible workflows currently supported
     A. raw data analysis
       1. convert raw data files into ntuples
       2. Uses single analyzer module: BernCRTAna_module.cc
       3. Run with "lar -c analyze_BernCRT.fcl -s data.artroot -T ntuple.root"
     B. noise monitoring and analysis
       1. sepecialized version of (A) for data where SiPM HV is off
       2. Uses single analyzer module: CrtNoiseMonTool_module.cc
       3. Run with "lar -c crt_noise_mon.fcl -s data.artroot -T ntuple.root"
       4. Handy ROOT analysis macro at ./Macros/noise_ana.cxx
     C. calibration and analysis
       1. determine active channels, pedestal values, gains
       2. Calibration algroithms live in dedicated class, CrtCal.h
       3. CrtCal applied in CrtCalAnalysis_module.cc which outputs ntuple
       4. Set daisychain to be analyzed in crt_cal_ana.fcl
       5. Run with "lar -c crt_cal_ana.fcl -s cal_date_region.artroot -T cal_date_regon_ntuple.root"
       6. For convenince, after all daisychains are processed, merge the
          seperate calntuples into one file. E.g.
          "hadd cal_date_merge.root cal_date_region1_in.root cal_date_region1_out.root ..."
       7. You can find a convenint ROOT macro for calibration analysis in ./Macros/CalAna.C
     D. process cosmics runs into CRTData data products
       1. first step in injecting our data into icaruscode
          reconstruction chains and analysis
       2. this is the most complicated workflow as we currently
          do not have an event builder for the CRT data
       3. The code supports performing a auto calibration, but we will not
          consider this here.
       4. The first step is to calibrate and extract the data we need from
          the raw files, facilitated by AnaProducer_module.cc
       5. We require a recent calibration ntuple before starting. Point our module to
          this file in  cal_ana.fcl
       6. Process each raw data file for each CRT layer/region using dedicated fcl files
          "lar -c crt_ana_cal_<region>_<layer>.fcl -s raw_cosmics_region_layer.artroot
            -T cosmics_region_layer_ntuple.root"
       7. We now need to merge the different pre-process cosmics ntuples into one file
          using "hadd cosmics_merge.root cosmics_region1_layer1.root ... cosmics_region2_layer2.root"
       8. Now time order the resulting fie. Use an executable for this,
          "./CRTMergePreProcessTrees.exe cosmics_merge.root cosmics_merge_timeordered.root"
       9. We will use this last file to create new art::events based on time slices of the file
          using the module, CRTEventProducer_module.cc
       10. The input file and time slice width are set in crt_event_prod.fcl
       11. To run: "lar -c prod_crt_event.fcl -o output.artroot"
       12. If everything went smoothly, output.root should be ready for further
           processing/analysis using existing icaruscode tools.

