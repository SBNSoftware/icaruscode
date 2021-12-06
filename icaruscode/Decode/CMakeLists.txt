cet_enable_asserts()

cet_find_library(LIBWDA NAMES wda PATHS ENV LIBWDA_LIB NO_DEFAULT_PATH)

art_make(
          EXCLUDE
                        PMTconfigurationExtraction_module.cc
                        DumpPMTconfiguration_module.cc
                        DumpArtDAQfragments_module.cc
                        DaqDecoderICARUSPMT_module.cc
          MODULE_LIBRARIES
                        icarus_signal_processing
                        icarus_signal_processing_Detection
                        icarus_signal_processing_Filters
                        sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_ICARUS 
                        artdaq_core::artdaq-core_Utilities
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        lardata_Utilities
                        larevt_Filters
                        lardataobj_RawData
                        lardata_ArtDataHelper
                        ${ROOT_BASIC_LIB_LIST}
                        LIBWDA
                        art::Framework_Core
                        art::Framework_Principal
                        art::Framework_Services_Registry
                        art_root_io::tfile_support
                        art_root_io::TFileService_service
                        art::Utilities
                        art::Persistency_Common canvas
                        art::Persistency_Provenance canvas
                        art::Utilities canvas
                        messagefacility::MF_MessageLogger
                        
                        fhiclcpp::fhiclcpp
                        cetlib::cetlib cetlib_except
                        ${TBB}
        )

simple_plugin(DaqDecoderICARUSPMT module
  sbnobj_Common_PMT_Data
  icaruscode_Decode_DecoderTools_Dumpers
  icaruscode_Decode_DataProducts # should become sbnobj_Common_Trigger
  icaruscode_Decode_DecoderTools
  sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_Common
  sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays
  lardataobj_RawData
  artdaq_core::artdaq-core_Data
  art_root_io::TFileService_service
  art_root_io::tfile_support
  art::Framework_Services_Registry
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib_except::cetlib_except
  ROOT::Tree
  )

simple_plugin(PMTconfigurationExtraction module
  icaruscode_Decode_DecoderTools
  sbnobj_Common_PMT_Data
  art::Framework_Services_Registry
  fhiclcpp
  cetlib_except::cetlib_except
  )

simple_plugin(DumpPMTconfiguration module
  sbnobj_Common_PMT_Data
  messagefacility::MF_MessageLogger
  fhiclcpp
  cetlib_except::cetlib_except
  )

simple_plugin(DumpArtDAQfragments module
  icaruscode_Decode_DecoderTools_Dumpers
  artdaq_core::artdaq-core_Data
  messagefacility::MF_MessageLogger
  fhiclcpp
  cetlib_except::cetlib_except
  )

install_headers()
install_fhicl()
install_source()

# Add our tools directory
add_subdirectory(DecoderTools)
add_subdirectory(ChannelMapping)
add_subdirectory(DataProducts)
add_subdirectory(fcl)

