cet_enable_asserts()

art_make(
          EXCLUDE
                        TriggerConfigurationExtraction_module.cc
                        PMTconfigurationExtraction_module.cc
                        DumpTriggerConfiguration_module.cc
                        DumpPMTconfiguration_module.cc
                        DumpArtDAQfragments_module.cc
                        DumpTrigger_module.cc
                        DaqDecoderICARUSPMT_module.cc
          MODULE_LIBRARIES
                        icarus_signal_processing
                        icarus_signal_processing_Detection
                        icarus_signal_processing_Filters
                        icaruscode_TPC_Utilities
                        sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_ICARUS
                        artdaq_core::artdaq-core_Utilities
                        larcorealg::Geometry
                        lardataobj::RawData
                        larcore::Geometry_Geometry_service
                        lardata::Utilities
                        larevt::Filters
                        lardataobj::RawData
                        lardata::ArtDataHelper
)

simple_plugin(DaqDecoderICARUSPMT module
	sbnobj::Common_PMT_Data
	sbnobj::Common_Trigger
	icaruscode_Decode_DecoderTools_Dumpers
	icaruscode_Decode_DecoderTools
	sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_Common
	sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays
	lardataobj::RawData
	larcore::Geometry_Geometry_service
	artdaq_core::artdaq-core_Data
	art_root_io::TFileService_service
	art_root_io::tfile_support
	art::Framework_Services_Registry
	messagefacility::MF_MessageLogger
	fhiclcpp::fhiclcpp
	cetlib_except::cetlib_except
	ROOT::Tree
	lardata::Utilities
	)

simple_plugin(PMTconfigurationExtraction module
	icaruscode::Decode_DecoderTools
	sbnobj::Common_PMT_Data
	sbnobj::Common_Trigger
	art::Framework_Services_Registry
	fhiclcpp::fhiclcpp
	cetlib_except::cetlib_except
	)

simple_plugin(TriggerConfigurationExtraction module
  icaruscode_Decode_DecoderTools
  icaruscode_Decode_DataProducts
  art::Framework_Services_Registry
  fhiclcpp::fhiclcpp
  cetlib_except::cetlib_except
  )

simple_plugin(DumpPMTconfiguration module
	sbnobj::Common_PMT_Data
        lardataobj::RawData
	messagefacility::MF_MessageLogger
	fhiclcpp::fhiclcpp
	cetlib_except::cetlib_except
	)


simple_plugin(DumpTriggerConfiguration module
  icaruscode_Decode_DataProducts
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib_except::cetlib_except
  )

simple_plugin(DumpArtDAQfragments module
	icaruscode::Decode_DecoderTools_Dumpers
        larcorealg::headers
	artdaq_core::artdaq-core_Data
	messagefacility::MF_MessageLogger
	fhiclcpp::fhiclcpp
	cetlib_except::cetlib_except
	)

simple_plugin(DumpTrigger module
  sbnobj::Common_Trigger
  lardataobj::Simulation
  lardataobj::RawData
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib_except::cetlib_except
  lardataalg::DetectorInfo
  )

install_headers()
install_fhicl()
install_source()

# Add our tools directory
add_subdirectory(DecoderTools)
add_subdirectory(ChannelMapping)
add_subdirectory(DataProducts)
#add_subdirectory(fcl)

