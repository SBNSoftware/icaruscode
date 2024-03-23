cet_enable_asserts()

art_make_library(
          EXCLUDE
                        TriggerConfigurationExtraction_module.cc
                        PMTconfigurationExtraction_module.cc
                        DumpTriggerConfiguration_module.cc
                        DumpPMTconfiguration_module.cc
                        DumpArtDAQfragments_module.cc
                        DumpTrigger_module.cc
                        DaqDecoderICARUSPMT_module.cc
)

set(          MODULE_LIBRARIES
                        icarus_signal_processing::Detection
                        icarus_signal_processing::Filters
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
                        icaruscode::Decode_DecoderTools_Dumpers
                        Eigen3::Eigen
)

cet_build_plugin(DaqDecoderICARUSTPC art::module LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(DaqDecoderICARUSTPCwROI art::module LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(DaqDecoderICARUSTrigger art::module LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(FilterNoiseICARUS art::module LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(MCDecoderICARUSTPCwROI art::module LIBRARIES ${MODULE_LIBRARIES})

cet_build_plugin(DaqDecoderICARUSPMT art::module LIBRARIES
	icaruscode::Timing
	icaruscode::Decode_DecoderTools_Dumpers
	icaruscode::Decode_DecoderTools
	icaruscode::PMT_Algorithms
	sbnobj::Common_PMT_Data
	sbnobj::Common_Trigger
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

cet_build_plugin(PMTconfigurationExtraction art::module LIBRARIES
	icaruscode::Decode_DecoderTools
	sbnobj::Common_PMT_Data
	sbnobj::Common_Trigger
	art::Framework_Services_Registry
	fhiclcpp::fhiclcpp
	cetlib_except::cetlib_except
	)

cet_build_plugin(TriggerConfigurationExtraction art::module LIBRARIES
  icaruscode_Decode_DecoderTools
  icaruscode_Decode_DataProducts
  art::Framework_Services_Registry
  fhiclcpp::fhiclcpp
  cetlib_except::cetlib_except
  )

cet_build_plugin(DumpPMTconfiguration art::module LIBRARIES
	sbnobj::Common_PMT_Data
        lardataobj::RawData
	messagefacility::MF_MessageLogger
	fhiclcpp::fhiclcpp
	cetlib_except::cetlib_except
	)


cet_build_plugin(DumpTriggerConfiguration art::module LIBRARIES
  icaruscode_Decode_DataProducts
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib_except::cetlib_except
  )

cet_build_plugin(DumpArtDAQfragments art::module LIBRARIES
	icaruscode::Decode_DecoderTools_Dumpers
        larcorealg::headers
	artdaq_core::artdaq-core_Data
	messagefacility::MF_MessageLogger
	fhiclcpp::fhiclcpp
	cetlib_except::cetlib_except
	)

cet_build_plugin(DumpTrigger art::module LIBRARIES
  sbnobj::Common_Trigger
  lardataobj::Simulation
  lardataobj::RawData
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib_except::cetlib_except
  lardataalg::DetectorInfo
  )

cet_build_plugin(CopyDaqICARUSTPC art::module LIBRARIES
  sbnobj::Common_Trigger
  artdaq_core::artdaq-core_Data
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib_except::cetlib_except
  )

install_headers()
install_fhicl()
install_source()

# Add our tools directory
add_subdirectory(DecoderTools)
add_subdirectory(ChannelMapping)
add_subdirectory(DataProducts)

