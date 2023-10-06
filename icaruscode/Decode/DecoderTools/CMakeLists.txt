cet_enable_asserts()

add_subdirectory("Dumpers")

art_make_library(SUBDIRS details
        EXCLUDE
                triggerPacketParser.cxx
        LIBRARIES
                lardataalg::DetectorInfo
                icaruscode::Utilities
                icaruscode::Decode_DataProducts
                sbnobj::Common_PMT_Data
                messagefacility::MF_MessageLogger
                fhiclcpp::fhiclcpp
                cetlib_except::cetlib_except
                cetlib::cetlib
                ROOT::RIO
                ROOT::Core
                art::Framework_Services_Registry
                FFTW3::FFTW3
)

set(        TOOL_LIBRARIES
                sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_ICARUS
                artdaq_core::artdaq-core_Utilities
                icarus_signal_processing::icarus_signal_processing
                icarus_signal_processing::Detection
                icarus_signal_processing::Filters
                icaruscode::TPC_Utilities_SignalShapingICARUSService_service
                icaruscode::Decode_DecoderTools
                icaruscode::Decode_DecoderTools_Dumpers
                icaruscode::Utilities
                larcorealg::Geometry
                larcore::Geometry_Geometry_service
                lardata::Utilities
                lardataalg::DetectorInfo
                larevt::Filters
                lardataobj::RawData
                lardata::ArtDataHelper
                ROOT::Core
                Eigen3::Eigen
                wda::wda
                art::Framework_Core
                art::Framework_Principal
                art::Framework_Services_Registry
                art_root_io::tfile_support
                art_root_io::TFileService_service
                art::Utilities
                art::Persistency_Common canvas::canvas
                art::Persistency_Provenance
                art::Utilities
                messagefacility::MF_MessageLogger
                messagefacility::headers
                fhiclcpp::fhiclcpp
                cetlib::cetlib cetlib_except::cetlib_except
                Boost::program_options
        )

cet_build_plugin(FakeParticle art::tool LIBRARIES ${TOOL_LIBRARIES})
cet_build_plugin(TPCDecoderFilter1D art::tool LIBRARIES ${TOOL_LIBRARIES})
cet_build_plugin(TPCDecoderFilter2D art::tool LIBRARIES ${TOOL_LIBRARIES})
cet_build_plugin(TPCDecoder art::tool LIBRARIES ${TOOL_LIBRARIES})
cet_build_plugin(TPCNoiseFilter1D art::tool LIBRARIES ${TOOL_LIBRARIES})
cet_build_plugin(TPCNoiseFilterCanny art::tool LIBRARIES ${TOOL_LIBRARIES})
cet_build_plugin(TriggerDecoder art::tool LIBRARIES ${TOOL_LIBRARIES})
cet_build_plugin(TriggerDecoderV2 art::tool LIBRARIES ${TOOL_LIBRARIES})
cet_build_plugin(TriggerDecoderV3 art::tool LIBRARIES ${TOOL_LIBRARIES})
cet_build_plugin(TriggerDecoderDispatcher art::tool LIBRARIES ${TOOL_LIBRARIES})



cet_make_exec(NAME triggerPacketParser
  LIBRARIES
    icaruscode_Decode_DecoderTools
    Boost::program_options
)

install_headers()
install_fhicl()
install_source()
