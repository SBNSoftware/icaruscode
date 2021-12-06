cet_enable_asserts()

##include_directories($ENV{LIBWDA_FQ_DIR}/include)

add_subdirectory("Dumpers")

art_make( SUBDIRS details
          LIB_LIBRARIES
                        icaruscode_Utilities
                        sbnobj_Common_PMT_Data
                        messagefacility::MF_MessageLogger
                        fhiclcpp::fhiclcpp
                        cetlib_except::cetlib_except
                        cetlib::cetlib
                        ROOT::RIO
                        ROOT::Core
          TOOL_LIBRARIES
                        sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_ICARUS 
                        artdaq_core::artdaq-core_Utilities
                        icarus_signal_processing
                        icarus_signal_processing_Detection
                        icarus_signal_processing_Filters
                        icaruscode_TPC_Utilities_SignalShapingICARUSService_service
                        icaruscode_Decode_DecoderTools
                        icaruscode_Decode_DecoderTools_Dumpers
                        icaruscode_Utilities
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        lardata_Utilities
                        lardataalg_DetectorInfo
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
        )

install_headers()
install_fhicl()
install_source()
