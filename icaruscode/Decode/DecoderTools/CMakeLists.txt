cet_enable_asserts()

##include_directories($ENV{LIBWDA_FQ_DIR}/include)

add_subdirectory("Dumpers")

art_make( SUBDIRS details
          EXCLUDE
                        triggerPacketParser.cxx
          LIB_LIBRARIES
                        icaruscode_Utilities
                        sbnobj_Common_PMT_Data
                        icaruscode_Decode_DataProducts
                        canvas
                        ${MF_MESSAGELOGGER}
                        ${FHICLCPP}
                        cetlib_except
                        cetlib
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
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_ROOT_IO_TFILE_SUPPORT}
                        ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                        ${ART_UTILITIES}
                        art_Persistency_Common canvas
                        art_Persistency_Provenance canvas
                        art_Utilities canvas
                        ${MF_MESSAGELOGGER}
                        ${MF_UTILITIES}
                        ${FHICLCPP}
                        cetlib cetlib_except
        )

cet_make_exec(NAME triggerPacketParser
  LIBRARIES PRIVATE
    icaruscode_Decode_DecoderTools
    Boost::program_options
)

install_headers()
install_fhicl()
install_source()
