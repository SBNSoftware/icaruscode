cet_enable_asserts()

cet_find_library(LIBWDA NAMES wda PATHS ENV LIBWDA_LIB NO_DEFAULT_PATH)

include_directories($ENV{LIBWDA_FQ_DIR}/include)

add_subdirectory("Dumpers")

art_make(
          EXCLUDE
                       "TriggerDecoder_tool.cc"
          LIB_LIBRARIES
                        icaruscode_Utilities
                        sbnobj_Common_PMT_Data
                        ${MF_MESSAGELOGGER}
                        ${FHICLCPP}
                        cetlib_except
                        cetlib
                        ROOT::RIO
                        ROOT::Core
          TOOL_LIBRARIES
                        sbndaq-artdaq-core_Overlays_ICARUS 
                        artdaq-core_Utilities
                        icarus_signal_processing
                        icarus_signal_processing_Detection
                        icarus_signal_processing_Filters
                        icaruscode_TPC_Utilities_SignalShapingICARUSService_service
                        icaruscode_Decode_DecoderTools
                        icaruscode_Utilities
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        lardata_Utilities
                        lardataalg_DetectorInfo
                        larevt_Filters
                        lardataobj_RawData
                        lardata_ArtDataHelper
                        ${ROOT_BASIC_LIB_LIST}
                        ${LIBWDA}
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

simple_plugin( "TriggerDecoder" tool
  icaruscode_Decode_DecoderTools_Dumpers
  icaruscode_Decode_DecoderTools
  lardataobj_RawData
  sbndaq-artdaq-core_Overlays_ICARUS
  artdaq-core_Utilities
  ${ART_PERSISTENCY_PROVENANCE}
  ${ART_FRAMEWORK_CORE}
  ${ART_FRAMEWORK_PRINCIPAL}
  ${ART_UTILITIES}
  ${CANVAS}
  ${MF_MESSAGELOGGER}
  ${FHICLCPP}
  ${CETLIB}
  )

install_headers()
install_fhicl()
install_source()
