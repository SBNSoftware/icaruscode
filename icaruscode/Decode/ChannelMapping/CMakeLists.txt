cet_enable_asserts()

art_make(
          EXCLUDE "PMTChannelMapDumper.cxx"
          LIB_LIBRARIES
                        art_Utilities
                        canvas
                        ${FHICLCPP}
                        ${MF_MESSAGELOGGER}
                        cetlib_except
                        cetlib
          SERVICE_LIBRARIES
                        icaruscode_Decode_ChannelMapping
                        lardata_Utilities
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
                        ${TBB}
          TOOL_LIBRARIES
                        lardata_Utilities
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
                        ${SQLITE3}
                        ${TBB}
        )

art_make_exec(NAME "PMTChannelMapDumper"
  LIBRARIES
    icaruscode_Decode_ChannelMapping
    ${MF_MESSAGELOGGER}
    ${FHICLCPP}
    cetlib
    cetlib_except
    Boost::filesystem
  )

install_headers()
install_fhicl()
install_source()

