cet_enable_asserts()

add_subdirectory(DataProducts)

art_make(
	LIB_LIBRARIES
		art_Framework_Core
		art_Framework_Principal
		${ART_FRAMEWORK_SERVICES_BASIC}
		${ART_FRAMEWORK_SERVICES_SYSTEM_FILECATALOGMETADATA_SERVICE}
        art_Utilities
        canvas
        ${FHICLCPP}
        ${MF_MESSAGELOGGER}
        cetlib_except
        cetlib 
        LIBWDA
        icaruscode_Timing_DataProducts
	SERVICE_LIBRARIES
		icaruscode_Timing
		icaruscode_Timing_DataProducts
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
install_headers()
install_fhicl()
install_source()