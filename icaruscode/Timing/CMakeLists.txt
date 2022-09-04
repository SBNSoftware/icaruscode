cet_enable_asserts()

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
	SERVICE_LIBRARIES
		icaruscode_Timing
		art_Framework_Core
		art_Framework_Principal
		art_Framework_Services_Registry
		art_Persistency_Common canvas
		art_Persistency_Provenance canvas
		${ART_FRAMEWORK_SERVICES_BASIC}
		${ART_FRAMEWORK_SERVICES_SYSTEM_FILECATALOGMETADATA_SERVICE}
		art_Utilities canvas
		messagefacility::MF_MessageLogger
		cetlib
		cetlib_except
		${ROOT_BASIC_LIB_LIST}
		LIBWDA

)
install_headers()
install_fhicl()
install_source()