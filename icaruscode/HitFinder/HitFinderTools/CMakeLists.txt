cet_enable_asserts()

set( hitfinder_tool_lib_list
			larcorealg_Geometry
			lardataobj_RecoBase
			larcore_Geometry_Geometry_service
			lardata_Utilities
			nutools_RandomUtils_NuRandomService_service
  			${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			art_Persistency_Common
			art_Persistency_Provenance
			${ART_FRAMEWORK_SERVICES_BASIC}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
            ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			${ART_FRAMEWORK_SERVICES_SYSTEM_FILECATALOGMETADATA_SERVICE}
			art_Utilities
			canvas
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
			${FHICLCPP}
			${CETLIB}
			${Boost_FILESYSTEM_LIBRARY}
			${Boost_SYSTEM_LIBRARY}
            ${CLHEP}
			${ROOT_BASIC_LIB_LIST}
    )

art_make(TOOL_LIBRARIES ${hitfinder_tool_lib_list}
        )

install_headers()
install_source()
install_fhicl()
