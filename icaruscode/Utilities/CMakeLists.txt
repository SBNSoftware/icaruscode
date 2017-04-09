cet_enable_asserts()

set( icarus_util_lib_list lardata_Utilities_LArFFT_service
			larcore_Geometry
			larevt_CalibrationDBI_IOVData
			larevt_CalibrationDBI_Providers
			larreco_HitFinder
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


simple_plugin( SignalShapingServiceICARUS  "service"
               ${icarus_util_lib_list}
	       BASENAME_ONLY
        )

install_headers()
install_source()
install_fhicl()
