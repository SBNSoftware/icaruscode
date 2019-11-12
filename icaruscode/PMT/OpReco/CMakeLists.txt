
art_make(
	  MODULE_LIBRARIES larcorealg_Geometry
	  		larcore_Geometry_Geometry_service
			lardataobj_RecoBase
			lardataobj_Simulation
			larcoreobj_SummaryData
			larana_OpticalDetector_OpHitFinder
			nusimdata_SimulationBase
			nurandom_RandomUtils_NuRandomService_service
			larsim_PhotonPropagation_PhotonVisibilityService_service
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
			${ART_ROOT_IO_TFILESERVICE_SERVICE}
			art_Persistency_Common
			art_Persistency_Provenance
			art_Utilities
			canvas
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
	                ${FHICLCPP}
			${CETLIB}
			${CLHEP}
                        ${ROOT_GEOM}
                        ${ROOT_XMLIO}
                        ${ROOT_GDML}
			${ROOT_FFTW}
			${ROOT_BASIC_LIB_LIST}
        )

install_headers()
install_source()
add_subdirectory(FlashFinder)
add_subdirectory(fcl)
add_subdirectory(driver)
