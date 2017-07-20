add_subdirectory(tools)


art_make( BASENAME_ONLY
          EXCLUDE SimWireICARUS_module.cc 
          MODULE_LIBRARIES larcore_Geometry
                        larcore_Geometry_Geometry_service
                        larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
                        lardata_Utilities
                        larevt_Filters
                        lardataobj_RawData
                        SignalShapingServiceICARUS_service
			larevt_CalibrationDBI_Providers
                        nutools_RandomUtils_NuRandomService_service
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
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
                        ${ROOT_BASIC_LIB_LIST}
        )

simple_plugin(SimWireICARUS "module"
          BASENAME_ONLY
              larcorealg_Geometry
              larcore_Geometry_Geometry_service
              larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
              lardata_Utilities
              larevt_Filters
              lardataobj_RawData
              lardata_DetectorInfoServices_DetectorClocksServiceStandard_service # BUG note that we link to a specific service implementation
              SignalShapingServiceICARUS_service
              larevt_CalibrationDBI_Providers
              nutools_RandomUtils_NuRandomService_service
              ${ART_FRAMEWORK_CORE}
              ${ART_FRAMEWORK_PRINCIPAL}
              ${ART_FRAMEWORK_SERVICES_REGISTRY}
              ${ART_FRAMEWORK_SERVICES_OPTIONAL}
              ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
              ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
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
              ${ROOT_BASIC_LIB_LIST}
        )

#install_headers()
install_fhicl()
install_source()
