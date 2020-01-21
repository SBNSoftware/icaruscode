add_subdirectory(tools)


art_make( 
          EXCLUDE SimWireICARUS_module.cc 
          MODULE_LIBRARIES larcorealg_Geometry
                           larcore_Geometry_Geometry_service
                           larsim_Simulation 
                           nug4_ParticleNavigation lardataobj_Simulation
                           lardata_Utilities
                           larevt_Filters
                           lardataobj_RawData
                           icaruscode_Utilities_SignalShapingServiceICARUS_service
			                  larevt_CalibrationDBI_Providers
                           nurandom_RandomUtils_NuRandomService_service
                           ${ART_FRAMEWORK_CORE}
                           ${ART_FRAMEWORK_PRINCIPAL}
                           ${ART_FRAMEWORK_SERVICES_REGISTRY}
                           ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
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
                           ${ROOT_BASIC_LIB_LIST}
        )

simple_plugin(SimWireICARUS "module"
                           larcorealg_Geometry
                           larcore_Geometry_Geometry_service
                           larsim_Simulation 
                           nug4_ParticleNavigation lardataobj_Simulation
                           lardata_Utilities
                           larevt_Filters
                           lardataobj_RawData
                           icaruscode_Utilities_SignalShapingServiceICARUS_service
                           larevt_CalibrationDBI_Providers
                           nurandom_RandomUtils_NuRandomService_service
                           ${ART_FRAMEWORK_CORE}
                           ${ART_FRAMEWORK_PRINCIPAL}
                           ${ART_FRAMEWORK_SERVICES_REGISTRY}
                           ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
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
                           ${ROOT_BASIC_LIB_LIST}
        )

#install_headers()
install_fhicl()
install_source()
