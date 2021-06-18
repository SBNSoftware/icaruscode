art_make( 
          MODULE_LIBRARIES
                        icarusalg_Utilities
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        larsim_Simulation 
                        nug4_ParticleNavigation
                        lardataobj_Simulation
                        lardata_Utilities
                        larevt_Filters
                        lardataalg_MCDumpers
                        lardataobj_RawData
                        lardataobj_RecoBase
                        lardata_RecoObjects
                        nusimdata_SimulationBase
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ROOT_CORE}
                        art_Persistency_Common 
                        art_Persistency_Provenance 
                        art_Utilities 
                        ${MF_MESSAGELOGGER}
                        ${MF_UTILITIES}
                        ${FHICLCPP}
                        ${CETLIB}
                        ${CETLIB_EXCEPT}
                        ${ROOT_GEOM}
                        ${ROOT_XMLIO}
                        ${ROOT_GDML}
                        ${ROOT_GENVECTOR}
                        ${ROOT_BASIC_LIB_LIST}
        )

# install_headers()
install_fhicl()
install_source()

