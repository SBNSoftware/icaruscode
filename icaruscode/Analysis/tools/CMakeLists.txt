cet_enable_asserts()
add_subdirectory(MCTruth)

art_make( TOOL_LIBRARIES lardataobj_RecoBase
                         lardataobj_Simulation
                         icaruscode_RawDigitFilter_Algorithms
                         larcorealg_Geometry
                         larcore_Geometry_Geometry_service
                         lardata_RecoObjects
                         lardata_Utilities
                         nusimdata_SimulationBase
                         ${ART_FRAMEWORK_CORE}
                         ${ART_FRAMEWORK_PRINCIPAL}
                         ${ART_FRAMEWORK_SERVICES_REGISTRY}
                         ${ART_FRAMEWORK_SERVICES_BASIC}
                         ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                         ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                         ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                         art_Persistency_Common
                         art_Persistency_Provenance
                         art_Utilities
                         canvas
                         ${MF_MESSAGELOGGER}
                         ${MF_UTILITIES}
                         ${FHICLCPP}
                         ${CETLIB}
                         cetlib_except
                         ${Boost_FILESYSTEM_LIBRARY}
                         ${Boost_SYSTEM_LIBRARY}
                         ${CLHEP}
                         ${ROOT_BASIC_LIB_LIST}
                         ${ROOT_MINUIT}
        )

install_headers()
install_source()
install_fhicl()
