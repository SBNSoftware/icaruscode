cet_enable_asserts()

art_make( LIB_LIBRARIES lardataobj_RecoBase
                        larcorealg_Geometry
                        nusimdata_SimulationBase
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
        )

install_headers()
install_source()
install_fhicl()
