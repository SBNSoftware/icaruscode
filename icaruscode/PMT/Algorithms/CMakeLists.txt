# Leave this out for the time being. One day the grid will be updated enough...
# cet_set_compiler_flags(EXTRA_CXX_FLAGS -mavx)

art_make(
  LIB_LIBRARIES 
    larcorealg_Geometry
    lardataobj_RawData
#    lardataobj_RecoBase
#    lardata_Utilities
#    larcorealg_Geometry
#    lardata_ArtDataHelper
#    nusimdata_SimulationBase
     canvas
     ${MF_MESSAGELOGGER}
     ${MF_UTILITIES}
     ${FHICLCPP}
     ${CETLIB}
     cetlib_except
     ${ROOT_FFTW}
     ${ROOT_GENVECTOR}
     ${ROOT_BASIC_LIB_LIST}
     ${CLHEP}
#    ${Boost_SYSTEM_LIBRARY}
  )

install_headers()
install_fhicl()
install_source()
