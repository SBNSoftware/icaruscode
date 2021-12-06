# Leave this out for the time being. One day the grid will be updated enough...
# cet_set_compiler_flags(EXTRA_CXX_FLAGS -mavx)

art_make(
  LIB_LIBRARIES 
    larcorealg_Geometry
    lardataobj_RawData
    icarusalg_Utilities
#    lardataobj_RecoBase
#    lardata_Utilities
#    larcorealg_Geometry
#    lardata_ArtDataHelper
#    nusimdata::SimulationBase
     canvas::canvas
     messagefacility::MF_MessageLogger
     
     fhiclcpp::fhiclcpp
     cetlib::cetlib
     cetlib_except::cetlib_except
     ROOT::FFTW
     ${ROOT_GENVECTOR}
     ${ROOT_BASIC_LIB_LIST}
     CLHEP::CLHEP
#    
  )

install_headers()
install_fhicl()
install_source()
