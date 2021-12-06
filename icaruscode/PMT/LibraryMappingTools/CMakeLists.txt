art_make_library(
  SOURCE
    "ICARUSPhotonMappingTransformations.cxx"
  LIBRARIES
    larsim_PhotonPropagation_LibraryMappingTools
    larcore_Geometry_Geometry_service
    larcorealg_Geometry
    art::Framework_Services_Registry
    art::Utilities
    canvas::canvas
    cetlib_except::cetlib_except
    messagefacility::MF_MessageLogger
    fhiclcpp::fhiclcpp
    ${ROOT_GENVECTOR}
    ${ROOT_BASIC_LIB_LIST}
    
  )

simple_plugin(ICARUSPhotonMappingTransformations "tool"
  icaruscode_PMT_LibraryMappingTools
  larsim_PhotonPropagation_LibraryMappingTools
  )

install_headers()
install_fhicl()
install_source()
