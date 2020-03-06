art_make_library(
  SOURCE
    "ICARUSPhotonMappingTransformations.cxx"
  LIBRARIES
    larsim_PhotonPropagation_LibraryMappingTools
    larcore_Geometry_Geometry_service
    larcorealg_Geometry
    ${ART_FRAMEWORK_SERVICES_REGISTRY}
    ${ART_UTILITIES}
    ${CANVAS}
    ${CETLIB_EXCEPT}
    ${MF_MESSAGELOGGER}
    ${FHICLCPP}
    ${ROOT_GENVECTOR}
    ${ROOT_BASIC_LIB_LIST}
    ${Boost_SYSTEM_LIBRARY}
  )

simple_plugin(ICARUSPhotonMappingTransformations "tool"
  icaruscode_PMT_LibraryMappingTools
  larsim_PhotonPropagation_LibraryMappingTools
  )

install_headers()
install_fhicl()
install_source()
