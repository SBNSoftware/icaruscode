add_subdirectory(Data)
add_subdirectory(Algorithms)
add_subdirectory(Utilities)

art_make(
  EXCLUDE "TriggerEfficiencyPlots_module.cc"
  MODULE_LIBRARIES
    icaruscode_PMT_Trigger_Algorithms
    icaruscode_PMT_Trigger_Data
    lardataalg_DetectorInfo
    lardataobj_RawData
    ${ART_FRAMEWORK_SERVICES_REGISTRY}
    ${MF_MESSAGELOGGER}
    ${FHICLCPP}
    ${ROOT_CORE}
    ${ROOT_PHYSICS}
  TOOL_LIBRARIES
    icaruscode_PMT_Trigger_Algorithms
    icaruscode_PMT_Trigger_Data
    lardataalg_DetectorInfo
    lardataobj_RawData
    ${MF_MESSAGELOGGER}
    ${FHICLCPP}
  )

simple_plugin(TriggerEfficiencyPlots module
  icaruscode_PMT_Trigger_Utilities
  icaruscode_PMT_Trigger_Data
  larcore_Geometry_Geometry_service
  lardataalg_MCDumpers
  larcorealg_Geometry
  lardataobj_RawData
  nusimdata_SimulationBase
  ${ART_ROOT_IO_TFILESERVICE_SERVICE}
  ${ART_ROOT_IO_TFILE_SUPPORT}
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ${ART_FRAMEWORK_PRINCIPAL}
  ${ART_FRAMEWORK_CORE}
  ${ART_FRAMEWORK_PRINCIPAL}
  ${CANVAS}
  ${MF_MESSAGELOGGER}
  ${FHICLCPP}
  ${CETLIB_EXCEPT}
  ${ROOT_GENVECTOR}
  ${ROOT_HIST}
  ${ROOT_TREE}
  ${ROOT_CORE}
  )

install_headers()
install_source()
install_fhicl()
