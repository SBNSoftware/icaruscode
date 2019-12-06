add_subdirectory(Data)
add_subdirectory(Algorithms)

art_make(
  MODULE_LIBRARIES
    icaruscode_PMT_Trigger_Algorithms
    icaruscode_PMT_Trigger_Data
    lardataalg_DetectorInfo
    lardataobj_RawData
    ${ART_FRAMEWORK_SERVICES_REGISTRY}
    ${MF_MESSAGELOGGER}
    ${FHICLCPP}
  TOOL_LIBRARIES
    icaruscode_PMT_Trigger_Algorithms
    icaruscode_PMT_Trigger_Data
    lardataalg_DetectorInfo
    lardataobj_RawData
    ${MF_MESSAGELOGGER}
    ${FHICLCPP}
  )

install_headers()
install_source()
install_fhicl()
