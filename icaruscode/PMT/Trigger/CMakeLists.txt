add_subdirectory(Data)
add_subdirectory(Algorithms)

art_make(
  TOOL_LIBRARIES
    icaruscode_PMT_Trigger_Algorithms
    icaruscode_PMT_Trigger_Data
    lardataalg_DetectorInfo
    lardataobj_RawData
    MF_MessageLogger
    fhiclcpp
  )

install_headers()
install_source()
install_fhicl()
