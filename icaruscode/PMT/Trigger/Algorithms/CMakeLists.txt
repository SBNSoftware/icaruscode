cet_make(
  LIBRARIES
    icaruscode_PMT_Trigger_Data
    icaruscode_PMT_Algorithms
    lardataalg_DetectorInfo
    lardataobj_RawData
    MF_MessageLogger
    fhiclcpp
  )

install_headers()
install_source()
install_fhicl()
