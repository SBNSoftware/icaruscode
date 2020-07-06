cet_make(
  SUBDIRS
    "details"
  LIBRARIES
    icaruscode_PMT_Trigger_Data
    icaruscode_PMT_Algorithms
    lardataalg_DetectorInfo
    lardataalg_MCDumpers
    lardataobj_RawData
    nusimdata_SimulationBase
    MF_MessageLogger
    fhiclcpp
  )

install_headers(SUBDIRS "details")
install_source(SUBDIRS "details")
install_fhicl(SUBDIRS "details")
