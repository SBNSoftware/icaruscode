cet_test(TriggerInfo_t_test
  LIBRARIES
    sbnobj::ICARUS_PMT_Trigger_Data
    larcorealg::headers
    lardataalg::headers
  USE_BOOST_UNIT
  )

cet_test(AdderSignalSimulation_test
  LIBRARIES
    icaruscode::PMT_Trigger_Algorithms
    # icarusalg::Utilities
    larcorealg::headers
    lardataalg::headers
    lardataobj::RawData
    larcoreobj::SimpleTypesAndConstants
  USE_BOOST_UNIT
  )

