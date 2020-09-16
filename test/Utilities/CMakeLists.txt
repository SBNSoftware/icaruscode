cet_test(rounding_test USE_BOOST_UNIT)
cet_test(ChangeMonitor_test USE_BOOST_UNIT)

cet_test(FastAndPoorGauss_test
  LIBRARIES
    ${ROOT_MATHCORE}
    ${ROOT_HIST}
    ${ROOT_RIO}
    ${ROOT_CORE}
  USE_BOOST_UNIT
  )

cet_test(SampledFunction_test USE_BOOST_UNIT)

cet_test(NonRandomCounter_test
  LIBRARIES
    icaruscode_Utilities
    ${CLHEP}
  USE_BOOST_UNIT
  )
