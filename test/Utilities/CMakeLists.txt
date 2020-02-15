cet_test(SampledFunction_test USE_BOOST_UNIT)

cet_test(NonRandomCounter_test
  LIBRARIES
    icaruscode_Utilities
    ${CLHEP}
  USE_BOOST_UNIT
  )
