cet_test(ParsingToolkit_test
  LIBRARIES
    icaruscode_PMT_Algorithms
  USE_BOOST_UNIT
  )

cet_test(KeyValueParser_test
  LIBRARIES
    icaruscode_Decode_DecoderTools
    icaruscode_PMT_Algorithms
  USE_BOOST_UNIT
  )
