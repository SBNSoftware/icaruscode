cet_make( 
  LIBRARIES
    lardataobj_RawData
    cetlib_except
  NO_DICTIONARY
  )

art_dictionary(DICTIONARY_LIBRARIES icaruscode_PMT_Trigger_Data)

install_headers()
install_source()
