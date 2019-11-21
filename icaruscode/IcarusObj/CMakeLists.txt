cet_make(NO_DICTIONARY
         LIBRARIES
#           canvas
#           MF_MessageLogger
           cetlib_except
           ${ROOT_BASIC_LIB_LIST}
           ${ROOT_EG}
           ${ROOT_CORE})
art_dictionary(DICTIONARY_LIBRARIES icaruscode_IcarusObj)
install_headers()
install_source()
install_fhicl()
