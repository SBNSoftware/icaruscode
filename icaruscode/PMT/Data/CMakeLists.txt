art_make_library(
  LIBRARIES
    lardataobj::headers
  NO_DICTIONARY
  )

art_dictionary(DICTIONARY_LIBRARIES icaruscode_PMT_Data )

install_headers()
install_source()
