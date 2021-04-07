art_make_library(
  SOURCE
    "FragmentDumper.cxx"
  LIBRARIES
    artdaq-core_Data
    ${CETLIB_EXCEPT}
  )
