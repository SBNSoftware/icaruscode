art_make_library(
  SOURCE
    "FragmentDumper.cxx"
  LIBRARIES
    artdaq_core::artdaq-core_Data
    ${CETLIB_EXCEPT}
  )
