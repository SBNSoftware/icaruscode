if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
cet_add_compiler_flags(CXX -Wno-maybe-uninitialized)
endif()

art_make(
  LIB_LIBRARIES
    lardataobj_RawData
    art_root_io::tfile_support
    art::Framework_Principal
    art::Framework_Core
    canvas
    messagefacility::MF_MessageLogger
    cetlib_except
    ROOT::Tree
    ROOT::RIO
    ROOT::Core
  )

install_headers()
install_source()
install_fhicl()
