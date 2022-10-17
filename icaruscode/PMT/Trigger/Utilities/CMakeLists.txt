if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	cet_add_compiler_flags(CXX -Wno-maybe-uninitialized)
endif()

art_make(
	LIB_LIBRARIES
                larcorealg::CoreUtils
		art_root_io::tfile_support
		messagefacility::MF_MessageLogger
		cetlib_except::cetlib_except
		ROOT::RIO
		ROOT::Core
	)

install_headers()
install_source()
install_fhicl()
