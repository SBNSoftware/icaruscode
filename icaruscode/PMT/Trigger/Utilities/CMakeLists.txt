if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	cet_add_compiler_flags(CXX -Wno-maybe-uninitialized)
endif()

art_make_library(
	LIBRARIES
		lardataobj::RawData
		larcorealg::CoreUtils
		art_root_io::tfile_support
		art::Framework_Core
		messagefacility::MF_MessageLogger
		cetlib_except::cetlib_except
		ROOT::RIO
		ROOT::Core
	)

install_headers()
install_source()
install_fhicl()
