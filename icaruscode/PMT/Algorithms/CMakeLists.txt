#Leave this out for the time being. One day the grid will be updated enough ...
#cet_set_compiler_flags(EXTRA_CXX_FLAGS -mavx)

include_directories(${CPPGSL_INC})

art_make(
	LIB_LIBRARIES
		larcorealg::Geometry
		lardataobj::RawData
		icarusalg::Utilities
                lardataalg::DetectorInfo
		canvas::canvas
		messagefacility::MF_MessageLogger
		messagefacility::headers
		fhiclcpp::fhiclcpp
		cetlib::cetlib
		cetlib_except::cetlib_except
		ROOT::FFTW
		ROOT::GenVector
		CLHEP::Random
	)

install_headers()
install_fhicl()
install_source()
