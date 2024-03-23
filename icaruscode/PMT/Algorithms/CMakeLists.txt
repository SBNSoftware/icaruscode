#Leave this out for the time being. One day the grid will be updated enough ...
#cet_set_compiler_flags(EXTRA_CXX_FLAGS -mavx)

file(GLOB lib_srcs *.cxx)
art_make_library(
        SOURCE ${lib_srcs}
	LIBRARIES
		icaruscode::Decode_DecoderTools
		icarusalg::Utilities
		larcorealg::Geometry
		lardataobj::RawData
		larcoreobj::headers
		Microsoft.GSL::GSL
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
