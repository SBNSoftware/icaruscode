
art_make_library(
	LIBRARIES
		lardataobj::RawData
		lardataobj::RecoBase
		lardata::Utilities
		fhiclcpp::fhiclcpp
		cetlib::cetlib
		cetlib_except::cetlib_except
		larcorealg::Geometry
		larcore::Geometry_Geometry_service
		lardata::Utilities
		larevt::Filters
		lardataobj::RawData
		larevt::CalibrationDBI_IOVData
		larevt::CalibrationDBI_Providers
		lardataobj::RecoBase
		icaruscode::TPC_Utilities_SignalShapingICARUSService_service
		icarus_signal_processing::icarus_signal_processing
		art::Framework_Core
		art::Framework_Principal
		art::Framework_Services_Registry
		art_root_io::tfile_support ROOT::Core
		art_root_io::TFileService_service
		art::Persistency_Common canvas::canvas
		art::Persistency_Provenance
		art::Utilities
		messagefacility::MF_MessageLogger
		messagefacility::headers
		fhiclcpp::fhiclcpp
		cetlib::cetlib
		cetlib_except::cetlib_except
		ROOT::Geom
		ROOT::XMLIO
		ROOT::Gdml
		ROOT::FFTW
		FFTW3::FFTW3
		Eigen3::Eigen
)

target_compile_definitions(icaruscode_TPC_SignalProcessing_RawDigitFilter_Algorithms PUBLIC EIGEN_FFTW_DEFAULT)

set(	TOOL_LIBRARIES
		icaruscode_TPC_SignalProcessing_RawDigitFilter_Algorithms
		icaruscode_TPC_Utilities_SignalShapingICARUSService_service
		lardataobj::RawData
		lardataobj::RecoBase
		lardata::Utilities
		fhiclcpp::fhiclcpp
		cetlib::cetlib cetlib_except::cetlib_except
		larcorealg::Geometry
		larcore::Geometry_Geometry_service
		larevt::CalibrationDBI_IOVData
		larevt::CalibrationDBI_Providers
		art::Framework_Core
		art::Framework_Principal
		art::Framework_Services_Registry
		art_root_io::tfile_support ROOT::Core
		art_root_io::TFileService_service
		art::Persistency_Common canvas::canvas
		art::Persistency_Provenance
		art::Utilities
		messagefacility::MF_MessageLogger
		messagefacility::headers
		fhiclcpp::fhiclcpp
		cetlib::cetlib cetlib_except::cetlib_except
		ROOT::Geom
		ROOT::FFTW
		FFTW3::FFTW3

)
cet_build_plugin(MorphologicalFilter art::tool LIBRARIES ${TOOL_LIBRARIES})
cet_build_plugin(RawDigitFilterAlg art::tool LIBRARIES ${TOOL_LIBRARIES})


install_headers()
install_fhicl()
install_source()
