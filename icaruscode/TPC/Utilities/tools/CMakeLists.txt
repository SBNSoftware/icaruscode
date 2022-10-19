cet_enable_asserts()
include_directories(${FFTW3_INCLUDE_DIRS})
link_directories(${FFTW3_LIBRARY_DIRS})


add_definitions(-DEIGEN_FFTW_DEFAULT)

art_make(
	TOOL_LIBRARIES
		larcorealg::Geometry
		larevt::CalibrationDBI_IOVData
		larevt::CalibrationDBI_Providers
		lardataobj::RecoBase
		larcore::Geometry_Geometry_service
		lardata::Utilities
		lardataalg::DetectorInfo
		${ICARUS_FFTW_LIBRARIES}
		art::Framework_Core
		art::Framework_Principal
		art::Framework_Services_Registry
		art::Persistency_Common
		art::Persistency_Provenance
		icarus_signal_processing::icarus_signal_processing
		art_root_io::tfile_support ROOT::Core
		art::Framework_Services_Optional_RandomNumberGenerator_service
		art_root_io::TFileService_service
		art::Framework_Services_System_FileCatalogMetadata_service
		art::Utilities
		canvas::canvas
		messagefacility::MF_MessageLogger
		messagefacility::headers
		fhiclcpp::fhiclcpp
		cetlib::cetlib
		Boost::filesystem
		Boost::system
		CLHEP::Random
	)

install_headers()
install_source()
install_fhicl()
