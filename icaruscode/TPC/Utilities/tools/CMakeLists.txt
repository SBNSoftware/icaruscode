cet_enable_asserts()


art_make_library()
set(	TOOL_LIBRARIES
		larcorealg::Geometry
		larevt::CalibrationDBI_IOVData
		larevt::CalibrationDBI_Providers
		lardataobj::RecoBase
		larcore::Geometry_Geometry_service
		lardata::Utilities
		lardataalg::DetectorInfo
		FFTW3::FFTW3
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

cet_build_plugin(ElectronicsResponseBesselApprox art::tool LIBRARIES ${TOOL_LIBRARIES})
cet_build_plugin(ElectronicsResponse art::tool LIBRARIES LIBRARIES ${TOOL_LIBRARIES})
cet_build_plugin(FieldResponse art::tool LIBRARIES ${TOOL_LIBRARIES})
cet_build_plugin(Filter art::tool LIBRARIES ${TOOL_LIBRARIES})
cet_build_plugin(Response art::tool LIBRARIES ${TOOL_LIBRARIES})


install_headers()
install_source()
install_fhicl()
