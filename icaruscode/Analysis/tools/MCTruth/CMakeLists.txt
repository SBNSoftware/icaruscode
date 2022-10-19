cet_enable_asserts()

art_make(MODULE_LIBRARIES lardataobj::RecoBase
		larcorealg::Geometry
		larcore::Geometry_Geometry_service
		lardata::RecoObjects
		lardata::Utilities
		larsim::MCCheater_BackTrackerService_service
		larsim::MCCheater_ParticleInventoryService_service
		icarusalg::gallery_MCTruthBase
		nug4::ParticleNavigation
		art::Framework_Core
		art::Framework_Principal
		art::Framework_Services_Registry
		art_root_io::tfile_support ROOT::Core
		art_root_io::TFileService_service
		art::Framework_Services_Optional_RandomNumberGenerator_service
		art::Persistency_Common
		art::Persistency_Provenance
		art::Utilities
		canvas::canvas
		messagefacility::MF_MessageLogger
		messagefacility::headers
		fhiclcpp::fhiclcpp
		cetlib::cetlib
		cetlib_except::cetlib_except
		Boost::filesystem
		Boost::system
		CLHEP::Random

	TOOL_LIBRARIES lardataobj::RecoBase
		larcorealg::Geometry
		larcore::Geometry_Geometry_service
		lardata::RecoObjects
		lardata::Utilities
		larsim::MCCheater_BackTrackerService_service
		larsim::MCCheater_ParticleInventoryService_service
		icarusalg::gallery_MCTruthBase
		nug4::ParticleNavigation
		art::Framework_Core
		art::Framework_Principal
		art::Framework_Services_Registry
		art_root_io::tfile_support ROOT::Core
		art_root_io::TFileService_service
		art::Framework_Services_Optional_RandomNumberGenerator_service
		art::Persistency_Common
		art::Persistency_Provenance
		art::Utilities
		canvas::canvas
		messagefacility::MF_MessageLogger
		messagefacility::headers
		fhiclcpp::fhiclcpp
		cetlib::cetlib
		cetlib_except::cetlib_except
		Boost::filesystem
		Boost::system
		CLHEP::Random
	)

install_headers()
install_source()
install_fhicl()
