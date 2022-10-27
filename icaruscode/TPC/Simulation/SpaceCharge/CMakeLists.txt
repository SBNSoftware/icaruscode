art_make_library(
	LIBRARIES
		lardataobj::RecoBase
		lardata::Utilities
		lardataalg::DetectorInfo
		larevt::SpaceCharge
		larcorealg::Geometry
		larsim::MCCheater_BackTracker
		larcore::Geometry_Geometry_service
		messagefacility::MF_MessageLogger
		cetlib::cetlib cetlib_except::cetlib_except
		art::Framework_Core
		art::Framework_Principal
		art::Framework_Services_Registry
		art_root_io::tfile_support ROOT::Core
		art::Framework_Services_Optional_RandomNumberGenerator_service
		art_root_io::TFileService_service
		art::Persistency_Common canvas::canvas
		art::Persistency_Provenance
		art::Utilities
		fhiclcpp::fhiclcpp
		Boost::system
	)
install_headers()
install_fhicl()
install_source()
