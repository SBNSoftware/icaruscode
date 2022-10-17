if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	cet_add_compiler_flags(CXX -Wno-maybe-uninitialized)
endif()
art_make(LIBRARY_NAME icaruscode_CRTUtils
	LIB_LIBRARIES larcorealg::Geometry
	larcore::Geometry_Geometry_service
	larsim::Simulation
	lardataobj::Simulation
	larsim::MCCheater_BackTrackerService_service
	larsim::MCCheater_ParticleInventoryService_service
	lardata::Utilities
	larevt::Filters
	lardataobj::RawData
	lardataobj::RecoBase
	lardata::RecoObjects
	lardataalg::DetectorInfo
	larpandora::LArPandoraInterface
	nusimdata::SimulationBase
	cetlib::cetlib cetlib_except::cetlib_except
	art::Framework_Core
	art::Framework_Principal
	art::Framework_Services_Registry
	art_root_io::tfile_support ROOT::Core
	art_root_io::TFileService_service
	art::Framework_Services_Optional_RandomNumberGenerator_service
	art::Persistency_Common
	art::Persistency_Provenance
	art::Utilities canvas::canvas
	messagefacility::MF_MessageLogger
	CLHEP::CLHEP
	fhiclcpp::fhiclcpp
	ROOT::Geom
	ROOT::XMLIO
	ROOT::Gdml
	ROOT::GenVector
	Boost::system
	sbnobj::ICARUS_CRT
	sbnobj::Common_CRT
	lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
	nurandom::RandomUtils_NuRandomService_service
	)

install_headers()
install_fhicl()
install_source()

