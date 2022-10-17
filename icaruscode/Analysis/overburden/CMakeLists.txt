art_make(
MODULE_LIBRARIES larcorealg::Geometry
	larcore::Geometry_Geometry_service
	larsim::Simulation
        lardataobj::Simulation
	larsim::MCCheater_BackTrackerService_service
	larsim::MCCheater_ParticleInventoryService_service
	lardata::Utilities
	larevt::Filters
	lardataobj::RawData
	lardataobj::RecoBase
	lardataobj::MCBase
	larreco::RecoAlg
	lardata::RecoObjects
	larpandora::LArPandoraInterface
	nusimdata::SimulationBase
	art::Framework_Core
	art::Framework_Principal
	art::Framework_Services_Registry
	art_root_io::tfile_support ROOT::Core
	art_root_io::TFileService_service
	art::Persistency_Common
	art::Persistency_Provenance
	art::Utilities canvas::canvas
	messagefacility::MF_MessageLogger
	fhiclcpp::fhiclcpp
	ROOT::Geom
	ROOT::XMLIO
	ROOT::Gdml
	icaruscode::RecoUtils
	)

#install_headers()
install_fhicl()
install_source()
