simple_plugin(SpaceChargeServiceICARUS "service"
	icaruscode::TPC_Simulation_SpaceCharge
	lardataalg::DetectorInfo
	art::Framework_Core
	art::Framework_Principal
	art::Framework_Services_Registry
	art::Persistency_Common
	art::Persistency_Provenance
	art_root_io::tfile_support ROOT::Core
	art_root_io::TFileService_service
	art::Framework_Services_System_FileCatalogMetadata_service
	art::Utilities canvas::canvas
	messagefacility::MF_MessageLogger								
	cetlib::cetlib cetlib_except::cetlib_except
)
install_headers()
install_fhicl()
install_source()
