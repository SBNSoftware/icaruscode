
art_make_library(
	LIBRARIES
		icarusalg::Geometry
		larcore::Geometry_Geometry_service
		larcorealg::Geometry
		lardataobj::RecoBase
		canvas::canvas
		cetlib::cetlib
		cetlib_except::cetlib_except
		art::Framework_Core
		art::Framework_Principal
		art::Framework_Services_Registry
		art_root_io::tfile_support ROOT::Core
		art_root_io::TFileService_service
		art::Persistency_Common
		art::Persistency_Provenance
		art::Utilities
		messagefacility::MF_MessageLogger
		fhiclcpp::fhiclcpp
		ROOT::Geom
		ROOT::XMLIO
		ROOT::Gdml
)
set(	MODULE_LIBRARIES
		icaruscode::PMT_OpReco_FlashFinder
		icarusalg::Geometry
		larcore::Geometry_Geometry_service
		larcorealg::Geometry
		lardataobj::RecoBase
		lardata::Utilities
		larcore::Geometry_AuxDetGeometry_service
		canvas::canvas
		cetlib::cetlib
		cetlib_except::cetlib_except
		art::Framework_Core
		art::Framework_Principal
		art::Framework_Services_Registry
		art_root_io::tfile_support
		art_root_io::TFileService_service
		art::Persistency_Common
		art::Persistency_Provenance
		art::Utilities
		messagefacility::MF_MessageLogger
		fhiclcpp::fhiclcpp
		ROOT::Geom
		ROOT::XMLIO
		ROOT::Gdml
		ROOT::Core
)
cet_build_plugin(ICARUSFlashFinder art::module LIBRARIES ${MODULE_LIBRARIES})


install_headers()
install_source()

