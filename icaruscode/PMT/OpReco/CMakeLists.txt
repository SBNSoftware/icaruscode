simple_plugin(ICARUSOpHitFinder "module"
  icaruscode_PMT_OpReco_Algorithms
  larana_OpticalDetector_OpHitFinder
  larcore_Geometry_Geometry_service
  larcorealg_Geometry
  lardataobj_RecoBase
  lardataobj_Simulation
  lardataobj_RawData
  art::Framework_Services_Registry
  messagefacility::MF_MessageLogger
  )


art_make(
		EXCLUDE
			"ICARUSOpHitFinder_module.cc"
	MODULE_LIBRARIES larcorealg::Geometry
			larcore::Geometry_Geometry_service
			lardataobj::RecoBase
			lardataobj::Simulation
			larcoreobj::SummaryData
			larana::OpticalDetector_OpHitFinder
                        larreco::Calorimetry
                        lardata::Utilities
			nusimdata::SimulationBase
			nurandom::RandomUtils_NuRandomService_service
			larsim::PhotonPropagation_PhotonVisibilityService_service
			art::Framework_Core
			art::Framework_Principal
			art::Framework_Services_Registry
			art_root_io::tfile_support
			art_root_io::TFileService_service
			art::Persistency_Common
			art::Persistency_Provenance
			art::Utilities
			canvas::canvas
			messagefacility::MF_MessageLogger
			messagefacility::headers
			fhiclcpp::fhiclcpp
			cetlib::cetlib
			CLHEP::Random
			ROOT::Geom
			ROOT::XMLIO
			ROOT::Gdml
			ROOT::FFTW
			ROOT::Core
)

install_headers()
install_source()
add_subdirectory(FlashFinder)
add_subdirectory(Algorithms)
add_subdirectory(fcl)
add_subdirectory(driver)
