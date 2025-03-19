cet_build_plugin(ICARUSOpHitFinder art::module
  LIBRARIES
  icaruscode_PMT_OpReco_Algorithms
  larana::OpticalDetector_OpHitFinder
  larcore::Geometry_Geometry_service
  larcorealg::Geometry
  lardataobj::RecoBase
  lardataobj::Simulation
  lardataobj::RawData
  larreco::Calorimetry
  art::Framework_Services_Registry
  messagefacility::MF_MessageLogger
  ROOT::Core
  lardata::Utilities
  )


set(	MODULE_LIBRARIES
            icarusalg::Utilities
			sbnobj::Common_Trigger
    	    icaruscode::Decode_DataProducts
            larcorealg::Geometry
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
            ROOT::Tree
)
cet_build_plugin(FakeFlash art::module LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(FakePhotoS art::module LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(FullOpHitFinder art::module LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(ICARUSFlashAssAna art::module LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(ICARUSMCOpFlash art::module LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(ICARUSMCOpHit art::module LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(ICARUSOpFlashAna art::module LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(ICARUSOpHitAna art::module LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(ICARUSOpHitTuple art::module LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(ICARUSParticleAna art::module LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(TPCPMTBarycenterMatchProducer art::module LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(ICARUSBeamStructureAna art::module LIBRARIES ${MODULE_LIBRARIES})


install_headers()
install_source()
add_subdirectory(FlashFinder)
add_subdirectory(Algorithms)
add_subdirectory(fcl)
add_subdirectory(driver)
