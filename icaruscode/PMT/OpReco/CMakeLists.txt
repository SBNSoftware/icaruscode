
art_make(
	  MODULE_LIBRARIES larcorealg_Geometry
	  		larcore_Geometry_Geometry_service
			lardataobj_RecoBase
			lardataobj_Simulation
			larcoreobj_SummaryData
			larana_OpticalDetector_OpHitFinder
			nusimdata::SimulationBase
			nurandom::RandomUtils_NuRandomService_service
			larsim_PhotonPropagation_PhotonVisibilityService_service
            art::Framework_Core
			art::Framework_Principal
			art::Framework_Services_Registry
			art_root_io::tfile_support ROOT::Core
			art_root_io::TFileService_service
			art::Persistency_Common
			art::Persistency_Provenance
			art::Utilities
			canvas::canvas
			messagefacility::MF_MessageLogger
			
	        fhiclcpp::fhiclcpp
			cetlib::cetlib
			CLHEP::CLHEP
            ROOT::Geom
            ROOT::XMLIO
            ROOT::Gdml
			ROOT::FFTW
			${ROOT_BASIC_LIB_LIST}

        )

install_headers()
install_source()
add_subdirectory(FlashFinder)
add_subdirectory(fcl)
add_subdirectory(driver)
