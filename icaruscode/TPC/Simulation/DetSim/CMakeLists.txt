add_subdirectory(tools)

add_definitions(-DEIGEN_FFTW_DEFAULT)

art_make( 
          EXCLUDE SimWireICARUS_module.cc 
          MODULE_LIBRARIES larcorealg_Geometry
                           larcore_Geometry_Geometry_service
                           larsim_Simulation 
                           nug4::ParticleNavigation lardataobj_Simulation
                           lardata_Utilities
                           larevt_Filters
                           lardataobj_RawData
                           icaruscode_TPC_Utilities_SignalShapingICARUSService_service
		           larevt_CalibrationDBI_Providers
                           nurandom::RandomUtils_NuRandomService_service
                           art::Framework_Core
                           art::Framework_Principal
                           art::Framework_Services_Registry
                           art_root_io::tfile_support ROOT::Core
                           art::Framework_Services_Optional_RandomNumberGenerator_service
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
                           ${ROOT_BASIC_LIB_LIST}
        )

simple_plugin(SimWireICARUS "module"
                           larcorealg_Geometry
                           larcore_Geometry_Geometry_service
                           larsim_Simulation 
                           nug4::ParticleNavigation lardataobj_Simulation
                           lardata_Utilities
                           lardataalg_DetectorInfo
                           larevt_Filters
                           lardataobj_RawData
                           icaruscode_TPC_Utilities_SignalShapingICARUSService_service
                           larevt_CalibrationDBI_Providers
                           nurandom::RandomUtils_NuRandomService_service
                           art::Framework_Core
                           art::Framework_Principal
                           art::Framework_Services_Registry
                           art_root_io::tfile_support ROOT::Core
                           art::Framework_Services_Optional_RandomNumberGenerator_service
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
                           ${ROOT_BASIC_LIB_LIST}
        )

#install_headers()
install_fhicl()
install_source()
