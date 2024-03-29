add_subdirectory(tools)
add_subdirectory(overburden)
add_subdirectory(trigger)

set(         MODULE_LIBRARIES
                larcorealg::Geometry
		larcore::Geometry_Geometry_service
		larsim::Simulation
		nug4::ParticleNavigation
		lardataobj::Simulation
		lardata::Utilities
		larevt::Filters
		larsim::MCCheater_BackTrackerService_service
		larsim::MCCheater_ParticleInventoryService_service
		lardataobj::RawData
		larevt::CalibrationDBI_Providers
		nurandom::RandomUtils_NuRandomService_service
		sbnobj::Common_Analysis
		art::Framework_Core
		art::Framework_Principal
		art::Framework_Services_Registry
		art_root_io::tfile_support
		art::Framework_Services_Optional_RandomNumberGenerator_service
		art_root_io::TFileService_service
		art::Persistency_Common
		art::Persistency_Provenance
		art::Utilities
		canvas::canvas
		messagefacility::MF_MessageLogger
		messagefacility::headers
		fhiclcpp::fhiclcpp
		cetlib::cetlib
		cetlib_except::cetlib_except
		CLHEP::Random
		ROOT::Geom
		ROOT::XMLIO
		ROOT::Gdml
		ROOT::Core
		ROOT::Tree
		icaruscode::RecoUtils
		Eigen3::Eigen
	)

cet_build_plugin(AnalysisTree art::module LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(HitEfficiencyAna art::module LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(ICARUSPurityDQM art::module LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(RawDigitAna art::module LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(TPCPurityInfoAna art::module LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(TPCPurityMonitor art::module LIBRARIES ${MODULE_LIBRARIES})
set_property(SOURCE TPCPurityMonitor_module.cc APPEND PROPERTY COMPILE_DEFINITIONS EIGEN_FFTW_DEFAULT)
cet_build_plugin(ThroughgoingmuonAnalyzer art::module LIBRARIES ${MODULE_LIBRARIES})
set_property(SOURCE ThroughgoingmuonAnalyzer_module.cc APPEND PROPERTY COMPILE_DEFINITIONS EIGEN_FFTW_DEFAULT)
cet_build_plugin(TrackHitAna art::module LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(WireAna art::module LIBRARIES ${MODULE_LIBRARIES})

#install_headers()
install_fhicl()
install_source()
