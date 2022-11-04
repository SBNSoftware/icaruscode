cet_enable_asserts()
add_subdirectory(MCTruth)

set(         TOOL_LIBRARIES
		lardataobj::RecoBase
		lardataobj::Simulation
		lardataalg::DetectorInfo
		icaruscode::TPC_SignalProcessing_RawDigitFilter_Algorithms
		larcorealg::Geometry
		larcore::Geometry_Geometry_service
		lardata::RecoObjects
		lardata::Utilities
		nusimdata::SimulationBase
		larsim::Simulation
		art::Framework_Core
		art::Framework_Principal
		art::Framework_Services_Registry
		art_root_io::tfile_support ROOT::Core
		art_root_io::TFileService_service
		art::Framework_Services_Optional_RandomNumberGenerator_service
		art::Persistency_Common
		art::Persistency_Provenance
		art::Utilities
		canvas::canvas
		messagefacility::MF_MessageLogger
		messagefacility::headers
		fhiclcpp::fhiclcpp
		cetlib::cetlib
		cetlib_except::cetlib_except
		Boost::filesystem
		Boost::system
		CLHEP::Random
		ROOT::Minuit ROOT::Tree
		FFTW3::FFTW3
		)

cet_build_plugin(BasicHitAnalysis art::tool LIBRARIES ${TOOL_LIBRARIES})
cet_build_plugin(BasicRawDigitAnalysis art::tool LIBRARIES ${TOOL_LIBRARIES})
cet_build_plugin(BasicTrackAnalysis art::tool LIBRARIES ${TOOL_LIBRARIES})
cet_build_plugin(BasicWireAnalysis art::tool LIBRARIES ${TOOL_LIBRARIES})
cet_build_plugin(HitEfficiencyAnalysis art::tool LIBRARIES ${TOOL_LIBRARIES})
cet_build_plugin(HitFinderAnalysis art::tool LIBRARIES ${TOOL_LIBRARIES} SOURCE HitFinderAnalysis_tool.cc)
set_property(SOURCE HitFinderAnalysis_tool.cc APPEND PROPERTY COMPILE_DEFINITIONS EIGEN_FFTW_DEFAULT)
cet_build_plugin(SpacePointAnalysisMC art::tool LIBRARIES ${TOOL_LIBRARIES} SOURCE SpacePointAnalysisMC_tool.cc)
set_property(SOURCE SpacePointAnalysisMC_tool.cc APPEND PROPERTY COMPILE_DEFINITIONS EIGEN_FFTW_DEFAULT)
cet_build_plugin(SpacePointAnalysis art::tool LIBRARIES ${TOOL_LIBRARIES} SOURCE SpacePointAnalysis_tool.cc)
set_property(SOURCE SpacePointAnalysis_tool.cc APPEND PROPERTY COMPILE_DEFINITIONS EIGEN_FFTW_DEFAULT)
cet_build_plugin(TrackHitEfficiencyAnalysis art::tool LIBRARIES ${TOOL_LIBRARIES} SOURCE TrackHitEfficiencyAnalysis_tool.cc)
set_property(SOURCE TrackHitEfficiencyAnalysis_tool.cc APPEND PROPERTY COMPILE_DEFINITIONS EIGEN_FFTW_DEFAULT)

install_headers()
install_source()
install_fhicl()
