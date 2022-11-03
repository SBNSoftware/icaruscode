cet_enable_asserts()
add_subdirectory(MCTruth)

add_definitions(-DEIGEN_FFTW_DEFAULT)

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
cet_build_plugin(HitFinderAnalysis art::tool LIBRARIES ${TOOL_LIBRARIES})
cet_build_plugin(SpacePointAnalysisMC art::tool LIBRARIES ${TOOL_LIBRARIES})
cet_build_plugin(SpacePointAnalysis art::tool LIBRARIES ${TOOL_LIBRARIES})
cet_build_plugin(TrackHitEfficiencyAnalysis art::tool LIBRARIES ${TOOL_LIBRARIES})

install_headers()
install_source()
install_fhicl()
