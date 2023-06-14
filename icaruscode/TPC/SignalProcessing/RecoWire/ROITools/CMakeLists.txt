cet_enable_asserts()

set(
	TOOL_LIBRARIES
		larcorealg::Geometry
		icarus_signal_processing::icarus_signal_processing
		icarus_signal_processing::Detection
		icarus_signal_processing::Filters
		larevt::CalibrationDBI_IOVData
		larevt::CalibrationDBI_Providers
		lardataobj::RecoBase
		larcore::Geometry_Geometry_service
		lardata::Utilities
		lardataalg::DetectorInfo
		icaruscode::TPC_Utilities_SignalShapingICARUSService_service
		nurandom::RandomUtils_NuRandomService_service
		art::Framework_Core
		art::Framework_Principal
		art::Framework_Services_Registry
		art::Persistency_Common
		art::Persistency_Provenance
		ROOT::Tree
		art_root_io::tfile_support ROOT::Core
		art::Framework_Services_Optional_RandomNumberGenerator_service
		art_root_io::TFileService_service
		art::Framework_Services_System_FileCatalogMetadata_service
		art::Utilities
		canvas::canvas
		messagefacility::MF_MessageLogger
		messagefacility::headers
		fhiclcpp::fhiclcpp
		cetlib::cetlib
		Boost::filesystem
		Boost::system
		CLHEP::Random
)

cet_build_plugin(ROICannyEdgeDetection art::tool LIBRARIES ${TOOL_LIBRARIES})
cet_build_plugin(ROIFromDecoder art::tool LIBRARIES ${TOOL_LIBRARIES})
cet_build_plugin(ROIMorphological2D art::tool LIBRARIES ${TOOL_LIBRARIES})
cet_build_plugin(ROIWavelets art::tool LIBRARIES ${TOOL_LIBRARIES})

install_headers()
install_source()
install_fhicl()
