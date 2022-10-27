cet_enable_asserts()

art_make_library()
set(		TOOL_LIBRARIES
			icaruscode::TPC_Utilities_SignalShapingICARUSService_service
			larcorealg::Geometry
			larevt::CalibrationDBI_IOVData
			larevt::CalibrationDBI_Providers
			lardataobj::RecoBase
			larcore::Geometry_Geometry_service
			lardata::Utilities
			nurandom::RandomUtils_NuRandomService_service
			art::Framework_Core
			art::Framework_Principal
			art::Framework_Services_Registry
			art::Persistency_Common
			art::Persistency_Provenance
			art_root_io::tfile_support
			ROOT::Core
			art::Framework_Services_Optional_RandomNumberGenerator_service
			art_root_io::TFileService_service
			art::Framework_Services_System_FileCatalogMetadata_service
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
)
cet_build_plugin(OpHitFinderStandard art::tool LIBRARIES ${TOOL_LIBRARIES})
cet_build_plugin(OpHitFinder art::tool LIBRARIES ${TOOL_LIBRARIES})

install_headers()
install_source()
install_fhicl()
