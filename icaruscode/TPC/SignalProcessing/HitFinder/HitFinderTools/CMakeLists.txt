cet_enable_asserts()

set(hitfinder_tool_lib_list
			larcorealg::Geometry
			lardataobj::RecoBase
			larcore::Geometry_Geometry_service
			lardata::Utilities
			larreco::HitFinder
			nurandom::RandomUtils_NuRandomService_service
			art::Framework_Core
			art::Framework_Principal
			art::Framework_Services_Registry
			art::Persistency_Common
			art::Persistency_Provenance		
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

art_make(TOOL_LIBRARIES ${hitfinder_tool_lib_list})

install_headers()
install_source()
install_fhicl()
