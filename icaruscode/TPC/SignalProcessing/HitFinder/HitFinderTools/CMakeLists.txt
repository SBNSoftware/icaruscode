cet_enable_asserts()

set( hitfinder_tool_lib_list
			larcorealg_Geometry
			lardataobj_RecoBase
			larcore_Geometry_Geometry_service
			lardata_Utilities
			nurandom::RandomUtils_NuRandomService_service
  			art::Framework_Core
			art::Framework_Principal
			art::Framework_Services_Registry
			art::Persistency_Common
			art::Persistency_Provenance
			${ART_FRAMEWORK_SERVICES_BASIC}
			art_root_io::tfile_support ROOT::Core
            art::Framework_Services_Optional_RandomNumberGenerator_service
			art_root_io::TFileService_service
			${ART_FRAMEWORK_SERVICES_SYSTEM_FILECATALOGMETADATA_SERVICE}
			art::Utilities
			canvas::canvas
			messagefacility::MF_MessageLogger
			
			fhiclcpp::fhiclcpp
			cetlib::cetlib
			Boost::filesystem
			
            CLHEP::CLHEP
			${ROOT_BASIC_LIB_LIST}
    )

art_make(TOOL_LIBRARIES ${hitfinder_tool_lib_list}
        )

install_headers()
install_source()
install_fhicl()
