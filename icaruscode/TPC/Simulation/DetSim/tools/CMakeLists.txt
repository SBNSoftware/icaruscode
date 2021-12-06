cet_enable_asserts()

add_definitions(-DEIGEN_FFTW_DEFAULT)

art_make( 
          LIB_LIBRARIES      art::Utilities
                             canvas::canvas
                             fhiclcpp::fhiclcpp
                             messagefacility::MF_MessageLogger
							 ${ART_FRAMEWORK_SERVICES_BASIC}
							 art::Framework_Services_Registry
							 art_root_io::tfile_support ROOT::Core
							 cetlib_except::cetlib_except
                             cetlib::cetlib
							 
		  SERVICE_LIBRARIES  icaruscode_TPC_Simulation_DetSim_tools
		  					 icaruscode_Decode_ChannelMapping
						     lardata_Utilities
						     lardata_ArtDataHelper
						     ${ROOT_BASIC_LIB_LIST}
						     LIBWDA
						     art::Framework_Core
						     art::Framework_Principal
						     art::Framework_Services_Registry
						     art_root_io::tfile_support
						     art_root_io::TFileService_service
						     art::Utilities
						     art::Persistency_Common canvas
						     art::Persistency_Provenance canvas
						     art::Utilities canvas
						     messagefacility::MF_MessageLogger
						     
						     fhiclcpp::fhiclcpp
						     cetlib::cetlib cetlib_except

		  TOOL_LIBRARIES     larcorealg_Geometry
			                 larevt_CalibrationDBI_IOVData
			                 larevt_CalibrationDBI_Providers
			                 larreco_HitFinder
			                 lardataobj_RecoBase
			                 larcore_Geometry_Geometry_service
			                 lardata_Utilities
                             lardataalg_DetectorInfo
			                 nurandom::RandomUtils_NuRandomService_service
                             ${ICARUS_FFTW_LIBRARIES}
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

install_headers()
install_source()
install_fhicl()
