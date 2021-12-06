cet_enable_asserts()

art_make( MODULE_LIBRARIES lardataobj_RecoBase
                           larcorealg_Geometry
                           larcore_Geometry_Geometry_service
                           lardata_RecoObjects
                           lardata_Utilities
                           larsim_MCCheater_BackTrackerService_service
                           larsim_MCCheater_ParticleInventoryService_service
                           icarusalg_gallery_MCTruthBase
                           nug4::ParticleNavigation
                           art::Framework_Core
                           art::Framework_Principal
                           art::Framework_Services_Registry
                           ${ART_FRAMEWORK_SERVICES_BASIC}
                           art_root_io::tfile_support ROOT::Core
                           art_root_io::TFileService_service
                           art::Framework_Services_Optional_RandomNumberGenerator_service
                           art::Persistency_Common
                           art::Persistency_Provenance
                           art::Utilities
                           canvas::canvas
                           messagefacility::MF_MessageLogger
                           
                           fhiclcpp::fhiclcpp
                           cetlib::cetlib
                           cetlib_except::cetlib_except
                           Boost::filesystem
                           
                           CLHEP::CLHEP
                           ${ROOT_BASIC_LIB_LIST}

          TOOL_LIBRARIES   lardataobj_RecoBase
                           larcorealg_Geometry
                           larcore_Geometry_Geometry_service
                           lardata_RecoObjects
                           lardata_Utilities
                           larsim_MCCheater_BackTrackerService_service
                           larsim_MCCheater_ParticleInventoryService_service
                           icarusalg_gallery_MCTruthBase
                           nug4::ParticleNavigation
                           art::Framework_Core
                           art::Framework_Principal
                           art::Framework_Services_Registry
                           ${ART_FRAMEWORK_SERVICES_BASIC}
                           art_root_io::tfile_support ROOT::Core
                           art_root_io::TFileService_service
                           art::Framework_Services_Optional_RandomNumberGenerator_service
                           art::Persistency_Common
                           art::Persistency_Provenance
                           art::Utilities
                           canvas::canvas
                           messagefacility::MF_MessageLogger
                           
                           fhiclcpp::fhiclcpp
                           cetlib::cetlib
                           cetlib_except::cetlib_except
                           Boost::filesystem
                           
                           CLHEP::CLHEP
                           ${ROOT_BASIC_LIB_LIST}
        )

install_headers()
install_source()
install_fhicl()
