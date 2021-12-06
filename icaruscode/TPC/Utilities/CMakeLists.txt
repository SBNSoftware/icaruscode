add_subdirectory(tools)

add_definitions(-DEIGEN_FFTW_DEFAULT)

set( icarus_util_lib_list
                        larcorealg_Geometry
                        larevt_CalibrationDBI_IOVData
                        larevt_CalibrationDBI_Providers
                        larreco_HitFinder
                        lardataobj_RecoBase
                        larcore_Geometry_Geometry_service
                        lardata_Utilities
                        nurandom::RandomUtils_NuRandomService_service
                        ${ICARUS_FFTW_LIBRARIES}
                        art::Framework_Core
                        art::Framework_Principal
                        art::Framework_Services_Registry
                        art::Persistency_Common
                        art_Framework_IO
                        art_Framework_IO_detail
                        art::Persistency_Provenance
                        ${ART_FRAMEWORK_SERVICES_BASIC}
                        art_root_io::tfile_support ROOT::Core
                        art::Framework_Services_Optional_RandomNumberGenerator_service
                        art_root_io::TFileService_service
                        ${ART_FRAMEWORK_SERVICES_SYSTEM_FILECATALOGMETADATA_SERVICE}
                        art_Framework_Services_System_TriggerNamesService_service
                        art::Utilities
                        canvas::canvas
                        messagefacility::MF_MessageLogger
                        
                        fhiclcpp::fhiclcpp
                        cetlib::cetlib
                        Boost::filesystem
                        
                        CLHEP::CLHEP
                        ${ROOT_BASIC_LIB_LIST}
    )

art_make( 
#          LIBRARY_NAME  icaruscode_TPC_Utilities
#          LIB_LIBRARIES ${icarus_util_lib_list}
#
          SERVICE_LIBRARIES 
                        ${icarus_util_lib_list}

          MODULE_LIBRARIES
                        art_root_io::TFileService_service
                        art_root_io::tfile_support ROOT::Core
                        art::Framework_Services_Registry
                        art::Framework_Principal
                        art::Framework_Core
                        canvas::canvas
                        messagefacility::MF_MessageLogger
                        hep_concurrency
                        fhiclcpp::fhiclcpp
                        cetlib::cetlib
                        ROOT::Core
       )


install_headers()
install_source()
install_fhicl()
