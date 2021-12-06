if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
cet_add_compiler_flags(CXX -Wno-maybe-uninitialized)
endif()

# Add the subfolder containing the deconvolution tools
add_subdirectory(DeconTools)
add_subdirectory(ROITools)
add_subdirectory(SimTestPulse)

add_definitions(-DEIGEN_FFTW_DEFAULT)

art_make( 
          LIB_LIBRARIES 
            lardataobj_RawData
            lardataobj_RecoBase
            lardata_Utilities
            fhiclcpp::fhiclcpp
            cetlib::cetlib
            ${ROOT_BASIC_LIB_LIST}
          MODULE_LIBRARIES 
            icarus_signal_processing
            icarus_signal_processing_Detection
            icarus_signal_processing_Filters
            larcorealg_Geometry
            larcore_Geometry_Geometry_service
            lardata_Utilities
            lardataalg_DetectorInfo
            larevt_Filters
            lardataobj_RawData
            larevt_CalibrationDBI_IOVData
            larevt_CalibrationDBI_Providers
            lardataobj_RecoBase
            lardata_ArtDataHelper
            icaruscode_TPC_Utilities_SignalShapingICARUSService_service
            art::Framework_Core
            art::Framework_Principal
            art::Framework_Services_Registry
            art_root_io::tfile_support ROOT::Core
            art_root_io::TFileService_service
            art::Persistency_Common
            art::Persistency_Provenance
            art::Utilities
            canvas::canvas
            messagefacility::MF_MessageLogger
            
            fhiclcpp::fhiclcpp
            cetlib::cetlib
            ROOT::FFTW
            ${ROOT_BASIC_LIB_LIST}
            ${TBB}
        )

install_headers()
install_fhicl()
install_source()
