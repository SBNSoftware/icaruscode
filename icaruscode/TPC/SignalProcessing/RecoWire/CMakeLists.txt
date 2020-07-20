# Add the subfolder containing the deconvolution tools
add_subdirectory(DeconTools)
add_subdirectory(SimTestPulse)

add_definitions(-DEIGEN_FFTW_DEFAULT)

art_make( 
          LIB_LIBRARIES 
            lardataobj_RawData
            lardataobj_RecoBase
            lardata_Utilities
            ${FHICLCPP}
            ${CETLIB}
            ${ROOT_BASIC_LIB_LIST}
          MODULE_LIBRARIES 
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
            ${ART_FRAMEWORK_CORE}
            ${ART_FRAMEWORK_PRINCIPAL}
            ${ART_FRAMEWORK_SERVICES_REGISTRY}
            ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
            ${ART_ROOT_IO_TFILESERVICE_SERVICE}
            art_Persistency_Common
            art_Persistency_Provenance
            art_Utilities
            canvas
            ${MF_MESSAGELOGGER}
            ${MF_UTILITIES}
            ${FHICLCPP}
            ${CETLIB}
            ${ROOT_FFTW}
            ${ROOT_BASIC_LIB_LIST}
            ${TBB}
        )

install_headers()
install_fhicl()
install_source()
