add_subdirectory(Algorithms)

add_definitions(-DEIGEN_FFTW_DEFAULT)

art_make(
          LIB_LIBRARIES 
                        lardataobj_RawData
                        lardataobj_RecoBase
                        lardata_Utilities
                        ${FHICLCPP}
                        cetlib cetlib_except
                        ${ROOT_BASIC_LIB_LIST}
          MODULE_LIBRARIES
                        larcorealg_Geometry
                        icaruscode_TPC_SignalProcessing_RawDigitFilter_Algorithms
                        larcore_Geometry_Geometry_service
                        lardata_Utilities
                        larevt_Filters
                        lardataobj_RawData
                        lardataalg_DetectorInfo
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
                        art_Persistency_Common canvas
                        art_Persistency_Provenance canvas
                        art_Utilities canvas
                        ${MF_MESSAGELOGGER}
                        ${MF_UTILITIES}
                        ${FHICLCPP}
                        cetlib cetlib_except
                        ${ROOT_GEOM}
                        ${ROOT_XMLIO}
                        ${ROOT_GDML}
			${ROOT_FFTW}
			${ROOT_BASIC_LIB_LIST}
			${TBB}
        )

install_headers()
install_fhicl()
install_source()
