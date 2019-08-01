art_make( 
          LIB_LIBRARIES 
                          lardataobj_RawData
                          lardataobj_RecoBase
                          lardata_Utilities
                          ${FHICLCPP}
                          cetlib cetlib_except
                          ${ROOT_BASIC_LIB_LIST}
                          larcorealg_Geometry
                          larcore_Geometry_Geometry_service
                          lardata_Utilities
                          larevt_Filters
                          lardataobj_RawData
                          larevt_CalibrationDBI_IOVData
                          larevt_CalibrationDBI_Providers
                          lardataobj_RecoBase
                          icaruscode_Utilities_SignalShapingServiceICARUS_service
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
          TOOL_LIBRARIES  icaruscode_TPC_SignalProcessing_RawDigitFilter_Algorithms
                          lardataobj_RawData
                          lardataobj_RecoBase
                          lardata_Utilities
                          ${FHICLCPP}
                          cetlib cetlib_except
                          larcorealg_Geometry
                          larcore_Geometry_Geometry_service
                          larevt_CalibrationDBI_IOVData
                          larevt_CalibrationDBI_Providers
                          icaruscode_Utilities_SignalShapingServiceICARUS_service
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
                          ${ROOT_FFTW}
                          ${ROOT_BASIC_LIB_LIST}
       )

install_headers()
install_fhicl()
install_source()
