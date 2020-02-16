

art_make(   LIB_LIBRARIES
                lardataobj_RawData
                lardataobj_RecoBase
                lardata_Utilities
                ${FHICLCPP}
                ${CETLIB}
                ${ROOT_BASIC_LIB_LIST}
            MODULE_LIBRARIES
                icaruscode_TPC_SignalProcessing_RecoWire_SimTestPulse
                icaruscode_Utilities_SignalShapingICARUSService_service
                icaruscode_Utilities
                lardata_Utilities_LArFFT_service
                lardata_Utilities
                lardata_DetectorInfoServices_DetectorClocksServiceStandard_service
                lardataobj_Simulation
                larcoreobj_SummaryData
                larcorealg_Geometry
                larevt_CalibrationDBI_IOVData
                larevt_CalibrationDBI_Providers
                larreco_HitFinder
                lardataobj_RecoBase
                larcorealg_Geometry
                larcore_Geometry_Geometry_service
                lardata_Utilities
                nurandom_RandomUtils_NuRandomService_service
                lardata_Utilities
                lardataobj_RawData
                larevt_CalibrationDBI_IOVData
                larevt_CalibrationDBI_Providers
                lardataobj_RecoBase
                lardata_ArtDataHelper
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
)

install_fhicl()

