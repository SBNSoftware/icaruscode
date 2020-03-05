add_subdirectory(HitFinderTools)

art_make( 
          LIB_LIBRARIES 
			           lardataobj_RawData
                       lardataobj_RecoBase
	                   lardata_Utilities
	                   ${FHICLCPP}
			           ${CETLIB}
			           ${ROOT_BASIC_LIB_LIST}
	  MODULE_LIBRARIES larcorealg_Geometry
	  		           larcore_Geometry_Geometry_service
	                   lardata_Utilities
			           larevt_Filters
			           lardataobj_RawData
			           larevt_CalibrationDBI_IOVData
			           larevt_CalibrationDBI_Providers
			           lardataobj_RecoBase
			           lardata_ArtDataHelper
			           larreco_RecoAlg
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
                       ${ROOT_GEOM}
                       ${ROOT_XMLIO}
                       ${ROOT_GDML}
			           ${ROOT_FFTW}
			           ${ROOT_BASIC_LIB_LIST}
        )

install_headers()
install_fhicl()
install_source()

