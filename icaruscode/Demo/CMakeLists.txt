art_make( 
	EXCLUDE CRTSimAnalysis_module.cc
	MODULE_LIBRARIES lardataobj_RecoBase
	  		   ${ROOT_BASIC_LIB_LIST}
                           icaruscode_RecoUtils
        )

simple_plugin( CRTSimAnalysis module
                        icaruscode_CRTData
                        icaruscode_CRT
                        icaruscode_CRTProducts
                        larcore_Geometry_Geometry_service
                        larcorealg_Geometry
                        nusimdata_SimulationBase
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
                        cetlib cetlib_except
                        ${ROOT_BASIC_LIB_LIST}
                        ${ROOT_GEOM}
                        ${ROOT_XMLIO}
                        ${ROOT_GDML}              

)

install_fhicl()
install_source()
