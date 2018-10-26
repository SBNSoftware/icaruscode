
art_make(
	  MODULE_LIBRARIES
			larcorealg_Geometry
	                ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_IO_ROOT}
			${ART_FRAMEWORK_MODULES}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			${BOOST_LIB}
                        ${FHICLCPP}
	  LIB_LIBRARIES larcorealg_Geometry
                        ${MF_MESSAGELOGGER}
                        ${CETLIB}
			${CETLIB_EXCEPT}
                        ${FHICLCPP}
                        ${ROOT_GEOM}
                        ${ROOT_BASIC_LIB_LIST}
	  SERVICE_LIBRARIES icaruscode_Geometry
	  		icaruscode_Geometry
			${MF_UTILITIES}
	                ${FHICLCPP}
			${CETLIB}
			${CETLIB_EXCEPT}
                        ${ROOT_GEOM}
                        ${ROOT_XMLIO}
                        ${ROOT_GDML}
			${ROOT_BASIC_LIB_LIST}
        )							   



add_subdirectory(gdml)

install_headers()
install_fhicl()
install_source()


