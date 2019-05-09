include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

set( ROOTLIB -L$ENV{ROOTSYS}/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl)

link_libraries( ${LIB_NAME} -L$ENV{BOOST_LIB} -lboost_system ${ROOTLIB} )

art_make( 
	  LIBRARY_NAME  icarus_FlashFinder
	  LIB_LIBRARIES
	  	   #larcore_Geometry
		   icaruscode_Geometry
		   larcore_Geometry_Geometry_service
                   #UBOpReadoutMap_service
		   lardataobj_RecoBase
		   #${UBOONECODE_LIB}
		   ${LARDATA_LIB}
		   canvas
		   cetlib
		   cetlib_except
		   ${ART_FRAMEWORK_ART}
		   ${ART_FRAMEWORK_CORE}
		   ${ART_FRAMEWORK_PRINCIPAL}
		   ${ART_FRAMEWORK_UTILITIES}
                   ${ART_FRAMEWORK_BASIC}
		   ${ART_FRAMEWORK_SERVICES_REGISTRY}
		   ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
		   ${ART_ROOT_IO_TFILESERVICE_SERVICE}
		   art_Persistency_Common
		   art_Persistency_Provenance
		   art_Utilities
		   ${MF_MESSAGELOGGER}
		   ${MF_UTILITIES}
		   ${FHICLCPP}
		   ${ROOT_GEOM}
		   ${ROOT_XMLIO}
		   ${ROOT_GDML}
		   ${ROOT_BASIC_LIB_LIST}
          MODULE_LIBRARIES
		   icarus_FlashFinder
	  	   #larcore_Geometry
		   icaruscode_Geometry
		   larcore_Geometry_Geometry_service
		   lardataobj_RecoBase
		   #${UBOONECODE_LIB}
		   ${LARDATA_LIB}
		   canvas
		   cetlib
		   cetlib_except
		   ${ART_FRAMEWORK_ART}
		   ${ART_FRAMEWORK_CORE}
		   ${ART_FRAMEWORK_PRINCIPAL}
                   ${ART_FRAMEWORK_BASIC}
		   ${ART_FRAMEWORK_SERVICES_REGISTRY}
		   ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
		   ${ART_ROOT_IO_TFILESERVICE_SERVICE}
		   art_Persistency_Common
		   art_Persistency_Provenance
		   art_Utilities
		   ${MF_MESSAGELOGGER}
		   ${MF_UTILITIES}
		   ${FHICLCPP}
		   ${ROOT_GEOM}
		   ${ROOT_XMLIO}
		   ${ROOT_GDML}
		   ${ROOT_BASIC_LIB_LIST}

        )							   

install_headers()
install_source()

