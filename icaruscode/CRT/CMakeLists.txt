include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
include_directories( $ENV{NURANDOM_INC} )
include_directories( $ENV{IFDHX_FQ_DIR}/inc )

add_subdirectory(CRTProducts)
add_subdirectory(CRTUtils)
add_subdirectory(CRTDecoder)

art_make(
    NO_PLUGINS
    EXCLUDE
        CRTChannelMapAlg.cxx
        CRTGeometryHelper_service.cc
        CRTDetSim_module.cc
        CRTSimHitProducer_module.cc
        CRTTrueHitProducer_module.cc
        CRTTzeroProducer_module.cc
        CRTTrackProducer_module.cc
        CRTSimAnalysis_module.cc
        CRTDataAnalysis_module.cc
	      CrtOpHitMatchAnalysis_module.cc
        CRTTruthMatchAnalysis_module.cc
        CRTAutoVeto_module.cc
        FlashResAna_module.cc
    LIBRARY_NAME
        icaruscode_CRTData
    LIB_LIBRARIES
        larcorealg_Geometry
        larcore_Geometry_Geometry_service
        ${ART_FRAMEWORK_CORE}
        ${ART_FRAMEWORK_IO_SOURCES}
        ${ART_FRAMEWORK_PRINCIPAL}
        ${ART_ROOT_IO_TFILESERVICE_SERVICE}
        canvas
        ${ART_FRAMEWORK_SERVICES_REGISTRY}
        ${MF_MESSAGELOGGER}
        ${MF_UTILITIES}
        ${FHICLCPP}
        ${ROOT_BASIC_LIB_LIST}
        ${ROOT_GEOM}
        ${CETLIB}
)

art_make_library(
    LIBRARY_NAME
        icaruscode_CRT
    SOURCE
        CRTChannelMapAlg.cxx
    LIBRARIES larcorealg_Geometry
        icaruscode_CRTData
        cetlib_except
        ${MF_MESSAGELOGGER}
        ${MF_UTILITIES}
        ${FHICLCPP}
        ${ROOT_BASIC_LIB_LIST}
        ${ROOT_GEOM}
        ${CETLIB}
)

simple_plugin(CRTGeometryHelper service
              larcorealg_Geometry
              icaruscode_CRT
              ${ART_ROOT_IO_TFILESERVICE_SERVICE}
              ${ART_FRAMEWORK_SERVICES_REGISTRY}
              ${MF_MESSAGELOGGER}
              ${MF_UTILITIES}
              ${ROOT_BASIC_LIB_LIST}
)

simple_plugin(CRTDetSim module
              larcorealg_Geometry
              icaruscode_CRT
              icaruscode_CRTData
              icaruscode_CRTProducts
              icaruscode_CRTUtils
              ${ART_ROOT_IO_TFILESERVICE_SERVICE}
              lardataalg_DetectorInfo
              nurandom_RandomUtils_NuRandomService_service
              ${ART_FRAMEWORK_SERVICES_REGISTRY}
              ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
              ${MF_MESSAGELOGGER}
              ${MF_UTILITIES}
              ${ROOT_BASIC_LIB_LIST}
              ${CLHEP}
)

simple_plugin(CRTSimHitProducer module
              larcorealg_Geometry
              icaruscode_CRT
              icaruscode_CRTData
              icaruscode_CRTProducts
              icaruscode_CRTUtils
              ${ART_ROOT_IO_TFILESERVICE_SERVICE}
              lardataalg_DetectorInfo
              nurandom_RandomUtils_NuRandomService_service
              ${ART_FRAMEWORK_SERVICES_REGISTRY}
              ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
              ${MF_MESSAGELOGGER}
              ${MF_UTILITIES}
              ${ROOT_BASIC_LIB_LIST}
              ${CLHEP}
)

simple_plugin( CRTTrueHitProducer module
        larcorealg_Geometry
        icaruscode_CRT
        icaruscode_CRTData
        icaruscode_CRTProducts
        icaruscode_CRTUtils
        ${ART_FRAMEWORK_CORE}
        ${ART_FRAMEWORK_IO_SOURCES}
        ${ART_FRAMEWORK_PRINCIPAL}
        ${ART_ROOT_IO_TFILESERVICE_SERVICE}
        art_Persistency_Provenance
        canvas
        cetlib_except
        lardata_DetectorInfoServices_DetectorClocksServiceStandard_service
        ${ART_FRAMEWORK_SERVICES_REGISTRY}
        ${MF_MESSAGELOGGER}
        ${MF_UTILITIES}
        ${ROOT_BASIC_LIB_LIST}
        ${CLHEP}
        ${CETLIB}
)

simple_plugin(CRTTzeroProducer module
              larcorealg_Geometry
              icaruscode_CRT
              icaruscode_CRTData
              icaruscode_CRTProducts
              ${ART_ROOT_IO_TFILESERVICE_SERVICE}
              lardataalg_DetectorInfo
              nurandom_RandomUtils_NuRandomService_service
              ${ART_FRAMEWORK_SERVICES_REGISTRY}
              ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
              ${MF_MESSAGELOGGER}
              ${MF_UTILITIES}
              ${ROOT_BASIC_LIB_LIST}
              ${CLHEP}
)

simple_plugin(CRTTrackProducer module
              larcorealg_Geometry
              icaruscode_CRT
              icaruscode_CRTData
              icaruscode_CRTProducts
              icaruscode_CRTUtils
              ${ART_ROOT_IO_TFILESERVICE_SERVICE}
              lardataalg_DetectorInfo
              nurandom_RandomUtils_NuRandomService_service
              ${ART_FRAMEWORK_SERVICES_REGISTRY}
              ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
              ${MF_MESSAGELOGGER}
              ${MF_UTILITIES}
              ${ROOT_BASIC_LIB_LIST}
              ${CLHEP}
)

simple_plugin(CRTSimAnalysis module
              icaruscode_CRTData
              icaruscode_CRT
              icaruscode_CRTProducts
              icaruscode_CRTUtils
              larcore_Geometry_Geometry_service
              larcorealg_Geometry
              nusimdata_SimulationBase
              ${ART_FRAMEWORK_SERVICES_REGISTRY}
              ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
              ${ART_ROOT_IO_TFILESERVICE_SERVICE}
              ${MF_MESSAGELOGGER}
              ${MF_UTILITIES}
              ${ROOT_BASIC_LIB_LIST}
              ${ROOT_GEOM}
              ${ROOT_XMLIO}
              ${ROOT_GDML}
)

simple_plugin(CRTDataAnalysis module
              icaruscode_CRTData
              icaruscode_CRT
              icaruscode_CRTProducts
              larcore_Geometry_Geometry_service
              larcorealg_Geometry
              nusimdata_SimulationBase
              ${ART_FRAMEWORK_SERVICES_REGISTRY}
              ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
              ${ART_ROOT_IO_TFILESERVICE_SERVICE}
              ${MF_MESSAGELOGGER}
              ${MF_UTILITIES}
              ${ROOT_BASIC_LIB_LIST}
              ${ROOT_GEOM}
              ${ROOT_XMLIO}
              ${ROOT_GDML}
)

simple_plugin(CrtOpHitMatchAnalysis module
              icaruscode_CRTData
              icaruscode_CRT
              icaruscode_CRTProducts
              icaruscode_CRTUtiles
              larcore_Geometry_Geometry_service
              larcorealg_Geometry
              nusimdata_SimulationBase
              lardataobj_RecoBase
              ${ART_FRAMEWORK_SERVICES_REGISTRY}
              ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
              ${ART_ROOT_IO_TFILESERVICE_SERVICE}
              lardataalg_DetectorInfo
              ${MF_MESSAGELOGGER}
              ${MF_UTILITIES}
              ${ROOT_BASIC_LIB_LIST}
)

simple_plugin(CRTTruthMatchAnalysis module
              icaruscode_CRTData
              icaruscode_CRT
              icaruscode_CRTProducts
              icaruscode_CRTUtils
              larcore_Geometry_Geometry_service
              larcorealg_Geometry
              nusimdata_SimulationBase
              ${ART_FRAMEWORK_SERVICES_REGISTRY}
              ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
              ${ART_ROOT_IO_TFILESERVICE_SERVICE}
              ${MF_MESSAGELOGGER}
              ${MF_UTILITIES}
              ${ROOT_BASIC_LIB_LIST}
              ${ROOT_GEOM}
              ${ROOT_XMLIO}
              ${ROOT_GDML}
)

simple_plugin( CRTAutoVeto module
                        icaruscode_CRTData
                        icaruscode_CRT
                        icaruscode_CRTProducts
                        icaruscode_CRTUtils
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

simple_plugin( FlashResAna module
                        larcore_Geometry_Geometry_service
                        larcorealg_Geometry
                        nusimdata_SimulationBase
                        lardataobj_RecoBase
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


install_headers()
install_fhicl()
install_source()
