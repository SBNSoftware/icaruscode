set(EXCLUDED_PLUGINS WireModMakeHists_plugin.cc ROICounter_plugin.cc HitComp_plugin.cc HitTTreeMaker_plugin.cc WireModifier_module.cc)

add_subdirectory(Utility)

art_make( 
	  LIB_LIBRARIES
			${ART_FRAMEWORK_CORE}
			art_Utilities canvas
			${MF_MESSAGELOGGER}
			
			cetlib cetlib_except
                        ${Boost_SYSTEM_LIBRARY}
			${ROOT_BASIC_LIB_LIST}
			${ROOT_HIST}
          MODULE_LIBRARIES 
                        art_root_io::TFileService_service
                        lardataobj::AnalysisBase
			lardataobj::RawData
			lardataobj::RecoBase
			lardata::RecoObjects 
			larreco::RecoAlg
			lardata::Utilities
			${LARRECO_LIB}
			${LARDATA_LIB}
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			${MF_MESSAGELOGGER}
	                ${FHICLCPP}
			${CLHEP}
                        ${ROOT_GEOM}
                        ${ROOT_XMLIO}
                        ${ROOT_GDML}
			${ROOT_BASIC_LIB_LIST}
			${ROOT_HIST}
            EXCLUDE     ${EXCLUDED_PLUGINS}
                        #WireMod_Utilities.hh
                        #WireMod_Utilities.cc
        )

find_library( sbncode::sbncode_WireMod_Utility_LIBRARY NAMES sbncode::sbncode_WireMod_Utility )

simple_plugin(WireModMakeHists "plugin"
              ${ART_FRAMEWORK_CORE}
              ${MF_MESSAGELOGGER}
              ${Boost_SYSTEM_LIBRARY}
              ${ROOT_BASIC_LIB_LIST}
              ${ROOT_HIST}
              art_root_io::TFileService_service
              lardataobj::AnalysisBase
              lardataobj::RawData
              lardataobj::RecoBase
              lardata::RecoObjects
              larreco::RecoAlg
              ${LARRECO_LIB}
              ${LARDATA_LIB}
              ${ART_FRAMEWORK_CORE}
              ${ART_FRAMEWORK_PRINCIPAL}
              ${ART_FRAMEWORK_SERVICES_REGISTRY}
              ${ART_FRAMEWORK_SERVICES_OPTIONAL}
              ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
              ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
              ${MF_MESSAGELOGGER}
              ${FHICLCPP}
              ${CLHEP}
              ${ROOT_GEOM}
              ${ROOT_XMLIO}
              ${ROOT_GDML}
              ${ROOT_BASIC_LIB_LIST}
              ${ROOT_HIST}
              BASENAME_ONLY
)

simple_plugin(HitComp "plugin"
              ${ART_FRAMEWORK_CORE}
              ${MF_MESSAGELOGGER}
              ${Boost_SYSTEM_LIBRARY}
              ${ROOT_BASIC_LIB_LIST}
              ${ROOT_HIST}
              art_root_io::TFileService_service
              lardataobj::AnalysisBase
              lardataobj::RawData
              lardataobj::RecoBase
              lardata::RecoObjects
              larreco::RecoAlg
              ${LARRECO_LIB}
              ${LARDATA_LIB}
              ${ART_FRAMEWORK_CORE}
              ${ART_FRAMEWORK_PRINCIPAL}
              ${ART_FRAMEWORK_SERVICES_REGISTRY}
              ${ART_FRAMEWORK_SERVICES_OPTIONAL}
              ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
              ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
              ${MF_MESSAGELOGGER}
              ${FHICLCPP}
              ${CLHEP}
              ${ROOT_GEOM}
              ${ROOT_XMLIO}
              ${ROOT_GDML}
              ${ROOT_BASIC_LIB_LIST}
              ${ROOT_HIST}
              BASENAME_ONLY
)

simple_plugin(HitTTreeMaker "plugin"
              ${ART_FRAMEWORK_CORE}
              ${MF_MESSAGELOGGER}
              ${Boost_SYSTEM_LIBRARY}
              ${ROOT_BASIC_LIB_LIST}
              ${ROOT_HIST}
              art_root_io::TFileService_service
              lardataobj::AnalysisBase
              lardataobj::RawData
              lardataobj::RecoBase
              lardata::RecoObjects
              larreco::RecoAlg
              ${LARRECO_LIB}
              ${LARDATA_LIB}
              ${ART_FRAMEWORK_CORE}
              ${ART_FRAMEWORK_PRINCIPAL}
              ${ART_FRAMEWORK_SERVICES_REGISTRY}
              ${ART_FRAMEWORK_SERVICES_OPTIONAL}
              ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
              ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
              ${MF_MESSAGELOGGER}
              ${FHICLCPP}
              ${CLHEP}
              ${ROOT_GEOM}
              ${ROOT_XMLIO}
              ${ROOT_GDML}
              ${ROOT_BASIC_LIB_LIST}
              ${ROOT_HIST}
              BASENAME_ONLY
)

simple_plugin(ROICounter "plugin"
              ${ART_FRAMEWORK_CORE}
              ${MF_MESSAGELOGGER}
              ${Boost_SYSTEM_LIBRARY}
              ${ROOT_BASIC_LIB_LIST}
              ${ROOT_HIST}
              art_root_io::TFileService_service
              lardataobj::AnalysisBase
              lardataobj::RawData
              lardataobj::RecoBase
              lardata::RecoObjects
              larreco::RecoAlg
              ${LARRECO_LIB}
              ${LARDATA_LIB}
              ${ART_FRAMEWORK_CORE}
              ${ART_FRAMEWORK_PRINCIPAL}
              ${ART_FRAMEWORK_SERVICES_REGISTRY}
              ${ART_FRAMEWORK_SERVICES_OPTIONAL}
              ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
              ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
              ${MF_MESSAGELOGGER}
              ${FHICLCPP}
              ${CLHEP}
              ${ROOT_GEOM}
              ${ROOT_XMLIO}
              ${ROOT_GDML}
              ${ROOT_BASIC_LIB_LIST}
              ${ROOT_HIST}
              BASENAME_ONLY
)

simple_plugin(WireModifier "module"
              ${ART_FRAMEWORK_CORE}
              ${MF_MESSAGELOGGER}
              ${Boost_SYSTEM_LIBRARY}
              ${ROOT_BASIC_LIB_LIST}
              ${ROOT_HIST}
              art_root_io::TFileService_service
              lardataobj::AnalysisBase
              lardataobj::RawData
              lardataobj::RecoBase
              lardata::RecoObjects
              larreco::RecoAlg
              ${LARRECO_LIB}
              ${LARDATA_LIB}
              ${ART_FRAMEWORK_CORE}
              ${ART_FRAMEWORK_PRINCIPAL}
              ${ART_FRAMEWORK_SERVICES_REGISTRY}
              ${ART_FRAMEWORK_SERVICES_OPTIONAL}
              ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
              ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
              ${MF_MESSAGELOGGER}
              ${FHICLCPP}
              ${CLHEP}
              ${ROOT_GEOM}
              ${ROOT_XMLIO}
              ${ROOT_GDML}
              ${ROOT_BASIC_LIB_LIST}
              ${ROOT_HIST}
              larcore::Geometry_Geometry_service
              larcorealg::Geometry
              lardataalg::DetectorInfo
              lardata::Utilities
              lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
              sbnobj::ICARUS_TPC
              BASENAME_ONLY
)

install_headers()
install_fhicl()
install_source()

