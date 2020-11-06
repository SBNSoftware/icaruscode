# Removing AVX instructions because assuming the grid nodes
# being less than 6 year old proved to be pretentious.
# AVX is not so relevant after the last commit anyway
# cet_set_compiler_flags(EXTRA_CXX_FLAGS -mavx)

add_subdirectory(CaloTools)

simple_plugin(PMTCalibration "module"
                icaruscode_PMT_Algorithms
                icaruscode_PMT_Calo_CaloTools
                lardataobj_RawData
                lardataobj_Simulation
                lardataobj_RecoBase
                nurandom_RandomUtils_NuRandomService_service
                ${ART_FRAMEWORK_CORE}
                ${ART_FRAMEWORK_PRINCIPAL}
                ${ART_FRAMEWORK_SERVICES_REGISTRY}
                ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                art_Persistency_Provenance
                ${CANVAS}
                ${MF_MESSAGELOGGER}
                ${FHICLCPP}
                ${CETLIB}
                cetlib_except
                ${ROOT_BASIC_LIB_LIST}
                ${Boost_SYSTEM_LIBRARY}
              )

simple_plugin(PMTPulseAna "module"
                icaruscode_PMT_Algorithms
                icaruscode_PMT_Calo_CaloTools
                lardataobj_RawData
                lardataobj_Simulation
                nurandom_RandomUtils_NuRandomService_service
                ${ART_FRAMEWORK_CORE}
                ${ART_FRAMEWORK_PRINCIPAL}
                ${ART_FRAMEWORK_SERVICES_REGISTRY}
                ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                art_Persistency_Provenance
                ${CANVAS}
                ${MF_MESSAGELOGGER}
                ${FHICLCPP}
                ${CETLIB}
                cetlib_except
                ${ROOT_BASIC_LIB_LIST}
                ${Boost_SYSTEM_LIBRARY}
              )

simple_plugin(PMTDump "module"
                icaruscode_PMT_Algorithms
                icaruscode_PMT_Calo_CaloTools
                lardataobj_RawData
                lardataobj_Simulation
                nurandom_RandomUtils_NuRandomService_service
                ${ART_FRAMEWORK_CORE}
                ${ART_FRAMEWORK_PRINCIPAL}
                ${ART_FRAMEWORK_SERVICES_REGISTRY}
                ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                art_Persistency_Provenance
                ${CANVAS}
                ${MF_MESSAGELOGGER}
                ${FHICLCPP}
                ${CETLIB}
                cetlib_except
                ${ROOT_BASIC_LIB_LIST}
                ${Boost_SYSTEM_LIBRARY}
              )

install_headers()
install_fhicl()
install_source()
