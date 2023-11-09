cet_enable_asserts()

#find_library(ICARUS_OVERLAY_LIB sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_ICARUS)
#if (NOT ICARUS_OVERLAY_LIB)
#  message(FATAL_ERROR "sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_ICARUS not found!")
#endif()

art_make_library(
        EXCLUDE
                  ProduceCompressedFragment_module.cc
                  ValidateFragmentCompression_plugin.cc
                  TestProduceCompressed_module.cc
        LIBRARIES
                  artdaq_core::artdaq-core_Data
                  ${CETLIB_EXCEPT}
                  sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_ICARUS
                )

set(MODULE_LIBRARIES
                   sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_ICARUS
                   sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays
                   artdaq_core::artdaq-core_Utilities
                   art_root_io::TFileService_service
                   ${ROOT_BASIC_LIB_LIST}
                   ${ART_FRAMEWORK_CORE}
                   ${ART_FRAMEWORK_PRINCIPAL}
                   ${ART_FRAMEWORK_SERVICES_REGISTRY}
                   ${ART_ROOT_IO_TFILE_SUPPORT}
                   ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                   ${ART_UTILITIES}
                   ${MF_MESSAGELOGGER}
                   ${MF_UTILITIES}
                   ${FHICLCPP}
                   )

simple_plugin(ValidateFragmentCompression "plugin"
                sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_ICARUS
                artdaq_core::artdaq-core_Utilities
                art_root_io::TFileService_service
                ${ROOT_BASIC_LIB_LIST}
                ${ROOT_HIST}
                ${ART_FRAMEWORK_CORE}
                ${ART_FRAMEWORK_PRINCIPAL}
                ${ART_FRAMEWORK_SERVICES_REGISTRY}
                ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                ${ART_ROOT_IO_TFILE_SUPPORT}
                ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                ${ART_UTILITIES}
                ${Boost_SYSTEM_LIBRARY}
                ${MF_MESSAGELOGGER}
                ${MF_UTILITIES}
                ${FHICLCPP}
                ${CLHEP}
                ROOT::Gdml
                ROOT::XMLIO
                ROOT::Geom
                ROOT::Tree
                ROOT::Core
                messagefacility::headers
                messagefacility::MF_MessageLogger
                art_root_io::tfile_support
                art::Framework_Services_Registry
             )

cet_build_plugin(ProduceCompressedFragment art::module LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(TestProduceCompressed art::module LIBRARIES ${MODULE_LIBRARIES})

install_headers()
install_source()
install_fhicl()
