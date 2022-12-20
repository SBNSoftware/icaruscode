cet_enable_asserts()

art_make_library(
        EXCLUDE
                  ProduceCompressedFragment_module.cc
        LIBRARIES
                  artdaq_core::artdaq-core_Data
                  ${CETLIB_EXCEPT}
                  sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_ICARUS
                )

simple_plugin(ProduceCompressedFragment module
                sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_ICARUS
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

install_headers()
install_source()
install_fhicl()
