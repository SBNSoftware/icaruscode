cet_enable_asserts()

art_make_library(
        EXCLUDE
                  ValidateFragmentCompression_plugin.cc
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
                   icaruscode::Decode_DecoderTools_Dumpers
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

cet_build_plugin(ICARUSProduceCompressed art::module LIBRARIES ${MODULE_LIBRARIES})

### These plugins and modules are to be used when producing or validating the reprocessing
### of TPC Fragments with compression and turning them into compressed fragments
### this requires a version of sbndaq_artdaq_core which has the updated PhysCrateOverlays
### remove these comments if you are building in such an enviroment

#find_library(ICARUS_OVERLAY_LIB sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_ICARUS)
#if (NOT ICARUS_OVERLAY_LIB)
#  message(FATAL_ERROR "sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_ICARUS not found!")
#endif()

#simple_plugin(ValidateFragmentCompression "plugin"
#                sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_ICARUS
#                artdaq_core::artdaq-core_Utilities
#                art_root_io::TFileService_service
#                ${ROOT_BASIC_LIB_LIST}
#                ${ROOT_HIST}
#                ${ART_FRAMEWORK_CORE}
#                ${ART_FRAMEWORK_PRINCIPAL}
#                ${ART_FRAMEWORK_SERVICES_REGISTRY}
#                ${ART_FRAMEWORK_SERVICES_OPTIONAL}
#                ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
#                ${ART_ROOT_IO_TFILE_SUPPORT}
#                ${ART_ROOT_IO_TFILESERVICE_SERVICE}
#                ${ART_UTILITIES}
#                ${Boost_SYSTEM_LIBRARY}
#                ${MF_MESSAGELOGGER}
#                ${MF_UTILITIES}
#                ${FHICLCPP}
#                ${CLHEP}
#                ROOT::Gdml
#                ROOT::XMLIO
#                ROOT::Geom
#                ROOT::Tree
#                ROOT::Core
#                messagefacility::headers
#                messagefacility::MF_MessageLogger
#                art_root_io::tfile_support
#                art::Framework_Services_Registry
#             )

#cet_build_plugin(ProduceCompressed art::module LIBRARIES ${MODULE_LIBRARIES})

#simple_plugin(CheckRawDigits "plugin"
#                ${MODULE_LIBRARIES})

install_headers()
install_source()
install_fhicl()
