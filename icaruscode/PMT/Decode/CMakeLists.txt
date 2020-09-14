cet_enable_asserts()

art_make(
          MODULE_LIBRARIES
                        sbndaq-artdaq-core_Overlays_ICARUS
                        artdaq-core_Utilities
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        lardata_Utilities
                        larevt_Filters
                        lardataobj_RawData
                        lardata_ArtDataHelper
                        larana_OpticalDetector_OpHitFinder
                        ${ROOT_BASIC_LIB_LIST}
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_ROOT_IO_TFILE_SUPPORT}
                        ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                        ${ART_UTILITIES}
                        art_Persistency_Common canvas
                        art_Persistency_Provenance canvas
                        art_Utilities canvas
                        ${MF_MESSAGELOGGER}
                        ${MF_UTILITIES}
                        ${FHICLCPP}
                        cetlib cetlib_except
        )

install_headers()
install_fhicl()
install_source()

# Add our tools directory
#add_subdirectory(DecoderTools)
