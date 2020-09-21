cet_enable_asserts()

add_definitions(-DEIGEN_FFTW_DEFAULT)

art_make(
          LIB_LIBRARIES
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        lardata_Utilities
                        larevt_Filters
                        lardataobj_RawData
                        lardata_ArtDataHelper
                        ${ROOT_BASIC_LIB_LIST}
                        ${ROOT_SPECTRUM}
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
                        ${ICARUS_FFTW_LIBRARIES}
        )

install_headers()
install_fhicl()
install_source()
