cet_enable_asserts()

art_make(
          EXCLUDE
                        PMTconfigurationExtraction_module.cc
                        DumpPMTconfiguration_module.cc
          MODULE_LIBRARIES
                        sbndaq-artdaq-core_Overlays_ICARUS 
                        artdaq-core_Utilities
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        lardata_Utilities
                        larevt_Filters
                        lardataobj_RawData
                        lardata_ArtDataHelper
                        ${ROOT_BASIC_LIB_LIST}
                        ${LIBWDA}
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
                        ${TBB}
        )

simple_plugin(PMTconfigurationExtraction module
  icaruscode_Decode_DecoderTools
  sbnobj_Common_PMT_Data
  art_Framework_Services_Registry
  fhiclcpp
  cetlib_except
  )

simple_plugin(DumpPMTconfiguration module
  sbnobj_Common_PMT_Data
  ${MF_MESSAGELOGGER}
  fhiclcpp
  cetlib_except
  )

install_headers()
install_fhicl()
install_source()

# Add our tools directory
add_subdirectory(DecoderTools)
add_subdirectory(ChannelMapping)
add_subdirectory(fcl)

