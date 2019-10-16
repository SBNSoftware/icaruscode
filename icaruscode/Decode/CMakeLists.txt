#find_ups_product(artdaq_core)
#find_ups_product(sbndaq_artdaq_core)

simple_plugin( ICARUSTPCDecoder module
#  tpcAnalysis_MODE
  sbndaq-artdaq-core_Overlays_ICARUS 
  artdaq-core_Utilities
  ${LARDATAOBJ}
  lardataobj_RawData
  ${ART_UTILITIES}
  ${FHICLCPP}
  ${MF_MESSAGELOGGER}
  ${MF_UTILITIES}
  ${ART_FRAMEWORK_CORE}
  ${ROOT_BASIC_LIB_LIST}
                           ${ART_FRAMEWORK_SERVICES_REGISTRY}
                           ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                           ${ART_ROOT_IO_TFILESERVICE_SERVICE}
)
install_headers()
install_fhicl()
install_source()
