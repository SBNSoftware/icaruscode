
art_make(
  EXCLUDE
    "SaveConfigurationIntoTFile_module.cc"
  LIB_LIBRARIES
    ${CETLIB_EXCEPT}
    ${CLHEP}
  )

simple_plugin(SaveConfigurationIntoTFile "module"
              ${ART_ROOT_IO_TFILESERVICE_SERVICE}
              ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
              ${ART_FRAMEWORK_SERVICES_REGISTRY}
              ${ART_FRAMEWORK_PRINCIPAL}
              ${ART_FRAMEWORK_CORE}
              canvas
              ${MF_MESSAGELOGGER}
              hep_concurrency
              ${FHICLCPP}
              ${CETLIB}
              ${ROOT_CORE}
)


install_headers()
install_source()
install_fhicl()
