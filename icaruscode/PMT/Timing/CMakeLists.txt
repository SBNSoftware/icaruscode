cet_enable_asserts()

simple_plugin(TimingCorrectionExtraction "module"
  lardataobj_RawData
  ${ART_FRAMEWORK_CORE}
  ${ART_FRAMEWORK_PRINCIPAL}
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  art_Persistency_Provenance
  ${CANVAS}
  ${MF_MESSAGELOGGER}
  ${FHICLCPP}
  ${CETLIB}
  cetlib_except
  ${ROOT_BASIC_LIB_LIST}
  LIBWDA
  ${Boost_SYSTEM_LIBRARY}
  ${ART_UTILITIES}
  art_Persistency_Common canvas
  art_Persistency_Provenance canvas
  art_Utilities canvas
  ${MF_MESSAGELOGGER}
  ${FHICLCPP}
  cetlib cetlib_except
)

install_headers()
install_fhicl()
install_source()