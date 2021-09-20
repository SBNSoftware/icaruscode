add_subdirectory(Objects)
add_subdirectory(scripts)

art_make(
  EXCLUDE
    "TimedTrackSelector_module.cc"
  )

simple_plugin(TimedTrackSelector module
  lardataobj_AnalysisBase
  lardataobj_RecoBase
  ${MF_MESSAGELOGGER}
  ROOT::Core
  )

install_headers()
install_source()
install_fhicl()
