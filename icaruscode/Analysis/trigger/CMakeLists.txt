link_directories( $ENV{ROOTSYS}/lib)
add_subdirectory(Objects)
add_subdirectory(scripts)

art_make(
  EXCLUDE
    "TimedTrackSelector_module.cc"
    "TimeTrackTreeStorage_module.cc"
  )

simple_plugin(TimedTrackSelector module
  lardataobj_AnalysisBase
  lardataobj_RecoBase
  ${MF_MESSAGELOGGER}
  ROOT::Core
  )

simple_plugin(TimeTrackTreeStorage module				   
  lardataobj_AnalysisBase  
  lardataobj_RecoBase
  lardataobj_Simulation
  lardataobj_RawData
  ${MF_MESSAGELOGGER}
  ROOT::Core
  ${ROOT_EVE_LIB_LIST}
  ${ROOT_X3d}
  ${ROOT_BASIC_LIB_LIST}
  art_root_io_tfile_support
  art_root_io
  art_root_io_dict
  ${ART_FRAMEWORK_CORE}
  ${ART_FRAMEWORK_SERVICES}
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ${ART_ROOT_IO_TFILESERVICE_SERVICE}
  ${ART_FRAMEWORK_PRINCIPAL}
  art_Persistency_Common
  art_Utilities canvas
)
		
install_headers()
install_source()
install_fhicl()
