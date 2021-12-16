#link_directories( $ENV{ROOTSYS}/lib)
add_subdirectory(Objects)
add_subdirectory(scripts)

art_make(
  SUBDIRS details
  EXCLUDE
    "TimedTrackSelector_module.cc"
    "BeamGateInfoFromTracks_module.cc"
    "TimeTrackTreeStorage_module.cc"
  LIB_LIBRARIES
    lardataobj_RawData
    art::Framework_Principal
    art::Framework_Core
    canvas
    messagefacility::MF_MessageLogger
    ROOT::Tree
  )

simple_plugin(TimedTrackSelector module
  lardataobj_AnalysisBase
  lardataobj_RecoBase
  ${MF_MESSAGELOGGER}
  ROOT::Core
  )

simple_plugin(BeamGateInfoFromTracks  module
  lardataobj_Simulation
  lardataobj_AnalysisBase
  lardataobj_RecoBase
  )

simple_plugin(TimeTrackTreeStorage module
  icaruscode_Analysis_trigger
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
