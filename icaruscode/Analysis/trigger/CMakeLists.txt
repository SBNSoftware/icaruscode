add_subdirectory(Objects)
add_subdirectory(scripts)

art_make(
  SUBDIRS details
  EXCLUDE
    "TimedTrackSelector_module.cc"
    "BeamGateInfoFromTracks_module.cc"
    "TimeTrackTreeStorage_module.cc"
  LIB_LIBRARIES
    larcorealg::Geometry
    lardataobj::RawData
    art::Framework_Principal
    art::Framework_Core
    canvas::canvas
    ROOT::Tree
  )

simple_plugin(TimedTrackSelector module
  lardataobj::AnalysisBase
  lardataobj::RecoBase
  messagefacility::MF_MessageLogger
  ROOT::Core
  )

simple_plugin(BeamGateInfoFromTracks  module
  lardataobj::Simulation
  lardataobj::AnalysisBase
  lardataobj::RecoBase
  lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
  )

simple_plugin(TimeTrackTreeStorage module
  icaruscode_Analysis_trigger
  icaruscode_PMT_Algorithms
  icaruscode_Decode_DataProducts
  larcore::Geometry_Geometry_service
  larcorealg::Geometry
  lardataobj::AnalysisBase
  lardataobj::RecoBase
  lardataobj::RawData
  lardataobj::Simulation
  ROOT::Tree
  ROOT::Core
  art_root_io::TFileService_service
  art::Framework_Services_Registry
  art::Framework_Core
  art::Framework_Principal
  canvas::canvas
  fhiclcpp::fhiclcpp
  messagefacility::MF_MessageLogger
)

install_headers()
install_source()
install_fhicl()
