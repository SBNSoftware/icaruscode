add_subdirectory(Objects)
add_subdirectory(ArtObjects)
add_subdirectory(scripts)

file(GLOB lib_srcs *.cxx details/*.cxx)

art_make_library(
  SOURCE ${lib_srcs}
  LIBRARIES
    larcorealg::Geometry
    lardataobj::RawData
    art::Framework_Principal
    art::Framework_Core
    canvas::canvas
    ROOT::Tree
  )

cet_build_plugin(TimedTrackSelector art::module
  LIBRARIES
  icarusalg::Utilities
  lardataobj::AnalysisBase
  lardataobj::RecoBase
  messagefacility::MF_MessageLogger
  ROOT::Core
  )

cet_build_plugin(TimedTrackSelectorCRT art::module
  LIBRARIES
  icarusalg::Utilities
  lardataobj::AnalysisBase
  lardataobj::RecoBase
  messagefacility::MF_MessageLogger
  ROOT::Core
  )

cet_build_plugin(BeamGateInfoFromTracks  art::module
  LIBRARIES
  lardataobj::Simulation
  lardataobj::AnalysisBase
  lardataobj::RecoBase
  lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
  )

cet_build_plugin(BeamGateInfoFromTracksCRT  art::module
  LIBRARIES
  lardataobj::Simulation
  lardataobj::AnalysisBase
  lardataobj::RecoBase
  lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
  )

cet_build_plugin(TimeTrackTreeStorage art::module
  LIBRARIES
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

cet_build_plugin(TimeTrackTreeStorageCRT art::module
  LIBRARIES
  icaruscode::Analysis_trigger
  icaruscode_CRTUtils
  icaruscode::PMT_Algorithms
  icaruscode::Decode_DataProducts
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
