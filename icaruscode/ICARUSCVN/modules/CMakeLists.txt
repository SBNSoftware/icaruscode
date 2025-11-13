cet_build_plugin(ICARUSCVNMapper art::EDProducer
  LIBRARIES PRIVATE
  larcorealg::Geometry
  larcore::Geometry_Geometry_service
  larsim::Simulation lardataobj::Simulation
  larsim::MCCheater_BackTrackerService_service
  larsim::MCCheater_ParticleInventoryService_service
  lardata::Utilities
  larevt::Filters
  lardataobj::RawData
  lardataobj::RecoBase
  larreco::Calorimetry
  larreco::RecoAlg
  lardata::RecoObjects  
  larpandora::LArPandoraInterface
  sbnobj::Common_CRT
  nusimdata::SimulationBase
  art::Framework_Core
  art::Framework_Principal
  art::Framework_Services_Registry
  art_root_io::tfile_support ROOT::Core
  art_root_io::TFileService_service
  art::Persistency_Common canvas::canvas
  art::Persistency_Provenance canvas::canvas
  art::Utilities canvas::canvas
  messagefacility::MF_MessageLogger
                           
  fhiclcpp::fhiclcpp
  ROOT::Geom
  ROOT::XMLIO
  ROOT::Gdml
  ${ROOT_BASIC_LIB_LIST}
  icaruscode_icaruscvn_module_helpers
)

cet_build_plugin(ICARUSCVNZlibMaker art::EDAnalyzer
  LIBRARIES PRIVATE
  larcorealg::Geometry
  larcore::Geometry_Geometry_service
  larsim::Simulation lardataobj::Simulation
  larsim::MCCheater_BackTrackerService_service
  larsim::MCCheater_ParticleInventoryService_service
  lardata::Utilities
  larevt::Filters
  lardataobj::RawData
  lardataobj::RecoBase
  larreco::Calorimetry
  larreco::RecoAlg
  lardata::RecoObjects  
  larpandora::LArPandoraInterface
  sbnobj::Common_CRT
  nusimdata::SimulationBase
  art::Framework_Core
  art::Framework_Principal
  art::Framework_Services_Registry
  art_root_io::tfile_support ROOT::Core
  art_root_io::TFileService_service
  art::Persistency_Common canvas::canvas
  art::Persistency_Provenance canvas::canvas
  art::Utilities canvas::canvas
  messagefacility::MF_MessageLogger
                           
  fhiclcpp::fhiclcpp
  ROOT::Geom
  ROOT::XMLIO
  ROOT::Gdml
  ${ROOT_BASIC_LIB_LIST}
  icaruscode_icaruscvn_module_helpers
)

cet_build_plugin(ICARUSTFNetHandler art::tool
  LIBRARIES PRIVATE
  larrecodnn::CVN_func
  larrecodnn::CVN_interfaces
  larrecodnn::ImagePatternAlgs_Tensorflow_TF
  icaruscode::ICARUSCVN_tf
  messagefacility::MF_MessageLogger
  fhiclcpp::types
  fhiclcpp::fhiclcpp
  ROOT::Hist
  icaruscode_icaruscvn_module_helpers
  icaruscode::ICARUSCVN_tf
)

cet_build_plugin(ICARUSCVNEvaluator art::EDProducer
  LIBRARIES PRIVATE
  larrecodnn::CVN_func
  larrecodnn::CVN_interfaces
  art_plugin_support::toolMaker
  art::Framework_Core
  art::Framework_Principal
  art_root_io::TFileService_service
  canvas::canvas
  icaruscode_icaruscvn_module_helpers
)

install_headers()
install_fhicl()
install_source()

