
cet_make_library(
  SOURCE
  ICARUSPixelMapProducer.cxx
  ICARUSICVNMapper.cxx
  ICARUSICVNZlibMaker.cxx
  ICARUSITFNetHandler.h
  LIBRARY_NAME
  icaruscode_icaruscvn_module_helpers
  LIBRARIES
  PUBLIC
  larrecodnn::CVN_interfaces
  larrecodnn::CVN_func
  lardataobj::RecoBase
  lardataobj::Simulation
  art::Framework_Services_Registry
  art::Framework_Core
  art::Framework_Principal
  art_root_io::tfile_support
  art_root_io::TFileService_service
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  art_plugin_types::module
  Boost::filesystem
  art_root_io::tfile_support
  nusimdata::SimulationBase
  ROOT::Hist
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  larcorealg::Geometry
  
  icaruscode::ICARUSCVN_tf
  ${ZLIB_LIBRARIES}
  PRIVATE
  larsim::MCCheater_BackTrackerService_service
  larsim::MCCheater_ParticleInventoryService_service
)

art_dictionary(DICTIONARY_LIBRARIES icaruscode_icaruscvn_module_helpers)

install_headers()
install_source()
