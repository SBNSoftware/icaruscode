set(MODULE_LIBRARIES
   larcorealg::Geometry
   larcore::Geometry_Geometry_service
   larsim::Simulation lardataobj_Simulation
   larsim::MCCheater_BackTrackerService_service
   larsim::MCCheater_ParticleInventoryService_service
   lardata::Utilities
   larevt::Filters
   lardataobj::RawData
   lardataobj::RecoBase
   larreco::Calorimetry
   larreco::RecoAlg
   larrecodnn::CVN_func
   larrecodnn::CVN_interfaces
   lardata::RecoObjects  
   larpandora::LArPandoraInterface
   icaruscode::CRT
   icaruscode::RecoUtils
   icaruscode::TPC_CVN_interfaces
   sbnobj::Common_CRT
   nusimdata::SimulationBase
   art::Framework_Core
   art::Framework_Principal
   art::Framework_Services_Registry
   art_root_io::tfile_support ROOT::Core
   art_root_io::TFileService_service
   art::Persistency_Common canvas
   art::Persistency_Provenance canvas
   art::Utilities canvas
   messagefacility::MF_MessageLogger

   fhiclcpp::fhiclcpp
   ROOT::Geom
   ROOT::XMLIO
   ROOT::Gdml
   ${ROOT_BASIC_LIB_LIST}
)

cet_build_plugin(CVNMapperICARUS art::module
   SOURCE CVNMapperICARUS_module.cc
   LIBRARIES ${MODULE_LIBRARIES}
)

cet_build_plugin(CVNZlibMakerICARUS art::module
   SOURCE CVNZlibMakerICARUS_module.cc
   LIBRARIES ${MODULE_LIBRARIES}
)

install_headers()
install_fhicl()
install_source()
