cet_make_library( SOURCE
                  ICVNMapperICARUS.cxx
                  ICVNZlibMakerICARUS.cxx
                  PixelMapProducerICARUS.cxx
                  LIBRARIES 
                  larcorealg::Geometry
                  larcore::Geometry_Geometry_service
                  larrecodnn::CVN_interfaces
                  larsim::Simulation
                  lardataobj::Simulation
                  larsim::MCCheater_BackTrackerService_service
                  larsim::MCCheater_ParticleInventoryService_service
                  lardata::Utilities
                  larevt::Filters
                  lardataobj::RawData
                  lardataobj::RecoBase
                  lardata::RecoObjects  
                  lardataalg::DetectorInfo
                  larpandora::LArPandoraInterface
                  nusimdata::SimulationBase
                  cetlib::cetlib
                  art::Framework_Core
                  art::Framework_Principal
                  art::Framework_Services_Registry
                  art_root_io::tfile_support ROOT::Core
                  art_root_io::TFileService_service
                  art::Persistency_Common
                  art::Persistency_Provenance
                  art::Utilities
                  canvas::canvas
                  messagefacility::MF_MessageLogger
                  
                  fhiclcpp::fhiclcpp
                  ROOT::Geom
                  ROOT::XMLIO
                  ROOT::Gdml
                  ${ROOT_BASIC_LIB_LIST}
)

install_headers()
install_fhicl()
install_source()