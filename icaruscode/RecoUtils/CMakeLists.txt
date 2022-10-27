art_make_library(
     LIBRARIES larcorealg::Geometry
          larcore::Geometry_Geometry_service
          larsim::Simulation lardataobj::Simulation
          larsim::MCCheater_BackTrackerService_service
          larsim::MCCheater_ParticleInventoryService_service
          lardata::Utilities
          larevt::Filters
          lardataobj::RawData
          lardataobj::RecoBase
          lardata::RecoObjects
          larpandora::LArPandoraInterface
          nusimdata::SimulationBase
          cetlib::cetlib cetlib_except::cetlib_except
          art::Framework_Core
          art::Framework_Principal
          art::Framework_Services_Registry
          art_root_io::tfile_support ROOT::Core
          art_root_io::TFileService_service
          art::Persistency_Common canvas::canvas
          art::Persistency_Provenance
          art::Utilities
          messagefacility::MF_MessageLogger
          messagefacility::headers
          fhiclcpp::fhiclcpp
          ROOT::Geom
          ROOT::XMLIO
          ROOT::Gdml
          Boost::system
)

install_headers()
install_fhicl()
install_source()
