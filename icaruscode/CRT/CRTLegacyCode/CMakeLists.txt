cet_build_plugin(CRTTzeroProducer art::module
     LIBRARIES
        larcorealg::Geometry
        sbnobj::ICARUS_CRT
        sbnobj::Common_CRT
        art_root_io::TFileService_service
        lardataalg::DetectorInfo
        nurandom::RandomUtils_NuRandomService_service
        art::Framework_Services_Registry
        art::Framework_Services_Optional_RandomNumberGenerator_service
        messagefacility::MF_MessageLogger
        messagefacility::headers
        CLHEP::CLHEP
        lardata::Utilities
        )

cet_build_plugin(CRTT0Matching art::module
     LIBRARIES
          sbnobj::Common_CRT
          icaruscode_CRTUtils
          larcore::Geometry_Geometry_service
          larsim::Simulation lardataobj::Simulation
          larsim::MCCheater_BackTrackerService_service
          larsim::MCCheater_ParticleInventoryService_service
          lardata::Utilities
          larevt::Filters
          lardataobj::RawData
          lardataobj::RecoBase
          lardataobj::AnalysisBase
          lardata::RecoObjects
          larpandora::LArPandoraInterface
          larcorealg::Geometry
          nusimdata::SimulationBase
          art::Framework_Core
          art::Framework_Principal
          art::Framework_Services_Registry
          art_root_io::tfile_support
          art_root_io::TFileService_service
          art::Persistency_Common canvas::canvas
          art::Persistency_Provenance
          art::Utilities
          messagefacility::MF_MessageLogger
          ROOT::Core
          ROOT::Geom
          ROOT::XMLIO
          ROOT::Gdml
          ROOT::Tree
          ROOT::Spectrum
          ROOT::RooFit
          ROOT::RooFitCore
)
   
   
cet_build_plugin(CRTT0MatchingAna art::module
     LIBRARIES
          sbnobj::ICARUS_CRT
          sbnobj::Common_CRT
          icaruscode_CRTUtils
          larcorealg::Geometry
          larcore::Geometry_Geometry_service
          larsim::Simulation lardataobj::Simulation
          larsim::MCCheater_BackTrackerService_service
          larsim::MCCheater_ParticleInventoryService_service
          lardata::Utilities
          larevt::Filters
          lardataobj::RawData
          lardataobj::RecoBase
          lardataobj::AnalysisBase
          lardata::RecoObjects
          larpandora::LArPandoraInterface
          larcorealg::Geometry
          nusimdata::SimulationBase
          art::Persistency_Common canvas::canvas
          art::Persistency_Provenance
          art::Utilities
          ROOT::Core
          ROOT::Tree
          ROOT::Geom
          ROOT::XMLIO
          ROOT::Gdml
          ROOT::Spectrum
          ROOT::RooFit
          ROOT::RooFitCore
)


install_headers()
install_fhicl()
install_source()
