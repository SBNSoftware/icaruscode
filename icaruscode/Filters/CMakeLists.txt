art_make_library()
set(  MODULE_LIBRARIES
    sbnobj::Common_Trigger
    icarusalg::Utilities
    sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_ICARUS
    artdaq_core::artdaq-core_Utilities
    larcorealg::Geometry
    larcore::Geometry_Geometry_service
    larsim::Simulation
    nug4::ParticleNavigation
    icaruscode::Decode_DataProducts
    lardataobj::Simulation
    lardata::Utilities
    larevt::Filters
    lardataalg::MCDumpers
    lardataobj::RawData
    lardataobj::RecoBase
    lardata::RecoObjects
    nusimdata::SimulationBase
    art::Framework_Core
    art::Framework_Principal
    art::Framework_Services_Registry
    art_root_io::tfile_support
    art_root_io::TFileService_service
    ROOT::Core
    art::Persistency_Common
    art::Persistency_Provenance
    art::Utilities
    messagefacility::MF_MessageLogger
    messagefacility::headers
    fhiclcpp::fhiclcpp
    cetlib::cetlib
    cetlib_except::cetlib_except
    ROOT::EG
    ROOT::Geom
    ROOT::XMLIO
    ROOT::Gdml
    ROOT::GenVector
    ROOT::Tree
  )
cet_build_plugin(FilterNeutrinosActiveVolume art::module LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(FilterNumberTPCHits art::module LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(FilterDataIntegrity art::module LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(FilterOpFlash art::module LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(FilterParticlesActiveVolume art::module LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(TriggerTypeFilter art::module LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(FilterCRTPMTMatching art::module LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(FilterDirts art::module LIBRARIES ${MODULE_LIBRARIES})

# install_headers()
install_fhicl()
install_source()

