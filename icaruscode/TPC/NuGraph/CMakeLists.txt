cet_enable_asserts()

set( nugraph_tool_lib_list
       lardataobj::RecoBase
       larrecodnn::NuGraphBaseTools
       TorchScatter::TorchScatter
       art::Framework_Core
       hep_concurrency::hep_concurrency
       art_plugin_types::tool
)

include_directories($ENV{HEP_HPC_INC})

cet_build_plugin(ICARUSNuGraphLoader art::tool
  LIBRARIES PRIVATE
  ${nugraph_tool_lib_list}
)

simple_plugin(ICARUSFilteredNuSliceHitsProducer "module"
    LIBRARIES PRIVATE
    art::Framework_Core
    lardataobj::RecoBase
    sbnobj::Common_Reco
)

simple_plugin(ICARUSNuSliceHitsProducer "module"
    LIBRARIES PRIVATE
    art::Framework_Core
    lardataobj::RecoBase
    sbnobj::Common_Reco
)

simple_plugin(ICARUSTrueNuSliceHitsProducer "module"
    LIBRARIES PRIVATE
    art::Framework_Core
    nusimdata::SimulationBase
    lardataobj::RecoBase
    lardataalg::DetectorInfo
    sbnobj::Common_Reco
    larsim::MCCheater_BackTrackerService_service           
    larsim::MCCheater_ParticleInventoryService_service
    lardata::DetectorClocksService
)

simple_plugin(ICARUSHDF5Maker "module"
    LIBRARIES PRIVATE
    art::Framework_Core
    nusimdata::SimulationBase
    larsim::MCCheater_BackTrackerService_service           
    larsim::MCCheater_ParticleInventoryService_service
    hep_hpc_hdf5
    lardataobj::RecoBase
    lardataalg::DetectorInfo
    lardata::DetectorClocksService
    larcorealg::Geometry
    larcore::WireReadout
    ${HDF5_LIBRARIES}
)

cet_build_plugin(ICARUSNuGraphAnalyzer art::EDAnalyzer
    LIBRARIES PRIVATE
    art_root_io::TFileService_service
    ROOT::Tree
    lardata::RecoBaseProxy
)

add_subdirectory(scripts)

install_headers()
install_fhicl()
install_source()

