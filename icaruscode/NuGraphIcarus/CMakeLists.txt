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

cet_build_plugin(IcarusNuGraphLoader art::tool
  LIBRARIES PRIVATE
  ${nugraph_tool_lib_list}
)

simple_plugin(IcarusNuSliceHitsProducer "module"
    LIBRARIES PRIVATE
    art::Framework_Core
    lardataobj::RecoBase
    sbnobj::Common_Reco
)

simple_plugin(IcarusTrueNuSliceHitsProducer "module"
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

simple_plugin(IcarusHDF5Maker "module"
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

add_subdirectory(scripts)

install_headers()
install_fhicl()
install_source()

