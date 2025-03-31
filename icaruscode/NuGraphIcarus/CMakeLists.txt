cet_enable_asserts()

set( nugraph_tool_lib_list
       lardataobj::RecoBase
       lardataobj::AnalysisBase
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
    nusimdata::SimulationBase
    lardataobj::RecoBase
    lardataalg::DetectorInfo
    lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
    larcore::Geometry_Geometry_service
    sbnobj::Common_Reco
)

simple_plugin(IcarusTrueNuSliceHitsProducer "module"
    LIBRARIES PRIVATE
    art::Framework_Core
    nusimdata::SimulationBase
    lardataobj::RecoBase
    lardataalg::DetectorInfo
    lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
    larcore::Geometry_Geometry_service
    sbnobj::Common_Reco
    larsim::MCCheater_BackTrackerService_service           
    larsim::MCCheater_ParticleInventoryService_service
)

simple_plugin(IcarusHDF5Maker "module"
    LIBRARIES PRIVATE
    art::Framework_Core
    nusimdata::SimulationBase
    larsim::MCCheater_BackTrackerService_service           
    larsim::MCCheater_ParticleInventoryService_service
    hep_hpc_hdf5
    #hep::hpc_Utilities
    lardata::RecoBaseProxy
    lardataobj::RecoBase
    lardataalg::DetectorInfo
    lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
    larcore::Geometry_Geometry_service
    ${HDF5_LIBRARIES}
)

install_headers()
install_fhicl()
install_source()

