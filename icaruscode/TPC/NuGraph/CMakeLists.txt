set(CMAKE_EXPORT_COMPILE_COMMANDS ON) # for clangd autocompletion
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

cet_build_plugin(ICARUSNuGraphMultiLoader art::tool
  LIBRARIES PRIVATE
  ${nugraph_tool_lib_list}
)

cet_build_plugin(ICARUSNuGraphInference art::EDProducer
    LIBRARIES PRIVATE
    ${nugraph_tool_lib_list}
    IMPL_TARGET_VAR ICARUSNuGraphInference_module
)

find_path(DELAUNATOR_INC delaunator-header-only.hpp HINTS ENV DELAUNATOR_INC REQUIRED DOC "Header location for delaunator_cpp")
target_include_directories(${ICARUSNuGraphInference_module} PRIVATE ${DELAUNATOR_INC})

simple_plugin(ICARUSNCCSlicesProducer "module"
    LIBRARIES PRIVATE
    art::Framework_Core
    lardataobj::RecoBase
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
