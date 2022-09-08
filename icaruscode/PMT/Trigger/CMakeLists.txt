add_subdirectory(Algorithms)
add_subdirectory(Utilities)
add_subdirectory(scripts)

art_make(
  EXCLUDE
    "WriteBeamGateInfo_module.cc"
    "FixBeamGateInfo_module.cc"
    "DiscriminatePMTwaveformsByChannel_module.cc"
    "MajorityTriggerSimulation_module.cc"
    "SlidingWindowTrigger_module.cc"
    "SlidingWindowTriggerSimulation_module.cc"
    "TriggerSimulationOnGates_module.cc"
    "TriggerEfficiencyPlots_module.cc"
    "MajorityTriggerEfficiencyPlots_module.cc"
    "SlidingWindowTriggerEfficiencyPlots_module.cc"
    "TriggerEmulationTree_module.cc"
    "MakeTriggerSimulationTree_module.cc"
  LIB_LIBRARIES
    icaruscode_PMT_Trigger_Algorithms
    icaruscode_PMT_Trigger_Utilities
    icarusalg_Utilities
    sbnobj_ICARUS_PMT_Trigger_Data
    larcore_Geometry_Geometry_service
    lardataalg_MCDumpers
    larcorealg_Geometry
    lardataobj_RawData
    nusimdata_SimulationBase
    ${ART_ROOT_IO_TFILESERVICE_SERVICE}
    ${ART_ROOT_IO_TFILE_SUPPORT}
    ${ART_FRAMEWORK_SERVICES_REGISTRY}
    ${ART_FRAMEWORK_PRINCIPAL}
    ${ART_FRAMEWORK_CORE}
    ${ART_PERSISTENCY_PROVENANCE}
    ${ART_UTILITIES}
    ${CANVAS}
    ${MF_MESSAGELOGGER}
    ${FHICLCPP}
    ${CETLIB_EXCEPT}
    ROOT::GenVector
    ROOT::Hist
    ROOT::Tree
    ROOT::Core
  MODULE_LIBRARIES
    icaruscode_PMT_Trigger_Algorithms
    sbnobj_ICARUS_PMT_Trigger_Data
    lardataalg_DetectorInfo
    lardataobj_RawData
    ${ART_FRAMEWORK_SERVICES_REGISTRY}
    ${MF_MESSAGELOGGER}
    ${FHICLCPP}
    ROOT::Core
    ROOT::Physics
  TOOL_LIBRARIES
    icaruscode_PMT_Trigger_Algorithms
    sbnobj_ICARUS_PMT_Trigger_Data
    lardataalg_DetectorInfo
    lardataobj_RawData
    ${MF_MESSAGELOGGER}
    ${FHICLCPP}
  )

simple_plugin(WriteBeamGateInfo module
  lardataobj_Simulation
  )

simple_plugin(FixBeamGateInfo module
  lardataobj_Simulation
  )

simple_plugin(DiscriminatePMTwaveformsByChannel module
  icaruscode_PMT_Trigger_Algorithms
  sbnobj_ICARUS_PMT_Trigger_Data
  sbnobj_ICARUS_PMT_Data
  sbnobj_Common_PMT_Data
  icaruscode_PMT_Trigger_Utilities
  icaruscode_PMT_Algorithms
  lardataobj_RawData
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ${ART_FRAMEWORK_PRINCIPAL}
  ${ART_FRAMEWORK_CORE}
  ${CANVAS}
  ${MF_MESSAGELOGGER}
  ${FHICLCPP}
  ${CETLIB_EXCEPT}
  ROOT::Core
  )

simple_plugin(SlidingWindowTrigger module
  icaruscode_PMT_Trigger_Algorithms
  sbnobj_ICARUS_PMT_Trigger_Data
  icaruscode_PMT_Trigger_Utilities
  icaruscode_PMT_Algorithms
  larcore_Geometry_Geometry_service
  larcorealg_Geometry
  lardataobj_RawData
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ${ART_FRAMEWORK_PRINCIPAL}
  ${ART_FRAMEWORK_CORE}
  ${CANVAS}
  ${MF_MESSAGELOGGER}
  ${FHICLCPP}
  ${CETLIB_EXCEPT}
  ROOT::Core
  )

simple_plugin(MajorityTriggerSimulation module
  icaruscode_PMT_Trigger_Algorithms
  icarusalg_Utilities
  sbnobj_ICARUS_PMT_Trigger_Data
  icaruscode_PMT_Trigger_Utilities
  larcore_Geometry_Geometry_service
  larcorealg_Geometry
  lardataobj_RawData
  ${ART_ROOT_IO_TFILESERVICE_SERVICE}
  ${ART_ROOT_IO_TFILE_SUPPORT}
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ${ART_FRAMEWORK_PRINCIPAL}
  ${ART_FRAMEWORK_CORE}
  ${CANVAS}
  ${MF_MESSAGELOGGER}
  ${FHICLCPP}
  ${CETLIB_EXCEPT}
  ROOT::Hist
  ROOT::Core
  )

simple_plugin(SlidingWindowTriggerSimulation module
  icaruscode_PMT_Trigger_Algorithms
  icarusalg_Utilities
  sbnobj_ICARUS_PMT_Trigger_Data
  icaruscode_PMT_Trigger_Utilities
  larcore_Geometry_Geometry_service
  larcorealg_Geometry
  lardataobj_RawData
  ${ART_ROOT_IO_TFILESERVICE_SERVICE}
  ${ART_ROOT_IO_TFILE_SUPPORT}
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ${ART_FRAMEWORK_PRINCIPAL}
  ${ART_FRAMEWORK_CORE}
  ${CANVAS}
  ${MF_MESSAGELOGGER}
  ${FHICLCPP}
  ${CETLIB_EXCEPT}
  ROOT::Hist
  ROOT::Core
  )

simple_plugin(TriggerSimulationOnGates module
  icaruscode_PMT_Trigger_Algorithms
  icarusalg_Utilities
  sbnobj_ICARUS_PMT_Trigger_Data
  icaruscode_PMT_Trigger_Utilities
  larcore_Geometry_Geometry_service
  larcorealg_Geometry
  lardataobj_Simulation
  lardataobj_RawData
  art_root_io::TFileService_service
  art_root_io::tfile_support
  art::Framework_Services_Registry
  art::Framework_Principal
  art::Framework_Core
  canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib_except
  ROOT::Hist
  ROOT::Core
  )

simple_plugin(TriggerEfficiencyPlots module
  icaruscode_PMT_Trigger_Utilities
  sbnobj_ICARUS_PMT_Trigger_Data
  larcore_Geometry_Geometry_service
  lardataalg_MCDumpers
  larcorealg_Geometry
  lardataobj_RawData
  nusimdata_SimulationBase
  ${ART_ROOT_IO_TFILESERVICE_SERVICE}
  ${ART_ROOT_IO_TFILE_SUPPORT}
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ${ART_FRAMEWORK_PRINCIPAL}
  ${ART_FRAMEWORK_CORE}
  ${CANVAS}
  ${MF_MESSAGELOGGER}
  ${FHICLCPP}
  ${CETLIB_EXCEPT}
  ROOT::GenVector
  ROOT::Hist
  ROOT::Tree
  ROOT::Core
  )

simple_plugin(MakeTriggerSimulationTree module
  icaruscode_PMT_Trigger_Algorithms
  icaruscode_PMT_Trigger_Utilities
  sbnobj_ICARUS_PMT_Trigger_Data
  larcore_Geometry_Geometry_service
  lardataalg_MCDumpers
  larcorealg_Geometry
  nusimdata_SimulationBase
  ${ART_ROOT_IO_TFILESERVICE_SERVICE}
  ${ART_ROOT_IO_TFILE_SUPPORT}
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ${ART_FRAMEWORK_PRINCIPAL}
  ${ART_FRAMEWORK_CORE}
  ${CANVAS}
  ${MF_MESSAGELOGGER}
  ${FHICLCPP}
  ${CETLIB_EXCEPT}
  ROOT::GenVector
  ROOT::Tree
  ROOT::Core
  )

simple_plugin(TriggerEmulationTree module
  icaruscode_PMT_Trigger_Algorithms
  icaruscode_PMT_Trigger_Utilities
  art_root_io::TFileService_service
  art_root_io::tfile_support
  art::Framework_Services_Registry
  art::Framework_Principal
  art::Framework_Core
  canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib_except
  ROOT::Tree
  ROOT::Core
  )

# trigger efficiency plot modules with the same build configuration:
foreach(TriggerEfficiencyPlotsModule
  MajorityTriggerEfficiencyPlots
  SlidingWindowTriggerEfficiencyPlots
)

  simple_plugin(${TriggerEfficiencyPlotsModule} module
    icaruscode_PMT_Trigger_Algorithms
    icaruscode_PMT_Trigger
    icaruscode_PMT_Trigger_Utilities
    icarusalg_Utilities
    sbnobj_ICARUS_PMT_Trigger_Data
    larcore_Geometry_Geometry_service
    lardataalg_MCDumpers
    larcorealg_Geometry
    lardataobj_RawData
    nusimdata_SimulationBase
    ${ART_ROOT_IO_TFILESERVICE_SERVICE}
    ${ART_ROOT_IO_TFILE_SUPPORT}
    ${ART_FRAMEWORK_SERVICES_REGISTRY}
    ${ART_FRAMEWORK_PRINCIPAL}
    ${ART_FRAMEWORK_CORE}
    ${ART_FRAMEWORK_PRINCIPAL}
    ${CANVAS}
    ${MF_MESSAGELOGGER}
    ${FHICLCPP}
    ${CETLIB_EXCEPT}
    ROOT::GenVector
    ROOT::Hist
    ROOT::Tree
    ROOT::Core
    )

endforeach()

install_headers()
install_source()
install_fhicl()
