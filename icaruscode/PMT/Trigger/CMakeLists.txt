add_subdirectory(Algorithms)
add_subdirectory(Utilities)
add_subdirectory(scripts)

art_make_library(
	EXCLUDE
    "WriteBeamGateInfo_module.cc"
    "FixBeamGateInfo_module.cc"
    "DiscriminatePMTwaveformsByChannel_module.cc"
    "DiscriminatedAdderSignal_module.cc"
    "MajorityTriggerSimulation_module.cc"
    "SlidingWindowTrigger_module.cc"
    "SlidingWindowTriggerSimulation_module.cc"
    "TriggerSimulationOnGates_module.cc"
    "TriggerEfficiencyPlots_module.cc"
    "MajorityTriggerEfficiencyPlots_module.cc"
    "SlidingWindowTriggerEfficiencyPlots_module.cc"
    "TriggerEmulationTree_module.cc"
    "MakeTriggerSimulationTree_module.cc"
	LIBRARIES
		icaruscode::PMT_Trigger_Algorithms
		icaruscode::PMT_Trigger_Utilities
		sbnobj::ICARUS_PMT_Trigger_Data
		larcore::Geometry_Geometry_service
		lardataalg::MCDumpers
		larcorealg::Geometry
		lardataobj::RawData
		lardata::Utilities
		lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
		nusimdata::SimulationBase
		art_root_io::TFileService_service
		art_root_io::tfile_support
		art::Framework_Services_Registry
		art::Framework_Principal
		art::Framework_Core
		art::Persistency_Provenance
		art::Utilities
		canvas::canvas
		messagefacility::MF_MessageLogger
		fhiclcpp::fhiclcpp
		cetlib_except::cetlib_except
		ROOT::GenVector
		ROOT::Hist
		ROOT::Tree
		ROOT::Core
)
set(	MODULE_LIBRARIES
		icaruscode::PMT_Trigger_Algorithms
		sbnobj::ICARUS_PMT_Trigger_Data
		lardataalg::DetectorInfo
		lardataobj::RawData
		lardata::Utilities
		larcore::Geometry_Geometry_service
		art::Framework_Services_Registry
		messagefacility::MF_MessageLogger
		fhiclcpp::fhiclcpp
		ROOT::Core
		ROOT::Physics
		TOOL_LIBRARIES
		icaruscode::PMT_Trigger_Algorithms
		sbnobj::ICARUS_PMT_Trigger_Data
		lardataalg::DetectorInfo
		lardataobj::RawData
		lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
		messagefacility::MF_MessageLogger
		fhiclcpp::fhiclcpp
		larcore::Geometry_Geometry_service
)
set(	TOOL_LIBRARIES
		icaruscode_PMT_Trigger_Algorithms
		sbnobj::ICARUS_PMT_Trigger_Data
		lardataalg::DetectorInfo
		lardataobj::RawData
	)

cet_build_plugin(DiscriminatePMTwaveforms art::module LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(DumpTriggerGateData art::module LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(DynamicTriggerGateBuilderTool art::tool LIBRARIES ${TOOL_LIBRARIES})
cet_build_plugin(ExtractEnergyDepositionSummary art::module LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(FixedTriggerGateBuilderTool art::tool LIBRARIES ${TOOL_LIBRARIES})
cet_build_plugin(LVDSgates art::module LIBRARIES ${MODULE_LIBRARIES})


cet_build_plugin(WriteBeamGateInfo art::module LIBRARIES
  lardataobj::Simulation
  lardata::Utilities
  larcorealg::headers
  lardataalg::headers
  )

cet_build_plugin(FixBeamGateInfo art::module LIBRARIES
  lardataobj::Simulation
  larcorealg::headers
  lardataalg::headers
  ) 

cet_build_plugin(DiscriminatePMTwaveformsByChannel art::module
                LIBRARIES
  icaruscode::PMT_Trigger_Algorithms
  sbnobj::ICARUS_PMT_Trigger_Data
  sbnobj::ICARUS_PMT_Data
  sbnobj::Common_PMT_Data
  icaruscode::PMT_Trigger_Utilities
  icaruscode::PMT_Algorithms
  lardataobj::RawData
  lardata::Utilities
  lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
  larcore::Geometry_Geometry_service
  art::Framework_Services_Registry
  art::Framework_Principal
  art::Framework_Core
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib_except::cetlib_except
  ROOT::Core
  )

cet_build_plugin(DiscriminatedAdderSignal art::module
  LIBRARIES
    icaruscode::PMT_Trigger_Algorithms
    icaruscode::PMT_Trigger_Utilities
    icaruscode::PMT_Algorithms
    icaruscode::IcarusObj
    sbnobj::ICARUS_PMT_Trigger_Data
    sbnobj::ICARUS_PMT_Data
    lardata::headers
    larcore::headers
    lardataalg::DetectorInfo
    lardataalg::headers
    larcorealg::Geometry
    larcorealg::CoreUtils
    lardataobj::RawData
    art::Framework_Services_Registry
  )

cet_build_plugin(SlidingWindowTrigger art::module
                LIBRARIES
		icaruscode::PMT_Trigger_Algorithms
		sbnobj::ICARUS_PMT_Trigger_Data
		icaruscode::PMT_Trigger_Utilities
		icaruscode::PMT_Trigger_Utilities
		icaruscode::PMT_Algorithms
		larcore::Geometry_Geometry_service
		larcorealg::Geometry
		lardata::Utilities
		lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
		lardataobj::RawData
		art::Framework_Services_Registry
		art::Framework_Principal
		art::Framework_Core
		canvas::canvas
		messagefacility::MF_MessageLogger
		fhiclcpp::fhiclcpp
		cetlib_except::cetlib_except
		ROOT::Core
		) 

cet_build_plugin(MajorityTriggerSimulation art::module
                LIBRARIES
		icaruscode::PMT_Trigger_Algorithms
		sbnobj::ICARUS_PMT_Trigger_Data
		icaruscode::PMT_Trigger_Utilities
		larcore::Geometry_Geometry_service
		larcorealg::Geometry
		lardata::Utilities
		lardataobj::RawData
		lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
		art_root_io::TFileService_service
		art_root_io::tfile_support
		art::Framework_Services_Registry
		ROOT::Hist
		ROOT::Core
		) 

cet_build_plugin(SlidingWindowTriggerSimulation art::module
                LIBRARIES
		icaruscode::PMT_Trigger_Algorithms
		sbnobj::ICARUS_PMT_Trigger_Data
		icaruscode::PMT_Trigger_Utilities
		larcore::Geometry_Geometry_service
		larcorealg::Geometry
		lardataobj::RawData
		lardata::Utilities
		lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
		art_root_io::TFileService_service
		art_root_io::tfile_support
		art::Framework_Services_Registry
		art::Framework_Principal
		art::Framework_Core
		canvas::canvas
		messagefacility::MF_MessageLogger
		fhiclcpp::fhiclcpp
		cetlib_except::cetlib_except
		ROOT::Hist
		ROOT::Core
		) 

cet_build_plugin(TriggerSimulationOnGates art::module LIBRARIES
  icaruscode::PMT_Trigger_Algorithms
  icarusalg::Utilities
  sbnobj::ICARUS_PMT_Trigger_Data
  icaruscode::PMT_Trigger_Utilities
  larcore::Geometry_Geometry_service
  larcorealg::Geometry
  lardataobj::Simulation
  lardataobj::RawData
  art_root_io::TFileService_service
  art_root_io::tfile_support
  art::Framework_Services_Registry
  art::Framework_Principal
  art::Framework_Core
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib_except::cetlib_except
  ROOT::Hist
  ROOT::Core
  ) 

cet_build_plugin(TriggerEfficiencyPlots art::module
                LIBRARIES
		icarusalg::Utilities
		icaruscode::PMT_Trigger_Utilities
		sbnobj::ICARUS_PMT_Trigger_Data
		larcore::Geometry_Geometry_service
		lardataalg::MCDumpers
		larcorealg::Geometry
		lardataobj::RawData
		lardata::Utilities
		lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
		nusimdata::SimulationBase
		art_root_io::TFileService_service
		art_root_io::tfile_support
		art::Framework_Services_Registry
		art::Framework_Principal
		art::Framework_Core
		canvas::canvas
		messagefacility::MF_MessageLogger
		fhiclcpp::fhiclcpp
		cetlib_except::cetlib_except
		ROOT::GenVector
		ROOT::Hist
		ROOT::Tree
		ROOT::Core
		) 

cet_build_plugin(MakeTriggerSimulationTree art::module
                LIBRARIES
		icaruscode::PMT_Trigger_Algorithms
		icaruscode::PMT_Trigger_Utilities
		sbnobj::ICARUS_PMT_Trigger_Data
		larcore::Geometry_Geometry_service
		lardataalg::MCDumpers
		lardata::Utilities
		lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
		larcorealg::Geometry
		nusimdata::SimulationBase
		art_root_io::TFileService_service
		art_root_io::tfile_support
		art::Framework_Services_Registry
		art::Framework_Principal
		art::Framework_Core
		canvas::canvas
		messagefacility::MF_MessageLogger
		fhiclcpp::fhiclcpp
		cetlib_except::cetlib_except
		ROOT::GenVector
		ROOT::Tree
		ROOT::Core
		) 

cet_build_plugin(TriggerEmulationTree art::module LIBRARIES
  icaruscode_PMT_Trigger_Algorithms
  icaruscode_PMT_Trigger_Utilities
  art_root_io::TFileService_service
  art_root_io::tfile_support
  art::Framework_Services_Registry
  art::Framework_Principal
  art::Framework_Core
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib_except::cetlib_except
  ROOT::Tree
  ROOT::Core
  ) 

# trigger efficiency plot modules with the same build configuration:
foreach(TriggerEfficiencyPlotsModule
	MajorityTriggerEfficiencyPlots
	SlidingWindowTriggerEfficiencyPlots
	)

	cet_build_plugin(${TriggerEfficiencyPlotsModule} art::module
                LIBRARIES
		icaruscode::PMT_Trigger_Algorithms
		icaruscode::PMT_Trigger
		icaruscode::PMT_Trigger_Utilities
		sbnobj::ICARUS_PMT_Trigger_Data
		larcore::Geometry_Geometry_service
		lardataalg::MCDumpers
		lardata::Utilities
		lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
		larcorealg::Geometry
		lardataobj::RawData
		nusimdata::SimulationBase
		art_root_io::TFileService_service
		art_root_io::tfile_support
		art::Framework_Services_Registry
		art::Framework_Principal
		art::Framework_Core
		art::Framework_Principal
		canvas::canvas
		messagefacility::MF_MessageLogger
		fhiclcpp::fhiclcpp
		cetlib_except::cetlib_except
		ROOT::GenVector
		ROOT::Hist
		ROOT::Tree
		ROOT::Core
		) 

endforeach()

install_headers()
install_source()
install_fhicl()
