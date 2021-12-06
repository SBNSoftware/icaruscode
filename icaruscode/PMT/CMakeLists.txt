add_subdirectory(OpticalTools)
add_subdirectory(OpReco)
add_subdirectory(Algorithms)
add_subdirectory(LibraryMappingTools)
add_subdirectory(Trigger)

# Removing AVX instructions because assuming the grid nodes
# being less than 6 year old proved to be pretentious.
# AVX is not so relevant after the last commit anyway
# cet_set_compiler_flags(EXTRA_CXX_FLAGS -mavx)

simple_plugin(SimPMTIcarus "module"
  icaruscode_PMT_Algorithms
  lardataobj_RawData
  lardataobj_Simulation
  nurandom::RandomUtils_NuRandomService_service
  art::Framework_Core
  art::Framework_Principal
  art::Framework_Services_Registry
  art::Framework_Services_Optional_RandomNumberGenerator_service
  art_root_io::TFileService_service
  art_root_io::tfile_support
  art::Persistency_Provenance
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib::cetlib
  cetlib_except::cetlib_except
  ${ROOT_BASIC_LIB_LIST}
  
  )

simple_plugin(PMTWaveformBaselines "module"
  icaruscode_PMT_Algorithms
  sbnobj_ICARUS_PMT_Data
  larcore_Geometry_Geometry_service
  larcorealg_Geometry
  lardataobj_RawData
  art_root_io::TFileService_service
  art_root_io::tfile_support
  art::Framework_Services_Registry
  messagefacility::MF_MessageLogger
  ROOT::Hist
  ROOT::Core
  )

simple_plugin(PMTWaveformBaselinesFromReadoutConfiguration "module"
  sbnobj_ICARUS_PMT_Data
  sbnobj_Common_PMT_Data
  larcore_Geometry_Geometry_service
  larcorealg_Geometry
  lardataobj_RawData
  art_root_io::TFileService_service
  art_root_io::tfile_support
  art::Framework_Services_Registry
  messagefacility::MF_MessageLogger
  ROOT::Hist
  ROOT::Core
  )


simple_plugin(AsymGaussPulseFunctionTool "tool"
  icaruscode_PMT_Algorithms
  )

simple_plugin(AsymExpPulseFunctionTool "tool"
  icaruscode_PMT_Algorithms
  )

simple_plugin(CustomPulseFunctionTool "tool"
  icaruscode_PMT_Algorithms
  ROOT::Hist
  )


art_make(
  EXCLUDE
    "SimPMTIcarus_module.cc"
    "PMTWaveformBaselines_module.cc"
    "PMTWaveformBaselinesFromReadoutConfiguration_module.cc"
    "AsymGaussPulseFunctionTool_tool.cc"
    "AsymExpPulseFunctionTool_tool.cc"
    "CustomPulseFunctionTool_tool.cc"
          LIB_LIBRARIES 
            lardataobj_RawData
            lardataobj_RecoBase
            lardata_Utilities
            fhiclcpp::fhiclcpp
			cetlib::cetlib
			${ROOT_BASIC_LIB_LIST}
			
	        MODULE_LIBRARIES larcorealg_Geometry
	  		larcore_Geometry_Geometry_service
	        lardata_Utilities
			larevt_Filters
			lardataobj_RawData
			larevt_CalibrationDBI_IOVData
			larevt_CalibrationDBI_Providers
			lardataobj_RecoBase
			lardata_ArtDataHelper
			nusimdata::SimulationBase
            nurandom::RandomUtils_NuRandomService_service
			CLHEP::CLHEP
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
			cetlib::cetlib
			cetlib_except::cetlib_except
            ROOT::Geom
            ROOT::XMLIO
            ROOT::Gdml
			ROOT::FFTW
			${ROOT_BASIC_LIB_LIST}
			
        )

install_headers()
install_fhicl()
install_source()
install_scripts(
  SUBDIRS scripts
  EXTRAS
    "scripts/MergePhotonLibrary.C"
    "scripts/RecompressROOTtrees.C"
  )

