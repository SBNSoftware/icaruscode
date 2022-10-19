add_subdirectory(Data)
add_subdirectory(OpticalTools)
add_subdirectory(OpReco)
add_subdirectory(Algorithms)
add_subdirectory(LibraryMappingTools)
add_subdirectory(Trigger)

# Removing AVX instructions because assuming the grid nodes
# being less than 6 year old proved to be pretentious.
# AVX is not so relevant after the last commit anyway
# cet_set_compiler_flags(EXTRA_CXX_FLAGS -mavx)

#include_directories(${CPPGSL_INC})

simple_plugin(SimPMTIcarus "module"
			icaruscode_PMT_Algorithms
			lardataobj::RawData
			lardataobj::Simulation
			larcore::Geometry_Geometry_service
			lardata::Utilities
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
			Boost::system
			)

simple_plugin(OpDetWaveformMetaMaker "module"
			icaruscode_PMT_Algorithms
			lardataalg::DetectorInfo
			lardataobj::RawData
			lardata::Utilities
			larcore::Geometry_Geometry_service
			)

simple_plugin(PMTWaveformBaselines "module"
			icaruscode_PMT_Algorithms
			sbnobj::ICARUS_PMT_Data
			larcore::Geometry_Geometry_service
			larcorealg::Geometry
			lardataobj::RawData
			lardata::Utilities
			art_root_io::TFileService_service
			art_root_io::tfile_support
			art::Framework_Services_Registry
			messagefacility::MF_MessageLogger
			ROOT::Hist
			ROOT::Core
			)

simple_plugin(PMTWaveformBaselinesFromReadoutConfiguration "module"
			sbnobj::ICARUS_PMT_Data
			sbnobj::Common_PMT_Data
			larcore::Geometry_Geometry_service
			larcorealg::Geometry
			lardataobj::RawData
			lardata::Utilities
			art_root_io::TFileService_service
			art_root_io::tfile_support
			art::Framework_Services_Registry
			messagefacility::MF_MessageLogger
			ROOT::Hist
			ROOT::Core
			)

simple_plugin(PMTWaveformBaselinesFromChannelData "module"
			icaruscode_PMT_Algorithms
			sbnobj::ICARUS_PMT_Data
			sbnobj::Common_PMT_Data
			larcore::Geometry_Geometry_service
			larcorealg::Geometry
			lardataobj::RawData
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

simple_plugin(SampledWaveformFunctionTool "tool"
	icaruscode_PMT_Algorithms
	icaruscode_Decode_DecoderTools
  )

simple_plugin(CustomPulseFunctionTool "tool"
	icaruscode_PMT_Algorithms
	ROOT::Hist
	)

art_make(
  EXCLUDE
    "SimPMTIcarus_module.cc"
    "OpDetWaveformMetaMaker_module.cc"
    "PMTWaveformBaselines_module.cc"
    "PMTWaveformBaselinesFromReadoutConfiguration_module.cc"
    "PMTWaveformBaselinesFromChannelData_module.cc"
    "AsymGaussPulseFunctionTool_tool.cc"
    "AsymExpPulseFunctionTool_tool.cc"
    "SampledWaveformFunctionTool_tool.cc"
    "CustomPulseFunctionTool_tool.cc"
          LIB_LIBRARIES
            lardataobj::RawData
            lardataobj::RecoBase
            lardata::Utilities
          MODULE_LIBRARIES
            larcorealg::Geometry
            larana::OpticalDetector
            larcore::Geometry_Geometry_service
            lardata::Utilities
            larevt::Filters
            lardataobj::RawData
            larevt::CalibrationDBI_IOVData
            larevt::CalibrationDBI_Providers
            lardataobj::RecoBase
            lardata::ArtDataHelper
            nusimdata::SimulationBase
            nurandom::RandomUtils_NuRandomService_service
            nug4::ParticleNavigation
            art::Framework_Core
            art::Framework_Principal
            art::Framework_Services_Registry
            art::Framework_Services_Optional_RandomNumberGenerator_service
            art_root_io::TFileService_service
            art_root_io::tfile_support
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
            ROOT::Tree
            Boost::system
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

