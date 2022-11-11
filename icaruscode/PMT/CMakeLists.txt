add_subdirectory(Data)
add_subdirectory(OpticalTools)
add_subdirectory(OpReco)
add_subdirectory(Algorithms)
add_subdirectory(LibraryMappingTools)
add_subdirectory(Trigger)
#add_subdirectory(Timing)

# Removing AVX instructions because assuming the grid nodes
# being less than 6 year old proved to be pretentious.
# AVX is not so relevant after the last commit anyway
# cet_set_compiler_flags(EXTRA_CXX_FLAGS -mavx)

cet_build_plugin(SimPMTIcarus art::module
                    LIBRARIES
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

cet_build_plugin(OpDetWaveformMetaMaker art::module
                    LIBRARIES
			icaruscode_PMT_Algorithms
			lardataalg::DetectorInfo
			lardataobj::RawData
			lardata::Utilities
			larcore::Geometry_Geometry_service
			)

cet_build_plugin(PMTWaveformBaselines art::module
                    LIBRARIES
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

cet_build_plugin(PMTWaveformBaselinesFromReadoutConfiguration art::module
                    LIBRARIES
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

cet_build_plugin(PMTWaveformBaselinesFromChannelData art::module
                    LIBRARIES
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

cet_build_plugin(AsymGaussPulseFunctionTool art::tool
        LIBRARIES
	icaruscode_PMT_Algorithms
	)

cet_build_plugin(AsymExpPulseFunctionTool art::tool
        LIBRARIES
	icaruscode_PMT_Algorithms
	)

cet_build_plugin(SampledWaveformFunctionTool art::tool
        LIBRARIES
	icaruscode_PMT_Algorithms
	icaruscode_Decode_DecoderTools
  )

cet_build_plugin(CustomPulseFunctionTool art::tool
        LIBRARIES
	icaruscode_PMT_Algorithms
	ROOT::Hist
	)

art_make_library(
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
          LIBRARIES
            lardataobj::RawData
            lardataobj::RecoBase
            lardata::Utilities)
set(      MODULE_LIBRARIES
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
cet_build_plugin(OpHitFinderICARUS art::module LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(PMTStartCalibTime art::module LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(PMTcoordinates art::module LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(PhotonPropogationICARUS art::module LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(TrigInfo art::module LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(wvfAnaicarus art::module LIBRARIES ${MODULE_LIBRARIES})


install_headers()
install_fhicl()
install_source()
install_scripts(
	SUBDIRS scripts
	EXTRAS
		"scripts/MergePhotonLibrary.C"
		"scripts/RecompressROOTtrees.C"
	)

