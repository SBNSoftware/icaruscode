if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	cet_add_compiler_flags(CXX -Wno-maybe-uninitialized)
endif()

#Add the subfolder containing the deconvolution tools
add_subdirectory(DeconTools)
add_subdirectory(ROITools)
add_subdirectory(SimTestPulse)

art_make_library(
		LIBRARIES
			lardataobj::RawData
			lardataobj::RecoBase
			lardata::Utilities
			fhiclcpp::fhiclcpp
			cetlib::cetlib
)
set(		MODULE_LIBRARIES
			icaruscode_TPC_Utilities
			icaruscode_IcarusObj
			icarus_signal_processing::icarus_signal_processing
			icarus_signal_processing::Detection
			icarus_signal_processing::Filters
			larcorealg::Geometry
			larcore::Geometry_Geometry_service
			lardata::Utilities
			lardataalg::DetectorInfo
			larevt::Filters
			lardataobj::Simulation
			lardataobj::RawData
			larevt::CalibrationDBI_IOVData
			larevt::CalibrationDBI_Providers
			lardataobj::RecoBase
			lardata::ArtDataHelper
			icaruscode::TPC_Utilities_SignalShapingICARUSService_service
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
			messagefacility::headers
			fhiclcpp::fhiclcpp
			cetlib::cetlib
			ROOT::FFTW
			FFTW3::FFTW3
		)
cet_build_plugin(Decon1DROI art::module LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(ROIConverter art::module LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(ROIFinder art::module LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(SimChannelROI art::module LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(RecoWireICARUSRaw art::module LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(RecoWireICARUS art::module LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(RecoWireROIICARUS art::module LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(RecoWireROI art::module LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(WaveformIntegrity art::module LIBRARIES ${MODULE_LIBRARIES})

install_headers()
install_fhicl()
install_source()
