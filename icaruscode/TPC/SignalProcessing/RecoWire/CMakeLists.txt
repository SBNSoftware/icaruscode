if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	cet_add_compiler_flags(CXX -Wno-maybe-uninitialized)
endif()

#Add the subfolder containing the deconvolution tools
add_subdirectory(DeconTools)
add_subdirectory(ROITools)
add_subdirectory(SimTestPulse)

add_definitions(-DEIGEN_FFTW_DEFAULT)
include_directories(${FFTW3_INCLUDE_DIRS})

art_make(
		LIB_LIBRARIES
			lardataobj::RawData
			lardataobj::RecoBase
			lardata::Utilities
			fhiclcpp::fhiclcpp
			cetlib::cetlib

		MODULE_LIBRARIES
			icarus_signal_processing::icarus_signal_processing
			icarus_signal_processing::Detection
			icarus_signal_processing::Filters
			larcorealg::Geometry
			larcore::Geometry_Geometry_service
			lardata::Utilities
			lardataalg::DetectorInfo
			larevt::Filters
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
			${TBB}
		)

install_headers()
install_fhicl()
install_source()
