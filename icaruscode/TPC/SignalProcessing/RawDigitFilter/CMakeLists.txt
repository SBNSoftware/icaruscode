add_subdirectory(Algorithms)

art_make_library(
		LIBRARIES
			lardataobj::RawData
			lardataobj::RecoBase
			lardata::Utilities
			fhiclcpp::fhiclcpp
			cetlib::cetlib cetlib_except::cetlib_except)

set(		MODULE_LIBRARIES
			larcorealg::Geometry
			icaruscode::TPC_SignalProcessing_RawDigitFilter_Algorithms
			larcore::Geometry_Geometry_service
			lardata::Utilities
			larevt::Filters
			lardataobj::RawData
			lardataalg::DetectorInfo
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
			canvas::canvas
			art::Persistency_Provenance
			art::Utilities
			messagefacility::MF_MessageLogger
			messagefacility::headers
			fhiclcpp::fhiclcpp
			cetlib::cetlib cetlib_except::cetlib_except
			ROOT::Geom
			ROOT::XMLIO
			ROOT::Gdml
			ROOT::FFTW
			ROOT::Core
			ROOT::Tree
			FFTW3::FFTW3
	)
cet_build_plugin(RawDigitFilterICARUS art::module LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(RawDigitSmoother art::module LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(TPCNoiseMC art::module LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(TPCNoise art::module LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(mtRawDigitFilterICARUS art::module LIBRARIES ${MODULE_LIBRARIES})


install_headers()
install_fhicl()
install_source()
