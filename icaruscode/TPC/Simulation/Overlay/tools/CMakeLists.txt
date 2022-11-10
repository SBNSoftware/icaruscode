cet_enable_asserts()

art_make_library()
set(	TOOL_LIBRARIES
			larcorealg::Geometry
			larevt::CalibrationDBI_IOVData
			larevt::CalibrationDBI_Providers
			larreco::HitFinder
			lardataobj::RecoBase
			larcore::Geometry_Geometry_service
			lardata::Utilities
			lardataalg::DetectorInfo
			nurandom::RandomUtils_NuRandomService_service
			FFTW3::FFTW3
			art::Framework_Core
			art::Framework_Principal
			art::Framework_Services_Registry
			art::Persistency_Common
			art::Persistency_Provenance
			icarus_signal_processing::icarus_signal_processing
			art_root_io::tfile_support ROOT::Core
			art::Framework_Services_Optional_RandomNumberGenerator_service
			art_root_io::TFileService_service
			art::Framework_Services_System_FileCatalogMetadata_service
			art::Utilities
			canvas::canvas
			messagefacility::MF_MessageLogger
			messagefacility::headers
			fhiclcpp::fhiclcpp
			cetlib::cetlib
			Boost::filesystem
			Boost::system
			CLHEP::Random
			Eigen3::Eigen
	)
cet_build_plugin(Overlay1D art::tool LIBRARIES ${TOOL_LIBRARIES} SOURCE Overlay1D_tool.cc)
set_property(SOURCE Overlay1D_tool.cc APPEND PROPERTY COMPILE_DEFINITIONS EIGEN_FFTW_DEFAULT)

install_headers()
install_source()
install_fhicl()
