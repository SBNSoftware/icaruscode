cet_enable_asserts()

add_subdirectory(tools)
include_directories(${FFTW3_INCLUDE_DIRS})
link_directories(${FFTW3_LIBRARY_DIRS})


add_definitions(-DEIGEN_FFTW_DEFAULT)

set(icarus_util_lib_list
	larcorealg::Geometry
	larevt::CalibrationDBI_IOVData
	larevt::CalibrationDBI_Providers
	larreco::HitFinder
	lardataobj::RecoBase
	lardataobj::RawData
	larcore::Geometry_Geometry_service
	lardata::Utilities
	nurandom::RandomUtils_NuRandomService_service
	${ICARUS_FFTW_LIBRARIES}
	art::Framework_Core
	art::Framework_Principal
	art::Framework_Services_Registry
	art::Persistency_Common
	art::Framework_IO
	art::Framework_IO_detail
	art::Persistency_Provenance
	icaruscode_IcarusObj
	art_root_io::tfile_support ROOT::Core
	art::Framework_Services_Optional_RandomNumberGenerator_service
	art_root_io::TFileService_service
	art::Framework_Services_System_FileCatalogMetadata_service
	art::Framework_Services_System_TriggerNamesService_service
	art::Utilities
	canvas::canvas
	messagefacility::MF_MessageLogger
	messagefacility::headers
	fhiclcpp::fhiclcpp
	cetlib::cetlib
	Boost::filesystem
	Boost::system
	CLHEP::Random
        larcoreobj::headers
)

art_make(
# LIBRARY_NAME icaruscode_TPC_Utilities
        LIB_LIBRARIES ${icarus_util_lib_list}

	SERVICE_LIBRARIES
		${icarus_util_lib_list}

	MODULE_LIBRARIES
		art_root_io::TFileService_service
		art_root_io::tfile_support ROOT::Core
		art::Framework_Services_Registry
		art::Framework_Principal
		art::Framework_Core
		canvas::canvas
		messagefacility::MF_MessageLogger
		hep_concurrency::hep_concurrency
		fhiclcpp::fhiclcpp
		cetlib::cetlib
		ROOT::Core
	)


install_headers()
install_source()
install_fhicl()
