cet_enable_asserts()

add_subdirectory(tools)


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
	FFTW3::FFTW3
	art::Framework_Core
	art::Framework_Principal
	art::Framework_Services_Registry
	art::Persistency_Common
	art::Framework_IO
	art::Framework_IO_detail
	art::Persistency_Provenance
	sbnobj::ICARUS_TPC	
	art_root_io::tfile_support ROOT::Core
	art_root_io::RootDB
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

art_make_library(
        LIBRARIES ${icarus_util_lib_list}
)

cet_build_plugin(FileCatalogMetadataICARUS art::service LIBRARIES ${icarus_util_lib_list})
cet_build_plugin(SignalShapingICARUSService art::service LIBRARIES ${icarus_util_lib_list})
cet_build_plugin(SIOVChannelStatusICARUSService art::service LIBRARIES ${icarus_util_lib_list})
cet_build_plugin(TFileMetadataICARUS art::service LIBRARIES ${icarus_util_lib_list})


install_headers()
install_source()
install_fhicl()
