cet_enable_asserts()

art_make_library(
	EXCLUDE "PMTChannelMapDumper.cxx"
	LIBRARIES
		art::Utilities
		canvas::canvas
		fhiclcpp::fhiclcpp
		messagefacility::MF_MessageLogger
		cetlib_except::cetlib_except
		cetlib::cetlib
		larcorealg::CoreUtils
		
		# for SQLite algorithms (includes duplicates above)
		icarusalg::Utilities
		larevt::CalibrationDBI_Providers
		larcorealg::CoreUtils
		messagefacility::MF_MessageLogger
		fhiclcpp::fhiclcpp
		cetlib::cetlib
		cetlib_except::cetlib_except
		SQLite::SQLite3
		
)

set(	SERVICE_LIBRARIES
		icaruscode::Decode_ChannelMapping
		lardata::Utilities
		lardata::ArtDataHelper
		larcorealg::CoreUtils
		wda::wda
		art::Framework_Core
		art::Framework_Principal
		art::Framework_Services_Registry
		art_root_io::tfile_support
		art_root_io::TFileService_service
		art::Utilities
		art::Persistency_Common
		art::Persistency_Provenance
		art::Utilities canvas::canvas
		messagefacility::MF_MessageLogger
		messagefacility::headers
		fhiclcpp::fhiclcpp
		cetlib::cetlib cetlib_except::cetlib_except
		larevt::CalibrationDBI_IOVData
		larevt::CalibrationDBI_Providers
)
set(	TOOL_LIBRARIES
		larcorealg::CoreUtils
		lardata::Utilities
		lardata::ArtDataHelper
		wda::wda
		art::Framework_Core
		art::Framework_Principal
		art::Framework_Services_Registry
		art_root_io::tfile_support
		art_root_io::TFileService_service
		art::Utilities
		art::Persistency_Common canvas::canvas
		art::Persistency_Provenance
		art::Utilities
		messagefacility::MF_MessageLogger
		messagefacility::headers
		fhiclcpp::fhiclcpp
		cetlib::cetlib cetlib_except::cetlib_except
		SQLite::SQLite3
		larevt::CalibrationDBI_IOVData
		larevt::CalibrationDBI_Providers
	)

cet_build_plugin(ChannelMapPostGres art::tool LIBRARIES ${TOOL_LIBRARIES})
# cet_build_plugin(ChannelMapSQLite art::tool LIBRARIES ${TOOL_LIBRARIES})
cet_build_plugin(ICARUSChannelMap art::service LIBRARIES ${SERVICE_LIBRARIES})

cet_build_plugin(ICARUSChannelMapSQLite art::service
	LIBRARIES
		icaruscode::Decode_ChannelMapping
		art::Framework_Principal
		messagefacility::MF_MessageLogger
		messagefacility::headers
	)

cet_make_exec(NAME "PMTChannelMapDumper"
	LIBRARIES
		icaruscode::Decode_ChannelMapping
		larevt::CalibrationDBI_IOVData
		messagefacility::MF_MessageLogger
		fhiclcpp::fhiclcpp
		cetlib::cetlib
		cetlib_except::cetlib_except
		Boost::filesystem
	) 

install_headers()
install_fhicl()
install_source()

add_subdirectory("Legacy")
