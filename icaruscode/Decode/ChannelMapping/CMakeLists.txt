cet_enable_asserts()

art_make_library(
  EXCLUDE "ChannelMapDumper.cxx"
  LIBRARIES
    icarusalg::Utilities
    larcorealg::CoreUtils
    art::Framework_Services_Registry
    fhiclcpp::fhiclcpp
    messagefacility::MF_MessageLogger
    cetlib::cetlib
    cetlib_except::cetlib_except
    
    # for PostgreSQL algorithm and service provider (includes duplicates above)
    icarusalg::Utilities
    larcorealg::CoreUtils
    messagefacility::MF_MessageLogger
    fhiclcpp::fhiclcpp
    cetlib_except::cetlib_except
    wda::wda
    
    # for SQLite algorithm and service provider (includes duplicates above)
    icarusalg::Utilities
    larevt::CalibrationDBI_Providers
    larcorealg::CoreUtils
    messagefacility::MF_MessageLogger
    fhiclcpp::fhiclcpp
    cetlib::cetlib
    cetlib_except::cetlib_except
    SQLite::SQLite3
    
)

cet_build_plugin(ICARUSChannelMapSQLite art::service
  LIBRARIES
    icaruscode::Decode_ChannelMapping
    art::Framework_Principal
    messagefacility::MF_MessageLogger
    messagefacility::headers
  )

add_subdirectory("Legacy")

cet_make_exec(NAME "ChannelMapDumper"
  LIBRARIES
    icaruscode::Decode_ChannelMapping_Legacy
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
