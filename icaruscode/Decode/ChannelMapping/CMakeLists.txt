cet_enable_asserts()

art_make(
          EXCLUDE "PMTChannelMapDumper.cxx"
          LIB_LIBRARIES
                        art::Utilities
                        canvas::canvas
                        fhiclcpp::fhiclcpp
                        messagefacility::MF_MessageLogger
                        cetlib_except::cetlib_except
                        cetlib::cetlib
          SERVICE_LIBRARIES
                        icaruscode_Decode_ChannelMapping
                        lardata_Utilities
                        lardata_ArtDataHelper
                        ${ROOT_BASIC_LIB_LIST}
                        LIBWDA
                        art::Framework_Core
                        art::Framework_Principal
                        art::Framework_Services_Registry
                        art_root_io::tfile_support
                        art_root_io::TFileService_service
                        art::Utilities
                        art::Persistency_Common canvas
                        art::Persistency_Provenance canvas
                        art::Utilities canvas
                        messagefacility::MF_MessageLogger
                        
                        fhiclcpp::fhiclcpp
                        cetlib::cetlib cetlib_except
          TOOL_LIBRARIES
                        lardata_Utilities
                        lardata_ArtDataHelper
                        ${ROOT_BASIC_LIB_LIST}
                        LIBWDA
                        art::Framework_Core
                        art::Framework_Principal
                        art::Framework_Services_Registry
                        art_root_io::tfile_support
                        art_root_io::TFileService_service
                        art::Utilities
                        art::Persistency_Common canvas
                        art::Persistency_Provenance canvas
                        art::Utilities canvas
                        messagefacility::MF_MessageLogger
                        
                        fhiclcpp::fhiclcpp
                        cetlib::cetlib cetlib_except
                        SQLite::SQLite3
        )

art_make_exec(NAME "PMTChannelMapDumper"
  LIBRARIES
    icaruscode_Decode_ChannelMapping
    messagefacility::MF_MessageLogger
    fhiclcpp::fhiclcpp
    cetlib::cetlib
    cetlib_except::cetlib_except
    Boost::filesystem
  )

install_headers()
install_fhicl()
install_source()

