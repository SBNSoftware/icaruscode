
art_make(
  EXCLUDE
    "SaveConfigurationIntoTFile_module.cc"
    "DummyFilter_module.cc"
    "DummyProducer_module.cc"
    "DuplicateEventTracker_service.cc"
  LIB_LIBRARIES
    art_root_io::RootDB
    SQLite::SQLite3
    cetlib_except::cetlib_except
    cetlib::cetlib
    fhiclcpp
    hep_concurrency # required by fhicl::ParameterSetRegistry
    ROOT::Tree
    ROOT::Core
    ROOT::RIO
)

simple_plugin(SaveConfigurationIntoTFile "module"
              art_root_io::TFileService_service
              art_root_io::tfile_support ROOT::Core
              art::Framework_Services_Registry
              art::Framework_Principal
              art::Framework_Core
              canvas::canvas
              messagefacility::MF_MessageLogger
              hep_concurrency
              fhiclcpp::fhiclcpp
              cetlib::cetlib
              ROOT::Core
)


simple_plugin(DuplicateEventTracker "service"
              icaruscode_Utilities
              art::Framework_Principal
              messagefacility::MF_MessageLogger
)

simple_plugin(DummyFilter "module")
simple_plugin(DummyProducer "module")


install_headers()
install_source()
install_fhicl()
