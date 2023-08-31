
art_make_library(
	LIBRARIES
		art_root_io::RootDB
		SQLite::SQLite3
		lardata::Utilities
		canvas::canvas
		cetlib_except::cetlib_except
		cetlib::cetlib
		fhiclcpp::fhiclcpp
		hep_concurrency::hep_concurrency # required by fhicl::ParameterSetRegistry
		ROOT::Tree
		ROOT::Core
		ROOT::RIO
                sbnobj::Common_CRT
	)

cet_build_plugin(SaveConfigurationIntoTFile art::module LIBRARIES
		art_root_io::TFileService_service
		art_root_io::tfile_support
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


cet_build_plugin(DuplicateEventTracker art::service LIBRARIES
		icaruscode::Utilities
		art::Framework_Principal
		messagefacility::MF_MessageLogger
	)

cet_build_plugin(DumpCRTPMTMatching art::module LIBRARIES
		icaruscode::Utilities
		icaruscode::IcarusObj
		sbnobj::Common_CRT
		lardataobj::RecoBase
		canvas::canvas
	)

cet_build_plugin(DummyFilter art::module)
cet_build_plugin(DummyProducer art::module)
cet_build_plugin(DummyAnalyzer art::module)
cet_build_plugin(FilterOnArtPathOutcome art::module)


install_headers()
install_source()
install_fhicl()
