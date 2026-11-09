cet_enable_asserts()

art_make_library(
	EXCLUDE "ICARUSChannelMap.cc"
	LIBRARIES
		icaruscode::Decode_ChannelMapping
		art::Utilities
		messagefacility::MF_MessageLogger
		fhiclcpp::fhiclcpp
		cetlib::cetlib
		cetlib_except::cetlib_except
)

cet_build_plugin(ChannelMapSQLite art::tool
	LIBRARIES PRIVATE
		icaruscode::Decode_ChannelMapping
	)

cet_build_plugin(ChannelMapPostGres art::tool
	LIBRARIES PRIVATE
		icaruscode::Decode_ChannelMapping
	)

cet_build_plugin(ICARUSChannelMap art::service
	LIBRARIES PUBLIC icaruscode::Decode_ChannelMapping_Legacy # service interface
	LIBRARIES PRIVATE
		icaruscode::Decode_ChannelMapping
		art::Framework_Principal
		messagefacility::MF_MessageLogger
		messagefacility::headers
	)

install_headers()
install_fhicl()
install_source()

