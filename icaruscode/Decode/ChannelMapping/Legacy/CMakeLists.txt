cet_enable_asserts()

cet_build_plugin(ChannelMapSQLite art::tool
	LIBRARIES
		icaruscode::Decode_ChannelMapping
	)

cet_build_plugin(ChannelMapPostGres art::tool
	LIBRARIES
		icaruscode::Decode_ChannelMapping
	)

install_headers()
install_fhicl()
install_source()

