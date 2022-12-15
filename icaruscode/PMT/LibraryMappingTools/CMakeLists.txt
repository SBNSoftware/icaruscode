art_make_library(
		SOURCE
			"ICARUSPhotonMappingTransformations.cxx"
		LIBRARIES
			larsim::PhotonMappingTransformations
			larcore::Geometry_Geometry_service
			larcorealg::Geometry
			lardataalg::DetectorInfo
			art::Framework_Services_Registry
			art::Utilities
			canvas::canvas
			cetlib_except::cetlib_except
			messagefacility::MF_MessageLogger
			fhiclcpp::fhiclcpp
			ROOT::GenVector
			Boost::system
		)

cet_build_plugin(ICARUSPhotonMappingTransformations art::tool LIBRARIES
		icaruscode::PMT_LibraryMappingTools
		larsim::PhotonMappingTransformations
		)

install_headers()
install_fhicl()
install_source()
