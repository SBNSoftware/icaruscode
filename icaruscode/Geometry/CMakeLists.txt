art_make_library()
set(	SERVICE_LIBRARIES
		icarusalg::Geometry
		larcorealg::Geometry
		larcore::Geometry_Geometry_service
		art::Utilities
		messagefacility::MF_MessageLogger
		fhiclcpp::fhiclcpp
		cetlib::cetlib
		cetlib_except::cetlib_except
		ROOT::Geom
		ROOT::XMLIO
		ROOT::Gdml
)
set(	TOOL_LIBRARIES
		icarusalg::Geometry
		larcorealg::Geometry
		art::Utilities
		messagefacility::headers
		fhiclcpp::fhiclcpp
		cetlib::cetlib
		cetlib_except::cetlib_except
		ROOT::Geom
		ROOT::XMLIO
		ROOT::Gdml

	)

cet_build_plugin(AuxDetGeoObjectSorterICARUS art::tool LIBRARIES PRIVATE ${TOOL_LIBRARIES})
cet_build_plugin(GeoObjectSorterPMTasTPC art::tool LIBRARIES PRIVATE ${TOOL_LIBRARIES})
cet_build_plugin(WireReadoutICARUS art::service LIBRARIES PRIVATE ${SERVICE_LIBRARIES})
