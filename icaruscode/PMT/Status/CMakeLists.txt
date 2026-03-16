cet_enable_asserts()

set( LIB_LIBRARIES
        art::Framework_Services_Registry
        messagefacility::MF_MessageLogger
        larcorealg::CoreUtils
        larevt::CalibrationDBI_IOVData
        larevt::CalibrationDBI_Providers
        fhiclcpp::fhiclcpp
        cetlib_except::cetlib_except
)

set( SERVICE_LIBRARIES
        icaruscode_PMT_Status
        art::Framework_Principal
        larcore::Geometry_Geometry_service
)

file(GLOB lib_srcs *.cxx)

art_make_library( SOURCE ${lib_srcs} LIBRARIES PUBLIC ${LIB_LIBRARIES})

cet_build_plugin( PMTChannelStatusService art::service LIBRARIES PUBLIC ${SERVICE_LIBRARIES})

install_headers()
install_fhicl()
install_source()
