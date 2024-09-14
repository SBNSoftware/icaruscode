find_package(gallery         REQUIRED)

#find_package(larwirecell REQUIRED )
#find_package(jsoncpp REQUIRED)
#find_package(spdlog REQUIRED)

include_directories( $ENV{JSONCPP_INC} )
link_directories( $ENV{JSONCPP_LIB} )
include_directories( $ENV{SPDLOG_INC} )
link_directories( $ENV{SPDLOG_LIB} )
include_directories( $ENV{WIRECELL_INC} )
link_directories( $ENV{WIRECELL_LIB} )
include_directories( $ENV{LARWIRECELL_INC} )
link_directories( $ENV{LARWIRECELL_LIB} )

cet_build_plugin(
  SimInfoOverlayFilter art::EDFilter
  LIBRARIES
  PRIVATE
  larcoreobj::SummaryData
  lardataobj::Simulation
  nusimdata::SimulationBase
  gallery::gallery
)

# This is what's needed for building the WireCellubsim WCT component
art_make_library(
  LIBRARY_NAME ICARUSWCDrifter
  SOURCE ICARUSDrifter.cxx
  LIBRARIES
        art::Framework_Services_Registry
        art::Utilities
        ${Boost_SYSTEM_LIBRARY}
        ${ROOT_BASIC_LIB_LIST}
        ${JSONCPP}
	larevt::CalibrationDBI_Providers
	lardata::ArtDataHelper
        ${FHICLCPP}
        cetlib_except::cetlib_except
	WireCellGen
	WireCellUtil
	WireCellAux
	WireCellIface
	WireCellPgraph
	WireCellRoot
	jsoncpp
)


install_headers()
install_source()
install_fhicl()
