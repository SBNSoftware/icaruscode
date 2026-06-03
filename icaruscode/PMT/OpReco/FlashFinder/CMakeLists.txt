
art_make_library(
  LIBRARIES
  INTERFACE
    fhiclcpp::fhiclcpp
  PUBLIC
    icarusalg::Geometry
    larcore::Geometry_Geometry_service
    larcorealg::Geometry
    art::Framework_Services_Registry
)

cet_build_plugin(ICARUSFlashFinder art::module
  LIBRARIES
    icaruscode::PMT_OpReco_FlashFinder
    icaruscode::PMT_OpReco_Algorithms
    lardataobj::RecoBase
  )


install_headers()
install_source()

