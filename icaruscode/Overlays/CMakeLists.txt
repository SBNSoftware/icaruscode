find_package(gallery         REQUIRED)

cet_build_plugin(
  SimInfoOverlayFilter art::EDFilter
  LIBRARIES
  PRIVATE
  larcoreobj::SummaryData
  lardataobj::Simulation
  nusimdata::SimulationBase
  gallery::gallery
)

install_headers()
install_source()
install_fhicl()