add_subdirectory(Algorithms)

# Build the module
art_make( MODULE_LIBRARIES
           icaruscode_TPC_Calorimetry_Algorithms
           larcorealg_Geometry
           larreco_Calorimetry
           lardataobj_RecoBase
           lardata_ArtDataHelper
           lardataalg_DetectorInfo
           larcorealg_Geometry
           lardataobj_AnalysisBase
           art::Framework_Services_Registry
           art::Utilities
           canvas::canvas
           fhiclcpp::fhiclcpp
           cetlib_except::cetlib_except
           ROOT::Core
           ROOT::Physics
           ROOT::Core
           ROOT::Hist
           ROOT::Physics
           messagefacility::MF_MessageLogger
         )

install_headers()
install_fhicl()
install_source()
