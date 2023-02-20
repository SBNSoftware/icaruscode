
# Build the module
art_make_library()
set( MODULE_LIBRARIES
           icaruscode::TPC_Calorimetry_Algorithms
           larcorealg::Geometry
           larreco::Calorimetry
           lardataobj::RecoBase
           lardata::ArtDataHelper
           lardataalg::DetectorInfo
           larcorealg::Geometry
           lardataobj::AnalysisBase
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
           )
set( TOOL_LIBRARIES
           larevt::CalibrationDBI_Providers
           larcorealg::Geometry
           larreco::Calorimetry
           lardataobj::RecoBase
           lardata::ArtDataHelper
           lardataalg::DetectorInfo
           larcorealg::Geometry
           lardataobj::AnalysisBase
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
           wda::wda
         )

cet_build_plugin(NormalizeDriftSQLite art::tool LIBRARIES ${TOOL_LIBRARIES})
cet_build_plugin(NormalizeDrift art::tool LIBRARIES ${TOOL_LIBRARIES})
cet_build_plugin(NormalizeTPCSQL art::tool LIBRARIES ${TOOL_LIBRARIES})
cet_build_plugin(NormalizeTPC art::tool LIBRARIES ${TOOL_LIBRARIES})
cet_build_plugin(NormalizeTPCLocal art::tool LIBRARIES ${TOOL_LIBRARIES})
cet_build_plugin(NormalizeWire art::tool LIBRARIES ${TOOL_LIBRARIES})
cet_build_plugin(NormalizeYZSQL art::tool LIBRARIES ${TOOL_LIBRARIES})
cet_build_plugin(NormalizeYZ art::tool LIBRARIES ${TOOL_LIBRARIES})

install_headers()
install_fhicl()
install_source()
