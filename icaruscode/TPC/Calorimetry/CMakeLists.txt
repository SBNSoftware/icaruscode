add_subdirectory(Algorithms)

# Build the module
art_make( MODULE_LIBRARIES
           icaruscode_TPC_Calorimetry_Algorithms
           larcorealg_Geometry
           larreco_Calorimetry
           lardataobj_RecoBase
           lardata_ArtDataHelper
           larcorealg_Geometry
           lardataobj_AnalysisBase
           ${ART_FRAMEWORK_SERVICES_REGISTRY}
           ${ART_UTILITIES}
           canvas
           ${FHICLCPP}
           cetlib_except
           ROOT::Core
           ROOT::Physics
           ROOT::Core
           ROOT::Hist
           ROOT::Physics
           ${MF_MESSAGELOGGER}
         )

install_headers()
install_fhicl()
install_source()

