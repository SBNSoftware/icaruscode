
art_make( 
             LIB_LIBRARIES    larsim_MCCheater_BackTrackerService_service
                              larsim_MCCheater_ParticleInventoryService_service
                              lardataobj_RecoBase
                              lardata_ArtDataHelper
                              larsim_Simulation 
                              nug4::ParticleNavigation lardataobj_Simulation
                              larevt_Filters
                              larcorealg_Geometry
                              larcore_Geometry_Geometry_service
                              lardata_RecoObjects
                              lardataobj_AnalysisBase
                              larreco_Calorimetry
                              lardata_Utilities
                              larreco_RecoAlg_ClusterRecoUtil
                              larreco_RecoAlg_CMTool_CMToolBase
                              larreco_RecoAlg_CMTool_CMTAlgMerge
                              larreco_RecoAlg_CMTool_CMTAlgPriority
                              larreco_RecoAlg_CMTool_CMTAlgMatch
                              larreco_RecoAlg_CMTool_CMToolApp
                              ROOT::MathCore
                              ${ROOT_GENVECTOR}
                              ROOT::Geom
                              ${ROOT_BASIC_LIB_LIST}
                              larreco_RecoAlg_TCAlg
                              larreco_RecoAlg_ImagePatternAlgs_DataProvider
                              larreco_RecoAlg_PMAlg
                              nusimdata::SimulationBase
                              art::Framework_Core
                              art::Framework_Principal
                              art::Framework_Services_Registry
                              ${ART_FRAMEWORK_SERVICES_BASIC}
                              ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                              
                              art::Framework_Services_Optional_RandomNumberGenerator_service
                              art::Persistency_Common
                              art::Persistency_Provenance
                              art::Utilities
                              canvas::canvas
                              messagefacility::MF_MessageLogger
                              fhiclcpp::fhiclcpp
                              cetlib::cetlib cetlib_except
                              CLHEP::CLHEP
                              ${ROOT_BASIC_LIB_LIST}
                              ${ROOT_MINUIT}
                              ${ROOT_MINUIT2}
             MODULE_LIBRARIES icaruscode_TPC_Tracking_MCS
                              larsim_MCCheater_BackTrackerService_service
                              larsim_MCCheater_ParticleInventoryService_service
                              lardataobj_RecoBase
                              lardata_ArtDataHelper
                              larsim_Simulation 
                              nug4::ParticleNavigation lardataobj_Simulation
                              larevt_Filters
                              larcorealg_Geometry
                              larcore_Geometry_Geometry_service
                              lardata_RecoObjects
                              lardataobj_AnalysisBase
                              larreco_Calorimetry
                              lardata_Utilities
                              larreco_RecoAlg_ClusterRecoUtil
                              larreco_RecoAlg_CMTool_CMToolBase
                              larreco_RecoAlg_CMTool_CMTAlgMerge
                              larreco_RecoAlg_CMTool_CMTAlgPriority
                              larreco_RecoAlg_CMTool_CMTAlgMatch
                              larreco_RecoAlg_CMTool_CMToolApp
                              ROOT::MathCore
                              ${ROOT_GENVECTOR}
                              ROOT::Geom
                              ${ROOT_BASIC_LIB_LIST}
                              larreco_RecoAlg_TCAlg
                              larreco_RecoAlg_ImagePatternAlgs_DataProvider
                              larreco_RecoAlg_PMAlg
                              nusimdata::SimulationBase
                              art::Framework_Core
                              art::Framework_Principal
                              art::Framework_Services_Registry
                              ${ART_FRAMEWORK_SERVICES_BASIC}
                              ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                              
                              art::Framework_Services_Optional_RandomNumberGenerator_service
                              art::Persistency_Common
                              art::Persistency_Provenance
                              art::Utilities
                              canvas::canvas
                              messagefacility::MF_MessageLogger
                              fhiclcpp::fhiclcpp
                              cetlib::cetlib cetlib_except
                              CLHEP::CLHEP
                              ${ROOT_BASIC_LIB_LIST}
                              ${ROOT_MINUIT}
                              ${ROOT_MINUIT2}
        )


install_headers()
install_fhicl()
install_source()

