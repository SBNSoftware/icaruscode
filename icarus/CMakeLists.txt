
set( IFDH_ART_DIR $ENV{IFDH_ART_DIR} )

add_subdirectory(Geometry)



#add_subdirectory(RawData)
#add_subdirectory(BeamData)
#add_subdirectory(BeamDAQ)
#add_subdirectory(BeamAna)
#add_subdirectory(Utilities)
#add_subdirectory(Database)
#add_subdirectory(Simulation)
#add_subdirectory(LArG4)
#add_subdirectory(PhotonPropagation)
#add_subdirectory(DetSim)
#add_subdirectory(EventGenerator)
#add_subdirectory(CalData)
#add_subdirectory(TriggerSim)
#add_subdirectory(AnalysisTree)
#add_subdirectory(Calibrations)
#add_subdirectory(OpticalDetectorSim)
#add_subdirectory(OpticalDetectorAna)
#add_subdirectory(SpaceCharge)
#add_subdirectory(SpaceChargeServices)
#add_subdirectory(TPCNeutrinoIDFilter)
#add_subdirectory(DataScanner)
#add_subdirectory(LiteMaker)
#add_subdirectory(MuCS)
#add_subdirectory(T0Reco)
#add_subdirectory(EventWeight)
#add_subdirectory(UBFlashFinder)

#add_subdirectory(DataOverlay)
#add_subdirectory(PatternFilter)

# Respect the order specified in lar_build, which is:
#
#     "SimpleTypesAndConstants",
#     "larcore_Geometry",
#     "larcoreobj_SummaryData",
#     "lardataobj_RawData",
#     "lardataobj_RecoBase",
#     "lardataobj_OpticalDetectorData",
#     "lardataobj_AnalysisBase",
#     "lardata_Utilities",
#     "larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation",
#     "larsim_MCCheater",
#     "larreco_Calorimetry",
#     "lardata_RecoObjects larreco_Deprecated",
#     "larreco_RecoAlg",
#     "larsim_LArG4",
#     "larsim_PhotonPropagation",
#     "larsim_DetSim",
#     "larsim_EventGenerator",
#     "larevt_Filters",
#     "CalData",
#     "larreco_ClusterFinder",
#     "larreco_Genfit",
#     "larreco_HitFinder",
#     "larreco_TrackFinder",
#     "larreco_VertexFinder",
#     "larreco_EventFinder",
#     "larsim_TriggerAlgo",
#     "larana_OpticalDetector",
#     "larreco_Calorimetry", 
#     "larreco_ShowerFinder",
#     "lareventdisplay_EventDisplay",
#     "larana_ParticleIdentification",
#     "larexamples_AnalysisExample",
#     "LArPandoraAlgorithms",
#     "larpandora_LArPandoraInterface"
