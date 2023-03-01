
# where should the scripts/..xml file be installed?  Perhaps in bin?

set( TOOL_LIBRARIES
	larpandora::LArPandora
	larpandora::GeometryComponents
	larpandoracontent::LArPandoraContent
	lardata::AssociationUtil
	lardata::Utilities
	lardataobj::AnalysisBase
	lardataobj::RecoBase
	lardataobj::Simulation
	larcoreobj::SimpleTypesAndConstants
	nusimdata::SimulationBase
	art::Persistency_Common
	canvas::canvas
	PandoraPFA::PandoraSDK
	larpandora::GetDetectorType
	larreco::ClusterFinder
	larreco::ClusterParamsImportWrapper
	larreco::RecoAlg_ClusterRecoUtil
	larevt::ChannelStatusProvider
	larevt::ChannelStatusService
	larsim::MCCheater_ParticleInventoryService_service
	lardata::DetectorClocksService
	lardata::DetectorPropertiesService
	larcore::Geometry_Geometry_service
	larcorealg::Geometry
	art_root_io::TFileService_service
	art::Framework_Principal
	art::Framework_Services_Registry
	messagefacility::MF_MessageLogger
	cetlib::cetlib
	cetlib_except::cetlib_except
	PandoraPFA::PandoraSDK
)

cet_build_plugin(LArPandoraHitCollectionToolICARUS art::tool LIBRARIES ${TOOL_LIBRARIES})

install_fhicl()
install_source()

add_subdirectory(scripts)

