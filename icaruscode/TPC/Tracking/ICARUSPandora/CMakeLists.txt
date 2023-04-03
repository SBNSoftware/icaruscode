
# where should the scripts/..xml file be installed?  Perhaps in bin?

set( TOOL_LIBRARIES
	larpandora::LArPandoraInterface
	lardataobj::RecoBase
	canvas::canvas
	art::Framework_Principal
	messagefacility::MF_MessageLogger
	cetlib::cetlib
)

cet_build_plugin(LArPandoraHitCollectionToolICARUS art::tool LIBRARIES ${TOOL_LIBRARIES})

install_fhicl()
install_source()

add_subdirectory(scripts)

