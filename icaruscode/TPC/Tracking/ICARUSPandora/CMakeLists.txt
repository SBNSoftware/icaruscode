
# where should the scripts/..xml file be installed?  Perhaps in bin?

art_make_library(
                LIBRARIES
                        lardataobj::RecoBase
                        fhiclcpp::fhiclcpp
                        cetlib::cetlib
)
set(            MODULE_LIBRARIES
                        larcorealg::Geometry
                        larcorealg::Geometry
                        larcore::Geometry_Geometry_service
                        lardata::Utilities
                        larevt::Filters
                        lardataobj::RecoBase
                        lardata::ArtDataHelper
                        art::Framework_Core
                        art::Framework_Principal
                        art::Framework_Services_Registry
                        art_root_io::tfile_support ROOT::Core
                        art_root_io::TFileService_service
                        art::Persistency_Common
                        art::Persistency_Provenance
                        art::Utilities
                        canvas::canvas
                        messagefacility::MF_MessageLogger
                        messagefacility::headers
                        fhiclcpp::fhiclcpp
                        cetlib::cetlib
                        ROOT::Geom
                        ROOT::XMLIO
                        ROOT::Gdml
                        ROOT::FFTW

        )
cet_build_plugin(HARPS art::module LIBRARIES ${MODULE_LIBRARIES})

set( TOOL_LIBRARIES
	larpandora::LArPandoraInterface
	lardataobj::RecoBase
	canvas::canvas
	art::Framework_Principal
	messagefacility::MF_MessageLogger
	cetlib::cetlib
)

cet_build_plugin(LArPandoraHitCollectionToolICARUS art::tool LIBRARIES ${TOOL_LIBRARIES})

install_headers()
install_fhicl()
install_source()

add_subdirectory(scripts)

