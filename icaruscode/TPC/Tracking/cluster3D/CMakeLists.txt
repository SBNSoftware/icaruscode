
# where should the scripts/..xml file be installed?  Perhaps in bin?

art_make(LIB_LIBRARIES
           lardataalg_DetectorInfo
           lardataobj_RecoBase
           larcorealg_Geometry
           lardata_ArtDataHelper
           lardata_RecoObjects
           art::Framework_Core
           art::Framework_Principal
           art::Framework_Services_Registry
           art_root_io::tfile_support
           art::Utilities
           ROOT::Core
           ROOT::Gpad
           ROOT::Hist
           ROOT::Matrix
           ROOT::Physics
           ROOT::Tree
           art_root_io::TFileService_service
           art::Persistency_Provenance
           canvas::canvas
           messagefacility::MF_MessageLogger
           fhiclcpp::fhiclcpp
           cetlib::cetlib
           cetlib_except::cetlib_except
          TOOL_LIBRARIES larreco_RecoAlg_Cluster3DAlgs
        )

install_headers()
install_fhicl()
install_source()
