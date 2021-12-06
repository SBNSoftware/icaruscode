add_subdirectory(HitFinderTools)

art_make( 
          LIB_LIBRARIES 
			           lardataobj_RawData
                       lardataobj_RecoBase
	                   lardata_Utilities
	                   fhiclcpp::fhiclcpp
			           cetlib::cetlib
			           ${ROOT_BASIC_LIB_LIST}
	  MODULE_LIBRARIES larcorealg_Geometry
	  		           larcore_Geometry_Geometry_service
	                   lardata_Utilities
			           larevt_Filters
			           lardataobj_RawData
			           larevt_CalibrationDBI_IOVData
			           larevt_CalibrationDBI_Providers
			           lardataobj_RecoBase
			           lardata_ArtDataHelper
			           larreco_RecoAlg
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
			           
	                   fhiclcpp::fhiclcpp
			           cetlib::cetlib
                       ROOT::Geom
                       ROOT::XMLIO
                       ROOT::Gdml
			           ROOT::FFTW
			           ${ROOT_BASIC_LIB_LIST}
        )

install_headers()
install_fhicl()
install_source()

