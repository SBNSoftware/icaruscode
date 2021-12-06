add_subdirectory(Algorithms)

add_definitions(-DEIGEN_FFTW_DEFAULT)

art_make(
          LIB_LIBRARIES 
                        lardataobj_RawData
                        lardataobj_RecoBase
                        lardata_Utilities
                        fhiclcpp::fhiclcpp
                        cetlib::cetlib cetlib_except
                        ${ROOT_BASIC_LIB_LIST}
          MODULE_LIBRARIES
                        larcorealg_Geometry
                        icaruscode_TPC_SignalProcessing_RawDigitFilter_Algorithms
                        larcore_Geometry_Geometry_service
                        lardata_Utilities
                        larevt_Filters
                        lardataobj_RawData
                        lardataalg_DetectorInfo
                        larevt_CalibrationDBI_IOVData
                        larevt_CalibrationDBI_Providers
                        lardataobj_RecoBase
                        lardata_ArtDataHelper
                        icaruscode_TPC_Utilities_SignalShapingICARUSService_service
                        art::Framework_Core
                        art::Framework_Principal
                        art::Framework_Services_Registry
                        art_root_io::tfile_support ROOT::Core
                        art_root_io::TFileService_service
                        art::Persistency_Common canvas
                        art::Persistency_Provenance canvas
                        art::Utilities canvas
                        messagefacility::MF_MessageLogger
                        
                        fhiclcpp::fhiclcpp
                        cetlib::cetlib cetlib_except
                        ROOT::Geom
                        ROOT::XMLIO
                        ROOT::Gdml
			ROOT::FFTW
			${ROOT_BASIC_LIB_LIST}
			${TBB}
        )

install_headers()
install_fhicl()
install_source()
