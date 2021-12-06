
add_definitions(-DEIGEN_FFTW_DEFAULT)

art_make( 
          LIB_LIBRARIES 
                          lardataobj_RawData
                          lardataobj_RecoBase
                          lardata_Utilities
                          fhiclcpp::fhiclcpp
                          cetlib::cetlib cetlib_except
                          ${ROOT_BASIC_LIB_LIST}
                          larcorealg_Geometry
                          larcore_Geometry_Geometry_service
                          lardata_Utilities
                          larevt_Filters
                          lardataobj_RawData
                          larevt_CalibrationDBI_IOVData
                          larevt_CalibrationDBI_Providers
                          lardataobj_RecoBase
                          icaruscode_TPC_Utilities_SignalShapingICARUSService_service
                          ${ICARUS_FFTW_LIBRARIES}
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

          TOOL_LIBRARIES  icaruscode_TPC_SignalProcessing_RawDigitFilter_Algorithms
                          icaruscode_TPC_Utilities_SignalShapingICARUSService_service
                          lardataobj_RawData
                          lardataobj_RecoBase
                          lardata_Utilities
                          fhiclcpp::fhiclcpp
                          cetlib::cetlib cetlib_except
                          larcorealg_Geometry
                          larcore_Geometry_Geometry_service
                          larevt_CalibrationDBI_IOVData
                          larevt_CalibrationDBI_Providers
                          ${ICARUS_FFTW_LIBRARIES}
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
                          ROOT::FFTW
                          ${ROOT_BASIC_LIB_LIST}
       )

install_headers()
install_fhicl()
install_source()
