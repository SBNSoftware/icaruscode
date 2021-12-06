include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
include_directories( $ENV{NURANDOM_INC} )
include_directories( $ENV{IFDHX_FQ_DIR}/inc )

add_subdirectory(CRTUtils)
add_subdirectory(CRTDecoder)

art_make(
    NO_PLUGINS
    EXCLUDE
        CRTChannelMapAlg.cxx
        CRTGeometryHelper_service.cc
        CRTDetSim_module.cc
        CRTSimHitProducer_module.cc
        CRTTrueHitProducer_module.cc
        CRTTzeroProducer_module.cc
        CRTTrackProducer_module.cc
        CRTSimAnalysis_module.cc
        CRTDataAnalysis_module.cc
        CrtOpHitMatchAnalysis_module.cc
        CRTTruthMatchAnalysis_module.cc
        CRTAutoVeto_module.cc
        FlashResAna_module.cc
        PhotBackground_module.cc
    LIBRARY_NAME
        icaruscode_CRTData
    LIB_LIBRARIES
        larcorealg_Geometry
        larcore_Geometry_Geometry_service
        art::Framework_Core
        art::Framework_IO_Sources
        art::Framework_Principal
        art_root_io::TFileService_service
        canvas::canvas
        art::Framework_Services_Registry
        messagefacility::MF_MessageLogger
        
        fhiclcpp::fhiclcpp
        ${ROOT_BASIC_LIB_LIST}
        ROOT::Geom
        cetlib::cetlib
)

art_make_library(
    LIBRARY_NAME
        icaruscode_CRT
    SOURCE
        CRTChannelMapAlg.cxx
    LIBRARIES larcorealg_Geometry
        icaruscode_CRTData
        cetlib_except::cetlib_except
        messagefacility::MF_MessageLogger
        
        fhiclcpp::fhiclcpp
        ${ROOT_BASIC_LIB_LIST}
        ROOT::Geom
        cetlib::cetlib
)

simple_plugin(CRTGeometryHelper service
              larcorealg_Geometry
              icaruscode_CRT
              art_root_io::TFileService_service
              art::Framework_Services_Registry
              messagefacility::MF_MessageLogger
              
              ${ROOT_BASIC_LIB_LIST}
)

simple_plugin(CRTDetSim module
              larcorealg_Geometry
              icaruscode_CRT
              icaruscode_CRTData
              sbnobj_Common_CRT
              icaruscode_CRTUtils
              art_root_io::TFileService_service
              lardataalg_DetectorInfo
              nurandom::RandomUtils_NuRandomService_service
              art::Framework_Services_Registry
              art::Framework_Services_Optional_RandomNumberGenerator_service
              messagefacility::MF_MessageLogger
              
              ${ROOT_BASIC_LIB_LIST}
              CLHEP::CLHEP
)

simple_plugin(CRTSimHitProducer module
              larcorealg_Geometry
              icaruscode_CRT
              icaruscode_CRTData
              sbnobj_Common_CRT
              icaruscode_CRTUtils
              art_root_io::TFileService_service
              lardataalg_DetectorInfo
              nurandom::RandomUtils_NuRandomService_service
              art::Framework_Services_Registry
              art::Framework_Services_Optional_RandomNumberGenerator_service
              messagefacility::MF_MessageLogger
              
              ${ROOT_BASIC_LIB_LIST}
              CLHEP::CLHEP
)

simple_plugin( CRTTrueHitProducer module
        larcorealg_Geometry
        icaruscode_CRT
        icaruscode_CRTData
        sbnobj_Common_CRT
        icaruscode_CRTUtils
        art::Framework_Core
        art::Framework_IO_Sources
        art::Framework_Principal
        art_root_io::TFileService_service
        art::Persistency_Provenance
        canvas::canvas
        cetlib_except::cetlib_except
        lardata_DetectorInfoServices_DetectorClocksServiceStandard_service
        art::Framework_Services_Registry
        messagefacility::MF_MessageLogger
        
        ${ROOT_BASIC_LIB_LIST}
        CLHEP::CLHEP
        cetlib::cetlib
)

simple_plugin(CRTTzeroProducer module
              larcorealg_Geometry
              icaruscode_CRT
              icaruscode_CRTData
              sbnobj_Common_CRT
              art_root_io::TFileService_service
              lardataalg_DetectorInfo
              nurandom::RandomUtils_NuRandomService_service
              art::Framework_Services_Registry
              art::Framework_Services_Optional_RandomNumberGenerator_service
              messagefacility::MF_MessageLogger
              
              ${ROOT_BASIC_LIB_LIST}
              CLHEP::CLHEP
)

simple_plugin(CRTTrackProducer module
              larcorealg_Geometry
              icaruscode_CRT
              icaruscode_CRTData
              sbnobj_Common_CRT
              icaruscode_CRTUtils
              art_root_io::TFileService_service
              lardataalg_DetectorInfo
              nurandom::RandomUtils_NuRandomService_service
              art::Framework_Services_Registry
              art::Framework_Services_Optional_RandomNumberGenerator_service
              messagefacility::MF_MessageLogger
              
              ${ROOT_BASIC_LIB_LIST}
              CLHEP::CLHEP
)

simple_plugin(CRTSimAnalysis module
              icaruscode_CRTData
              icaruscode_CRT
              sbnobj_Common_CRT
              icaruscode_CRTUtils
              larcore_Geometry_Geometry_service
              larcorealg_Geometry
              nusimdata::SimulationBase
              art::Framework_Services_Registry
              art_root_io::tfile_support ROOT::Core
              art_root_io::TFileService_service
              messagefacility::MF_MessageLogger
              
              ${ROOT_BASIC_LIB_LIST}
              ROOT::Geom
              ROOT::XMLIO
              ROOT::Gdml
)

simple_plugin(CRTDataAnalysis module
              icaruscode_CRTData
              icaruscode_CRT
              sbnobj_Common_CRT
              icaruscode_CRTUtils
              larcore_Geometry_Geometry_service
              larcorealg_Geometry
              nusimdata::SimulationBase
              art::Framework_Services_Registry
              art_root_io::tfile_support ROOT::Core
              art_root_io::TFileService_service
              messagefacility::MF_MessageLogger
              
              ${ROOT_BASIC_LIB_LIST}
              ROOT::Geom
              ROOT::XMLIO
              ROOT::Gdml
)

simple_plugin(CrtOpHitMatchAnalysis module
              icaruscode_CRTData
              icaruscode_CRT
              sbnobj_Common_CRT
              icaruscode_CRTUtils
              larcore_Geometry_Geometry_service
              larcorealg_Geometry
              nusimdata::SimulationBase
              lardataobj_RecoBase
              art::Framework_Services_Registry
              art_root_io::tfile_support ROOT::Core
              art_root_io::TFileService_service
              lardataalg_DetectorInfo
              messagefacility::MF_MessageLogger
              
              ${ROOT_BASIC_LIB_LIST}
)

simple_plugin(CRTTruthMatchAnalysis module
              icaruscode_CRTData
              icaruscode_CRT
              sbnobj_Common_CRT
              icaruscode_CRTUtils
              larcore_Geometry_Geometry_service
              larcorealg_Geometry
              nusimdata::SimulationBase
              art::Framework_Services_Registry
              art_root_io::tfile_support ROOT::Core
              art_root_io::TFileService_service
              messagefacility::MF_MessageLogger
              
              ${ROOT_BASIC_LIB_LIST}
              ROOT::Geom
              ROOT::XMLIO
              ROOT::Gdml
)

simple_plugin( CRTAutoVeto module
                        icaruscode_CRTData
                        icaruscode_CRT
                        sbnobj_Common_CRT
                        icaruscode_CRTUtils
                        larcore_Geometry_Geometry_service
                        larcorealg_Geometry
                        nusimdata::SimulationBase
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
                        
                        cetlib::cetlib cetlib_except
                        ${ROOT_BASIC_LIB_LIST}
                        ROOT::Geom
                        ROOT::XMLIO
                        ROOT::Gdml              
)

simple_plugin( FlashResAna module
                        larcore_Geometry_Geometry_service
                        larcorealg_Geometry
                        nusimdata::SimulationBase
                        lardataobj_RecoBase
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
                        
                        cetlib::cetlib cetlib_except
                        ${ROOT_BASIC_LIB_LIST}
                        ROOT::Geom
                        ROOT::XMLIO
                        ROOT::Gdml              
)

simple_plugin( PhotBackground module
                        larcore_Geometry_Geometry_service
                        larcorealg_Geometry
                        nusimdata::SimulationBase
                        icaruscode_CRTData
                        icaruscode_CRT
                        sbnobj_Common_CRT
                        icaruscode_CRTUtils
                        art::Framework_Core
                        art::Framework_Principal
                        art::Framework_Services_Registry
                        art_root_io::tfile_support ROOT::Core
                        art_root_io::TFileService_service
                        art::Persistency_Common
                        art::Persistency_Provenance
                        art::Utilities
                        ${ROOT_BASIC_LIB_LIST}
)

install_headers()
install_fhicl()
install_source()
