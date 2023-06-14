
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
        CRTT0Matching_module.cc
        CRTTPCMatchingAna_module.cc	
        CRTPMTMatchingAna_module.cc
        CRTTPCTruthEff_module.cc
    LIBRARY_NAME
        icaruscode_CRTData
    SOURCE
        CRTGeoObjectSorter.cxx
    LIBRARIES
        larcorealg::Geometry
        larcore::Geometry_Geometry_service
        canvas::canvas
        cetlib_except::cetlib_except
        messagefacility::MF_MessageLogger
        messagefacility::headers
        fhiclcpp::fhiclcpp
        art::Framework_Services_Registry
        CLHEP::CLHEP
        ROOT::Geom
        cetlib::cetlib
)

art_make_library(
        LIBRARY_NAME
                icaruscode_CRT
        SOURCE
                CRTChannelMapAlg.cxx
        LIBRARIES
                larcorealg::Geometry
                larcore::Geometry_Geometry_service
                sbnobj::ICARUS_CRT
                icaruscode_CRTData
                cetlib_except::cetlib_except
                messagefacility::MF_MessageLogger
                messagefacility::headers
                fhiclcpp::fhiclcpp
                art::Framework_Services_Registry
                CLHEP::CLHEP
                ROOT::Geom
                cetlib::cetlib
        )

cet_build_plugin(CRTGeometryHelper art::service
     LIBRARIES
        larcorealg::Geometry
        larcore::Geometry_Geometry_service
        icaruscode_CRT
        art_root_io::TFileService_service
        art::Framework_Services_Registry
        messagefacility::MF_MessageLogger
        messagefacility::headers
        )

cet_build_plugin(CRTDetSim art::module
     LIBRARIES
        larcorealg::Geometry
        icaruscode_CRT
        sbnobj::ICARUS_CRT
        sbnobj::Common_CRT
        icaruscode_CRTUtils
        art_root_io::TFileService_service
        lardataalg::DetectorInfo
        nurandom::RandomUtils_NuRandomService_service
        art::Framework_Services_Registry
        art::Framework_Services_Optional_RandomNumberGenerator_service
        messagefacility::MF_MessageLogger
        messagefacility::headers
        CLHEP::CLHEP
        )

cet_build_plugin(CRTSimHitProducer art::module
     LIBRARIES
        larcorealg::Geometry
        icaruscode_CRT
        sbnobj::ICARUS_CRT
        sbnobj::Common_CRT
        icaruscode_CRTUtils
        art_root_io::TFileService_service
        lardataalg::DetectorInfo
        nurandom::RandomUtils_NuRandomService_service
        art::Framework_Services_Registry
        art::Framework_Services_Optional_RandomNumberGenerator_service
        messagefacility::MF_MessageLogger
        messagefacility::headers
        CLHEP::CLHEP
        )

cet_build_plugin( CRTTrueHitProducer art::module
     LIBRARIES
        larcorealg::Geometry
        icaruscode_CRT
        sbnobj::ICARUS_CRT
        sbnobj::Common_CRT
        icaruscode_CRTUtils
        art::Framework_Core
        art::Framework_IO_Sources
        art::Framework_Principal
        art_root_io::TFileService_service
        art::Persistency_Provenance
        canvas::canvas
        cetlib_except::cetlib_except
        lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
        art::Framework_Services_Registry
        messagefacility::MF_MessageLogger
        messagefacility::headers
        CLHEP::CLHEP
        cetlib::cetlib
        )

cet_build_plugin(CRTTzeroProducer art::module
     LIBRARIES
        larcorealg::Geometry
        icaruscode_CRT
        sbnobj::ICARUS_CRT
        sbnobj::Common_CRT
        art_root_io::TFileService_service
        lardataalg::DetectorInfo
        nurandom::RandomUtils_NuRandomService_service
        art::Framework_Services_Registry
        art::Framework_Services_Optional_RandomNumberGenerator_service
        messagefacility::MF_MessageLogger
        messagefacility::headers
        CLHEP::CLHEP
        lardata::Utilities
        )

cet_build_plugin(CRTTrackProducer art::module
     LIBRARIES
        larcorealg::Geometry
        icaruscode_CRT
        sbnobj::ICARUS_CRT
        sbnobj::Common_CRT
        icaruscode_CRTUtils
        art_root_io::TFileService_service
        ROOT::Tree
        lardataalg::DetectorInfo
        nurandom::RandomUtils_NuRandomService_service
        art::Framework_Services_Registry
        art::Framework_Services_Optional_RandomNumberGenerator_service
        messagefacility::MF_MessageLogger
        messagefacility::headers
        CLHEP::CLHEP
        )

cet_build_plugin(CRTSimAnalysis art::module
     LIBRARIES
        sbnobj::ICARUS_CRT
        icaruscode_CRT
        sbnobj::Common_CRT
        icaruscode_CRTUtils
        larcore::Geometry_Geometry_service
        larcorealg::Geometry
        nusimdata::SimulationBase
        art::Framework_Services_Registry
        art_root_io::tfile_support ROOT::Core
        art_root_io::TFileService_service
        messagefacility::MF_MessageLogger
        messagefacility::headers
        ROOT::Core
        ROOT::Geom
        ROOT::XMLIO
        ROOT::Gdml
        ROOT::Tree
        )

cet_build_plugin(CRTDataAnalysis art::module
     LIBRARIES
        artdaq_core::artdaq-core_Data
        sbnobj::ICARUS_CRT
        icaruscode_CRT
        sbnobj::Common_CRT
        icaruscode_CRTUtils
        larcore::Geometry_Geometry_service
        larcorealg::Geometry
        nusimdata::SimulationBase
        art::Framework_Services_Registry
        art_root_io::tfile_support ROOT::Core
        art_root_io::TFileService_service
        messagefacility::MF_MessageLogger
        messagefacility::headers
        ROOT::Core
        ROOT::Tree
        ROOT::Geom
        ROOT::XMLIO
        ROOT::Gdml
        )

cet_build_plugin(CRTCalibrationAnalysis art::module
     LIBRARIES
        icaruscode_CRTData
        icaruscode_CRT
        sbnobj::Common_CRT
        icaruscode_CRTUtils
        larcore::Geometry_Geometry_service
        larcorealg::Geometry
        nusimdata::SimulationBase
        art::Framework_Services_Registry
        art_root_io::tfile_support ROOT::Core
        art_root_io::TFileService_service
        messagefacility::MF_MessageLogger
        messagefacility::headers
        ROOT::Core
        ROOT::Tree
        ROOT::Geom
        ROOT::XMLIO
        ROOT::Gdml
)


cet_build_plugin(CrtOpHitMatchAnalysis art::module
     LIBRARIES
        sbnobj::ICARUS_CRT
        icaruscode_CRT
        sbnobj::Common_CRT
        icaruscode_CRTUtils
        larcore::Geometry_Geometry_service
        larcorealg::Geometry
        nusimdata::SimulationBase
        lardataobj::RecoBase
        art::Framework_Services_Registry
        art_root_io::tfile_support ROOT::Core
        art_root_io::TFileService_service
        lardataalg::DetectorInfo
        messagefacility::MF_MessageLogger
        messagefacility::headers
        ROOT::Core
        ROOT::Tree
        )

cet_build_plugin(CRTPMTMatchingAna art::module
     LIBRARIES
        sbnobj::ICARUS_CRT
        icaruscode_CRT
        sbnobj::Common_CRT
        icaruscode_CRTUtils
        icaruscode::Decode_DataProducts
        larcore::Geometry_Geometry_service
        larcorealg::Geometry
        nusimdata::SimulationBase
        lardataobj::RecoBase
        lardataalg::DetectorInfo
        art::Framework_Services_Registry
        art_root_io::tfile_support ROOT::Core
        art_root_io::TFileService_service
        lardataalg::DetectorInfo
        messagefacility::MF_MessageLogger
        messagefacility::headers
        ROOT::Core
        ROOT::Tree
)

cet_build_plugin(CRTPMTMatchingProducer art::module
     LIBRARIES
        icaruscode::CRT
        sbnobj::Common_CRT
        icaruscode::CRTUtils
        larcore::Geometry_Geometry_service
        larcorealg::Geometry
        larcoreobj::headers
        lardataobj::RecoBase
        art::Framework_Services_Registry
        messagefacility::MF_MessageLogger
)

cet_build_plugin(CRTTruthMatchAnalysis art::module
     LIBRARIES
        sbnobj::ICARUS_CRT
        icaruscode_CRT
        sbnobj::Common_CRT
        icaruscode_CRTUtils
        larcore::Geometry_Geometry_service
        larcorealg::Geometry
        nusimdata::SimulationBase
        art::Framework_Services_Registry
        art_root_io::tfile_support
        art_root_io::TFileService_service
        messagefacility::MF_MessageLogger
        messagefacility::headers
        ROOT::Core
        ROOT::Tree
        ROOT::Geom
        ROOT::XMLIO
        ROOT::Gdml
        )

cet_build_plugin( CRTAutoVeto art::module
     LIBRARIES
        sbnobj::ICARUS_CRT
        icaruscode_CRT
        sbnobj::Common_CRT
        icaruscode_CRTUtils
        larcore::Geometry_Geometry_service
        larcorealg::Geometry
        nusimdata::SimulationBase
        art::Framework_Core
        art::Framework_Principal
        art::Framework_Services_Registry
        art_root_io::tfile_support
        art_root_io::TFileService_service
        art::Persistency_Common
        art::Persistency_Provenance
        art::Utilities
        canvas::canvas
        messagefacility::MF_MessageLogger
        messagefacility::headers
        cetlib::cetlib cetlib_except::cetlib_except
        ROOT::Core
        ROOT::Tree
        ROOT::Geom
        ROOT::XMLIO
        ROOT::Gdml
        )

cet_build_plugin( FlashResAna art::module
     LIBRARIES
        larcore::Geometry_Geometry_service
        larcorealg::Geometry
        nusimdata::SimulationBase
        lardataobj::RecoBase
        art::Framework_Core
        art::Framework_Principal
        art::Framework_Services_Registry
        art_root_io::tfile_support
        art_root_io::TFileService_service
        art::Persistency_Common
        art::Persistency_Provenance
        art::Utilities
        canvas::canvas
        messagefacility::MF_MessageLogger
        messagefacility::headers
        cetlib::cetlib cetlib_except::cetlib_except
        ROOT::Tree
        ROOT::Core
        ROOT::Geom
        ROOT::XMLIO
        ROOT::Gdml
        )

cet_build_plugin( PhotBackground art::module
     LIBRARIES
        larcore::Geometry_Geometry_service
        larcorealg::Geometry
        nusimdata::SimulationBase
        sbnobj::ICARUS_CRT
        icaruscode_CRT
        sbnobj::Common_CRT
        icaruscode_CRTUtils
        art::Framework_Core
        art::Framework_Principal
        art::Framework_Services_Registry
        art_root_io::tfile_support
        art_root_io::TFileService_service
        art::Persistency_Common
        art::Persistency_Provenance
        art::Utilities
        ROOT::Core
        ROOT::Geom
        ROOT::XMLIO
        ROOT::Gdml
        ROOT::Tree
        )

cet_build_plugin(CRTT0Matching art::module
     LIBRARIES
              icaruscode_CRTData
              icaruscode_CRT
              sbnobj::Common_CRT
              icaruscode_CRTUtils
              larcore::Geometry_Geometry_service
              larsim::Simulation lardataobj::Simulation
              larsim::MCCheater_BackTrackerService_service
              larsim::MCCheater_ParticleInventoryService_service
              lardata::Utilities
              larevt::Filters
              lardataobj::RawData
              lardataobj::RecoBase
              lardataobj::AnalysisBase
              lardata::RecoObjects
              larpandora::LArPandoraInterface
              larcorealg::Geometry
              nusimdata::SimulationBase
              art::Framework_Core
              art::Framework_Principal
              art::Framework_Services_Registry
              art_root_io::tfile_support
              art_root_io::TFileService_service
              art::Persistency_Common canvas::canvas
              art::Persistency_Provenance
              art::Utilities
              messagefacility::MF_MessageLogger
              ROOT::Core
              ROOT::Geom
              ROOT::XMLIO
              ROOT::Gdml
              ROOT::Tree
              ROOT::Spectrum
              ROOT::RooFit
              ROOT::RooFitCore
)


cet_build_plugin(CRTT0MatchingAna art::module
     LIBRARIES
              sbnobj::ICARUS_CRT
              icaruscode_CRT
              sbnobj::Common_CRT
              icaruscode_CRTUtils
              larcorealg::Geometry
              larcore::Geometry_Geometry_service
              larsim::Simulation lardataobj::Simulation
              larsim::MCCheater_BackTrackerService_service
              larsim::MCCheater_ParticleInventoryService_service
              lardata::Utilities
              larevt::Filters
              lardataobj::RawData
              lardataobj::RecoBase
              lardataobj::AnalysisBase
              lardata::RecoObjects
              larpandora::LArPandoraInterface
              larcorealg::Geometry
              nusimdata::SimulationBase
              art::Persistency_Common canvas::canvas
              art::Persistency_Provenance
              art::Utilities
              ROOT::Core
              ROOT::Tree
              ROOT::Geom
              ROOT::XMLIO
              ROOT::Gdml
              ROOT::Spectrum
              ROOT::RooFit
              ROOT::RooFitCore
)

cet_build_plugin(CRTTPCMatchingAna art::module
     LIBRARIES
              sbnobj::ICARUS_CRT
              icaruscode_CRT
              sbnobj::Common_CRT
              icaruscode_CRTUtils
              larcorealg::Geometry
              larcore::Geometry_Geometry_service
              larsim::Simulation lardataobj::Simulation
              larsim::MCCheater_BackTrackerService_service
              larsim::MCCheater_ParticleInventoryService_service
              lardata::Utilities
              larevt::Filters
              lardataobj::RawData
              lardataobj::RecoBase
              lardataobj::AnalysisBase
              lardata::RecoObjects
              larpandora::LArPandoraInterface
              larcorealg::Geometry
              nusimdata::SimulationBase
              art::Persistency_Common canvas::canvas
              art::Persistency_Provenance
              art::Utilities
              ROOT::Core
              ROOT::Geom
              ROOT::XMLIO
              ROOT::Gdml
              ROOT::Spectrum
              ROOT::RooFit
              ROOT::RooFitCore
)
simple_plugin(CRTTPCTruthEff module
              sbnobj::ICARUS_CRT
              icaruscode_CRT
              sbnobj::Common_CRT
              icaruscode_CRTUtils
              larcorealg::Geometry
              larcore::Geometry_Geometry_service
              larsim::Simulation lardataobj::Simulation
              larsim::MCCheater_BackTrackerService_service
              larsim::MCCheater_ParticleInventoryService_service
              lardata::Utilities
              larevt::Filters
              lardataobj::RawData
              lardataobj::RecoBase
              lardataobj::AnalysisBase
              lardata::RecoObjects
              larpandora::LArPandoraInterface
              larcorealg::Geometry
              nusimdata::SimulationBase
              art::Persistency_Common canvas::canvas
              art::Persistency_Provenance
              art::Utilities
              ROOT::Core
              ROOT::Geom
              ROOT::XMLIO
              ROOT::Gdml
              ROOT::Spectrum
              ROOT::RooFit
              ROOT::RooFitCore

)


install_headers()
install_fhicl()
install_source()
