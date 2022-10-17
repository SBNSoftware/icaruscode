include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories($ENV{NURANDOM_INC})
include_directories($ENV{IFDHX_FQ_DIR}/inc)
include_directories(${TRACE_INC})

add_subdirectory(CRTUtils)
add_subdirectory(CRTDecoder)

art_make_library(
	LIBRARY_NAME
		icaruscode_CRT
	SOURCE
		CRTChannelMapAlg.cxx
                CRTGeoObjectSorter.cxx
	LIBRARIES larcorealg::Geometry
		sbnobj::ICARUS_CRT
		cetlib_except::cetlib_except
		messagefacility::MF_MessageLogger
		messagefacility::headers
		fhiclcpp::fhiclcpp
	        art::Framework_Services_Registry
	        CLHEP::CLHEP	
		ROOT::Geom
		cetlib::cetlib
	)

simple_plugin(CRTGeometryHelper service
	larcorealg::Geometry
	larcore::Geometry_Geometry_service
	icaruscode::CRT
	art_root_io::TFileService_service
	art::Framework_Services_Registry
	messagefacility::MF_MessageLogger
	messagefacility::headers
	
	)

simple_plugin(CRTDetSim module
	larcorealg::Geometry
	icaruscode::CRT
	sbnobj::ICARUS_CRT
	sbnobj::Common_CRT
	icaruscode::CRTUtils
	art_root_io::TFileService_service
	lardataalg::DetectorInfo
	nurandom::RandomUtils_NuRandomService_service
	art::Framework_Services_Registry
	art::Framework_Services_Optional_RandomNumberGenerator_service
	messagefacility::MF_MessageLogger
	messagefacility::headers
	CLHEP::CLHEP
	)

simple_plugin(CRTSimHitProducer module
	larcorealg::Geometry
	icaruscode::CRT
	sbnobj::ICARUS_CRT
	sbnobj::Common_CRT
	icaruscode::CRTUtils
	art_root_io::TFileService_service
	lardataalg::DetectorInfo
	nurandom::RandomUtils_NuRandomService_service
	art::Framework_Services_Registry
	art::Framework_Services_Optional_RandomNumberGenerator_service
	messagefacility::MF_MessageLogger
	messagefacility::headers
	CLHEP::CLHEP
	)

simple_plugin(CRTTrueHitProducer module
	larcorealg::Geometry
	icaruscode::CRT
	sbnobj::ICARUS_CRT
	sbnobj::Common_CRT
	icaruscode::CRTUtils
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

simple_plugin(CRTTzeroProducer module
	larcorealg::Geometry
	icaruscode::CRT
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

simple_plugin(CRTTrackProducer module
	larcorealg::Geometry
	icaruscode::CRT
	sbnobj::ICARUS_CRT
	sbnobj::Common_CRT
	icaruscode::CRTUtils
	art_root_io::TFileService_service
	lardataalg::DetectorInfo
	nurandom::RandomUtils_NuRandomService_service
	art::Framework_Services_Registry
	art::Framework_Services_Optional_RandomNumberGenerator_service
	messagefacility::MF_MessageLogger
	messagefacility::headers
	
	CLHEP::CLHEP
	)

simple_plugin(CRTSimAnalysis module
	sbnobj::ICARUS_CRT
	icaruscode::CRT
	sbnobj::Common_CRT
	icaruscode::CRTUtils
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
	)

simple_plugin(CRTDataAnalysis module
        artdaq_core::artdaq-core_Data
	sbnobj::ICARUS_CRT
	icaruscode::CRT
	sbnobj::Common_CRT
	icaruscode::CRTUtils
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
	)

simple_plugin(CrtOpHitMatchAnalysis module
	sbnobj::ICARUS_CRT
	icaruscode::CRT
	sbnobj::Common_CRT
	icaruscode::CRTUtils
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
	
	)

simple_plugin(CRTPMTMatchingAna module
              sbnobj::ICARUS_CRT
	      icaruscode::CRT
              sbnobj::Common_CRT
              icaruscode::CRTUtils
	      larcore::Geometry_Geometry_service
              larcorealg::Geometry
              nusimdata::SimulationBase
              lardataobj::RecoBase
              ${ART_FRAMEWORK_SERVICES_REGISTRY}
              ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
              ${ART_ROOT_IO_TFILESERVICE_SERVICE}
              lardataalg::DetectorInfo
              ${MF_MESSAGELOGGER}
              ${MF_UTILITIES}
              ${ROOT_BASIC_LIB_LIST}
)

simple_plugin(CRTPMTMatchingAna module
              icaruscode_CRTData
	      icaruscode_CRT
              sbnobj_Common_CRT
              icaruscode_CRTUtils
	      larcore_Geometry_Geometry_service
              larcorealg_Geometry
              nusimdata_SimulationBase
              lardataobj_RecoBase
              ${ART_FRAMEWORK_SERVICES_REGISTRY}
              ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
              ${ART_ROOT_IO_TFILESERVICE_SERVICE}
              lardataalg_DetectorInfo
              ${MF_MESSAGELOGGER}
              ${MF_UTILITIES}
              ${ROOT_BASIC_LIB_LIST}
)

simple_plugin(CRTTruthMatchAnalysis module
	sbnobj::ICARUS_CRT
	icaruscode::CRT
	sbnobj::Common_CRT
	icaruscode::CRTUtils
	larcore::Geometry_Geometry_service
	larcorealg::Geometry
	nusimdata::SimulationBase
	art::Framework_Services_Registry
	art_root_io::tfile_support
	art_root_io::TFileService_service
	messagefacility::MF_MessageLogger
	messagefacility::headers
	ROOT::Core
	ROOT::Geom
	ROOT::XMLIO
	ROOT::Gdml
	)

simple_plugin(CRTAutoVeto module
	sbnobj::ICARUS_CRT
	icaruscode::CRT
	sbnobj::Common_CRT
	icaruscode::CRTUtils
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
	ROOT::Geom
	ROOT::XMLIO
	ROOT::Gdml
	)

simple_plugin(FlashResAna module
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

simple_plugin(PhotBackground module
	larcore::Geometry_Geometry_service
	larcorealg::Geometry
	nusimdata::SimulationBase
	sbnobj::ICARUS_CRT
	icaruscode::CRT
	sbnobj::Common_CRT
	icaruscode::CRTUtils
	art::Framework_Core
	art::Framework_Principal
	art::Framework_Services_Registry
	art_root_io::tfile_support
	art_root_io::TFileService_service
	art::Persistency_Common
	art::Persistency_Provenance
	art::Utilities
	ROOT::Core
	)

simple_plugin(CRTT0Matching module
              sbnobj::ICARUS_CRT
              icaruscode::CRT
              sbnobj::Common_CRT
              icaruscode::CRTUtils
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
              ${ART_FRAMEWORK_SERVICES_REGISTRY}
              ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
              ${ART_ROOT_IO_TFILESERVICE_SERVICE}
              art::Persistency_Common canvas::canvas
              art::Persistency_Provenance
              art::Utilities
              ${MF_MESSAGELOGGER}
              ${MF_UTILITIES}
              ${ROOT_BASIC_LIB_LIST}
              ${ROOT_GEOM}
              ${ROOT_XMLIO}
              ${ROOT_GDML}
              ${ROOT_SPECTRUM}
              ${ROOT_ROOFIT}
              ${ROOT_ROOFITCORE}
)

simple_plugin(CRTT0MatchingAna module
              sbnobj::ICARUS_CRT
              icaruscode::CRT
              sbnobj::Common_CRT
              icaruscode::CRTUtils
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
              ${ART_FRAMEWORK_SERVICES_REGISTRY}
              ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
              ${ART_ROOT_IO_TFILESERVICE_SERVICE}
              art::Persistency_Common canvas::canvas
              art::Persistency_Provenance
              art::Utilities
              ${MF_MESSAGELOGGER}
              ${MF_UTILITIES}
              ${ROOT_BASIC_LIB_LIST}
              ${ROOT_GEOM}
              ${ROOT_XMLIO}
              ${ROOT_GDML}
              ${ROOT_SPECTRUM}
              ${ROOT_ROOFIT}
              ${ROOT_ROOFITCORE}
)

simple_plugin(CRTTPCMatchingAna module
	      sbnobj::ICARUS_CRT
              icaruscode::CRT
              sbnobj::Common_CRT
              icaruscode::CRTUtils
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
              ${ART_FRAMEWORK_SERVICES_REGISTRY}
              ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
              ${ART_ROOT_IO_TFILESERVICE_SERVICE}
              art::Persistency_Common canvas::canvas
              art::Persistency_Provenance
              art::Utilities
              ${MF_MESSAGELOGGER}
              ${MF_UTILITIES}
              ${ROOT_BASIC_LIB_LIST}
              ${ROOT_GEOM}
              ${ROOT_XMLIO}
              ${ROOT_GDML}
              ${ROOT_SPECTRUM}
              ${ROOT_ROOFIT}
              ${ROOT_ROOFITCORE}
)

simple_plugin(CRTTPCMatchingAna module
              icaruscode_CRTData
              icaruscode_CRT
              sbnobj_Common_CRT
              icaruscode_CRTUtils
              larcorealg_Geometry
              larcore_Geometry_Geometry_service
              larsim_Simulation lardataobj_Simulation
              larsim_MCCheater_BackTrackerService_service
              larsim_MCCheater_ParticleInventoryService_service
              lardata_Utilities
              larevt_Filters
              lardataobj_RawData
              lardataobj_RecoBase
              lardataobj_AnalysisBase
              lardata_RecoObjects
              larpandora_LArPandoraInterface
              larcorealg_Geometry
              nusimdata_SimulationBase
              ${ART_FRAMEWORK_SERVICES_REGISTRY}
              ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
              ${ART_ROOT_IO_TFILESERVICE_SERVICE}
              art_Persistency_Common canvas
              art_Persistency_Provenance canvas
              art_Utilities canvas
              ${MF_MESSAGELOGGER}
              ${MF_UTILITIES}
              ${ROOT_BASIC_LIB_LIST}
              ${ROOT_GEOM}
              ${ROOT_XMLIO}
              ${ROOT_GDML}
              ${ROOT_SPECTRUM}
              ${ROOT_ROOFIT}
              ${ROOT_ROOFITCORE}
)

install_headers()
install_fhicl()
install_source()
