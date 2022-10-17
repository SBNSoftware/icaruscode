include_directories(${TRACE_INC})

art_make(
	EXCLUDE
		CrtCal.cc
		CrtCalTree.cxx
	CRTTiming.cc
	CRTRawTree.cc
	CRTPreProcessTree.cc
	CRTMergePreProcessTrees.cxx
	
	NO_PLUGINS
	BASENAME_ONLY
	LIBRARY_NAME sbndaq-artdaq_ArtModules_Common
	LIB_LIBRARIES sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays
		sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_Common
                canvas::canvas
)

art_make_library(
	LIBRARY_NAME
		CRT_CAL
	SOURCE
		CrtCal.cc
	LIBRARIES
		cetlib_except::cetlib_except
                ROOT::Graf
		ROOT::Spectrum
		ROOT::Geom
		ROOT::Core
                ROOT::Gpad
		cetlib::cetlib
)

art_make_library(
	LIBRARY_NAME
		CRT_CAL_TREE
	SOURCE
		CrtCalTree.cxx
	LIBRARIES
		cetlib_except::cetlib_except
		ROOT::Core
		ROOT::Spectrum
		ROOT::Geom
		cetlib::cetlib
	)

art_make_library(
	LIBRARY_NAME
		CRT_RAW_TREE
	SOURCE
		CRTRawTree.cc
	LIBRARIES
		cetlib_except::cetlib_except
		ROOT::Core
		ROOT::Spectrum
		ROOT::Geom
		cetlib::cetlib
	)

art_make_library(
	LIBRARY_NAME
		CRT_TIMING
	SOURCE
		CRTTiming.cc
	LIBRARIES
		cetlib_except::cetlib_except
		ROOT::Core
		ROOT::Spectrum
		ROOT::Geom
		cetlib::cetlib ROOT::Tree
		CRT_RAW_TREE
		CRT_PREPROCESS_TREE
	)

art_make_library(
	LIBRARY_NAME
		CRT_PREPROCESS_TREE
	SOURCE
		CRTPreProcessTree.cc
	LIBRARIES
		cetlib_except::cetlib_except
		ROOT::Core
		ROOT::Spectrum
		ROOT::Geom
		cetlib::cetlib
	)

art_make_library(
	LIBRARY_NAME
		CRT_MERGE_PREPROCESS_TREE
	SOURCE
		CRTMergePreProcessTrees.cxx
	LIBRARIES
		cetlib_except::cetlib_except
		ROOT::Core
		ROOT::Spectrum
		ROOT::Geom
                ROOT::Tree
		cetlib::cetlib
		CRT_TIMING
		CRT_PREPROCESS_TREE
)

install_headers()
install_source()

cet_build_plugin( BernCRTAna art::module
  LIBRARIES
  art::Framework_Services_Registry
  art_root_io::tfile_support
  art_root_io::TFileService_service
  ROOT::Core ROOT::Tree
  canvas::canvas
        ${ARTDAQ-CORE_UTILITIES}
  sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_Common 
  artdaq_core::artdaq-core_Data
)

simple_plugin( CrtNoiseMonTool "module"
  art::Framework_Services_Registry
  art_root_io::tfile_support
  art_root_io::TFileService_service
  #art::Framework_Services_Optional
  #art::Framework_Services_Optional_TFileService_service
  ROOT::Core ROOT::Tree ROOT::Hist
  ${ROOT_BASIC_LIB_LIST}
        ${ARTDAQ-CORE_UTILITIES}
  sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_Common 
  artdaq_core::artdaq-core_Data
)

simple_plugin( CrtCalAnalysis "module"
  art::Framework_Services_Registry
  art_root_io::tfile_support
  art_root_io::TFileService_service
  #art::Framework_Services_Optional
  #art::Framework_Services_Optional_TFileService_service
  ROOT::Core ROOT::Tree ROOT::Hist
  ${ROOT_BASIC_LIB_LIST}
        ${ARTDAQ-CORE_UTILITIES}
  sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_Common 
  artdaq_core::artdaq-core_Data
  CRT_CAL
)

simple_plugin( AnaProducer "module"
  art::Framework_Services_Registry
  art_root_io::tfile_support
  art_root_io::TFileService_service
  #art::Framework_Services_Optional
  #art::Framework_Services_Optional_TFileService_service
  ROOT::Core ROOT::Tree ROOT::Hist
  ${ROOT_BASIC_LIB_LIST}
        ${ARTDAQ-CORE_UTILITIES}
  sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_Common 
  artdaq_core::artdaq-core_Data
  CRT_CAL
  CRT_CAL_TREE
  CRT_RAW_TREE
  CRT_TIMING
)

simple_plugin( CRTEventProducer "module"
  art::Framework_Services_Registry
  art_root_io::tfile_support
  art_root_io::TFileService_service
  #art::Framework_Services_Optional
  #art::Framework_Services_Optional_TFileService_service
  #${ROOT_BASIC_LIB_LIST}
  CRT_PREPROCESS_TREE
  #CRT_CAL_TREE
  sbnobj::Common_CRT
  lardata::headers
  larcore::headers
  larcorealg::headers
  larcoreobj::headers
  lardataalg::headers
)

simple_plugin( DecoderICARUSCRT "module"
  icaruscode_Utilities
  art::Framework_Services_Registry
  art_root_io::tfile_support
  art_root_io::TFileService_service
  ROOT::Core ROOT::Tree ROOT::Hist
  ${ROOT_BASIC_LIB_LIST}
  ${ARTDAQ-CORE_UTILITIES}
  sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_Common
  artdaq_core::artdaq-core_Data
  sbnobj::Common_CRT
)

install_fhicl()
