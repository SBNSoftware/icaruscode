art_make(
 EXCLUDE
        CrtCal.cc
        CrtCalTree.cxx
	CRTTiming.cc
	CRTRawTree.cc
	CRTPreProcessTree.cc
	CRTMergePreProcessTrees.cxx
 NO_PLUGINS # WAnt different link lists for modues
 BASENAME_ONLY
 LIBRARY_NAME sbndaq-artdaq_ArtModules_Common
 LIB_LIBRARIES sbndaq-artdaq-core_Overlays
               sbndaq-artdaq-core_Overlays_Common
)

art_make_library(
    LIBRARY_NAME
        CRT_CAL
    SOURCE
        CrtCal.cc
    LIBRARIES 
        cetlib_except
        ${ROOT_BASIC_LIB_LIST}
        ${ROOT_SPECTRUM}
        ${ROOT_GEOM}
        ${ROOT_CORE}
        ${CETLIB}
)

art_make_library(
    LIBRARY_NAME
        CRT_CAL_TREE
    SOURCE
        CrtCalTree.cxx
    LIBRARIES 
        cetlib_except
        ${ROOT_CORE}
        ${ROOT_BASIC_LIB_LIST}
        ${ROOT_SPECTRUM}
        ${ROOT_GEOM}
        ${CETLIB}
)

art_make_library(
    LIBRARY_NAME
        CRT_RAW_TREE
    SOURCE
        CRTRawTree.cc
    LIBRARIES 
        cetlib_except
        ${ROOT_CORE}
        ${ROOT_BASIC_LIB_LIST}
        ${ROOT_SPECTRUM}
        ${ROOT_GEOM}
        ${CETLIB}
)

art_make_library(
    LIBRARY_NAME
        CRT_TIMING
    SOURCE
        CRTTiming.cc
    LIBRARIES 
        cetlib_except
        ${ROOT_CORE}
        ${ROOT_BASIC_LIB_LIST}
        ${ROOT_SPECTRUM}
        ${ROOT_GEOM}
        ${CETLIB}
	CRT_RAW_TREE
	CRT_PREPROCESS_TREE
)

art_make_library(
    LIBRARY_NAME
        CRT_PREPROCESS_TREE
    SOURCE
        CRTPreProcessTree.cc
    LIBRARIES 
        cetlib_except
        ${ROOT_CORE}
        ${ROOT_BASIC_LIB_LIST}
        ${ROOT_SPECTRUM}
        ${ROOT_GEOM}
        ${CETLIB}
)

art_make_library(
    LIBRARY_NAME
        CRT_MERGE_PREPROCESS_TREE
    SOURCE
        CRTMergePreProcessTrees.cxx
    LIBRARIES 
        cetlib_except
        ${ROOT_CORE}
        ${ROOT_BASIC_LIB_LIST}
        ${ROOT_SPECTRUM}
        ${ROOT_GEOM}
        ${CETLIB}
	CRT_TIMING
	CRT_PREPROCESS_TREE
)


install_headers()
install_source()

simple_plugin( CrtNoiseMonTool "module"
  art_Framework_Services_Registry
  art_root_io_tfile_support
  art_root_io_TFileService_service
  #art_Framework_Services_Optional
  #art_Framework_Services_Optional_TFileService_service
  ${ROOT_BASIC_LIB_LIST}
        ${ARTDAQ-CORE_UTILITIES}
  sbndaq-artdaq-core_Overlays_Common 
)

simple_plugin( CrtCalAnalysis "module"
  art_Framework_Services_Registry
  art_root_io_tfile_support
  art_root_io_TFileService_service
  #art_Framework_Services_Optional
  #art_Framework_Services_Optional_TFileService_service
  ${ROOT_BASIC_LIB_LIST}
        ${ARTDAQ-CORE_UTILITIES}
  sbndaq-artdaq-core_Overlays_Common 
  CRT_CAL
)

simple_plugin( AnaProducer "module"
  art_Framework_Services_Registry
  art_root_io_tfile_support
  art_root_io_TFileService_service
  #art_Framework_Services_Optional
  #art_Framework_Services_Optional_TFileService_service
  ${ROOT_BASIC_LIB_LIST}
        ${ARTDAQ-CORE_UTILITIES}
  sbndaq-artdaq-core_Overlays_Common 
  CRT_CAL
  CRT_CAL_TREE
  CRT_RAW_TREE
  CRT_TIMING
)

simple_plugin( CRTEventProducer "module"
  art_Framework_Services_Registry
  art_root_io_tfile_support
  art_root_io_TFileService_service
  #art_Framework_Services_Optional
  #art_Framework_Services_Optional_TFileService_service
  #${ROOT_BASIC_LIB_LIST}
  CRT_PREPROCESS_TREE
  #CRT_CAL_TREE
)

install_fhicl()
