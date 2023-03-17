art_make_library(
	EXCLUDE
		CrtCal.cc
		CrtCalTree.cxx
		CRTTiming.cc
		CRTRawTree.cc
		CRTPreProcessTree.cc
		CRTMergePreProcessTrees.cxx
	
	NO_PLUGINS
	LIBRARY_NAME sbndaq-artdaq_ArtModules_Common
	LIBRARIES sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays
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
  artdaq_core::artdaq-core_Utilities
  sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_Common
  artdaq_core::artdaq-core_Data
)

cet_build_plugin( CrtNoiseMonTool art::module
  LIBRARIES
  art::Framework_Services_Registry
  art_root_io::tfile_support
  art_root_io::TFileService_service
  ROOT::Core ROOT::Tree ROOT::Hist
  artdaq_core::artdaq-core_Utilities
  sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_Common
  artdaq_core::artdaq-core_Data
)

cet_build_plugin( CrtCalAnalysis art::module
  LIBRARIES
  art::Framework_Services_Registry
  art_root_io::tfile_support
  art_root_io::TFileService_service
  ROOT::Core ROOT::Tree ROOT::Hist
  artdaq_core::artdaq-core_Utilities
  sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_Common
  artdaq_core::artdaq-core_Data
  CRT_CAL
)

cet_build_plugin( AnaProducer art::module
  LIBRARIES
  art::Framework_Services_Registry
  art_root_io::tfile_support
  art_root_io::TFileService_service
  ROOT::Core ROOT::Tree ROOT::Hist
  artdaq_core::artdaq-core_Utilities
  sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_Common
  artdaq_core::artdaq-core_Data
  CRT_CAL
  CRT_CAL_TREE
  CRT_RAW_TREE
  CRT_TIMING
)

cet_build_plugin( CRTEventProducer art::module
  LIBRARIES
  art::Framework_Services_Registry
  art_root_io::tfile_support
  art_root_io::TFileService_service
  CRT_PREPROCESS_TREE
  #CRT_CAL_TREE
  sbnobj::Common_CRT
  lardata::headers
  larcore::headers
  larcorealg::headers
  larcoreobj::headers
  lardataalg::headers
  lardataobj::RawData
)

cet_build_plugin( DecoderICARUSCRT art::module
  LIBRARIES
  icaruscode_Utilities
  art::Framework_Services_Registry
  art_root_io::tfile_support
  art_root_io::TFileService_service
  ROOT::Core ROOT::Tree ROOT::Hist
  artdaq_core::artdaq-core_Utilities
  sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_Common
  artdaq_core::artdaq-core_Data
  sbnobj::Common_CRT
)
cet_build_plugin( DecoderICARUSCRTNew art::module
  LIBRARIES
  icaruscode_Utilities
  art::Framework_Services_Registry
  art_root_io::tfile_support
  art_root_io::TFileService_service
  ROOT::Core ROOT::Tree ROOT::Hist
  artdaq_core::artdaq-core_Utilities
  sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_Common
  artdaq_core::artdaq-core_Data
  sbnobj::Common_CRT
)

install_fhicl()
