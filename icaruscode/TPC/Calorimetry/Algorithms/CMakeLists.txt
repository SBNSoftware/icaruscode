include_directories(${FFTW3_INCLUDE_DIRS})
link_directories(${FFTW3_LIBRARY_DIRS})

art_make(
	LIB_LIBRARIES
		lardataobj::RawData
		lardataobj::RecoBase
		lardata::Utilities
		fhiclcpp::fhiclcpp
		cetlib::cetlib cetlib_except::cetlib_except
		larcorealg::Geometry
		larcore::Geometry_Geometry_service
		lardata::Utilities
		larevt::Filters
		lardataobj::RawData
		larevt::CalibrationDBI_IOVData
		larevt::CalibrationDBI_Providers
		lardataalg::DetectorInfo
		lardataobj::RecoBase
		icaruscode::TPC_Utilities_SignalShapingICARUSService_service
		${ICARUS_FFTW_LIBRARIES}
		art::Framework_Core
		art::Framework_Principal
		art::Framework_Services_Registry
		art_root_io::tfile_support ROOT::Core
		art_root_io::TFileService_service
		art::Persistency_Common canvas::canvas
		art::Persistency_Provenance canvas::canvas
		art::Utilities canvas::canvas
		messagefacility::MF_MessageLogger
		messagefacility::headers
		fhiclcpp::fhiclcpp
		cetlib::cetlib cetlib_except::cetlib_except
		ROOT::Geom
		ROOT::XMLIO
		ROOT::Gdml
		ROOT::FFTW

	TOOL_LIBRARIES
		icaruscode::TPC_Calorimetry_Algorithms
		icaruscode::TPC_Utilities_SignalShapingICARUSService_service
		lardataobj::RawData
		lardataobj::RecoBase
		lardata::Utilities
		fhiclcpp::fhiclcpp
		cetlib::cetlib cetlib_except::cetlib_except
		larcorealg::Geometry
		larcore::Geometry_Geometry_service
		larevt::CalibrationDBI_IOVData
		larevt::CalibrationDBI_Providers
		${ICARUS_FFTW_LIBRARIES}
		art::Framework_Core
		art::Framework_Principal
		art::Framework_Services_Registry
		art_root_io::tfile_support ROOT::Core
		art_root_io::TFileService_service
		art::Persistency_Common canvas::canvas
		art::Persistency_Provenance canvas::canvas
		art::Utilities canvas::canvas
		messagefacility::MF_MessageLogger
		messagefacility::headers
		fhiclcpp::fhiclcpp
		cetlib::cetlib cetlib_except::cetlib_except
		ROOT::Geom
		ROOT::FFTW
	)

install_headers()
install_fhicl()
install_source()
