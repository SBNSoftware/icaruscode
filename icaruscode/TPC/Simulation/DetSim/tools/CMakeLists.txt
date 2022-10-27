cet_enable_asserts()
include_directories(${FFTW3_INCLUDE_DIRS})

add_definitions(-DEIGEN_FFTW_DEFAULT)

art_make_library(
        SOURCE
          CoherentNoiseFactorProvider.cxx
        LIBRARIES
          art::Framework_Services_Registry
          larcore::Geometry_Geometry_service
          messagefacility::MF_MessageLogger
)

foreach ( TOOL
CorrelatedNoise
NoiseFromHist
NoNoise
RandomNoise
SBNDataNoiseBoard
SBNDataNoise
SBNNoise
)
  cet_build_plugin( ${TOOL} art::tool LIBRARIES
  larcorealg::Geometry
  larevt::CalibrationDBI_IOVData larevt::CalibrationDBI_Providers
  larreco::HitFinder lardataobj::RecoBase larcore::Geometry_Geometry_service
  lardata::Utilities lardataalg::DetectorInfo
  nurandom::RandomUtils_NuRandomService_service PkgConfig::FFTW3
  art::Framework_Core art::Framework_Principal
  art::Framework_Services_Registry art::Persistency_Common
  art::Persistency_Provenance
  icarus_signal_processing::icarus_signal_processing
  art_root_io::tfile_support ROOT::Core
  art::Framework_Services_Optional_RandomNumberGenerator_service
  art_root_io::TFileService_service
  art::Framework_Services_System_FileCatalogMetadata_service art::Utilities
  canvas::canvas messagefacility::MF_MessageLogger messagefacility::headers
  fhiclcpp::fhiclcpp cetlib::cetlib Boost::filesystem Boost::system
  CLHEP::Random Eigen3::Eigen messagefacility::MF_MessageLogger
  )
endforeach()

cet_build_plugin(CoherentNoiseFactor art::service LIBRARIES
  larcorealg::Geometry
  larevt::CalibrationDBI_IOVData larevt::CalibrationDBI_Providers
  larreco::HitFinder lardataobj::RecoBase larcore::Geometry_Geometry_service
  lardata::Utilities lardataalg::DetectorInfo
  nurandom::RandomUtils_NuRandomService_service PkgConfig::FFTW3
  art::Framework_Core art::Framework_Principal
  art::Framework_Services_Registry art::Persistency_Common
  art::Persistency_Provenance
  icarus_signal_processing::icarus_signal_processing
  art_root_io::tfile_support ROOT::Core
  art::Framework_Services_Optional_RandomNumberGenerator_service
  art_root_io::TFileService_service
  art::Framework_Services_System_FileCatalogMetadata_service art::Utilities
  canvas::canvas messagefacility::MF_MessageLogger messagefacility::headers
  fhiclcpp::fhiclcpp cetlib::cetlib Boost::filesystem Boost::system
  CLHEP::Random Eigen3::Eigen icaruscode_TPC_Simulation_DetSim_tools
  messagefacility::MF_MessageLogger )

install_headers()
install_source()
install_fhicl()
