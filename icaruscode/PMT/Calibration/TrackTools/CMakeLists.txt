cet_make_library(
  SOURCE
    StoppingTrackSelector.cxx
  LIBRARIES
    PUBLIC
      sbnobj::Common_Calibration_dict   # sbn::TrackInfo
      sbnobj::Common_CRT                # sbn::crt::CRTHitT0TaggingInfo
      sbncode::Utilities
      icarusalg::Utilities              # lar::util::TrackTimeInterval
      lardataobj::AnalysisBase          # anab::Calorimetry, anab::T0
      lardataobj::RecoBase
      larcore::Geometry_Geometry_service
      larcorealg::Geometry
      lardata::DetectorInfoServices_DetectorPropertiesServiceStandard_service
      lardataalg::DetectorInfo
      art::Framework_Principal
      art::Framework_Services_Registry
      art::Utilities
      canvas::canvas
      fhiclcpp::fhiclcpp
      messagefacility::MF_MessageLogger
      cetlib_except::cetlib_except
)

install_headers()
install_source()
