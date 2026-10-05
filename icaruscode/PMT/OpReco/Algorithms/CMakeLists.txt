art_make_library(
  LIBRARIES
    larana::OpticalDetector_OpHitFinder
    lardataalg::UtilitiesHeaders
    lardataalg::DetectorInfo
    lardataobj::RecoBase
    sbnobj::ICARUS_PMT_Data
    lardataobj::RawData
    messagefacility::MF_MessageLogger
    fhiclcpp::fhiclcpp
    cetlib_except::cetlib_except
  )

install_headers()
install_source()
