art_make(
  LIB_LIBRARIES
    larana::OpticalDetector_OpHitFinder
    sbnobj::ICARUS_PMT_Data
    lardataobj::RecoBase
    lardataobj::RawData
    messagefacility::MF_MessageLogger
    fhiclcpp::fhiclcpp
    cetlib_except::cetlib_except
  )

install_headers()
install_source()
