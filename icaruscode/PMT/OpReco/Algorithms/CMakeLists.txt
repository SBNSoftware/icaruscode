art_make(
  LIB_LIBRARIES
    larana_OpticalDetector_OpHitFinder
    sbnobj_ICARUS_PMT_Data
    lardataobj_RecoBase
    lardataobj_RawData
    messagefacility::MF_MessageLogger
    fhiclcpp::fhiclcpp
    cetlib_except::cetlib_except
  )

install_headers()
install_source()
