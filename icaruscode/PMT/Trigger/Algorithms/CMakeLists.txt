art_make_library(
  SUBDIRS
    "details"
  LIBRARIES
    icaruscode::PMT_Algorithms
    icarusalg::Utilities
    icaruscode::IcarusObj
    sbnobj::ICARUS_PMT_Trigger_Data
    sbnobj::ICARUS_PMT_Data
    sbnalg::Utilities
    lardataalg::DetectorInfo
    lardataalg::MCDumpers
    lardataalg::headers                # lardataalg/Utilities
    larcorealg::Geometry
    larcorealg::CoreUtils
    lardataobj::RawData
    lardataobj::Simulation
    larcoreobj::headers                # larcoreobj/SimpleTypesAndConstants
    nusimdata::SimulationBase
    messagefacility::MF_MessageLogger
    canvas::canvas                     # art::EventID, art::InputTag
    cetlib_except::cetlib_except
    fhiclcpp::fhiclcpp
    range-v3::range-v3
    FFTW3f::FFTW3f
    ROOT::Hist                         # TGraph.h
    ROOT::Tree
    ROOT::RIO                          # TDirectory.h
    ROOT::Core                         # Rtypes.h
  )

install_headers(SUBDIRS "details")
install_source(SUBDIRS "details")
install_fhicl(SUBDIRS "details")
