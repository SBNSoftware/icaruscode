///////////////////////////////////////////////////////////////////////////////
/// Class: CRTDetSim
/// Module Type: producer
/// File: CRTDetSim_module.cc
///
/// Based on LArIAT TOFSimDigits.cc (Author: Lucas Mendes Santos)
/// with modifications for SBND (Author: mastbaum@uchicago.edu)
//
/// Author: chilge@rams.colostate.edu
///////////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"

#include "lardataalg/DetectorInfo/ElecClock.h"
#include "icaruscode/CRT/CRTProducts/CRTChannelData.h"
#include "icaruscode/CRT/CRTProducts/CRTData.hh"

#include <string>
#include <utility>

namespace icarus {
namespace crt {

class CRTDetSim : public art::EDProducer {
public:
  explicit CRTDetSim(fhicl::ParameterSet const & p);

  CRTDetSim(CRTDetSim const &) = delete;
  CRTDetSim(CRTDetSim &&) = delete;
  CRTDetSim& operator = (CRTDetSim const &) = delete;
  CRTDetSim& operator = (CRTDetSim &&) = delete;
  void reconfigure(fhicl::ParameterSet const & p);

  void produce(art::Event & e) override;
  std::string fG4ModuleLabel;

private:
  char GetAuxDetType(geo::AuxDetGeo const& adgeo);
  std::string GetAuxDetRegion(geo::AuxDetGeo const& adgeo);
  uint8_t GetStackNum(geo::AuxDetGeo const& adgeo);
  //bool TimeOrderCRTData(icarus::crt::CRTData::ChannelData crtdat1, icarus::crt::CRTData::ChannelData crtdat2);
  /**
   * Get the channel trigger time relative to the start of the MC event.
   *
   * @param engine The random number generator engine
   * @param clock The clock to count ticks on
   * @param t0 The starting time (which delay is added to)
   * @param npe Number of observed photoelectrons
   * @param r Distance between the energy deposit and strip readout end [mm]
   * @return Trigger clock ticks at this true hit time
   */
  uint32_t GetChannelTriggerTicks(CLHEP::HepRandomEngine* engine,
                                detinfo::ElecClock& clock,
                                float t0, float npeMean, float r);

  bool   fVerbose;
  double fGlobalT0Offset;  //!< Time delay fit: Gaussian normalization
  double fTDelayNorm;  //!< Time delay fit: Gaussian normalization
  double fTDelayShift;  //!< Time delay fit: Gaussian x shift
  double fTDelaySigma;  //!< Time delay fit: Gaussian width
  double fTDelayOffset;  //!< Time delay fit: Gaussian baseline offset
  double fTDelayRMSGausNorm;  //!< Time delay RMS fit: Gaussian normalization
  double fTDelayRMSGausShift;  //!< Time delay RMS fit: Gaussian x shift
  double fTDelayRMSGausSigma;  //!< Time delay RMS fit: Gaussian width
  double fTDelayRMSExpNorm;  //!< Time delay RMS fit: Exponential normalization
  double fTDelayRMSExpShift;  //!< Time delay RMS fit: Exponential x shift
  double fTDelayRMSExpScale;  //!< Time delay RMS fit: Exponential scale
  double fQ0;  //!< Average energy deposited for mips in 1cm for charge scaling [GeV]
  double fQPed;  //!< ADC offset for the single-peak peak mean [ADC]
  double fQSlope;  //!< Slope in mean ADC / Npe [ADC]
  double fQRMS;  //!< ADC single-pe spectrum width [ADC]
  double fQThresholdC;  //!< ADC charge threshold for CERN system [ADC]
  double fQThresholdM;  //!< ADC charge threshold for MINOS system [ADC]
  double fQThresholdD;  //!< ADC charge threshold for DC system [ADC]
  double fTResInterpolator;  //!< Interpolator time resolution [ns]
  double fPropDelay;  //!< Delay in pulse arrival time [ns/m]
  double fPropDelayError;  //!< Delay in pulse arrival time, uncertainty [ns/m]
  double fStripCoincidenceWindow;  //!< Time window for two-fiber coincidence [ns]
  bool fApplyCoincidenceC; //!< Whether or not to apply coincidence between hits in adjacent layers
  bool fApplyCoincidenceM; //!< Whether or not to apply coincidence between hits in adjacent layers
  bool fApplyCoincidenceD; //!< Whether or not to apply coincidence between hits in adjacent layers
  double fLayerCoincidenceWindowC;  //!< Time window for two layer coincidence in a CERN module [ns]
  double fLayerCoincidenceWindowM;  //!< Time window for two layer coincidence between MINOS modules [ns]
  double fLayerCoincidenceWindowD;  //!< Time window for two layer coincidence in a DC module [ns]
  //double fAbsLenEffC;  //!< Effective abs. length for transverse Npe scaling in CERN scintillator [cm]
  //double fAbsLenEffM;  //!< Effective abs. length for transverse Npe scaling in MINOS scintillator [cm]
  //double fAbsLenEffD;  //!< Effective abs. length for transverse Npe scaling in DC scintillator [cm]
  bool fUseEdep;  //!< Use the true G4 energy deposited, assume mip if false.
  double fDeadTime; //!< Dead Time inherent in the front-end electronics
  double fBiasTime; //!< Hard cut off for follow-up hits after primary trigger to bias ADC level
};

}  // namespace crt
}  // namespace icarus

