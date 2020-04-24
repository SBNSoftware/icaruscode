#ifndef IC_CRTDETSIMALG_H
#define IC_CRTDETSIMALG_H

//art includes
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"
#include "nurandom/RandomUtils/NuRandomService.h"

//larsoft includes
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataalg/DetectorInfo/ElecClock.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/AuxDetGeo.h"
#include "larcore/Geometry/AuxDetGeometry.h"
#include "larcorealg/CoreUtils/NumericUtils.h"

//CLHEP includes
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandPoisson.h"

//C++ includes
#include <cmath>
#include <map>
#include <set>
#include <vector>
#include <string>
#include <utility>

//ROOT includes
#include "TGeoManager.h"
#include "TGeoNode.h"

//CRT includes
#include "icaruscode/CRT/CRTProducts/CRTChannelData.h"
#include "icaruscode/CRT/CRTProducts/CRTData.hh"
#include "icaruscode/CRT/CRTUtils/CRTCommonUtils.h"

using std::vector;
using std::pair;

namespace icarus {
 namespace crt {
    class CRTDetSimAlg;
 }
}

struct Tagger {
    char type;
    int modid;
    std::string reg; //crt region where FEB is located
    std::set<int> layerid; //keep track of layers hit accross whole event window
    std::map<int,int> chanlayer; //map chan # to layer
    std::vector<icarus::crt::CRTChannelData> data; //time and charge info for each channel > thresh
    std::vector<int> ichan; //index of the auxDetIDE vector provided by DetSimProducer 
};//Tagger


class icarus::crt::CRTDetSimAlg {

 public:

    CRTDetSimAlg(fhicl::ParameterSet const & p, CLHEP::HepRandomEngine& fRandEngine);
    void reconfigure(fhicl::ParameterSet const & p);

    void ClearTaggers();
    //given a vector of AuxDetIDEs, fill a map of tagger objects with intermediate ChannelData + aux info
    void FillTaggers(const uint32_t adid, const uint32_t adsid, const std::unique_ptr<vector<sim::AuxDetIDE>>& ides, const int nide);
    vector<pair<CRTData, vector<int>>> CreateData();//vector<art::Ptr<sim::AuxDetSimChannel>> channels);


 private:

    //fhicl configurable vars
    bool   fVerbose;
    bool   fUltraVerbose;
    double fGlobalT0Offset;    //!< Time delay fit: Gaussian normalization
    double fTDelayNorm;        //!< Time delay fit: Gaussian normalization
    double fTDelayShift;       //!< Time delay fit: Gaussian x shift
    double fTDelaySigma;       //!< Time delay fit: Gaussian width
    double fTDelayOffset;      //!< Time delay fit: Gaussian baseline offset
    double fTDelayRMSGausNorm; //!< Time delay RMS fit: Gaussian normalization
    double fTDelayRMSGausShift;//!< Time delay RMS fit: Gaussian x shift
    double fTDelayRMSGausSigma;//!< Time delay RMS fit: Gaussian width
    double fTDelayRMSExpNorm;  //!< Time delay RMS fit: Exponential normalization
    double fTDelayRMSExpShift; //!< Time delay RMS fit: Exponential x shift
    double fTDelayRMSExpScale; //!< Time delay RMS fit: Exponential scale
    double fQ0;                //!< Average energy deposited for mips in 1cm for charge scaling [GeV]
    double fQPed;              //!< ADC offset for the single-peak peak mean [ADC]
    double fQSlope;            //!< Slope in mean ADC / Npe [ADC]
    double fQRMS;              //!< ADC single-pe spectrum width [ADC]
    double fQThresholdC;       //!< ADC charge threshold for CERN system [ADC]
    double fQThresholdM;       //!< ADC charge threshold for MINOS system [ADC]
    double fQThresholdD;       //!< ADC charge threshold for DC system [ADC]
    double fTResInterpolator;  //!< Interpolator time resolution [ns]
    double fPropDelay;         //!< Delay in pulse arrival time [ns/m]
    double fPropDelayError;    //!< Delay in pulse arrival time, uncertainty [ns/m]
    double fStripCoincidenceWindow;  //!< Time window for two-fiber coincidence [ns]
    bool   fApplyStripCoinC;   //!< Whether or not to apply coincence between fibers in a strip (c modules only)
    bool   fApplyCoincidenceC; //!< Whether or not to apply coincidence between hits in adjacent layers
    bool   fApplyCoincidenceM; //!< Whether or not to apply coincidence between hits in adjacent layers
    bool   fApplyCoincidenceD; //!< Whether or not to apply coincidence between hits in adjacent layers
    double fLayerCoincidenceWindowC;  //!< Time window for two layer coincidence in a CERN module [ns]
    double fLayerCoincidenceWindowM;  //!< Time window for two layer coincidence between MINOS modules [ns]
    double fLayerCoincidenceWindowD;  //!< Time window for two layer coincidence in a DC module [ns]
    bool   fUseEdep;  //!< Use the true G4 energy deposited, assume mip if false.
    double fDeadTime; //!< Dead Time inherent in the front-end electronics
    double fBiasTime; //!< Hard cut off for follow-up hits after primary trigger to bias ADC level

    //bookkeeping / stats vars
    int fNsim_m, fNsim_d, fNsim_c; //number of strips in each subsystem with deposited energy
    int fNchandat_m, fNchandat_d, fNchandat_c; //number of SiPM channel signals above threshold
    int fNmissthr_c, fNmissthr_d, fNmissthr_m; //number of channel signals below threshold
    int fNmiss_strcoin_c; //number of channel signals missed due to no fiber-fiber coincidence in a cern strip
    int fNdual_m; //number of energy deposits producing signals above threshold at both ends of a minos strip

    std::map<int,int> fRegCounts;
    std::set<int> fRegions;

    CLHEP::HepRandomEngine& fRandEngine;
    std::map<int,vector<pair<int,int>>> fFebMap;

    // A list of hit taggers, before any coincidence requirement
    std::map<int, Tagger> fTaggers;

    pair<double,double> GetTransAtten(const double pos);
    double GetLongAtten(const double dist);

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
    //uint32_t
    double GetChannelTriggerTicks(detinfo::ElecClock& clock,
                                  float t0, float npeMean, float r);

};

#endif
