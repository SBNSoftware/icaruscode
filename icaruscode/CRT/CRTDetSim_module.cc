////////////////////////////////////////////////////////////////////////////////
/// \file CRTDetSim_module.cc
///
/// Based on LArIAT TOFSimDigits.cc (Author: Lucas Mendes Santos)
/// with modifications for SBND (Author: mastbaum@uchicago.edu)
//
/// Author: chilge@rams.colostate.edu
////////////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "nutools/RandomUtils/NuRandomService.h"

#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataalg/DetectorInfo/ElecClock.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/AuxDetGeo.h"
#include "larcore/Geometry/AuxDetGeometry.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/CoreUtils/NumericUtils.h"

#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandPoisson.h"

#include "TFile.h"
#include "TNtuple.h"
#include "TGeoManager.h"
#include "TGeoNode.h"
#include "icaruscode/CRT/CRTData.hh"
#include "icaruscode/CRT/CRTDetSim.h"

#include <cmath>
#include <memory>
#include <string>

namespace icarus {
namespace crt {

void CRTDetSim::reconfigure(fhicl::ParameterSet const & p) {
  fG4ModuleLabel = p.get<std::string>("G4ModuleLabel");

  fGlobalT0Offset = p.get<double>("GlobalT0Offset");
  fTDelayNorm = p.get<double>("TDelayNorm");
  fTDelayShift = p.get<double>("TDelayShift");
  fTDelaySigma = p.get<double>("TDelaySigma");
  fTDelayOffset = p.get<double>("TDelayOffset");
  fTDelayRMSGausNorm = p.get<double>("TDelayRMSGausNorm");
  fTDelayRMSGausShift = p.get<double>("TDelayRMSGausShift");
  fTDelayRMSGausSigma = p.get<double>("TDelayRMSGausSigma");
  fTDelayRMSExpNorm = p.get<double>("TDelayRMSExpNorm");
  fTDelayRMSExpShift = p.get<double>("TDelayRMSExpShift");
  fTDelayRMSExpScale = p.get<double>("TDelayRMSExpScale");
  fPropDelay = p.get<double>("PropDelay");
  fPropDelayError = p.get<double>("PropDelayError");
  fTResInterpolator = p.get<double>("TResInterpolator");
  fNpeScaleNorm = p.get<double>("NpeScaleNorm");
  fNpeScaleShift = p.get<double>("NpeScaleShift");
  fUseEdep = p.get<bool>("UseEdep");
  fQ0 = p.get<double>("Q0");
  fQPed = p.get<double>("QPed");
  fQSlope = p.get<double>("QSlope");
  fQRMS = p.get<double>("QRMS");
  fQThreshold = p.get<double>("QThreshold");
  fStripCoincidenceWindow = p.get<double>("StripCoincidenceWindow");
  fLayerCoincidenceWindow = p.get<double>("LayerCoincidenceWindow");
  fAbsLenEff = p.get<double>("AbsLenEff");
}

// constructor
CRTDetSim::CRTDetSim(fhicl::ParameterSet const & p) {
  art::ServiceHandle<rndm::NuRandomService> seeds;
  seeds->createEngine(*this, "HepJamesRandom", "crt", p, "Seed");

  this->reconfigure(p);

  produces<std::vector<icarus::crt::CRTData> >();
}

//function takes reference to AuxDetGeo object and gives parent subsystem
char CRTDetSim::GetAuxDetType(geo::AuxDetGeo const& adgeo)
{
  std::string base = "volAuxDet_Module_010_";
  std::string volName(adgeo.TotalVolume()->GetName());
  std::string reg  = volName.substr(base.size(),volName.size());

  if(reg == "Top")        return 'c';
  if(reg == "SlopeLeft")  return 'c';
  if(reg == "SlopeRight") return 'c';
  if(reg == "SlopeFront") return 'c';
  if(reg == "SlopeBack")  return 'c';
  if(reg == "Left")       return 'm';
  if(reg == "Right")      return 'm';
  if(reg == "Front")      return 'm';
  if(reg == "Back")       return 'm';
  if(reg == "Bottom")     return 'd';
  std::cout << " Region not found for module " << volName << std::endl;
  return 'e';
}


//function for simulating time response
//  takes true hit time, LY(PE) observed, and longitudinal distance from readout
//  uses 12 FHiCL configurable parameters
//  returns simulated time in units of clock ticks
uint32_t CRTDetSim::GetChannelTriggerTicks(CLHEP::HepRandomEngine* engine,
                                         detinfo::ElecClock& clock,
                                         float t0, float npeMean, float r) {
  // Hit timing, with smearing and NPE dependence
  double tDelayMean = \
    fTDelayNorm *
      exp(-0.5 * pow((npeMean - fTDelayShift) / fTDelaySigma, 2)) +
    fTDelayOffset;

  double tDelayRMS = \
    fTDelayRMSGausNorm *
      exp(-pow(npeMean - fTDelayRMSGausShift, 2) / fTDelayRMSGausSigma) +
    fTDelayRMSExpNorm *
      exp(-(npeMean - fTDelayRMSExpShift) / fTDelayRMSExpScale);

  double tDelay = CLHEP::RandGauss::shoot(engine, tDelayMean, tDelayRMS);

  // Time resolution of the interpolator
  tDelay += CLHEP::RandGauss::shoot(engine, 0, fTResInterpolator);

  // Propagation time
  double tProp = CLHEP::RandGauss::shoot(fPropDelay, fPropDelayError) * r;

  double t = t0 + tProp + tDelay;

  // Get clock ticks
  clock.SetTime(t / 1e3);  // SetTime takes microseconds

  mf::LogInfo("CRT")
    << "CRT TIMING: t0=" << t0
    << ", tDelayMean=" << tDelayMean << ", tDelayRMS=" << tDelayRMS
    << ", tDelay=" << tDelay << ", tDelay(interp)="
    << tDelay << ", tProp=" << tProp << ", t=" << t << ", ticks=" << clock.Ticks() << "\n";

  return clock.Ticks();
}


struct Tagger {
  std::set<unsigned> layerID;
  std::vector<icarus::crt::CRTData> data;
};

//module producer
void CRTDetSim::produce(art::Event & e) {
  // A list of hit taggers, before any coincidence requirement
  std::map<int, Tagger> taggers;

  // Services: Geometry, DetectorClocks, RandomNumberGenerator
  art::ServiceHandle<geo::Geometry> geoService;
  art::ServiceHandle<detinfo::DetectorClocksService> detClocks;
  detinfo::ElecClock trigClock = detClocks->provider()->TriggerClock();

  art::ServiceHandle<art::RandomNumberGenerator> rng;
  CLHEP::HepRandomEngine* engine = &rng->getEngine("crt");

  // Handle for (truth) AuxDetSimChannels
  art::Handle<std::vector<sim::AuxDetSimChannel> > channels;
  e.getByLabel(fG4ModuleLabel, channels);

  // Loop through truth AD channels
  for (auto& adsc : *channels) {

    const geo::AuxDetGeo& adGeo = \
        geoService->AuxDet(adsc.AuxDetID());

    //check stripID is consistent with number of sensitive volumes
    if( adGeo.NSensitiveVolume() < adsc.AuxDetSensitiveID()+1){
        std::cout << "adsID out of bounds! Skipping..." << "\n"
                  << "   " << adGeo.Name()  << " / modID "   << adsc.AuxDetID()
                  << " / stripID " << adsc.AuxDetSensitiveID() 
        << std::endl;
        continue;
    }

    const geo::AuxDetSensitiveGeo& adsGeo = \
        adGeo.SensitiveVolume(adsc.AuxDetSensitiveID());

    //parent subsystem of module (c,m,d), e if not found
    char auxDetType = GetAuxDetType(adGeo);

    // Simulate the CRT response for each hit
    for (auto ide : adsc.AuxDetIDEs()) {

      // Finally, what is the distance from the hit (centroid of the entry
      // and exit points) to the readout end?
      double x = (ide.entryX + ide.exitX) / 2;
      double y = (ide.entryY + ide.exitY) / 2;
      double z = (ide.entryZ + ide.exitZ) / 2;
      double world[3] = {x, y, z};
      double svHitPosLocal[3];
      double modHitPosLocal[3];
      adsGeo.WorldToLocal(world, svHitPosLocal); //position in strip frame (origin at center)
      adGeo.WorldToLocal(world, modHitPosLocal); //position in module fram (origin at center)

      //for CERN and DC modules, check which layer of strips hit is in for dual layer coincidence
      unsigned layid = 0;
      if (auxDetType=='c'||auxDetType=='d') layid = (modHitPosLocal[1] > 0);
      //else layid = 0;

      //longitudinal distance (cm) along the strip for fiber atten. calculation
      double distToReadout = abs( adsGeo.HalfLength() - svHitPosLocal[2]);
      double distToReadout2 = abs(-adsGeo.HalfHeight() - svHitPosLocal[2]);

      // The expected number of PE, using a quadratic model for the distance
      // dependence, and scaling linearly with deposited energy.
      double qr = fUseEdep ? 1.0 * ide.energyDeposited / fQ0 : 1.0;

      double npeExpected = \
        fNpeScaleNorm / pow(distToReadout - fNpeScaleShift, 2) * qr;
      double npeExpected2 = \
        fNpeScaleNorm / pow(distToReadout2 - fNpeScaleShift, 2) * qr;

      // Put PE on channels weighted by transverse distance across the strip,
      // using an exponential model
      double d0=0., d1=0.;

      switch(auxDetType){
          case 'c' : 
              d0 = abs(-adsGeo.HalfWidth1() + 6 - svHitPosLocal[0]);  // L (fibers 6cm from edge)
              d1 = abs(adsGeo.HalfWidth1() - 6  - svHitPosLocal[0]);  // R
              break;
          case 'm' : 
              d0 = abs(svHitPosLocal[1]); //one fiber centered
              d1 = 0.;
              break;
          case 'd' : 
              d0 = abs(svHitPosLocal[0]); //one fiber centered
              d1 = 0.;
              break;
      }

      //FIX ME (ONLY DOES SOMETHING FOR C MODULES!
      double abs0 = exp(-d0 / fAbsLenEff);
      double abs1 = exp(-d1 / fAbsLenEff);
      double npeExp0 = npeExpected * abs0 / (abs0 + abs1);
      double npeExp1 = npeExpected * abs1 / (abs0 + abs1);
      double npeExp0Dual = npeExpected2 * abs0 / (abs0 + abs1);

      // Observed PE (Poisson-fluctuated)
      long npe0 = CLHEP::RandPoisson::shoot(engine, npeExp0);
      long npe1 = CLHEP::RandPoisson::shoot(engine, npeExp1);
      long npe0Dual = CLHEP::RandPoisson::shoot(engine, npeExp0Dual);

      // Time relative to trigger, accounting for propagation delay and 'walk'
      // for the fixed-threshold discriminator
      double tTrue = (ide.entryT + ide.exitT) / 2 + fGlobalT0Offset;
      uint32_t t0 = \
        GetChannelTriggerTicks(engine, trigClock, tTrue, npe0, distToReadout);
      uint32_t t1 = \
        GetChannelTriggerTicks(engine, trigClock, tTrue, npe1, distToReadout);
      uint32_t t0Dual = \
        GetChannelTriggerTicks(engine, trigClock, tTrue, npe0Dual, distToReadout2);

      // Time relative to PPS: Random for now! (FIXME)
      uint32_t ppsTicks = \
        CLHEP::RandFlat::shootInt(engine, trigClock.Frequency() * 1e6);

      // SiPM and ADC response: Npe to ADC counts
      short q0 = \
        CLHEP::RandGauss::shoot(engine, fQPed + fQSlope * npe0, fQRMS * sqrt(npe0));
      short q1 = \
        CLHEP::RandGauss::shoot(engine, fQPed + fQSlope * npe1, fQRMS * sqrt(npe1));
      short q0Dual = \
        CLHEP::RandGauss::shoot(engine, fQPed + fQSlope * npe0Dual, fQRMS * sqrt(npe0Dual));

      // Adjacent channels on a strip are numbered sequentially.
      //
      // In the AuxDetChannelMapAlg methods, channels are identified by an
      // AuxDet name (retrievable given the hit AuxDet ID) which specifies a
      // module, and a channel number from 0 to 32.

      std::string modBase = "volAuxDet_Module_";
      std::string stripBase = "volAuxDetSensitive_Module_xxx_Strip_";

      std::string adName = adGeo.TotalVolume()->GetName();
      std::string adsName = adsGeo.TotalVolume()->GetName();

      uint32_t moduleID = atoi(adName.substr( modBase.size(), 3).c_str());
      uint32_t stripID = atoi(adsName.substr( stripBase.size(), 2).c_str());

      //uint32_t moduleID = adsc.AuxDetID();
      //uint32_t stripID = adsc.AuxDetSensitiveID();
      uint32_t channel0ID=0, channel1ID=0;

      switch (auxDetType){
          case 'c' : //channels 4052-7956
              channel0ID = 32 * moduleID + 2 * stripID + 0 - 1132;
              channel1ID = 32 * moduleID + 2 * stripID + 1 - 1132;
              break;
          case 'd' : //channels 3156-4051
              channel0ID = 64 * moduleID + stripID - 6316;
              break;
          case 'm' : //channels 0-1577, 1588-3155
              channel0ID = 32 * (moduleID/3) + stripID/2 + 10*(moduleID % 3);
              channel1ID = 32 * (moduleID/3) + stripID/2 + 10*(moduleID % 3) + 1578;
              break;

      }

      // Apply ADC threshold and strip-level coincidence (both fibers fire)
      if (auxDetType=='c' && q0 > fQThreshold && q1 > fQThreshold && util::absDiff(t0, t1) < fStripCoincidenceWindow) {
              Tagger& tagger = taggers[moduleID];
              tagger.layerID.insert(layid);
              tagger.data.push_back(icarus::crt::CRTData(channel0ID, t0, ppsTicks, q0));
              tagger.data.push_back(icarus::crt::CRTData(channel1ID, t1, ppsTicks, q1));
      }//if fiber-fiber coincidence

      if (auxDetType=='d' && q0 > fQThreshold) {
              Tagger& tagger = taggers[moduleID];
              tagger.layerID.insert(layid);
              tagger.data.push_back(icarus::crt::CRTData(channel0ID, t0, ppsTicks, q0));
      }//if fiber-fiber coincidence

      if (auxDetType=='m' && (q0 > fQThreshold || q0Dual > fQThreshold)) {
              Tagger& tagger = taggers[moduleID];
              tagger.layerID.insert(layid);
              if(q0 > fQThreshold) tagger.data.push_back(icarus::crt::CRTData(channel0ID, t0, ppsTicks, q0));
              if(q0Dual > fQThreshold) tagger.data.push_back(icarus::crt::CRTData(channel1ID, t0Dual, ppsTicks, q0Dual));
      }//if

      mf::LogInfo("CRT")
        << "CRT HIT VOL " << (adGeo.TotalVolume())->GetName() << " with " << adGeo.NSensitiveVolume() << " AuxDetSensitive volumes" << "\n"
        << "CRT HIT SENSITIVE VOL " << (adsGeo.TotalVolume())->GetName() << "\n"
        << "CRT HIT AuxDetID " <<  adsc.AuxDetID() << " / AuxDetSensitiveID " << adsc.AuxDetSensitiveID() << "\n"
        << "CRT HIT POS " << x << " " << y << " " << z << "\n"
        << "CRT STRIP POS " << svHitPosLocal[0] << " " << svHitPosLocal[1] << " " << svHitPosLocal[2] << "\n"
        << "CRT MODULE POS " << modHitPosLocal[0] << " " << modHitPosLocal[1] << " "<< modHitPosLocal[2] << " " << "\n"
        << "CRT layer ID: " << layid << "\n"
        << "CRT distToReadout: " << distToReadout << "\n"
        << "CRT q0: " << q0 << ", q1: " << q1 << ", t0: " << t0 << ", t1: " << t1 << ", dt: " << util::absDiff(t0,t1) << "\n";
    }//for AuxDetIDEs 
  }//for AuxDetChannels

  // Apply coincidence trigger requirement
  std::unique_ptr<std::vector<icarus::crt::CRTData> > triggeredCRTHits(
      new std::vector<icarus::crt::CRTData>);

  // Logic: For normal taggers, require at least one hit in each perpendicular
  // plane. For the bottom tagger, any hit triggers read out.
  for (auto trg : taggers) {
    if ((trg.first >= 163 && trg.first <= 283 && trg.second.layerID.size() > 1 )
           || trg.first < 163 || trg.first > 283){
      for (auto d : trg.second.data) {
        triggeredCRTHits->push_back(d);
      }
    }
  }

  mf::LogInfo("CRT") << "CRT TRIGGERED HITS: " << triggeredCRTHits->size() << "\n";

  e.put(std::move(triggeredCRTHits));
}

DEFINE_ART_MODULE(CRTDetSim)

}  // namespace crt
}  // namespace icarus

