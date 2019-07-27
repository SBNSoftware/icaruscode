///////////////////////////////////////////////////////////////////////////////
/// Class: CRTDetSim
/// Module Type: producer
/// \file CRTDetSim_module.cc
///
/// Based on LArIAT TOFSimDigits.cc (Author: Lucas Mendes Santos)
/// with modifications for SBND (Author: mastbaum@uchicago.edu)
/// then modified for ICARUS
///
/// Author: Chris.Hilgenberg@colostate.edu
///
/// Extracts position, time, and energy deposited, aint with
///   trackID associated with deposit from AuxDetSimChannel objects
///   produced by G4 stage. The trackID is used for truth matching.
/// Simulated CRT data products are intended to match the known
///   front-end electronics behavior and basic DAQ format anticipated.
///   Each hit strip is assigned a channel and associated with 
//      1) T0 - hit time relative to t0 (beam signal)
//      2) T1 - hit time relative to PPS from GPS
//      3) ADC from hit
///   Effects included are
///     * time 
///         - light propegation delay in scintillator and WLS fiber
///         - smearing in the arrival time at the SiPM due to scattering
///         - amplitude dependant smearing due to the discriminator
///     * position
///         - AuxDetSimChannel provides true entry and exit point in scintillator strip
///         - take arithmetic mean as true "hit" position
///         - "hit" position defines transverse distance (hit to fiber) and
///           intitudinal distance (length of WLS fiber between fiber entry and SiPM)
///     * energy deposited
///         - take true energy deposited in scintillator strip
///         - convert true energy to photons yielded (taken from measurements on
///               normally incident MIP muons
///         - using known attenuation lengths for bulk scinillator and WLS fiber,
///               apply attenuation correction using transverse and intitudinal
///               propegation distances respectively and a simple exponential model
///         - for the attenuated light arriving at the SiPM, correct for counting
///               statistics sampling from Poisson
///         - convert N photons seen by SiPM into ADCs using known pedestal and 
///               gain values for each channel (assumed all the same for now) 
///               ADC_i = gain_i * nphotons + ped_i
///         - smear ADC using guasian with ADC_i as mean and 
///               width = width_pedestal+sqrt(nphotons)
///     * front-end electonics deadtime
///         - sort vector of data products by time (T0)
///         - for C and D type modules, intermodule coincidence is applied
///            - one strip from each of 2 layers must be above threshold
///         - for M modules, only one channel is required to be above threshold
///         - the earliest channel that was part of the trigger provides the T0
///            defining the readout window (set in FHiCL)
///         - for each channel in the readout window with an entry above threshold,
///            the data (charge and time) is added the "track and hold list"
///         - channels in track and hold do not accept new data during the readout window
///         - after the readout window has passed, no channels can receive data until the 
///           (FHiCL congifurable) deadtime has passed. Channels are then reset and
///           the front-end board is again able to trigger
///
/////////////////////////////////////////////////////////////////////////////////////////

//art includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "nurandom/RandomUtils/NuRandomService.h"

//larsoft includes
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataalg/DetectorInfo/ElecClock.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/AuxDetGeo.h"
#include "larcore/Geometry/AuxDetGeometry.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/CoreUtils/NumericUtils.h"

//CLHEP includes
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandPoisson.h"

//ROOT includes
#include "TFile.h"
#include "TNtuple.h"
#include "TGeoManager.h"
#include "TGeoNode.h"

//C++ includes
#include <cmath>
#include <memory>
#include <string>
#include <utility>
#include <map> 
#include <vector>
#include <fstream>
#include <iostream>

//CRT includes
#include "icaruscode/CRT/CRTProducts/CRTChannelData.h"
#include "icaruscode/CRT/CRTProducts/CRTData.hh"

using std::string;

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
    string fG4ModuleLabel;
    
private:
    char GetAuxDetType(geo::AuxDetGeo const& adgeo);
    string GetAuxDetRegion(geo::AuxDetGeo const& adgeo);
    uint8_t GetStackNum(geo::AuxDetGeo const& adgeo);
    void FillFebMap(std::map< int,std::vector< std::pair<int,int> > >& m);
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
    //uint32_t
    double GetChannelTriggerTicks(detinfo::ElecClock& clock,
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
    
    CLHEP::HepRandomEngine& fRandEngine;
};

//getting parameter values from FHiCL
void CRTDetSim::reconfigure(fhicl::ParameterSet const & p) {
  fVerbose = p.get<bool>("Verbose");
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
  fUseEdep = p.get<bool>("UseEdep");
  fQ0 = p.get<double>("Q0");
  fQPed = p.get<double>("QPed");
  fQSlope = p.get<double>("QSlope");
  fQRMS = p.get<double>("QRMS");
  fQThresholdC = p.get<double>("QThresholdC");
  fQThresholdM = p.get<double>("QThresholdM");
  fQThresholdD = p.get<double>("QThresholdD");
  fStripCoincidenceWindow = p.get<double>("StripCoincidenceWindow");
  fApplyCoincidenceC = p.get<bool>("ApplyCoincidenceC");
  fApplyCoincidenceM = p.get<bool>("ApplyCoincidenceM");
  fApplyCoincidenceD = p.get<bool>("ApplyCoincidenceD");
  fLayerCoincidenceWindowC = p.get<double>("LayerCoincidenceWindowC");
  fLayerCoincidenceWindowM = p.get<double>("LayerCoincidenceWindowM");
  fLayerCoincidenceWindowD = p.get<double>("LayerCoincidenceWindowD");
  //fAbsLenEffC = p.get<double>("AbsLenEffC");
  //fAbsLenEffM = p.get<double>("AbsLenEffM");
  //fAbsLenEffD = p.get<double>("AbsLenEffD");
  fDeadTime = p.get<double>("DeadTime");
  fBiasTime = p.get<double>("BiasTime");
}

// constructor
CRTDetSim::CRTDetSim(fhicl::ParameterSet const & p) : EDProducer(p),
        fRandEngine(art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this, "HepJamesRandom", "crt", p, "Seed"))
{
  this->reconfigure(p);

  produces<std::vector<icarus::crt::CRTData> >();
}

//read from file and fill map modID->FEB(s), FEB channel subset(s). can be 2 if MINOS module (not cut)
void CRTDetSim::FillFebMap(std::map< int,std::vector< std::pair<int,int> > >& m){

    std::string geodir="/icarus/app/users/chilgenb/dev_areas/v08_22_00_prof/srcs/icaruscode/icaruscode/Geometry/gdml";
    std::ifstream fin;
    fin.open(geodir+"feb_map.txt",std::ios::in);
    if(fin.good()) std::cout << "opened file 'feb_map.txt' for reading..." << std::endl;
    else std::cout << "could not open file 'feb_map.txt' for reading!" << std::endl;
    std::vector<std::string> row;
    std::string line, word;
    while(getline(fin,line)) {
        row.clear();
        std::stringstream s(line);
        int mod;
        while (std::getline(s, word, ',')) {
            row.push_back(word);
        }
        mod = std::stoi(row[0]);
        m[mod].push_back(std::make_pair(std::stoi(row[1]),std::stoi(row[2])));
        if(row.size()>3)
            m[mod].push_back(std::make_pair(std::stoi(row[3]),std::stoi(row[4])));
    }
    std::cout << "filled febMap with " << m.size() << " entries" << std::endl;
    fin.close();
}

//function takes reference to AuxDetGeo object and gives parent subsystem
char CRTDetSim::GetAuxDetType(geo::AuxDetGeo const& adgeo)
{
  std::string volName(adgeo.TotalVolume()->GetName());
  if (volName.find("MINOS") != std::string::npos) return 'm';
  if (volName.find("CERN")  != std::string::npos) return 'c';
  if (volName.find("DC")    != std::string::npos) return 'd';

  mf::LogError("CRT") << "AuxDetType not found!" << '\n';
  return 'e';
}

//function takes reference to AuxDetGeo object and gives crt region
std::string CRTDetSim::GetAuxDetRegion(geo::AuxDetGeo const& adgeo)
{
  char type = CRTDetSim::GetAuxDetType(adgeo);
  std::string base = "volAuxDet_", region="";
  switch ( type ) {
    case 'c' : base+= "CERN"; break;
    case 'd' : base+= "DC"; break;
    case 'm' : base+= "MINOS"; break;
  }
  base+="_module_###_";
  std::string volName(adgeo.TotalVolume()->GetName());

  //module name has 2 possible formats
  //  volAuxDet_<subsystem>_module_###_<region>
  //  volAuxDet_<subsystem>_module_###_cut###_<region>

  region = volName.substr(base.length(),volName.length());
  if( region.find("_")==string::npos)
    return region;

  else
	return region.substr(region.find("_")+1,region.length());
}

int GetAuxDetRegionNum(std::string reg)
{
    if(reg == "Top")        return 30;
    if(reg == "RimWest")    return 31;
    if(reg == "RimEast")    return 32;
    if(reg == "RimSouth")   return 33;
    if(reg == "RimNorth")   return 34;
    if(reg == "WestSouth")  return 40;
    if(reg == "WestCenter") return 41;
    if(reg == "WestNorth")  return 42;
    if(reg == "EastSouth")  return 43;
    if(reg == "EastCenter") return 44;
    if(reg == "EastNorth")  return 45;
    if(reg == "South")      return 46;
    if(reg == "North")      return 47;
    if(reg == "Bottom")     return 50;
    mf::LogError("CRT") << "region not found!" << '\n';
    return INT_MAX;
}

//function for simulating time response
//  takes true hit time, LY(PE) observed, and intitudinal distance from readout
//  uses 12 FHiCL configurable parameters
//  returns simulated time in units of clock ticks
double CRTDetSim::GetChannelTriggerTicks(detinfo::ElecClock& clock,
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

  double tDelay = CLHEP::RandGauss::shoot(&fRandEngine, tDelayMean, tDelayRMS);

  // Time resolution of the interpolator
  tDelay += CLHEP::RandGauss::shoot(&fRandEngine, 0, fTResInterpolator);

  // Propagation time
  double tProp = CLHEP::RandGauss::shoot(fPropDelay, fPropDelayError) * r;

  double t = t0 + tProp + tDelay;

  // Get clock ticks
  clock.SetTime(t / 1e3);  // SetTime takes microseconds

  //if (fVerbose) mf::LogInfo("CRT")
   // << "CRT TIMING: t0=" << t0
   // << ", tDelayMean=" << tDelayMean << ", tDelayRMS=" << tDelayRMS
   // << ", tDelay=" << tDelay << ", tDelay(interp)="
   // << tDelay << ", tProp=" << tProp << ", t=" << t << "\n";//, ticks=" << clock.Ticks() << "\n"; 

  return t;//clock.Ticks();
}

bool TimeOrderCRTData(icarus::crt::CRTChannelData crtdat1, icarus::crt::CRTChannelData crtdat2) {
    return ( crtdat1.T0() < crtdat2.T0() );
}

struct Tagger {
  char type;
  std::string reg; //crt region where FEB is located
  std::set<int> layerid; //keep track of layers hit accross whole event window
  std::map<int,int> chanlayer; //map chan # to layer
  //std::pair<int,int> macPair; //which two FEBs provided coincidence (applies to m mods only)
  std::vector<icarus::crt::CRTChannelData> data; //time and charge info for each channel > thresh
};


//module producer
void CRTDetSim::produce(art::Event & e) {
  // A list of hit taggers, before any coincidence requirement
  std::map<int, Tagger> taggers;
  std::map< int,std::vector< std::pair<int,int> > > febMap;
  FillFebMap(febMap);

  // Services: Geometry, DetectorClocks, RandomNumberGenerator
  art::ServiceHandle<geo::Geometry> geoService;
  art::ServiceHandle<detinfo::DetectorClocksService> detClocks;
  detinfo::ElecClock trigClock = detClocks->provider()->TriggerClock();

  // Handle for (truth) AuxDetSimChannels
  art::Handle<std::vector<sim::AuxDetSimChannel> > channels;
  e.getByLabel(fG4ModuleLabel, channels);

  int eve = e.id().event();

  int nsim_m=0, nsim_d=0, nsim_c=0;
  int nchandat_m=0, nchandat_d=0, nchandat_c=0;
  int nmissthr_c = 0, nmissthr_d = 0, nmissthr_m = 0;
  int nmiss_strcoin_c = 0;

  std::map<int,int> regCounts;
  std::set<int> regions;

  // Loop through truth AuxDetChannel objects
  // Each channel contains the module and strip IDs and
  //  collection of position/time/deposited E for all MCParticles
  //  (AuxDetIDE) for all tracks in the event
  for (auto& adsc : *channels) {

    const int adid = adsc.AuxDetID(); //CRT module ID number (from gdml)
    const int adsid = adsc.AuxDetSensitiveID(); //CRT strip ID number (from gdml)
    //if (adsid == 0 ) continue; //skip AuxDetSensitiveID=0 (bug in AuxDetSimChannels)

    const geo::AuxDetGeo& adGeo = geoService->AuxDet(adid); //pointer to module object

    //check stripID is consistent with number of sensitive volumes
    if( (int)adGeo.NSensitiveVolume() < adsid){
        mf::LogError("CRT") << "adsID out of bounds! Skipping..." << "\n"
                  << "   " << adGeo.Name()  << " / modID "   << adid
                  << " / stripID " << adsid << '\n';
        continue;
    }

    const geo::AuxDetSensitiveGeo& adsGeo = adGeo.SensitiveVolume(adsid); //pointer to strip object
    const char auxDetType = GetAuxDetType(adGeo); //CRT module type (c, d, or m)
    if (auxDetType=='e') mf::LogError("CRT") << "COULD NOT GET AD TYPE!" << '\n';
    const std::string region = GetAuxDetRegion(adGeo); //CRT region

    int layid = INT_MAX; //set to 0 or 1 if layerid determined
    int mac5=INT_MAX, mac5dual=INT_MAX;

    // Find the path to the strip geo node, to locate it in the hierarchy
    std::set<std::string> volNames = { adsGeo.TotalVolume()->GetName() };
    std::vector<std::vector<TGeoNode const*> > paths = geoService->FindAllVolumePaths(volNames);

    std::string path = "";
    for (size_t inode=0; inode<paths.at(0).size(); inode++) {
      path += paths.at(0).at(inode)->GetName();
      if (inode < paths.at(0).size() - 1) {
        path += "/";
      }
    }

    TGeoManager* manager = geoService->ROOTGeoManager();
    manager->cd(path.c_str());
    TGeoNode* nodeStrip = manager->GetCurrentNode();
    TGeoNode* nodeInner = manager->GetMother(1);
    TGeoNode* nodeModule = manager->GetMother(2);
    double origin[3] = {0, 0, 0};

    // Module position in parent (tagger) frame
    double modulePosMother[3]; //position in CRT region volume
    nodeModule->LocalToMaster(origin, modulePosMother);

    // strip position in module frame
    double stripPosMother[3];
    double stripPosModule[3];
    nodeStrip->LocalToMaster(origin, stripPosMother);
    nodeInner->LocalToMaster(stripPosMother,stripPosModule);

    // Determine layid
    //  for C modules, two diff. lay thick
    //    1cm for top (y>0) and 1.5cm for bottom (y<0)
    if (auxDetType == 'c' || auxDetType == 'd') 
        layid = (stripPosModule[1] > 0);
       
    if (auxDetType == 'm') {
      //lateral stacks (6 total, 3 per side)
      if ( region.find("West")!=std::string::npos || region.find("East")!=std::string::npos ) {
          layid = ( abs(modulePosMother[0]>0) );
      }
      //longitudinal walls
      if ( region=="South" || region=="North" ) {
          layid = ( modulePosMother[2]> 0 );
      }
    }

    if(layid==INT_MAX) mf::LogError("CRT") << "layid NOT SET!!!" << '\n'
                               << "   ADType: " << auxDetType << '\n'
                               << "   ADRegion: " << region << '\n';

    // Simulate the CRT response for each hit in this strip
    for (auto ide : adsc.AuxDetIDEs()) {

      //if(ide.trackID!=1) continue;

      if (auxDetType=='c') nsim_c++;
      if (auxDetType=='d') nsim_d++;
      if (auxDetType=='m') nsim_m++;

      //track ID of MC paritcle depositing energy for truth matching
      std::vector<int> trkid;
      trkid.push_back(ide.trackID);

      // What is the distance from the hit (centroid of the entry
      // and exit points) to the readout end?
      double x = (ide.entryX + ide.exitX) / 2;
      double y = (ide.entryY + ide.exitY) / 2;
      double z = (ide.entryZ + ide.exitZ) / 2;
      double world[3] = {x, y, z};
      double svHitPosLocal[3];
      double modHitPosLocal[3];
      adsGeo.WorldToLocal(world, svHitPosLocal); //position in strip frame  (origin at center)
      adGeo.WorldToLocal(world, modHitPosLocal); //position in module frame (origin at center)

      if ( abs(svHitPosLocal[0])>adsGeo.HalfWidth1()+0.001 || 
           abs(svHitPosLocal[1])>adsGeo.HalfHeight()+0.001 ||
           abs(svHitPosLocal[2])>adsGeo.HalfLength()+0.001) 
         mf::LogWarning("CRT") << "HIT POINT OUTSIDE OF SENSITIVE VOLUME!" << '\n'
                            << "  AD: " << adid << " , ADS: " << adsid << '\n'
                            << "  Local position (x,y,z): ( " << svHitPosLocal[0]
                            << " , " << svHitPosLocal[1] << " , " << svHitPosLocal[2] << " )" << '\n';

      // The expected number of PE, using a quadratic model for the distance
      // dependence, and scaling linearly with deposited energy.
      double qr = fUseEdep ? ide.energyDeposited / fQ0 : 1.0;
      if (auxDetType == 'c'&& layid==0) qr *= 1.5; //c bottom layer strips 50% thicker

      //longitudinal distance (m) along the strip for fiber atten. calculation
      //assuming SiPM is on +z end (also -z for m modules)
      double distToReadout = abs( adsGeo.HalfLength() - svHitPosLocal[2])*0.01; 
      double distToReadout2 = abs(-adsGeo.HalfLength() - svHitPosLocal[2])*0.01; 

      //coefficients for quadratic fit to MINOS test data w/S14
      //obtained for normally incident cosmic muons
      //p2_m * x^2 + p1_m * x + p0_m, x is distance from readout [m]
      double p0_m = 36.5425; //initial light yield (pe) before any attenuation in m scintillator
      double p1_m = -6.3895;
      double p2_m =  0.3742;

      //coefficiencts for transverse attenuation in CERN modules 
      //from early beta-source tranverse scan data
      double at0_c = 0.682976;
      double at1_c = -0.0204477;
      double at2_c = -0.000707564;
      double at3_c = 0.000636617;
      double at4_c = 0.000147957;
      double at5_c = -3.89078e-05;

      double at0_r = 0.139941;
      double at1_r = 0.168238;
      double at2_r = -0.0198199;
      double at3_r = 0.000781752;

      double at0_l = 8.78875;
      double at1_l = 3.54602;
      double at2_l = 0.595592;
      double at3_l = 0.0449169;
      double at4_l = 0.00127892;

      //scale to LY from normally incident MIP muon  (PE)
      double npeExpected = \
        (p2_m * pow(distToReadout,2) + p1_m * distToReadout + p0_m) * qr;
      double npeExpected2 = \
        (p2_m * pow(distToReadout2,2) + p1_m * distToReadout2 + p0_m) * qr;

      // Put PE on channels weighted by transverse distance across the strip,
      // using an exponential model
      double abs0=0.0, abs1=0.0, arg=0.0; 

      switch(auxDetType){
          case 'c' :
              //hit between both fibers 
              if ( abs(svHitPosLocal[0]) <= 5.5 ) {
                arg=svHitPosLocal[0];
                abs0 = at5_c*pow(arg,5) + at4_c*pow(arg,4) + at3_c*pow(arg,3) \
		  + at2_c*pow(arg,2) + at1_c*arg + at0_c;
                abs1 = -1*at5_c*pow(arg,5) + at4_c*pow(arg,4) - at3_c*pow(arg,3) \
                  + at2_c*pow(arg,2) - at1_c*arg + at0_c;
                break;
              }
              //hit to right of both fibers
	      if ( svHitPosLocal[0] > 5.5 ) {
                arg=svHitPosLocal[0];
		abs0 = at3_r*pow(arg,3) + at2_r*pow(arg,2) + at1_r*arg + at0_r;
		abs1 = at4_l*pow(arg,4) - at3_l*pow(arg,3) \
                  + at2_l*pow(arg,2) - at1_l*arg + at0_l;
                break;
	      }
              //hit to left of both fibers
              if ( svHitPosLocal[0] < -5.5 ) {
		arg=svHitPosLocal[0];
		abs0 = at4_l*pow(arg,4) + at3_l*pow(arg,3) \
                  + at2_l*pow(arg,2) + at1_l*arg + at0_l;
                abs1 = -1*at3_r*pow(arg,3) + at2_r*pow(arg,2) - at1_r*arg + at0_r;
	      }
              break;
          case 'm' : 
              abs0 = 1.0; abs1 = 1.0;
              break;
          case 'd' : 
              abs0 = 1.0; abs1 = 1.0;
              break;
      }

      double npeExp0 = npeExpected * abs0;// / (abs0 + abs1);
      double npeExp1 = npeExpected * abs1;// / (abs0 + abs1);
      double npeExp0Dual = npeExpected2 * abs0;// / (abs0 + abs1);

      if (npeExp0<0||npeExp1<0||npeExp0Dual<0) mf::LogError("CRT") << "NEGATIVE PE!!!!!" << '\n';

      // Observed PE (Poisson-fluctuated)
      int npe0 = CLHEP::RandPoisson::shoot(&fRandEngine, npeExp0);
      int npe1 = CLHEP::RandPoisson::shoot(&fRandEngine, npeExp1);
      int npe0Dual = CLHEP::RandPoisson::shoot(&fRandEngine, npeExp0Dual);

      // Time relative to trigger [ns], accounting for propagation delay and 'walk'
      // for the fixed-threshold discriminator
      double tTrue = (ide.entryT + ide.exitT) / 2 + fGlobalT0Offset;
      double t0 = \
        GetChannelTriggerTicks(trigClock, tTrue, npe0, distToReadout*100);
      double t1 = \
        GetChannelTriggerTicks(trigClock, tTrue, npe1, distToReadout*100);
      double t0Dual = \
        GetChannelTriggerTicks(trigClock, tTrue, npe0Dual, distToReadout2*100);

      // Time relative to PPS: Random for now! (FIXME)
      int ppsTicks = \
        CLHEP::RandFlat::shootInt(&fRandEngine, trigClock.Frequency() * 1e6);

      // SiPM and ADC response: Npe to ADC counts
      int q0 = \
        CLHEP::RandGauss::shoot(&fRandEngine, fQPed + fQSlope * npe0, fQRMS * sqrt(npe0));
      int q1 = \
        CLHEP::RandGauss::shoot(&fRandEngine, fQPed + fQSlope * npe1, fQRMS * sqrt(npe1));
      int q0Dual = \
        CLHEP::RandGauss::shoot(&fRandEngine, fQPed + fQSlope * npe0Dual, fQRMS * sqrt(npe0Dual));

      if (q0<0||q1<0||q0Dual<0) mf::LogError("CRT") << "NEGATIVE ADC!!!!!" << '\n';

      // Adjacent channels on a strip are numbered sequentially.
      //
      // In the AuxDetChannelMapAlg methods, channels are identified by an
      // AuxDet name (retrievable given the hit AuxDet ID) which specifies a
      // module, and a channel number from 0 to 32.

      int channel0ID=0, channel1ID=0;

      switch (auxDetType){
          case 'c' :
              mac5 = febMap[adid][0].first;
              channel0ID = 2 * adsid + 0;
              channel1ID = 2 * adsid + 1;
              break;
          case 'd' : 
              mac5 = febMap[adid][0].first;
              channel0ID = adsid;
              break;
          case 'm' :
              mac5 = febMap[adid][0].first;
              channel0ID = adsid/2 + 10*(febMap[adid][0].second-1);
              if (febMap[adid].size()==2) 
                  mac5dual = febMap[adid][1].first;
              break;

      }

      if (mac5==INT_MAX) mf::LogError("CRT") << "mac addrs not set!" << '\n';

      // Apply ADC threshold and strip-level coincidence (both fibers fire)
      if (auxDetType=='c' && q0 > fQThresholdC && q1 > fQThresholdC && util::absDiff(t0, t1) < fStripCoincidenceWindow) {
              Tagger& tagger = taggers[mac5];
              tagger.layerid.insert(layid);
              tagger.chanlayer[channel0ID] = layid;
              tagger.chanlayer[channel1ID] = layid;
              tagger.reg = region;
              tagger.type = 'c';
              tagger.data.push_back(icarus::crt::CRTChannelData(channel0ID,t0,ppsTicks,q0,trkid));
              tagger.data.push_back(icarus::crt::CRTChannelData(channel1ID,t1,ppsTicks,q1,trkid));
              nchandat_c+=2;
      }//if fiber-fiber coincidence

      if (auxDetType=='d' && q0 > fQThresholdD) {
              Tagger& tagger = taggers[mac5];
              tagger.layerid.insert(layid);
              tagger.chanlayer[channel0ID] = layid;
              tagger.reg = region;
              tagger.type = 'd';
              tagger.data.push_back(icarus::crt::CRTChannelData(channel0ID,t0,ppsTicks,q0,trkid));
              nchandat_d++;
      }//if one strip above threshold

      if (auxDetType=='m') {
              if(q0 > fQThresholdM) {
                Tagger& tagger = taggers[mac5];
                tagger.layerid.insert(layid);
                tagger.chanlayer[channel0ID] = layid;
                tagger.reg = region;
                tagger.type = 'm';
                tagger.data.push_back(icarus::crt::CRTChannelData(channel0ID,t0,ppsTicks,q0,trkid));
                nchandat_m++;
              }
              if(q0Dual > fQThresholdM) { 
                Tagger& tagger = taggers[mac5dual];
                tagger.layerid.insert(layid);
                tagger.chanlayer[channel0ID] = layid;
                tagger.reg = region;
                tagger.type = 'm';
                tagger.data.push_back(icarus::crt::CRTChannelData(channel0ID,t0Dual,ppsTicks,q0Dual,trkid));
                nchandat_m++;
              }
      }//if one strip above threshold at either end

      //counting losses
      if (auxDetType == 'c') {
          if (q0 < fQThresholdC || q1 < fQThresholdC) nmissthr_c++;
          if ( util::absDiff(t0,t1) >= fStripCoincidenceWindow ) nmiss_strcoin_c++;
      }
      if (auxDetType == 'd' && q0 < fQThresholdD) nmissthr_d++;
      if (auxDetType == 'm' && ( q0 < fQThresholdM || q0Dual < fQThresholdM)) nmissthr_m++;

      //print detsim info (if enabled)
      if (fVerbose&&
         ( (auxDetType=='c' && q0>fQThresholdC && q1>fQThresholdC) ||
           (auxDetType=='d' && q0>fQThresholdD ) ||
           (auxDetType=='m' && (q0>fQThresholdM || q0Dual>fQThresholdM)) ))
        std::cout << '\n'
        << "CRT HIT VOL " << (adGeo.TotalVolume())->GetName() << "\n"
        << "CRT HIT SENSITIVE VOL " << (adsGeo.TotalVolume())->GetName() << "\n"
        << "CRT HIT AuxDetID " <<  adsc.AuxDetID() << " / AuxDetSensitiveID " << adsc.AuxDetSensitiveID() << "\n"
        << "CRT module type: " << auxDetType << " , CRT region: " << region << '\n'
        << "CRT channel: " << channel0ID << " , mac5: " << mac5 << '\n'
        << "CRT HIT POS " << x << " " << y << " " << z << "\n"
        << "CRT STRIP POS " << svHitPosLocal[0] << " " << svHitPosLocal[1] << " " << svHitPosLocal[2] << "\n"
        << "CRT MODULE POS " << modHitPosLocal[0] << " " << modHitPosLocal[1] << " "<< modHitPosLocal[2] << " " << "\n"
        << "CRT layer ID: " << layid << "\n"
        << "CRT distToReadout: " << distToReadout << ", distToReadout2: " << distToReadout2 << '\n'
        << "CRT abs0: " << abs0 << " , abs1: " << abs1 << '\n'
        << "CRT npeExpected: " << npeExpected << " , npeExpected2: " << npeExpected2 << '\n'
        << "CRT npeExp0: " << npeExp0 << " , npeExp1: " << npeExp1 << " , npeExp0Dual: " << npeExp0Dual << '\n'
        << "CRT charge q0: " << q0 << ", q1: " << q1 << '\n'
        << "CRT timing: tTrue: " << tTrue << ", t0: " << t0 << ", t1: " << t1 << ", dt: " << util::absDiff(t0,t1) << '\n'
        << " recoT-trueT = " << t0-tTrue << std::endl; 
       
    }//for AuxDetIDEs 
  }//for AuxDetChannels

  if(fVerbose) std::cout << "outside of AD loop" << std::endl;

  // Apply coincidence trigger requirement
  std::unique_ptr<std::vector<icarus::crt::CRTData> > triggeredCRTHits(
      new std::vector<icarus::crt::CRTData>);

  int nmiss_lock_c=0, nmiss_lock_d=0, nmiss_lock_m=0;
  int nmiss_dead_c=0, nmiss_dead_d=0, nmiss_dead_m=0;
  int nmiss_opencoin_c = 0, nmiss_opencoin_d = 0;
  int nmiss_coin_c = 0;
  int nmiss_coin_d = 0;
  int nmiss_coin_m = 0;
  int event = 0;
  int nhit_m=0, nhit_c=0, nhit_d=0;
  int neve_m=0, neve_c=0, neve_d=0;

  // Loop over all FEBs (key for taggers) with a hit and check coincidence requirement.
  // For each FEB, find channel providing trigger and determine if
  //  other hits are in concidence with the trigger (keep) 
  //  or if hits occur during R/O (dead time) (lost)
  //  or if hits are part of a different event (keep for now)
  // First apply dead time correction, biasing effect if configured to do so.
  // Front-end logic: For CERN or DC modules require at least one hit in each X-X layer.
  if (fVerbose) std::cout << '\n' << "about to loop over taggers (size " << taggers.size() << " )" << std::endl;

  for (auto trg : taggers) {

      event = 0;
      icarus::crt::CRTChannelData *chanTrigData(nullptr);
      icarus::crt::CRTChannelData *chanTmpData(nullptr);
      std::set<int> trackNHold = {};
      std::set<int> layerNHold = {};
      std::pair<int,int> macPair = std::make_pair(trg.first,trg.first); //FEB-FEB validation pair
      std::pair<int,int> tpair; 
      bool minosPairFound = false;
      std::vector<icarus::crt::CRTChannelData> passingData;
      double ttrig=0.0, ttmp=0.0;
      size_t trigIndex = 0;

      //check "open" coincidence (just check if coincdence possible w/hit in both layers) 
      if (trg.second.type=='c' && fApplyCoincidenceC && trg.second.layerid.size()<2) {
          nmiss_opencoin_c++;
          continue;
      }
      if (trg.second.type=='d' && fApplyCoincidenceD && trg.second.layerid.size()<2) {
          nmiss_opencoin_d++;
          continue;
      }

      //time order ChannalData objects in this FEB by T0
      std::sort((trg.second.data).begin(),(trg.second.data).end(),TimeOrderCRTData);

      if (fVerbose) std::cout << "processing data for FEB " << trg.first << " with "
                              << trg.second.data.size() << " entries..." << '\n'
                              << "    type: " <<  trg.second.type << '\n'
                              << "    region: " <<  trg.second.reg << '\n'
                              << "    layerID: " << *trg.second.layerid.begin() << '\n' << std::endl;
 
      //FIX ME! DOESN'T WORK IF ONLY ONE DATA ENTRY!!!!!!
      //outer (primary) loop over all data products for this FEB
      for ( size_t i=0; i< trg.second.data.size(); i++ ) {

        //get data for earliest entry
        if(i==0) {
          chanTrigData = &(trg.second.data[0]);
          ttrig = chanTrigData->T0(); //time stamp [ns]
          trackNHold.insert(chanTrigData->Channel());
          layerNHold.insert(trg.second.chanlayer[chanTrigData->Channel()]);
          passingData.push_back(*chanTrigData);
          ttmp = ttrig;
        }

        else {
          chanTmpData = &(trg.second.data[i]);
          ttmp = chanTmpData->T0(); 
        }

        //check that time sorting works
        if ( ttmp < ttrig ) mf::LogError("CRT") << "SORTING OF DATA PRODUCTS FAILED!!!"<< "\n";

        //for C and D modules only and coin. enabled, if assumed trigger channel has no coincidence
        // set trigger channel to tmp channel and try again
        if ( layerNHold.size()==1 &&
             ( (trg.second.type=='c' && fApplyCoincidenceC && (ttmp-ttrig>fLayerCoincidenceWindowC || trg.second.data.size()==1 )) ||
               (trg.second.type=='d' && fApplyCoincidenceD && (ttmp-ttrig>fLayerCoincidenceWindowD || trg.second.data.size()==1 )) ) ) {
             trigIndex++;
             chanTrigData = &(trg.second.data[trigIndex]);
             i = trigIndex+1;
             ttrig = chanTrigData->T0();
             trackNHold.clear();
             layerNHold.clear();
             passingData.clear();
             trackNHold.insert(chanTrigData->Channel());
             layerNHold.insert(trg.second.chanlayer[chanTrigData->Channel()]);
             passingData.push_back(*chanTrigData);
             if(trg.second.type=='c') nmiss_coin_c++;
             if(trg.second.type=='d') nmiss_coin_d++;
             continue;
        }

        //check if coincidence condtion met
        //for c and d modules, just need time stamps within tagger obj
        //for m modules, need to check coincidence with other tagger objs
        if (trg.second.type=='m' && !minosPairFound && fApplyCoincidenceM) {
            for (auto trg2 : taggers) {

                if(fVerbose) std::cout << "checking coincidence with FEB " << trg2.first << std::endl;

                if( trg2.second.type!='m' || //is other mod 'm' type
                  trg.first == trg2.first || //other mod not same as this one
                  trg.second.reg != trg2.second.reg || //other mod is in same region
                  *trg2.second.layerid.begin() == *trg.second.layerid.begin()) //modules are in adjacent layers
                    continue;

                //time sort all data for this FEB (all times for this event)
                std::sort((trg2.second.data).begin(),(trg2.second.data).end(),TimeOrderCRTData);

                //find entry within coincidence window starting with this FEB's
                //triggering channel in coincidence candidate's FEB
                for ( size_t j=0; j< trg2.second.data.size(); j++ ) {
                    double t2tmp = trg2.second.data[j].T0(); //in us
		    if ( util::absDiff(t2tmp,ttrig) < fLayerCoincidenceWindowM) {
                           minosPairFound = true;
                           macPair = std::make_pair(trg.first,trg2.first);
                        if (macPair.first==macPair.second+50 || macPair.first==macPair.second+50)
                            mf::LogError("CRT") << "invalid macPair!!! " << macPair.first << ", " << macPair.second << "\n";
                               
                        tpair = std::make_pair(chanTrigData->Channel(),trg2.second.data[j].Channel());
                        break;
                    }
                }
                //we found a valid pair so move on to next step
                if (minosPairFound) {
                    if(fVerbose) std::cout << "MINOS pair found! (" << macPair.first 
                                           << "," << macPair.second << ")" << std::endl;
                    break;
                }
            }//inner loop over febs (taggers)

            //if no coincidence pairs found, reinitialize and move to next FEB
            if(!minosPairFound) {
                if(fVerbose) std::cout << "MINOS pair NOT found! Skipping to next FEB..." << std::endl;
                if(trg.second.data.size()==1) continue;
                trigIndex++;
                chanTrigData = &(trg.second.data[trigIndex]);
                i = trigIndex+1;
                ttrig = chanTrigData->T0();
                trackNHold.clear();
                layerNHold.clear();
                passingData.clear();
                trackNHold.insert(chanTrigData->Channel());
                layerNHold.insert(trg.second.chanlayer[chanTrigData->Channel()]);
                passingData.push_back(*chanTrigData);
                nmiss_coin_m++;
                continue;
            }
        }//if minos module and no pair yet found

        if(fVerbose) std::cout << "done checking coinceidence...moving on to latency effects..." << std::endl;

        int adctmp = 0;
        std::vector<int> combined_trackids;
        //currently assuming bias time is same as track and hold window (FIX ME!)
        if (i>0 && ((trg.second.type=='c' && ttmp < ttrig + fLayerCoincidenceWindowC) || 
            (trg.second.type=='d' && ttmp < ttrig + fLayerCoincidenceWindowD) ||
            (trg.second.type=='m' && ttmp < ttrig + fLayerCoincidenceWindowM)) )
        {

            //if channel not locked
            if ((trackNHold.insert(chanTmpData->Channel())).second) {

                passingData.push_back(*chanTmpData);
                if (layerNHold.insert(trg.second.chanlayer[chanTmpData->Channel()]).second
                   && trg.second.type != 'm') 
                {
                
                    tpair=std::make_pair(chanTrigData->Channel(),chanTmpData->Channel());

                    if (trg.second.type=='c'&&
                      ((tpair.first<16&&tpair.second<16)||(tpair.first>15&&tpair.second>15)) )
                        mf::LogError("CRT")<< "incorrect CERN trigger pair!!!" << '\n'
                                           << "  " << tpair.first << ", " << tpair.second << "\n";
                }
            } //end if channel not locked

            else if (ttmp < ttrig + fBiasTime) {
                adctmp = (passingData.back()).ADC();
                adctmp += chanTmpData->ADC();
                (passingData.back()).SetADC(adctmp);

                for ( auto const& ids : passingData.back().TrackID())
                    combined_trackids.push_back(ids);
                combined_trackids.push_back(chanTmpData->TrackID()[0]);
                passingData.back().SetTrackID(combined_trackids);
            } //channel is locked but hit close in time to bias pulse height

            else switch (trg.second.type) {
                case 'c' : nmiss_lock_c++; break;
                case 'd' : nmiss_lock_d++; break;
                case 'm' : nmiss_lock_m++; break;
            } //data is discarded (not writeable)

        }//if hits inside readout window

        else if ( i>0 && ttmp <= ttrig + fDeadTime ) {
            switch (trg.second.type) {
                case 'c' : nmiss_dead_c++; break;
                case 'd' : nmiss_dead_d++; break;
                case 'm' : nmiss_dead_m++; break;
            }
              //continue;
        } // hits occuring during digitization lost (dead time)

        //"read out" data for this event, first hit after dead time as next trigger channel
        else if ( ttmp > ttrig + fDeadTime) {// || i==trg.second.data.size()-1) {

            int regnum = GetAuxDetRegionNum(trg.second.reg);
            if( (regions.insert(regnum)).second) regCounts[regnum] = 1;
            else regCounts[regnum]++;
            if (fVerbose) std::cout << "creating CRTData product just after deadtime" << std::endl;
            triggeredCRTHits->push_back(
              icarus::crt::CRTData(eve, trg.first,event,ttrig,chanTrigData->Channel(),tpair,macPair,passingData));
            if (fVerbose) std::cout << " ...success!" << std::endl;
            event++;
            if (trg.second.type=='c') {neve_c++; nhit_c+=passingData.size(); }
            if (trg.second.type=='d') {neve_d++; nhit_d+=passingData.size(); }
            if (trg.second.type=='m') {neve_m++; nhit_m+=passingData.size(); }
            ttrig = ttmp;
            chanTrigData = chanTmpData;
            passingData.clear();
            trackNHold.clear();
            layerNHold.clear();
            passingData.push_back(*chanTrigData);
            trackNHold.insert(chanTrigData->Channel());
            layerNHold.insert(trg.second.chanlayer[chanTmpData->Channel()]);
            minosPairFound = false;
        }

        if (!(ttmp > ttrig + fDeadTime) && i==trg.second.data.size()-1) {
            if (fVerbose) std::cout << "creating CRTData product at end of FEB events..." << std::endl;
            triggeredCRTHits->push_back(
              icarus::crt::CRTData(eve,trg.first,event,ttrig,chanTrigData->Channel(),tpair,macPair,passingData));
            if (fVerbose) std::cout << " ...success!" << std::endl;
            event++;
            if (trg.second.type=='c') {neve_c++; nhit_c+=passingData.size(); }
            if (trg.second.type=='d') {neve_d++; nhit_d+=passingData.size(); }
            if (trg.second.type=='m') {neve_m++; nhit_m+=passingData.size(); }

        } 
      }//for data entries (hits)

      if(fVerbose) std::cout << " outside loop over FEB data entries...moving on to next FEB..." << std::endl;

  } // for taggers

  if (fVerbose) { 
    std::cout << "|---------------------- CRTDetSim Summary ----------------------|" << '\n'
     << " - - - EDeps from AuxDetIDE's - - -" << '\n'
     << "  CERN  sim hits: " << nsim_c << '\n'
     << "  MINOS sim hits: " << nsim_m << '\n'
     << "  DC    sim hits: " << nsim_d << '\n' << '\n'
     << " - - - Single Channel Threshold - - -" << '\n'
     << " Pass:" << '\n'
     << "   CERN hits  > thresh: " << nchandat_c << '\n'
     << "   MINOS hits > thresh: " << nchandat_m << '\n'
     << "   DC hits    > thresh: " << nchandat_d << '\n'
     << " Lost:" << '\n'
     << "   CERN hits  < thresh: " << nmissthr_c << '\n'
     << "   MINOS hits < thresh: " << nmissthr_m << '\n'
     << "   DC hits    < thresh: " << nmissthr_d << '\n' << '\n'
     << " - - - System Specifc Coincidence Loss- - -" << '\n'
     << "  CERN      fiber-fiber: " << nmiss_strcoin_c << '\n'
     << "  CERN open layer-layer: " << nmiss_opencoin_c << '\n'
     << "  DC   open layer-layer: " << nmiss_opencoin_d << '\n'
     << "  CERN      layer-layer: " << nmiss_coin_c << " (" << 100.0*nmiss_coin_c/nchandat_c << "%)" << '\n'
     << "  MINOS     layer-layer: " << nmiss_coin_m << " (" << 100.0*nmiss_coin_m/nchandat_m << "%)" << '\n'
     << "  DC        layer-layer: " << nmiss_coin_d << " (" << 100.0*nmiss_coin_d/nchandat_d << "%)" << '\n' << '\n'
     << " - - - Front-End Electronics Effects Losses - - -" << '\n'
     << "  CERN  trackNHold: " << nmiss_lock_c << " (" << 100.0*nmiss_lock_c/nchandat_c << "%)" << '\n' 
     << "  MINOS trackNHold: " << nmiss_lock_m << " (" << 100.0*nmiss_lock_m/nchandat_m << "%)" << '\n'
     << "  DC    trackNHold: " << nmiss_lock_d << " (" << 100.0*nmiss_lock_d/nchandat_d << "%)" << '\n'
     << "  CERN    deadTime: " << nmiss_dead_c << " (" << 100.0*nmiss_dead_c/nchandat_c << "%)" << '\n'
     << "  MINOS   deadTime: " << nmiss_dead_m << " (" << 100.0*nmiss_dead_m/nchandat_m << "%)" << '\n'
     << "  DC      deadTime: " << nmiss_dead_d << " (" << 100.0*nmiss_dead_d/nchandat_d << "%)" << '\n' << '\n'
     << " - - - Passing Hits Pushed To Event - - -" << '\n'
     << "  CERN    hits: " << nhit_c << " (" << 100.0*nhit_c/nchandat_c << "%)" << '\n'
     << "  MINOS   hits: " << nhit_m << " (" << 100.0*nhit_m/nchandat_m << "%)" << '\n'
     << "  DC      hits: " << nhit_d << " (" << 100.0*nhit_d/nchandat_d << "%)" << '\n'
     << "  CERN  events: " << neve_c << '\n'
     << "  MINOS events: " << neve_m << '\n'
     << "  DC    events: " << neve_d << '\n' 
     << "  Total pushes to event: " << triggeredCRTHits->size() << std::endl;
     
    std::map<int,int>::iterator it = regCounts.begin();
    std::cout << '\n' << "FEB events per CRT region: " << '\n' << std::endl;
     
    while ( it != regCounts.end() ) {
        std::cout << "reg: " << (*it).first << " , events: " << (*it).second << '\n' << std::endl;
        it++;
    }
  } //if verbose

  e.put(std::move(triggeredCRTHits));
}

DEFINE_ART_MODULE(CRTDetSim)

}  // namespace crt
}  // namespace icarus

