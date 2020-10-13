/**
 * @file   CRTDataAnalysis_module.cc
 * @brief  Access CRT data and reco products and compare to MCTruth info 
 * @author Chris Hilgenberg (Chris.Hilgenberg@colostate.edu)
 * 
 * The last revision of this code was done in October 2018 with LArSoft v07_06_01.
 */

// LArSoft includes
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/AuxDetGeometryCore.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/Exception.h"

// Utility libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"
#include "cetlib/pow.h" // cet::sum_of_squares()

// ROOT includes
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TGeoManager.h"
#include "TMath.h"
#include "TROOT.h"

// C++ includes
#include <map>
#include <vector>
#include <string>
#include <set>
#include <cmath>
#include <iostream>
#include <utility>
#include <array>

// CRT data products
#include "sbnobj/ICARUS/CRT/CRTData.hh"
#include "sbnobj/Common/CRT/CRTHit.hh"

using std::string;
using std::vector;
using std::map;
using std::set;
using std::pair;


namespace {
  //uint32_t ModToTypeCode(geo::AuxDetGeo const& adgeo); 
  //char ModToAuxDetType(geo::AuxDetGeo const& adgeo);
  //int GetAuxDetRegion(geo::AuxDetGeo const& adgeo);
  uint32_t MacToADReg(uint32_t mac);
  char MacToType(uint32_t mac);
  uint32_t MacToTypeCode(uint32_t mac);
  //std::pair<uint32_t,uint32_t> ADToMac(const map<int,vector<pair<int,int>>>& febMap, uint32_t adid);
  int MacToAuxDetID(const map<int,vector<pair<int,int>>>& febMap, int mac, int chan);
  int ChannelToAuxDetSensitiveID(int mac, int chan);
  int RegToTypeCode(int reg);

} // local namespace


namespace icarus {
namespace crt {

  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------
  // class definition
  /**
   * @brief Example analyzer
   * 
   * This class extracts information from the generated and reconstructed
   * particles.
   *
   * It produces histograms for the simulated particles in the input file:
   * - PDG ID (flavor) of all particles
   * - momentum of the primary particles selected to have a specific PDG ID
   * - length of the selected particle trajectory
   * 
   * It also produces two ROOT trees.
   *
   * The first ROOT tree contains information on the simulated
   * particles, including "dEdx", a binned histogram of collected
   * charge as function of track range.
   * 
   * Configuration parameters
   * =========================
   * 
   * - *SimulationLabel* (string, default: "largeant"): tag of the input data
   *   product with the detector simulation information (typically an instance
   *   of the LArG4 module)
   *
   */
  class CRTDataAnalysis : public art::EDAnalyzer
  {
  public:

    struct Config {
      
      // Save some typing:
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      

      fhicl::Atom<art::InputTag> CRTHitLabel {
        Name("CRTHitLabel"),
        Comment("tag of the input data product with reconstructed CRT hits")
        };

      fhicl::Atom<art::InputTag> CRTDAQLabel {
        Name("CRTDAQLabel"),
        Comment("tag of the input data product with calibrated CRT data")
        };
   

    }; // Config
    
    using Parameters = art::EDAnalyzer::Table<Config>;
    
    // -------------------------------------------------------------------
    // -------------------------------------------------------------------
    // Standard constructor for an ART module with configuration validation;
    // we don't need a special destructor here.

    /// Constructor: configures the module (see the Config structure above)
    explicit CRTDataAnalysis(Parameters const& config);

    virtual void beginJob() override;
    virtual void beginRun(const art::Run& run) override;
    virtual void analyze (const art::Event& event) override;

  private:

    void FillFebMap();

    // The parameters we'll read from the .fcl file.
    art::InputTag fCRTHitProducerLabel;        ///< The name of the producer that created hits
    art::InputTag fCRTDAQProducerLabel;

    static map<int, vector<pair<int,int>>> fFebMap;

    // The n-tuples we'll create.
    TTree* fDAQNtuple;
    TTree* fHitNtuple;

    // The comment lines with the @ symbols define groups in doxygen. 
    /// @name The variables that will go into both n-tuples.
    /// @{
    int fEvent;        ///< number of the event being processed
    int fRun;          ///< number of the run being processed
    int fSubRun;       ///< number of the sub-run being processed
    /// @}

    /// @name The variables that will go into the CosmicDisplay n-tuple.
    /// @{
    static const int LAR_PROP_DELAY = 1.0/(30.0/1.38); //[ns/cm]

    //CRT data product vars
    //static const int kDetMax = 64;
    int      fDetEvent;
    int      fNChan; ///< number of channels above threshold for this front-end board readout
    int      fEntry; ///< front-end board entry number (reset for each event)
    int      fFEBReg; ///< CRT region for this front-end board
    int      fMac5; ///< Mac5 address for this front-end board
    int      fDetSubSys;
    double   fT0;///< signal time w.r.t. global event time
    double   fT1;///< signal time w.r.t. PPS
    int      fADC[64];///< signal amplitude
    vector<vector<int>> fTrackID;///< track ID(s) of particle that produced the signal
    vector<vector<int>> fDetPDG; /// signal inducing particle(s)' PDG code

    //CRT hit product vars
    int       fHitEvent;
    float    fXHit; ///< reconstructed X position of CRT hit (cm)
    float    fYHit; ///< reconstructed Y position of CRT hit (cm)
    float    fZHit; ///< reconstructed Z position of CRT hit (cm)
    float    fXErrHit; ///< stat error of CRT hit reco X (cm)
    float    fYErrHit; ///< stat error of CRT hit reco Y (cm)
    float    fZErrHit; ///< stat error of CRT hit reco Z (cm)
    int32_t    fT0Hit; ///< hit time w.r.t. global event time
    int32_t    fT1Hit; ///< hit time w.r.t. PPS
    //double    fT0CorrHit;
    //double    fT1CorrHit;
    int       fHitReg; ///< region code of CRT hit
    int       fHitSubSys;
    int       fNHit; ///< number of CRT hits for this event
    int       fHitStrip;
    int       fHitMod;

    TH1F* fModMultHistC;   ///< true N C-modules hit / muon track
    TH1F* fModMultHistM;   ///< true N M-modules hit / muon track
    TH1F* fModMultHistD;   ///< true N D-modules hit / muon track

    // note that the following four variables are not used so are being commented out here
    //uint32_t fNmuTagC;     //N muon tracks producing >0 CRT triggers in C-subsystem
    //uint32_t fNmuTagM;     //N muon tracks producing >0 CRT triggers in M-subsystem
    //uint32_t fNmuTagD;     //N muon tracks producing >0 CRT triggers in D-subsystem
    //uint32_t fNmuTagTot;
    TH1F* fChanMultHistC;  //N FEB channels > threshold / muon track
    TH1F* fChanMultHistM;  
    TH1F* fChanMultHistD;
    TH1F* fFEBMultHistC;   //N FEBs w/trigger / muon track
    TH1F* fFEBMultHistM;
    TH1F* fFEBMultHistD;

    /// @}
    
    // Other variables that will be shared between different methods.
    geo::GeometryCore const* fGeometryService;   ///< pointer to Geometry provider
    int                      fTriggerOffset;     ///< (units of ticks) time of expected neutrino event
    
  }; // class CRTDataAnalysis


  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------
  // class implementation

  //-----------------------------------------------------------------------
  // Constructor
  // 
  // Note that config is a Table<Config>, and to access the Config
  // value we need to use an operator: "config()". In the same way,
  // each element in Config is an Atom<Type>, so to access the type we
  // again use the call operator, e.g. "SimulationLabel()".

  map<int,vector<pair<int,int>>> CRTDataAnalysis::fFebMap;
 
  CRTDataAnalysis::CRTDataAnalysis(Parameters const& config)
    : EDAnalyzer(config)
    , fCRTHitProducerLabel(config().CRTHitLabel())
    , fCRTDAQProducerLabel(config().CRTDAQLabel())
  {
    // Get a pointer to the geometry service provider.
    fGeometryService = lar::providerFrom<geo::Geometry>();
    // The same for detector TDC clock services.
    // Access to detector properties.
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
    fTriggerOffset = trigger_offset(clockData);
  }

  void CRTDataAnalysis::FillFebMap() { 
    if(!this->fFebMap.empty())
        return;
    std::string fullFileName;
    cet::search_path searchPath("FW_SEARCH_PATH");
    searchPath.find_file("feb_map.txt",fullFileName);
    std::ifstream fin;
    fin.open(fullFileName,std::ios::in);
    if(fin.good()) std::cout << "opened file 'feb_map.txt' for reading..." << std::endl;
    else //std::cout << "could not open file 'feb_map.txt' for reading!" << std::endl;
        throw cet::exception("CRTDataAnalysis::FillFebMap") << "Unable to find/open file 'feb_map.txt'" << std::endl;
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
        (this->fFebMap)[mod].push_back(std::make_pair(std::stoi(row[1]),std::stoi(row[2])));
        if(row.size()>3)
            (this->fFebMap)[mod].push_back(std::make_pair(std::stoi(row[3]),std::stoi(row[4])));
    }
    std::cout << "filled febMap with " << (this->fFebMap).size() << " entries" << std::endl;
    fin.close();
  }

  
  //-----------------------------------------------------------------------
  void CRTDataAnalysis::beginJob()
  {
    std::cout << " starting analysis job" << std::endl;

    // Access ART's TFileService, which will handle creating and writing
    // histograms and n-tuples for us. 
    art::ServiceHandle<art::TFileService> tfs;

    // Define our n-tuples
    fDAQNtuple        = tfs->make<TTree>("DAQTree",          "MyCRTDAQ");
    fHitNtuple        = tfs->make<TTree>("HitTree",          "MyCRTHit");

    // Construct truth matching histograms
    //fStripMultHistC   = tfs->make<TH1F>("StripMultC",";no. strips hit / module / #mu;",64,0,64);
    //fStripMultHistM   = tfs->make<TH1F>("StripMultM",";no. strips hit / module / #mu;",64,0,64);
    //fStripMultHistD   = tfs->make<TH1F>("StripMultD",";no. strips hit / module / #mu;",64,0,64);
    fModMultHistC     = tfs->make<TH1F>("ModMultC",";no. modules hit / #mu;",10,0,10);
    fModMultHistM     = tfs->make<TH1F>("ModMultM",";no. modules hit / #mu;",10,0,10);
    fModMultHistD     = tfs->make<TH1F>("ModMultD",";no. modules hit / #mu;",10,0,10);

    fChanMultHistC    =tfs->make<TH1F>("ChanMultC",";no. FEB channels > threshold / FEB / #mu;",64,0,64);
    fChanMultHistM    =tfs->make<TH1F>("ChanMultD",";no. FEB channels > threshold / FEB / #mu;",64,0,64);
    fChanMultHistD    =tfs->make<TH1F>("ChanMultM",";no. FEB channels > threshold / FEB / #mu;",64,0,64);
    fFEBMultHistC     =tfs->make<TH1F>("FEBMultC",";no. FEB triggers / #mu;",64,0,64);
    fFEBMultHistM     =tfs->make<TH1F>("FEBMultD",";no. FEB triggers / #mu;",64,0,64);
    fFEBMultHistD     =tfs->make<TH1F>("FEBMultM",";no. FEB triggers / #mu;",64,0,64);

    // Define the branches of our DetSim n-tuple 
    fDAQNtuple->Branch("event",                 &fDetEvent,          "event/I");
    fDAQNtuple->Branch("nChan",                 &fNChan,             "nChan/I");
    fDAQNtuple->Branch("t0",                    &fT0,                "t0/D");
    fDAQNtuple->Branch("t1",                    &fT1,                "t1/D");
    fDAQNtuple->Branch("adc",                   fADC);
    fDAQNtuple->Branch("entry",                 &fEntry,             "entry/I");
    fDAQNtuple->Branch("mac5",                  &fMac5,              "mac5/I");
    fDAQNtuple->Branch("region",                &fFEBReg,            "region/I");
    fDAQNtuple->Branch("subSys",                &fDetSubSys,         "subSys/I");

    // Define the branches of our SimHit n-tuple
    fHitNtuple->Branch("event",       &fHitEvent,    "event/I");
    fHitNtuple->Branch("nHit",        &fNHit,        "nHit/I");
    fHitNtuple->Branch("x",           &fXHit,        "x/F");
    fHitNtuple->Branch("y",           &fYHit,        "y/F");
    fHitNtuple->Branch("z",           &fZHit,        "z/F");
    fHitNtuple->Branch("xErr",        &fXErrHit,     "xErr/F");
    fHitNtuple->Branch("yErr",        &fYErrHit,     "yErr/F");
    fHitNtuple->Branch("zErr",        &fZErrHit,     "zErr/F");
    fHitNtuple->Branch("t0",          &fT0Hit,       "t0/I");
    fHitNtuple->Branch("t1",          &fT1Hit,       "t1/I");
    fHitNtuple->Branch("region",      &fHitReg,      "region/I");  
    fHitNtuple->Branch("subSys",      &fHitSubSys,   "subSys/I");
    fHitNtuple->Branch("modID",       &fHitMod,      "modID/I");
    fHitNtuple->Branch("stripID",     &fHitStrip,    "stripID/I");

}
   
  void CRTDataAnalysis::beginRun(const art::Run& /*run*/)
  {
  }

  //-----------------------------------------------------------------------
  void CRTDataAnalysis::analyze(const art::Event& event) 
  {
    MF_LOG_DEBUG("CRT") << "beginning analyis" << '\n';

    // Start by fetching some basic event information for our n-tuple.
    fEvent  = event.id().event(); 
    fRun    = event.run();
    fSubRun = event.subRun();

    FillFebMap();//febMap);

    art::Handle<vector<icarus::crt::CRTData>> crtDAQHandle;
    bool isCRTDAQ = event.getByLabel(fCRTDAQProducerLabel, crtDAQHandle);

    if (isCRTDAQ)  {
     std::cout << "about to loop over detsim entries" << std::endl;
     for ( auto const& febdat : (*crtDAQHandle) ) {
        fDetEvent       = fEvent;
        fMac5           = febdat.fMac5;
        fEntry          = febdat.fEntry;
        fFEBReg         = MacToADReg(fMac5);
        fNChan = 0;
        fDetSubSys = MacToTypeCode(fMac5);
        fT0 = febdat.fTs0;
        fT1 = febdat.fTs1;

        int maxchan =0;
        if(fDetSubSys!=2) maxchan=32;
        else maxchan = 64;
        for(int ch=0; ch<maxchan; ch++) {
            fADC[ch] = febdat.fAdc[ch]; 
        } 


        fDAQNtuple->Fill();

     } //for CRT FEB events

    }//if crtdetsim products present

    else 
	throw cet::exception("CRTSimAnalysis") << "CRTDAQ products not found!" << std::endl;

    art::Handle<std::vector<sbn::crt::CRTHit>> crtHitHandle;
    
    bool isCRTHit = event.getByLabel(fCRTHitProducerLabel, crtHitHandle);
    std::vector<int> ids;
    fNHit = 0;
    if (isCRTHit) {

        //art::fill_ptr_vector(crtSimHits,crtSimHitHandle);
        art::FindManyP<icarus::crt::CRTData> findManyData(crtHitHandle, event, fCRTHitProducerLabel);
        std::vector<art::Ptr<icarus::crt::CRTData>> data = findManyData.at(0);

        std::cout << "looping over sim hits..." << std::endl;
        for ( auto const& hit : *crtHitHandle )
        {
            fNHit++;
            fHitEvent = fEvent;
            fXHit    = hit.x_pos;
            fYHit    = hit.y_pos;
            fZHit    = hit.z_pos;
            fXErrHit = hit.x_err;
            fYErrHit = hit.y_err;
            fZErrHit = hit.z_err;
            fT0Hit   = hit.ts0_ns;
            fT1Hit   = hit.ts1_ns;

            int mactmp = hit.feb_id[0];
            fHitReg  = MacToADReg(mactmp);
            fHitSubSys = RegToTypeCode(fHitReg);

            auto ittmp = hit.pesmap.find(mactmp);
            if (ittmp==hit.pesmap.end()) {
                std::cout << "hitreg: " << fHitReg << std::endl;
                std::cout << "fHitSubSys: "<< fHitSubSys << std::endl;
                std::cout << "mactmp = " << mactmp << std::endl;
                std::cout << "could not find mac in pesmap!" << std::endl;
                continue;
            }
            int chantmp = (*ittmp).second[0].first;

            fHitMod  = MacToAuxDetID(this->fFebMap, mactmp, chantmp);
            fHitStrip = ChannelToAuxDetSensitiveID(mactmp, chantmp);

            fHitNtuple->Fill();
        }//for CRT Hits
    }//if CRT Hits

    else std::cout << "CRTHit products not found! (expected if decoder step)" << std::endl;


  } // CRTDataAnalysis::analyze()
  
  
  DEFINE_ART_MODULE(CRTDataAnalysis)

} // namespace crt
} // namespace icarus


// Back to our local namespace.
namespace {

  /*char ModToAuxDetType(geo::AuxDetGeo const& adgeo) {
    size_t nstrips = adgeo.NSensitiveVolume();
    if (nstrips==16) return 'c'; 
    if (nstrips==20) return 'm';
    if (nstrips==64) return 'd';
    return 'e';
  }

  uint32_t ModToTypeCode(geo::AuxDetGeo const& adgeo) {
    size_t nstrips = adgeo.NSensitiveVolume();
    if (nstrips==16) return 0; //'c'
    if (nstrips==20) return 1; //'m'
    if (nstrips==64) return 2; //'d'
    return UINT32_MAX;
  }*/

  /*int GetAuxDetRegion(geo::AuxDetGeo const& adgeo)
  {
    char type = ModToAuxDetType(adgeo);
    std::string base = "volAuxDet_";
    switch ( type ) {
      case 'c' : base+= "CERN"; break;
      case 'd' : base+= "DC"; break;
      case 'm' : base+= "MINOS"; break;
      case 'e' :
          std::cout << "error in GetAuxDetRegion: type error" << std::endl;
          return INT_MAX;
    }
    base+="_module_###_";
    std::string volName(adgeo.TotalVolume()->GetName());

    //module name has 2 possible formats
    //  volAuxDet_<subsystem>_module_###_<region>
    //  volAuxDet_<subsystem>_module_###_cut###_<region>

    std::string reg = volName.substr(base.length(),volName.length());
    if( reg.find("_")!=std::string::npos)
        reg = reg.substr(reg.find("_")+1,reg.length());

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

    return INT_MAX;
  }*/

  uint32_t MacToADReg(uint32_t mac) {

      if(mac>=107 && mac<=190) return 30; //top
      if(mac>=191 && mac<=204) return 31; //rim west
      if(mac>=205 && mac<=218) return 32; //rim east
      if(mac>=219 && mac<=224) return 33; //rim south
      if(mac>=225 && mac<=230) return 34; //rim north
      if(            mac<=12 ) return 40; //west side, south stack
      if(mac>=13  && mac<=24 ) return 41; //west side, center stack
      if(mac>=25  && mac<=36 ) return 42; //west side, north stack
      if(mac>=37  && mac<=48 ) return 43; //east side, south stack
      if(mac>=49  && mac<=60 ) return 44; //east side, center stack
      if(mac>=61  && mac<=72 ) return 45; //east side, north stack
      if(mac>=73  && mac<=84 ) return 46; //south
      if(mac>=85  && mac<=92 ) return 47; //north
      if(mac>=93 && mac<=106) return 50; //bottom

      return 0;
  }

  char MacToType(uint32_t mac) {

      uint32_t reg = MacToADReg(mac);

      if( reg>29 && reg<40 )
        return 'c';
      if( reg>39 && reg<50 )
        return 'm';
      if( reg==50 )
        return 'd';

      return 'e';
  }

  uint32_t MacToTypeCode(uint32_t mac) {

      char reg = MacToType(mac);

      if(reg=='c')
        return 0; //'c';
      if(reg=='m')
        return 1; //'m';
      if(reg=='d')
        return 2; //'d';

      return UINT32_MAX;//'e';
  }

  //for C- and D-modules, mac address is same as AD ID
  //three M-modules / FEB, each modules readout at both ends
  //  numbering convention is module from FEB i 
  //  is readout on the opposite end by FEB i+50
  //  return FEB i
  /*std::pair<uint32_t,uint32_t> ADToMac(const map<int,vector<pair<int,int>>>& febMap, uint32_t adid) {
      for(auto const& p : febMap) {
          if((uint32_t)p.first!=adid)
              continue;
          if(p.second.size()==2)
              return std::make_pair((uint32_t)p.second[0].first,(uint32_t)p.second[1].first);
          else
              return std::make_pair((uint32_t)p.second[0].first,(uint32_t)p.second[0].first);
      }
      return std::make_pair(UINT32_MAX,UINT32_MAX);
  }*/

  int MacToAuxDetID(const map<int,vector<pair<int,int>>>& febMap, int mac, int chan){
      char type = MacToType(mac);
      int pos=1;
      if(type=='m')
          pos = chan/10 + 1;

      for(auto const& p : febMap) {
          if(p.second[0].first == mac && p.second[0].second==pos)
              return (uint32_t)p.first;
          if(p.second.size()==2)
              if(p.second[1].first==mac && p.second[1].second==pos)
                  return (uint32_t)p.first;
      }

    return INT_MAX;
  }

  int ChannelToAuxDetSensitiveID(int mac, int chan) {
    int type = MacToTypeCode(mac);
    if (type==2) return chan; //d
    if (type==0) return chan/2; //c
    if (type==1) return (chan % 10)*2; //m

    return INT_MAX;
  }

  int RegToTypeCode(int reg){
      if(reg>=30&&reg<40)
          return 0;
      if(reg>=40&&reg<50)
          return 1;
      if(reg==50)
          return 2;
      std::cout << "ERROR in RegToTypeCode: unknown reg code!" << std::endl;
      return -1;
  }

}//local namespace
