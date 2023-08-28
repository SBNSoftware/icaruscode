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
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/AuxDetGeometryCore.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "sbnobj/Common/Trigger/ExtraTriggerInfo.h"
//#include "icaruscode/CRT/CRTUtils/CRTHitRecoAlg.h"

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
#include "icaruscode/CRT/CRTUtils/CRTCommonUtils.h"
#include "icaruscode/CRT/CRTUtils/CRTPMTMatchingUtils.h"
#include "icaruscode/Decode/DecoderTools/IDecoder.h"
#include "icaruscode/Decode/ChannelMapping/IICARUSChannelMap.h"

using std::string;
using std::vector;
using std::map;
using std::set;
using std::pair;

using namespace sbn::crt;

namespace icarus {
namespace crt {

  //-----------------------------------------------------------------------

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
   
      fhicl::Atom<art::InputTag> TriggerLabel {
        Name("TriggerLabel"),
	  Comment("Label for the Trigger fragment label")
	  };
      fhicl::Atom<art::InputTag> CRTPMTLabel {
	Name("CRTPMTLabel"),
	  Comment("Label for the CRTPMT Matched variables from the crtpmt data product")
	  };
      fhicl::Atom<double> QPed {
	Name("QPed"),
	  Comment("Pedestal offset [ADC]")
	  };
      fhicl::Atom<double> QSlope {
	Name("QSlope"),
	  Comment("Pedestal slope [ADC/photon]")
	  };

      fhicl::Atom<double> PEThresh {
	Name("PEThresh"),
	  Comment("threshold in photoelectrons above which charge amplitudes used in hit reco")
	  };

      fhicl::Atom<uint64_t> CrtWindow {
	Name("CrtWindow"),
	  Comment("window for looking data [ns]")
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

    //void reconfigure(fhicl::ParameterSet const & p);
  private:

    void FillFebMap();
    
    // Declare member data here.
    const icarusDB::IICARUSChannelMap* fChannelMap = nullptr;
    //    CRTHitRecoAlg hitAlg;
    void ClearVecs();
    // The parameters we'll read from the .fcl file.
    art::InputTag fTriggerLabel;
    art::InputTag fCRTHitProducerLabel;        ///< The name of the producer that created hits
    art::InputTag fCRTDAQProducerLabel;
    art::InputTag fCRTPMTProducerLabel;

    //    bool fVerbose;          ///< print info
    double fQPed;           ///< Pedestal offset of SiPMs [ADC]
    double fQSlope;         ///< Pedestal slope of SiPMs [ADC/photon]
    double fPEThresh;       ///< threshold[PE] above which charge amplitudes used in hit reco
    uint64_t fCrtWindow;    ///< Looking data window within trigger timestamp [ns]

    static map<int, vector<pair<int,int>>> fFebMap;

    // The n-tuples we'll create.
    TTree* fDAQNtuple;
    TTree* fHitNtuple;
    TTree* fCRTPMTNtuple;

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

    //add trigger data product vars
    unsigned int m_gate_type;
    std::string m_gate_name;
    uint64_t m_trigger_timestamp;
    uint64_t m_gate_start_timestamp;
    uint64_t m_trigger_gate_diff;
    uint64_t m_gate_crt_diff;
    uint64_t m_crt_global_trigger;
    Long64_t m_crtGT_trig_diff;
    //CRT data product vars
    //static const int kDetMax = 64;
    int      fDetEvent;
    int      fDetRun;
    int      fDetSubRun;
    int      fNChan; ///< number of channels above threshold for this front-end board readout
    int      fEntry; ///< front-end board entry number (reset for each event)
    int      fFEBReg; ///< CRT region for this front-end board
    int      fMac5; ///< Mac5 address for this front-end board
    int      fDetSubSys;
    uint64_t fT0;///< signal time w.r.t. PPS
    uint64_t fT1;///< signal time w.r.t. global event time
    int      fNMaxCh;/// Max number of channel
    int      fADC[64];///< signal amplitude
    float    fPE[64];///< signal amplitude
    int      fFlags;///< Flags
    vector<vector<int>> fTrackID;///< track ID(s) of particle that produced the signal
    vector<vector<int>> fDetPDG; /// signal inducing particle(s)' PDG code

    //CRT hit product vars
    int      fHitEvent;
    float    fXHit; ///< reconstructed X position of CRT hit (cm)
    float    fYHit; ///< reconstructed Y position of CRT hit (cm)
    float    fZHit; ///< reconstructed Z position of CRT hit (cm)
    float    fXErrHit; ///< stat error of CRT hit reco X (cm)
    float    fYErrHit; ///< stat error of CRT hit reco Y (cm)
    float    fZErrHit; ///< stat error of CRT hit reco Z (cm)
    uint64_t    fT0Hit; ///< hit time w.r.t. PPS
    Long64_t    fT1Hit; ///< hit time w.r.t. global trigger
    int       fHitReg; ///< region code of CRT hit
    int       fHitSubSys;
    int       fNHit; ///< number of CRT hits for this event
    int       fHitStrip;
    int       fHitMod;
    int       fNHitFeb;
    float     fHitTotPe;

    //CRT-PMT Matching vars
    int          fMatchEvent;///< Event number.
    int          fMatchRun;///< Run number.
    unsigned int fGateType;///< Beam gate type.
    int          fFlashID; ///< ID of the optical flash.
    double       fFlashTime_us;///< Time of the optical flash w.r.t. the global trigger in us.    
    double fFlashGateTime_ns;///< Time of the optical flash w.r.t. the beam gate opening in ns.
    double fFirstOpHitPeakTime;///< Time of the first optical hit peak time w.r.t. the global trigger [us]
    double fFirstOpHitStartTime; ///< Time of the first optical hit start time w.r.t. the global trigger [us]
    bool fFlashInGate;///< Flash within gate or not.
    bool fFlashInBeam;///< Flash within the beam window of the gate or not.
    double fFlashPE;///< Total reconstructed light in the flash [photoelectrons]
    double       fFlashYWidth;///< Flash spread along Y.
    double       fFlashZWidth;///< Flash spread along Z.
    double fFlashPos_x;///< Flash barycenter coordinates evaluated using ADCs as weights, X-position.
    double fFlashPos_y;///< Flash barycenter coordinates evaluated using ADCs as weights, Y-position.
    double fFlashPos_z;///< Flash barycenter coordinates evaluated using ADCs as weights, Z-position.
    MatchType fFlashClassification;///< Classication of the optical flash.
    std::vector<MatchedCRT> matchedCRTHits;///< Matched CRT Hits with the optical flash.
    // add contents of MatchedCRT struct to be put into branches, 
    //geo::Point_t CRTHitPos;
    int nMatchedCRTHits; ///< Number of Matched CRT hits to flash 
    vector<double> CRTHitPos_x;
    vector<double> CRTHitPos_y;
    vector<double> CRTHitPos_z;
    vector<double> fCRTPMTTimeDiff_ns;
    vector<double> fCRTTime_us;
    vector<int> fCRTSys;
    vector<int> fCRTRegion;     

    int fNtopCRTBefore;
    int fNtopCRTAfter;
    int fNsideCRTBefore;
    int fNsideCRTAfter;
    //std::vector<recob::OpHit>opHits;///< Optical hits of the flash.
    
    // Other variables that will be shared between different methods.
    geo::GeometryCore const* fGeometryService;   ///< pointer to Geometry provider
    int                      fTriggerOffset;     ///< (units of ticks) time of expected neutrino event
    CRTCommonUtils* fCrtutils;  
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
    , fTriggerLabel( config().TriggerLabel() )
    , fCRTHitProducerLabel(config().CRTHitLabel())
    , fCRTDAQProducerLabel(config().CRTDAQLabel())
    , fCRTPMTProducerLabel(config().CRTPMTLabel())
    , fQPed(config().QPed())
    , fQSlope(config().QSlope())
    , fPEThresh(config().PEThresh())
    , fCrtWindow(config().CrtWindow())
    , fCrtutils(new CRTCommonUtils())
      
  {
    // Get a pointer to the geometry service provider.
    fGeometryService = lar::providerFrom<geo::Geometry>();
    fChannelMap = art::ServiceHandle<icarusDB::IICARUSChannelMap const>{}.get();
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
    if(fin.good())  mf::LogError("CRTDataAnalysis") << "opened file 'feb_map.txt' for reading..." << std::endl;
    else // mf::LogError("CRTDataAnalysis") << "could not open file 'feb_map.txt' for reading!" << std::endl;
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
    // mf::LogError("CRTDataAnalysis") << "filled febMap with " << (this->fFebMap).size() << " entries" << std::endl;
    fin.close();
  }

  
  //-----------------------------------------------------------------------
  void CRTDataAnalysis::beginJob()
  {
     mf::LogError("CRTDataAnalysis") << " starting analysis job" << std::endl;

    // Access ART's TFileService, which will handle creating and writing
    // histograms and n-tuples for us. 
    art::ServiceHandle<art::TFileService> tfs;

    // Define our n-tuples
    fDAQNtuple        = tfs->make<TTree>("DAQTree",          "MyCRTDAQ");
    fHitNtuple        = tfs->make<TTree>("HitTree",          "MyCRTHit");
    fCRTPMTNtuple     = tfs->make<TTree>("CRTPMTTree",       "MyCRTPMTMatch");

    // Define the branches of our DetSim n-tuple 
    fDAQNtuple->Branch("event",                 &fDetEvent,          "event/I");
    fDAQNtuple->Branch("run",                   &fDetRun,            "run/I");
    fDAQNtuple->Branch("subrun",                &fDetSubRun,         "subrun/I");
    fDAQNtuple->Branch("nChan",                 &fNChan,             "nChan/I");
    fDAQNtuple->Branch("t0",                    &fT0,                "t0/l");
    fDAQNtuple->Branch("t1",                    &fT1,                "t1/l");
    fDAQNtuple->Branch("flags",                 &fFlags,             "flags/I");
    fDAQNtuple->Branch("nmaxch",                &fNMaxCh,            "nmaxch/I");
    fDAQNtuple->Branch("adc",                   fADC,                "adc[nmaxch]/I");
    fDAQNtuple->Branch("pe",                    fPE,                "pe[nmaxch]/F");
    fDAQNtuple->Branch("entry",                 &fEntry,             "entry/I");
    fDAQNtuple->Branch("mac5",                  &fMac5,              "mac5/I");
    fDAQNtuple->Branch("region",                &fFEBReg,            "region/I");
    fDAQNtuple->Branch("subSys",                &fDetSubSys,         "subSys/I");
    fDAQNtuple->Branch("gate_type", &m_gate_type, "gate_type/b");
    fDAQNtuple->Branch("gate_start_timestamp", &m_gate_start_timestamp, "gate_start_timestamp/l");

    // Define the branches of our SimHit n-tuple
    fHitNtuple->Branch("event",       &fHitEvent,    "event/I");
    fHitNtuple->Branch("nHit",        &fNHit,        "nHit/I");
    fHitNtuple->Branch("x",           &fXHit,        "x/F");
    fHitNtuple->Branch("y",           &fYHit,        "y/F");
    fHitNtuple->Branch("z",           &fZHit,        "z/F");
    fHitNtuple->Branch("xErr",        &fXErrHit,     "xErr/F");
    fHitNtuple->Branch("yErr",        &fYErrHit,     "yErr/F");
    fHitNtuple->Branch("zErr",        &fZErrHit,     "zErr/F");
    fHitNtuple->Branch("t0",          &fT0Hit,       "t0/l");
    fHitNtuple->Branch("t1",          &fT1Hit,       "t1/L");
    fHitNtuple->Branch("region",      &fHitReg,      "region/I");  
    //    fHitNtuple->Branch("tagger",      &ftagger,      "tagger/C");  
    fHitNtuple->Branch("subSys",      &fHitSubSys,   "subSys/I");
    fHitNtuple->Branch("modID",       &fHitMod,      "modID/I");
    fHitNtuple->Branch("stripID",     &fHitStrip,    "stripID/I");
    fHitNtuple->Branch("nFeb",        &fNHitFeb,     "nFeb/I");
    fHitNtuple->Branch("totPe",       &fHitTotPe,    "totPe/F");
    fHitNtuple->Branch("gate_type", &m_gate_type, "gate_type/b");
    fHitNtuple->Branch("gate_name", &m_gate_name);
    fHitNtuple->Branch("trigger_timestamp", &m_trigger_timestamp, "trigger_timestamp/l");
    fHitNtuple->Branch("gate_start_timestamp", &m_gate_start_timestamp, "gate_start_timestamp/l");
    fHitNtuple->Branch("trigger_gate_diff", &m_trigger_gate_diff, "trigger_gate_diff/l");
    fHitNtuple->Branch("gate_crt_diff",&m_gate_crt_diff, "gate_crt_diff/l");
    fHitNtuple->Branch("crt_global_trigger",&m_crt_global_trigger,"crt_global_trigger/l");
    fHitNtuple->Branch("crtGT_trig_diff",&m_crtGT_trig_diff,"crtGT_trig_diff/L");
    
    // Define the branches of our CRTPMTMatch ntuple
    fCRTPMTNtuple->Branch("event", &fMatchEvent, "event/I");
    fCRTPMTNtuple->Branch("run", &fMatchRun, "run/I");
    fCRTPMTNtuple->Branch("gate_type", &fGateType, "gate_type/b");
    fCRTPMTNtuple->Branch("fFlashID", &fFlashID);
    fCRTPMTNtuple->Branch("flashTime_us", &fFlashTime_us, "flashTime_us/D");
    fCRTPMTNtuple->Branch("flashGateTime_ns", &fFlashGateTime_ns, "flashGateTime_ns/D");
    fCRTPMTNtuple->Branch("firstOpHitPeakTime", &fFirstOpHitPeakTime);
    fCRTPMTNtuple->Branch("firstOpHitStartTime", &fFirstOpHitStartTime);
    fCRTPMTNtuple->Branch("flashInGate", &fFlashInGate, "flashInGate/O");
    fCRTPMTNtuple->Branch("flashInBeam", &fFlashInBeam, "flashInBeam/O");
    fCRTPMTNtuple->Branch("flashPE", &fFlashPE);
    fCRTPMTNtuple->Branch("fFlashPos_x", &fFlashPos_x, "flashPos_x/D");
    fCRTPMTNtuple->Branch("fFlashPos_y", &fFlashPos_y, "flashPos_y/D");
    fCRTPMTNtuple->Branch("fFlashPos_z", &fFlashPos_z, "flashPos_z/D");
    fCRTPMTNtuple->Branch("fFlashYWidth",&fFlashYWidth);
    fCRTPMTNtuple->Branch("fFlashZWidth",&fFlashZWidth);
    fCRTPMTNtuple->Branch("fFlashClassification", &fFlashClassification, "flashClassification/I");
    fCRTPMTNtuple->Branch("nMatchedCRTHits", &nMatchedCRTHits);
    fCRTPMTNtuple->Branch("CRTHitPos_x", &CRTHitPos_x);
    fCRTPMTNtuple->Branch("CRTHitPos_y", &CRTHitPos_y);
    fCRTPMTNtuple->Branch("CRTHitPos_z", &CRTHitPos_z);
    fCRTPMTNtuple->Branch("CRTPMTTimeDiff_ns", &fCRTPMTTimeDiff_ns);
    fCRTPMTNtuple->Branch("CRTTime_us", &fCRTTime_us);
    fCRTPMTNtuple->Branch("CRTSys", &fCRTSys);
    fCRTPMTNtuple->Branch("CRTRegion", &fCRTRegion);
    fCRTPMTNtuple->Branch("topCRTBefore", &fNtopCRTBefore);
    fCRTPMTNtuple->Branch("topCRTAfter", &fNtopCRTAfter);
    fCRTPMTNtuple->Branch("sideCRTBefore", &fNsideCRTBefore);
    fCRTPMTNtuple->Branch("sideCRTAfter", &fNsideCRTAfter);
}
  
  void CRTDataAnalysis::beginRun(const art::Run&)
  {
  }

  //-----------------------------------------------------------------------
  void CRTDataAnalysis::analyze(const art::Event& event) 
  {
    MF_LOG_DEBUG("CRTDataAnalysis") << "beginning analyis" << '\n';

    // Start by fetching some basic event information for our n-tuple.
    fEvent  = event.id().event(); 
    fRun    = event.run();
    fSubRun = event.subRun();
    FillFebMap();//febMap);

    //add trigger info
    m_gate_type = value(sbn::triggerSource::Unknown);
    if( !fTriggerLabel.empty() ) {

      art::Handle<sbn::ExtraTriggerInfo> trigger_handle;
      event.getByLabel( fTriggerLabel, trigger_handle );
      if( trigger_handle.isValid() ) {
	sbn::triggerSource bit = trigger_handle->sourceType;
        m_gate_type = value(bit);
        m_gate_name = bitName(bit);
        m_trigger_timestamp = trigger_handle->triggerTimestamp;
        m_gate_start_timestamp =  trigger_handle->beamGateTimestamp;
        m_trigger_gate_diff = trigger_handle->triggerTimestamp - trigger_handle->beamGateTimestamp;

      }
      else{
	mf::LogError("CRTDataAnalysis") << "No raw::Trigger associated to label: " << fTriggerLabel.label() << "\n" ;
      }
    }
    else {
       mf::LogError("CRTDataAnalysis")  << "Trigger Data product " << fTriggerLabel.label() << " not found!\n" ;
    }

    art::Handle<vector<icarus::crt::CRTData>> crtDAQHandle;
    vector<art::Ptr<icarus::crt::CRTData> > crtList;
    if ( event.getByLabel(fCRTDAQProducerLabel, crtDAQHandle))
      art::fill_ptr_vector(crtList, crtDAQHandle);
    
    vector<art::Ptr<icarus::crt::CRTData>> crtData;
    bool presel = false;

    for (size_t febdat_i=0; febdat_i<crtList.size(); febdat_i++) {
      
      uint8_t mac = crtList[febdat_i]->fMac5;
      int adid  = fCrtutils->MacToAuxDetID(mac,0);
      char type = fCrtutils->GetAuxDetType(adid);
      /*
      for(int chan=0; chan<32; chan++) {
	 mf::LogError("CRTDataAnalysis") << "\nfebP (mac5, channel, adc, type, adid) = (" << (int)crtList[febdat_i]->fMac5 << " , " << chan << " , "
		  << crtList[febdat_i]->fAdc[chan]  << " , " << type << " , " << adid << ")\n";
      }      
      */
      /// Looking for data within +/- 3ms within trigger time stamp
      /// Here t0 - trigger time -ve, only adding 1s makes the value +ve or -ve
      //    if (std::fabs(int64_t(crtList[febdat_i]->fTs0 - m_trigger_timestamp) + 1e9) > fCrtWindow) continue;
      if ( type == 'm'){
	for(int chan=0; chan<32; chan++) {
	  std::pair<double,double> const chg_cal = fChannelMap->getSideCRTCalibrationMap((int)crtList[febdat_i]->fMac5,chan);
	  float pe = (crtList[febdat_i]->fAdc[chan]-chg_cal.second)/chg_cal.first;
	  // In order to have Reset TS1 hits in CRTData from Side CRT, we have to explicitly include them
	  // The current threshold cut (6.5 PE) was applied to filter out noise, but this also filters out
	  // Reset events which are random trigger around the pedestal. These Reset hits are removed in 
	  // CRT Hit reconstruction. Top CRT has in internal triggering logic and threshold  that screens
	  // from the noise (hence presel = true for all the hits).
	  // Please revise this in the future if also T0 Reset hits need to be kept in CRTData. 
	  // To do so, include !0crtList[febdat_i]->IsReference_TS0()
	  if(pe<=fPEThresh && !crtList[febdat_i]->IsReference_TS1()) continue;
	  presel = true;
	}
      }else if ( type == 'c' ) {

	  presel = true;

      }else if ( type == 'd'){
	for(int chan=0; chan<64; chan++) {
	  float pe = (crtList[febdat_i]->fAdc[chan]-fQPed)/fQSlope;
	  if(pe<=fPEThresh) continue;
	  presel = true;
	}
      }
      if (presel) crtData.push_back(crtList[febdat_i]);
      presel = false;
    } // end of crtList
    
        
    mf::LogError("CRTDataAnalysis") << "about to loop over " << crtData.size() <<" crtData entries \n";
    for (size_t febdat_i=0; febdat_i<crtData.size(); febdat_i++) {
      
      
      fDetEvent       = fEvent;
      fDetRun         = fRun;
      fDetSubRun      = fSubRun;
      fMac5           = crtData[febdat_i]->fMac5;
      fEntry          = crtData[febdat_i]->fEntry;
      fFEBReg         = fCrtutils->AuxDetRegionNameToNum(fCrtutils->MacToRegion(fMac5));
      fNChan = 0;
      fDetSubSys = fCrtutils->MacToTypeCode(fMac5);
      fT0 = crtData[febdat_i]->fTs0;
      fT1 = crtData[febdat_i]->fTs1;
      fFlags = crtData[febdat_i]->fFlags;     
      int maxchan =0;
      if(fDetSubSys!=2) maxchan=32;
      else maxchan = 64;
      fNMaxCh = maxchan;
      for(int ch=0; ch<maxchan; ch++) {
	fADC[ch] = crtData[febdat_i]->fAdc[ch]; 
	std::pair<double,double> const chg_cal = fChannelMap->getSideCRTCalibrationMap((int)fMac5,ch);
	if (fDetSubSys == 0 || fDetSubSys == 1){
	  float pe = (fADC[ch]-chg_cal.second)/chg_cal.first;
	  if (pe < 0) continue;
	  fPE[ch] = pe;
	}else{
	  float pe = (fADC[ch]-fQPed)/fQSlope;
	  if (pe < 0) continue;
          fPE[ch] = pe;
	} 
	
      }
            
      fDAQNtuple->Fill();
      
    } //for CRT FEB events
    
  
    // Fill CRT Hit Tree
    art::Handle<std::vector<sbn::crt::CRTHit>> crtHitHandle;
    
    bool isCRTHit = event.getByLabel(fCRTHitProducerLabel, crtHitHandle);
    std::vector<int> ids;
    fNHit = 0;
    if (isCRTHit) {
      
       mf::LogError("CRTDataAnalysis") << "looping over reco hits..." << std::endl;
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
	   
	  fNHitFeb  = hit.feb_id.size();
	  fHitTotPe = hit.peshit;
	  int mactmp = hit.feb_id[0];
	  fHitReg  = fCrtutils->AuxDetRegionNameToNum(fCrtutils->MacToRegion(mactmp));
	  fHitSubSys =  fCrtutils->MacToTypeCode(mactmp);
	  
	  m_gate_crt_diff = m_gate_start_timestamp - hit.ts0_ns;
	  m_crt_global_trigger = hit.ts0_ns - hit.ts1_ns;
	  m_crtGT_trig_diff = m_crt_global_trigger - (m_trigger_timestamp%1'000'000'000);//'''						      
	  auto ittmp = hit.pesmap.find(mactmp);
	  if (ittmp==hit.pesmap.end()) {
	     mf::LogError("CRTDataAnalysis") << "hitreg: " << fHitReg << std::endl;
	     mf::LogError("CRTDataAnalysis") << "fHitSubSys: "<< fHitSubSys << std::endl;
	     mf::LogError("CRTDataAnalysis") << "mactmp = " << mactmp << std::endl;
	     mf::LogError("CRTDataAnalysis") << "could not find mac in pesmap!" << std::endl;
	    continue;
	  }	  
	  int chantmp = (*ittmp).second[0].first;
	  
	  fHitMod  = fCrtutils->MacToAuxDetID(mactmp, chantmp);
	  fHitStrip = fCrtutils->ChannelToAuxDetSensitiveID(mactmp, chantmp);
	  
	  fHitNtuple->Fill();
       }//for CRT Hits
    }//if CRT Hits
    
    else  mf::LogError("CRTDataAnalysis") << "CRTHit products not found! (expected if decoder step)" << std::endl;

       
    //Fill CRTPMT Match TTree
    art::Handle<vector<sbn::crt::CRTPMTMatching>> CRTPMTMatchingHandle;
    if ( event.getByLabel(fCRTPMTProducerLabel, CRTPMTMatchingHandle)){
      
      for (auto const& match: *CRTPMTMatchingHandle){
	int TopEn = 0, TopEx = 0, SideEn = 0, SideEx = 0;
	fMatchEvent = fEvent;
	fMatchRun = fRun;
	fGateType = m_gate_type;
	fFlashID = match.flashID;
	fFlashTime_us = match.flashTime;
	fFlashGateTime_ns = match.flashGateTime;
	fFirstOpHitPeakTime = match.firstOpHitPeakTime;
	fFirstOpHitStartTime = match.firstOpHitStartTime;
	fFlashInGate = match.flashInGate;
	fFlashInBeam = match.flashInBeam;
	fFlashPE = match.flashPE;
	fFlashPos_x = match.flashPosition.X();
	fFlashPos_y = match.flashPosition.Y();
	fFlashPos_z = match.flashPosition.Z();
	fFlashYWidth = match.flashYWidth;
	fFlashZWidth = match.flashZWidth;
	fFlashClassification = match.flashClassification;
	nMatchedCRTHits = match.matchedCRTHits.size();
	for(auto const& crthit: match.matchedCRTHits){
	  CRTHitPos_x.push_back(crthit.position.X());
	  CRTHitPos_y.push_back(crthit.position.Y());
	  CRTHitPos_z.push_back(crthit.position.Z());
	  fCRTPMTTimeDiff_ns.push_back(1e3*crthit.PMTTimeDiff);
	  fCRTTime_us.push_back(crthit.time);
	  fCRTSys.push_back(crthit.sys);
	  fCRTRegion.push_back(crthit.region);
	  int fMatchType = static_cast<int>(fFlashClassification);
	  if(fMatchType == 1 || fMatchType == 3 || fMatchType == 6 || fMatchType == 7 || fMatchType == 11) TopEn++;
	  if(fMatchType == 4 || fMatchType == 13) TopEx++;
	  if(fMatchType == 2 || fMatchType == 12) SideEn++;
	  if(fMatchType == 3 || fMatchType == 5 || fMatchType == 7 || fMatchType == 14) SideEx++;
	}
	fNtopCRTBefore = TopEn;
	fNtopCRTAfter = TopEx;
	fNsideCRTBefore = SideEn;
	fNsideCRTAfter = SideEx;
	fCRTPMTNtuple->Fill();
	ClearVecs();
      } // for match in handle 
    } // if valid label 
    else{
      mf::LogError("CRTDataAnalysis") << "not Valid CRTPMTProducer label!\n";
    }
    

  } // CRTDataAnalysis::analyze()
  
    void CRTDataAnalysis::ClearVecs(){
      CRTHitPos_x.clear();
      CRTHitPos_y.clear();
      CRTHitPos_z.clear();
      fCRTPMTTimeDiff_ns.clear();
      fCRTTime_us.clear();
      fCRTSys.clear();
      fCRTRegion.clear();
    }
  
  DEFINE_ART_MODULE(CRTDataAnalysis)

} // namespace crt
} // namespace icarus

