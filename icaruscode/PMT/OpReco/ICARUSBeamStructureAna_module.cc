////////////////////////////////////////////////////////////////////////
// Class:       ICARUSBeamStructureAna
// Plugin Type: analyzer (Unknown Unknown)
// File:        ICARUSBeamStructureAna_module.cc
//
// Generated at Mon Jul 22 13:50:03 2024 by Matteo Vicenzi using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

// framework libraries
#include "canvas/Utilities/InputTag.h"
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Utilities/Exception.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "canvas/Persistency/Common/Assns.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "art_root_io/TFileService.h"

// LArSoft libraries
#include "icaruscode/Decode/DataProducts/TriggerConfiguration.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "icaruscode/Decode/ChannelMapping/IICARUSChannelMap.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataalg/DetectorInfo/DetectorTimings.h"
#include "lardataalg/DetectorInfo/DetectorTimingTypes.h" // electronics_time
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/Simulation/BeamGateInfo.h"
#include "lardataobj/RawData/TriggerData.h"
#include "sbnobj/Common/Trigger/ExtraTriggerInfo.h"
#include "icaruscode/CRT/CRTUtils/CRTPMTMatchingUtils.h"
#include "sbnobj/Common/CRT/CRTPMTMatching.hh"
#include "icaruscode/IcarusObj/PMTBeamSignal.h"

// ROOT libraries
#include "TTree.h"
#include "TFile.h"

// C/C++ standard libraries
#include <string>
#include <vector>
#include <map>

// -----------------------------------------------------------------------------
namespace opana {   class ICARUSBeamStructureAna; }

class opana::ICARUSBeamStructureAna : public art::EDAnalyzer {
  public:
    struct Config {
      fhicl::Sequence<art::InputTag> FlashLabels {
        fhicl::Name("FlashLabels"),
        fhicl::Comment("Tags for the recob::OpFlash data products")
      };
      fhicl::Atom<art::InputTag> TriggerLabel {
        fhicl::Name("TriggerLabel"),
        fhicl::Comment("Tag for trigger info")
      };
      fhicl::Atom<art::InputTag> RWMLabel {
        fhicl::Name("RWMLabel"),
        fhicl::Comment("Tag for RWM info")
      };
      fhicl::Atom<art::InputTag> TriggerConfigLabel {
        fhicl::Name("TriggerConfigLabel"),
        fhicl::Comment("Trigger configuration label")
      };
      fhicl::Atom<art::InputTag> CRTPMTMatchingLabel {
        fhicl::Name("CRTPMTMatchingLabel"),
        fhicl::Comment("CRTPMT matching label")
      };			
    }; // struct Config 

    using Parameters = art::EDAnalyzer::Table<Config>;

    /// constructor
    explicit ICARUSBeamStructureAna(Parameters const& config);
    
    /// Return RWM-relative time from a trigger-relative time
    double getRWMRelativeTime(int channel, double t);
 
    /// Return side/wall from channel id
    int getSideByChannel(const int channel);
    
    /// Return the RWM-relative flash interaction time
    double getFlashBunchTime(std::vector<double> hit_rise_time); 
    
    /// Clear all data structures
    void clear();

    void analyze(art::Event const& e) override;
    void beginJob();
    void beginRun(const art::Run& run) override;

  private:

    art::ServiceHandle<art::TFileService> tfs;
    icarus::TriggerConfiguration fTriggerConfiguration;

    std::vector<art::InputTag> fFlashLabels;
    art::InputTag fTriggerLabel;
    art::InputTag fRWMLabel;
    art::InputTag fTriggerConfigurationLabel;
    art::InputTag fCRTPMTMatchingLabel;

    /// data members
    std::vector<TTree*> fOpFlashTrees;

    int m_run;
    int m_event;
    int m_timestamp; 

    // trigger info 
    unsigned int m_gate_type;
    int m_trigger_type=-1;
    std::string m_gate_name;
    uint64_t m_trigger_timestamp;
    uint64_t m_beam_gate_timestamp;
    double m_beam_us;
    double m_trigger_us;
    double m_beam_gate_width;
		
    // flash info
    int m_cryo; 
    int m_flash_id;
    double m_flash_time;
    double m_flash_time_rwm;
    double m_flash_z;
    double m_flash_y;
    double m_flash_pe;
    int m_flash_nhits;
    std::vector<int> m_channel_id;
    std::vector<double> m_hit_start_time;
    std::vector<double> m_hit_peak_time;
    std::vector<double> m_hit_rise_time;
    std::vector<double> m_hit_start_time_rwm;
    std::vector<double> m_hit_peak_time_rwm;
    std::vector<double> m_hit_rise_time_rwm;
    std::vector<double> m_hit_pe;

    //crt-pmt match
    int m_flash_classification;
    int m_flash_ncrthits;
    std::vector<double> m_crthit_x;
    std::vector<double> m_crthit_y;
    std::vector<double> m_crthit_z;
    std::vector<double> m_crttime_us;
    std::vector<double> m_crtpmttimediff_ns;
    std::vector<int> m_crtsys;
    std::vector<int> m_crtregion;

    // RWM times
    std::vector<icarus::timing::PMTBeamSignal> fRWMTimes;
};

// --------------------------------------------------------------------------
opana::ICARUSBeamStructureAna::ICARUSBeamStructureAna(Parameters const& config)
  : art::EDAnalyzer(config)
  , fFlashLabels( config().FlashLabels() )
  , fTriggerLabel( config().TriggerLabel() )
  , fRWMLabel ( config().RWMLabel() )
  , fTriggerConfigurationLabel( config().TriggerConfigLabel() )
  , fCRTPMTMatchingLabel( config().CRTPMTMatchingLabel() )
{}

// ---------------------------------------------------------------------------
void opana::ICARUSBeamStructureAna::beginJob() {

  if ( !fFlashLabels.empty() ) {
    
    for( auto const & label : fFlashLabels ) {
      
      // TTree for the flash in a given cryostat
      std::string name = "beamtiming_"+label.label();
      std::string info = "Beam timing from label "+label.label();

      TTree* ttree = tfs->make<TTree>(name.c_str(), info.c_str() );
      ttree->Branch("run",&m_run);
      ttree->Branch("event",&m_event);
      ttree->Branch("timestamp",&m_timestamp);
  
      ttree->Branch("gate_type", &m_gate_type);
      ttree->Branch("gate_name", &m_gate_name);
      ttree->Branch("beam_gate_us", &m_beam_us);
      ttree->Branch("trigger_us", &m_trigger_us);
      ttree->Branch("beam_gate_width",&m_beam_gate_width);
      ttree->Branch("trigger_type", &m_trigger_type, "trigger_type/I");
      ttree->Branch("trigger_timestamp", &m_trigger_timestamp, "trigger_timestamp/l");
      ttree->Branch("beam_gate_timestamp", &m_beam_gate_timestamp, "beam_gate_timestamp/l");
	
      ttree->Branch("cryo",&m_cryo);
      ttree->Branch("flash_id",&m_flash_id);
      ttree->Branch("flash_time",&m_flash_time);
      ttree->Branch("flash_time_rwm",&m_flash_time_rwm);
      ttree->Branch("flash_pe",&m_flash_pe);
      ttree->Branch("flash_z",&m_flash_z);
      ttree->Branch("flash_y",&m_flash_y);
      ttree->Branch("flash_nhits",&m_flash_nhits);
      ttree->Branch("channels",&m_channel_id);
      ttree->Branch("hit_start_time",&m_hit_start_time);
      ttree->Branch("hit_peak_time",&m_hit_peak_time);
      ttree->Branch("hit_rise_time",&m_hit_rise_time);
      ttree->Branch("hit_start_time_rwm",&m_hit_start_time_rwm);
      ttree->Branch("hit_peak_time_rwm",&m_hit_peak_time_rwm);
      ttree->Branch("hit_rise_time_rwm",&m_hit_rise_time);
      ttree->Branch("hit_pe",&m_hit_pe);

      ttree->Branch("flash_classification",&m_flash_classification);
      ttree->Branch("flash_ncrthits",&m_flash_ncrthits);
      ttree->Branch("crthit_x",&m_crthit_x);
      ttree->Branch("crthit_y",&m_crthit_y);
      ttree->Branch("crthit_z",&m_crthit_z);
      ttree->Branch("crthit_sys",&m_crtsys);
      ttree->Branch("crthit_region",&m_crtregion);	
      ttree->Branch("crttime_us",&m_crttime_us);
      ttree->Branch("crtpmttimediff_ns",&m_crtpmttimediff_ns);

      fOpFlashTrees.push_back( ttree );
    }
  } else {
    mf::LogError("ICARUSBeamStructureAna")
      << "No flash labels selected!!"; 
  }
}

// ------------------------------------------------------------------------------ 

void opana::ICARUSBeamStructureAna::beginRun(const art::Run& r)
{
  fTriggerConfiguration = r.getProduct<icarus::TriggerConfiguration>(fTriggerConfigurationLabel);
}

// -------------------------------------------------------------------------------

void opana::ICARUSBeamStructureAna::analyze(art::Event const& e)
{
  // ----
  // Event metadata information
  m_run = e.id().run();
  m_event = e.id().event();
  m_timestamp = e.time().timeHigh(); // precision to the second    
  
  // ----
  // Trigger metadata information
  detinfo::DetectorTimings const detTimings = detinfo::makeDetectorTimings(
                                              art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e));
  detinfo::timescales::electronics_time triggerTime = detTimings.TriggerTime();
  detinfo::timescales::electronics_time beamGateTime =  detTimings.BeamGateTime();

  if( !fTriggerLabel.empty() ) { 

    auto const& extraInfo = e.getProduct<sbn::ExtraTriggerInfo>(fTriggerLabel);
    
    if( extraInfo.isValid() ) {

      sbn::triggerSource bit = extraInfo.sourceType;
      m_gate_type = (unsigned int)bit; //1 BNB 2 NumI 3 offbeamBNB 4 offbeamNuMi
      m_gate_name = bitName(bit);
      m_trigger_type = value( extraInfo.triggerType ); //1 majority, 2 minbias

      // absolute timestamp
      m_trigger_timestamp = extraInfo.triggerTimestamp;
      m_beam_gate_timestamp =  extraInfo.beamGateTimestamp;
     
      // time in electronics time
      m_trigger_us = triggerTime.value();
      m_beam_us = beamGateTime.value();
      m_beam_gate_width = fTriggerConfiguration.getGateWidth(m_gate_type);

    } else { 
      mf::LogError("ICARUSBeamStructureAna") << "No raw::Trigger associated to label: " << fTriggerLabel.label() << "!" ; 
    }
  } else {
    mf::LogError("ICARUSBeamStructureAna") << "No trigger labels selected!!" ; 
  }
  
  // ----
  // RWM times

  fRWMTimes = e.getProduct<std::vector<icarus::timing::PMTBeamSignal>>(fRWMLabel);
  if ( fRWMTimes.empty() )
    mf::LogWarning("ICARUSBeamStructureAna") << "Data product std::vector<icarus::timing::PMTBeamSignal for " << fRWMLabel.label()
                                        << " is empty in " << m_gate_name << " event!";
  // ----
  // FLASH/CRT timing information

  if ( !fFlashLabels.empty() ) {

    for ( size_t iFlashLabel=0; iFlashLabel<fFlashLabels.size(); iFlashLabel++  ) {
      
      auto const label = fFlashLabels[iFlashLabel];
      m_cryo = iFlashLabel;

      auto const& flash_handle = e.getValidHandle<std::vector<recob::OpFlash>>(label);
      auto const& flashes = *flash_handle;

      // we want our flashes to be valid and not empty
      if( !flash_handle.isValid() ) {
        mf::LogError("ICARUSBeamStructureAna") << "Not found a recob::OpFlash with label '" << label.encode() << "'"; 
      } else if ( flashes.empty() ) {
        mf::LogWarning("ICARUSBeamStructureAna") << "No recob::OpFlash in collection with label '" << label.encode() << "'"; 
      } else {
        
        art::FindManyP<recob::OpHit> ophitsPtr( flash_handle, e, label );
        art::FindOneP<sbn::crt::CRTPMTMatching> matchPtr( flash_handle, e, fCRTPMTMatchingLabel );

	// loop all flashes
	size_t idx = 0;
        for ( auto const& flash : flashes ) {
          
 	  // Filling flash info...
          m_flash_id = idx;
          m_flash_time = flash.Time();
          m_flash_pe = flash.TotalPE();
          m_flash_z = flash.ZCenter();
          m_flash_y = flash.YCenter();
          
          // ----
          // CRT match info   

          auto const & match = matchPtr.at(idx);
          
	  // if there is no match, there is no product
	  // fill null parameters
	  if( !match ) {
	    m_flash_classification = 0;
	    m_flash_ncrthits = 0;		
	  } else {
	    m_flash_classification = static_cast<int>(match->flashClassification);
	    m_flash_ncrthits = match->matchedCRTHits.size();
					
	    for( auto const& crthit: match->matchedCRTHits ){
	      m_crthit_x.push_back(crthit.position.X());
              m_crthit_y.push_back(crthit.position.Y());
              m_crthit_z.push_back(crthit.position.Z());
              m_crttime_us.push_back(crthit.time);
              m_crtpmttimediff_ns.push_back(1e3*crthit.PMTTimeDiff);
              m_crtsys.push_back(crthit.sys);
              m_crtregion.push_back(crthit.region);
            }
          }
	
          // ----
          // OPHITS info

	  auto const & ophits = ophitsPtr.at(idx);
     
          m_channel_id.resize(360);     
          m_hit_start_time.resize(360);
          m_hit_peak_time.resize(360);
          m_hit_rise_time.resize(360);
          m_hit_start_time_rwm.resize(360);
          m_hit_peak_time_rwm.resize(360);
          m_hit_rise_time_rwm.resize(360);
          m_hit_pe.resize(360);
				
          // loop all hits in the flash: save only the first one
          for ( auto const hit : ophits ){
            
            const int ch = hit->OpChannel();
            double ts = hit->StartTime();
	    double tp = hit->PeakTime();
	    double tr = hit->RiseTime();
	    double pe = hit->PE();
	
            // select the first ophit (by time) in each channel
            if( ( m_hit_start_time[ch] == 0 ) || ( m_hit_start_time[ch] > ts )) {
              m_channel_id[ch] = ch;
              m_hit_start_time[ch] = ts;
              m_hit_peak_time[ch] = tp;
              m_hit_rise_time[ch] = ts + tr;
              m_hit_start_time_rwm[ch] = getRWMRelativeTime(ch, ts);
              m_hit_peak_time_rwm[ch] = getRWMRelativeTime(ch, tp);
              m_hit_rise_time_rwm[ch] = getRWMRelativeTime(ch, ts+tr);
              m_hit_pe[ch] = pe;
            }
          }

          // get the flash interaction time w.r.t. RWM
          // this is currently the mean between the first ophits on opposite walls
          m_flash_time_rwm = getFlashBunchTime(m_hit_rise_time_rwm);

          fOpFlashTrees[iFlashLabel]->Fill();
          clear();	  
        }
      }
    }
  }
}

// -----------------------------------------------------------------------------

int opana::ICARUSBeamStructureAna::getSideByChannel( const int channel ) {

  /* 
  Channels are numbered from east to west, from North (cryo side) to South (beam side)
  We look in the opposide direction wrt to the beam direction South->North: 
  - Left is the east wall of each cryostat;
  - Right is the west side of each cryostat;
  - [ 0:89 ] and [180:269] are on the left, 
    the return value of the function is 0;
  - [ 90-179 ] and [ 270:359 ] are on the right,
    the return value of the function is 1;
  */

  int side = channel / 90; // always round down
  return side % 2;
}

// -----------------------------------------------------------------------------

double opana::ICARUSBeamStructureAna::getRWMRelativeTime(int channel, double t) {

  if( fRWMTimes.empty() ) return 0;

  auto rwm = fRWMTimes.at(channel);
  if ( !rwm.isValid() ) return 0;

  double rwm_trigger = rwm.startTime; //rwm time w.r.t. trigger time [us]
  return (t - rwm_trigger);

}

// -----------------------------------------------------------------------------

double opana::ICARUSBeamStructureAna::getFlashBunchTime(std::vector<double> hit_rise_time_rwm) {

  float tfirst_left = std::numeric_limits<float>::max();
  float tfirst_right = std::numeric_limits<float>::max();
    
  // if no RWM info available, all pmt_start_time_rwm is zero
  // return zero as well for the flash
  if (fRWMTimes.empty()) return 0;
  
  for(size_t i=0; i<hit_rise_time_rwm.size(); i++ ){
  
    int cryo = getSideByChannel(i);
    float t = hit_rise_time_rwm[i]; // rise time w.r.t. rwm

    // if RWM is missing for some PMT channels, 
    // it might not be possible to use the first hits (might not have RMW time)
    // so return zero as in other bad cases
    if( !fRWMTimes[i].isValid() ) return 0; 
 
    if( (cryo == 0) && (t < tfirst_left) )
      tfirst_left = t;   
    if( (cryo == 1) && (t < tfirst_right) )
      tfirst_right = t;   
  } 

  return (tfirst_left + tfirst_right)/2.;

} 

// -----------------------------------------------------------------------------

void opana::ICARUSBeamStructureAna::clear() {
  
          m_channel_id.clear();
          m_hit_start_time.clear();
	  m_hit_peak_time.clear();
	  m_hit_rise_time.clear();
	  m_hit_start_time_rwm.clear();
	  m_hit_peak_time_rwm.clear();
	  m_hit_rise_time_rwm.clear();
          m_hit_pe.clear();
	      
          m_crthit_x.clear();
          m_crthit_y.clear();
          m_crthit_z.clear();
          m_crttime_us.clear();
          m_crtpmttimediff_ns.clear();
          m_crtsys.clear();
          m_crtregion.clear();
}

// -----------------------------------------------------------------------------

DEFINE_ART_MODULE(opana::ICARUSBeamStructureAna)
