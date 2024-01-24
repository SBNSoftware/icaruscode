////////////////////////////////////////////////////////////////////////
// Class:       HitTiming
// Plugin Type: analyzer (Unknown Unknown)
// File:        HitTiming_module.cc
//
// Generated at Tue Nov 21 10:35:03 2023 by Matteo Vicenzi using cetskelgen
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
#include "canvas/Persistency/Common/Assns.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "art_root_io/TFileService.h"

// LArSoft libraries
#include "icaruscode/Decode/DataProducts/TriggerConfiguration.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "icaruscode/IcarusObj/PMTWaveformTimeCorrection.h"
#include "icaruscode/Decode/ChannelMapping/IICARUSChannelMap.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataalg/DetectorInfo/DetectorTimings.h"
#include "lardataalg/DetectorInfo/DetectorTimingTypes.h" // electronics_time
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RawData/TriggerData.h"
#include "sbnobj/Common/Trigger/ExtraTriggerInfo.h"

// ROOT libraries
#include "TTree.h"
#include "TFile.h"

// C/C++ standard libraries
#include <string>
#include <vector>
#include <map>

// -----------------------------------------------------------------------------
namespace icarus {   class HitTiming; }

class icarus::HitTiming : public art::EDAnalyzer {
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
			fhicl::Atom<art::InputTag> CorrectionLabel {
				fhicl::Name("CorrectionLabel"),
				fhicl::Comment("Instance for correction info")
			};
			fhicl::Atom<art::InputTag> TriggerConfigLabel {
				fhicl::Name("TriggerConfigLabel"),
				fhicl::Comment("Trigger configuration label")
			};
			

		}; // struct Config 

		using Parameters = art::EDAnalyzer::Table<Config>;

		/// constructor
		explicit HitTiming(Parameters const& config);

		/// process the event
  		template<typename T> T Median(std::vector<T> data) const;
        	template<typename T> static size_t getMaxBin(std::vector<T> const& vv, size_t startElement, size_t endElement);
		template<typename T> static size_t getMinBin(std::vector<T> const& vv, size_t startElement, size_t endElement);
		template<typename T> static size_t getStartSample( std::vector<T> const& vv );
		
		void analyze(art::Event const& e) override;
		void beginJob();
		void beginRun(const art::Run& run) override;
		void endJob();

	private:

		art::ServiceHandle<art::TFileService> tfs;
  		icarus::TriggerConfiguration fTriggerConfiguration;

		std::vector<art::InputTag> fFlashLabels;
  		art::InputTag fTriggerLabel;
		art::InputTag fCorrectionLabel;
		art::InputTag fTriggerConfigurationLabel;

		/// data members
		std::vector<TTree*> fOpFlashTrees;
		TTree* fRWMTree;

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
			
		// special signal
		int m_n_ew;
		int m_n_rwm;
		int m_n_trig;
		std::vector<int> m_ew_channel;
		std::vector<double> m_ew_start;
		std::vector<double> m_ew_time;
		std::vector<int> m_rwm_channel;
		std::vector<double> m_rwm_start;
		std::vector<double> m_rwm_time;
		std::vector<int> m_trig_channel;
		std::vector<double> m_trig_corr;
	
		// flash info
		int m_cryo; 
		int m_flash_id;
		double m_flash_time;
		double m_flash_z;
		double m_flash_y;
		double m_flash_pe;
		int m_flash_nhits;
		std::vector<int> m_channel_id;
		std::vector<double> m_start_time;
		std::vector<double> m_peak_time;
		std::vector<double> m_rise_time;
		std::vector<double> m_hit_pe;
};

// --------------------------------------------------------------------------
icarus::HitTiming::HitTiming(Parameters const& config)
	: art::EDAnalyzer(config)
	  , fFlashLabels( config().FlashLabels() )
    	  , fTriggerLabel( config().TriggerLabel() )
	  , fCorrectionLabel( config().CorrectionLabel() ) 
	  , fTriggerConfigurationLabel( config().TriggerConfigLabel() )
{ 

}

// ---------------------------------------------------------------------------
void icarus::HitTiming::beginJob() {

	if ( !fFlashLabels.empty() ) {

		for( auto const & label : fFlashLabels ) {

			// TTree for the flash in a given cryostat
			std::string name = "hittiming_"+label.label();
			std::string info = "Time info from label "+label.label();

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
			ttree->Branch("flash_pe",&m_flash_pe);
			ttree->Branch("flash_z",&m_flash_z);
			ttree->Branch("flash_y",&m_flash_y);
			ttree->Branch("flash_nhits",&m_flash_nhits);
			ttree->Branch("channels",&m_channel_id);
			ttree->Branch("start_time",&m_start_time);
			ttree->Branch("peak_time",&m_peak_time);
			ttree->Branch("rise_time",&m_rise_time);
			ttree->Branch("hit_pe",&m_hit_pe);

			fOpFlashTrees.push_back( ttree );
		}
	}

	fRWMTree = tfs->make<TTree>("rwmtree","RWM and EW info");
	fRWMTree->Branch("run",&m_run);
	fRWMTree->Branch("event",&m_event);
	fRWMTree->Branch("timestamp",&m_timestamp);
	fRWMTree->Branch("n_ew",&m_n_ew);
	fRWMTree->Branch("ew_channel",&m_ew_channel);
	fRWMTree->Branch("ew_start", &m_ew_start);
	fRWMTree->Branch("ew_time", &m_ew_time);
	fRWMTree->Branch("n_rwm",&m_n_rwm);
	fRWMTree->Branch("rwm_channel",&m_rwm_channel);
	fRWMTree->Branch("rwm_start", &m_rwm_start);
	fRWMTree->Branch("rwm_time", &m_rwm_time);
	fRWMTree->Branch("n_trig",&m_n_trig);
	fRWMTree->Branch("trig_channel",&m_trig_channel);		
	fRWMTree->Branch("trig_corr", &m_trig_corr);

}

// ---------------------------------------------------------------------------

template<typename T> T icarus::HitTiming::Median( std::vector<T> data ) const {

    std::nth_element( data.begin(), data.begin() + data.size()/2, data.end() );
    return data[ data.size()/2 ];

}

// -----------------------------------------------------------------------------

template<typename T>
  size_t icarus::HitTiming::getMinBin( 
        std::vector<T> const& vv, size_t startElement, size_t endElement ){

    auto minel = 
        std::min_element( vv.begin()+startElement, vv.begin()+endElement );
    size_t minsample = std::distance( vv.begin()+startElement, minel );

    return minsample;
}

// -----------------------------------------------------------------------------

template<typename T>
  size_t icarus::HitTiming::getMaxBin( 
            std::vector<T> const& vv, size_t startElement, size_t endElement){

    auto maxel = 
        std::max_element( vv.begin()+startElement, vv.begin()+endElement );
    
    size_t maxsample = std::distance( vv.begin()+startElement, maxel );

    return maxsample;
} 

// -----------------------------------------------------------------------------

template<typename T>
  size_t icarus::HitTiming::getStartSample( std::vector<T> const& vv ){
    
    // NOTE: when changing this algorithm, also update the documentation
    // in the section "Signal timing extraction" of the class documentation

    // We are thinking in inverted polarity
    size_t minbin = getMinBin( vv, 0, vv.size() );

    //Search only a cropped region of the waveform backward from the min
    size_t maxbin =  minbin-20; //getMaxBin( wave, minbin-20, minbin );

    // Now we crawl betweem maxbin and minbin and we stop when:
      //( maxbin value - minbin value )*0.05 > (maxbin value - bin value )
    size_t startbin = maxbin;
    auto delta = vv[maxbin]-vv[minbin];
    for( size_t bin=maxbin; bin<minbin; bin++ ){
      auto val = vv[maxbin]-vv[bin];
      if( val >= 0.2*delta ){
        startbin = bin - 1;
        break;
      }
    }

    if( startbin < maxbin ){
      startbin=maxbin;
    }

    return startbin;
}


// ------------------------------------------------------------------------------ 

void icarus::HitTiming::beginRun(const art::Run& r)
{
  fTriggerConfiguration = r.getProduct<icarus::TriggerConfiguration>(fTriggerConfigurationLabel);
}

// -------------------------------------------------------------------------------

void icarus::HitTiming::analyze(art::Event const& e)
{
	m_run = e.id().run();
	m_event = e.id().event();
	m_timestamp = e.time().timeHigh(); // precision to the second    
  
	detinfo::DetectorTimings const detTimings = detinfo::makeDetectorTimings(art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e));
  	detinfo::timescales::electronics_time triggerTime = detTimings.TriggerTime();
  	detinfo::timescales::electronics_time beamGateTime =  detTimings.BeamGateTime();

  	// We work out the trigger information here 
  	if( !fTriggerLabel.empty() ) { 

      		// Trigger information
      		art::Handle<sbn::ExtraTriggerInfo> trigger_handle;
      		e.getByLabel( fTriggerLabel, trigger_handle );

      		if( trigger_handle.isValid() ) {

        		sbn::triggerSource bit = trigger_handle->sourceType;

        		m_gate_name = bitName(bit); //1 BNB 2 NumI 3 offbeamBNB 4 offbeamNuMi
			m_gate_type = (unsigned int)bit;
        		m_trigger_type = value( trigger_handle->triggerType ); //1 majority, 2 minbias
        
			// absolute timestamp
			m_trigger_timestamp = trigger_handle->triggerTimestamp;
        		m_beam_gate_timestamp =  trigger_handle->beamGateTimestamp;
     
			// time in electronics time
			m_trigger_us = triggerTime.value();
			m_beam_us = beamGateTime.value();
			m_beam_gate_width = fTriggerConfiguration.getGateWidth(m_gate_type);

      		}
}
	if ( !fFlashLabels.empty() ) {

		for ( size_t iFlashLabel=0; iFlashLabel<fFlashLabels.size(); iFlashLabel++  ) {

			auto const label = fFlashLabels[iFlashLabel];
			m_cryo = iFlashLabel;

			art::Handle<std::vector<recob::OpFlash>> flash_handle;
			e.getByLabel( label, flash_handle );

			// We want our flashes to be valid and not empty
			if( !flash_handle.isValid() ) {
				mf::LogError("HitTiming")
					<< "Not found a recob::OpFlash with label '" << label.encode() << "'"; 
			} else if ( flash_handle->empty() ) {
				mf::LogWarning("HitTiming")
					<< "No recob::OpFlash in collection with label '" << label.encode() << "'"; 
			}
			else {

				art::FindManyP<recob::OpHit> ophitsPtr( flash_handle, e, label );

				// loop all flashes
				for ( size_t idx=0; idx<flash_handle->size(); idx++ ) {

					m_channel_id.clear();
					m_start_time.clear();
					m_peak_time.clear();
					m_rise_time.clear();
					m_hit_pe.clear();

					m_flash_id = idx;
					auto const & flash = (*flash_handle)[idx];

					m_flash_time = flash.Time();
					m_flash_pe = flash.TotalPE();
					m_flash_z = flash.ZCenter();
					m_flash_y = flash.YCenter();

					auto const & ophits = ophitsPtr.at(idx);

					std::map<int,double> hitmap;
					std::map<int,double> peakmap;
					std::map<int,double> risemap;
					std::map<int,double> pemap;

					// loop all hits in the flash: save only the first one
					for ( auto const hit : ophits ){

						const int ch = hit->OpChannel();
						double ts = hit->StartTime();
						double tp = hit->PeakTime();
						double tr = hit->RiseTime();
						double pe = hit->PE();
	
						if ( hitmap.find(ch) != hitmap.end() ){
							if ( ts < hitmap[ch] ){
								hitmap[ch] = ts;
								peakmap[ch] = tp;
								risemap[ch] = tr;
								pemap[ch] = pe;
							}	
						}else{
							 hitmap.insert(std::make_pair(ch,ts));
							 peakmap.insert(std::make_pair(ch,tp));
							 risemap.insert(std::make_pair(ch,tr));
							 pemap.insert(std::make_pair(ch,pe));
						}
					}

					m_flash_nhits = hitmap.size();

					for(auto it = hitmap.begin(); it != hitmap.end(); it++){
						m_channel_id.push_back(it->first);
						m_start_time.push_back(it->second);
						m_peak_time.push_back(peakmap.at(it->first));
						m_rise_time.push_back(risemap.at(it->first));
						m_hit_pe.push_back(pemap.at(it->first));
						//std::cout<< "ch " << it->first << " t " << it->second << std::endl;
					}		  

					fOpFlashTrees[iFlashLabel]->Fill();
					hitmap.clear();
					peakmap.clear();
					risemap.clear();
					pemap.clear();
				}
			}
		}
	}

	m_ew_channel.clear();
	m_ew_start.clear();
	m_ew_time.clear();
	m_rwm_channel.clear();
	m_rwm_start.clear();
	m_rwm_time.clear();
	m_trig_channel.clear();
	m_trig_corr.clear();

	// get trig corrections
	auto const& wfCorrections = e.getProduct<std::vector<icarus::timing::PMTWaveformTimeCorrection>>(fCorrectionLabel);
	m_n_trig = wfCorrections.size();
			
	if( m_n_trig < 1 ) mf::LogError("HitTimin") << "Not found PMTWaveformTimeCorrections with label '" << fCorrectionLabel.label() << "'"; 
	
	for (unsigned int k=0; k<wfCorrections.size(); k++){
		auto wfc = wfCorrections.at(k);
		m_trig_channel.push_back(wfc.channelID);
		m_trig_corr.push_back(wfc.startTime);
	}

	// get trig corrections
	auto const& ewWaveforms = e.getProduct<std::vector<raw::OpDetWaveform>>("daqPMT:EW");
	auto const& rwmWaveforms = e.getProduct<std::vector<raw::OpDetWaveform>>("daqPMT:RWM");
	m_n_ew = ewWaveforms.size();
	m_n_rwm = rwmWaveforms.size();
	if( m_n_ew< 1 ) mf::LogError("HitTimin") << "Not found raw::OpDetWaveform with label 'daqPMT:EW'"; 
	if( m_n_rwm< 1 ) mf::LogError("HitTimin") << "Not found raw::OpDetWaveform with label 'daqPMT:RWM'"; 

        for( auto const & wave : ewWaveforms ){
            m_ew_channel.push_back( wave.ChannelNumber() );
	    detinfo::timescales::electronics_time tstart = util::quantities::points::microsecond{wave.TimeStamp()} ;
            m_ew_start.push_back( tstart.value() );
            
   	    double sample_rise = getStartSample( wave.Waveform() ); 
	    m_ew_time.push_back( tstart.value() + 0.002*sample_rise );
        }

        for( auto const & wave : rwmWaveforms ){
            
	    m_rwm_channel.push_back( wave.ChannelNumber() );
	    detinfo::timescales::electronics_time tstart = util::quantities::points::microsecond{wave.TimeStamp()} ;
            m_rwm_start.push_back( tstart.value() );
            
   	    double sample_rise = getStartSample( wave.Waveform() ); 
	    m_rwm_time.push_back( tstart.value() + 0.002*sample_rise );
        }

	fRWMTree->Fill();
}

// ----------------------------------------------------------------------------
void icarus::HitTiming::endJob()
{
}

// -----------------------------------------------------------------------------

DEFINE_ART_MODULE(icarus::HitTiming)
