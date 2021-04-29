////////////////////////////////////////////////////////////////////////
// 
// Class:       PMTMonitor
// Plugin Type: analyzer (art v3_05_00)
// File:        PMTMonitor_module.cc
//
// Generated at Mon Sep 21 15:21:37 2020 by Andrea Scarpelli
// 
// Base for the PMT Online Monitor. Reuse some code from the PMTDecoder_tool.cc
// and FullOpHitFinder_module.cc 
// 
// Use it as base for the new PMT Online Monitor
//
// mailto:ascarpell@bnl.gov
//
////////////////////////////////////////////////////////////////////////

#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/FileBlock.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

#include "artdaq-core/Data/Fragment.hh" // Fragment
#include "sbndaq-artdaq-core/Overlays/Common/CAENV1730Fragment.hh" // Fragment
#include "icaruscode/Decode/ChannelMapping/IICARUSChannelMap.h" // Channel map

// Base optical reconstruction 
#include "lardataobj/RawData/OpDetWaveform.h"
#include "larana/OpticalDetector/OpHitFinder/PMTPulseRecoBase.h"
#include "larana/OpticalDetector/OpHitFinder/AlgoThreshold.h"
#include "larana/OpticalDetector/OpHitFinder/AlgoSiPM.h"
#include "larana/OpticalDetector/OpHitFinder/AlgoSlidingWindow.h"
#include "larana/OpticalDetector/OpHitFinder/AlgoFixedWindow.h"
#include "larana/OpticalDetector/OpHitFinder/AlgoCFD.h"
#include "larana/OpticalDetector/OpHitFinder/PedAlgoEdges.h"
#include "larana/OpticalDetector/OpHitFinder/PedAlgoRollingMean.h"
#include "larana/OpticalDetector/OpHitFinder/PedAlgoUB.h"
#include "larana/OpticalDetector/OpHitFinder/PulseRecoManager.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "larana/OpticalDetector/OpHitFinder/OpHitAlg.h"

#include "canvas/Utilities/Exception.h"

#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"


#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"
#include "lardataobj/RawData/OpDetWaveform.h"

#include "TTree.h"


namespace pmtcalo {
  class PMTMonitor;
}

class pmtcalo::PMTMonitor : public art::EDAnalyzer {

public:


  // the struct hold channels metadata
  struct channelMetadata
  {
    raw::Channel_t channel_id;
    size_t fragment_id;
    size_t digitizer_id;
    raw::TimeStamp_t timestamp;
    uint16_t median;
    double rms;
    size_t npulses;
  };
  
  explicit PMTMonitor(fhicl::ParameterSet const& pset);

  PMTMonitor(PMTMonitor const&) = delete;
  PMTMonitor(PMTMonitor&&) = delete;
  PMTMonitor& operator=(PMTMonitor const&) = delete;
  PMTMonitor& operator=(PMTMonitor&&) = delete;

  virtual void beginJob() override;

  uint16_t findMedian(std::vector<uint16_t> a, size_t n);
  double   findMedian(std::vector<double> a, size_t n);

  void process_fragment(const artdaq::Fragment &artdaqFragment);

  void analyze(art::Event const& event) override;

  void clean();

private:

  const icarusDB::IICARUSChannelMap* fChannelMap;
  pmtana::PulseRecoManager  fPulseRecoMgr;
  pmtana::PMTPulseRecoBase* fThreshAlg;
  pmtana::PMTPedestalBase*  fPedAlg;

  std::string fInputTag;

  std::vector<channelMetadata> m_channel_metadata;
  std::map<int, double> m_temperature_fragment;

  int m_run;
  int m_event;
  std::vector<size_t> m_fragment_id;
  std::vector<size_t> m_digitizer_id;
  std::vector<raw::Channel_t> m_channel_id;
  std::vector<raw::TimeStamp_t> m_timestamp;
  std::vector<uint16_t> m_baseline;
  std::vector<double> m_rms;
  std::vector<size_t> m_npulses;
  std::vector<double> m_temperature;

  TTree *m_om_ttree;

  art::ServiceHandle<art::TFileService> tfs;

};


//------------------------------------------------------------------------------


pmtcalo::PMTMonitor::PMTMonitor(fhicl::ParameterSet const& pset)
  : art::EDAnalyzer(pset),
  fPulseRecoMgr()  // ,
{
   fInputTag = pset.get<std::string>("FragmentLabel", "");

   // Pulse finder algorithm
   fChannelMap = art::ServiceHandle<icarusDB::IICARUSChannelMap const>{}.get();


   // Initialize the hit finder algorithm
  auto const hit_alg_pset = pset.get< fhicl::ParameterSet >("HitAlgoPset");
  std::string threshAlgName = hit_alg_pset.get< std::string >("Name");
  if      (threshAlgName == "Threshold")
    fThreshAlg = new pmtana::AlgoThreshold(hit_alg_pset);
  else if (threshAlgName == "SiPM")
    fThreshAlg = new pmtana::AlgoSiPM(hit_alg_pset);
  else if (threshAlgName == "SlidingWindow")
    fThreshAlg = new pmtana::AlgoSlidingWindow(hit_alg_pset);
  else if (threshAlgName == "FixedWindow")
    fThreshAlg = new pmtana::AlgoFixedWindow(hit_alg_pset);
  else if (threshAlgName == "CFD" )
    fThreshAlg = new pmtana::AlgoCFD(hit_alg_pset);
  else throw art::Exception(art::errors::UnimplementedFeature)
                    << "Cannot find implementation for "
                    << threshAlgName << " algorithm.\n";

  auto const ped_alg_pset = pset.get< fhicl::ParameterSet >("PedAlgoPset");
  std::string pedAlgName = ped_alg_pset.get< std::string >("Name");
  if      (pedAlgName == "Edges")
    fPedAlg = new pmtana::PedAlgoEdges(ped_alg_pset);
  else if (pedAlgName == "RollingMean")
    fPedAlg = new pmtana::PedAlgoRollingMean(ped_alg_pset);
  else if (pedAlgName == "UB"   )
    fPedAlg = new pmtana::PedAlgoUB(ped_alg_pset);
  else throw art::Exception(art::errors::UnimplementedFeature)
                    << "Cannot find implementation for "
                    << pedAlgName << " algorithm.\n";

  fPulseRecoMgr.AddRecoAlgo(fThreshAlg);
  fPulseRecoMgr.SetDefaultPedAlgo(fPedAlg);

}


//------------------------------------------------------------------------------

// Function overloading, I'm not fancy enough to use templates

uint16_t pmtcalo::PMTMonitor::findMedian(std::vector<uint16_t> a, size_t n) {
      // First we sort the array
      sort(a.begin(), a.end());

      // check for even case
      if (n % 2 != 0)
          return (uint16_t)a[n / 2];

      return (uint16_t)(a[(n - 1) / 2] + a[n / 2]) / 2.0;
  }

double pmtcalo::PMTMonitor::findMedian(std::vector<double> a, size_t n) {
      // First we sort the array
      sort(a.begin(), a.end());

      // check for even case
      if (n % 2 != 0)
          return a[n / 2];

      return a[(n - 1) / 2] + a[n / 2] / 2.0;
  }

//------------------------------------------------------------------------------


void pmtcalo::PMTMonitor::process_fragment(const artdaq::Fragment &artdaqFragment)
{
    size_t fragment_id = artdaqFragment.fragmentID()& 0x0fff;

    
    // convert fragment to Nevis fragment
    sbndaq::CAENV1730Fragment         fragment(artdaqFragment);
    sbndaq::CAENV1730FragmentMetadata metafrag = *fragment.Metadata();
    sbndaq::CAENV1730Event            evt      = *fragment.Event();
    sbndaq::CAENV1730EventHeader      header   = evt.Header;

    size_t nChannelsPerBoard  = metafrag.nChannels; //fragment.nChannelsPerBoard();

    uint32_t ev_size_quad_bytes         = header.eventSize;
    uint32_t evt_header_size_quad_bytes = sizeof(sbndaq::CAENV1730EventHeader)/sizeof(uint32_t);
    uint32_t data_size_double_bytes     = 2*(ev_size_quad_bytes - evt_header_size_quad_bytes);
    uint32_t nSamplesPerChannel         = data_size_double_bytes/nChannelsPerBoard;

    raw::TimeStamp_t time_tag =  metafrag.timeStampSec;

    const uint16_t* data_begin = reinterpret_cast<const uint16_t*>(artdaqFragment.dataBeginBytes() + sizeof(sbndaq::CAENV1730EventHeader));
    const uint16_t* value_ptr  =  data_begin;
    uint16_t        value      = 0;
    size_t          ch_offset  = 0;

    // Recover the information for this fragment
    if (fChannelMap->hasPMTDigitizerID(fragment_id))
    {
        const icarusDB::DigitizerChannelChannelIDPairVec& digitizerChannelVec = fChannelMap->getChannelIDPairVec(fragment_id);

        // Allocate the vector outside the loop just since we'll resuse it over and over...
        std::vector<uint16_t> wvfm(nSamplesPerChannel);


        uint32_t channel_temperature = 0;

        for(const auto& digitizerChannelPair : digitizerChannelVec)
        {
            size_t         digitizerChannel = digitizerChannelPair.first;
            raw::Channel_t channelID        = digitizerChannelPair.second;


            channel_temperature += metafrag.chTemps[digitizerChannel];

	    	    // Actual mapping is from 0 to 15 with 15 active channel and a last inactive channel. Database mapping has channels numbered from 1 to 15 without the spare channel. Subtract 1 to find the correct association
            ch_offset = digitizerChannel * nSamplesPerChannel;

            uint16_t min = 17000;

            for(size_t i_t=0; i_t < nSamplesPerChannel; ++i_t)
            {
                value_ptr = data_begin + ch_offset + i_t; 
                value     = *(value_ptr);
  	            wvfm[i_t] = value;

                if( value < min ){
                  min=value;
                }
            }

            // If the waveform saturates then remove it 

            if( min < 1000 ){
              continue;
            }

            uint16_t median = findMedian(wvfm, nSamplesPerChannel);

            //fPulseRecoMgr wants a pmtana::Waveform_t type 
            pmtana::Waveform_t waveform(wvfm.begin(), wvfm.end());

            // Here we geopardize the same algorithm behind the OpHit reconstruction
            fPulseRecoMgr.Reconstruct(waveform);

            //_ped_mean_v = fPedAlg->Mean();
            std::vector<double>  pedestal_sigma = fPedAlg->Sigma();
            double rms = findMedian(pedestal_sigma, pedestal_sigma.size());

            auto const& pulses = fThreshAlg->GetPulses();
            size_t npulses = pulses.size();

            channelMetadata chmd = { channelID, fragment_id, digitizerChannel, time_tag, median, rms, npulses };
            m_channel_metadata.push_back( chmd );

        }

        // Make the mean of the temperature of the four boards on the digitizer
        m_temperature_fragment[fragment_id] = (double)channel_temperature/(double)nChannelsPerBoard;

    }
    else std::cout << "*** PMT could not find channel information for fragment: " << fragment_id << std::endl;

    return;
}



//------------------------------------------------------------------------------


void pmtcalo::PMTMonitor::beginJob()
{

  m_om_ttree = tfs->make<TTree>("pmtmonitor", "ttree to monitor pmt system");

  m_om_ttree->Branch("run", &m_run, "run/I" );
  m_om_ttree->Branch("event", &m_event, "event/I" );
  m_om_ttree->Branch("channel_id", &m_channel_id );
  m_om_ttree->Branch("fragment_id", &m_fragment_id );
  m_om_ttree->Branch("digitizer_id", &m_digitizer_id );
  m_om_ttree->Branch("timestamp", &m_timestamp );
  m_om_ttree->Branch("baseline", &m_baseline );
  m_om_ttree->Branch("rms", &m_rms );
  m_om_ttree->Branch("temperature", &m_temperature );
  m_om_ttree->Branch("npulses", &m_npulses);
 
}


//-----------------------------------------------------------------------------


void pmtcalo::PMTMonitor::analyze(art::Event const& event)
{
	
  m_run = event.id().run();
	m_event = event.id().event();
  
	// Protect for runs with no PMT info
	try
	{
		// Recover the data fragments for the PMT 
		auto const& daq_handle = event.getValidHandle<artdaq::Fragments>(fInputTag);
		// Make sure data available
        
        if (daq_handle.isValid() && daq_handle->size() > 0)
        {
            for (auto const &rawFrag: *daq_handle)
            {
            	process_fragment(rawFrag);
            }
        }
	}
	catch(...)
	{
		std::cout << "pmtcalo::PMTMonitor::analyze: Did not find daq data products to decode" << std::endl;
    return;
	}


  for( auto channel_metadata : m_channel_metadata )
  {

    m_fragment_id.push_back( channel_metadata.fragment_id );
    m_digitizer_id.push_back( channel_metadata.digitizer_id );
    m_channel_id.push_back( channel_metadata.channel_id );
    m_timestamp.push_back( channel_metadata.timestamp );
    m_baseline.push_back( channel_metadata.median );
    m_rms.push_back( channel_metadata.rms );
    m_npulses.push_back( channel_metadata.npulses );
    m_temperature.push_back( m_temperature_fragment[channel_metadata.fragment_id] );

  }
  
  
  m_om_ttree->Fill();

  clean();

 
} // end analyze


//-----------------------------------------------------------------------------


void pmtcalo::PMTMonitor::clean()
{

  m_channel_metadata.clear();
  m_temperature_fragment.clear();
  m_fragment_id.clear();
  m_digitizer_id.clear();
  m_channel_id.clear();
  m_timestamp.clear();
  m_baseline.clear();
  m_rms.clear();
  m_npulses.clear();
  m_temperature.clear();

}



DEFINE_ART_MODULE(pmtcalo::PMTMonitor)
