////////////////////////////////////////////////////////////////////////
// Class:       PMTMonitor
// Plugin Type: analyzer (art v3_05_00)
// File:        PMTMonitor_module.cc
//
// Generated at Mon Sep 21 15:21:37 2020 by Andrea Scarpelli
// 
// Base for the PMT Online Monitor
//
//  mailto:ascarpell@bnl.gov
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
    //float baseline;
    //float rms;
    //int npulses;
  };


  explicit PMTMonitor(fhicl::ParameterSet const& pset);

  PMTMonitor(PMTMonitor const&) = delete;
  PMTMonitor(PMTMonitor&&) = delete;
  PMTMonitor& operator=(PMTMonitor const&) = delete;
  PMTMonitor& operator=(PMTMonitor&&) = delete;

  virtual void beginJob() override;

  void process_fragment(const artdaq::Fragment &artdaqFragment);

  void analyze(art::Event const& event) override;

  //void clean();

private:

  //const icarusDB::IICARUSChannelMap* fChannelMap;

  //std::string fInputTag;

  //std::vector<channelMetadata> m_channel_metadata;

  int m_run;
  int m_event;
  //std::vector<size_t>   *m_fragment_id = NULL;
  //std::vector<size_t>   *m_digitizer_id = NULL;
  //std::vector<raw::Channel_t>   *m_channel_id = NULL;
  //std::vector<raw::TimeStamp_t>   *m_timestamp = NULL;
  //std::vector<int> *m_baseline = NULL;
  //std::vector<float> *m_rms = NULL;
  
  //std::vector<float> *m_temperature = NULL;
  //std::vector<int> *m_pulses_channel = NULL;

  TTree *m_om_ttree;


  art::ServiceHandle<art::TFileService> tfs;

};


//------------------------------------------------------------------------------


pmtcalo::PMTMonitor::PMTMonitor(fhicl::ParameterSet const& pset)
  : art::EDAnalyzer(pset)  // ,
{
   //fInputTag = pset.get<std::string>("FragmentsLabel", "daq:CAEN1730");

   // Pulse finder algorithm
   //fChannelMap = art::ServiceHandle<icarusDB::IICARUSChannelMap const>{}.get();
}


//------------------------------------------------------------------------------

/*
void pmtcalo::PMTMonitor::process_fragment(const artdaq::Fragment &artdaqFragment)
{
    size_t fragment_id = artdaqFragment.fragmentID();

    
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

        for(const auto& digitizerChannelPair : digitizerChannelVec)
        {
            size_t         digitizerChannel = digitizerChannelPair.first;
            raw::Channel_t channelID        = digitizerChannelPair.second;

	    	// Actual mapping is from 0 to 15 with 15 active channel and a last inactive channel. Database mapping has channels numbered from 1 to 15 without the spare channel. Subtract 1 to find the correct association
            ch_offset = digitizerChannel * nSamplesPerChannel;

            for(size_t i_t=0; i_t < nSamplesPerChannel; ++i_t)
            {
                value_ptr = data_begin + ch_offset + i_t; 
                value     = *(value_ptr);
  	            wvfm[i_t] = value;
            }

            channelMetadata chmd = { channelID, fragment_id, digitizerChannel, time_tag };
            m_channel_metadata.push_back( chmd );

            // time_tag, channelID, wvfm[i_t]
            //raw::TimeStamp_t time_tag = 0.0; //metafrag.timeStampSec;

            // Baseline 
            // RMS 

            // Temperature

            // N Pulses

        }


        // Make the mean of the temperature of the four boards on the digitizer

    }
    else std::cout << "*** PMT could not find channel information for fragment: " << fragment_id << std::endl;

    return;
}
*/


//------------------------------------------------------------------------------


void pmtcalo::PMTMonitor::beginJob()
{

  
  m_om_ttree = tfs->make<TTree>("monitor", "ttree to monitor pmt system");

  m_om_ttree->Branch("run", &m_run, "run/I" );
  m_om_ttree->Branch("event", &m_event, "event/I" );
  //m_om_ttree->Branch("channel_id", &m_channel_id );
  //m_om_ttree->Branch("channel_id", &m_fragment_id );
  //m_om_ttree->Branch("channel_id", &m_digitizer_id );
  //m_om_ttree->Branch("timestamp", &m_timestamp );
  //m_om_ttree->Branch("baseline", &m_baseline );
  //m_om_ttree->Branch("rms", &m_rms );
  //m_om_ttree->Branch("temperature", &m_temperature );
  //m_om_ttree->Branch("pulses_channel", &m_pulses_channel );
 
}


//-----------------------------------------------------------------------------


void pmtcalo::PMTMonitor::analyze(art::Event const& event)
{
	m_run = event.id().run();
	m_event = event.id().event();

  /*
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

    m_fragment_id->push_back( channel_metadata.fragment_id );
    m_digitizer_id->push_back( channel_metadata.digitizer_id );
    m_channel_id->push_back( channel_metadata.channel_id );
    m_timestamp->push_back( channel_metadata.timestamp );

  }
  */

  m_om_ttree->Write();

  

  //clean();

} // end analyze


//-----------------------------------------------------------------------------

/*
void pmtcalo::PMTMonitor::clean()
{

  m_channel_metadata.clear();
  m_fragment_id->clear();
  m_digitizer_id->clear();
  m_channel_id->clear();
  m_timestamp->clear();

}
*/


DEFINE_ART_MODULE(pmtcalo::PMTMonitor)
