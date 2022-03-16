////////////////////////////////////////////////////////////////////////
// Class:       ICARUSTimeChecksDQM
// Plugin Type: analyzer (Unknown Unknown)
// File:        ICARUSTimeChecksDQM_module.cc
//
// Generated at Wed Mar 16 04:45:48 2022 by Wesley Ketchum using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "sbndaq-artdaq-core/Overlays/ICARUS/ICARUSTriggerUDPFragment.hh"
#include "sbndaq-artdaq-core/Overlays/ICARUS/PhysCrateFragment.hh"

#include <map>
#include <iostream>

namespace icarus {
  class ICARUSTimeChecksDQM;
}


class icarus::ICARUSTimeChecksDQM : public art::EDAnalyzer {

typedef std::pair<uint16_t,uint16_t> tpc_id_t; //fragment id, board id

public:
  explicit ICARUSTimeChecksDQM(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ICARUSTimeChecksDQM(ICARUSTimeChecksDQM const&) = delete;
  ICARUSTimeChecksDQM(ICARUSTimeChecksDQM&&) = delete;
  ICARUSTimeChecksDQM& operator=(ICARUSTimeChecksDQM const&) = delete;
  ICARUSTimeChecksDQM& operator=(ICARUSTimeChecksDQM&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:

  // Declare member data here.
  std::vector<art::InputTag> const m_tpc_data_labels;
  art::InputTag              const m_trigger_data_label;
  bool                       const m_verbose;

  size_t const m_n_tpc_collections;

  std::map <tpc_id_t, uint32_t>  m_tpcid_to_timestamp_map;
  std::map <uint32_t, std::vector<tpc_id_t> > m_timestamp_to_tpcid_map;
  std::map <tpc_id_t, uint32_t>  m_tpcid_to_event_map;
  std::map <uint32_t, std::vector<tpc_id_t> > m_event_to_tpcid_map;
};


icarus::ICARUSTimeChecksDQM::ICARUSTimeChecksDQM(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
  , m_tpc_data_labels(p.get<std::vector<art::InputTag>>("TPCDataLabels"))
  , m_trigger_data_label(p.get<art::InputTag>("TriggerDataLabel"))
  , m_verbose(p.get<bool>("Verbose",true))
  , m_n_tpc_collections(m_tpc_data_labels.size())  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void icarus::ICARUSTimeChecksDQM::analyze(art::Event const& e)
{
  //clear maps
  m_tpcid_to_event_map.clear();
  m_tpcid_to_timestamp_map.clear();
  m_timestamp_to_tpcid_map.clear();
  m_event_to_tpcid_map.clear();

  //create our handles
  std::vector< art::Handle< std::vector<artdaq::Fragment> > > tpc_data_handles(m_n_tpc_collections);
  art::Handle< std::vector<artdaq::Fragment> > trigger_data_handle;

  //fill them by the labels
  for(size_t itpc=0; itpc<m_n_tpc_collections; ++itpc)
    e.getByLabel(m_tpc_data_labels[itpc],tpc_data_handles[itpc]);
  e.getByLabel(m_trigger_data_label,trigger_data_handle);

  //process trigger information
  auto const& trigger_data_vec(*trigger_data_handle);
  if(trigger_data_vec.size()!=1){
    MF_LOG_WARNING("ICARUSTimeChecksDQM")
      << "Trigger data vector of unexpected size "
      << trigger_data_vec.size();
  }

  //icarus::ICARUSTriggerUDPFragment trig_frag(trigger_data_vec.at(0));
  //char *buffer = const_cast<char*>(trig_frag.GetDataString().c_str());
  //icarus::ICARUSTriggerInfo datastream_info = icarus::parse_ICARUSTriggerString(buffer);

  std::string trig_data((char*)trigger_data_vec[0].dataBeginBytes(),
                        trigger_data_vec[0].dataSizeBytes());

  std::cout << "Get Trigger Info" << std::endl;

  //uint64_t wr_ts = 0;//datastream_info.getNanoseconds_since_UTC_epoch() + 2e9;
  //auto ev_num_trig = trig_frag.getEventNo();
  auto ev_num_frag_trig = trigger_data_vec[0].sequenceID();
  //int64_t time_trig_ns = trig_frag.getSeconds()*1e9 + trig_frag.getNanoSeconds();


  //process TPC information
  uint16_t tpc_frag_id,tpc_board_id;
  tpc_id_t tpc_id;
  uint32_t tpc_timestamp;
  uint32_t tpc_event;
  for(auto const& tpc_handle : tpc_data_handles){
    auto const& tpc_data_vec(*tpc_handle);

    for(auto const& frag : tpc_data_vec){
        tpc_frag_id = frag.fragmentID();
        icarus::PhysCrateFragment tpc_frag(frag);
        for(tpc_board_id=0; tpc_board_id<tpc_frag.nBoards(); ++tpc_board_id){
          tpc_id = std::make_pair(tpc_frag_id,tpc_board_id);
          tpc_timestamp = tpc_frag.BoardTimeStamp(tpc_board_id);
          tpc_event = tpc_frag.BoardEventNumber(tpc_board_id);

          m_tpcid_to_timestamp_map[tpc_id] = tpc_timestamp;
          m_timestamp_to_tpcid_map[tpc_timestamp].push_back(tpc_id);
          m_tpcid_to_event_map[tpc_id] = tpc_event;
          m_event_to_tpcid_map[tpc_event].push_back(tpc_id);
        }
    }
  }

  //print out information
  if(m_verbose){

    std::cout << "Trigger seqID=" << ev_num_frag_trig << std::endl;
    //std::cout << "Trigger timestamps (wr,hw): (" << wr_ts << "," << time_trig_ns << ")" << std::endl;
    std::cout << "Trigger string: " << trig_data << std::endl;

    std::cout << "Number of TPC boards: " << m_tpcid_to_event_map.size() << std::endl;

    std::cout << "Event numbers of TPC boards:";
    for( auto const& tpc_evs : m_event_to_tpcid_map){
      std::cout << "\n\t" << tpc_evs.first
                << " (" << tpc_evs.second.size() << " boards)";
    }
    std::cout << std::endl;
    std::cout << "Timestamps of TPC boards:";
    for( auto const& tpc_ts : m_timestamp_to_tpcid_map){
      std::cout << "\n\t" << tpc_ts.first
                << " (" << tpc_ts.second.size() << " boards)";
    }
    std::cout << std::endl;

  }
}

DEFINE_ART_MODULE(icarus::ICARUSTimeChecksDQM)
