/**
 * @file   WriteAdderChannelMap_module.cc
 * @brief  Writes the adder mapping as a data product.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   December 4, 2025
 */


// ICARUS libraries
#include "icaruscode/PMT/Trigger/Algorithms/AdderChannelMaps.h"
#include "icaruscode/Decode/ChannelMapping/IICARUSChannelMap.h"
#include "icaruscode/Decode/ChannelMapping/IICARUSChannelMapProvider.h"
#include "icaruscode/IcarusObj/ChannelToChannelMap.h"

// // LArSoft libraries
#include "larcorealg/CoreUtils/values.h" // util::const_values()

// framework libraries
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/SharedProducer.h"
#include "art/Framework/Principal/Run.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/Atom.h"

// C/C++ standard libraries
#include <iomanip> // std::setw()
#include <memory> // std::unique_ptr
#include <string>
#include <utility> // std::move()


//------------------------------------------------------------------------------
namespace icarus::trigger { class WriteAdderChannelMap; }
/**
 * @brief Builds and writes as data product the ICARUS adder channel mapping.
 * 
 * This module writes the adder channel mapping as a run-level event.
 * The data product is a simple list of adder channel numbers, each with a list
 * of PMT channels associated to it.
 * The mapping is built on every run extracting information from the ICARUS
 * channel mapping database (see `icarus::trigger::AdderChannelMapBuilder`).
 * 
 * 
 * Output
 * =======
 * 
 * In each run:
 * * `icarus::ChannelToChannelMap<raw::Channel_t>`: the mapping between adder
 *   and PMT channels.
 * 
 * 
 * Service requirements
 * =====================
 * 
 * The following services are _required_:
 * 
 * * _art_ message facility (mapping is dumped in the INFO stream)
 * * `IICARUSChannelMap` (ICARUS hardware channel mapping service)
 * 
 * 
 * Configuration parameters
 * =========================
 * 
 * A terse description of the parameters is printed by running
 * `lar --print-description WriteAdderChannelMap`.
 * 
 * * `LogCategory` (string, default: `"WriteAdderChannelMap"`): label
 *   for the category of messages in the console output; this is the label
 *   that can be used for filtering messages via MessageFacility service
 *   configuration.
 * 
 */
class icarus::trigger::WriteAdderChannelMap: public art::SharedProducer {
  
    public:
  
  /// Type of channel map serialized.
  using ChannelMap_t = icarus::ChannelToChannelMap<raw::Channel_t>;
  
  /// FHiCL configuration of the module.
  struct Config {
    
    fhicl::Atom<std::string> LogCategory {
      fhicl::Name{ "LogCategory" },
      fhicl::Comment{ "message facility stream for logging" },
      "WriteAdderChannelMap"
      };
    
  }; // Config
  
  using Parameters = art::SharedProducer::Table<Config>;
  
  
  WriteAdderChannelMap(Parameters const& params, art::ProcessingFrame const&);
  
  /// Builds and writes the map.
  virtual void beginRun(art::Run& run, art::ProcessingFrame const&) override;
  
  // required
  virtual void produce(art::Event&, art::ProcessingFrame const&) override {}
  
  
    private:
  
  std::string const fLogCategory; ///< Message facility stream name.
  
  icarusDB::IICARUSChannelMapProvider const& fPMTchannelMap;
  
}; // icarus::trigger::WriteAdderChannelMap


//------------------------------------------------------------------------------
//--- Implementation
//------------------------------------------------------------------------------
icarus::trigger::WriteAdderChannelMap::WriteAdderChannelMap
  (Parameters const& params, art::ProcessingFrame const&)
  : art::SharedProducer{ params }
  , fLogCategory{ params().LogCategory() }
  , fPMTchannelMap
    { *art::ServiceHandle<icarusDB::IICARUSChannelMap>()->provider() }
{
  async<art::InEvent>();
  produces<ChannelMap_t, art::InRun>();
}


//------------------------------------------------------------------------------
void icarus::trigger::WriteAdderChannelMap::beginRun
  (art::Run& run, art::ProcessingFrame const&)
{
  
  AdderChannelMapBuilder const channelMapBuilder{ fLogCategory };
  AdderChannelMap adderChannelMap
    = channelMapBuilder.build(fPMTchannelMap, run.run());
  
  ChannelMap_t map;
  mf::LogInfo log{ fLogCategory };
  log << "Mapping of " << adderChannelMap.size() << " adder channels:";
  for (auto const& adderInfo: util::const_values(adderChannelMap)) {
    log << "\n  " << adderInfo.channel << " => " << adderInfo.PMTs.size()
      << " PMT:";
    for (raw::Channel_t const PMTchannel: adderInfo.PMTs)
      log << std::setw(5) << PMTchannel;
    map.addChannel
      (adderInfo.channel.channel(), { begin(adderInfo.PMTs), end(adderInfo.PMTs) });
  } // for
  
  run.put(std::make_unique<ChannelMap_t>(std::move(map)), art::fullRun());
  
} // icarus::trigger::WriteAdderChannelMap::beginRun()


//------------------------------------------------------------------------------
DEFINE_ART_MODULE(icarus::trigger::WriteAdderChannelMap)


//------------------------------------------------------------------------------
