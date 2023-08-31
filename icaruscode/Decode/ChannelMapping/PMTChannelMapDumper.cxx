/**
 * @file   icaruscode/Decode/ChannelMapping/PMTChannelMapDumper.cxx
 * @brief  Utility dumping the content of PMT channel mapping on screen.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * 
 * This utility can be run with any configuration file including a configuration
 * for `IICARUSChannel` service.
 * 
 * It is using _art_ facilities for tool loading, but it does not run in _art_
 * environment. So it may break without warning and without solution.
 * 
 */


// ICARUS libraries
#include "icaruscode/Decode/ChannelMapping/ICARUSChannelMapProvider.h"
#include "icaruscode/Decode/ChannelMapping/RunPeriods.h"

// LArSoft and framework libraries
#include "larcorealg/CoreUtils/enumerate.h"
#include "larcorealg/TestUtils/unit_test_base.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// C/C++ standard libraries
#include <iomanip> // std::setw()
#include <iostream>
#include <algorithm>
#include <numeric> // std::iota()
#include <array>


// -----------------------------------------------------------------------------
template <std::size_t KeyNo = 0U>
struct SortByElement {
  
  template <typename TupleA, typename TupleB>
  bool operator() (TupleA const& A, TupleB const& B) const
    { return less(A, B); }
  
  template <typename Tuple>
  static decltype(auto) key(Tuple const& t) { return std::get<KeyNo>(t); }
  
  template <typename TupleA, typename TupleB>
  static bool less(TupleA const& A, TupleB const& B)
    { using std::less; return less{}(key(A), key(B)); }
  
}; // SortByElement<>


// -----------------------------------------------------------------------------
void dumpMapping(icarusDB::ICARUSChannelMapProvider const& channelMapping) {
  
  // hard-coded list of fragment ID; don't like it?
  // ask for an extension of the channel mapping features.
  std::array<unsigned int, 24U> FragmentIDs;
  std::iota(FragmentIDs.begin(), FragmentIDs.end(), 0x2000);
  
  mf::LogVerbatim("PMTchannelMappingDumper") << "Fragment IDs:";
  for (auto const [ iFragment, fragmentID]: util::enumerate(FragmentIDs)) {
    unsigned int const effFragmentID = fragmentID & 0xFF;
    
    if (!channelMapping.hasPMTDigitizerID(effFragmentID)) {
      mf::LogVerbatim("PMTchannelMappingDumper")
        << "[" << iFragment << "] " << std::hex << fragmentID << std::dec
        << " not found in the database (as "
        << std::hex << effFragmentID << std::dec << ")";
      continue;
    }
    
    icarusDB::DigitizerChannelChannelIDPairVec digitizerChannels
      = channelMapping.getChannelIDPairVec(effFragmentID);
    
    std::sort
      (digitizerChannels.begin(), digitizerChannels.end(), SortByElement<1U>{});
    
    
    mf::LogVerbatim log("PMTchannelMappingDumper");
    log
      << "[" << iFragment << "] " << std::hex << fragmentID << std::dec
      << " includes " << digitizerChannels.size()
      << " LArSoft channels between " << std::get<1U>( digitizerChannels.front() )
      << " and " << std::get<1U>( digitizerChannels.back() )
      << " [board channel index in brackets]:";
    constexpr unsigned int Cols = 8U;
    unsigned int n = 0;
    for(auto const [ digitizerChannel, channelID, laserChannel ]: digitizerChannels) {
      if (n-- == 0) { log << "\n     "; n = Cols - 1U; }
      log << " "  << std::setw(3) << channelID
        << " [" << std::setw(3) << digitizerChannel << "]";
    } // for channel
    
  } // for fragment
  
}


// -----------------------------------------------------------------------------
int main(int argc, char** argv) {
  
  using Environment
    = testing::TesterEnvironment<testing::BasicEnvironmentConfiguration>;
  
  testing::BasicEnvironmentConfiguration config("PMTchannelMappingDumper");

  //
  // parameter parsing
  //
  int iParam = 0;

  // first argument: configuration file (mandatory)
  if (++iParam < argc)
    config.SetConfigurationPath(argv[iParam]);
  else {
    std::cerr << "FHiCL configuration file path required as first argument!"
      << std::endl;
    return 1;
  }

  Environment const Env { config };
  
  icarusDB::ICARUSChannelMapProvider channelMapping
    { Env.ServiceParameters("IICARUSChannelMap") };
  
  for (icarusDB::RunPeriod const period: icarusDB::RunPeriods::All)
  {
    mf::LogVerbatim("PMTchannelMappingDumper")
      << std::string(80, '=') << "\nRun period #"
      << static_cast<unsigned int>(period) << ":\n";
    channelMapping.forPeriod(period);
    dumpMapping(channelMapping);
  }
  
  mf::LogVerbatim("PMTchannelMappingDumper")
    << std::string(80, '=')
    << "\nDumped " << icarusDB::RunPeriods::All.size() << " run periods."
    ;
  
  return 0;
} // main()


