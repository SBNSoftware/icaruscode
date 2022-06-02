/**
 * @file   PMTconfigurationExtraction_module.cc
 * @brief  Writes PMT configuration from FHiCL into data product.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   February 23, 2021
 */

// ICARUS/SBN libraries
#include "sbnobj/Common/PMT/Data/PMTconfiguration.h" // sbn::PMTconfiguration
#include "icaruscode/Decode/ChannelMapping/IICARUSChannelMap.h"
#include "icaruscode/PMT/Timing/PMTTimeCorrection.h" // data product holding the corrections

// framework libraries
#include "canvas/Persistency/Provenance/ProcessConfiguration.h"
#include "canvas/Persistency/Provenance/ProcessHistory.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"

#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/FileBlock.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "cetlib_except/exception.h"

// LArSoft libraries
#include "lardataobj/RawData/OpDetWaveform.h"

// C/C++ standard libraries
#include <memory> // std::unique_ptr<>
#include <optional>
#include <string>
#include <utility> // std::move()
#include <cassert>

// -----------------------------------------------------------------------------
namespace icarus { class TimingCorrectionExtraction; }
/**
 * @brief 
 * 
 * This module reads 
 * 
 * Input
 * ------
 * 
 * 
 * Output
 * -------
 * 
 * 
 * 
 * Configuration parameters
 * -------------------------
 * 
 * 
 * 
 * Multithreading
 * ---------------
 * 
 * 
 */
class icarus::TimingCorrectionExtraction: public art::EDProducer {
  
  /// Pointer to the online channel mapping service.
  icarusDB::IICARUSChannelMap const* fChannelMap = nullptr;

  bool fVerbose = false; ///< Whether to print the configuration we read.
  
  std::string fLogCategory; ///< Category tag for messages.
  
  public:
  
  /// Configuration of the module.
  struct Config {

    fhicl::Sequence<art::InputTag> InputLabels {
        fhicl::Name("InputLabels"), 
        fhicl::Comment("list of the input lables to be used")
    };

    fhicl::Atom<bool> Verbose {
      fhicl::Name("Verbose"),
      fhicl::Comment("print the times read and the associated channels"),
      false // default
    };
    
    fhicl::Atom<std::string> LogCategory {
      fhicl::Name("LogCategory"),
      fhicl::Comment("category tag used for messages to message facility"),
      "TimingCorrectionExtraction" // default
    };
    
  }; // struct Config

  using Parameters = art::EDProducer::Table<Config>;

  /// Constructor: just reads the configuration.
  TimingCorrectionExtraction(Parameters const& config);
    
  /// Mandatory method, unused.
  virtual void produce(art::Event&) override {}
  
}; // icarus::TimingCorrectionExtractor



// -----------------------------------------------------------------------------
icarus::TimingCorrectionExtraction( Parameters const& config ) 
: art::EDProducer(config)
, fChannelMap(config().AssignOfflineChannelIDs()
    ? art::ServiceHandle<icarusDB::IICARUSChannelMap const>{}.get()
    : nullptr
  )
, fVerbose(config().Verbose())
, fLogCategory(config().LogCategory())
{

    for ( auto const & tag : fInputLabels )
        consumes<std::vector<raw::OpDetWaveform>>(tag);

    produces<icarus::PMTTimeCorrection>();

};

// -----------------------------------------------------------------------------


// -----------------------------------------------------------------------------


// -----------------------------------------------------------------------------
icarus::TimingCorrectionExtraction::produce( art::Event& ) {


};