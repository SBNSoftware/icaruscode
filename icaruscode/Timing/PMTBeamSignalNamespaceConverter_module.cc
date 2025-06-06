////////////////////////////////////////////////////////////////////////
// Class:       PMTBeamSignalNameSpaceConverter_module
// Plugin Type: producer (Unknown Unknown)
// File:        PMTBeamSignalNameSpaceConverter_module.cc
//
// Generated at Thu Jun  5 14:59:32 2025 by Anna Heggestuen using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////
#include "sbnobj/Common/PMT/Data/PMTBeamSignal.hh"
#include "icaruscode/IcarusObj/Legacy/PMTBeamSignal.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>

namespace sbn::timing {
  class PMTBeamSignalNamespaceConverter;
}
/**
 * @brief Copy PMTBeamSignal data product from icaruscode/IcarusObj namespace to sbnobj/Common
 *
 * This producer module simply copies what is stored in the icaruscode/IcarusObj data product 
 * `std::vector<icarus::timing::PMTBeamSignal>` to a new sbnobj/Common data product 
 * `std::vector<sbn::timing::PMTBeamSignal>`. This is to help with backward compatibility for
 * files that are produced with the `icarus::timing::PMTBeamSignal` data product. Converting 
 * the icarus::timing namespace to sbn::timing allows us to fill this variable in the CAFs. 
 * 
 */

class sbn::timing::PMTBeamSignalNamespaceConverter : public art::EDProducer {
public:
  explicit PMTBeamSignalNamespaceConverter(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PMTBeamSignalNamespaceConverter(PMTBeamSignalNamespaceConverter const&) = delete;
  PMTBeamSignalNamespaceConverter(PMTBeamSignalNamespaceConverter&&) = delete;
  PMTBeamSignalNamespaceConverter& operator=(PMTBeamSignalNamespaceConverter const&) = delete;
  PMTBeamSignalNamespaceConverter& operator=(PMTBeamSignalNamespaceConverter&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

   /// RWM waveform instance label
  art::InputTag const fRWMlabel;

  // RWM times
  std::vector<icarus::timing::PMTBeamSignal> fRWMTimes;
  using BeamSignalCollection = std::vector<icarus::timing::PMTBeamSignal>;
};


sbn::timing::PMTBeamSignalNamespaceConverter::PMTBeamSignalNamespaceConverter(fhicl::ParameterSet const& p)
  : EDProducer{p},  
    fRWMlabel(p.get<art::InputTag>("RWMlabel"))
{
  // Call appropriate produces<>() functions here.
  produces<std::vector<sbn::timing::PMTBeamSignal>>("RWM"); 

  // Call appropriate consumes<>() for any products to be retrieved by this module.
  consumes<std::vector<icarus::timing::PMTBeamSignal>>(fRWMlabel);
}

void sbn::timing::PMTBeamSignalNamespaceConverter::produce(art::Event& e)
{
  //old `icarus::timing::PMTBeamSignal` data product
  fRWMTimes = e.getProduct<std::vector<icarus::timing::PMTBeamSignal>>(fRWMlabel);
  
  //new `sbn::timing::PMTBeamSignal` data product
  auto PMTBeamSignalColl = std::make_unique<std::vector<sbn::timing::PMTBeamSignal>>();
  
  if (fRWMTimes.empty())
    mf::LogTrace("ICARUSBeamStructureAna") << "Data product std::vector<icarus::timing::PMTBeamSignal> for '" << fRWMlabel.label()
                                           << "' is empty in event!";
  else{
      for (size_t i = 0; i < fRWMTimes.size(); i++){
        //std::cout << "i: " << i << ", fRWMTimes[i].specialChannel = " << fRWMTimes[i].specialChannel << "\n";
        sbn::timing::PMTBeamSignal newRWMTimes;// = fRWMTimes[i];
        newRWMTimes.specialChannel = fRWMTimes[i].specialChannel;
        newRWMTimes.digitizerLabel = fRWMTimes[i].digitizerLabel;
        newRWMTimes.crate = fRWMTimes[i].crate;
        newRWMTimes.sample = fRWMTimes[i].sample;
        newRWMTimes.startTimeAbs = fRWMTimes[i].startTimeAbs;
        newRWMTimes.startTime = fRWMTimes[i].startTime;
        PMTBeamSignalColl->push_back(newRWMTimes);
      }
  }
  e.put(std::move(PMTBeamSignalColl), "RWM");
} 

DEFINE_ART_MODULE(sbn::timing::PMTBeamSignalNamespaceConverter)
