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

namespace icarus::timing {
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

 * This module is set up to convert either the "RWM" or "EW" instance of the data product 
 * following what is set in the `SignalLabel` input parameter. 
 *
 * Input parameters
 * ------------------------
 * * `SignalLabel` (input tag): tag for the "RWM" or "EW" instance of the 
 *    std::vector<icarus::timing::PMTBeamSignal> data product.
 * 
 * Output products
 * -------------------------
 * This module produces a `std::vector<sbn::timing::PMTBeamSignal>` with 360 elements 
 * representing the "RWM" or "EW" time for the corresponding PMT channel.
 * 
 */

class icarus::timing::PMTBeamSignalNamespaceConverter : public art::EDProducer {
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

  // RWM or EW waveform instance label from fcl 
  art::InputTag const fSignalLabel;

  // Vector for input data product 
  std::vector<icarus::timing::PMTBeamSignal> fPMTBeamSignal;
};


icarus::timing::PMTBeamSignalNamespaceConverter::PMTBeamSignalNamespaceConverter(fhicl::ParameterSet const& p)
  : EDProducer{p},  
    fSignalLabel(p.get<art::InputTag>("SignalLabel"))
{
  // Call appropriate produces<>() functions here.
  produces<std::vector<sbn::timing::PMTBeamSignal>>(fSignalLabel.instance()); 

  // Call appropriate consumes<>() for any products to be retrieved by this module.
  consumes<std::vector<icarus::timing::PMTBeamSignal>>(fSignalLabel);
}
 
void icarus::timing::PMTBeamSignalNamespaceConverter::produce(art::Event& e)
{
  //old `icarus::timing::PMTBeamSignal` data product
  fPMTBeamSignal = e.getProduct<std::vector<icarus::timing::PMTBeamSignal>>(fSignalLabel);
  
  //new `sbn::timing::PMTBeamSignal` data product
  auto PMTBeamSignalColl = std::make_unique<std::vector<sbn::timing::PMTBeamSignal>>();
  
  if (fPMTBeamSignal.empty())
    mf::LogTrace("ICARUSBeamStructureAna") << "Data product std::vector<icarus::timing::PMTBeamSignal> for '" << fSignalLabel.label()
                                           << "' is empty in event!";
  else{
      for (size_t i = 0; i < fPMTBeamSignal.size(); i++){
        sbn::timing::PMTBeamSignal newPMTBeamSignal;
        newPMTBeamSignal.specialChannel = fPMTBeamSignal[i].specialChannel;
        newPMTBeamSignal.digitizerLabel = fPMTBeamSignal[i].digitizerLabel;
        newPMTBeamSignal.crate = fPMTBeamSignal[i].crate;
        newPMTBeamSignal.sample = fPMTBeamSignal[i].sample;
        newPMTBeamSignal.startTimeAbs = fPMTBeamSignal[i].startTimeAbs;
        newPMTBeamSignal.startTime = fPMTBeamSignal[i].startTime;
        PMTBeamSignalColl->push_back(newPMTBeamSignal);
      }
  }
  e.put(std::move(PMTBeamSignalColl), fSignalLabel.instance());
} 

DEFINE_ART_MODULE(icarus::timing::PMTBeamSignalNamespaceConverter)
