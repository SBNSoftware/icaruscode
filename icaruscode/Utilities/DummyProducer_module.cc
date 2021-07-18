/**
 * @file   DummyProducer_module.cc
 * @brief  An _art_ producer which does nothing.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   July 15, 2021
 *
 */

// framework libraries
#include "art/Framework/Core/SharedProducer.h"
#include "art/Framework/Core/ModuleMacros.h"


//------------------------------------------------------------------------------
/**
 * @brief Module which does nothing.
 * @see DummyFilter
 *
 * This module accepts any FHiCL configuration (which it ignores), acts like a
 * producer that does not consume anything and does not produce anything.
 * 
 * It is a "good" way to remove a producer module from a job without having to
 * mess with much configuration: just target it with a:
 *     
 *     physics.producers.redundantfilter.module_type: DummyProducer
 *     
 *
 * Configuration parameters
 * =========================
 *
 * Anything.
 *
 */
class DummyProducer: public art::SharedProducer {
  
    public:

  DummyProducer(fhicl::ParameterSet const& params, const art::ProcessingFrame&)
    : art::SharedProducer{ params }
    { async<art::InEvent>(); }
  
  virtual void produce(art::Event&, const art::ProcessingFrame&) override {}
  
}; // class DummyProducer


DEFINE_ART_MODULE(DummyProducer)

//------------------------------------------------------------------------------

