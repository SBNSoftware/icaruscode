/**
 * @file   DummyFilter_module.cc
 * @brief  An _art_ filter which does nothing.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   July 15, 2021
 *
 */

// framework libraries
#include "art/Framework/Core/SharedFilter.h"
#include "art/Framework/Core/ModuleMacros.h"


//------------------------------------------------------------------------------
/**
 * @brief Filter module which passes everything.
 * @see DummyProducer
 *
 * This module accepts any FHiCL configuration (which it ignores), acts like a
 * filter that passes everything, does not consume anything and does not produce
 * anything.
 * 
 * It is a "good" way to remove a filter module from a job without having to
 * mess with much configuration: just target it with a:
 *     
 *     physics.producers.redundantmod.module_type: DummyFilter
 *     
 *
 * Configuration parameters
 * =========================
 *
 * Anything.
 *
 */
class DummyFilter: public art::SharedFilter {
  
    public:

  DummyFilter(fhicl::ParameterSet const& params, const art::ProcessingFrame&)
    : art::SharedFilter{ params }
    { async<art::InEvent>(); }
  
  virtual bool filter(art::Event&, const art::ProcessingFrame&) override
    { return true; }
  
}; // class DummyFilter


DEFINE_ART_MODULE(DummyFilter)

//------------------------------------------------------------------------------

