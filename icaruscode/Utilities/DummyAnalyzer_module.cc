/**
 * @file   DummyAnalyzer_module.cc
 * @brief  An _art_ analyzer which does nothing.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   July 26, 2022
 *
 */

// framework libraries
#include "art/Framework/Core/SharedAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"


//------------------------------------------------------------------------------
/**
 * @brief Module which does nothing.
 * @see DummyFilter, DummyProducer
 *
 * This module accepts any FHiCL configuration (which it ignores), acts like an
 * analyzer that does not consume anything.
 * 
 * It is a "good" way to remove an analyzer module from a job without having to
 * mess with much configuration: just target it with a:
 *     
 *     physics.analyzer.redundantanalyzer.module_type: DummyAnalyzer
 *     
 *
 * Configuration parameters
 * =========================
 *
 * Anything.
 *
 */
class DummyAnalyzer: public art::SharedAnalyzer {
  
    public:

  DummyAnalyzer(fhicl::ParameterSet const& params, const art::ProcessingFrame&)
    : art::SharedAnalyzer{ params }
    { async<art::InEvent>(); }
  
  virtual void analyze(art::Event const&, const art::ProcessingFrame&) override
    {}
  
}; // class DummyAnalyzer


DEFINE_ART_MODULE(DummyAnalyzer)

//------------------------------------------------------------------------------

