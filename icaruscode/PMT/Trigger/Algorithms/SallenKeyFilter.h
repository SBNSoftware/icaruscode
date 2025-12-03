/**
 * @file   icaruscode/PMT/Trigger/Algorithms/SallenKeyFilter.h
 * @brief  Simulation of adder board output.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   August 25, 2025
 * @see    icaruscode/PMT/Trigger/Algorithms/SallenKeyFilter.cxx
 * 
 */

#ifndef ICARUSCODE_PMT_TRIGGER_ALGORITHMS_SALLENKEYFILTER_H
#define ICARUSCODE_PMT_TRIGGER_ALGORITHMS_SALLENKEYFILTER_H

// ICARUS libraries
#include "icaruscode/PMT/Trigger/Algorithms/FilteringCircuit.h"

// LArSoft libraries
#include "lardataalg/Utilities/quantities/electromagnetism.h" // farad, ohm

// framework libraries
#include "fhiclcpp/types/Atom.h"

// C/C++ standard libraries
#include <complex>


// -----------------------------------------------------------------------------
namespace icarus::trigger {
  class SallenKeyFilter;
  class SallenKeyFilterConfig;
}
/**
 * @brief Sallen-Key circuit simulation.
 *
 * The schematics of the ICARUS Sallen Key low-pass filter are:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *                            ..................................................
 *                            :                               [C2]             :
 *                            :                           | |                  :
 *                            :           ,---------------| |---------------.  :
 *                            :           |               | |               |  :
 *               [950 Ω]      :   [R1]    |    [R2]                         |  :
 *    V(in)                   :           |                     |`-._       |  :  [220 Ω]
 *       o---+---'\/\/\,--+---:--'\/\/\,--+--'\/\/\,--+---------| U1 `-._   |  :
 *       |   |            |   :                       |         | LT1818 >--+--:--'\/\/\,--o
 *       |   \            \   :                       |     ,---|    _,-'   |  :         V(out)
 *       |   /            /   :                     -----   |   |_,-'       |  :
 *       |   \ [3 Ω]      \   :                     -----   |               |  :
 *     __|   /            /   :                  [C1] |     |----'\/\/\,----'  :
 *    |  |   \     [50 Ω] \   :                       |     \     [R4]         :
 *    |  |   |            |   :                       |     /                  :
 *    |  |   |            |   :                       |     \ [R3]             :
 *    |  |   |            |   :                       |     /                  :
 *    |  o---+            |   :                       |     \                  :
 *    |                   |   :                       |     |                  :
 * ===========================:================================================:=================
 *                            :             [ ground ]                         :
 *                            :................................................:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * The bracketed part is simulated with the response:
 * @f[
 * \frac{V_{\text{out}}}{V_{\text{in}}} = \frac{k}{a_{2} s^{2} + a_{1} s + 1}
 * @f]
 * where @f$ k = 1 + R_{4}/R{3} @f$,
 * @f$ a_{1} = R_{1} C_{1} + R_{2} C_{1} + R_{1} C_{2} (1 - k)@f$
 * and @f$ a_{2} = R_{1} R{2} C_{1} C_{2} @f$
 *
 * See [SBN DocDB 42951](https://sbn-docdb.fnal.gov/cgi-bin/sso/ShowDocument?docid=42951)
 * for details, including a reference to the derivation of this formula.
 */
class icarus::trigger::SallenKeyFilter: public FilteringCircuit {
    public:
  
  using Data_t = FilteringCircuit::Data_t; // just to reiterate
  
  /// Filter configuration parameters.
  struct Parameters {
    using ohm = util::quantities::ohm_as<Data_t>;
    using farad = util::quantities::farad_as<Data_t>;
    
    ohm R1;
    ohm R2;
    ohm R3;
    ohm R4;
    farad C1;
    farad C2;
  }; // Parameters
  
  /// Constructor: imports the filter parameters.
  SallenKeyFilter(Parameters parameters);
  
    private:
  using Cache_t = Data_t;
  
  Parameters fParams; ///< Physical components.
  
  Cache_t fK, fA1, fA2; ///< Cached response constants.
  
  /// Response of the circuit to a (complex) frequency `s` [Hz].
  std::complex<Data_t> response(Frequency_t s) const;
  
  virtual WaveformAmplitudes_t doApply
    (Frequencies_t const& s, WaveformAmplitudes_t amplitudes) const
    override;
  
  virtual void doDumpConfig
    (std::ostream& out, details::Indenter& nextLine) const override;
  
}; // SallenKeyFilter


// -----------------------------------------------------------------------------
/// FHiCL configuration for `icarus::trigger::SallenKeyFilter`.
struct icarus::trigger::SallenKeyFilterConfig {
  
  using ohm = util::quantities::ohm;
  using farad = util::quantities::farad;
  
  using Name = fhicl::Name;
  using Comment = fhicl::Comment;
  
  fhicl::Atom<ohm> R1 {
    Name{ "R1" },
    Comment{ "Sallen-Key filter circuit \"R1\" resistance [Ω]" }
    };
  
  fhicl::Atom<ohm> R2 {
    Name{ "R2" },
    Comment{ "Sallen-Key filter circuit \"R2\" resistance [Ω]" }
    };
  
  fhicl::Atom<ohm> R3 {
    Name{ "R3" },
    Comment{ "Sallen-Key filter circuit \"R3\" resistance [Ω]" }
    };
  
  fhicl::Atom<ohm> R4 {
    Name{ "R4" },
    Comment{ "Sallen-Key filter circuit \"R4\" resistance [Ω]" }
    };
  
  fhicl::Atom<farad> C1 {
    Name{ "C1" },
    Comment{ "Sallen-Key filter circuit \"C1\" capacitance [F]" }
    };
  
  fhicl::Atom<farad> C2 {
    Name{ "C2" },
    Comment{ "Sallen-Key filter circuit \"C2\" capacitance [F]" }
    };
  
}; // icarus::trigger::SallenKeyFilterConfig



namespace icarus::trigger {
  
  /// Automatically-invoked conversion of FHiCL configuration
  /// into Sallen-Key algorithm parameters.
  SallenKeyFilter::Parameters convert(SallenKeyFilterConfig const& config);

} // namespace icarus::trigger


// -----------------------------------------------------------------------------

#endif // ICARUSCODE_PMT_TRIGGER_ALGORITHMS_SALLENKEYFILTER_H
