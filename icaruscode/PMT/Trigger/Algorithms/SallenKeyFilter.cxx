/**
 * @file   icaruscode/PMT/Trigger/Algorithms/SallenKeyFilter.cxx
 * @brief  Sallen-Key filter circuit response.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   August 25, 2025
 * @see    icaruscode/PMT/Trigger/Algorithms/SallenKeyFilter.h
 * 
 */

#undef NDEBUG // FIXME

// library header
#include "icaruscode/PMT/Trigger/Algorithms/SallenKeyFilter.h"

// ICARUS and framework libraries
#include "lardataalg/Utilities/quantities/spacetime.h" // second_as
#include "larcorealg/CoreUtils/enumerate.h"
#include "larcorealg/CoreUtils/counter.h"
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h" // util::pi()

// C/C++ standard libraries
#include <cstdint> // std::size_t
#include <ostream>
#include <utility> // std::move()


// -----------------------------------------------------------------------------
namespace {
  
  /// Returns the value of `t` in seconds.
  template <typename S, typename T>
  T toSecond(util::quantities::scaled_second<S, T> t)
    { return t.template convertInto<util::quantities::second_as<T>>().value(); }
  
} // local namespace


// -----------------------------------------------------------------------------
auto icarus::trigger::convert(SallenKeyFilterConfig const& config)
  -> SallenKeyFilter::Parameters
{
  return {
    /* .R1 = */ config.R1().asValueType<float>(),
    /* .R2 = */ config.R2().asValueType<float>(),
    /* .R3 = */ config.R3().asValueType<float>(),
    /* .R4 = */ config.R4().asValueType<float>(),
    /* .C1 = */ config.C1().asValueType<float>(),
    /* .C2 = */ config.C2().asValueType<float>()
    };
} // convert(SallenKeyFilterConfig)


// -----------------------------------------------------------------------------
// --- icarus::trigger::SallenKeyFilter
// -----------------------------------------------------------------------------
icarus::trigger::SallenKeyFilter::SallenKeyFilter
  (Parameters parameters)
  : fParams{ std::move(parameters) }
  , fK { 1.0f + fParams.R4 / fParams.R3 }
  , fA1{ toSecond(
      (fParams.R1 + fParams.R2) * fParams.C1
    + (1.0f - fK) * fParams.R1 * fParams.C2)
    }
  , fA2{ toSecond(fParams.R1*fParams.C1) * toSecond(fParams.R2*fParams.C2) }
  {}


// -----------------------------------------------------------------------------
auto icarus::trigger::SallenKeyFilter::response
  (Frequency_t s) const -> std::complex<Data_t>
{
  // s in hertz!
  return fK / ((fA2 * s + fA1) * s + static_cast<Data_t>(1.0));
  
} // icarus::trigger::SallenKeyFilter::response()


// -----------------------------------------------------------------------------
auto icarus::trigger::SallenKeyFilter::doApply
  (Frequencies_t const& s, WaveformAmplitudes_t amplitudes) const
  -> WaveformAmplitudes_t
{
  for (std::size_t const k: util::counter(size(s)))
    amplitudes[k] *= response(s[k]);
  
  return amplitudes;
} // icarus::trigger::SallenKeyFilter::doApply()


// -----------------------------------------------------------------------------
void icarus::trigger::SallenKeyFilter::doDumpConfig
  (std::ostream& out, details::Indenter& nextLine) const
{
  out
    << nextLine << "Sallen-Key filter component parameters:"
    << nextLine << "  R1 = " << fParams.R1
    << nextLine << "  R2 = " << fParams.R2
    << nextLine << "  R3 = " << fParams.R3
    << nextLine << "  R4 = " << fParams.R4
    << nextLine << "  C1 = " << fParams.C1
    << nextLine << "  C2 = " << fParams.C2
    << nextLine << "Constant values ( V(out)/V(in) = k / (1 + A1 s + A2 s^2) ):"
    << nextLine << "   k = " << fK
    << nextLine << "  A1 = " << fA1 << " s"
    << nextLine << "  A2 = " << fA2 << " s²"
    ;
} // icarus::trigger::SallenKeyFilter::doDumpConfig()


//------------------------------------------------------------------------------
