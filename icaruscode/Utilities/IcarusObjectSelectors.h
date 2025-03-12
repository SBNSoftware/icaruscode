/**
 * @file   icaruscode/Utilities/IcarusObjectSelectors.h
 * @brief  Selector implementations for some ICARUS enumerator data types.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   May 11, 2023
 * @see    icaruscode/Utilities/IcarusObjectSelectors.cxx
 */


#ifndef ICARUSCODE_UTILITIES_ICARUSOBJECTSELECTORS_H
#define ICARUSCODE_UTILITIES_ICARUSOBJECTSELECTORS_H


#include "sbnobj/Common/CRT/CRTPMTMatching.hh"
#include "lardataalg/Utilities/MultipleChoiceSelection.h"


namespace icarus::crt {
  
  /// Selector for `icarus::crt::MatchType`.
  extern util::MultipleChoiceSelection<sbn::crt::MatchType> const MatchTypeSelector;
  
} // icarus::crt


#endif // ICARUSCODE_UTILITIES_ICARUSOBJECTSELECTORS_H
