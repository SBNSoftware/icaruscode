/**
 * @file   icaruscode/Utilities/IcarusObjectSelectors.cxx
 * @brief  Selector implementations for some ICARUS enumerator data types.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   May 11, 2023
 * @see    icaruscode/Utilities/IcarusObjectSelectors.h
 */

#include "icaruscode/Utilities/IcarusObjectSelectors.h"

#include <string>


// ---  BEGIN  --- icarus::crt::MatchType --------------------------------------
using icarus::crt::MatchType;
using MatchTypeSelOption_t = util::MultipleChoiceSelection<MatchType>::Option_t;

util::MultipleChoiceSelection<MatchType> const icarus::crt::MatchTypeSelector{
    MatchTypeSelOption_t{ MatchType::noMatch          , "noMatch"           }
  , MatchTypeSelOption_t{ MatchType::enTop            , "enTop"             }
  , MatchTypeSelOption_t{ MatchType::enSide           , "enSide"            }
  , MatchTypeSelOption_t{ MatchType::enTop_exSide     , "enTop_exSide"      }
  , MatchTypeSelOption_t{ MatchType::exTop            , "exTop"             }
  , MatchTypeSelOption_t{ MatchType::exSide           , "exSide"            }
  , MatchTypeSelOption_t{ MatchType::enTop_mult       , "enTop_mult"        }
  , MatchTypeSelOption_t{ MatchType::enTop_exSide_mult, "enTop_exSide_mult" }
  , MatchTypeSelOption_t{ MatchType::others           , "others"            }
};

// ---  END  ----- icarus::crt::MatchType --------------------------------------
