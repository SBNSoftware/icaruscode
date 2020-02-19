/**
 * @file   icaruscode/PMT/Trigger/Utilities/ROOTutils.tcc
 * @brief  A bunch of diverse utilities and futilities related to ROOT.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   April 1, 2019
 * @see    `icaruscode/PMT/Trigger/Utilities/ROOTutils.h`
 */

#ifndef ICARUSCODE_PMT_TRIGGER_UTILITIES_ROOTUTILS_TCC
#define ICARUSCODE_PMT_TRIGGER_UTILITIES_ROOTUTILS_TCC


#ifndef ICARUSCODE_PMT_TRIGGER_UTILITIES_ROOTUTILS_H
# error("'ROOTutils.tcc' is not meant to be included directly." \
  "Please #include \"icaruscode/PMT/Trigger/Utilities/ROOTutils.h\" instead.")
#endif // ICARUSCODE_PMT_TRIGGER_UTILITIES_ROOTUTILS_H



// LArSoft libraries
#include "larcorealg/CoreUtils/enumerate.h"
#include "larcorealg/CoreUtils/StdUtils.h" // util::to_string()

// C/C++ libraries
#include <iterator> // std::next()
#include <vector>
#include <string>
#include <cassert>


// -----------------------------------------------------------------------------
template <typename Coll>
std::pair<std::vector<double>, std::vector<std::string>>
util::ROOT::makeVariableBinningAndLabels(Coll const& centralPoints) {
  assert(std::size(centralPoints) > 1U);
  
  std::pair<std::vector<double>, std::vector<std::string>> result;
  auto& [ binning, labels ] = result;
  
  auto const nPoints = util::size(centralPoints);
  binning.reserve(nPoints + 1U);
  labels.reserve(nPoints);
  
  auto iCenter = util::begin(centralPoints);
  auto const pend = util::end(centralPoints);
  
  auto const startHalfWidth =
    (static_cast<double>(*std::next(iCenter)) - static_cast<double>(*iCenter))
    / 2
    ;
  // lower boundary of current bin:
  double lower = static_cast<double>(*iCenter) - startHalfWidth;
  
  do {
    double const center = static_cast<double>(*iCenter);
    
    if (center <= lower) {
      // sorry, no solution to this except to use a better algorithm
      cet::exception e("makeVariableBinningAndLabels");
      e << "Logic error:"
        " the algorithm does not cope well with the specified set of points:";
      for (auto const& point: centralPoints)
        e << " " << util::to_string(point);
      e << " (trouble at " << util::to_string(*iCenter) << ").\n";
      throw e;
    } // if algorithm failure
    
    double const halfWidth = center - lower;
    
    labels.push_back(util::to_string(*iCenter));
    
    binning.push_back(lower);
    lower += halfWidth * 2;
    
  } while (++iCenter != pend);
  
  binning.push_back(lower); // at this point, it's "upper"
  return result;
} // util::ROOT::makeVariableBinningAndLabels()


// -----------------------------------------------------------------------------


#endif // ICARUSCODE_PMT_TRIGGER_UTILITIES_ROOTUTILS_TCC
