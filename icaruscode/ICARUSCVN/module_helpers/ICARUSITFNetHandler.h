/////////////////////////////////////////////////////////////////////////////////////////////////
/// \file    ICARUSTFNetHandler.h
/// \brief   ICARUSTFNetHandler for CVN
/// \author  Varuna Meddage by looking at the format of larrecodnn/CVN/interfaces/ITFNetHandler.h
////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef LCVN_ICARUSTFNETHANDLER_H
#define LCVN_ICARUSTFNETHANDLER_H

#include <memory>
#include <vector>

#include "larrecodnn/CVN/func/InteractionType.h"
#include "larrecodnn/CVN/func/PixelMap.h"

namespace lcvn {

  /// Wrapper for caffe::Net which handles construction and prediction
  class ICARUSITFNetHandler {
  public:
    virtual ~ICARUSITFNetHandler() noexcept = default;
    /// Return prediction arrays for PixelMap
    virtual std::vector<std::vector<float>> Predict(const PixelMap& pm, const int event, const std::string cryo) const = 0;
    virtual std::vector<std::vector<float>> PredictFromArray(const std::vector<unsigned char> &pa, const std::string event) const = 0;
  };

}

#endif // LCVN_ICARUSTFNETHANDLER_H
