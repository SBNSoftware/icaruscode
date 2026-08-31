#include "icaruscode/ICARUSCVN/module_helpers/ICARUSICVNMapper.cxx"
#include "icaruscode/ICARUSCVN/module_helpers/ICARUSICVNMapper.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/Simulation/SimChannel.h"

#include "icaruscode/ICARUSCVN/module_helpers/ICARUSPixelMapProducer.h"
#include "larrecodnn/CVN/func/PixelMap.h"

namespace lcvn {

  typedef ICARUSICVNMapper<lcvn::ICARUSPixelMapHitProducer, recob::Hit> ICARUSCVNMapper;
  template class ICARUSICVNMapper<lcvn::ICARUSPixelMapHitProducer, recob::Hit>;

DEFINE_ART_MODULE(lcvn::ICARUSCVNMapper)
}
