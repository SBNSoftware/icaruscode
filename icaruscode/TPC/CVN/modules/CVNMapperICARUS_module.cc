#include "icaruscode/TPC/CVN/module_helpers/ICVNMapperICARUS.h"
#include "icaruscode/TPC/CVN/module_helpers/ICVNMapperICARUS.cxx"
#include "larrecodnn/CVN/module_helpers/ICVNMapper.h"
#include "larrecodnn/CVN/module_helpers/ICVNMapper.cxx"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/Simulation/SimChannel.h"

#include "icaruscode/TPC/CVN/module_helpers/PixelMapProducerICARUS.h"

namespace cvn {

  typedef ICVNMapperICARUS<PixelMapHitProducerICARUS, recob::Hit> CVNMapperICARUS;
  template class ICVNMapperICARUS<PixelMapHitProducerICARUS, recob::Hit>;

DEFINE_ART_MODULE(CVNMapperICARUS)
}
