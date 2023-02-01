#include "icaruscode/TPC/CVN/interfaces/ICVNMapperICARUS.h"
#include "icaruscode/TPC/CVN/interfaces/ICVNMapperICARUS.cxx"
#include "icaruscode/TPC/CVN/interfaces/PixelMapProducerICARUS.h"
#include "larrecodnn/CVN/interfaces/ICVNMapper.h"
#include "larrecodnn/CVN/interfaces/ICVNMapper.cxx"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/Simulation/SimChannel.h"

namespace lcvn {

  typedef ICVNMapperICARUS<PixelMapHitProducerICARUS, recob::Hit> CVNMapperICARUS;
  template class ICVNMapperICARUS<PixelMapHitProducerICARUS, recob::Hit>;

}

DEFINE_ART_MODULE(lcvn::CVNMapperICARUS)