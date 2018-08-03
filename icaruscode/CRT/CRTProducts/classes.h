#include "canvas/Persistency/Common/Wrapper.h"
#include "icaruscode/CRT/CRTProducts/CRTChannelData.h"
#include "icaruscode/CRT/CRTProducts/CRTData.hh"
#include "icaruscode/CRT/CRTProducts/CRTHit.h"
#include <vector>
//#include <utility>
//#include <cstdint>

namespace {
  struct dictionary {
    art::Wrapper<icarus::crt::CRTChannelData> wc;
    std::vector<icarus::crt::CRTChannelData> vc;
    art::Wrapper<std::vector<icarus::crt::CRTChannelData>> wvc;

    //std::pair<unsigned short,unsigned short> p8;
    //art::Wrapper<std::pair<unsigned short,unsigned short>> wp8;

    art::Wrapper<icarus::crt::CRTData> w;
    std::vector<icarus::crt::CRTData> v;
    art::Wrapper<std::vector<icarus::crt::CRTData> > wv;

    art::Wrapper<icarus::crt::CRTHit> wh;
    std::vector<icarus::crt::CRTHit> vh;
    art::Wrapper< std::vector<icarus::crt::CRTHit> > wvh;

  };
}
