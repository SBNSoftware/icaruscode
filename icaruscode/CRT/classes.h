#include "canvas/Persistency/Common/Wrapper.h"
#include "icaruscode/CRT/CRTData.hh"
#include <vector>

namespace {
  struct dictionary {
    art::Wrapper<icarus::crt::CRTData> w;
    std::vector<icarus::crt::CRTData> v;
    art::Wrapper<std::vector<icarus::crt::CRTData> > wv;
  };
}


