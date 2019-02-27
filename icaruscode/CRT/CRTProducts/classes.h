#include "canvas/Persistency/Common/Wrapper.h"
#include "canvas/Persistency/Common/Assns.h"
#include "icaruscode/CRT/CRTProducts/CRTChannelData.h"
#include "icaruscode/CRT/CRTProducts/CRTData.hh"
#include "icaruscode/CRT/CRTProducts/CRTHit.hh"
#include <vector>
#include <map>
#include <utility>
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

    std::vector<pair<int,float> > vpif;
    art::Wrapper<std::vector<pair<int,float> > > wvpif;
    std::map<uint8_t,std::vector<std::pair<int,float> > > muvpif;
    art::Wrapper<std::map<uint8_t,std::vector<std::pair<int,float> > > > wmuvpif;
    art::Wrapper<icarus::crt::CRTHit> wh;
    std::vector<icarus::crt::CRTHit> vh;
    art::Wrapper< std::vector<icarus::crt::CRTHit> > wvh;
    art::Assns<icarus::crt::CRTData,icarus::crt::CRTHit,void> adh;
    art::Assns<icarus::crt::CRTHit,icarus::crt::CRTData,void> ahd;
    art::Wrapper<art::Assns<icarus::crt::CRTData,icarus::crt::CRTHit,void> > wadh;
    art::Wrapper<art::Assns<icarus::crt::CRTHit,icarus::crt::CRTData,void> > wahd;


  };
}
