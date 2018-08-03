#include "icaruscode/CRT/CRTProducts/CRTChannelData.h"

namespace icarus{
 namespace crt{

    CRTChannelData::CRTChannelData() : fChannel(0), fT0(0), fT1(0), fAdc(0) {}
    CRTChannelData::CRTChannelData(uint32_t chan, uint32_t time0, uint32_t time1, uint32_t q):
      fChannel(chan),
      fT0(time0),
      fT1(time1),
      fAdc(q) {}
    CRTChannelData::~CRTChannelData() {}
    uint32_t CRTChannelData::Channel() const { 
      return fChannel; 
    }
    uint32_t CRTChannelData::T0() const {
      return fT0;
    }
    uint32_t CRTChannelData::T1() const { 
      return fT1;
    }
    uint32_t CRTChannelData::ADC() const {
      return fAdc;
    }
 }
}

