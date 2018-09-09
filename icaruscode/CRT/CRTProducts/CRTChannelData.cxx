#include "icaruscode/CRT/CRTProducts/CRTChannelData.h"

namespace icarus{
 namespace crt{

    CRTChannelData::CRTChannelData() : fChannel(0), fT0(0), fT1(0), fAdc(0), fTrackID(0) {}
    CRTChannelData::CRTChannelData(uint32_t chan, int time0, int time1, uint32_t q, int trackid):
      fChannel(chan),
      fT0(time0),
      fT1(time1),
      fAdc(q),
      fTrackID(trackid) {}
    CRTChannelData::~CRTChannelData() {}
    uint32_t CRTChannelData::Channel() const { 
      return fChannel; 
    }
    int CRTChannelData::T0() const {
      return fT0;
    }
    int CRTChannelData::T1() const { 
      return fT1;
    }
    uint32_t CRTChannelData::ADC() const {
      return fAdc;
    }
    int CRTChannelData::TrackID() const {
      return fTrackID;
    }
    void CRTChannelData::SetADC(uint32_t adc)  {
        this->fAdc = adc;
    }
 }
}

