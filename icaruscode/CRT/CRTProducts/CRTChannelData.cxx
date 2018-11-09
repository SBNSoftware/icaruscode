#include "icaruscode/CRT/CRTProducts/CRTChannelData.h"

namespace icarus{
 namespace crt{

    CRTChannelData::CRTChannelData() : fChannel(0), fT0(0.0), fT1(0.0), fAdc(0), fTrackID{} {}
    CRTChannelData::CRTChannelData(int chan, double time0, double time1, int q, std::vector<int> trackid):
      fChannel(chan),
      fT0(time0),
      fT1(time1),
      fAdc(q),
      fTrackID(trackid) {}
    CRTChannelData::~CRTChannelData() {}
    int CRTChannelData::Channel() const { 
      return fChannel; 
    }
    double CRTChannelData::T0() const {
      return fT0;
    }
    double CRTChannelData::T1() const { 
      return fT1;
    }
    int CRTChannelData::ADC() const {
      return fAdc;
    }
    std::vector<int> CRTChannelData::TrackID() const {
      return fTrackID;
    }
    void CRTChannelData::SetADC(int adc)  {
        this->fAdc = adc;
    }
    void CRTChannelData::SetTrackID(std::vector<int> vec) {
        this->fTrackID = vec;
    }
 }
}

