#include "icaruscode/CRT/CRTProducts/CRTData.hh"

namespace icarus{
namespace crt{

  CRTData::CRTData(): fEvent(0), fMac5(0), fEntry(0), fTTrig(0.0), fChanTrig(0), fTrigPair(std::pair<int,int>()), 
          fMacPair(std::pair<int,int>()), fFEBData(std::vector<CRTChannelData>()){
  }
  CRTData::CRTData(int event, int mac5, int entry, double ttrig, int chantrig,
      std::pair<int,int> trigpair, std::pair<int,int> macpair, std::vector<CRTChannelData> febdata):
    fEvent(event),
    fMac5(mac5),
    fEntry(entry),
    fTTrig(ttrig),
    fChanTrig(chantrig),
    fTrigPair(trigpair),
    fMacPair(macpair),
    fFEBData(febdata) {}
  CRTData::~CRTData(){}
  int CRTData::Event() const {
    return fEvent;
  }
  int CRTData::Mac5() const {
    return fMac5;
  }
  int CRTData::Entry() const {
    return fEntry;
  }
  double CRTData::TTrig() const {
    return fTTrig;
  }
  int CRTData::ChanTrig() const {
    return fChanTrig;
  }
  std::pair<int,int> CRTData::TrigPair() const {
    return fTrigPair;
  }
  std::pair<int,int> CRTData::MacPair() const {
    return fMacPair;
  }

  std::vector<CRTChannelData> CRTData::ChanData() const {
    return fFEBData;
  }


} // namespace crt
} // namespace icarus
