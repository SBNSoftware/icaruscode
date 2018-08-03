#include "icaruscode/CRT/CRTProducts/CRTData.hh"

namespace icarus{
namespace crt{

  CRTData::CRTData(): fMac5(0), fEntry(0), fTTrig(0), fChanTrig(0), fTrigPair(std::pair<uint32_t,uint32_t>()), 
          fMacPair(std::pair<uint32_t,uint32_t>()), fFEBData(std::vector<CRTChannelData>()){
  }
  CRTData::CRTData(uint32_t mac5, uint32_t entry, uint32_t ttrig, uint32_t chantrig,
      std::pair<uint32_t,uint32_t> trigpair, std::pair<uint32_t,uint32_t> macpair, std::vector<CRTChannelData> febdata):
    fMac5(mac5),
    fEntry(entry),
    fTTrig(ttrig),
    fChanTrig(chantrig),
    fTrigPair(trigpair),
    fMacPair(macpair),
    fFEBData(febdata) {}
  CRTData::~CRTData(){}
  uint32_t CRTData::Mac5() const {
    return fMac5;
  }
  uint32_t CRTData::Entry() const {
    return fEntry;
  }
  uint32_t CRTData::TTrig() const {
    return fTTrig;
  }
  uint32_t CRTData::ChanTrig() const {
    return fChanTrig;
  }
  std::pair<uint32_t,uint32_t> CRTData::TrigPair() const {
    return fTrigPair;
  }
  std::pair<uint32_t,uint32_t> CRTData::MacPair() const {
    return fMacPair;
  }

  std::vector<CRTChannelData> CRTData::ChanData() const {
    return fFEBData;
  }


} // namespace crt
} // namespace icarus
