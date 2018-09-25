#ifndef CRTData_hh_
#define CRTData_hh_

#include <stdint.h> //uint32_t
#include <map>
#include <utility>
#include <vector>
#include "icaruscode/CRT/CRTProducts/CRTChannelData.h"

namespace icarus {
namespace crt {

  class CRTData {

  public:
    CRTData();
    CRTData(uint32_t event, uint32_t mac5, uint32_t entry, uint32_t ttrig, uint32_t chantrig,
        std::pair<uint32_t,uint32_t> tpair, std::pair<uint32_t,uint32_t> macpair, std::vector<CRTChannelData> febdata);
    virtual ~CRTData();

    uint32_t Event() const;
    uint32_t Entry() const;
    uint32_t Mac5() const;
    uint32_t TTrig() const;
    uint32_t ChanTrig() const;
    std::pair<uint32_t,uint32_t> TrigPair() const;
    std::pair<uint32_t,uint32_t> MacPair() const;
    std::vector<CRTChannelData> ChanData() const;

  private:
    uint32_t fEvent;
    uint32_t fMac5;
    uint32_t fEntry;
    uint32_t fTTrig;
    uint32_t  fChanTrig;
    std::pair<uint32_t,uint32_t> fTrigPair;
    std::pair<uint32_t,uint32_t> fMacPair;
    std::vector<CRTChannelData> fFEBData;
  };

} // namespace crt
} // namespace icarus

#endif
