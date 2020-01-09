#ifndef ICCRTData_hh_
#define ICCRTData_hh_

#include <stdint.h> //int
#include <map>
#include <utility>
#include <vector>
#include "icaruscode/CRT/CRTProducts/CRTChannelData.h"

namespace icarus {
namespace crt {

  class CRTData {

  public:
    CRTData();
    CRTData(int event, int mac5, int entry, double ttrig, int chantrig,
        std::pair<int,int> tpair, std::pair<int,int> macpair, std::vector<CRTChannelData> febdata);
    virtual ~CRTData();

    int Event() const;
    int Entry() const;
    int Mac5() const;
    double TTrig() const;
    int ChanTrig() const;
    std::pair<int,int> TrigPair() const;
    std::pair<int,int> MacPair() const;
    std::vector<CRTChannelData> ChanData() const;

  private:
    int fEvent;
    int fMac5;
    int fEntry;
    double fTTrig;
    int  fChanTrig;
    std::pair<int,int> fTrigPair;
    std::pair<int,int> fMacPair;
    std::vector<CRTChannelData> fFEBData;
  };

} // namespace crt
} // namespace icarus

#endif
