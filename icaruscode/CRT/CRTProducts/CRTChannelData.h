#ifndef CRTChannelData_h_
#define CRTChannelData_h_

#include <stdint.h>

namespace icarus {
namespace crt {

  class CRTChannelData {

    public:
      CRTChannelData();
      CRTChannelData(uint32_t chan, uint32_t time0, uint32_t time1, uint32_t q);
      virtual ~CRTChannelData();
      uint32_t Channel() const;
      uint32_t T0() const;
      uint32_t T1() const;
      uint32_t ADC() const;

    private:
      uint32_t fChannel;
      uint32_t fT0;
      uint32_t fT1;
      uint32_t fAdc;
  };

 }
}

#endif
