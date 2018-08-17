#ifndef CRTChannelData_h_
#define CRTChannelData_h_

#include <stdint.h>

namespace icarus {
namespace crt {

  class CRTChannelData {

    public:
      CRTChannelData();
      CRTChannelData(uint32_t chan, uint32_t time0, uint32_t time1, uint32_t q, uint32_t trackid);
      virtual ~CRTChannelData();
      uint32_t Channel() const;
      uint32_t T0() const;
      uint32_t T1() const;
      uint32_t ADC() const;
      uint32_t TrackID() const;
      void SetADC(uint32_t adc);

    private:
      uint32_t fChannel;
      uint32_t fT0;
      uint32_t fT1;
      uint32_t fAdc;
      uint32_t fTrackID;
  };

 }
}

#endif
