#ifndef CRTChannelData_h_
#define CRTChannelData_h_

#include <stdint.h>
#include <vector>

namespace icarus {
namespace crt {

  class CRTChannelData {

    public:
      CRTChannelData();
      CRTChannelData(int chan, double time0, double time1, int q, std::vector<int> trackid);
      virtual ~CRTChannelData();
      int Channel() const;
      double T0() const;
      double T1() const;
      int ADC() const;
      std::vector<int> TrackID() const;
      void SetADC(int adc);
      void SetTrackID(std::vector<int> vec);

    private:
      int fChannel;
      double fT0;
      double fT1;
      int fAdc;
      std::vector<int> fTrackID;
  };

 }
}

#endif
