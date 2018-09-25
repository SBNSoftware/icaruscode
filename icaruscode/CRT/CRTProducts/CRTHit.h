/**
 * \class CRTData
 *
 * \ingroup crt
 *
 * \brief CRT Hit Info
 *
 * \author $Author: David Lorca $
 *
 */
#ifndef CRTHit_hh_
#define CRTHit_hh_
#include <cstdint>
#include <utility> //needed for pair

namespace icarus{
namespace crt {

  class CRTHit {

  public:
     CRTHit();
     CRTHit(float x, float y, float z, \
            float xerr, float yerr, float zerr,  \
            float t0, float t0corr, float t1, float t1corr, \
            std::pair<uint16_t,uint16_t> macpair, uint32_t reg,  \
            int trackID, uint32_t stripID, uint32_t modID);
     virtual ~CRTHit();

     std::pair<uint16_t, uint16_t> MacPair() const;
     float X() const;
     float XErr() const;
     float Y() const;
     float YErr() const;
     float Z() const;
     float ZErr() const;
     float T0() const;
     float T0Corr() const;
     float T1() const;
     float T1Corr() const;
     uint32_t Region() const;     
     int TrackID() const;
     uint32_t Strip() const;
     uint32_t Module() const;

  private:

      float fX;
      float fY;
      float fZ;
      float fXErr;
      float fYErr;
      float fZErr;
      float fT0;
      float fT0Corr;
      float fT1;
      float fT1Corr;
      std::pair<uint16_t,uint16_t> fMacPair;
      uint32_t fReg;
      int fTrackID;
      uint32_t fStrip;
      uint32_t fModule;

  }; //class

 }   //crt
}  //icarus
#endif
