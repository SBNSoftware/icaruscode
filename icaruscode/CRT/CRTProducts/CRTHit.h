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
            std::pair<uint16_t,uint16_t> macpair);
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

  }; //class

 }   //crt
}  //icarus
#endif
