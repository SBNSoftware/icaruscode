#include "icaruscode/CRT/CRTProducts/CRTHit.h"

namespace icarus {
namespace crt {

  CRTHit::CRTHit():  \
    fX(0), fY(0), fZ(0),    \
    fXErr(0), fYErr(0), fZErr(0),    \
    fT0(0), fT0Corr(0), fT1(0), fT1Corr(0),  \
    fMacPair(std::make_pair(0,0)), fReg(0),  \
    fTrackID(0), fStrip(0), fModule(0) {}

  CRTHit::CRTHit(float x, float y, float z,    \
                float xerr, float yerr, float zerr,   \
                float t0, float t0corr, float t1, float t1corr,   \
                std::pair<uint16_t,uint16_t> macpair, uint32_t reg, \
                int trackID, uint32_t strip, uint32_t mod):   \
    fX(x), fY(y), fZ(z),   \
    fXErr(xerr), fYErr(yerr), fZErr(zerr),   \
    fT0(t0), fT0Corr(t0corr), fT1(t1), fT1Corr(t1corr),   \
    fMacPair(macpair), fReg(reg),  \
    fTrackID(trackID), fStrip(strip), fModule(mod) {}

  CRTHit::~CRTHit() {}

  std::pair<uint16_t,uint16_t> CRTHit::MacPair() const {
    return fMacPair;
  }
  float CRTHit::X() const {
    return fX;
  }
  float CRTHit::Y() const {
    return fY;
  }
  float CRTHit::Z() const {
    return fZ;
  }
  float CRTHit::XErr() const {
    return fXErr;
  }
  float CRTHit::YErr() const {
    return fYErr;
  }
  float CRTHit::ZErr() const {
    return fZErr;
  }
  float CRTHit::T0() const {
    return fT0;
  }
  float CRTHit::T0Corr() const {
    return fT0Corr;
  }
  float CRTHit::T1() const {
    return fT1;
  }
  float CRTHit::T1Corr() const {
    return fT1Corr;
  }

  uint32_t CRTHit::Region() const {
    return fReg;
  }

  int CRTHit::TrackID() const {
    return fTrackID;
  }

  uint32_t CRTHit::Strip() const {
    return fStrip;
  }

  uint32_t CRTHit::Module() const {
    return fModule;
  }
}
}
