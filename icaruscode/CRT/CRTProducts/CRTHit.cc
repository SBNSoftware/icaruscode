#include "icaruscode/CRT/CRTProducts/CRTHit.h"

namespace icarus {
namespace crt {

  CRTHit::CRTHit():  \
    fEvent(0),       \
    fX(0), fY(0), fZ(0),    \
    fXErr(0), fYErr(0), fZErr(0),    \
    fT0(0), fT0Corr(0), fT1(0), fT1Corr(0),  \
    fMacPair(std::make_pair(0,0)), fReg(0),  \
    fTrackID{}, fStrip(0), fModule(0) {}

  CRTHit::CRTHit(int event, double x, double y, double z,    \
                double xerr, double yerr, double zerr,   \
                double t0, double t0corr, double t1, double t1corr,   \
                std::pair<int,int> macpair, int reg, \
                std::set<int> trackID, int strip, int mod):   \
    fEvent(event),         \
    fX(x), fY(y), fZ(z),   \
    fXErr(xerr), fYErr(yerr), fZErr(zerr),   \
    fT0(t0), fT0Corr(t0corr), fT1(t1), fT1Corr(t1corr),   \
    fMacPair(macpair), fReg(reg),  \
    fTrackID(trackID), fStrip(strip), fModule(mod) {}

  CRTHit::~CRTHit() {}

  std::pair<int,int> CRTHit::MacPair() const {
    return fMacPair;
  }
  int    CRTHit::Event() const {
    return fEvent;
  }
  double CRTHit::X() const {
    return fX;
  }
  double CRTHit::Y() const {
    return fY;
  }
  double CRTHit::Z() const {
    return fZ;
  }
  double CRTHit::XErr() const {
    return fXErr;
  }
  double CRTHit::YErr() const {
    return fYErr;
  }
  double CRTHit::ZErr() const {
    return fZErr;
  }
  double CRTHit::T0() const {
    return fT0;
  }
  double CRTHit::T0Corr() const {
    return fT0Corr;
  }
  double CRTHit::T1() const {
    return fT1;
  }
  double CRTHit::T1Corr() const {
    return fT1Corr;
  }

  int CRTHit::Region() const {
    return fReg;
  }

  std::set<int> CRTHit::TrackID() const {
    return fTrackID;
  }

  int CRTHit::Strip() const {
    return fStrip;
  }

  int CRTHit::Module() const {
    return fModule;
  }
}
}
