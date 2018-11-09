/**
 * \class CRTHit
 *
 * \ingroup crt
 *
 * \brief CRT Hit Info
 *
 * \original author $Author: David Lorca $
 * \ported to icaruscode and modified by Chris.Hilgenberg@colostate.edu
 *
 */
#ifndef CRTHit_hh_
#define CRTHit_hh_
#include <cstdint>
#include <utility> //needed for pair
#include <set>

namespace icarus{
namespace crt {

  class CRTHit {

  public:
     CRTHit();
     CRTHit(int event, double x, double y, double z, \
            double xerr, double yerr, double zerr,  \
            double t0, double t0corr, double t1, double t1corr, \
            std::pair<int,int> macpair, int reg,  \
            std::set<int> trackID, int stripID, int modID);
     virtual ~CRTHit();

     std::pair<int, int> MacPair() const;
     int Event() const;
     double X() const;
     double XErr() const;
     double Y() const;
     double YErr() const;
     double Z() const;
     double ZErr() const;
     double T0() const;
     double T0Corr() const;
     double T1() const;
     double T1Corr() const;
     int Region() const;     
     //int TrackID() const;
     std::set<int> TrackID() const;
     int Strip() const;
     int Module() const;

  private:

      int    fEvent;
      double fX;
      double fY;
      double fZ;
      double fXErr;
      double fYErr;
      double fZErr;
      double fT0;
      double fT0Corr;
      double fT1;
      double fT1Corr;
      std::pair<int,int> fMacPair;
      int fReg;
      //int fTrackID;
      std::set<int> fTrackID;
      int fStrip;
      int fModule;

  }; //class

 }   //crt
}  //icarus
#endif
