#ifndef MCSFitResultGS_h
#define MCSFitResultGS_h

#include <vector>

namespace recob {
  /**
   * @file  icaruscode/TPC/Tracking/MCS/MCSFitResultGS.h
   * @class recob::MCSFitResultGS
   * @brief Class storing the result of fit of Multiple Coulomb Scattering angles between segments within a Track or Trajectory. It stores: 
   * the PID 
   * the best value of momentum [GeV/c]
   * the error on momentum [GeV/c]
   * the sigma3p [cm]
   * the track length on 1D drift direction [cm]
   * the track length on 2D wire plane [cm]
   * the track length on 3D space [cm]
   * the vector of segment number of hits
   * the vector of cumulative segment number of hits
   * the vector of segment lengths [cm]
   * the vector of cumulative segment lengths [cm]
   * the vector of linear scattering angles [rad]
   * the vector of polygonal scattering angles [rad]
   * the geometrical check if track is stopping
   *
   * @author  G. Chiello (Pisa, ICARUS) based on code from G. Cerati 
   * @date    2025
   * @version 2.0
   */
  class MCSFitResultGS {
  public:
    MCSFitResultGS() = default;
    MCSFitResultGS(
      int pid,
      float bestp, float errorp, float minp, float maxp,
      float alpha, float dalpha, float beta, float dbeta,
      std::vector<float> testp, std::vector<float> c2function,
      float c2atrange,
      float sigma3p, 
      float L1D, float L2D, float L3D,
      std::vector<float> seglens, std::vector<float> cumseglens, 
      std::vector<int> seghits, std::vector<int> cumseghits,
      std::vector<float> dthetalinexp, std::vector<float> dthetalin,
      std::vector<float> dthetapolyexp, std::vector<float> dthetapoly,
      bool stop, std::vector<int> isdelta) 
      : pid_(pid), 
        bestp_(bestp), errorp_(errorp), minp_(minp), maxp_(maxp), 
        alpha_(alpha), dalpha_(dalpha), beta_(beta), dbeta_(dbeta), 
        testp_(testp), c2function_(c2function), 
        c2atrange_(c2atrange),
        sigma3p_(sigma3p),
        L1D_(L1D), L2D_(L2D), L3D_(L3D),
        seglens_(seglens), cumseglens_(cumseglens),
        seghits_(seghits), cumseghits_(cumseghits),
        dthetalinexp_(dthetalinexp), dthetalin_(dthetalin),
        dthetapolyexp_(dthetapolyexp), dthetapoly_(dthetapoly),
        stop_(stop), isdelta_(isdelta)
    {}

    /// particle id hypothesis used in the fit
    int particleIdHyp() const { 
      return pid_; }

    /// momentum best value from fthe it [GeV/c]
    float bestMomentum() const { 
      return bestp_; }

    /// error on momentum from the fit [GeV/c]
    float errorMomentum() const { 
      return errorp_; }

    /// minimum momentum used for the fit [GeV/c]
    float minMomentum() const { 
      return minp_; }

    /// maximum momentum used for the fit [GeV/c]
    float maxMomentum() const { 
      return maxp_; }

    /// alpha parameter of the fit
    float alphaFit() const { 
      return alpha_; }

    /// error on alpha parameter of the fit
    float alphaFitErr() const { 
      return dalpha_; }

    /// beta parameter of the fit
    float betaFit() const { 
      return beta_; }

    /// error on beta parameter of the fit
    float betaFitErr() const { 
      return dbeta_; }

    /// vector of test momenta [GeV/c]
    std::vector<float> testMomentum() const { 
      return testp_; }

    /// vector of c2 function values
    std::vector<float> C2Function() const { 
      return c2function_; }

    /// c2 function value at range momentum
    float C2AtRange() const { 
      return c2atrange_; }

    /// RMS of delta3p distribution [cm]
    float sigma3P() const { 
      return sigma3p_; }

    /// track length [cm] projected in 1 dimension (drift direction)
    float length1D() const { 
      return L1D_; }

    /// track length [cm] projected in 2 dimension (wire plane frame)
    float length2D() const { 
      return L2D_; }

    /// track length [cm] projected in 3 dimension (xyz reference frame)
    float length3D() const { 
      return L3D_; }

    /// vector of lengths of the segments [cm]
    std::vector<float> segmentLengths() const { 
      return seglens_; }

    /// vector of cumulative lengths of the segments [cm]
    std::vector<float> segmentCumLengths() const { 
      return cumseglens_; }

    /// vector of number of hits in the segments 
    std::vector<int> segmentHits() const { 
      return seghits_; }

    /// vector of cumulative number of hits in the segments 
    std::vector<int> segmentCumHits() const { 
      return cumseghits_; }

    /// vector of expected linear angles between the segments used in the fit
    std::vector<float> expectedLinAngles() const { 
      return dthetalinexp_; }

    /// vector of measured linear angles between the segments used in the fit
    std::vector<float> measuredLinAngles() const { 
      return dthetalin_; }

    /// vector of expected polygonal angles between the segments used in the fit
    std::vector<float> expectedPolyAngles() const { 
      return dthetapolyexp_; }

    /// vector of measured polygonal angles between the segments used in the fit
    std::vector<float> measuredPolyAngles() const { 
      return dthetapoly_; }

    /// geometrical check if track is stopping inside the detector
    bool GeoStopCheck() const { 
      return stop_; }

    /// vector of check on delta rays
    std::vector<int> IsDelta() const {
      return isdelta_; }
  private:
    int pid_;
    float bestp_; float errorp_; float minp_; float maxp_;
    float alpha_; float dalpha_; float beta_; float dbeta_;
    std::vector<float> testp_; std::vector<float> c2function_;
    float c2atrange_;
    float sigma3p_;
    float L1D_; float L2D_; float L3D_;
    std::vector<float> seglens_; std::vector<float> cumseglens_;
    std::vector<int> seghits_; std::vector<int> cumseghits_;
    std::vector<float> dthetalinexp_; std::vector<float> dthetalin_; 
    std::vector<float> dthetapolyexp_; std::vector<float> dthetapoly_;
    bool stop_; std::vector<int> isdelta_;
  };
}

#endif