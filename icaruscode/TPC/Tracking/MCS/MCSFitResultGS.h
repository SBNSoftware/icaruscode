#ifndef MCSFitResultGS_h
#define MCSFitResultGS_h

#include <vector>

namespace recob {
  /**
   * @file  icaruscode/TPC/Tracking/MCS/MCSFitResultGS.h
   * @class recob::MCSFitResultGS
   * @brief Class storing the result of fit of Multiple Coulomb Scattering angles between segments within a Track or Trajectory. It stores: 
   * the PID 
   * the best value and error of C2-fit momentum [GeV/c]
   * the min and max C2-fit momenta hypothesis [GeV/c]
   * the best value and error C2-fit parameters 
   * the best value and error of likelihood-fit momentum [GeV/c]
   * the sigma3p [cm]
   * the track length on 1D drift direction [cm]
   * the track length on 2D wire plane [cm]
   * the track length on 3D space [cm]
   * the vector of segment number of hits
   * the vector of cumulative segment number of hits
   * the vector of segment lengths [cm]
   * the vector of cumulative segment lengths [cm]
   * the vector of linear-fit measured angles [rad]
   * the vector of polygonal measured angles [rad]
   * the validity check if fit is successful
   *
   * @author  G. Chiello (Pisa, ICARUS) and F. Varanini (Padova, ICARUS)
   * @date    2026
   * @version 1.0
   */
  class MCSFitResultGS {
  public:
    /// default empty class
    MCSFitResultGS()
      : bestp_(0), errorp_(0), minp_(0), maxp_(0),
        alpha_(0), dalpha_(0), beta_(0), dbeta_(0),
        bestpL_(0), errorpL_(0),
        nhits_(0), sigma3p_(0), delta3p_(), 
        L1D_(0), L2D_(0), L3D_(0),
        seglens_(), cumseglens_(),
        seghits_(), cumseghits_(),
        segsides_(),
        nsegused_(), nsegtot_(),
        dthetalin_(), dthetapoly_(),
        check_(false)
    {}

    MCSFitResultGS(bool check)
      : MCSFitResultGS() {
        check_ = check; }

    /// default general class
    MCSFitResultGS(
      float bestp, float errorp, float minp, float maxp,
      float alpha, float dalpha, float beta, float dbeta,
      float bestpL, float errorpL,
      int nhits, double sigma3p, std::vector<double> delta3p, 
      float L1D, float L2D, float L3D,
      std::vector<float> seglens, std::vector<float> cumseglens, 
      std::vector<int> seghits, std::vector<int> cumseghits,
      std::vector<int> segsides,
      int nsegused, int nsegtot,
      std::vector<float> dthetalin, std::vector<float> dthetapoly,
      bool check) 
      : bestp_(bestp), errorp_(errorp), minp_(minp), maxp_(maxp), 
        alpha_(alpha), dalpha_(dalpha), beta_(beta), dbeta_(dbeta), 
        bestpL_(bestpL), errorpL_(errorpL),
        nhits_(nhits), sigma3p_(sigma3p), delta3p_(delta3p), 
        L1D_(L1D), L2D_(L2D), L3D_(L3D),
        seglens_(seglens), cumseglens_(cumseglens),
        seghits_(seghits), cumseghits_(cumseghits),
        segsides_(segsides),
        nsegused_(nsegused), nsegtot_(nsegtot),
        dthetalin_(dthetalin), dthetapoly_(dthetapoly),
        check_(check)
    {}

    /// momentum best value from the likelihood fit [GeV/c]
    float bestMomentum_Likelihood() const { 
      return bestpL_; }

    /// error on momentum from the likelihood fit [GeV/c]
    float errorMomentum_Likelihood() const { 
      return errorpL_; }

    /// momentum best value from the fit [GeV/c]
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

    /// RMS of delta3p distribution [cm]
    double sigma3P() const { 
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

    /// vector of TPC side for the segments 
    std::vector<int> segmentSides() const { 
      return segsides_; }

    /// number of used segments 
    int nSegmentsUsed() const {
      return nsegused_; }

    /// number of total segments 
    int nSegmentsTotal() const {
      return nsegtot_; }

    /// vector of measured linear angles between the segments used in the fit
    std::vector<float> measuredLinAngles() const { 
      return dthetalin_; }

    /// vector of measured polygonal angles between the segments used in the fit
    std::vector<float> measuredPolyAngles() const { 
      return dthetapoly_; }

    /// validity check if fitter runned
    bool ValidityCheck() const { 
      return check_; }

  private:
    float bestp_; float errorp_; float minp_; float maxp_;
    float alpha_; float dalpha_; float beta_; float dbeta_;
    float bestpL_; float errorpL_;
    int nhits_; double sigma3p_; std::vector<double> delta3p_;
    float L1D_; float L2D_; float L3D_;
    std::vector<float> seglens_; std::vector<float> cumseglens_;
    std::vector<int> seghits_; std::vector<int> cumseghits_;
    std::vector<int> segsides_;
    int nsegused_; int nsegtot_;
    std::vector<float> dthetalin_; std::vector<float> dthetapoly_;
    bool check_; }; }

#endif