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
   * the C2 function value at range momentum 
   * the number of tails at range momentum
   * the sigma3p [cm]
   * the track length on 1D drift direction [cm]
   * the track length on 2D wire plane [cm]
   * the track length on 3D space [cm]
   * the vector of segment number of hits
   * the vector of cumulative segment number of hits
   * the vector of segment lengths [cm]
   * the vector of cumulative segment lengths [cm]
   * the vector of linear-fit measured angles [rad]
   * the vector of linear-fit expected angles [rad]
   * the vector of polygonal measured angles [rad]
   * the vector of polygonal expected angles [rad]
   * the vector of corresponding-angle momenta [GeV/c]
   * the geometrical check if track is stopping
   * the vector of delta rays flag
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
      : pid_(0),
        bestp_(0), errorp_(0), minp_(0), maxp_(0),
        alpha_(0), dalpha_(0), beta_(0), dbeta_(0),
        testp_(), c2function_(), tailssize_(),
        bestpL_(0), errorpL_(0),
        c2atrange_(0), tailsatrange_(0),
        nhits_(0), sigma3p_(0), delta3p_(0), 
        //sigma3p_1skip_(0), sigma3p_2skip_(0), sigma3p_3skip_(0), sigma3p_4skip_(0), 
        //delta3p_1skip_(), delta3p_2skip_(), delta3p_3skip_(), delta3p_4skip_(), 
        L1D_(0), L2D_(0), L3D_(0),
        seglens_(), cumseglens_(),
        seghits_(), cumseghits_(),
        nsegused_(), nsegtot_(),
        dthetalinexp_(), dthetalin_(),
        dthetapolyexp_(), dthetapoly_(),
        seglens2D_(), dthetamom_(),
        stop_(false), isdelta_(), isdeltaindex_(),
        check_(false)
    {}

    MCSFitResultGS(bool check)
      : MCSFitResultGS() {
        check_ = check; }

    /// default general class
    MCSFitResultGS(
      int pid,
      float bestp, float errorp, float minp, float maxp,
      float alpha, float dalpha, float beta, float dbeta,
      std::vector<float> testp, std::vector<float> c2function,
      std::vector<int> tailssize,
      float bestpL, float errorpL,
      float c2atrange, int tailsatrange, 
      int nhits, double sigma3p, std::vector<double> delta3p, 
      //double sigma3p_1skip, double sigma3p_2skip, double sigma3p_3skip, double sigma3p_4skip, 
      //std::vector<double> delta3p_1skip, std::vector<double> delta3p_2skip, std::vector<double> delta3p_3skip, std::vector<double> delta3p_4skip, 
      float L1D, float L2D, float L3D,
      std::vector<float> seglens, std::vector<float> cumseglens, 
      std::vector<int> seghits, std::vector<int> cumseghits,
      int nsegused, int nsegtot,
      std::vector<float> dthetalinexp, std::vector<float> dthetalin,
      std::vector<float> dthetapolyexp, std::vector<float> dthetapoly,
      std::vector<float> seglens2D, std::vector<float> dthetamom,
      bool stop, std::vector<int> isdelta, std::vector<int> isdeltaindex,
      bool check) 
      : pid_(pid), 
        bestp_(bestp), errorp_(errorp), minp_(minp), maxp_(maxp), 
        alpha_(alpha), dalpha_(dalpha), beta_(beta), dbeta_(dbeta), 
        testp_(testp), c2function_(c2function), tailssize_(tailssize),
        bestpL_(bestpL), errorpL_(errorpL),
        c2atrange_(c2atrange), tailsatrange_(tailsatrange),
        nhits_(nhits), sigma3p_(sigma3p), delta3p_(delta3p), 
        //sigma3p_1skip_(sigma3p_1skip), sigma3p_2skip_(sigma3p_2skip), sigma3p_3skip_(sigma3p_3skip), sigma3p_4skip_(sigma3p_4skip), 
        //delta3p_1skip_(delta3p_1skip), delta3p_2skip_(delta3p_2skip), delta3p_3skip_(delta3p_3skip), delta3p_4skip_(delta3p_4skip), 
        L1D_(L1D), L2D_(L2D), L3D_(L3D),
        seglens_(seglens), cumseglens_(cumseglens),
        seghits_(seghits), cumseghits_(cumseghits),
        nsegused_(nsegused), nsegtot_(nsegtot),
        dthetalinexp_(dthetalinexp), dthetalin_(dthetalin),
        dthetapolyexp_(dthetapolyexp), dthetapoly_(dthetapoly),
        seglens2D_(seglens2D), dthetamom_(dthetamom),
        stop_(stop), isdelta_(isdelta), isdeltaindex_(isdeltaindex),
        check_(check)
    {}

    /// particle id hypothesis used in the fit
    int particleIdHyp() const { 
      return pid_; }

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

    /// vector of test momenta [GeV/c]
    std::vector<float> testMomentum() const { 
      return testp_; }

    /// vector of c2 function values
    std::vector<float> C2Function() const { 
      return c2function_; }

    /// vector of number of tails
    std::vector<int> TailsSize() const { 
      return tailssize_; }

    /// c2 function value at range momentum
    float C2AtRange() const { 
      return c2atrange_; }

    /// number of tails at range momentum
    int TailsAtRange() const { 
      return tailsatrange_; }

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

    /// vector of 2D lengths of the segments [cm]
    std::vector<float> segmentLengths2D() const { 
      return seglens2D_; }

    /// vector of cumulative lengths of the segments [cm]
    std::vector<float> segmentCumLengths() const { 
      return cumseglens_; }

    /// vector of number of hits in the segments 
    std::vector<int> segmentHits() const { 
      return seghits_; }

    /// vector of cumulative number of hits in the segments 
    std::vector<int> segmentCumHits() const { 
      return cumseghits_; }

    /// number of used segments 
    int nSegmentsUsed() const {
      return nsegused_; }

    /// number of total segments 
    int nSegmentsTotal() const {
      return nsegtot_; }

    /// vector of expected linear angles between the segments used in the fit
    std::vector<float> expectedLinAngles() const { 
      return dthetalinexp_; }

    /// vector of measured linear angles between the segments used in the fit
    std::vector<float> measuredLinAngles() const { 
      return dthetalin_; }

    /// vector of measured linear angles between the segments used in the fit
    std::vector<float> momentumLinAngles() const { 
      return dthetamom_; }

    /// vector of expected polygonal angles between the segments used in the fit
    std::vector<float> expectedPolyAngles() const { 
      return dthetapolyexp_; }

    /// vector of measured polygonal angles between the segments used in the fit
    std::vector<float> measuredPolyAngles() const { 
      return dthetapoly_; }

    /// geometrical check if track is stopping inside the detector
    bool GeoStopCheck() const { 
      return stop_; }

    /// validity check if fitter runned
    bool ValidityCheck() const { 
      return check_; }

    /// vector of check on delta rays
    std::vector<int> IsDelta() const {
      return isdelta_; }

    std::vector<int> IsDeltaIndex() const {
      return isdeltaindex_; }

  private:
    int pid_;
    float bestp_; float errorp_; float minp_; float maxp_;
    float alpha_; float dalpha_; float beta_; float dbeta_;
    std::vector<float> testp_; std::vector<float> c2function_; std::vector<int> tailssize_; 
    float bestpL_; float errorpL_;
    float c2atrange_; int tailsatrange_;
    int nhits_; double sigma3p_; std::vector<double> delta3p_;
    //double sigma3p_1skip_; double sigma3p_2skip_; double sigma3p_3skip_; double sigma3p_4skip_; 
    //std::vector<double> delta3p_1skip_; std::vector<double> delta3p_2skip_; std::vector<double> delta3p_3skip_; std::vector<double> delta3p_4skip_; 
    float L1D_; float L2D_; float L3D_;
    std::vector<float> seglens_; std::vector<float> cumseglens_;
    std::vector<int> seghits_; std::vector<int> cumseghits_;
    int nsegused_; int nsegtot_;
    std::vector<float> dthetalinexp_; std::vector<float> dthetalin_; 
    std::vector<float> dthetapolyexp_; std::vector<float> dthetapoly_;
    std::vector<float> seglens2D_; std::vector<float> dthetamom_; 
    bool stop_; std::vector<int> isdelta_; std::vector<int> isdeltaindex_; 
    bool check_; }; }

#endif