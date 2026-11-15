#ifndef MCSFitResultGS_h
#define MCSFitResultGS_h

#include <vector>

namespace recob {
  /**
   * @file    icaruscode/TPC/Tracking/MCS/MCSFitResultGS.h
   * @class   recob::MCSFitResultGS
   * @brief   Class storing the result of Multiple Coulomb Scattering fit. It stores: 
   * the best value and error of C2-fit momentum [GeV/c]
   * the test momentum and the C2 function values
   * the best value and error of likelihood-fit momentum [GeV/c]
   * the test momentum and the likelihood values
   * the sigma3p [cm]
   * the track length on 1D drift direction [cm]
   * the track length on 2D wire plane [cm]
   * the track length on 3D space [cm]
   * the vector of segment number of hits
   * the vector of cumulative segment number of hits
   * the vector of segment lengths [cm]
   * the vector of cumulative segment lengths [cm]
   * the vector of segment sides
   * the vector of linear-fit measured angles [rad]
   * the vector of polygonal measured angles [rad]
   * the validity check if MCS fit is successful
   * @author  G. Chiello (Pisa, ICARUS) and F. Varanini (Padova, ICARUS)
   * @date    2026
   * @version 1.0
   */
  class MCSFitResultGS {
  public:
    /// default empty class
    MCSFitResultGS()
      : bestpC2_(0), errorpC2_(0), testpC2_(), valuesC2_(),
        bestpLL_(0), errorpLL_(0), testpLL_(), valuesLL_(),
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
      double bestpC2, double errorpC2, std::vector<double> testpC2, std::vector<double> valuesC2,
      double bestpLL, double errorpLL, std::vector<double> testpLL, std::vector<double> valuesLL,
      int nhits, double sigma3p, std::vector<double> delta3p, 
      double L1D, double L2D, double L3D,
      std::vector<double> seglens, std::vector<double> cumseglens, 
      std::vector<int> seghits, std::vector<int> cumseghits,
      std::vector<int> segsides,
      int nsegused, int nsegtot,
      std::vector<double> dthetalin, std::vector<double> dthetapoly,
      bool check) 
      : bestpC2_(bestpC2), errorpC2_(errorpC2), testpC2_(testpC2), valuesC2_(valuesC2), 
        bestpLL_(bestpLL), errorpLL_(errorpLL), testpLL_(testpLL), valuesLL_(valuesLL), 
        nhits_(nhits), sigma3p_(sigma3p), delta3p_(delta3p), 
        L1D_(L1D), L2D_(L2D), L3D_(L3D),
        seglens_(seglens), cumseglens_(cumseglens),
        seghits_(seghits), cumseghits_(cumseghits),
        segsides_(segsides),
        nsegused_(nsegused), nsegtot_(nsegtot),
        dthetalin_(dthetalin), dthetapoly_(dthetapoly),
        check_(check)
    {}

    /// momentum best value from the fit [GeV/c]
    double bestMomentumC2() const { 
      return bestpC2_; }

    /// error on momentum from the fit [GeV/c]
    double errorMomentumC2() const { 
      return errorpC2_; }

    /// vector of test momenta used in the fit [GeV/c]
    std::vector<double> testMomentumC2() const { 
      return testpC2_; }

    /// vector of C2 values computed in the fit
    std::vector<double> valuesC2() const { 
      return valuesC2_; }

    /// momentum best value from the likelihood fit [GeV/c]
    double bestMomentumLikelihood() const { 
      return bestpLL_; }

    /// error on momentum from the likelihood fit [GeV/c]
    double errorMomentumLikelihood() const { 
      return errorpLL_; }

    /// vector of test momenta used in the fit [GeV/c]
    std::vector<double> testMomentumLikelihood() const { 
      return testpLL_; }

    /// vector of C2 values computed in the fit
    std::vector<double> valuesLikelihood() const { 
      return valuesLL_; }

    /// RMS of delta3p distribution [cm]
    double sigma3P() const { 
      return sigma3p_; }

    /// track length [cm] projected in 1 dimension (drift direction)
    double length1D() const { 
      return L1D_; }

    /// track length [cm] projected in 2 dimension (wire plane frame)
    double length2D() const { 
      return L2D_; }

    /// track length [cm] projected in 3 dimension (xyz reference frame)
    double length3D() const { 
      return L3D_; }

    /// vector of lengths of the segments [cm]
    std::vector<double> segmentLengths() const { 
      return seglens_; }

    /// vector of cumulative lengths of the segments [cm]
    std::vector<double> segmentCumLengths() const { 
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
    std::vector<double> measuredLinAngles() const { 
      return dthetalin_; }

    /// vector of measured polygonal angles between the segments used in the fit
    std::vector<double> measuredPolyAngles() const { 
      return dthetapoly_; }

    /// validity check if fitter runned
    bool ValidityCheck() const { 
      return check_; }

  private:
    double bestpC2_; double errorpC2_; std::vector<double> testpC2_; std::vector<double> valuesC2_;
    double bestpLL_; double errorpLL_; std::vector<double> testpLL_; std::vector<double> valuesLL_;
    int nhits_; double sigma3p_; std::vector<double> delta3p_;
    double L1D_; double L2D_; double L3D_;
    std::vector<double> seglens_; std::vector<double> cumseglens_;
    std::vector<int> seghits_; std::vector<int> cumseghits_;
    std::vector<int> segsides_;
    int nsegused_; int nsegtot_;
    std::vector<double> dthetalin_; std::vector<double> dthetapoly_;
    bool check_; }; }

#endif