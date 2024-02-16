////////////////////////////////////////////////////////////////////////
// \file DetectorPropertiesICARUSClockOffsetMC.h
//
// \brief Goes with DetectorPropertiesICARUSClockOffsetMC.cxx ()
//
// \author brebel@fnal.gov
//
// Separation of service from Detector info class:
// jpaley@fnal.gov
//
// Update for ICARUS by Bruce Howard and Jaesung Kim
////////////////////////////////////////////////////////////////////////
#ifndef DETINFO_DETECTORPROPERTIESICARUS_H
#define DETINFO_DETECTORPROPERTIESICARUS_H

// LArSoft libraries
#include "larcorealg/CoreUtils/ProviderPack.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"
#include "lardataalg/DetectorInfo/DetectorClocks.h"
#include "lardataalg/DetectorInfo/DetectorProperties.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesData.h"
#include "lardataalg/DetectorInfo/LArProperties.h"

// framework libraries
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/fwd.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/Sequence.h"

// C/C++ standard libraries
#include <set>

/// General LArSoft Utilities
namespace detinfo {

  class DetectorPropertiesICARUSClockOffsetMC final : public DetectorProperties {
  public:
    /// List of service providers we depend on
    using providers_type = lar::ProviderPack<geo::GeometryCore, detinfo::LArProperties>;

    /// Structure for configuration parameters
    struct Configuration_t {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Sequence<double> Efield{
        Name("Efield"),
        Comment("electric field in front of each wire plane (the last one is "
                "the big one!) [kV/cm]")};

      fhicl::Atom<double> Electronlifetime{Name("Electronlifetime"),
                                           Comment("electron lifetime in liquid argon [us]")};
      fhicl::Atom<double> Temperature{Name("Temperature"), Comment("argon temperature [K]")};
      fhicl::Atom<double> ElectronsToADC{
        Name("ElectronsToADC"),
        Comment("conversion factor: (ADC counts)/(ionization electrons)")};
      fhicl::Atom<unsigned int> NumberTimeSamples{
        Name("NumberTimeSamples"),
        Comment("number of TPC readout TDC clock ticks per event")};
      fhicl::Atom<unsigned int> ReadOutWindowSize{
        Name("ReadOutWindowSize"),
        Comment("number of TPC readout TDC clock ticks per readout window")};

      // The following are not really "optional": the ones for the views which
      // are present are mandatory.
      fhicl::OptionalAtom<double> TimeOffsetU{
        Name("TimeOffsetU"),
        Comment("tick offset subtracted to to convert spacepoint coordinates "
                "to hit times on view U")};
      fhicl::OptionalAtom<double> TimeOffsetV{
        Name("TimeOffsetV"),
        Comment("tick offset subtracted to to convert spacepoint coordinates "
                "to hit times on view V")};
      fhicl::OptionalAtom<double> TimeOffsetZ{
        Name("TimeOffsetZ"),
        Comment("tick offset subtracted to to convert spacepoint coordinates "
                "to hit times on view Z")};
      fhicl::OptionalAtom<double> TimeOffsetY{
        Name("TimeOffsetY"),
        Comment("tick offset subtracted to to convert spacepoint coordinates "
                "to hit times on view Y")};
      fhicl::OptionalAtom<double> TimeOffsetX{
        Name("TimeOffsetX"),
        Comment("tick offset subtracted to to convert spacepoint coordinates "
                "to hit times on view X")};

      fhicl::Atom<double> SternheimerA{
        Name("SternheimerA"),
        Comment("parameter a of Sternheimer correction delta = 2log(10) x - "
                "cbar + { a (x1-x)^k } theta(x1-x), x = log10(p/m)")};
      fhicl::Atom<double> SternheimerK{
        Name("SternheimerK"),
        Comment("parameter k of Sternheimer correction delta = 2log(10) x - "
                "cbar + { a (x_1-x)^k } theta(x1-x), x = log10(p/m)")};
      fhicl::Atom<double> SternheimerX0{
        Name("SternheimerX0"),
        Comment("minimum x = log10(p/m) for the application of Sternheimer "
                "correction")};
      fhicl::Atom<double> SternheimerX1{
        Name("SternheimerX1"),
        Comment("parameter x_1 of Sternheimer correction delta = 2log(10) x - "
                "cbar + { a (x_1-x)^k } theta(x1-x), x = log10(p/m)")};
      fhicl::Atom<double> SternheimerCbar{
        Name("SternheimerCbar"),
        Comment("parameter cbar of Sternheimer correction delta = 2log(10) x - "
                "cbar + { a (x_1-x)^k } theta(x1-x), x = log10(p/m)")};
      fhicl::Atom<double> DriftVelFudgeFactor{
        Name("DriftVelFudgeFactor"),
        Comment("Allows a scaling factor to fudge the drift velocity "
                "calculation (as suggested by DriftVel Stancari")};

      fhicl::Atom<bool> UseIcarusMicrobooneDriftModel{
        Name("UseIcarusMicrobooneDriftModel"),
        Comment("Allows user to decide to use the ICARUS+MicroBooNE drift "
                "model for velocity calculation as in arXiv:2008.09765"),
        false};

      fhicl::Atom<bool> IncludeInterPlanePitchInXTickOffsets{
        Name("IncludeInterPlanePitchInXTickOffsets"),
        Comment("Historically, ConvertTicksToX has allowed for the drift time "
                "between the wire planes. This is appropriate for "
                "recob::RawDigits, and recob::Wires from the 1D unfolding, "
                "but is not appropriate for recob::Wires from WireCell. "
                "The default value is 'true', retaining the 'classic' behaviour"),
        true};

      fhicl::Atom<bool> SimpleBoundary{Name("SimpleBoundaryProcess"), Comment("")};

      fhicl::Atom<double> ModBoxAlpha{
        Name("ModBoxAlpha"),
        Comment("alpha parameter in the Modified Box recombination model."),
        util::kModBoxA};

      fhicl::Atom<double> ModBoxBeta{
        Name("ModBoxBeta"),
        Comment("beta parameter in the Modified Box recombination model."),
        util::kModBoxB};

    }; // Configuration_t

    DetectorPropertiesICARUSClockOffsetMC(fhicl::ParameterSet const& pset,
                                          const geo::GeometryCore* geo,
                                          const detinfo::LArProperties* lp,
                                          std::set<std::string> const& ignore_params = {});

    DetectorPropertiesICARUSClockOffsetMC(DetectorPropertiesICARUSClockOffsetMC const&) = delete;
    virtual ~DetectorPropertiesICARUSClockOffsetMC() = default;

    void SetNumberTimeSamples(unsigned int nsamp) { fNumberTimeSamples = nsamp; }

    // Accessors.

    double Efield(unsigned int planegap = 0) const override; ///< kV/cm

    double DriftVelocity(double efield = 0.,
                         double temperature = 0.) const override; ///< cm/us

    /// dQ/dX in electrons/cm, returns dE/dX in MeV/cm.
    double BirksCorrection(double dQdX) const override;
    double BirksCorrection(double dQdX, double EField) const override;
    double ModBoxCorrection(double dQdX) const override;
    double ModBoxCorrection(double dQdX, double EField) const override;

    double ElectronLifetime() const override
    {
      return fElectronlifetime; //< microseconds
    }

    /**
     * @brief Returns argon density at a given temperature
     * @param temperature the temperature in kelvin
     * @return argon density in g/cm^3
     *
     * Density is nearly a linear function of temperature.
     * See the NIST tables for details
     * Slope is between -6.2 and -6.1, intercept is 1928 kg/m^3.
     * This parameterization will be good to better than 0.5%.
     */
    double Density(double temperature = 0.) const override; ///< g/cm^3

    /// In kelvin.
    double Temperature() const override { return fTemperature; }

    /**
     * @brief Restricted mean energy loss (dE/dx)
     * @param mom  momentum of incident particle [GeV/c]
     * @param mass mass of incident particle [GeV/c^2]
     * @param tcut maximum kinetic energy of delta rays [MeV]; 0 for unlimited
     * @return the restricted mean energy loss (dE/dx) in units of MeV/cm
     *
     * Returned value is always positive.
     * For unrestricted mean energy loss, set tcut = 0 (special case),
     * or tcut large.
     *
     * Based on Bethe-Bloch formula as contained in particle data book.
     * Material parameters are from the configuration.
     */
    double Eloss(double mom, double mass, double tcut) const override;

    /**
     * @brief Energy loss fluctuation (@f$ \sigma_{E}^2 / x @f$)
     * @param mom  momentum of incident particle in [GeV/c]
     * @param mass mass of incident particle [GeV/c^2]
     * @return energy loss fluctuation in MeV^2/cm
     *
     * Based on Bichsel formula referred to but not given in PDG.
     */
    double ElossVar(double mom, double mass) const override;

    double ElectronsToADC() const override { return fElectronsToADC; }
    unsigned int NumberTimeSamples() const override { return fNumberTimeSamples; }
    unsigned int ReadOutWindowSize() const override { return fReadOutWindowSize; }
    double TimeOffsetU() const override { return fTimeOffsetU; };
    double TimeOffsetV() const override { return fTimeOffsetV; };
    double TimeOffsetZ() const override { return fTimeOffsetZ; };
    double TimeOffsetY() const override { return fTimeOffsetY; };

    bool SimpleBoundary() const override { return fSimpleBoundary; }

    DetectorPropertiesData DataFor(detinfo::DetectorClocksData const& clock_data) const override;

  private:
    /**
     * @brief Configures the provider, first validating the configuration
     * @param p configuration parameter set
     * @param ignore_params parameters to be ignored (optional)
     *
     * This method will validate the parameter set (except for the parameters
     * it's explicitly told to ignore) and extract the useful information out
     * of it.
     */
    void ValidateAndConfigure(fhicl::ParameterSet const& p,
                              std::set<std::string> const& ignore_params);

    std::string CheckTimeOffsets(std::set<geo::View_t> const& requested_views) const;

    /// Parameters for Sternheimer density effect corrections
    struct SternheimerParameters_t {
      double a;    ///< parameter a
      double k;    ///< parameter k
      double x0;   ///< parameter x0
      double x1;   ///< parameter x1
      double cbar; ///< parameter Cbar
    };

    // service providers we depend on;
    // in principle could be replaced by a single providerpack_type.
    const detinfo::LArProperties* fLP;
    const geo::GeometryCore* fGeo;

    std::vector<double> fEfield;     ///< kV/cm (per inter-plane volume) !
    double fElectronlifetime;        ///< microseconds
    double fTemperature;             ///< kelvin
    double fElectronsToADC;          ///< conversion factor for # of ionization electrons
                                     ///< to 1 ADC count
    unsigned int fNumberTimeSamples; ///< number of clock ticks per event
    unsigned int fReadOutWindowSize; ///< number of clock ticks per readout window
    double fTimeOffsetU;             ///< time offset to convert spacepoint coordinates to
                                     ///< hit times on view U
    double fTimeOffsetV;             ///< time offset to convert spacepoint coordinates to
                                     ///< hit times on view V
    double fTimeOffsetZ;             ///< time offset to convert spacepoint coordinates to
                                     ///< hit times on view Z
    double fTimeOffsetY;             ///< time offset to convert spacepoint coordinates to
                                     ///< hit times on view Y
    double fTimeOffsetX;             ///< time offset to convert spacepoint coordinates to
                                     ///< hit times on view X
    double fDriftVelFudgeFactor;     ///< Scaling factor to allow "fudging" of drift
                                     ///< velocity

    bool fUseIcarusMicrobooneDriftModel; ///< if true, use alternative ICARUS-MicroBooNE drift
                                         ///< model instead of Walkowiak-based one

    /// Historically, ConvertTicksToX has allowed for the drift time between
    /// the wire planes. This is appropriate for recob::RawDigits, and
    /// recob::Wires from the 1D unfolding, but is not appropriate for
    /// recob::Wires from WireCell.
    bool fIncludeInterPlanePitchInXTickOffsets;

    SternheimerParameters_t fSternheimerParameters; ///< Sternheimer parameters

    std::vector<std::vector<double>> fDriftDirection;

    bool fSimpleBoundary;

    double fModBoxA;
    double fModBoxB;

  }; // class DetectorPropertiesICARUSClockOffsetMC
} // namespace detinfo

#endif // DETINFO_DETECTORPROPERTIESICARUS_H
