/**
 * @file lardataalg/DetectorInfo/DetectorPropertiesICARUSClockOffsetMC.cxx
 * @brief DetectorPropertiesStandard but adding in a shift
 *        for the trigger time to be used with ICARUS MC *when* 
 *        shifting of sim products themselves has not taken place.
 * @author Bruce Howard and Jaesung Kim, building on the original module (DetectorPropertiesStandard) by Jonathan Paley (jpaley@fnal.gov)
 */

// Framework includes

// icaruscode includes
#include "icaruscode/Utilities/DetectorPropertiesICARUSClockOffsetMC.h"

// LArSoft includes
#include "larcorealg/CoreUtils/ProviderUtil.h" // lar::IgnorableProviderConfigKeys()
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"

#include "cetlib/pow.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Art includes
#include "fhiclcpp/types/Table.h"

// C/C++ libraries
#include <sstream> // std::ostringstream

namespace detinfo {

  //--------------------------------------------------------------------
  DetectorPropertiesICARUSClockOffsetMC::DetectorPropertiesICARUSClockOffsetMC(fhicl::ParameterSet const& pset,
                                                                               const geo::GeometryCore* geo,
                                                                               const detinfo::LArProperties* lp,
                                                                               std::set<std::string> const& ignore_params)
    : fLP(lp), fGeo(geo)
  {
    ValidateAndConfigure(pset, ignore_params);
  }

  //--------------------------------------------------------------------
  void DetectorPropertiesICARUSClockOffsetMC::ValidateAndConfigure(fhicl::ParameterSet const& p,
                                                                   std::set<std::string> const& ignore_params)
  {
    {
      mf::LogInfo debug("setupProvider<DetectorPropertiesICARUSClockOffsetMC>");

      debug << "Asked to ignore " << ignore_params.size() << " keys:";
      for (auto const& key : ignore_params)
        debug << " '" << key << "'";
    }

    std::set<std::string> ignorable_keys = lar::IgnorableProviderConfigKeys();
    ignorable_keys.insert(ignore_params.begin(), ignore_params.end());

    // parses and validates the parameter set:
    fhicl::Table<Configuration_t> const config{p, ignorable_keys};

    fEfield = config().Efield();
    fElectronlifetime = config().Electronlifetime();
    fTemperature = config().Temperature();
    fElectronsToADC = config().ElectronsToADC();
    fNumberTimeSamples = config().NumberTimeSamples();
    fReadOutWindowSize = config().ReadOutWindowSize();

    std::set<geo::View_t> present_views;
    if (config().TimeOffsetU(fTimeOffsetU)) present_views.insert(geo::kU);
    if (config().TimeOffsetV(fTimeOffsetV)) present_views.insert(geo::kV);
    if (config().TimeOffsetZ(fTimeOffsetZ)) present_views.insert(geo::kZ);
    if (config().TimeOffsetY(fTimeOffsetY)) present_views.insert(geo::kY);
    if (config().TimeOffsetX(fTimeOffsetX)) present_views.insert(geo::kX);

    std::string const errors = CheckTimeOffsets(present_views);
    if (!errors.empty()) {
      throw cet::exception("DetectorPropertiesICARUSClockOffsetMC") << "Detected configuration errors: \n"
                                                                    << errors;
    }

    fSternheimerParameters.a = config().SternheimerA();
    fSternheimerParameters.k = config().SternheimerK();
    fSternheimerParameters.x0 = config().SternheimerX0();
    fSternheimerParameters.x1 = config().SternheimerX1();
    fSternheimerParameters.cbar = config().SternheimerCbar();

    fDriftVelFudgeFactor = config().DriftVelFudgeFactor();

    fUseIcarusMicrobooneDriftModel = config().UseIcarusMicrobooneDriftModel();

    fIncludeInterPlanePitchInXTickOffsets = config().IncludeInterPlanePitchInXTickOffsets();

    fSimpleBoundary = config().SimpleBoundary();

    fModBoxA = config().ModBoxAlpha();
    fModBoxB = config().ModBoxBeta();
  }

  //------------------------------------------------------------------------------------//
  double DetectorPropertiesICARUSClockOffsetMC::Efield(unsigned int const planegap) const
  {
    if (planegap >= fEfield.size())
      throw cet::exception("DetectorPropertiesICARUSClockOffsetMC")
        << "requesting Electric field in a plane gap that is not defined\n";

    return fEfield[planegap];
  }

  //------------------------------------------------
  double DetectorPropertiesICARUSClockOffsetMC::Density(double temperature) const
  {
    // Default temperature use internal value.
    if (temperature == 0.) temperature = Temperature();

    return -0.00615 * temperature + 1.928;
  }

  //----------------------------------------------------------------------------------
  // Restricted mean energy loss (dE/dx) in units of MeV/cm.
  //
  // For unrestricted mean energy loss, set tcut = 0, or tcut large.
  //
  // Arguments:
  //
  // mom  - Momentum of incident particle in GeV/c.
  // mass - Mass of incident particle in GeV/c^2.
  // tcut - Maximum kinetic energy of delta rays (MeV).
  //
  // Returned value is positive.
  //
  // Based on Bethe-Bloch formula as contained in particle data book.
  // Material parameters (stored in larproperties.fcl) are taken from
  // pdg web site http://pdg.lbl.gov/AtomicNuclearProperties/
  //
  double DetectorPropertiesICARUSClockOffsetMC::Eloss(double mom, double mass, double tcut) const
  {
    // Some constants.
    constexpr double K = 0.307075;     // 4 pi N_A r_e^2 m_e c^2 (MeV cm^2/mol).
    constexpr double me = 0.510998918; // Electron mass (MeV/c^2).

    // Calculate kinematic quantities.
    double const bg = mom / mass;            // beta*gamma.
    double const gamma = sqrt(1. + bg * bg); // gamma.
    double const beta = bg / gamma;          // beta (velocity).
    double const mer = 0.001 * me / mass;    // electron mass / mass of incident particle.
    double const tmax =
      2. * me * bg * bg / (1. + 2. * gamma * mer + mer * mer); // Maximum delta ray energy (MeV).

    // Make sure tcut does not exceed tmax.
    if (tcut == 0. || tcut > tmax) tcut = tmax;

    // Calculate density effect correction (delta).
    double const x = std::log10(bg);
    double delta = 0.;
    if (x >= fSternheimerParameters.x0) {
      delta = 2. * std::log(10.) * x - fSternheimerParameters.cbar;
      if (x < fSternheimerParameters.x1)
        delta += fSternheimerParameters.a *
                 std::pow(fSternheimerParameters.x1 - x, fSternheimerParameters.k);
    }

    // Calculate stopping number.
    double B =
      0.5 * std::log(2. * me * bg * bg * tcut / (1.e-12 * cet::square(fLP->ExcitationEnergy()))) -
      0.5 * beta * beta * (1. + tcut / tmax) - 0.5 * delta;

    // Don't let the stopping number become negative.
    if (B < 1.) B = 1.;

    // Calculate dE/dx.
    return Density() * K * fLP->AtomicNumber() * B / (fLP->AtomicMass() * beta * beta);
  }

  //----------------------------------------------------------------------------------
  double DetectorPropertiesICARUSClockOffsetMC::ElossVar(double const mom, double const mass) const
  {
    // Some constants.
    constexpr double K = 0.307075;     // 4 pi N_A r_e^2 m_e c^2 (MeV cm^2/mol).
    constexpr double me = 0.510998918; // Electron mass (MeV/c^2).

    // Calculate kinematic quantities.
    double const bg = mom / mass;          // beta*gamma.
    double const gamma2 = 1. + bg * bg;    // gamma^2.
    double const beta2 = bg * bg / gamma2; // beta^2.
    return gamma2 * (1. - 0.5 * beta2) * me * (fLP->AtomicNumber() / fLP->AtomicMass()) * K *
           Density();
  }

  //------------------------------------------------------------------------------------//
  double DetectorPropertiesICARUSClockOffsetMC::DriftVelocity(double efield, double temperature) const
  {
    // Drift Velocity as a function of Electric Field and LAr Temperature
    // from : W. Walkowiak, NIM A 449 (2000) 288-294
    //
    // Option to use MicroBooNE+ICARUS model (as in arXiv:2008.09765) provided as
    // well, with temperature depenence as prescribed by Mike Mooney based on
    // looking at the Walkowiak data.
    //
    // Efield should have units of kV/cm
    // Temperature should have units of Kelvin

    // Default Efield, use internal value.
    if (efield == 0.) efield = Efield();

    if (efield > 4.0)
      mf::LogWarning("DetectorPropertiesStandard")
        << "DriftVelocity Warning! : E-field value of " << efield
        << " kV/cm is outside of range covered by drift"
        << " velocity parameterization. Returned value"
        << " may not be correct";

    // Default temperature use internal value.
    if (temperature == 0.) temperature = Temperature();

    if (temperature < 87.0 || temperature > 94.0)
      mf::LogWarning("DetectorPropertiesStandard")
        << "DriftVelocity Warning! : Temperature value of " << temperature
        << " K is outside of range covered by drift velocity"
        << " parameterization. Returned value may not be"
        << " correct";

    double vd;

    if (!fUseIcarusMicrobooneDriftModel) {
      double const tshift = -87.203 + temperature;
      double const xFit = 0.0938163 - 0.0052563 * tshift - 0.0001470 * tshift * tshift;
      double const uFit = 5.18406 + 0.01448 * tshift - 0.003497 * tshift * tshift -
                          0.000516 * tshift * tshift * tshift;

      // Icarus Parameter Set, use as default
      constexpr double P1 = -0.04640; // K^-1
      constexpr double P2 = 0.01712;  // K^-1
      constexpr double P3 = 1.88125;  // (kV/cm)^-1
      constexpr double P4 = 0.99408;  // kV/cm
      constexpr double P5 = 0.01172;  // (kV/cm)^-P6
      constexpr double P6 = 4.20214;
      constexpr double T0 = 105.749; // K

      // Walkowiak Parameter Set
      constexpr double P1W = -0.01481; // K^-1
      constexpr double P2W = -0.0075;  // K^-1
      constexpr double P3W = 0.141;    // (kV/cm)^-1
      constexpr double P4W = 12.4;     // kV/cm
      constexpr double P5W = 1.627;    // (kV/cm)^-P6
      constexpr double P6W = 0.317;
      constexpr double T0W = 90.371; // K

      // From Craig Thorne . . . currently not documented
      // smooth transition from linear at small fields to
      //     icarus fit at most fields to Walkowiak at very high fields
      if (efield < xFit)
        vd = efield * uFit;
      else if (efield < 0.619) {
        vd = ((P1 * (temperature - T0) + 1) *
                (P3 * efield * std::log(1 + P4 / efield) + P5 * std::pow(efield, P6)) +
              P2 * (temperature - T0));
      }
      else if (efield < 0.699) {
        vd = 12.5 * (efield - 0.619) *
               ((P1W * (temperature - T0W) + 1) *
                  (P3W * efield * std::log(1 + P4W / efield) + P5W * std::pow(efield, P6W)) +
                P2W * (temperature - T0W)) +
             12.5 * (0.699 - efield) *
               ((P1 * (temperature - T0) + 1) *
                  (P3 * efield * std::log(1 + P4 / efield) + P5 * std::pow(efield, P6)) +
                P2 * (temperature - T0));
      }
      else {
        vd = ((P1W * (temperature - T0W) + 1) *
                (P3W * efield * std::log(1 + P4W / efield) + P5W * std::pow(efield, P6W)) +
              P2W * (temperature - T0W));
      }
    }

    // MicroBooNE+ICARUS model (arXiv:2008.09765) with temperature dependence given by
    // Mike Mooney based on looking at Walkowiak data (NIM A 449 (2000) 288-294)
    if (fUseIcarusMicrobooneDriftModel) {
      constexpr double P0 = 0.;
      constexpr double P1 = 5.53416;
      constexpr double P2 = -6.53093;
      constexpr double P3 = 3.20752;
      constexpr double P4 = 0.389696;
      constexpr double P5 = -0.556184;
      vd = (1.0 - 0.0184 * (temperature - 89.0)) *
           (P0 + P1 * cet::pow<1>(efield) + P2 * cet::pow<2>(efield) + P3 * cet::pow<3>(efield) +
            P4 * cet::pow<4>(efield) + P5 * cet::pow<5>(efield));
    }

    vd *= fDriftVelFudgeFactor / 10.;

    return vd; // in cm/us
  }

  //----------------------------------------------------------------------------------
  // The below function assumes that the user has applied the lifetime
  // correction and effective pitch between the wires (usually after 3D
  // reconstruction). Using with mean wire pitch will not give correct results.
  // parameters:
  //  dQdX in electrons/cm, charge (amplitude or integral obtained) divided by
  //         effective pitch for a given 3D track.
  //  Electric Field in the drift region in KV/cm/
  // returns dEdX in MeV/cm
  double DetectorPropertiesICARUSClockOffsetMC::BirksCorrection(double dQdx) const
  {
    return BirksCorrection(dQdx, Efield());
  }
  double DetectorPropertiesICARUSClockOffsetMC::BirksCorrection(double dQdx, double E_field) const
  {
    // Correction for charge quenching using parameterization from
    // S.Amoruso et al., NIM A 523 (2004) 275

    constexpr double A3t = util::kRecombA;
    double K3t = util::kRecombk;                                    // in KV/cm*(g/cm^2)/MeV
    double const rho = Density();                                   // LAr density in g/cm^3
    constexpr double Wion = 1000. / util::kGeVToElectrons;          // 23.6 eV = 1e, Wion in MeV/e
    K3t /= rho;                                                     // KV/MeV
    double const dEdx = dQdx / (A3t / Wion - K3t / E_field * dQdx); // MeV/cm

    return dEdx;
  }

  //----------------------------------------------------------------------------------
  // Modified Box model correction
  double DetectorPropertiesICARUSClockOffsetMC::ModBoxCorrection(double dQdx) const
  {
    return ModBoxCorrection(dQdx, Efield());
  }
  double DetectorPropertiesICARUSClockOffsetMC::ModBoxCorrection(double dQdx, double E_field) const
  {
    // Modified Box model correction has better behavior than the Birks
    // correction at high values of dQ/dx.
    double const rho = Density();                          // LAr density in g/cm^3
    constexpr double Wion = 1000. / util::kGeVToElectrons; // 23.6 eV = 1e, Wion in MeV/e
    double const Beta = fModBoxB / (rho * E_field);
    double const Alpha = fModBoxA;
    double const dEdx = (exp(Beta * Wion * dQdx) - Alpha) / Beta;

    return dEdx;
  }

  //--------------------------------------------------------------------
  //  x<--> ticks conversion methods
  //
  //  Ben Jones April 2012,
  //  based on code by Herb Greenlee in SpacePointService
  //

  //--------------------------------------------------------------------
  // Recalculate x<-->ticks conversion parameters from detector constants

  DetectorPropertiesData DetectorPropertiesICARUSClockOffsetMC::DataFor(
    detinfo::DetectorClocksData const& clock_data) const
  {
    double const samplingRate = sampling_rate(clock_data);
    double const efield = Efield();
    double const temperature = Temperature();
    double const driftVelocity = DriftVelocity(efield, temperature);
    double const x_ticks_coefficient = 0.001 * driftVelocity * samplingRate;

    double extraOffset = 0.;
    if ( clock_data.TriggerTime() > std::numeric_limits<double>::lowest() + std::numeric_limits<double>::epsilon() &&
         clock_data.TriggerTime() < std::numeric_limits<double>::max() - std::numeric_limits<double>::epsilon()
    ){
      extraOffset = clock_data.TriggerTime()-clock_data.BeamGateTime();
    }

    double const triggerOffset = trigger_offset(clock_data) + clock_data.TPCClock().Ticks(extraOffset);

    std::vector<std::vector<std::vector<double>>> x_ticks_offsets(fGeo->Ncryostats());
    std::vector<std::vector<double>> drift_direction(fGeo->Ncryostats());

    for (size_t cstat = 0; cstat < fGeo->Ncryostats(); ++cstat) {
      auto const& cryostat = fGeo->Cryostat(geo::CryostatID(cstat));
      x_ticks_offsets[cstat].resize(cryostat.NTPC());
      drift_direction[cstat].resize(cryostat.NTPC());

      for (size_t tpc = 0; tpc < cryostat.NTPC(); ++tpc) {
        const geo::TPCGeo& tpcgeom = cryostat.TPC(tpc);

        const double dir((tpcgeom.DriftDirection() == geo::kNegX) ? +1.0 : -1.0);
        drift_direction[cstat][tpc] = dir;

        int nplane = tpcgeom.Nplanes();
        x_ticks_offsets[cstat][tpc].resize(nplane, 0.);
        for (int plane = 0; plane < nplane; ++plane) {
          const geo::PlaneGeo& pgeom = tpcgeom.Plane(plane);

          // Calculate geometric time offset.
          // only works if xyz[0]<=0
          auto const xyz = tpcgeom.Plane(0).GetCenter();

          x_ticks_offsets[cstat][tpc][plane] =
            -xyz.X() / (dir * x_ticks_coefficient) + triggerOffset;

          if (fIncludeInterPlanePitchInXTickOffsets) {
            // Get field in gap between planes
            double efieldgap[3];
            double driftVelocitygap[3];
            double x_ticks_coefficient_gap[3];
            for (int igap = 0; igap < 3; ++igap) {
              efieldgap[igap] = Efield(igap);
              driftVelocitygap[igap] = DriftVelocity(efieldgap[igap], temperature);
              x_ticks_coefficient_gap[igap] = 0.001 * driftVelocitygap[igap] * samplingRate;
            }

            if (nplane == 3) {
              /*
                |    ---------- plane = 2 (collection)
                |                      Coeff[2]
                |    ---------- plane = 1 (2nd induction)
                |                      Coeff[1]
                |    ---------- plane = 0 (1st induction) x = xyz[0]
                |                      Coeff[0]
                |    ---------- x = 0
                V     For plane = 0, t offset is -xyz[0]/Coeff[0]
                x   */
              for (int ip = 0; ip < plane; ++ip) {
                x_ticks_offsets[cstat][tpc][plane] +=
                  tpcgeom.PlanePitch(ip, ip + 1) / x_ticks_coefficient_gap[ip + 1];
              }
            }
            else if (nplane == 2) { ///< special case for ArgoNeuT
              /*
                |    ---------- plane = 1 (collection)
                |                      Coeff[2]
                |    ---------- plane = 0 (2nd induction) x = xyz[0]
                |    ---------- x = 0, Coeff[1]
                V    ---------- first induction plane
                x                      Coeff[0]
                For plane = 0, t offset is pitch/Coeff[1] -
                (pitch+xyz[0])/Coeff[0] = -xyz[0]/Coeff[0] -
                pitch*(1/Coeff[0]-1/Coeff[1])
              */
              for (int ip = 0; ip < plane; ++ip) {
                x_ticks_offsets[cstat][tpc][plane] +=
                  tpcgeom.PlanePitch(ip, ip + 1) / x_ticks_coefficient_gap[ip + 2];
              }
              x_ticks_offsets[cstat][tpc][plane] -=
                tpcgeom.PlanePitch() * (1 / x_ticks_coefficient - 1 / x_ticks_coefficient_gap[1]);
            }

          } // end if fIncludeInterPlanePitchInXTickOffsets

          // Add view dependent offset
          // FIXME the offset should be plane-dependent
          geo::View_t view = pgeom.View();
          switch (view) {
          case geo::kU: x_ticks_offsets[cstat][tpc][plane] += fTimeOffsetU; break;
          case geo::kV: x_ticks_offsets[cstat][tpc][plane] += fTimeOffsetV; break;
          case geo::kZ: x_ticks_offsets[cstat][tpc][plane] += fTimeOffsetZ; break;
          case geo::kY: x_ticks_offsets[cstat][tpc][plane] += fTimeOffsetY; break;
          case geo::kX: x_ticks_offsets[cstat][tpc][plane] += fTimeOffsetX; break;
          default: throw cet::exception(__FUNCTION__) << "Bad view = " << view << "\n";
          } // switch
        }
      }
    }

    return DetectorPropertiesData{
      *this, x_ticks_coefficient, move(x_ticks_offsets), move(drift_direction)};
  }

  std::string DetectorPropertiesICARUSClockOffsetMC::CheckTimeOffsets(
    std::set<geo::View_t> const& requested_views) const
  {
    auto const& present_views = fGeo->Views();

    auto view_diff = [&present_views, &requested_views](geo::View_t const view) {
      return static_cast<int>(present_views.count(view)) -
             static_cast<int>(requested_views.count(view));
    };

    // It is not an error to specify an offset if the view does not
    // exist.  However, if a view does exist, and an offset does not,
    // then that will end the job.
    std::ostringstream errors;
    if (auto diff = view_diff(geo::kU); diff > 0) { errors << "TimeOffsetU missing for view U.\n"; }
    if (auto diff = view_diff(geo::kV); diff > 0) { errors << "TimeOffsetV missing for view V.\n"; }
    if (auto diff = view_diff(geo::kZ); diff > 0) { errors << "TimeOffsetZ missing for view Z.\n"; }
    if (auto diff = view_diff(geo::kY); diff > 0) { errors << "TimeOffsetY missing for view Y.\n"; }
    if (auto diff = view_diff(geo::kX); diff > 0) { errors << "TimeOffsetX missing for view X.\n"; }
    return errors.str();
  }
} // namespace
