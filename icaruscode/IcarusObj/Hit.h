/** ****************************************************************************
 * @file lardataobj/RecoBase/Hit.h
 * @brief Declaration of signal hit object.
 * @author mitchell.soderberg@yale.edu
 * @see  lardataobj/RecoBase/Hit.cxx
 *
 * Changes:
 * 20141212 Gianluca Petrillo (petrillo@fnal.gov)
 *   data architecture revision changes (v13 -> v14):
 *   - fRawDigit and RawDigit() removed
 *   - fWire and Wire() removed
 *   - fSignalType and SignalType() removed
 *   - fChannel and Channel() added
 *   - constructors now take directly a RawDigit, not its art::Ptr
 * 20150129 Gianluca Petrillo (petrillo@fnal.gov)
 *   data architecture revision changes (v14 -> v15):
 *   - removed fHitSignal
 *
 * ****************************************************************************/

#ifndef ICARUSOBJ_RECOBASE_HIT_H
#define ICARUSOBJ_RECOBASE_HIT_H

// C/C++ standard librraies
#include <iosfwd>

// LArSoft libraries
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h" // geo::View_t, geo::SignalType, geo::WireID

namespace icarus {

  /**
   * @brief 2D representation of charge deposited in the TDC/wire plane
   *
   * Hits are assumed to be made from deconvoluted, unipolar, calibrated
   * signals.
   * They identify a charge deposit in a specific location and time;
   * the location is absolute and unique in the detector, while the time is
   * relative to the start of sampling (tick time).
   *
   * The version 14 of icarus::Hit introduces the following changes:
   * - StartTime() becomes StartTick(), StopTime() becomes StopTick()
   * - Charge(true) is now PeakAmplitude(), Charge(false) (default) is Integral()
   */
  class Hit {
  public:
    /// Default constructor: a hit with no signal
    Hit();

  private:
    raw::ChannelID_t fChannel; ///< ID of the readout channel the hit was extracted from
    raw::TDCtick_t fStartTick; ///< initial tdc tick for hit
    raw::TDCtick_t fEndTick;   ///< final tdc tick for hit
    float fPeakTime;           ///< time of the signal peak, in tick units
    float fSigmaPeakTime;      ///< uncertainty for the signal peak, in tick units
    float fRMS;                ///< RMS of the hit shape, in tick units
    float fPeakAmplitude;      ///< the estimated amplitude of the hit at its peak, in ADC units
    float
      fSigmaPeakAmplitude; ///< uncertainty on estimated amplitude of the hit at its peak, in ADC units
    float fSummedADC;      ///< the sum of calibrated ADC counts of the hit
    float
      fIntegral; ///< the integral under the calibrated signal waveform of the hit, in tick x ADC units
    float
      fSigmaIntegral; ///< the uncertainty of integral under the calibrated signal waveform of the hit, in ADC units
    short int fMultiplicity; ///< how many hits could this one be shared with
    short int fLocalIndex; ///< index of this hit among the Multiplicity() hits in the signal window
    float fGoodnessOfFit;  ///< how well do we believe we know this hit?
    int fNDF;              ///< degrees of freedom in the determination of the hit shape
    geo::View_t fView;     ///< view for the plane of the hit
    geo::SigType_t fSignalType; ///< signal type for the plane of the hit
    geo::WireID fWireID;        ///< WireID for the hit (Cryostat, TPC, Plane, Wire)

    friend class HitCreator; // helper to create hits

  public:
    /**
       * @brief Constructor: directly sets all the fields
       * @param channel         ID of the readout channel the hit was extracted from
       * @param start_tick      initial tdc tick for hit
       * @param end_tick        final tdc tick for hit
       * @param peak_time       tdc for the peak charge deposition
       * @param sigma_peak_time uncertainty for tdc of the peak
       * @param rms             RMS of the hit shape
       * @param peak_amplitude  the estimated amplitude of the hit at its peak
       * @param sigma_peak_amplitude  the uncertainty on the estimated amplitude of the hit at its peak
       * @param summedADC       the sum of calibrated ADC counts of the hit
       * @param hit_integral    the integral under the calibrated signal waveform of the hit
       * @param hit_sigma_integral uncertainty on the integral under the calibrated signal waveform of the hit
       * @param multiplicity    how many hits could this one be shared with
       * @param local_index     index of this hit among the Multiplicity() hits in the signal window
       * @param goodness_of_fit how well do we believe we know this hit?
       * @param dof             number of degrees of freedom in the determination of hit shape
       * @param view            view for the plane of the hit
       * @param signal_type     signal type for the plane of the hit
       * @param wireID          WireID for the hit (Cryostat  TPC  Plane  Wire)
       *
       * The tick parameters are real numbers, since they can in principle come
       * from some processing.
       */
    Hit(raw::ChannelID_t channel,
        raw::TDCtick_t start_tick,
        raw::TDCtick_t end_tick,
        float peak_time,
        float sigma_peak_time,
        float rms,
        float peak_amplitude,
        float sigma_peak_amplitude,
        float summedADC,
        float hit_integral,
        float hit_sigma_integral,
        short int multiplicity,
        short int local_index,
        float goodness_of_fit,
        int dof,
        geo::View_t view,
        geo::SigType_t signal_type,
        geo::WireID wireID);

    /// @{
    /// @name Accessors

    /// Initial tdc tick for hit
    raw::TDCtick_t StartTick() const;

    /// Final tdc tick for hit
    raw::TDCtick_t EndTick() const;

    /// Time of the signal peak, in tick units
    float PeakTime() const;

    /// Uncertainty for the signal peak, in tick units
    float SigmaPeakTime() const;

    /// RMS of the hit shape, in tick units
    float RMS() const;

    /// The estimated amplitude of the hit at its peak, in ADC units
    float PeakAmplitude() const;

    /// Uncertainty on estimated amplitude of the hit at its peak, in ADC units
    float SigmaPeakAmplitude() const;

    /// The sum of calibrated ADC counts of the hit (0. by default)
    float SummedADC() const;

    /// Integral under the calibrated signal waveform of the hit, in tick x ADC units
    float Integral() const;

    ///< Uncertainty of integral under the calibrated signal waveform of the hit, in ADC units
    float SigmaIntegral() const;

    /// How many hits could this one be shared with
    short int Multiplicity() const;

    ///< Index of this hit among the Multiplicity() hits in the signal window (-1 by default)
    short int LocalIndex() const;

    ///< How well do we believe we know this hit?
    float GoodnessOfFit() const;

    ///< Degrees of freedom in the determination of the hit signal shape (-1 by default)
    int DegreesOfFreedom() const;

    /// ID of the readout channel the hit was extracted from
    raw::ChannelID_t Channel() const;

    /// View for the plane of the hit
    geo::View_t View() const;

    /// Signal type for the plane of the hit
    geo::SigType_t SignalType() const;

    ///< ID of the wire the hit is on (Cryostat, TPC, Plane, Wire)
    geo::WireID const& WireID() const;

    /// @}

    //@{
    /**
       * @brief Returns a time sigmas RMS away from the peak time
       * @param sigmas the number of RMS units to move away
       * @return the shifted time in TDC ticks
       *
       * PeakTimePlusRMS() returns PeakTime() + sigmas x RMS();
       * PeakTimeMinusRMS() returns PeakTime() - sigmas x RMS().
       *
       * @note StartTime() of icarus::Hit version <=13 was defined by
       *   GausHitFinder to be PeakTimePlusRMS(-1.), and EndTime() was
       *   PeakTimePlusRMS(+1.).
       */
    float PeakTimePlusRMS(float sigmas = +1.) const;
    float PeakTimeMinusRMS(float sigmas = +1.) const;
    //@}

    /**
       * @brief Returns the distance of the specified time from peak, in RMS units
       * @param time the time, in TDC tick units
       * @return the distance of the specified time from peak, in RMS units
       *
       * This returns (ticks - PeakTime()) / RMS().
       * There is no protection in case RMS is 0!
       */
    float TimeDistanceAsRMS(float time) const;

    friend std::ostream& operator<<(std::ostream& o, const Hit& a);
    friend bool operator<(const Hit& a, const Hit& b);

  }; // class Hit
} // namespace icarus

inline raw::TDCtick_t icarus::Hit::StartTick() const
{
  return fStartTick;
}
inline raw::TDCtick_t icarus::Hit::EndTick() const
{
  return fEndTick;
}
inline float icarus::Hit::PeakTime() const
{
  return fPeakTime;
}
inline float icarus::Hit::SigmaPeakTime() const
{
  return fSigmaPeakTime;
}
inline float icarus::Hit::RMS() const
{
  return fRMS;
}
inline float icarus::Hit::PeakAmplitude() const
{
  return fPeakAmplitude;
}
inline float icarus::Hit::SigmaPeakAmplitude() const
{
  return fSigmaPeakAmplitude;
}
inline float icarus::Hit::SummedADC() const
{
  return fSummedADC;
}
inline float icarus::Hit::Integral() const
{
  return fIntegral;
}
inline float icarus::Hit::SigmaIntegral() const
{
  return fSigmaIntegral;
}
inline short int icarus::Hit::Multiplicity() const
{
  return fMultiplicity;
}
inline short int icarus::Hit::LocalIndex() const
{
  return fLocalIndex;
}
inline float icarus::Hit::GoodnessOfFit() const
{
  return fGoodnessOfFit;
}
inline int icarus::Hit::DegreesOfFreedom() const
{
  return fNDF;
}
inline raw::ChannelID_t icarus::Hit::Channel() const
{
  return fChannel;
}
inline geo::SigType_t icarus::Hit::SignalType() const
{
  return fSignalType;
}
inline geo::View_t icarus::Hit::View() const
{
  return fView;
}
inline geo::WireID const& icarus::Hit::WireID() const
{
  return fWireID;
}

inline float icarus::Hit::PeakTimePlusRMS(float sigmas /* = +1. */) const
{
  return PeakTime() + sigmas * RMS();
}

inline float icarus::Hit::PeakTimeMinusRMS(float sigmas /* = +1. */) const
{
  return PeakTimePlusRMS(-sigmas);
}

inline float icarus::Hit::TimeDistanceAsRMS(float time) const
{
  return (time - PeakTime()) / RMS();
}

#endif // LARDATAOBJ_RECOBASE_HIT_H
