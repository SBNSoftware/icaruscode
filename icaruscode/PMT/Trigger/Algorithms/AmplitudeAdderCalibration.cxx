/**
 * @file   icaruscode/PMT/Trigger/Algorithms/AmplitudeAdderCalibration.cxx
 * @brief  Calibration for adder board output based on signal amplitude.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   September 5, 2025
 * @see    icaruscode/PMT/Trigger/Algorithms/AmplitudeAdderCalibration.h
 */

#undef NDEBUG // FIXME

// library header
#include "icaruscode/PMT/Trigger/Algorithms/AmplitudeAdderCalibration.h"

// C/C++ standard libraries
#include <algorithm> // std::sort(), std::lower_bound()
#include <iomanip> // std::setw()
#include <utility> // std::move()
#include <memory> // std::make_unique()


// -----------------------------------------------------------------------------
namespace {
  template <typename T> constexpr T sqr(T v) noexcept { return v*v; }
}


// -----------------------------------------------------------------------------

icarus::trigger::AmplitudeAdderCalibration::ChannelCalibration_t::ChannelCalibration_t(
  AdderChannelID channel,
  double slope, millivolt offset,
  double slopeErr /* = 0.0 */, millivolt offsetErr /* = 0_mV */,
  millivolt covariance /* = 0_mV */,
  microseconds timeOffset /* = 0_us */,
  microseconds timeOffsetErr /* = 0_us */
)
  : channel{ channel }
  , parameters{ offset.value(), slope }
  , covariance
    { sqr(offsetErr.value()), covariance.value(), covariance.value(), sqr(slopeErr) }
  , timeOffset{ timeOffset }
  , timeOffsetErr{ timeOffsetErr }
{}


// -----------------------------------------------------------------------------
double icarus::trigger::AmplitudeAdderCalibration::ChannelCalibration_t::factor
  (millivolt amplitude) const
{
  // computeAmplitude(amplitude) / amplitude
  return slope() + offset() / amplitude;
}


// -----------------------------------------------------------------------------
double icarus::trigger::AmplitudeAdderCalibration::ChannelCalibration_t::factorUncertainty
  (millivolt amplitude) const
{
  // computeUncertainty(amplitude) / amplitude
  return slopeVar() / sqr(amplitude.value())
    + offsetSlopeCovar() / amplitude.value() + offsetVar();
}


// -----------------------------------------------------------------------------
double icarus::trigger::AmplitudeAdderCalibration::ChannelCalibration_t::factorRelUncertainty
  (millivolt amplitude) const
{
  return computeUncertainty(amplitude) / computeAmplitude(amplitude);
}


// -----------------------------------------------------------------------------
auto icarus::trigger::AmplitudeAdderCalibration::ChannelCalibration_t::computeAmplitude
  (millivolt amplitude) const -> millivolt
{
  return offset() + slope() * amplitude;
}


// -----------------------------------------------------------------------------
auto icarus::trigger::AmplitudeAdderCalibration::ChannelCalibration_t::computeUncertainty
  (millivolt amplitude) const -> millivolt
{
  return millivolt{ std::sqrt(computeVariance(amplitude)) };
}


// -----------------------------------------------------------------------------
double icarus::trigger::AmplitudeAdderCalibration::ChannelCalibration_t::computeVariance
  (millivolt amplitude) const
{
  // not very robust against rounding errors
  return covariance[0][0] + covariance[0][1] * amplitude.value()
    + covariance[1][1] * sqr(amplitude.value());
}


// -----------------------------------------------------------------------------
// ---  Interface implementations
// -----------------------------------------------------------------------------
// NOTE: these classes could be moved later to their own header/source
//       (especially if they are going to be toolized).

// C/C++ standard libraries
#include <algorithm> // std::lower_bound()


// -----------------------------------------------------------------------------
icarus::trigger::AmplitudeAdderCalibration::RunCalibration::RunCalibration
  (int run, CalibrationConstants_t calibration, std::string logCategory)
  : Base_t{ run, std::move(logCategory) }, fCalibration{ std::move(calibration) }
{
  std::sort(
    fCalibration.begin(), fCalibration.end(),
    details::ChannelCalibration_t_sorter{}
    );
}


// -----------------------------------------------------------------------------
double icarus::trigger::AmplitudeAdderCalibration::RunCalibration::doCalibrationFactor
  (AdderChannelID channel, WaveformSamples_t const& waveform) const
{
  // find the calibration constants for the channel
  ChannelCalibration_t const& calib = fetchChannel(channel);
  
  // compute the peak amplitude (assume 0 baseline)
  Voltage_t const peak = computePeakAmplitude(waveform);
  mfLogTrace() << "Waveform amplitude for channel " << channel
    << " is " << peak << " mV";
  
  return calib.factor(millivolt{ peak });
  
} // icarus::trigger::AmplitudeAdderCalibration::RunCalibration::doCalibrationFactor()


// -----------------------------------------------------------------------------
auto icarus::trigger::AmplitudeAdderCalibration::RunCalibration::doTimeOffset
  (AdderChannelID channel, WaveformSamples_t const&) const -> microseconds
{
  return fetchChannel(channel).timeOffset;
} // icarus::trigger::AmplitudeAdderCalibration::RunCalibration::doCalibrationFactor()


// -----------------------------------------------------------------------------
auto icarus::trigger::AmplitudeAdderCalibration::RunCalibration::findChannel
  (AdderChannelID channel) const -> ChannelCalibration_t const*
{
  auto const it = std::lower_bound(
    fCalibration.begin(), fCalibration.end(), channel,
    details::ChannelCalibration_t_sorter{}
    );
  return
    ((it != fCalibration.end()) && (it->channel == channel))? &*it: nullptr;
} // icarus::trigger::AmplitudeAdderCalibration::RunCalibration::findChannel()


// -----------------------------------------------------------------------------
auto icarus::trigger::AmplitudeAdderCalibration::RunCalibration::fetchChannel
  (AdderChannelID channel) const -> ChannelCalibration_t const&
{
  // find the calibration constants for the channel
  ChannelCalibration_t const* calib = findChannel(channel);
  if (!calib) throw UnknownChannelError{ channel };
  return *calib;
} // icarus::trigger::AmplitudeAdderCalibration::RunCalibration::fetchChannel()


// -----------------------------------------------------------------------------
icarus::trigger::AmplitudeAdderCalibration::AmplitudeAdderCalibration
  (Config const& config)
  : AmplitudeAdderCalibration{ config.CalibrationCoefficients() }
{}


// -----------------------------------------------------------------------------
icarus::trigger::AmplitudeAdderCalibration::AmplitudeAdderCalibration
  (CalibrationConstants_t calibration)
  : fCalibration{ std::move(calibration) }
{
  std::sort(
    fCalibration.begin(), fCalibration.end(),
    details::ChannelCalibration_t_sorter{}
    );
} // icarus::trigger::AmplitudeAdderCalibration::AmplitudeAdderCalibration()


// -----------------------------------------------------------------------------
auto icarus::trigger::AmplitudeAdderCalibration::doCalibrationForRun
  (unsigned int run) const
  -> std::unique_ptr<AdderCalibrationDatabase::RunCalibration>
{
  return std::make_unique<RunCalibration>(run, fCalibration);
}


// -----------------------------------------------------------------------------
void icarus::trigger::AmplitudeAdderCalibration::doDumpConfig
  (std::ostream& out, details::Indenter& nextLine) const
{
  // no not end the last line.
  out << nextLine << "calibration based on amplitude (AmplitudeAdderCalibration):"
    << nextLine << "linear calibration for " << fCalibration.size()
    << " adder channels (valid for all runs):"
    ;
  auto nextLineL1 = nextLine.nested("  ");
  nextLineL1.useLine(); // start from the next line
  for (ChannelCalibration_t const& calib: fCalibration) {
    out << nextLineL1 << "channel " << calib.channel
      << ": slope=(" << std::setw(12) << calib.slope()
      << " +/- " << std::setw(12) << calib.slopeErr()
      << "), offset=(" << std::setw(12) << calib.offset()
      << " +/- " << std::setw(12) << calib.offsetErr()
      << "); time delay=(" << std::setw(12) << calib.timeOffset
      << " +/- " << std::setw(12) << calib.timeOffsetErr << ")"
      ;
  }
  
} // doDumpConfig()


// -----------------------------------------------------------------------------
auto icarus::trigger::convert
  (AmplitudeAdderCalibration::ChannelCalibrationConfig const& config)
  -> AmplitudeAdderCalibration::ChannelCalibration_t
{
  return {
    /* channel       = */ config.channel(),
    /* slope         = */ config.slope(),
    /* offset        = */ config.offset(),
    /* slopeErr      = */ config.slopeErr(),
    /* offsetErr     = */ config.offsetErr(),
    /* covariance    = */ config.slopeOffsetCov(),
    /* timeOffset    = */ config.timeOffset(),
    /* timeOffsetErr = */ config.timeOffsetErr()
    };
} // namespace icarus::trigger


// -----------------------------------------------------------------------------
