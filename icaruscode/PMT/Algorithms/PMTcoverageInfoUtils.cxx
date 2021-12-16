/**
 * @file   icaruscode/PMT/Algorithms/PMTcoverageInfoUtils.cxx
 * @brief  Writes a collection of sbn::PMTcoverageInfo from PMT waveforms.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   November 22, 2021
 * @see    icaruscode/PMT/Algorithms/PMTcoverageInfoUtils.h
 */


// library header
#include "icaruscode/PMT/Algorithms/PMTcoverageInfoUtils.h"

// LArSoft libraries
#include "lardataalg/DetectorInfo/DetectorTimings.h"
#include "lardataobj/RawData/OpDetWaveform.h"


// -----------------------------------------------------------------------------
// ---  sbn::PMTcoverageInfoMaker
// -----------------------------------------------------------------------------
sbn::PMTcoverageInfoMaker::PMTcoverageInfoMaker
  (detinfo::DetectorTimings const& detTimings)
  : fOpDetTickPeriod{ detTimings.OpticalClockPeriod() }
  , fTriggerTime{ detTimings.TriggerTime() }
  , fBeamGateTime{ detTimings.BeamGateTime() }
  {}


// -----------------------------------------------------------------------------
sbn::PMTcoverageInfoMaker::PMTcoverageInfoMaker(microseconds opDetTickPeriod)
  : fOpDetTickPeriod{ opDetTickPeriod }
  {}


// -----------------------------------------------------------------------------
sbn::PMTcoverageInfo sbn::PMTcoverageInfoMaker::make
  (raw::OpDetWaveform const& waveform) const
{
  
  using detinfo::timescales::electronics_time;
  
  raw::Channel_t const channel = waveform.ChannelNumber();
  electronics_time const startTime { waveform.TimeStamp() };
  electronics_time const endTime
    = startTime + waveform.Waveform().size() * fOpDetTickPeriod;
  
  sbn::PMTcoverageInfo info {
      channel                   // channel
    , startTime.value()         // startTime
    , endTime.value()           // endTime
    /* the following are left default:
    // flags
    */
    };
  
  auto const setFlag = [&info]
    (sbn::PMTcoverageInfo::Flags_t::Flag_t flag, bool value)
    { if (value) info.flags.set(flag); else info.flags.unset(flag); };
  
  auto const isInWaveform = [startTime,endTime](electronics_time t)
    { return (t >= startTime) && (t < endTime); };
  
  if (fTriggerTime) {
    setFlag
      (sbn::PMTcoverageInfo::bits::WithTrigger, isInWaveform(*fTriggerTime));
  }
  
  if (fBeamGateTime) {
    setFlag
      (sbn::PMTcoverageInfo::bits::WithBeamGate, isInWaveform(*fBeamGateTime));
  }
  
  return info;
  
} // sbn::PMTcoverageInfoMaker::make()


// -----------------------------------------------------------------------------
// ---  functions
// -----------------------------------------------------------------------------
sbn::PMTcoverageInfo sbn::makePMTcoverageInfo(
  raw::OpDetWaveform const& waveform,
  detinfo::DetectorTimings const& detTimings
) {
  return sbn::PMTcoverageInfoMaker{ detTimings }.make(waveform);
}


sbn::PMTcoverageInfo sbn::makePMTcoverageInfo(
  raw::OpDetWaveform const& waveform,
  util::quantities::intervals::microseconds opDetTickPeriod
) {
  return sbn::PMTcoverageInfoMaker{ opDetTickPeriod }.make(waveform);
}


// -----------------------------------------------------------------------------
