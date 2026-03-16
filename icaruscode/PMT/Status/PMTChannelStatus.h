/**
 * @file   icaruscode/PMT/Status/PMTChannelStatus.h
 * @brief  Interface for PMT channel status provider.
 * @author Matteo Vicenzi (mvicenzi@bnl.gov)
 */

#ifndef ICARUSCODE_PMT_STATUS_PMTCHANNELSTATUS_H
#define ICARUSCODE_PMT_STATUS_PMTCHANNELSTATUS_H

// LArSoft libraries
#include "larcorealg/CoreUtils/UncopiableAndUnmovableClass.h"

// Framework libraries
#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"

// C++ standard libraries
#include <set>
#include <string>


namespace icarusDB {

  /// Possible status values for a PMT channel, matching the database integers.
  enum PMTChannelStatusValue : int {
    kOFF = 0, ///< Channel is off (not powered).
    kON  = 1, ///< Channel is on and good.
    kBad = 2  ///< Channel is powered but flagged as bad.
  };


  /**
   * @brief Interface for PMT channel status information.
   *
   * Currently, the class provides interface for the following information:
   * - channel status: kON, kOFF, or kBad
   * - convenience predicates: isOn(), isOff(), isBad()
   * - channel sets by status: getOnChannels(), getOffChannels(), getBadChannels()
   *
   */
  class PMTChannelStatus : private lar::UncopiableAndUnmovableClass {

  public:

    using ChannelSet_t = std::set<unsigned int>;

    virtual ~PMTChannelStatus() noexcept = default;

    /// Returns the status of the specified channel.
    virtual PMTChannelStatusValue getChannelStatus(unsigned int channel) const = 0;

    /// Returns whether the specified channel is on and good.
    virtual bool isOn(unsigned int channel) const
      { return getChannelStatus(channel) == kON; }

    /// Returns whether the specified channel is off (not powered).
    virtual bool isOff(unsigned int channel) const
      { return getChannelStatus(channel) == kOFF; }

    /// Returns whether the specified channel is powered but bad.
    virtual bool isBad(unsigned int channel) const
      { return getChannelStatus(channel) == kBad; }

    /// Returns the set of channels with status kON.
    virtual ChannelSet_t getOnChannels()  const = 0;

    /// Returns the set of channels with status kOFF.
    virtual ChannelSet_t getOffChannels() const = 0;

    /// Returns the set of channels with status kBad.
    virtual ChannelSet_t getBadChannels() const = 0;

    /// Returns the nominal voltage [V] of the specified channel.
    virtual double getChannelVoltage(unsigned int channel) const = 0;

    /// Returns the database tag currently in use.
    virtual std::string getDatabaseTag() const = 0;

  }; // class PMTChannelStatus

} // namespace icarusDB


DECLARE_ART_SERVICE_INTERFACE(icarusDB::PMTChannelStatus, SHARED)


#endif // ICARUSCODE_PMT_STATUS_PMTCHANNELSTATUS_H
