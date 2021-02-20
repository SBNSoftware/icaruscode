/**
 * @file   icaruscode/Decode/DecoderTools/PMTconfigurationExtractor.h
 * @brief  Utility to extract PMT readout configuration from data.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   February 18, 2021
 * @file   icaruscode/Decode/DecoderTools/PMTconfigurationExtractor.cxx
 */

#ifndef ICARUSCODE_DECODE_DECODERTOOLS_PMTCONFIGURATIONEXTRACTOR_H
#define ICARUSCODE_DECODE_DECODERTOOLS_PMTCONFIGURATIONEXTRACTOR_H


// ICARUS libraries
#include "sbnobj/ICARUS/PMT/Data/PMTconfiguration.h"
#include "sbnobj/ICARUS/PMT/Data/V1730Configuration.h"

// framework libraries
#include "fhiclcpp/ParameterSet.h"

// C/C++ standard libraries
#include <string>


// -----------------------------------------------------------------------------
namespace icarus { class PMTconfigurationExtractor; }
/**
 * @brief Class to extract PMT readout board configuration.
 * 
 * This is an example of PMT readout board configuration taken from ICARUS run
 * 4774 (`destinations` and `metrics` have been omitted):
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * icaruspmteetop01: {
 *   daq: {
 *     fragment_receiver: {
 *       BaselineCh1: 14768
 *       BaselineCh10: 14843
 *       BaselineCh11: 14846
 *       BaselineCh12: 14840
 *       BaselineCh13: 14981
 *       BaselineCh14: 14896
 *       BaselineCh15: 14966
 *       BaselineCh16: 14800
 *       BaselineCh2: 14836
 *       BaselineCh3: 14806
 *       BaselineCh4: 14954
 *       BaselineCh5: 14766
 *       BaselineCh6: 14835
 *       BaselineCh7: 14861
 *       BaselineCh8: 14857
 *       BaselineCh9: 14949
 *       BoardChainNumber: 0
 *       CalibrateOnConfig: true
 *       ChargePedstalBitCh1: 14400
 *       CircularBufferSize: 5e8
 *       CombineReadoutWindows: false
 *       GetNextFragmentBunchSize: 20
 *       GetNextSleep: 100
 *       InterruptEventNumber: 1
 *       InterruptLevel: 0
 *       InterruptStatusID: 0
 *       LVDSLogicValueG1: 3
 *       LVDSLogicValueG2: 3
 *       LVDSLogicValueG3: 3
 *       LVDSLogicValueG4: 3
 *       LVDSLogicValueG5: 3
 *       LVDSLogicValueG6: 3
 *       LVDSLogicValueG7: 3
 *       LVDSLogicValueG8: 3
 *       LVDSOutWidthC1: 20
 *       LVDSOutWidthC10: 20
 *       LVDSOutWidthC11: 20
 *       LVDSOutWidthC12: 20
 *       LVDSOutWidthC13: 20
 *       LVDSOutWidthC14: 20
 *       LVDSOutWidthC15: 20
 *       LVDSOutWidthC16: 20
 *       LVDSOutWidthC2: 20
 *       LVDSOutWidthC3: 20
 *       LVDSOutWidthC4: 20
 *       LVDSOutWidthC5: 20
 *       LVDSOutWidthC6: 20
 *       LVDSOutWidthC7: 20
 *       LVDSOutWidthC8: 20
 *       LVDSOutWidthC9: 20
 *       LockTempCalibration: false
 *       MajorityCoincidenceWindow: 1
 *       MajorityLevel: 0
 *       MajorityTimeWindow: 4
 *       ModeLVDS: 1
 *       SWTrigger: false
 *       SelfTrigBit: 16
 *       SelfTriggerMask: 255
 *       SelfTriggerMode: 0
 *       TimeOffsetNanoSec: 0
 *       UseTimeTagForTimeStamp: false
 *       Verbosity: 1
 *       acqMode: 0
 *       allowTriggerOverlap: false
 *       analogMode: 1
 *       boardId: 0
 *       board_id: 18
 *       channelEnable0: true
 *       channelEnable1: true
 *       channelEnable10: true
 *       channelEnable11: true
 *       channelEnable12: true
 *       channelEnable13: true
 *       channelEnable14: true
 *       channelEnable15: true
 *       channelEnable2: true
 *       channelEnable3: true
 *       channelEnable4: true
 *       channelEnable5: true
 *       channelEnable6: true
 *       channelEnable7: true
 *       channelEnable8: true
 *       channelEnable9: true
 *       channelPedestal0: 6554
 *       channelPedestal1: 6554
 *       channelPedestal10: 6554
 *       channelPedestal11: 6554
 *       channelPedestal12: 6554
 *       channelPedestal13: 6554
 *       channelPedestal14: 6554
 *       channelPedestal15: 6554
 *       channelPedestal2: 6554
 *       channelPedestal3: 6554
 *       channelPedestal4: 6554
 *       channelPedestal5: 6554
 *       channelPedestal6: 6554
 *       channelPedestal7: 6554
 *       channelPedestal8: 6554
 *       channelPedestal9: 6554
 *       dacValue: 32768
 *       data_buffer_depth_fragments: 10000
 *       data_buffer_depth_mb: 2000
 *       debugLevel: 7
 *       destinations: {} # ...
 *       dynamicRange: 0
 *       enableReadout: 1
 *       eventCounterWarning: 1
 *       eventsPerInterrupt: 1
 *       extTrgMode: 3
 *       fragment_id: 18
 *       fragment_type: "CAENV1730"
 *       generator: "CAENV1730Readout"
 *       ioLevel: 1
 *       irqWaitTime: 1
 *       link: 2
 *       maxEventsPerTransfer: 1
 *       max_fragment_size_bytes: 1e7
 *       memoryAlmostFull: 2
 *       multicast_interface_ip: "192.168.191.0"
 *       nChannels: 16
 *       outputSignalMode: 0
 *       postPercent: 70
 *       readoutMode: 0
 *       receive_requests: true
 *       recordLength: 25000
 *       request_address: "227.128.1.129"
 *       request_mode: "sequence"
 *       request_port: 3502
 *       routing_table_config: { use_routing_master: false }
 *       runSyncMode: 0
 *       separate_data_thread: true
 *       swTrgMode: 0
 *       testPattern: 0
 *       triggerPolarity: 0
 *       triggerPulseWidth: 20
 *       triggerThreshold0: 14368
 *       triggerThreshold1: 14436
 *       triggerThreshold10: 14446
 *       triggerThreshold11: 14440
 *       triggerThreshold12: 14581
 *       triggerThreshold13: 14596
 *       triggerThreshold14: 14566
 *       triggerThreshold15: 14400
 *       triggerThreshold2: 14406
 *       triggerThreshold3: 14554
 *       triggerThreshold4: 14366
 *       triggerThreshold5: 14435
 *       triggerThreshold6: 14461
 *       triggerThreshold7: 14457
 *       triggerThreshold8: 14549
 *       triggerThreshold9: 14443
 *       usePedestals: false
 *     }
 *     metrics: {} # ...
 *   }
 * }
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * 
 * 
 * 
 * 
 */
class icarus::PMTconfigurationExtractor {
  
    public:
  
  /**
   * @brief Extracts all supported PMT configuration from `config`.
   * @param config a FHiCL parameter set with component configuration
   * @return an object with the supported PMT configuration
   * 
   * All PMT-related configuration that is known to this code is extracted and
   * returned.
   */
  icarus::PMTconfiguration extract(fhicl::ParameterSet const& config) const;
  
    private:
  
  /**
   * @brief Extracts PMT readout board configuration from `pset`.
   * @param pset information source (FHiCL configuration)
   * @param boardName name of the board we are given the configuration of
   * @return the requested configuration
   * 
   * 
   */
  icarus::V1730Configuration extractV1730configuration
    (fhicl::ParameterSet const& pset, std::string const& boardName)
    const;
  
  /**
   * @brief Returns the specified V1730 readout channel configuration.
   * @param boardPSet readout board configuration
   * @param channelNo number of channel on board (0-15)
   * @return the configuration of the specified channel
   */
  icarus::V1730channelConfiguration extractChannelConfiguration
    (fhicl::ParameterSet const& boardPSet, unsigned short int channelNo) const;
  
  /**
   * @brief Returns the specified PMT readout board configuration.
   * @param pset parameter set including `key`
   * @param key key of the PMT readout configuration candidate
   * @return the configuration, or an empty object if key does not represent one
   */
  std::optional<fhicl::ParameterSet> readBoardConfig
    (fhicl::ParameterSet const& pset, std::string const& key) const;
  
}; // icarus::PMTconfigurationExtractor


// ---------------------------------------------------------------------------


#endif // ICARUSCODE_DECODE_DECODERTOOLS_PMTCONFIGURATIONEXTRACTOR_H
