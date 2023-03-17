#ifndef ICCRTDataRaw_hh_
#define ICCRTDataRaw_hh_

#include <cstdint>

namespace icarus::crt::CRTDataFlagsNew {
  static constexpr uint32_t TS0Present   = 0x0001; ///< TS0 counter was restarted in time
  static constexpr uint32_t TS1Present   = 0x0002; ///< TS1 counter was restarted in time
  static constexpr uint32_t TS0Reference = 0x0004; ///< This hit is a T0 reference signal
  static constexpr uint32_t TS1Reference = 0x0008; ///< This hit is a T1 reference signal
}

namespace icarus::crt {

   struct CRTDataRaw {

      uint8_t  fMac5    { 0 }; ///< FEB address following numbering convention common for both data and MC.
      uint32_t fEntry;         ///< Hit index for given FEB in an event (starts from 0 for each event).
      uint32_t fCounterT0 { 0 }; ///< Value of the T0 Counter.
      uint32_t fCounterT1 { 0 }; ///< Value of the T1 Counter.
      uint64_t fFullSecond { 0 }; ///> Full Second value of the Hit in Unix Server time [ns].
      uint64_t fTs0     { 0 }; ///< Absolute hit timestamp [ns].
      uint64_t fTs1     { 0 }; ///< Trigger time, not well defined as of Apr 14, 2021.
      uint16_t fAdc[64] {};    ///< ADC readout for each channel. CAEN (Bern) CRT FEBs use only indices 0-31.

      uint32_t fFlags    { 0 }; ///< DAQ uses 8 bits, remaining bits reserved for future use
      uint64_t fThisPollStart         { 0 }; ///< Absolute timestamp [ns] of the start of the data transfer from FEB
      uint64_t fLastPollStart         { 0 }; ///< Absolute timestamp [ns] of the start of the previous data transfer from FEB
      uint32_t fHitsInPoll            { 0 }; ///< Total number of hits in the data transfer (including hits omitted by FEB)

      uint32_t fCoinc                 { 0 }; ///< Top CRT coincidence bit map

      uint64_t fLastAcceptedTimestamp { 0 }; ///< Timestamp of the previous measured hit
      uint16_t fLostHits              { 0 }; ///< Number of hits omitted by FEB since the previous measured hit


      bool IsOverflow_TS0()  const { return !(fFlags & CRTDataFlagsNew::TS0Present); }
      bool IsOverflow_TS1()  const { return !(fFlags & CRTDataFlagsNew::TS1Present); }
      bool IsReference_TS0() const { return   fFlags & CRTDataFlagsNew::TS0Reference; }
      bool IsReference_TS1() const { return   fFlags & CRTDataFlagsNew::TS1Reference; }
      
  };

} // namespace icarus::crt


#endif


