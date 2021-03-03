/**
 * This class provides a tool to read out CRT data.
 * The main method, getCRTData reads an event and
 * returns a vector containing all CRT hits in a
 * unified format.
 * This class automatically recognizes data
 * collected in pull and push modes, and with any
 * version of raw fragment type.
 */

#ifndef icaruscode_CRT_CRTDecoder_BERNCRTTranslator_hh
#define icaruscode_CRT_CRTDecoder_BERNCRTTranslator_hh

#include "art/Framework/Principal/Event.h"

#include <vector>


namespace icarus::crt {

  class BernCRTTranslator;
  std::ostream & operator << (std::ostream &, BernCRTTranslator const &);

}

class icarus::crt::BernCRTTranslator {

public:
  uint8_t   mac5                    = 0; //last 8 bits of FEB MAC5 address
  uint64_t  run_start_time          = 0; //absolute timestamp [ns]
  uint64_t  this_poll_start         = 0; //absolute timestamp [ns]
  uint64_t  this_poll_end           = 0; //absolute timestamp [ns]
  uint64_t  last_poll_start         = 0; //absolute timestamp [ns]
  uint64_t  last_poll_end           = 0; //absolute timestamp [ns]
  int32_t   system_clock_deviation  = 0; //system clock deviation [ns] w.r.t. steady clock, synchronised at the beginning of the run
  uint32_t  hits_in_poll            = 0; //includes lost hits
  uint16_t  hits_in_fragment        = 0; //size of the fragment hit array

  uint8_t  flags                    = 0; //see corresponding methods
  uint16_t lostcpu                  = 0; //hits lost before the trigger hit due to too long polling time
  uint16_t lostfpga                 = 0; //hits lost after the trigger hit due to dead time
  uint32_t ts0                      = 0; //number of [ns] since last T0 signal
  uint32_t ts1                      = 0; //number of [ns] since last T1 signal
  uint16_t adc[32]                  = {0}; //ADC readout for each of 32 FEB channels
  uint32_t coinc                    = 0; //used by TOP CRT

  //Data added by fragment generator
  uint64_t  feb_hit_number          = 0; //hit counter for individual FEB, including hits lost in FEB or fragment generator
  uint64_t  timestamp               = 0; //absolute hit timestamp
  uint64_t  sequence_id             = 0; //*fragment* sequence_id. Note a fragment may contain multiple hits
  uint64_t  last_accepted_timestamp = 0; //timestamp of previous accepted hit
  uint16_t  lost_hits               = 0; //number of lost hits from the previous one

  bool IsOverflow_TS0()  const { return !(flags&1); } //TODO double check the logic
  bool IsOverflow_TS1()  const { return !(flags&2); }
  bool IsReference_TS0() const { return   flags&4; }
  bool IsReference_TS1() const { return   flags&8; }
  
  static std::vector<BernCRTTranslator> getCRTData(art::Event const & evt);


private:
  static BernCRTTranslator analyze_BernCRTZMQFragment(artdaq::Fragment & frag); 
  static BernCRTTranslator analyze_BernCRTFragment(artdaq::Fragment & frag); 
  static std::vector<BernCRTTranslator> analyze_BernCRTFragmentV2(artdaq::Fragment & frag); 

};


#endif /* icaruscode_CRT_CRTDecoder_BERNCRTTranslator_hh */
