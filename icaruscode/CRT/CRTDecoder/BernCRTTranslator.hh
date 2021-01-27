#ifndef icaruscode_CRT_CRTDecoder_BERNCRTTranslator_hh
#define icaruscode_CRT_CRTDecoder_BERNCRTTranslator_hh

#include "art/Framework/Principal/Event.h"

//#include <iostream>
//#include <sstream>
//#include <string>

#include <vector>


namespace icarus::crt {

  class BernCRTTranslator;
//  std::ostream & operator << (std::ostream &, BernCRTTranslator const &);

}

class icarus::crt::BernCRTTranslator {

public:
  uint8_t   _mac5                    = 0;
  uint64_t  _run_start_time          = 0;
  uint64_t  _this_poll_start         = 0;
  uint64_t  _this_poll_end           = 0;
  uint64_t  _last_poll_start         = 0;
  uint64_t  _last_poll_end           = 0;
  int32_t   _system_clock_deviation  = 0; //system clock deviation w.r.t. steady clock, synchronised at the beginning of the run
  uint32_t  _hits_in_poll            = 0; //includes lost ones
  uint16_t  _hits_in_fragment        = 0; //size of the fragment hit array

  uint8_t  flags    = 0;
  uint16_t lostcpu  = 0;
  uint16_t lostfpga = 0;
  uint32_t ts0      = 0;
  uint32_t ts1      = 0;
  uint16_t adc[32]  = {0};
  uint32_t coinc    = 0;

  //Data added by fragment generator
  uint64_t  feb_hit_number          = 0; //hit counter for individual FEB, including hits lost in FEB or fragment generator
  uint64_t  timestamp               = 0; //absolute timestamp
  uint64_t  last_accepted_timestamp = 0; //timestamp of previous accepted hit
  uint16_t  lost_hits               = 0; //number of lost hits from the previous one
  
  
  static std::vector<BernCRTTranslator> getCRTData(art::Event const & evt);
//  static std::vector<int> getCRTData(art::Event const  & evt);


private:


};


#endif /* icaruscode_CRT_CRTDecoder_BERNCRTTranslator_hh */
