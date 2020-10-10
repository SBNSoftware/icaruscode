#ifndef ICCRTData_hh_
#define ICCRTData_hh_

#include <cstdint>

namespace icarus {
namespace crt {

  struct CRTData {

      uint8_t  fMac5;
      uint32_t fEntry;
      uint64_t fTs0;
      uint64_t fTs1;
      uint16_t fAdc[64]; //only use indices 0-31 for C or M modules

      CRTData() {}
  };

} // namespace crt
} // namespace icarus

#endif
