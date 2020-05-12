#ifndef ICCRTData_hh_
#define ICCRTData_hh_

namespace icarus {
namespace crt {

  struct CRTData {

      uint8_t  fMac5;
      uint32_t fEntry;
      uint32_t fTs0;
      uint32_t fTs1;
      uint16_t fAdc[32];

      CRTData() {}
  };

} // namespace crt
} // namespace icarus

#endif
