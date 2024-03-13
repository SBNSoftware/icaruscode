/* Author: Matthew Strait <mstrait@fnal.gov> */

#include "CRTFragment.hh"

#include <ostream>

namespace CRT
{
  // Return the module number for this fragment.  A CRT fragment consists
  // of a set of hits sharing a time stamp from one module.
  uint16_t Fragment::module_num() const
  {
    return header()->module_num;
  }

  // Return the number of hits claimed in the header of this fragment
  size_t Fragment::num_hits() const
  {
    return header()->nhit;
  }

  // Return the value of the 50MHz counter that holds the raw internal CRT time
  uint64_t Fragment::raw_backend_time() const
  {
    return header()->raw_backend_time;
  }

  // Return the value of the 50MHz counter that holds the global ProtoDUNE time
  uint64_t Fragment::fifty_mhz_time() const
  {
    return header()->fifty_mhz_time;
  }

  // Return the channel number of the ith hit.  That hit must exist.
  uint8_t Fragment::channel(const int i) const
  {
    return hit(i)->channel;
  }

  // Return the ADC value of the ith hit.  That hit must exist.
  int16_t Fragment::adc(const int i) const
  {
    return hit(i)->adc;
  }

  // Print the header to stdout, even if it is bad, but not if it
  // isn't all there.
  void Fragment::print_header() const
  {
    if(size() < sizeof(header_t)){
      std::cerr << "CRT fragment smaller (" << size() << "B) than header (" 
                << sizeof(header_t) << "B), can't print\n";
      return;
    }

    std::cout << "CRT header: Magic = '" << header()->magic << "'\n"
              << "            n hit = " << +(header()->nhit) << "\n"
              << "            module = " << header()->module_num << "\n"
              << "            50Mhz time = " << header()->fifty_mhz_time 
              << " (0x" << std::hex << header()->fifty_mhz_time << std::dec << ")\n"
              << "            raw time   = " << header()->raw_backend_time << " (0x" 
              << std::hex << header()->raw_backend_time << std::dec << ")\n";
  }

  // Print the given hit to stdout, even if it is bad, but not if it
  // isn't all there.
  void Fragment::print_hit(const int i) const
  {
    // Is pointer arithmetic valid in the case that we're checking for?  Not
    // sure what the standard says, but I think that practically this should
    // work.
    if((uint8_t *)hit(i) + sizeof(hit_t) > thefrag.dataEndBytes()){
      std::cerr << "Hit " << i << " would be past end of fragment, can't print\n";
      return;
    }

    std::cout << "CRT hit " << i << ": Magic = '" << hit(i)->magic << "'\n"
              << "           channel = " << +(hit(i)->channel) << "\n"
              << "           ADC     = " << hit(i)->adc << "\n";
  }

  // Print all the hits
  void Fragment::print_hits() const
  {
    for(int i = 0; i < header()->nhit; i++)
      print_hit(i);
    puts("");
  }

  // Returns true if the header contains sensible values.  Otherwise,
  // prints a complaint and returns false.  Assumes header is complete
  // (check good_size() first).
  bool Fragment::good_header() const
  {
    const header_t * const h = header();
    if(h->magic != 'M'){
      std::cerr << "CRT header has wrong magic: " << h->magic << "\n";
      return false;
    }
    if(h->nhit == 0){
      std::cerr << "CRT event has no hits\n";
      return false;
    }
    if(h->nhit > 64){
      std::cerr << "CRT event has more hits (" << h->nhit << ") than channels (64)\n";
      return false;
    }
    return true;
  }

  // Returns true if the hit contains sensible values.  Otherwise,
  // prints a complaint and returns false.  Assumes hit exists and
  // is complete (check good_size() first).
  bool Fragment::good_hit(const int i) const
  {
    const hit_t * const h = hit(i);
    if(h->magic != 'H'){
      std::cerr << "CRT hit has wrong magic: " << h->magic << "\n";
      return false;
    }
    if(h->channel >= 64){
      std::cerr << "CRT hit has bad channel " << h->channel << " >= 64\n";
      return false;
    }
    if(h->adc >= 4096){
      // It is a 12-bit ADC.  This number probably represents the raw
      // ADC value before pedestal subtraction, but in any case, the
      // pedestal is positive, so the value still can't exceed 4095.
      std::cerr << "CRT hit has bad ADC value " << h->adc << " >= 4096\n";
      return false;
    }
    return true;
  }

  // Return the size of the CRT fragment in bytes.
  unsigned int Fragment::size() const
  {
    return thefrag.dataEndBytes() - thefrag.dataBeginBytes();
  }

  // Returns true if the fragment is as big as the header says it is,
  // and false if it isn't, or it doesn't even have a full header.
  // i.e. if this is false, you're going to seg fault (or wish you had)
  // if you read the fragment.
  bool Fragment::good_size() const
  {
    if(size() < sizeof(header_t)){
      std::cerr << "CRT fragment isn't as big (" << size() << "B) as header (" << sizeof(header_t) << "B)\n";
      return false;
    }

    const unsigned int expect_size =
      (sizeof(header_t) + header()->nhit * sizeof(hit_t)
       + sizeof(artdaq::RawDataType) - 1)
       /sizeof(artdaq::RawDataType)
       *sizeof(artdaq::RawDataType);

    if(size() != expect_size){
      std::cerr << "CRT fragment: N hit (" << header()->nhit << " -> " << expect_size << "B) mismatches size " 
                << size() << "B\n";
      for(char * c = (char *) thefrag.dataBeginBytes();
                 c < (char *)thefrag.dataEndBytes();
                 c++){
        std::cerr << (unsigned char)*c << "/" 
                  << (isprint((unsigned char)*c)?(unsigned char)*c:'.');
        if((c - (char*)thefrag.dataBeginBytes())%0x08 == 0x07)
          std::cerr << " ";
        if((c - (char*)thefrag.dataBeginBytes())%0x10 == 0x0f)
          std::cerr << "\n";
      }
      std::cerr << "\n";
      return false;
    }
    return true;
  }

  // Return true if the fragment contains a complete and sensible event.
  bool Fragment::good_event() const
  {
    if(!good_size()) return false;

    if(!good_header()) return false;

    for(unsigned int i = 0; i < header()->nhit; i++)
      if(!good_hit(i))
        return false;

    return true;
  }

  // Return a pointer to hit 'i'.  Not range checked.
  const CRT::Fragment::hit_t * Fragment::hit(const int i) const
  {
    return reinterpret_cast<const hit_t *>
      (thefrag.dataBeginBytes() + sizeof(header_t) + i*sizeof(hit_t));
  }

  // Return a pointer to the header
  const CRT::Fragment::header_t * Fragment::header() const
  {
    return reinterpret_cast<const header_t *>(thefrag.dataBeginBytes());
  }
}
