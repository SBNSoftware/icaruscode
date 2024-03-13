/* Author: Matthew Strait <mstrait@fnal.gov> */

#ifndef artdaq_demo_Overlays_CRTFragment_hh
#define artdaq_demo_Overlays_CRTFragment_hh

//#include "Fragment.hh"
#include "artdaq-core/Data/Fragment.hh"
#include <ostream>

namespace CRT
{
  class Fragment;
}

class CRT::Fragment
{
public:

  struct header_t{
    uint8_t magic; // must be 'M'
    uint8_t nhit;
    uint16_t module_num;

    // The global ProtoDUNE-SP timestamp
    uint64_t fifty_mhz_time;

    // The raw timestamp held by the CRT's internal 32-bit counter.  You don't
    // want to look at this unless you are a CRT expert.  Should start at zero
    // at the start of the run, increment one-to-one with fifty_mhz_time, and
    // roll over every 86 seconds.
    uint32_t raw_backend_time;

    // Must tell GCC not to add padding, even at the cost of performance,
    // because this struct is describing exactly how the input data is
    // laid out.
  } __attribute__((packed));

  struct hit_t{
    uint8_t magic; // must be 'H'
    uint8_t channel;
    int16_t adc;

    // This time, disabling padding doesn't make a difference in struct
    // layout, but I am including it for consistency.
  } __attribute__((packed));

  // Return the module number for this fragment.  A CRT fragment consists
  // of a set of hits sharing a time stamp from one module.
  uint16_t module_num() const;

  // Return the number of hits claimed in the header of this fragment
  size_t num_hits() const;

  // Return the value of the 50MHz counter that holds the raw internal CRT time
  uint64_t raw_backend_time() const;

  // Return the value of the 50MHz counter that holds the global ProtoDUNE time
  uint64_t fifty_mhz_time() const;

  // Return the channel number of the ith hit.  That hit must exist.
  uint8_t channel(const int i) const;

  // Return the ADC value of the ith hit.  That hit must exist.
  int16_t adc(const int i) const;

  // Print the header to stdout, even if it is bad, but not if it
  // isn't all there.
  void print_header() const;

  // Print the given hit to stdout, even if it is bad, but not if it
  // isn't all there.
  void print_hit(const int i) const;

  // Print all the hits
  void print_hits() const;

  // Returns true if the header contains sensible values.  Otherwise,
  // prints a complaint and returns false.  Assumes header is complete
  // (check good_size() first).
  bool good_header() const;

  // Returns true if the hit contains sensible values.  Otherwise,
  // prints a complaint and returns false.  Assumes hit exists and
  // is complete (check good_size() first).
  bool good_hit(const int i) const;

  // Return the size of the CRT fragment in bytes.
  unsigned int size() const;

  // Returns true if the fragment is as big as the header says it is,
  // and false if it isn't, or it doesn't even have a full header.
  // i.e. if this is false, you're going to seg fault (or wish you had)
  // if you read the fragment.
  bool good_size() const;

  // Return true if the fragment contains a complete and sensible event.
  bool good_event() const;

  // Return a pointer to hit 'i'.  Not range checked.
  const hit_t * hit(const int i) const;

  // Return a pointer to the header
  const header_t * header() const;

  explicit Fragment(artdaq::Fragment const& f) : thefrag(f) {}

private:
  artdaq::Fragment const& thefrag;
};

#endif /* artdaq_demo_Overlays_CRTFragment_hh */
