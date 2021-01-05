/**
 * @file   icaruscode/Decode/Dumpers/PMTdataDumpers.h
 * @brief  Functions to dump the content of PMT (V1730) data fragments to console
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   December 22, 2020
 * @see    icaruscode/Decode/Dumpers/PMTdataDumpers.cxx
 */

#ifndef ICARUSCODE_DECODE_DUMPERS_PMTDATADUMPERS_H
#define ICARUSCODE_DECODE_DUMPERS_PMTDATADUMPERS_H


// SBN DAQ libraries
#include "sbndaq-artdaq-core/Overlays/Common/CAENV1730Fragment.hh"

// C++ standard libraries
#include <iosfwd>


namespace sbndaq {
  
  std::ostream& operator<< (std::ostream& out, sbndaq::CAENV1730FragmentMetadata const& meta);
  std::ostream& operator<< (std::ostream& out, sbndaq::CAENV1730EventHeader const& header);
  std::ostream& operator<< (std::ostream& out, sbndaq::CAENV1730Event const& event);
  std::ostream& operator<< (std::ostream& out, sbndaq::CAENV1730Fragment const& frag);
  
} // namespace sbndaq


#endif // ICARUSCODE_DECODE_DUMPERS_PMTDATADUMPERS_H
