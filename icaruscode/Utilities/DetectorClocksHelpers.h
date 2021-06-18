/**
 * @file   icaruscode/Utilities/DetectorClocksHelpers.h
 * @brief  Simple functions to streamline the creation of `DetectorClocksData`.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   September 16, 2020
 */

#ifndef ICARUSCODE_UTILITIES_DETECTORCLOCKSHELPERS_H
#define ICARUSCODE_UTILITIES_DETECTORCLOCKSHELPERS_H

// LArSoft libraries
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataalg/DetectorInfo/DetectorTimings.h"
#include "lardataalg/DetectorInfo/DetectorClocksData.h"

// framework libraries
#include "art/Framework/Services/Registry/ServiceHandle.h"


// -----------------------------------------------------------------------------
// forward declarations
namespace art { class Event; }

// -----------------------------------------------------------------------------
namespace icarus::ns::util {
  
  // --- BEGIN -- DetectorClocksData helpers -----------------------------------
  /// @name DetectorClocksData helpers
  /// @{
  /**
   * @brief Returns a `detinfo::DetectorClocksData` from
   *        `DetectorClocksService`.
   * @param event pointer to an _art_ event
   * @return a `detinfo::DetectorClocksData` object (see the details below)
   * 
   * The `detinfo::DetectorClocksService` is queried to obtain a
   * `detinfo::DetectorClocksData` object. The returned object is a static copy
   * whose values will never update, even if the `DetectorClocksService` state
   * changes.
   * 
   * If the `event` pointer is valid (i.e. not `nullptr`), the event is passed
   * to the service `DataFor()` method for a complete timing record.
   * Otherwise, `DataForJob()` is used, and the information might be incomplete
   * or become outdated by the time a new event, run or file is accessed.
   * 
   * @note This function is _art_ dependent and accesses `DetectorClocksService`
   *       service.
   */
  detinfo::DetectorClocksData makeDetClockData(art::Event const* event)
    {
      auto const& detClocks
        = *(art::ServiceHandle<detinfo::DetectorClocksService const>());
      return event? detClocks.DataFor(*event): detClocks.DataForJob();
    } // makeDetClockData()

  /// Returns detector clock data for the specified event.
  /// @see `makeDetClockData(art::Event const*)`
  detinfo::DetectorClocksData makeDetClockData(art::Event const& event)
    { return makeDetClockData(&event); }
  
  /// Returns generic detector clock data for the job.
  /// @see `makeDetClockData(art::Event const*)`
  detinfo::DetectorClocksData makeDetClockData()
    { return makeDetClockData(nullptr); }
  
  /// @}
  // --- END -- DetectorClocksData helpers -------------------------------------

  
  // --- BEGIN -- DetectorTimings helpers --------------------------------------
  /// @name DetectorTimings helpers
  /// @{
  
  /**
   * @brief Returns a `detinfo::DetectorTimings` from `DetectorClocksService`.
   * @param event pointer to an _art_ event
   * @return a `detinfo::DetectorTimings` object (see the details below)
   * @see `makeDetClockData(art::Event const*)`
   * 
   * A `detinfo::DetectorTimings` object is created out of information obtained
   * from `detinfo::DetectorClocksService`.
   *
   * All the considerations described in `makeDetClockData(art::Event const*)`
   * also hold for this function.
   */
  detinfo::DetectorTimings makeDetTimings(art::Event const* event)
    { return detinfo::DetectorTimings{ makeDetClockData(event) }; }
  
  /// Returns detector clock data for the specified event.
  /// @see `makeDetTimings(art::Event const*)`
  detinfo::DetectorTimings makeDetTimings(art::Event const& event)
    { return makeDetTimings(&event); }
  
  /// Returns generic detector clock data for the job.
  /// @see `makeDetTimings(art::Event const*)`
  detinfo::DetectorTimings makeDetTimings()
    { return makeDetTimings(nullptr); }
  
  
  /// @}
  // --- END -- DetectorTimings helpers ----------------------------------------

  
} // namespace icarus::ns::util


#endif // ICARUSCODE_UTILITIES_DETECTORCLOCKSHELPERS_H
