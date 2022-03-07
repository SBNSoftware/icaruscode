/**
 * @file   icaruscode/PMT/Trigger/Algorithms/WindowTopologyAlg.h
 * @brief  Assembles the topology of trigger windows.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   March 25, 2021
 * @see    icaruscode/PMT/Trigger/Algorithms/WindowTopologyAlg.cxx
 */


#ifndef ICARUSCODE_PMT_TRIGGER_ALGORITHMS_WINDOWTOPOLOGYALG_H
#define ICARUSCODE_PMT_TRIGGER_ALGORITHMS_WINDOWTOPOLOGYALG_H


// ICARUS libraries
#include "icaruscode/PMT/Trigger/Algorithms/WindowChannelMap.h"
#include "icaruscode/PMT/Trigger/Utilities/TrackedOpticalTriggerGate.h"
#include "icaruscode/IcarusObj/OpDetWaveformMeta.h" // sbn::OpDetWaveformMeta
#include "icarusalg/Utilities/mfLoggingClass.h"

// framework libraries
#include "cetlib_except/exception.h"

// C/C++ standard libraries
#include <optional>
#include <vector>
#include <string>


// -----------------------------------------------------------------------------
// forward declaration
namespace geo {
  class GeometryCore;
  class CryostatGeo;
}

// -----------------------------------------------------------------------------
namespace icarus::trigger {
  class WindowTopologyAlg;
  class WindowTopologyVerification;
  class WindowTopologyManager;
}
// -----------------------------------------------------------------------------
/**
 * @brief Algorithm to create trigger window topology information.
 * 
 * The algorithm defines the windows based on trigger gates and their channels,
 * and establishes their topology cryostat by cryostat.
 * 
 */
class icarus::trigger::WindowTopologyAlg
  : protected icarus::ns::util::mfLoggingClass
{
  
    public:
  
  /// Type of trigger gate extracted from the input event.
  using InputTriggerGate_t
    = icarus::trigger::TrackedOpticalTriggerGate<sbn::OpDetWaveformMeta>;
  
  /// A list of trigger gates from input.
  using TriggerGates_t = std::vector<InputTriggerGate_t>;

  /// Type of lists of gates, one list per cryostat (outer index: cryostat no).
  using TriggerGatesPerCryostat_t = std::vector<TriggerGates_t>;

  
  /**
   * @brief Constructor.
   * @param geom geometry service provider
   * @param logCategory category tag for messages from this algorithm
   */
  WindowTopologyAlg(
    geo::GeometryCore const& geom,
    std::string const& logCategory = "WindowTopologyAlg"
    );
  
  /**
   * @brief Returns the topology of the windows described by the gates.
   * @param gates a set of gates, grouped by cryostat
   * @return topology of the windows described by the gates
   * 
   * The input `gates` is a collection with as many entries as there are
   * cryostats in the detector, the first entry collecting all gates from the
   * first cryostat (`C:0`), the next all the gates from the second one (`C:2`),
   * and so on.
   * 
   * Each gate covers one or more channels: each window is made out of all the
   * channels in a gate. The opening content of the gates is not considered.
   */
  WindowChannelMap createFromGates
    (TriggerGatesPerCryostat_t const& gates) const;
  
  
  /**
   * @brief Returns the topology of the windows described by the gates.
   * @param gates a set of gates
   * @return topology of the windows described by the gates
   * 
   * The input `gates` is a collection with all relevant gates in the detector.
   * 
   * Each gate covers one or more channels: each window is made out of all the
   * channels in a gate. The opening content of the gates is not considered.
   * 
   * Gates are first split by the cryostat they belong to (cross-cryostat gates
   * are forbidden), and then the window topology is extracted.
   */
  WindowChannelMap createFromGates(TriggerGates_t const& gates) const;
  
  
  /**
   * @brief Returns the topology of the windows described by the gates.
   * @param gates a set of gates within a single cryostat
   * @param cryo the cryostat the gates (and their channels) belong to
   * @return topology of the windows described by the gates
   * 
   * The input `gates` is a collection with all the gates from the specified
   * cryostat `cryo`.
   * 
   * Each gate covers one or more channels: each window is made out of all the
   * channels in a gate. The opening content of the gates is not considered.
   */
  WindowChannelMap createFromCryostatGates
    (TriggerGates_t const& windowChannels, geo::CryostatGeo const& cryo) const;


  /**
   * @brief Returns the topology of the windows described by the gates.
   * @param gates a set of gates within a single cryostat
   * @param cryoID ID of the cryostat the gates (and their channels) belong to
   * @return topology of the windows described by the gates
   * @see `createFromCryostatGates(TriggerGates_t const&, geo::CryostatGeo const&)`
   * 
   * The input `gates` is a collection with all the gates from the specified
   * cryostat with ID `cryoID`.
   * 
   * See `createFromCryostatGates(WindowChannelsPerCryostat_t const&, geo::CryostatGeo const&)`
   * for more information.
   */
  WindowChannelMap createFromCryostatGates
    (TriggerGates_t const& windowChannels, geo::CryostatID cryoID) const;
  
  
    private:
  /// All channels in a gate.
  using WindowChannels_t = std::vector<raw::Channel_t>;
  
  /// All channels in many gates (one list per gate).
  using WindowChannelColl_t = std::vector<WindowChannels_t>;

  
  geo::GeometryCore const* const fGeom; ///< Geometry service provider.

  
  /// Convenience function: creates and returns a `WindowChannelMap` from the
  /// arguments (`args`); in the meanwhile, it also dumps it into debug stream.
  template <typename... Args>
  WindowChannelMap emplaceAndDumpMap(Args&&... args) const;
  
  
  /**
   * @brief Extracts topology information from a set of windows.
   * @param windowChannels the windows, specified by composition in channel ID
   * @param cryo the cryostat the windows belong to
   * @param geom geometry service provider
   * @param firstWindowIndex (default: `0`) the index assigned to the first
   *                         window
   * @return window topology information, one entry per window
   */
  static std::vector<WindowChannelMap::WindowInfo_t> createWindowsFromCryostat(
    WindowChannelColl_t const& windowChannels,
    geo::CryostatGeo const& cryo,
    geo::GeometryCore const& geom,
    std::size_t firstWindowIndex = 0U
    );

  /**
   * @brief Extracts the channel ID from a collection of gates.
   * @param gates a collection of gates
   * @return a collection of channel lists
   * 
   * The returned collection has one entry per input gate in `gates`: each entry
   * is the list of channels that the corresponding gate covers.
   */
  static WindowChannelColl_t extractGateChannels(TriggerGates_t const& gates);

  /**
   * @brief Returns the window in `windowList` closest to the `targetWindow`.
   * @param windowList list of (pointers to) the candidate windows
   * @param targetWindow the reference window
   * @return a pointer to the closest among `windowList`, `nullptr` if none
   * 
   * A pointer to the window with the smallest distance from `targetWindow` is
   * returned. The distance is computed between the centers of the windows as
   * reported by the windows themselves.
   * 
   * The only case where no window (`nullptr`) is returned is when the window
   * candidate list `windowList` is empty.
   */
  static WindowChannelMap::WindowInfo_t const* findClosestWindow(
    std::vector<WindowChannelMap::WindowInfo_t*> const& windowList,
    WindowChannelMap::WindowInfo_t const* targetWindow
    );
  
}; // icarus::trigger::WindowTopologyAlg


//------------------------------------------------------------------------------
/**
 * @brief Algorithm verifying the topology from trigger gates.
 * 
 * The algorithm needs to be provided a reference topology information.
 * It will then analyze any set of trigger gates to confirm that their grouping
 * matches the one defined in that topology.
 */
class icarus::trigger::WindowTopologyVerification
  : protected icarus::ns::util::mfLoggingClass
{
  
  /// The reference topology to check against.
  std::optional<icarus::trigger::WindowChannelMap> fWindowMap;
  
    public:
  
  /// Type of sets of trigger gates.
  using TriggerGates_t = icarus::trigger::WindowTopologyAlg::TriggerGates_t;
  
  /// Type of sets of trigger gates, grouped by cryostat.
  using TriggerGatesPerCryostat_t
    = icarus::trigger::WindowTopologyAlg::TriggerGatesPerCryostat_t;
  
  
  /**
   * @brief Constructor.
   * @param topology the window topology to be verified against
   * @param logCategory category tag for messages from this algorithm
   */
  WindowTopologyVerification(
    icarus::trigger::WindowChannelMap topology,
    std::string const& logCategory = "WindowTopologyVerification"
    );
  
  /**
   * @brief Constructor.
   * @param logCategory category tag for messages from this algorithm
   * 
   * Topology must be set up with `setTopology()` before any verification can
   * occur.
   */
  WindowTopologyVerification
    (std::string const& logCategory = "WindowTopologyVerification");
  
  
  /// Returns whether a window topology is set (see `setTopology()`).
  bool hasTopology() const;
  
  /// Resets the topology to be verified against.
  void setTopology(icarus::trigger::WindowChannelMap topology);
  
  /// Forgets the current topology. Set a new one with `setTopology()`.
  void clearTopology();
  
  /// Returns the current window topology.
  /// @throw std::bad_optional_access if no topology has been set
  WindowChannelMap const& getTopology() const;
  
  //@{
  /**
   * @brief Verifies that `gates` match the topology current set up.
   * @param gates the set of gates to match to the topology
   * @return a composite error message (empty on success)
   * @throw cet::exception (category: `WindowTopologyVerification`) if no window
   *                       topology is set up yet
   * 
   * Gates can be specified either as a collection, or as a set of collections,
   * one per cryostat.
   */
  std::string verify(TriggerGatesPerCryostat_t const& gates) const;
  std::string verify(TriggerGates_t const& gates) const;
  //@}
  
  //@{
  /**
   * @brief Verifies that `gates` match the topology current set up.
   * @param gates the set of gates to match to the topology
   * @throw cet::exception (category: `WindowTopologyVerification`) if no window
   *                       topology is set up yet
   * @throw cet::exception (category: `WindowTopologyVerification`) if
   *                       verification fails
   * @see `verify()`
   */
  void operator() (TriggerGatesPerCryostat_t const& gates) const;
  void operator() (TriggerGates_t const& gates) const;
  //@}
  
  
    private:
  
  /// Type of trigger gate used for input to the algorithm.
  using InputTriggerGate_t
    = icarus::trigger::WindowTopologyAlg::InputTriggerGate_t;
  
  /// Checks that the specified `gate` matched the window with index `iWindow`.
  /// @return message describing the failure, empty on success
  std::string verifyGate
    (std::size_t iWindow, InputTriggerGate_t const& gate) const;
  
  /// Implementation of `operator()`.
  template <typename Gates>
  void verifyOrThrow(Gates const& gates) const;
  
}; // icarus::trigger::WindowTopologyVerification


//------------------------------------------------------------------------------
/**
 * @brief Class to extract and verify a window topology from trigger gates.
 * 
 * This class is meant to provide a convenient interface for the creation
 * (`WindowTopologyAlg`) and following verification
 * (`WindowTopologyVerification`) of a window topology.
 * The first time the class is given a trigger gate set, it extracts the
 * topology from it. The following times, it will verify that the gates match
 * that topology, and throw an exception otherwise.
 * 
 * The class can be reset to restart the cycle, and the topology can be obtained
 * directly with the `operator->`.
 * 
 * The main methods support different types of input, collectively known as
 * `Gates`. The complete set of supported input depends on what is supported
 * by the underlying classes. Currently, supported input includes (flat)
 * collections of gates and collections of gates grouped by cryostat.
 */
class icarus::trigger::WindowTopologyManager
  : protected icarus::ns::util::mfLoggingClass
{
  
  /// Verification algorithm; holds the current topology information.
  WindowTopologyVerification fVerify;
  
  geo::GeometryCore const* const fGeom; ///< Geometry service provider.
  
    public:
  
  /// Type of sets of trigger gates.
  using TriggerGates_t = icarus::trigger::WindowTopologyAlg::TriggerGates_t;
  
  /// Type of sets of trigger gates, grouped by cryostat.
  using TriggerGatesPerCryostat_t
    = icarus::trigger::WindowTopologyAlg::TriggerGatesPerCryostat_t;
  
  
  /**
   * @brief Constructor.
   * @param topology the window topology to be verified against
   * @param logCategory category tag for messages from this algorithm
   */
  WindowTopologyManager(
    geo::GeometryCore const& geom,
    std::string const& logCategory = "WindowTopologyManager"
    );
  
  /**
   * @brief Extracts topology or verifies it against `gates` .
   * @tparam Gates type of gate collection
   * @param gates the set of gates to match to the topology
   * @return a composite error message (empty on success)
   * 
   */
  template <typename Gates>
  std::string setOrVerify(Gates const& gates);

  /**
   * @brief Extracts topology or verifies it against `gates` .
   * @tparam Gates type of gate collection
   * @param gates the set of gates to match to the topology
   * @return whether this method extracted the topology
   * @throw cet::exception (category: `WindowTopologyVerification`) if
   *                       verification fails
   */
  template <typename Gates>
  bool operator() (Gates const& gates);

  /// Access to the stored topology.
  /// @throw std::bad_optional_access if no topology has been set
  WindowChannelMap const& operator* () const;
  
  /// Access to the stored topology.
  /// @throw std::bad_optional_access if no topology has been set
  WindowChannelMap const* operator-> () const;
  
  
  
    private:
  
  /// Helper: creates the topology from `gates`.
  template <typename Gates>
  void extractTopology(Gates const& gates);
  
}; // icarus::trigger::WindowTopologyManager


// -----------------------------------------------------------------------------
// ---  inline implementation
// -----------------------------------------------------------------------------
// --- icarus::trigger::WindowTopologyVerification
// -----------------------------------------------------------------------------
inline icarus::trigger::WindowTopologyVerification::WindowTopologyVerification(
  icarus::trigger::WindowChannelMap topology,
  std::string const& logCategory /* = "WindowTopologyVerification" */
  )
  : icarus::ns::util::mfLoggingClass(logCategory)
  , fWindowMap(std::move(topology))
{}


// -----------------------------------------------------------------------------
inline icarus::trigger::WindowTopologyVerification::WindowTopologyVerification(
  std::string const& logCategory /* = "WindowTopologyVerification" */
  )
  : icarus::ns::util::mfLoggingClass(logCategory)
{}


// -----------------------------------------------------------------------------
inline bool icarus::trigger::WindowTopologyVerification::hasTopology() const
  { return fWindowMap.has_value(); }


// -----------------------------------------------------------------------------
inline void icarus::trigger::WindowTopologyVerification::setTopology
  (icarus::trigger::WindowChannelMap topology)
  { fWindowMap = std::move(topology); }


//------------------------------------------------------------------------------
inline auto icarus::trigger::WindowTopologyVerification::getTopology() const
  -> WindowChannelMap const&
  { return fWindowMap.value(); }


//------------------------------------------------------------------------------
inline void icarus::trigger::WindowTopologyVerification::clearTopology()
  { return fWindowMap.reset(); }


//------------------------------------------------------------------------------
inline void icarus::trigger::WindowTopologyVerification::operator()
  (TriggerGates_t const& gates) const
  { verifyOrThrow(gates); }


//------------------------------------------------------------------------------
inline void icarus::trigger::WindowTopologyVerification::operator()
  (TriggerGatesPerCryostat_t const& gates) const
  { verifyOrThrow(gates); }


//------------------------------------------------------------------------------
//--- icarus::trigger::WindowTopologyManager
//------------------------------------------------------------------------------
inline icarus::trigger::WindowTopologyManager::WindowTopologyManager(
  geo::GeometryCore const& geom,
  std::string const& logCategory /* = "WindowTopologyManager" */
  )
  : icarus::ns::util::mfLoggingClass(logCategory)
  , fVerify{ logCategory }
  , fGeom(&geom)
{}


//------------------------------------------------------------------------------
inline auto icarus::trigger::WindowTopologyManager::operator* () const
  -> WindowChannelMap const&
  { return fVerify.getTopology(); }


//------------------------------------------------------------------------------
inline auto icarus::trigger::WindowTopologyManager::operator-> () const
  -> WindowChannelMap const*
  { return &(fVerify.getTopology()); }


//------------------------------------------------------------------------------
//--- Template implementation
//------------------------------------------------------------------------------
//--- icarus::trigger::WindowTopologyVerification
// -----------------------------------------------------------------------------
template <typename Gates>
void icarus::trigger::WindowTopologyVerification::verifyOrThrow
  (Gates const& gates) const
{
  
  std::string const errorMsg = verify(gates);
  if (errorMsg.empty()) return;
  
  // put together the exception message and throw it.
  throw cet::exception("WindowTopologyVerification")
    << "Some channels from trigger gates do not match the previous window allocation:\n"
    << errorMsg
    << "\n"
    << "Window allocation: "
    << fWindowMap.value()
    ;
  
} // icarus::trigger::WindowTopologyVerification::operator()


//------------------------------------------------------------------------------
//--- icarus::trigger::WindowTopologyManager
//------------------------------------------------------------------------------
template <typename Gates>
std::string icarus::trigger::WindowTopologyManager::setOrVerify
  (Gates const& gates)
{
  if (fVerify.hasTopology()) {
    return fVerify.verify(gates);
  }
  else {
    extractTopology(gates);
    return {};
  }
} // icarus::trigger::WindowTopologyManager::setOrVerify()


//------------------------------------------------------------------------------
template <typename Gates>
bool icarus::trigger::WindowTopologyManager::operator() (Gates const& gates) {
  if (fVerify.hasTopology()) {
    fVerify(gates);
    return false;
  }
  else {
    extractTopology(gates);
    return true;
  }
} // icarus::trigger::WindowTopologyManager::operator()


//------------------------------------------------------------------------------
template <typename Gates>
void icarus::trigger::WindowTopologyManager::extractTopology(Gates const& gates)
{
  icarus::trigger::WindowTopologyAlg const topoMaker
    { *fGeom, logCategory() + "_Extractor" };
  fVerify.setTopology(topoMaker.createFromGates(gates));
} // icarus::trigger::WindowTopologyManager::extractTopology()


//------------------------------------------------------------------------------

#endif // ICARUSCODE_PMT_TRIGGER_ALGORITHMS_WINDOWTOPOLOGYALG_H

