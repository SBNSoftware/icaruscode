/**
 * @file   icaruscode/PMT/OpReco/Algorithms/OpHitTimeSelector.h
 * @brief  Utility to extract the time of a reconstructed optical hit.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   April 11, 2023
 * 
 * This is a header-only library (so far).
 */

#ifndef ICARUSCODE_PMT_OPRECO_ALGORITHMS_OPHITTIMESELECTOR_H
#define ICARUSCODE_PMT_OPRECO_ALGORITHMS_OPHITTIMESELECTOR_H


// LArSoft libraries
#include "lardataalg/DetectorInfo/DetectorTimingTypes.h" // electronics_time
#include "lardataalg/Utilities/MultipleChoiceSelection.h"
#include "lardataalg/Utilities/quantities/spacetime.h" // microsecond
#include "lardataobj/RecoBase/OpHit.h"

// C/C++ standard libraries
#include "stdexcept" // std::logic_error
#include "string"
#include "utility" // std::exchange()


// -----------------------------------------------------------------------------
namespace recob { class OpHitTimeSelector; }
/**
 * @brief Extracts a configured time from reconstructed optical hits.
 * 
 * This class provides three functionalities:
 *  * `Type`: an enumerator listing all the optical hit times supported by the
 *    extractor;
 *  * `timeOf()` (also `operator()`): returns the type of time configured in the
 *    object;
 *  * `timeOf()` (`static` version): returns the type of time specified as
 *    argument.
 * 
 * Note that the exact definition of each of the times depends on the algorithm
 * used for the reconstruction of optical hits.
 * 
 * Example of usage:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * recob::OpHitTimeSelector const hitTime { recob::OpHitTimeSelector::Start };
 * 
 * for (recob::OpHit const& opHit: OpHits) {
 *   
 *   double const time = hitTime(opHit);
 * 
 *   // ... do something with it
 *   
 * }
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * 
 * The supported time types are defined in `Type` enumerator.
 * 
 */
class recob::OpHitTimeSelector {
  
    public:
  
  /// Type of optical hit reconstructed time.
  enum class Type {
      Peak   ///< `recob::OpHit::PeakTime()`: representing the maximum light.
    , Start  ///< `recob::OpHit::StartTime()`: representing the earliest light.
    , Rise   ///< `recob::OpHit::StartTime()` + `recob::OpHit::RiseTime()`: fine reconstruction of earliest light.
  }; // Type
  
  /// Type used for default construction.
  static constexpr Type DefaultType { Type::Peak };
  
  /// Helper to convert a type from and to a string.
  static util::MultipleChoiceSelection<Type> const TypeSelector;
  
  
  // --- BEGIN -- Constructors -------------------------------------------------
  
  /// Initializes the extractor with the specified `timeType`
  /// @param timeType (default: `DefaultType`) the type of time to extract
  OpHitTimeSelector(Type timeType = DefaultType);
  
  /// Initializes the extractor with the type specified by `timeTypeName`.
  /// @param timeType name of the type to extract (see `TypeSelector`)
  OpHitTimeSelector(std::string const& timeTypeName);
  
  // --- END ---- Constructors -------------------------------------------------
  
  
  // --- BEGIN --- Time type extraction management -----------------------------
  /// @name Time type extraction management
  /// @{
  
  /// Returns the type of time being extracted.
  Type getType() const { return fType; }
  
  /// Returns the type of time being extracted.
  std::string getTypeName() const;
  
  /// Changes the time type to be extracted to `newType`.
  /// @return the type that was being extracted before the change
  Type selectType(Type newType);
  
  /// Changes the time type to be extracted to the one with name `newTypeName`.
  /// @return the type that was being extracted before the change
  Type selectType(std::string const& newTypeName)
    { return selectType(parseTypeName(newTypeName)); }
  
  /// @}
  // --- END ----- Time type extraction management -----------------------------
  
  
  // --- BEGIN --- Extraction --------------------------------------------------
  /// @name Time extraction
  /// @{
  
  /// Returns the configured time of the hit, in the same scale as
  /// `recob::OpHit::PeakTime()`.
  /// @see `timeOf(recob::OpHit const&)`
  double operator() (recob::OpHit const& opHit) const { return timeOf(opHit); }
  
  /// Returns the configured time of the hit, in the same scale as
  /// `recob::OpHit::PeakTime()`.
  double timeOf(recob::OpHit const& opHit) const
    { return (*fTimeProc)(opHit); }
  
  /// Returns the configured time of the hit, in the same scale as
  /// `recob::OpHit::PeakTimeAbs()`.
  double absTimeOf(recob::OpHit const& opHit) const
    { return timeOf(opHit) + relToAbs(opHit); }
  
  /// Returns the configured time of the hit, on the
  /// @ref DetectorClocksElectronicsTime "electronics time scale".
  /// @note Hit absolute time is assumed to be electronics time scale already.
  detinfo::timescales::electronics_time elecTimeOf
    (recob::OpHit const& opHit) const
    { return { util::quantities::microsecond{ absTimeOf(opHit) } }; }
  

  /// Returns the `opHit` time of the specified `type`, in the same scale as
  /// `recob::OpHit::PeakTime()`.
  static double timeOf(recob::OpHit const& opHit, Type type)
    { return selectTimeProc(type)(opHit); }
  
  /// Returns the `opHit` time of the specified `type`, in the same scale as
  /// `recob::OpHit::PeakTimeAbs()`.
  static double absTimeOf(recob::OpHit const& opHit, Type type)
    { return timeOf(opHit, type) + relToAbs(opHit); }
  
  /// Returns the `opHit` time of the specified `type`, on the
  /// @ref DetectorClocksElectronicsTime "electronics time scale".
  /// @note Hit absolute time is assumed to be electronics time scale already.
  static detinfo::timescales::electronics_time elecTimeOf
    (recob::OpHit const& opHit, Type type)
    { return { util::quantities::microsecond{ absTimeOf(opHit, type) } }; }
  
  /// @}
  // --- END ----- Extraction --------------------------------------------------
  
    private:
  
  /// Type of method of this class used to extract the time.
  using TimeProc_t = double (*)(recob::OpHit const&);
  
  
  Type fType; ///< Type of time being extracted (for book-keeping).
  
  TimeProc_t fTimeProc; ///< Pointer to the method extracting the right time.
  
  
  /// Returns the offset to add to a relative time to make it absolute.
  /// @note Assumes both `recob::OpHit::PeakTime()` and
  ///       `recob::OpHit::PeakTimeAbs()` to be correctly set.
  static double relToAbs(recob::OpHit const& opHit);
  
  
  // --- BEGIN --- Time extraction methods -------------------------------------
  static double extractPeakTime(recob::OpHit const& opHit)
    { return opHit.PeakTime(); }
  
  static double extractStartTime(recob::OpHit const& opHit)
    { return opHit.StartTime(); }
  
  static double extractRiseTime(recob::OpHit const& opHit)
    { return opHit.StartTime() + opHit.RiseTime(); }
  
  // --- END ----- Time extraction methods -------------------------------------
  
  
  /// Returns the method to extract the specified time `type`.
  static TimeProc_t selectTimeProc(Type type);
  
  /// Returns the type associated to the specified `name` (or alias).
  static Type parseTypeName(std::string const& name)
    { return TypeSelector.parse(name).value(); }
  
}; // recob::OpHitTimeSelector


// -----------------------------------------------------------------------------
// ---  Inline implementation
// -----------------------------------------------------------------------------
inline util::MultipleChoiceSelection<recob::OpHitTimeSelector::Type> const
recob::OpHitTimeSelector::TypeSelector {
    { Type::Peak,  "Peak",  "PeakTime" }
  , { Type::Start, "Start", "StartTime" }
  , { Type::Rise,  "Rise",  "RiseTime" }
};


// -----------------------------------------------------------------------------
inline recob::OpHitTimeSelector::OpHitTimeSelector
  (Type timeType /* = DefaultType */)
  : fType{ timeType }, fTimeProc{ selectTimeProc(fType) }
{}


// -----------------------------------------------------------------------------
inline recob::OpHitTimeSelector::OpHitTimeSelector
  (std::string const& timeTypeName)
  : OpHitTimeSelector{ parseTypeName(timeTypeName) }
{}


// -----------------------------------------------------------------------------
inline std::string recob::OpHitTimeSelector::getTypeName() const {
  return TypeSelector.get(fType).name();
}


// -----------------------------------------------------------------------------
inline auto recob::OpHitTimeSelector::selectType(Type newType) -> Type {
  fTimeProc = selectTimeProc(newType);
  return std::exchange(fType, newType);
}


// -----------------------------------------------------------------------------
inline double recob::OpHitTimeSelector::relToAbs(recob::OpHit const& opHit) {
  return opHit.PeakTimeAbs() - opHit.PeakTime(); 
}


// -----------------------------------------------------------------------------
inline auto recob::OpHitTimeSelector::selectTimeProc(Type type) -> TimeProc_t {
  switch (type) {
    case Type::Peak:  return &recob::OpHitTimeSelector::extractPeakTime;
    case Type::Start: return &recob::OpHitTimeSelector::extractStartTime;
    case Type::Rise:  return &recob::OpHitTimeSelector::extractRiseTime;
    default: throw std::logic_error{ "Unsupported time type" };
  } // switch
}


// -----------------------------------------------------------------------------


#endif // ICARUSCODE_PMT_OPRECO_ALGORITHMS_OPHITTIMESELECTOR_H


// this section was envisioned to be set into a different file;
// I am leaving it here until I am forced to depend on FHiCL
// (right now I am avoiding that by using templates).

#ifndef ICARUSCODE_PMT_OPRECO_ALGORITHMS_OPHITTIMESELECTOR_FHICL_H
#define ICARUSCODE_PMT_OPRECO_ALGORITHMS_OPHITTIMESELECTOR_FHICL_H


namespace recob {
  
  /**
   * @brief Configuration helper for `recob::OpHitTimeSelector::Type` parameter.
   * @param typeName type name
   * @param paramName (default: empty) name used for error message
   * @return the `Type` specified by `configAtom`
   * @throw std::runtime_error if the string is not valid
   * 
   * In case of invalid `name`, an exception is thrown mentioning the
   * `paramName` (if non-empty) and the list of allowable values.
   */
  recob::OpHitTimeSelector::Type opHitTimeType
    (std::string const& typeName, std::string const& paramName = "") {
    try {
      return recob::OpHitTimeSelector::TypeSelector.parse(typeName).value();
    }
    catch (util::MultipleChoiceSelectionBase::UnknownOptionError const& e)
    {
      throw std::runtime_error{
        "Invalid value for "
        + (paramName.empty()? "optical hit time": ("'" + paramName + "'"))
        + " parameter: '" + e.label() + "'; valid options: "
        + recob::OpHitTimeSelector::TypeSelector.optionListString() + ".\n"
        };
    }
  } // opHitTimeType(std::string)
  
  /**
   * @brief Configuration helper for `recob::OpHitTimeSelector::Type` parameter.
   * @tparam ConfigAtom configuration class type holding the time type name
   * @param configAtom the configuration object holding the time type name
   * @return the `Type` specified by `configAtom`
   * @throw std::runtime_error if the string is not valid
   * @see `opHitTimeType(std::string const&, std::string const&)`
   * 
   * A version of `opHitTimeType()` serving FHiCL-like parameters.
   * 
   * 
   * Requirements
   * -------------
   * 
   * The type `ConfigAtom` is expected to have an interface like
   * `fhicl::Atom<std::string>`, and in particular:
   * 
   * * `operator() const` returning the value of the parameter in an object
   *   convertible to string (note: `fhicl::OptionalAtom` won't do).
   * * `name() const` returning the name of the configuration parameter.
   */
  template <typename ConfigAtom>
  recob::OpHitTimeSelector::Type opHitTimeType(ConfigAtom const& configAtom)
    { return opHitTimeType(configAtom(), configAtom.name()); }
  
} // namespace recob


#endif // ICARUSCODE_PMT_OPRECO_ALGORITHMS_OPHITTIMESELECTOR_FHICL_H

