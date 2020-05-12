/**
 * @file   icaruscode/Utilities/WeakCurrentType.h
 * @brief  A C++ type to describe the type of weak current.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   August 8, 2019
 * @see    `icaruscode/Utilities/WeakCurrentType.cxx`
 */

#ifndef ICARUSCODE_UTILITIES_WEAKCURRENTTYPE_H
#define ICARUSCODE_UTILITIES_WEAKCURRENTTYPE_H


// framework libraries
#include "cetlib_except/exception.h"

// C/C++ standard library
#include <string>


// --- BEGIN -- Weak current types ---------------------------------------------
/**
 * @defgroup ICARUS_WeakCurrentTypes Weak current types
 * 
 * Data structures and constants to represent weak current types.
 */
/// @{

//------------------------------------------------------------------------------
namespace icarus { class WeakCurrentType; }
/**
 * @brief Represents a type of weak current.
 * 
 * This type can be initialised with a string (see `parse()`) or a mnemonic
 * value (`CurrentType`).
 * 
 * Example of usage:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * icarus::WeakCurrentType const current { "neutral" };
 * 
 * std::cout << "Current type: " << std::string(current) << std::endl;
 * if (current == icarus::NeutralCurrentType) {
 *   std::cout << "Well, yeah." << std::endl;
 * }
 * 
 * // use WeakCurrentType::CurrentType in switch statements:
 * std::string boson;
 * switch (current) {
 *   case WeakCurrentType::CC:  boson = "W"; break;
 *   case WeakCurrentType::NC:  boson = "Z"; break;
 *   case WeakCurrentType::any: boson = "V"; break;
 * } // switch
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */
class icarus::WeakCurrentType {
  
    public:
  
  /// Type of weak current.
  enum CurrentType {
    CC, ///< Charged current.
    NC, ///< Neutral current.
    any ///< Any interaction.
  }; // enum CurrentType
  
  
  // --- BEGIN -- Constructors -------------------------------------------------
  
  /// Default constructor: matches any current type.
  constexpr WeakCurrentType() = default;
  
  /// Constructor: assigns the specified current type.
  constexpr WeakCurrentType(CurrentType type): fType(type) {}
  
  /// Constructor: assigns the value interpreting the specification `spec`.
  /// @see `WeakCurrentType::parse()`
  explicit WeakCurrentType(std::string const& spec)
    : WeakCurrentType(parse(spec)) {}
  
  // --- END -- Constructors ---------------------------------------------------
  
  
  
  // --- BEGIN -- String operations --------------------------------------------
  /// @name String operations
  /// @{
  
  /// Returns a string with the name of the current (`"charged"`, `"neutral"`,
  /// `"any"`).
  std::string name() const;
  
  /// Returns a string with the short name of the current
  /// (`"CC"`, `"NC"`, `"any"`).
  std::string shortName() const;
  
  /// Converts to the `name()` of the current.
  operator std::string() const { return name(); }
  
  /// @}
  // --- END -- String operations ----------------------------------------------
  
  
  // --- BEGIN -- Comparisons --------------------------------------------------
  /// @name Comparisons
  /// @{
  
  /// Returns whether this current type is the same as the `other`.
  constexpr bool operator== (WeakCurrentType const& other) const
    { return fType == other.fType; }
  
  /// Returns whether this current type is different than the `other`.
  constexpr bool operator!= (WeakCurrentType const& other) const
    { return fType != other.fType; }
  
  /// @}
  // --- END -- Comparisons ----------------------------------------------------
  
  
  // this allows using `switch` on `WeakCurrentType` objects
  operator int() const { return static_cast<int>(fType); }
  
  
  /**
   * @brief Converts a string into a `WeakCurrentType`.
   * @param spec the specification string to be converted
   * @throws cet::exception (category `"WeakCurrentType"`) if `spec` is not
   *           supported
   * 
   * The case of `spec` is irrelevant.
   * Accepted values include the shortened name (`"CC"`; see `shortName()`),
   * and the full name (`"charged"`; see `name()`).
   * Also the empty string is converted to a `CurrentType` of `any`.
   */
  static CurrentType parse(std::string const& spec);
  
  
    private:
  
  CurrentType fType = any; ///< Type of current stored.
  
  /// Returns a upper-case copy of `s`.
  static std::string to_upper(std::string const& s);
  
}; // class icarus::WeakCurrentType


//------------------------------------------------------------------------------
namespace icarus {
  
  /// Constant value for a weak neutral current type.
  constexpr WeakCurrentType NeutralCurrentType { WeakCurrentType::NC };

  /// Constant value for a weak charged current type.
  constexpr WeakCurrentType ChargedCurrentType { WeakCurrentType::CC };

  /// Constant value for any weak current type.
  constexpr WeakCurrentType AnyWeakCurrentType { WeakCurrentType::any };

} // namespace icarus


// --- END -- Weak current types -----------------------------------------------
/// @}

//------------------------------------------------------------------------------


#endif // ICARUSCODE_UTILITIES_WEAKCURRENTTYPE_H
