/**
 * @file   icaruscode/Decode/DecoderTools/details/KeyValuesData.h
 * @brief  Simple parsed data format.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   May 9, 2021
 * 
 * This library is header only.
 */

#ifndef ICARUSCODE_DECODE_DECODERTOOLS_DETAILS_KEYVALUESDATA_H
#define ICARUSCODE_DECODE_DECODERTOOLS_DETAILS_KEYVALUESDATA_H


// C++ standard libraries
#include <iosfwd> // std::ostream
#include <vector>
#include <string_view>
#include <string>
#include <optional>
#include <stdexcept> // std::runtime_error
#include <utility> // std::move()
#include <charconv> // std::from_chars()
#include <type_traits> // std::is_floating_point_v, std::enable_if_t
#include <cstddef> // std::size_t


// -----------------------------------------------------------------------------
namespace icarus {
  
  namespace details {
    template <typename T, typename Enable = void>
    struct KeyValuesConverter;
  }
  
  struct KeyValuesData;
  
  std::ostream& operator<< (std::ostream& out, KeyValuesData const& data);
  
} // namespace icarus


// -----------------------------------------------------------------------------
/**
 * @class icarus::KeyValuesData
 * @brief Collection of items with key/values structure.
 * 
 * This data type is a collection of `Item` objects that contain unparsed
 * strings. Each `Item` has a key and a sequence of values. The values can be
 * queried as strings or as other data types. Conversions are performed by
 * `icarus::details::KeyValuesConverter`, which can be specialized with
 * the needed conversion logic. Only one type of conversion to any given type
 * is supported.
 * 
 * Each converter object specialization for a type `T` should support a call
 * with argument `std::string` returning a `std::optional<T>`.
 * 
 */
struct icarus::KeyValuesData {
  
  // --- BEGIN --- Exception definition ----------------------------------------
  /// @name Exceptions
  /// @{
  
  struct Error;
  struct DuplicateKey;
  struct ConversionFailed;
  struct ItemNotFound;
  struct ValueNotAvailable;
  
  /// @}
  // --- END ----- Exception definition ----------------------------------------
  
  /// Representation of a single item of data: a key and several values.
  struct Item {
    std::string key;
    std::vector<std::string> values;
    
    
    /// Constructs a new item assigning it a key (which should not be changed).
    Item(std::string key): key(std::move(key)) {}
    
    
    // --- BEGIN -- Setting ----------------------------------------------------
    /// @name Setting interface
    /// @{
    
    //@{
    /// Appends a value to this key.
    Item& addValue(std::string value)
      { values.push_back(std::move(value)); return *this; }
    Item& addValue(std::string_view value)
      { return addValue(std::string{ value }); }
    //@}
    
    /// Appends a sequence of values to this key.
    template <typename BIter, typename EIter>
    Item& addValues(BIter begin, EIter end)
      { while (begin != end) addValue(*begin++); return *this; }
    
    /// @}
    // --- END ---- Setting ----------------------------------------------------
    
    
    // --- BEGIN -- Query ------------------------------------------------------
    /// @name Query interface
    /// @{
    
    /// Returns the number of values currently present.
    std::size_t nValues() const noexcept { return values.size(); }
    
    /// Converts the value to `T`, throws on failure.
    template <typename T>
    T getNumber(std::size_t index) const;
    
    /// Converts the value to `T`, throws on conversion failures
    /// unless `ignoreFormatErrors` is `true`.
    template <typename T>
    std::optional<T> getOptionalNumber
      (std::size_t index, bool ignoreFormatErrors = false) const;
    
    /// @}
    // --- END ---- Query ------------------------------------------------------
    
    
    /// Lexicographic order by key (case-sensitive).
    bool operator< (Item const& other) const noexcept
      { return key < other.key; }
    
    
      private:
    
    /// Conversion functions.
    template <typename T>
    std::optional<T> convertStringInto(std::string const& valueStr) const
      { return icarus::details::KeyValuesConverter<T>{}(valueStr); }
    
  }; // struct Item
  
  
  // --- BEGIN -- Setter interface ---------------------------------------------
  /// @name Setter interface
  /// @{
  
  /// @brief Creates and registers a new item with the specified `key`.
  /// @return the newly created item for modifications
  Item& makeItem(std::string key);
  
  /// @brief Creates or retrieves an item with the specified `key`.
  /// @return the newly created or existing item for modifications
  Item& makeOrFetchItem(std::string const& key);
  
  /// Returns the item with specified `key`, `nullptr` if none.
  Item* findItem(std::string const& key) noexcept;
  
  /// @}
  // --- END ---- Setter interface ---------------------------------------------
  
  
  // --- BEGIN -- Query interface ----------------------------------------------
  /// @name Query interface
  /// @{
  
  /// Returns the item with specified `key`, `nullptr` if none.
  Item const* findItem(std::string const& key) const noexcept;
  
  /// Returns the item with specified `key`, throws `std::out_of_range` if none.
  Item const& getItem(std::string const& key) const;
  
  /// Returns whether an item with the specified key is present.
  bool hasItem(std::string const& key) const noexcept;
  
  /// Returns whether there is no item in data.
  bool empty() const noexcept;
  
  /// Returns the number of items in the data.
  std::size_t size() const noexcept;
  
  /// Returns a forward-iterable list of references to items.
  decltype(auto) items() const noexcept;
  
  /// @}
  // --- END ---- Query interface ----------------------------------------------
  
    private:
  
  std::vector<Item> fItems; ///< Collection of data items.
  
  
  /// Creates, registers and return a new item (assumed not to exist yet).
  Item& makeItemImpl(std::string key);
  
}; // icarus::KeyValuesData


// -----------------------------------------------------------------------------
// ---  Exception class definitions
// -----------------------------------------------------------------------------
struct icarus::KeyValuesData::Error: public std::runtime_error {
  
  Error(std::string msg): std::runtime_error(std::move(msg)) {}
  
}; // icarus::KeyValuesData::Error()


// -----------------------------------------------------------------------------
struct icarus::KeyValuesData::DuplicateKey: public Error {
  
  DuplicateKey(std::string const& msg)
    : Error("KeyValuesData::DuplicateKey: '" + msg + '\'')
    {}
  
}; // icarus::KeyValuesData::DuplicateKey()


// -----------------------------------------------------------------------------
struct icarus::KeyValuesData::ConversionFailed: public Error {
  
  ConversionFailed(std::string const& s, std::string const& tname = "")
    : Error{
      "conversion of '" + s + "'"
       + (tname.empty()? "": (" to type '" + tname + "'")) + " failed"
      }
    {}
  
  template <typename T>
  static ConversionFailed makeFor(std::string const& s)
    { return { s, typeid(T).name() }; }
  
}; // icarus::KeyValuesData::ConversionFailed()


// -----------------------------------------------------------------------------
struct icarus::KeyValuesData::ItemNotFound: public Error {
  
  ItemNotFound(std::string const& key): Error("item not found: '" + key + '\'')
    {}
  
}; // icarus::KeyValuesData::ItemNotFound()


// -----------------------------------------------------------------------------
struct icarus::KeyValuesData::ValueNotAvailable: public Error {
  
  ValueNotAvailable(std::size_t index)
    : Error("item value #" + std::to_string(index) + " not available")
    {}
  
}; // icarus::KeyValuesData::ValueNotAvailable()


// -----------------------------------------------------------------------------
// --- Template implementation
// -----------------------------------------------------------------------------
// ---  icarus::KeyValuesData::Item
// -----------------------------------------------------------------------------
template <typename T>
T icarus::KeyValuesData::Item::getNumber(std::size_t index) const {
  
  if (index >= values.size()) throw ValueNotAvailable(index);
  
  auto const& valueStr = values[index];
  auto const number = convertStringInto<T>(valueStr);
  return number? number.value(): throw ConversionFailed::makeFor<T>(valueStr);
  
} // icarus::KeyValuesData<>::Item::getNumber<>()


// -----------------------------------------------------------------------------
template <typename T>
std::optional<T> icarus::KeyValuesData::Item::getOptionalNumber
  (std::size_t index, bool ignoreFormatErrors /* = false */) const
{
  if (index < values.size()) return std::nullopt;

  auto const& valueStr = values[index];
  auto const number = convertStringInto<T>(valueStr);
  return (number || ignoreFormatErrors)
    ? number: throw ConversionFailed::makeFor<T>(valueStr);

} // icarus::KeyValuesData<>::Item::getOptionalNumber()


// -----------------------------------------------------------------------------
// ---  icarus::details::KeyValuesConverter
// -----------------------------------------------------------------------------
template <typename T, typename /* = void */>
struct icarus::details::KeyValuesConverter {
  
  //@{
  /// Convert a string `s` into a type `T`;
  /// may return `std::nullopt` on "non-fatal" failure.
  std::optional<T> operator() (std::string const& s) const
    { return convert(s); }
  
  static std::optional<T> convert(std::string const& s)
    {
      T number {}; // useless initialization to avoid GCC complains
      char const *b = s.data(), *e = b + s.length();
      return (std::from_chars(b, e, number).ptr == e)
        ? std::make_optional(number): std::nullopt;
    } // convert()
  
  //@}
}; // icarus::details::KeyValuesConverter<>


// ---  BEGIN  --- WORKAROUND --------------------------------------------------
/*
 * Compilers like GCC 9.3.0 and Clang 7.0 are not fully C++17 compliant yet.
 * `std::from_chars()` is not provided for floating points types.
 * The first compiler versions supporting them are GCC 11.1 and Clang 12.0.1.
 * I am providing a specialization to cover from that, until the compilers
 * are updated.
 */
#ifdef __GNUC__
#  if (__GNUC__ >= 11)
     // GCC 11 should support std::from_chars() for floating point types;
     // remove this #if branch, and if Clang's is also already removed,
     // remove the workaround too
#    error "Redundant workaround on std::from_chars() for GCC"
#  else
#    define ICARUSCODE_DECODE_DECODERTOOLS_DETAILS_KEYVALUESDATA_NEEDS_FROMCHARS_FLOAT
#  endif // __GNUC__
#endif // __GNUC__

#ifdef __clang_major__
#  if (__clang_major__ >= 12) || ((__clang_major__ == 12) && ((__clang_minor__ >= 1) || (__clang_patchlevel__ >= 1)))
     // Clang 12.0.1 should support std::from_chars() for floating point types;
     // remove this #if branch, and if GCC's is also already removed,
     // remove the workaround too
#    error "Redundant workaround on std::from_chars() for GCC"
#  else
#    define ICARUSCODE_DECODE_DECODERTOOLS_DETAILS_KEYVALUESDATA_NEEDS_FROMCHARS_FLOAT
#  endif // __clang_major__
#endif // __clang_major__

#ifdef ICARUSCODE_DECODE_DECODERTOOLS_DETAILS_KEYVALUESDATA_NEEDS_FROMCHARS_FLOAT

#include <sstream>

template <typename T>
struct icarus::details::KeyValuesConverter
  <T, std::enable_if_t<std::is_floating_point_v<T>>>
{
  
  std::optional<T> operator() (std::string const& s) const
    { return convert(s); }
  
  static std::optional<T> convert(std::string const& s)
    {
      T number {}; // useless initialization to avoid GCC complains
      std::istringstream sstr{ s };
      sstr >> number;
      // check that no non-space character is left in the stream
      return (sstr && (sstr >> std::ws).eof())
        ? std::make_optional(number): std::nullopt;
    } // convert()
  
}; // icarus::details::KeyValuesConverter<floating point>

#endif
// ---  END ------ WORKAROUND --------------------------------------------------

// -----------------------------------------------------------------------------

#endif // ICARUSCODE_DECODE_DECODERTOOLS_DETAILS_KEYVALUESDATA_H
