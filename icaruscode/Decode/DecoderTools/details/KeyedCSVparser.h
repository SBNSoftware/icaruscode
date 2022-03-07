/**
 * @file   icaruscode/Decode/DecoderTools/details/KeyedCSVparser.h
 * @brief  Simple parser for comma-separated text.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   May 9, 2021
 * @see    icaruscode/Decode/DecoderTools/details/KeyedCSVparser.cxx
 */

#ifndef ICARUSCODE_DECODE_DECODERTOOLS_DETAILS_KEYEDCVSPARSER_H
#define ICARUSCODE_DECODE_DECODERTOOLS_DETAILS_KEYEDCVSPARSER_H


// C++ standard libraries
#include <iosfwd> // std::ostream
#include <string_view>
#include <vector>
#include <string>
#include <optional>
#include <regex>
#include <initializer_list>
#include <stdexcept> // std::runtime_error
#include <utility> // std::move(), std::pair
#include <limits>
#include <type_traits> // std::is_constructible_v, std::is_arithmetic_v, ...
#include <charconv> // std::from_chars()
#include <cstddef> // std::size_t


// -----------------------------------------------------------------------------
namespace icarus::details {
  class KeyedCSVparser;
  struct KeyValuesData;
  
  std::ostream& operator<< (std::ostream& out, KeyValuesData const& data);
} // namespace icarus::details


// -----------------------------------------------------------------------------
/**
 * @brief A set of key/values items.
 * 
 * This class collects `Item` objects, each being a string key and a sequence
 * of zero or more values. An specific item can be accessed by its key
 * (`findItem()`, `getItem()`) or all items may be iterated through (`items()`).
 * 
 * 
 * Initialization and updates
 * ---------------------------
 * 
 * A `KeyValuesData` object always starts empty (default constructor).
 * A new item is also always created empty (`makeItem()`) and with a set key.
 * 
 * After an empty item (i.e. an item with a key but no values) is created,
 * values can be added using the `Item` subclass interface (`addValue()`).
 * Item values and keys can be modified by changing `Item` data members
 * directly; note however that the mean to do that is the reference that
 * the `makeItem()` call returned, as there is no other way to access an item
 * in a writeable way.
 * 
 * This object will refuse to create a new item with the same key as an existing
 * one. However, the key of the item may be changed to any value after
 * `makeItem()` is called, although the interface does not encourage that.
 * If this introduces a duplicate key, the query functions will systematically
 * retrieve only one of the items with the repeated key (which one and whether
 * always the same one are undefined).
 * 
 * Example:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * icarus::details::KeyValuesData data;
 * data.makeItem("TriggerType").addValue("S5");
 * data.makeItem("Triggers");
 * data.makeItem("TriggerWindows").addValue("0C0B");
 * data.makeItem("TPChits")
 *   .addValue("12").addValue("130").addValue("0").addValue("0");
 * data.makeItem("TPChitTimes")
 *   .addValue("3").addValue("-1.1").addValue("-0.3").addValue("0.1");
 * data.makeItem("PMThits").addValue("8");
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * creates a `data` object with four items, a `TriggerType` one with one value,
 * a `Triggers` one with no values, a `TriggerWindows` with a single value
 * (expressed as a hexadecimal number), a `TPChits` one with four values,
 * a `TPChitTimes` with four values (the first meant to be the number of
 * remaining ones) and a `PMThits` with one.
 * 
 * 
 * Query
 * ------
 * 
 * The interface is quite terse.
 * 
 * General query methods reports whether there is no item in the object
 * (`empty()`) and how many items there are (`size()`).
 * 
 * A item with a known key can be retrieved (`getItem()`, `findItem()`), or
 * its existence may be tested (`hasItem()`).
 * 
 * Finally, all items can be iterated (`items()`). In this case, the items
 * are presented in the creation order.
 * 
 * The `Item` interface is documented in that subclass.
 * 
 * Example using the `data` object from the previous example:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * std::string triggerType = data.getItem("TriggerType").values()[0];
 * std::vector<int> triggers = data.getItem("Triggers").getVector<int>();
 * std::uint32_t triggerWindowBits
 *  = data.getItem("TriggerWindows").getNumber<std::uint32_t>(0, 16); // base 16
 * std::vector<int> TPChits = data.getItem("TPChits").getVector<int>();
 * std::vector<float> TPCtimes
 *  = data.getItem("TPChitTimes").getSizedVector<float>();
 * std::vector<int> CRThits;
 * if (auto const* item = data.findItem("CRThits"))
 *   CRThits = item->getVector<int>();
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * 
 */
struct icarus::details::KeyValuesData {
  
  struct Error;
  struct ErrorOnKey;
  struct DuplicateKey;
  struct ConversionFailed;
  struct ItemNotFound;
  struct ValueNotAvailable;
  struct MissingSize;
  struct WrongSize;
  
  /**
   * @brief Representation of a single item of data: a key and several values.
   * 
   * Values can be added directly accessing the `values` data member, or
   * with `addValue()` method call.
   * 
   * Access to the values happens by index, with `nValues()` indices starting
   * from `0` available. Direct access to `values()` is expected, and additional
   * shortcuts are available:
   * * `getAs()`, `getOptionalAs()` to convert to an arbitrary type; the
   *     standard conversion uses `from_chars()`, and specialization is
   *     possible
   * * `getNumber()`, `getOptionalNumber()` to convert to a number
   * * `getVector()` to convert to a vector of values (each like in `getAs()`)
   * * `getSizedVector()` to convert to a vector of values; the first value is
   *   the size of the actual values.
   */
  struct Item: public std::pair<std::string, std::vector<std::string>> {
    
    using pair_t = std::pair<std::string, std::vector<std::string>>;
    
    /**
     * @brief Converter class for interpreting strings as type `T`.
     * @tparam T type to convert strings into
     * 
     * This class is designed to allow handling special types via
     * specialization. The standard implementation uses `from_chars()`
     * (including the one from `std`, but also any pulled in by C++ name lookup)
     * to perform the conversion, and returns a empty optional on failure.
     * Specializations may decide to throw an exception instead
     * (type `ConversionFailed` or derivatives are suggested).
     */
    template <typename T, typename = void>
    struct StringConverter {
      unsigned int base;
      StringConverter(unsigned int base = 10): base(base) {}
      std::optional<T> operator() (std::string const& valueStr) const;
    }; // StringConverter
    
    std::string const& key() const { return first; } ///< Key of the item.
    /// All item values, as strings.
    std::vector<std::string>& values() { return second; }
    std::vector<std::string> const& values() const { return second; }
    
    /// Constructor: sets the key (mandatory).
    Item(std::string key): pair_t{ std::move(key), {} } {}
    
    /// Adds a string value to the list of values.
    /// @return this same object (allows queueing calls in the same statement)
    Item& addValue(std::string value)
      { values().push_back(std::move(value)); return *this; }
    
    /// The current number of values in this item.
    std::size_t nValues() const noexcept { return values().size(); }
    
    /**
     * @brief Returns the requested value, converted into type `T`
     * @tparam T type to convert the value into
     * @tparam Conv type of a functor for conversion of the value into type `T`
     * @param index the index of the requested value
     * @param converter a functor for conversion of the value into type `T`
     * @return the requested value as an object of type `T`
     * @throw ValueNotAvailable if no value is available with that index
     * @throw ConversionFailed if the value could not be converted to type `T`
     * 
     * Conversion is performed via `converter` object, functor taking a string
     * and returning an object of type `std::optional<T>`. The functor can
     * decline the conversion by returning an empty `std::optional`, or directly
     * throw an exception on error.
     */
    template <typename T, typename Conv>
    T getAs(std::size_t index, Conv converter) const;
    
    /**
     * @brief Returns the requested value, converted into type `T`
     * @tparam T type to convert the value into
     * @param index the index of the requested value
     * @return the requested value as an object of type `T`
     * @throw ValueNotAvailable if no value is available with that index
     * @throw ConversionFailed if the value could not be converted to type `T`
     * 
     * Conversion is performed via an helper class `StringConverter` which can
     * be specialized if needed, and that uses `from_chars()` for conversion.
     */
    template <typename T>
    T getAs(std::size_t index) const;
    
    
    /**
     * @brief Returns the requested value, converted into type `T`
     * @tparam T type to convert the value into
     * @tparam IgnoreFormatErrors (default: `true`) how to treat conversion
     *                            errors
     * @tparam Conv type of a functor for conversion of the value into type `T`
     * @param index the index of the requested value
     * @param converter a functor for conversion of the value into type `T`
     * @return the requested value, or an empty optional on failure
     * @throw ConversionFailed if the value could not be converted to type `T`
     * 
     * Conversion is performed via `converter` object, functor taking a string
     * and returning an object of type `std::optional<T>`. The functor can
     * decline the conversion by returning an empty `std::optional`, or directly
     * throw an exception on error.
     * 
     * If no value is available for the specified `index`, an empty optional
     * is returned.
     * 
     * An exception is thrown on conversion failures unless `IgnoreFormatErrors`
     * is `true`, in which case an empty optional is also returned.
     */
    template <typename T, bool IgnoreFormatErrors = false, typename Conv>
    std::optional<T> getOptionalAs(std::size_t index, Conv converter) const;
    
    /**
     * @brief Returns the requested value, converted into type `T`
     * @tparam T type to convert the value into
     * @tparam IgnoreFormatErrors (default: `true`) how to treat conversion
     *                            errors
     * @param index the index of the requested value
     * @param ignoreFormatErrors (default: `false`) ignore conversion errors
     * @return the requested value, or an empty optional on failure
     * @throw ConversionFailed if the value could not be converted to type `T`
     * 
     * Conversion is performed via `converter` object, functor taking a string
     * and returning an object of type `std::optional<T>`. The functor can
     * decline the conversion by returning an empty `std::optional`, or directly
     * throw an exception on error.
     * 
     * If no value is available for the specified `index`, an empty optional
     * is returned.
     * 
     * An exception is thrown on conversion failures unless `IgnoreFormatErrors`
     * is `true`, in which case an empty optional is also returned.
     */
    template <typename T, bool IgnoreFormatErrors = false>
    std::optional<T> getOptionalAs(std::size_t index) const;
    
    
    /**
     * @brief Returns the requested value, converted into a number of type `T`
     * @tparam T type of number to convert the value into
     * @param index the index of the requested value
     * @param base (default: `10`) numerical base of the input number
     * @return the requested value as a number of type `T`
     * @throw ValueNotAvailable if no value is available with that index
     * @throw ConversionFailed if the value could not be converted to type `T`
     * 
     * See `getAs()` for details.
     * 
     * Note that the number must have no base prefix (e.g. `"F5"` for
     * hexadecimal rather than `"0xF5"`).
     */
    template <typename T>
    T getNumber(std::size_t index, unsigned int base = 10) const;
    
    /**
     * @brief Returns the requested value, converted into a number of type `T`
     * @tparam T type of number to convert the value into
     * @param index the index of the requested value
     * @param base (default: `10`) numerical base of the input number
     * @return the requested value, or an empty optional on failure
     * @throw ConversionFailed if the value could not be converted to type `T`
     * 
     * See `getOptionalAs()` for details.
     * 
     * Note that the number must have no base prefix (e.g. `"F5"` for
     * hexadecimal rather than `"0xF5"`).
     */
    template <typename T>
    std::optional<T> getOptionalNumber
      (std::size_t index, unsigned int base = 10) const;
    
    
    /**
     * @brief Returns all the values, each converted into type `T`
     * @tparam T type to convert the value into
     * @tparam Conv type of a functor for conversion of the value into type `T`
     * @param converter a functor for conversion of the value into type `T`
     * @return a vector with all the converted values
     * @throw ConversionFailed if any value could not be converted to type `T`
     * 
     * Conversion of each element is performed by `getAs<T, Conv>()`.
     * 
     * An exception is thrown on any conversion failure.
     */
    template <typename T, typename Conv = StringConverter<T>>
    std::vector<T> getVector(Conv converter = {}) const;
    
    
    /**
     * @brief Returns all the values, each converted into type `T`
     * @tparam T type to convert the value into
     * @tparam Conv type of a functor for conversion of the value into type `T`
     * @param converter a functor for conversion of the value into type `T`
     * @return a vector with all the converted values
     * @throw MissingSize on any error converting the first value to a size
     * @throw WrongSize if the actual number of values does not match the size
     * @throw ConversionFailed if any value could not be converted to type `T`
     * 
     * The first value (mandatory) is converted to represent the size of the
     * vector. That is used as verification when converting all the other
     * elements: if there is the wrong number of elements, an exception is
     * thrown.
     * 
     * Conversion of each element is performed by `getAs<T, Conv>()`.
     * 
     * An exception is also thrown on conversion failure of any of the values.
     */
    template <typename T, typename Conv = StringConverter<T>>
    std::vector<T> getSizedVector(Conv converter = Conv{}) const;
    
    
    /// Sorting operator (by key, lexicographic).
    bool operator< (Item const& other) const noexcept
      { return key() < other.key(); }
    
      private:

    template <typename T, typename Iter, typename Conv>
    std::vector<T> convertVector(Iter begin, Iter end, Conv converter) const;

    template <typename T>
    static std::optional<T> convertStringInto(std::string const& valueStr);
    
  }; // struct Item
  
  Item& makeItem(std::string key);
  
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
  
    private:
  
  std::vector<Item> fItems; ///< Collection of data items.
  
}; // struct icarus::details::KeyValuesData


// -----------------------------------------------------------------------------
/**
 * @brief Parser to fill a `KeyValuesData` structure out of a character buffer.
 * 
 * It currently supports only single-line buffer.
 * 
 * The parser operates one "line" at a time, returning a `KeyValuesData` with
 * the values assigned to each detected key. No data type is implied: all
 * elements are treated as strings, either a key or a value.
 * The parser separates the elements according to a separator, strips them of
 * trailing and heading spaces, then it decides whether each element is a value
 * to be assigned to the last key found, or a new key.
 * Keys are elements that have letters in them, values are anything else.
 * This simple (and arguable) criterion can be broken with specific parser
 * configuration: a pattern can be specified that when matched to an element
 * will make it a key; the pattern can also set the number of values that key
 * will require.
 * 
 * For example:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * icarus::details::KeyedCSVparser parser;
 * parser.addPatterns({
 *       { "TriggerType", 1U } // expect one value (even if contains letters)
 *     , { "TriggerWindows", 1U } // expect one value (even if contains letters)
 *     , { "TPChitTimes", icarus::details::KeyedCSVparser::FixedSize }
 *          // the first value is an integer, count of how many other values
 *   });
 * 
 * icarus::details::KeyValuesData data = parser(
 *   "TriggerType, S5, Triggers, TriggerWindows, 0C0B,"
 *   " TPChits, 12, 130, 0, 0, TPChitTimes, 3, -1.1, -0.3, 0.1, PMThits, 8"
 *   );
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * will return `data` with 6 items.
 */
class icarus::details::KeyedCSVparser {
  
    public:
  
  using ParsedData_t = icarus::details::KeyValuesData;
  
  using Error = KeyValuesData::Error; ///< Base of all errors by KeyedCSVparser.
  using ErrorOnKey = KeyValuesData::ErrorOnKey; ///< Base of some errors.
  struct ParserError; ///< Generic error: base of all errors by KeyedCSVparser.
  struct InvalidFormat; ///< Parsing format is not understood.
  /// Expected number of values is missing.
  using MissingSize = KeyValuesData::MissingSize;
  struct MissingValues; ///< Expected values are missing.
  
  /// Mnemonic size value used in `addPattern()` calls.
  static constexpr unsigned int FixedSize
    = std::numeric_limits<unsigned int>::max();
  /// Mnemonic size value used in `addPattern()` calls.
  static constexpr unsigned int DynamicSize = FixedSize - 1U;
  
  
  /// Constructor: specifies the separator character.
  KeyedCSVparser(char sep = ','): fSep(sep) {}
  
  //@{
  /// Parses the buffer `s` and returns a data structure with the content.
  ParsedData_t parse(std::string_view const& s) const;
  ParsedData_t parse(std::string const& s) const;
  template <typename BIter, typename EIter>
  ParsedData_t parse(BIter b, EIter e) const;
  
  ParsedData_t operator() (std::string_view const& s) const { return parse(s); }
  ParsedData_t operator() (std::string const& s) const { return parse(s); }
  template <typename BIter, typename EIter>
  ParsedData_t operator() (BIter b, EIter e) const { return parse(b, e); }
  //@}
  
  //@{
  /// Parses the buffer `s` and fills `data` with it.
  void parse(std::string_view const& s, ParsedData_t& data) const;
  //@}
  
  /**
   * @name Know patterns
   * 
   * The parser normally treats as a value everything that does not start with a
   * letter.
   * Known patterns may override this behaviour: if a token matches a known
   * pattern, it is considered a key and it is possible to specify the expected
   * number of values.
   * 
   * The number of values can be:
   * * a number: exactly that number of values are required; an exception will
   *   be thrown if not enough tokens are available;
   * * `FixedSize`: the next token must be a non-negative integer specifying how
   *   many other values to add (read this with `Item::getSizedVector()`);
   *   an exception will be thrown if not enough tokens are available;
   * * `DynamicSize`: the standard algorithm is used and values are added as
   *   long as they don't look like keys; the token matching the pattern is
   *   interpreted as a key though.
   * 
   * Patterns are considered in the order they were added.
   */
  /// @{
  
  //@{
  /**
   * @brief Adds a single known pattern.
   * @param pattern the regular expression matching the key for this pattern
   * @param values the number of values for this pattern
   * @return this parser (`addPattern()` calls may be chained)
   */
  KeyedCSVparser& addPattern(std::regex pattern, unsigned int values)
    { fPatterns.emplace_back(std::move(pattern), values); return *this; }
  KeyedCSVparser& addPattern(std::string const& pattern, unsigned int values)
    { return addPattern(std::regex{ pattern }, values); }
  //@}
  
  //@{
  /**
   * @brief Adds known patterns.
   * @param patterns sequence of patterns to be added
   * @return this parser (`addPatterns()` calls may be chained)
   * 
   * Each pattern is a pair key regex/number of values, like in `addPattern()`.
   */
  KeyedCSVparser& addPatterns
    (std::initializer_list<std::pair<std::regex, unsigned int>> patterns);
  KeyedCSVparser& addPatterns
    (std::initializer_list<std::pair<std::string, unsigned int>> patterns);
  //@}
  
  /// @}
  
    private:
  using Buffer_t = std::string_view;
  using SubBuffer_t = std::string_view;
  
  char const fSep = ','; ///< Character used as token separator.
  
  /// List of known patterns for matching keys, and how many values they hold.
  std::vector<std::pair<std::regex, unsigned int>> fPatterns;
  
  /// Returns the length of the next toke, up to the next separator (excluded).
  std::size_t findTokenLength(Buffer_t const& buffer) const noexcept;

  /// Returns the value of the next token, stripped.
  SubBuffer_t peekToken(Buffer_t const& buffer) const noexcept;
  
  /// Extracts the next token from the `buffer` and returns its value, stripped.
  SubBuffer_t extractToken(Buffer_t& buffer) const noexcept;
  
  /// Is content of `buffer` a key (as opposed to a value)?
  bool isKey(SubBuffer_t const& buffer) const noexcept;
  
  
  template <typename String>
  static Buffer_t makeBuffer(String const& s) noexcept;

  static Buffer_t& moveBufferHead(Buffer_t& buffer, std::size_t size) noexcept;
  
  static SubBuffer_t strip(SubBuffer_t s) noexcept;
  static SubBuffer_t stripLeft(SubBuffer_t s) noexcept;
  static SubBuffer_t stripRight(SubBuffer_t s) noexcept;
  static SubBuffer_t stripRightChar(SubBuffer_t s, char c) noexcept;
  
  template <char... Chars>
  static SubBuffer_t stripRightChars(SubBuffer_t s) noexcept;
  
}; // icarus::details::KeyedCSVparser



// -----------------------------------------------------------------------------
// ---  Exception class definitions
// -----------------------------------------------------------------------------
struct icarus::details::KeyValuesData::Error: public std::runtime_error {
  
  Error(std::string msg): std::runtime_error(std::move(msg)) {}
  
}; // icarus::details::KeyValuesData::Error()


// -----------------------------------------------------------------------------
struct icarus::details::KeyValuesData::ErrorOnKey: public Error {
  
  ErrorOnKey(std::string const& key, std::string const& msg)
    : Error("Key '" + key + "': " + msg) {}
  
}; // icarus::details::KeyValuesData::Error()


// -----------------------------------------------------------------------------
struct icarus::details::KeyValuesData::DuplicateKey: public ErrorOnKey {
  
  DuplicateKey(std::string const& key)
    : ErrorOnKey(key, "duplicate key")
    {}
  
}; // icarus::details::KeyValuesData::DuplicateKey()


// -----------------------------------------------------------------------------
struct icarus::details::KeyValuesData::ConversionFailed: public ErrorOnKey {
  
  ConversionFailed(
    std::string const& key, std::string const& s, std::string const& tname = ""
    )
    : ErrorOnKey{
      key,
      "conversion of '" + s + "'"
       + (tname.empty()? "": (" to type '" + tname + "'")) + " failed"
      }
    {}
  
  template <typename T>
  static ConversionFailed makeFor(std::string const& key, std::string const& s)
    { return { key, s, typeid(T).name() }; }
  
}; // icarus::details::KeyValuesData::ConversionFailed()


// -----------------------------------------------------------------------------
struct icarus::details::KeyValuesData::ItemNotFound: public ErrorOnKey {
  
  ItemNotFound(std::string const& key): ErrorOnKey(key, "key not found") {}
  
}; // icarus::details::KeyValuesData::ItemNotFound()


// -----------------------------------------------------------------------------
struct icarus::details::KeyValuesData::ValueNotAvailable: public ErrorOnKey {
  
  ValueNotAvailable(std::string const& key, std::size_t index)
    : ErrorOnKey(key, "item value #" + std::to_string(index) + " not available")
    {}
  
}; // icarus::details::KeyValuesData::ValueNotAvailable()


// -----------------------------------------------------------------------------
struct icarus::details::KeyValuesData::MissingSize: public ErrorOnKey {
  
  MissingSize(std::string const& key)
    : ErrorOnKey
      (key, "is required to have a size as first value, but it has no values")
    {}
  
  MissingSize(std::string const& key, std::string const& valueStr)
    : ErrorOnKey(
      key,
      " first value '" + valueStr + "' can't be converted into a vector size"
      )
    {}
  
}; // icarus::details::KeyValuesData::MissingSize


// -----------------------------------------------------------------------------
struct icarus::details::KeyValuesData::WrongSize: public ErrorOnKey {
  
  WrongSize(std::string const& key, std::size_t expected, std::size_t actual)
    : ErrorOnKey(key,
      std::to_string(expected) + " values (except the size) were expected, "
      + std::to_string(actual) + " are present instead"
      )
    {}
  
}; // icarus::details::KeyValuesData::WrongSize


// -----------------------------------------------------------------------------
struct icarus::details::KeyedCSVparser::ParserError: public Error {
  
  ParserError(std::string msg): KeyValuesData::Error(std::move(msg)) {}
  
}; // icarus::details::KeyedCSVparser::ParseError


// -----------------------------------------------------------------------------
struct icarus::details::KeyedCSVparser::InvalidFormat: public ParserError {
  
  InvalidFormat(std::string const& msg): ParserError("Format error: " + msg) {}
  
}; // icarus::details::KeyedCSVparser::InvalidFormat


// -----------------------------------------------------------------------------
struct icarus::details::KeyedCSVparser::MissingValues: public ErrorOnKey {
  
  MissingValues(std::string const& key, unsigned int values)
    : ErrorOnKey(key,
      "data ended while expecting " + std::to_string(values) + " more values"
      )
    {}
  
}; // icarus::details::KeyedCSVparser::MissingValues


// -----------------------------------------------------------------------------
// --- Inline implementation
// -----------------------------------------------------------------------------
inline auto icarus::details::KeyedCSVparser::parse
  (std::string_view const& s) const -> ParsedData_t
  { ParsedData_t data; parse(s, data); return data; }


// -----------------------------------------------------------------------------
// --- Template implementation
// -----------------------------------------------------------------------------
template <typename T, typename Enable>
std::optional<T>
icarus::details::KeyValuesData::Item::StringConverter<T, Enable>::operator()
  (std::string const& valueStr) const
{
  using std::from_chars;
  T number {}; // useless initialization to avoid GCC complains
  char const *b = valueStr.data(), *e = b + valueStr.length();
  if constexpr(std::is_integral_v<T>) {
    return (from_chars(b, e, number, base).ptr == e)
      ? std::make_optional(number): std::nullopt;
  }
  else if constexpr(std::is_arithmetic_v<T>) {
#if defined(__GNUC__)
# if __GNUC__ < 11
    char* str_end;
    number = std::strtod(b, &str_end);
    return (str_end == e)? std::make_optional(number): std::nullopt;
# else
#   error("Code update required: GCC 11 finally supports `std::from_chars(float)`")
    // just remove all this #if stuff
# endif
#else
    return (from_chars(b, e, number).ptr == e)
      ? std::make_optional(number): std::nullopt;
#endif
  }
  else if constexpr(std::is_constructible_v<T, std::string>){
    return std::make_optional(T{ valueStr });
  }
  else return std::nullopt;
  
} // icarus::details::KeyValuesData::Item::StringConverter::operator()


// -----------------------------------------------------------------------------
template <typename T, typename Conv>
T icarus::details::KeyValuesData::Item::getAs
  (std::size_t index, Conv converter) const
{
  
  if (index >= values().size()) throw ValueNotAvailable(key(), index);
  
  auto const& valueStr = values()[index];
  auto const number = converter(valueStr);
  return number? *number: throw ConversionFailed::makeFor<T>(key(), valueStr);
  
} // icarus::details::KeyValuesData::Item::getAs<>()


// -----------------------------------------------------------------------------
template <typename T>
T icarus::details::KeyValuesData::Item::getAs(std::size_t index) const
  { return getAs<T>(index, StringConverter<T>{}); }


// -----------------------------------------------------------------------------
template <typename T, bool IgnoreFormatErrors /* = false */, typename Conv>
std::optional<T> icarus::details::KeyValuesData::Item::getOptionalAs
  (std::size_t index, Conv converter) const
{
  if (index < values().size()) return std::nullopt;

  auto const& valueStr = values()[index];
  auto const number = converter(valueStr);
  return (number || IgnoreFormatErrors)
    ? number: throw ConversionFailed::makeFor<T>(key(), valueStr);

} // icarus::details::KeyValuesData::Item::getOptionalAs()


// -----------------------------------------------------------------------------
template <typename T, bool IgnoreFormatErrors /* = false */>
std::optional<T> icarus::details::KeyValuesData::Item::getOptionalAs
  (std::size_t index) const
  { return getOptionalAs<T, IgnoreFormatErrors>(index, StringConverter<T>{}); }


// -----------------------------------------------------------------------------
template <typename T>
T icarus::details::KeyValuesData::Item::getNumber
  (std::size_t index, unsigned int base /* = 10 */) const
  { return getAs<T>(index, StringConverter<T>{ base }); }


// -----------------------------------------------------------------------------
template <typename T>
std::optional<T> icarus::details::KeyValuesData::Item::getOptionalNumber
  (std::size_t index, unsigned int base /* = 10 */) const
  { return getOptionalAs<T>(index, StringConverter<T>{ base }); }


// -----------------------------------------------------------------------------
template <typename T, typename Conv>
std::vector<T> icarus::details::KeyValuesData::Item::getVector
  (Conv converter /* = {} */) const
{
  return
    convertVector<T>(values().begin(), values().end(), std::move(converter));
}


// -----------------------------------------------------------------------------
template <typename T, typename Conv>
std::vector<T> icarus::details::KeyValuesData::Item::getSizedVector
  (Conv converter /* = {} */) const
{
  
  if (values().empty()) throw MissingSize(key());
  
  std::size_t const n = getNumber<std::size_t>(0U);
  if (n != values().size() - 1)
    throw WrongSize(key(), n, values().size() - 1);
  
  return convertVector<T>
    (std::next(values().begin()), values().end(), std::move(converter));
} // icarus::details::KeyValuesData::Item::getSizedVector()


// -----------------------------------------------------------------------------
template <typename T, typename Iter, typename Conv>
std::vector<T> icarus::details::KeyValuesData::Item::convertVector
  (Iter begin, Iter end, Conv converter) const
{
  std::vector<T> data;
  data.reserve(std::distance(begin, end));
  Iter it = begin;
  while (it != end) {
    std::string const& valueStr = *it++;
    std::optional const number = converter(valueStr);
    if (!number) throw ConversionFailed::makeFor<T>(key(), valueStr);
    data.push_back(*number);
  }
  return data;
} // icarus::details::KeyValuesData::Item::getVector()


// -----------------------------------------------------------------------------
template <typename T>
std::optional<T> icarus::details::KeyValuesData::Item::convertStringInto
  (std::string const& valueStr)
  { return StringConverter<T>{}.convert(valueStr); }


// -----------------------------------------------------------------------------
inline decltype(auto) icarus::details::KeyValuesData::items() const noexcept
  { return fItems; }


// -----------------------------------------------------------------------------
template <typename BIter, typename EIter>
auto icarus::details::KeyedCSVparser::parse(BIter b, EIter e) const
  -> ParsedData_t
  { return parse(std::string_view{ &*b, std::distance(b, e) }); }


// -----------------------------------------------------------------------------


#endif // ICARUSCODE_DECODE_DECODERTOOLS_DETAILS_KEYEDCVSPARSER_H
