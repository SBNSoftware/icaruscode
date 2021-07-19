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
#include <stdexcept> // std::runtime_error
#include <utility> // std::move()
#include <charconv> // std::from_chars()
#include <cstddef> // std::size_t


// -----------------------------------------------------------------------------
namespace icarus::details {
  class KeyedCSVparser;
  struct KeyValuesData;
  
  std::ostream& operator<< (std::ostream& out, KeyValuesData const& data);
} // namespace icarus::details


// -----------------------------------------------------------------------------
struct icarus::details::KeyValuesData {
  
  struct Error;
  struct DuplicateKey;
  struct ConversionFailed;
  struct ItemNotFound;
  struct ValueNotAvailable;
  
  /// Representation of a single item of data: a key and several values.
  struct Item {
    std::string key;
    std::vector<std::string> values;
    
    Item(std::string key): key(std::move(key)) {}
    
    Item& addValue(std::string value)
      { values.push_back(std::move(value)); return *this; }
    
    std::size_t nValues() const noexcept { return values.size(); }
    
    bool operator< (Item const& other) const noexcept
      { return key < other.key; }
    
    /// Converts the value to `T`, throws on failure.
    template <typename T>
    T getNumber(std::size_t index) const;
    
    /// Converts the value to `T`, throws on conversion failures
    /// unless `ignoreFormatErrors` is `true`.
    template <typename T>
    std::optional<T> getOptionalNumber
      (std::size_t index, bool ignoreFormatErrors = false) const;
    
      private:
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
 */
class icarus::details::KeyedCSVparser {
  
    public:
  
  using ParsedData_t = icarus::details::KeyValuesData;
  
  using Error = KeyValuesData::Error; ///< Base of all errors by KeyedCSVparser.
  struct ParserError; ///< Generic error: base of all errors by KeyedCSVparser.
  struct InvalidFormat; ///< Parsing format is not understood.
  
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
  
    private:
  using Buffer_t = std::string_view;
  using SubBuffer_t = std::string_view;
  
  char const fSep = ','; ///< Character used as token separator.
  
  
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
struct icarus::details::KeyValuesData::DuplicateKey: public Error {
  
  DuplicateKey(std::string const& msg)
    : Error("KeyValuesData::DuplicateKey: '" + msg + '\'')
    {}
  
}; // icarus::details::KeyValuesData::DuplicateKey()


// -----------------------------------------------------------------------------
struct icarus::details::KeyValuesData::ConversionFailed: public Error {
  
  ConversionFailed(std::string const& s, std::string const& tname = "")
    : Error{
      "conversion of '" + s + "'"
       + (tname.empty()? "": (" to type '" + tname + "'")) + " failed"
      }
    {}
  
  template <typename T>
  static ConversionFailed makeFor(std::string const& s)
    { return { s, typeid(T).name() }; }
  
}; // icarus::details::KeyValuesData::ConversionFailed()


// -----------------------------------------------------------------------------
struct icarus::details::KeyValuesData::ItemNotFound: public Error {
  
  ItemNotFound(std::string const& key): Error("item not found: '" + key + '\'')
    {}
  
}; // icarus::details::KeyValuesData::ItemNotFound()


// -----------------------------------------------------------------------------
struct icarus::details::KeyValuesData::ValueNotAvailable: public Error {
  
  ValueNotAvailable(std::size_t index)
    : Error("item value #" + std::to_string(index) + " not available")
    {}
  
}; // icarus::details::KeyValuesData::ValueNotAvailable()



// -----------------------------------------------------------------------------
struct icarus::details::KeyedCSVparser::ParserError: public Error {
  
  ParserError(std::string msg): KeyValuesData::Error(std::move(msg)) {}
  
}; // icarus::details::KeyedCSVparser::Error()


// -----------------------------------------------------------------------------
struct icarus::details::KeyedCSVparser::InvalidFormat: public ParserError {
  
  InvalidFormat(std::string const& msg): ParserError("Format error: " + msg) {}
  
}; // icarus::details::KeyedCSVparser::InvalidFormat()


// -----------------------------------------------------------------------------
// --- Inline implementation
// -----------------------------------------------------------------------------
inline auto icarus::details::KeyedCSVparser::parse
  (std::string_view const& s) const -> ParsedData_t
  { ParsedData_t data; parse(s, data); return data; }


// -----------------------------------------------------------------------------
// --- Template implementation
// -----------------------------------------------------------------------------
template <typename T>
T icarus::details::KeyValuesData::Item::getNumber(std::size_t index) const {
  
  if (index >= values.size()) throw ValueNotAvailable(index);
  
  auto const& valueStr = values[index];
  auto const number = convertStringInto<T>(valueStr);
  return number? number.value(): throw ConversionFailed::makeFor<T>(valueStr);
  
} // icarus::details::KeyValuesData::Item::getNumber<>()


// -----------------------------------------------------------------------------
template <typename T>
std::optional<T> icarus::details::KeyValuesData::Item::getOptionalNumber
  (std::size_t index, bool ignoreFormatErrors /* = false */) const
{
  if (index < values.size()) return std::nullopt;

  auto const& valueStr = values[index];
  auto const number = convertStringInto<T>(valueStr);
  return (number || ignoreFormatErrors)
    ? number: throw ConversionFailed::makeFor<T>(valueStr);

} // icarus::details::KeyValuesData::Item::getOptionalNumber()


// -----------------------------------------------------------------------------
template <typename T>
std::optional<T> icarus::details::KeyValuesData::Item::convertStringInto
  (std::string const& valueStr)
{
  T number {}; // useless initialization to avoid GCC complains
  char const *b = valueStr.data(), *e = b + valueStr.length();
  return (std::from_chars(b, e, number).ptr == e)
    ? std::make_optional(number): std::nullopt;
} // icarus::details::KeyValuesData::Item::convertStringInto()


// -----------------------------------------------------------------------------
template <typename BIter, typename EIter>
auto icarus::details::KeyedCSVparser::parse(BIter b, EIter e) const
  -> ParsedData_t
  { return parse(std::string_view{ &*b, std::distance(b, e) }); }


// -----------------------------------------------------------------------------


#endif // ICARUSCODE_DECODE_DECODERTOOLS_DETAILS_KEYEDCVSPARSER_H
