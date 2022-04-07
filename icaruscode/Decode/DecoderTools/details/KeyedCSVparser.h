/**
 * @file   icaruscode/Decode/DecoderTools/details/KeyedCSVparser.h
 * @brief  Simple parser for comma-separated text.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   May 9, 2021
 * @see    icaruscode/Decode/DecoderTools/details/KeyedCSVparser.cxx
 */

#ifndef ICARUSCODE_DECODE_DECODERTOOLS_DETAILS_KEYEDCVSPARSER_H
#define ICARUSCODE_DECODE_DECODERTOOLS_DETAILS_KEYEDCVSPARSER_H

// ICARUS libraries
#include "icaruscode/Decode/DecoderTools/details/KeyValuesData.h"

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
namespace icarus::details { class KeyedCSVparser; }
/**
 * @class icarus::details::KeyedCSVparser
 * @brief Parser to fill a `KeyValuesData` structure out of a character buffer.
 * 
 * It currently supports only single-line buffer.
 */
class icarus::details::KeyedCSVparser {
  
    public:
  
  using ParsedData_t = icarus::KeyValuesData;
  
  /// Base of all errors by KeyedCSVparser.
  using Error = icarus::KeyValuesData::Error;
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
struct icarus::details::KeyedCSVparser::ParserError: public Error {
  
  ParserError(std::string msg): Error(std::move(msg)) {}
  
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
template <typename BIter, typename EIter>
auto icarus::details::KeyedCSVparser::parse(BIter b, EIter e) const
  -> ParsedData_t
  { return parse(std::string_view{ &*b, std::distance(b, e) }); }


// -----------------------------------------------------------------------------


#endif // ICARUSCODE_DECODE_DECODERTOOLS_DETAILS_KEYEDCVSPARSER_H
