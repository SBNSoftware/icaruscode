/**
 * @file   icaruscode/PMT/Algorithms/KeyValueParser.h
 * @brief  Simple parser for "key: value" text.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   May 13, 2022
 * @see    icaruscode/PMT/Algorithms/KeyValueParser.cxx
 */

#ifndef ICARUSCODE_PMT_ALGORITHMS_KEYVALUEPARSER_H
#define ICARUSCODE_PMT_ALGORITHMS_KEYVALUEPARSER_H

// ICARUS libraries
#include "icaruscode/PMT/Algorithms/ParsingToolkit.h"
#include "icaruscode/Decode/DecoderTools/details/KeyValuesData.h"

// C++ standard libraries
#include <iosfwd> // std::istream
#include <string_view>
#include <vector>
#include <string>


// -----------------------------------------------------------------------------
namespace icarus::details { class KeyValueParser; }
/**
 * @class icarus::details::KeyValueParser
 * @brief Parser to fill a `KeyValuesData` structure out of a character buffer.
 * 
 * The parser processes an input stream or string.
 * Parsing produces a data structure of type `icarus::KeyValuesData`; it is not
 * possible to reconstruct the input from that structure, since quotations,
 * escapes, formatting and comments are lost.
 * 
 * The supported format is:
 * * one key/value pair per line
 * * all keys and values are stripped of trailing and ending blank characters
 * * key is everything before the key separator (by default, `:` to
 *   create/overwrite, `:+` to create/append values)
 * * array values are separated by spaces
 * * quotation supports both single and double quotes
 * * comments start with a word beginning with `#` and extend to the end of the
 *   line
 * * escape character (backslash by default) drops the special meaning of
 *   comment, key and quote symbols
 * 
 * 
 * 
 * 
 * 
 * Example of format
 * ------------------
 * 
 * This is the text from the unit test:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * 
 * # SPR test input file
 * 
 * Description: "
 * This is a test for the key-values parser with default settings.
 * It is expected to be used to describe the Single Photoelectron Response.
 * "
 * 
 * Contact: Gianluca Petrillo (petrillo@slac.stanford.edu)
 * 
 * Gain: 9.7e6  # from amplitude 4 mV
 * Tick: '2 ns'
 * Samples: 0.0 1.0 2.5 \
 *          4.5 3.0 2.5
 * Samples:+1.8 1.6 1.2 0.8 0.8 0.7 0.7 0.6
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * 
 */
class icarus::details::KeyValueParser {
  
    public:
  
  /// Base of all errors by KeyValueParser.
  using Error = icarus::KeyValuesData::Error;
  struct ParserError; ///< Generic error: base of all errors by KeyValueParser.
  
  using ParsedData_t = icarus::KeyValuesData; ///< Type of returned data.
  
  /// Parameters of the format.
  struct FormatParams_t {
    std::string newKey { ":" }; ///< Sequence starting the values.
    std::string addKey { ":+" }; ///< Sequence appending the values.
  };
  
  static const FormatParams_t DefaultFormatParameters;
  
  /// Creates a parser with the specified parsing parameters.
  KeyValueParser(
    FormatParams_t formatParams = DefaultFormatParameters,
    icarus::ParsingToolkit::Params_t parserParams
      = icarus::ParsingToolkit::DefaultParameters
    );
  

  //@{
  /// Parses the `stream` and returns a data structure with the content.
  ParsedData_t parse(std::istream& stream) const;
  ParsedData_t parse(std::istream&& stream) const { return parse(stream); }
  
  ParsedData_t operator() (std::istream& stream) const { return parse(stream); }
  //@}
  
  /// Runs the parser on a string.
  ParsedData_t parse(std::string const& s) const;
  
  
    private:
  
  enum class keyType { unsupported, create, add };

  icarus::ParsingToolkit fPTK; ///< Parsing toolkit (and its parameters).
  FormatParams_t fFmt; ///< Parser format parameters.
  
  std::vector<std::string> fKeys; ///< Sorted keys cache.
  
  
  /**
   * @brief Modifies `tokens` placing the key/value separator in its own token.
   * @param[in,out] tokens list of tokens
   * @return the type of separator found
   * 
   * On success, the first token is the key, the second the separator that was
   * found, and the rest are values. The list of `tokens` is modified in place.
   * In case of failure, where the separator is not found where expected,
   * `tokens` are unchanged.
   * 
   * The separator is expected to be either unquoted and unescaped in the first
   * of the tokens, or otherwise at the beginning of the second token.
   */
  keyType highlightSeparator(std::vector<std::string_view>& tokens) const;
  
  
  /// Returns the type of `key`.
  template <typename Key>
  keyType keySepType(Key const& key) const;
  
}; // icarus::details::KeyValueParser


// -----------------------------------------------------------------------------
// ---  Exception class definitions
// -----------------------------------------------------------------------------
struct icarus::details::KeyValueParser::ParserError: public Error {
  
  ParserError
    (unsigned int iLine, std::string const& line, std::string const& msg);
  
}; // icarus::details::KeyValueParser::Error()


// -----------------------------------------------------------------------------

#endif // ICARUSCODE_PMT_ALGORITHMS_KEYVALUEPARSER_H
