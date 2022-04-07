/**
 * @file   icaruscode/PMT/Algorithms/KeyValueParser.cxx
 * @brief  Simple parser for "key: value" text.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   May 13, 2022
 * @see    icaruscode/PMT/Algorithms/KeyValueParser.h
 */

// library header
#include "icaruscode/PMT/Algorithms/KeyValueParser.h"

// C++ standard libraries
#include <algorithm> // std::copy(), std::sort()
#include <istream>
#include <sstream>
#include <initializer_list>
#include <iterator> // std::advance(), std::next(), std::move_iterator()
#include <string>
#include <utility> // std::move()
#include <cstdint> // std::size_t


// -----------------------------------------------------------------------------
namespace {
  
  /// Returns a copy of `keys` collection sorted by length then
  /// lexicographically.
  template <typename Keys>
  std::vector<std::string> sortKeys(std::initializer_list<Keys> keys) {
    using std::size;
    std::vector<std::string> sorted;
    sorted.reserve(size(keys));
    std::copy(begin(keys), end(keys), back_inserter(sorted));
    
    auto const byOpeningLength = [](auto const& a, auto const& b)
      {
        std::size_t const al = a.length(), bl = b.length();
        return (al != bl)? (al > bl): (a < b);
      };
    std::sort(sorted.begin(), sorted.end(), byOpeningLength);
    return sorted;
  } // sortKeys()
  
} // local namespace


// -----------------------------------------------------------------------------
// ---  icarus::details::KeyValueParser
// -----------------------------------------------------------------------------
const icarus::details::KeyValueParser::FormatParams_t
icarus::details::KeyValueParser::DefaultFormatParameters;


// -----------------------------------------------------------------------------
  /// Creates a parser with the specified parsing parameters.
icarus::details::KeyValueParser::KeyValueParser(
  FormatParams_t formatParams /* = DefaultFormatParameters */,
  icarus::ParsingToolkit::Params_t parserParams
    /* = icarus::ParsingToolkit::DefaultParameters */
  )
  : fPTK{ std::move(parserParams) }
  , fFmt{ std::move(formatParams) }
  , fKeys{ sortKeys({ fFmt.newKey, fFmt.addKey }) }
{}


// -----------------------------------------------------------------------------
auto icarus::details::KeyValueParser::parse(std::istream& s) const
  -> ParsedData_t
{
  
  ParsedData_t data;
  
  unsigned int iSrcLine = 0U; // count of lines in the source
  while (s) {
    
    auto [ line, nMultiLines ] = fPTK.readMultiline(s);
    iSrcLine += nMultiLines;
    
    if (line.empty()) continue;
    
    std::vector<std::string_view> tokens = fPTK.splitWords(line);
    
    fPTK.removeCommentLine(tokens);
    if (tokens.empty()) continue;
    
    keyType const kType = highlightSeparator(tokens);
    if (kType == keyType::unsupported) {
      auto iKey = fKeys.cbegin(), kend = fKeys.cend();
      std::string l { "'" + *iKey + "'" };
      while (++iKey != kend) l += ", '" + *iKey + "'";
      throw ParserError
        { iSrcLine - nMultiLines + 1, line, "no key separator (" + l + ")" };
    }
    assert(tokens.size() >= 2U);
    
    std::vector<std::string> words
      = fPTK.removeEscapes(fPTK.removeQuotations(tokens));
    
    auto iWord = words.begin();
    ParsedData_t::Item& item = data.makeOrFetchItem(std::move(*iWord));
    std::advance(iWord, 2); // skip to after the separator
    
    switch (kType) {
      case keyType::create:
        item.values.clear();
        [[fallthrough]];
      case keyType::add:
        item.addValues
          (std::move_iterator{ iWord }, std::move_iterator{ words.end() });
        break;
      default:
        throw ParserError{ iSrcLine - nMultiLines + 1, line,
          "LOGIC ERROR: '" + std::string{ tokens[1] }
          + "' should have been a key separator"
          };
    } // switch
    
  } // while
  
  return data;
  
} // icarus::KeyValueParser::parse()


// -----------------------------------------------------------------------------
auto icarus::details::KeyValueParser::parse
  (std::string const& s) const -> ParsedData_t
  { return parse(std::istringstream{ s }); }


// -----------------------------------------------------------------------------
auto icarus::details::KeyValueParser::highlightSeparator
  (std::vector<std::string_view>& tokens) const -> keyType
{
  
  // need to find the separator in the unquoted, unescaped parts of first token,
  // or the separator must be the start of the second token
  
  if (tokens.empty()) return keyType::unsupported;
  
  //
  // separator in the first word
  //
  std::string_view const firstKey
    = fPTK.findFirstUnquoted(tokens.front(), fKeys);
  
  if (!firstKey.empty()) {
    auto const [ pre, sep, post ] = fPTK.splitOn(tokens.front(), firstKey);
    // if pre is empty, it's still considered a (empty) key (questionable...)
    // if post is empty, we omit it
    tokens.front() = sep;
    if (!post.empty()) tokens.insert(std::next(tokens.begin()), post);
    tokens.insert(tokens.begin(), pre);
    return keySepType(sep);
  }
  
  //
  // separator is the second word start
  //
  if (tokens.size() < 2U)
    return keyType::unsupported; // ah, actually: no, it wasn't.
  
  for (std::string const& sep: fKeys) {
    
    if (tokens[1].compare(0, sep.length(), sep) != 0) continue;
    
    // so now we have [0] key [1] sep[+first value] [2...] values;
    // we just need to insert the "first value", if any, as its own token;
    // we also take care to leave the views pointing to the token rather than to
    // a member of 'fKeys'
    if (!tokens[1].empty()) {
      std::string_view firstValue = tokens[1];
      firstValue.remove_prefix(sep.length());
      tokens.insert(std::next(tokens.begin(), 2), firstValue);
      tokens[1].remove_suffix(firstValue.length());
    }
    
    return keySepType(sep);
  } // for
  
  //
  // no separator at all
  //
  return keyType::unsupported;
  
} // icarus::details::KeyValueParser::highlightSeparator()


// -----------------------------------------------------------------------------
template <typename Key>
auto icarus::details::KeyValueParser::keySepType(Key const& key) const
  -> keyType
{
  if (key == fFmt.newKey) return keyType::create;
  if (key == fFmt.addKey) return keyType::add;
  return keyType::unsupported;
} // icarus::details::KeyValueParser::keySepType()


// -----------------------------------------------------------------------------
icarus::details::KeyValueParser::ParserError::ParserError
  (unsigned int iLine, std::string const& line, std::string const& msg)
  : Error{
    "KeyValueParser::ParserError on line " + std::to_string(iLine)
    + " ('" + line + "'): " + msg
    }
  {}

// -----------------------------------------------------------------------------

