/**
 * @file icaruscode/Decode/DecoderTools/details/KeyedCSVparser.cxx
 * @brief  Simple parser for comma-separated text (implementation).
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   May 9, 2021
 * @see icaruscode/Decode/DecoderTools/details/KeyedCSVparser.h
 */

// library header
#include "icaruscode/Decode/DecoderTools/details/KeyedCSVparser.h"

// C++ standard libraries
#include <ostream>
#include <cassert>
#include <cctype> // std::isspace()


// -----------------------------------------------------------------------------
// ---  icarus::details::KeyValuesData

// -----------------------------------------------------------------------------
auto icarus::details::KeyValuesData::makeItem(std::string key) -> Item& {
  
  if (hasItem(key)) throw DuplicateKey{ key };
  
  fItems.emplace_back(std::move(key));
  
  return fItems.back();
} // icarus::details::KeyValuesData::makeItem()


// -----------------------------------------------------------------------------
auto icarus::details::KeyValuesData::findItem
  (std::string const& key) const noexcept -> Item const*
{
  for (auto const& item: fItems) if (key == item.key()) return &item;
  return nullptr;
} // icarus::details::KeyValuesData::findItem()


// -----------------------------------------------------------------------------
auto icarus::details::KeyValuesData::getItem(std::string const& key) const
  -> Item const&
{
  if (auto item = findItem(key); item) return *item;
  throw ItemNotFound(key);
} // icarus::details::KeyValuesData::getItem()


// -----------------------------------------------------------------------------
bool icarus::details::KeyValuesData::hasItem
  (std::string const& key) const noexcept
{
  return findItem(key);
} // icarus::details::KeyValuesData::hasItem()


// -----------------------------------------------------------------------------
bool icarus::details::KeyValuesData::empty() const noexcept
  { return fItems.empty(); }


// -----------------------------------------------------------------------------
std::size_t icarus::details::KeyValuesData::size() const noexcept
  { return fItems.size(); }


// -----------------------------------------------------------------------------
std::ostream& icarus::details::operator<<
  (std::ostream& out, KeyValuesData const& data)
{
  out << data.size() << " items";
  if (!data.empty()) {
    out << ':';
    for (auto const& item: data.items()) {
      out << "\n  '" << item.key() << "' (" << item.values().size() << ")";
      if (item.values().empty()) continue;
      out << ':';
      for (auto const& value: item.values()) out << " '" << value << '\'';
    } // for
  } // if has data
  return out << "\n";
} // icarus::details::operator<< (KeyValuesData)


// -----------------------------------------------------------------------------
// ---  icarus::details::KeyedCSVparser
// -----------------------------------------------------------------------------
void icarus::details::KeyedCSVparser::parse
  (std::string_view const& s, ParsedData_t& data) const
{
  
  auto stream = s;
  
  ParsedData_t::Item* currentItem = nullptr;
  
  // this many tokens will be assigned to the current key:
  int forcedValues = -1; // 0 would force the first entry to be a key
  
  while (!stream.empty()) {
    
    auto const token = extractToken(stream);
    
    std::string tokenStr { cbegin(token), cend(token) };
    
    bool bKey = false;
    do {
      
      // if there are values pending, this is not a key, period;
      // if all required values have been assigned, the next token is a key.
      if (forcedValues >= 0) {
        bKey = (forcedValues == 0); // if no more forced values, next is key
        --forcedValues;
        // if we know the token is a key, we still need to check for matching
        // patterns to assign required values
        if (!bKey) break;
      }
      
      // the token may still be a key (if `bKey` is true, it is for sure: we can
      // decide that a non-key (!bKey) is actually a key, but not the opposite)
      for (auto const& [ pattern, values ]: fPatterns) {
        if (!std::regex_match(begin(token), end(token), pattern)) continue;
        bKey = true; // matching a pattern implies this is a key
        std::string const& key = tokenStr;
        // how many values to expect:
        switch (values) {
          case FixedSize: // read the next token immediately as fixed size
            {
              if (stream.empty()) throw MissingSize(key);
              
              auto const sizeToken = peekToken(stream);
              if (empty(sizeToken)) throw MissingSize(key);
              
              // the value is loaded in `forcedValues` and already excludes
              // the size token just read
              char const *b = begin(sizeToken), *e = end(sizeToken);
              if (std::from_chars(b, e, forcedValues).ptr != e)
                throw MissingSize(key, std::string{ sizeToken });
              
              ++forcedValues; // the size will be forced in the values anyway
              
            } // FixedSize
            break;
          case DynamicSize:
            // nothing to do, the normal algorithm rules will follow
            break;
          default:
            forcedValues = values;
            break;
        } // switch
        break;
      } // for pattern
      if (bKey) break;
      
      // let the "standard" pattern decide
      bKey = isKey(token);
      
    } while (false);
    
    if (bKey) currentItem = &(data.makeItem(std::move(tokenStr)));
    else {
      if (!currentItem) {
        throw InvalidFormat(
         "values started without a key ('" + tokenStr + "' is not a valid key)."
         );
      }
      currentItem->addValue(std::move(tokenStr));
    }
    
  } // while
  
  if (forcedValues > 0) {
    assert(currentItem);
    throw MissingValues(currentItem->key(), forcedValues);
  }
  
} // icarus::KeyedCSVparser::parse()


// -----------------------------------------------------------------------------
auto icarus::details::KeyedCSVparser::addPatterns
  (std::initializer_list<std::pair<std::regex, unsigned int>> patterns)
  -> KeyedCSVparser&
{
  for (auto& pattern: patterns) fPatterns.emplace_back(std::move(pattern));
  return *this;
} // icarus::details::KeyedCSVparser::addPatterns()


// -----------------------------------------------------------------------------
auto icarus::details::KeyedCSVparser::addPatterns
  (std::initializer_list<std::pair<std::string, unsigned int>> patterns)
  -> KeyedCSVparser&
{
  for (auto& pattern: patterns)
    fPatterns.emplace_back(std::regex{ pattern.first }, pattern.second);
  return *this;
} // icarus::details::KeyedCSVparser::addPatterns()


// -----------------------------------------------------------------------------
auto icarus::details::KeyedCSVparser::parse
  (std::string const& s) const -> ParsedData_t
  { return parse(std::string_view{ s.data(), s.size() }); }


// -----------------------------------------------------------------------------
std::size_t icarus::details::KeyedCSVparser::findTokenLength
  (Buffer_t const& buffer) const noexcept
{
  
  auto const start = cbegin(buffer), bend = cend(buffer);
  auto finish = start;
  while (finish != bend) {
    if (*finish == fSep) break;
    ++finish;
  } // for
  
  return std::distance(start, finish);
} // icarus::details::KeyedCSVparser::findTokenLength()


// -----------------------------------------------------------------------------
auto icarus::details::KeyedCSVparser::peekToken
  (Buffer_t const& buffer) const noexcept -> SubBuffer_t
{
  return strip({ cbegin(buffer), findTokenLength(buffer) });
} // icarus::details::KeyedCSVparser::peekToken()


// -----------------------------------------------------------------------------
auto icarus::details::KeyedCSVparser::extractToken
  (Buffer_t& buffer) const noexcept -> SubBuffer_t
{
#if 1
  auto const start = cbegin(buffer), bend = cend(buffer);
  std::size_t const length = findTokenLength(buffer);
  moveBufferHead(buffer, length + ((start + length == bend)? 0: 1));
  return strip({ start, length });
#else
  
  auto const start = cbegin(buffer), bend = cend(buffer);
  auto finish = start;
  while (finish != bend) {
    if (*finish == fSep) break;
    ++finish;
  } // for
  
  // update the start of the buffer
  std::size_t const tokenLength = std::distance(start, finish);
  moveBufferHead(buffer, tokenLength + ((finish == bend)? 0: 1));
  
  return strip({ start, tokenLength });
#endif
} // icarus::details::KeyedCSVparser::extractToken()


// -----------------------------------------------------------------------------
bool icarus::details::KeyedCSVparser::isKey
  (SubBuffer_t const& buffer) const noexcept
{
  
  return !buffer.empty() && std::isalpha(buffer.front());
  
} // icarus::details::KeyedCSVparser::isKey()


// -----------------------------------------------------------------------------
template <typename String>
auto icarus::details::KeyedCSVparser::makeBuffer(String const& s) noexcept
  -> Buffer_t
  { return { data(s), size(s) }; } // C++20: use begin/end constructor



// -----------------------------------------------------------------------------
auto icarus::details::KeyedCSVparser::moveBufferHead
  (Buffer_t& buffer, std::size_t size) noexcept -> Buffer_t&
{
  
  size = std::min(size, buffer.size());
  return buffer = { buffer.data() + size, buffer.size() - size };
  
} // details::KeyedCSVparser::eatBufferHead()


// -----------------------------------------------------------------------------
auto icarus::details::KeyedCSVparser::strip(SubBuffer_t s) noexcept
  -> SubBuffer_t
  { return stripRight(stripLeft(stripRightChars<'\n', '\r', '\0'>(s))); }


// -----------------------------------------------------------------------------
auto icarus::details::KeyedCSVparser::stripLeft(SubBuffer_t s) noexcept
  -> SubBuffer_t
{
  
  while (!s.empty()) {
    if (!std::isspace(s.front())) break;
    s.remove_prefix(1);
  }
  return s;
  
} // icarus::details::KeyedCSVparser::stripLeft()


// -----------------------------------------------------------------------------
auto icarus::details::KeyedCSVparser::stripRight(SubBuffer_t s) noexcept
  -> SubBuffer_t
{
  
  while (!s.empty()) {
    if (!std::isspace(s.back())) break;
    s.remove_suffix(1);
  }
  return s;
  
} // icarus::details::KeyedCSVparser::stripRight()


// -----------------------------------------------------------------------------
auto icarus::details::KeyedCSVparser::stripRightChar
  (SubBuffer_t s, char c) noexcept -> SubBuffer_t
{
  
  while (!s.empty()) {
    if (s.back() != c) break;
    s.remove_suffix(1);
  }
  return s;
  
} // icarus::details::KeyedCSVparser::stripRightChar()


// ----------------------------------------------------------------------------
template <char... Chars>
auto icarus::details::KeyedCSVparser::stripRightChars
  (SubBuffer_t s) noexcept -> SubBuffer_t
{
  while (true) {
    auto ns = s;
    for (char c: { Chars... }) ns = stripRightChar(ns, c);
    if (ns == s) return ns;
    s = ns;
  } // while(true)
  
} // icarus::details::KeyedCSVparser::stripRightChars()


// -----------------------------------------------------------------------------

