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
#include <cctype> // std::isspace()


// -----------------------------------------------------------------------------
// ---  icarus::details::KeyedCSVparser
// -----------------------------------------------------------------------------
void icarus::details::KeyedCSVparser::parse
  (std::string_view const& s, ParsedData_t& data) const
{
  
  auto stream = s;
  
  ParsedData_t::Item* currentItem = nullptr;
  
  while (!stream.empty()) {
    
    auto const token = extractToken(stream);
    
    std::string tokenStr { cbegin(token), cend(token) };
    
    if (isKey(token)) currentItem = &(data.makeItem(std::move(tokenStr)));
    else {
      if (!currentItem) {
        throw InvalidFormat(
         "values started without a key ('" + tokenStr + "' is not a valid key)."
         );
      }
      currentItem->addValue(std::move(tokenStr));
    }
    
  } // while
  
} // icarus::KeyedCSVparser::parse()


// -----------------------------------------------------------------------------
auto icarus::details::KeyedCSVparser::parse
  (std::string const& s) const -> ParsedData_t
  { return parse(std::string_view{ s.data(), s.size() }); }


// -----------------------------------------------------------------------------
auto icarus::details::KeyedCSVparser::extractToken
  (Buffer_t& buffer) const noexcept -> SubBuffer_t
{
  
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

