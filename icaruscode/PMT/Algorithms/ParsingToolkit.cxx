/**
 * @file   icaruscode/PMT/Algorithms/ParsingToolkit.cxx
 * @brief  Simple text parsing utilities.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   May 13, 2022
 * @see    icaruscode/PMT/Algorithms/ParsingToolkit.h
 */

// library header
#include "icaruscode/PMT/Algorithms/ParsingToolkit.h"

// C/C++ standard libraries
#include <algorithm> // std::sort(), std::unique()
#include <string>
#include <string_view>


// -----------------------------------------------------------------------------
// ---  icarus::ParsingToolkit
// -----------------------------------------------------------------------------
icarus::ParsingToolkit::Params_t const
icarus::ParsingToolkit::DefaultParameters;


// -----------------------------------------------------------------------------
std::pair<std::string, unsigned int> icarus::ParsingToolkit::readMultiline
  (std::istream& in) const
{
  
  std::string fullLine;
  std::string openQuoteLine;
  unsigned int nLines = 0U;
  while (in) {
    
    std::string line;
    std::getline(in, line, in.widen(fParams.EOL));
    bool const isEOF = in.eof();
    if (!isEOF || !line.empty()) ++nLines;
    openQuoteLine.append(line);
    
    if (isQuotationUnclosed(make_view(openQuoteLine))) {
      if (isCharacterEscaped(line.begin(), line.end())) {
        fullLine.append(openQuoteLine);
        throw Error{ "Parser error: escaped end-of-line inside a quotation:\n"
          + fullLine + "\n" };
      }
      // if the newline is quoted, it's preserved
      if (!isEOF) openQuoteLine += fParams.EOL;
      continue;
    }
    fullLine.append(openQuoteLine);
    openQuoteLine.clear();
    
    if (!isCharacterEscaped(fullLine.begin(), fullLine.end())) break;
    
    fullLine.pop_back(); // remove the escape character
    
  } // while
  fullLine.append(openQuoteLine); // usually empty
  
  return { std::move(fullLine), nLines };
} // icarus::ParsingToolkit::readMultiline()


// -----------------------------------------------------------------------------
auto icarus::ParsingToolkit::findQuotationStart(std::string_view sv) const
  -> std::pair<std::string_view, QuotSpec_t const*>
{
  
  while (!sv.empty()) {
    
    // look for a character that could start a quotation opening
    std::size_t const startPos = sv.find_first_of(fQuoteStarts);
    
    // no such character found:
    if (startPos == std::string_view::npos) break;
    
    // if the character is escaped, this is not a quotation opening:
    if (isCharacterEscaped(sv.begin(), sv.begin() + startPos)) {
      sv.remove_prefix(std::min(startPos + 1, sv.length()));
      continue;
    }
    
    sv.remove_prefix(std::min(startPos, sv.length()));
    
    // try all the opening quotes
    // (may be optimized by grouping them by first character)
    for (auto const& qSpec: fParams.quotes) {
      
//       if (sv.starts_with(qSpec.first)) return { sv, &qSpec }; // C++20
      if (sv.compare(0, qSpec.first.length(), qSpec.first) == 0)
        return { sv, &qSpec };
      
    } // for quotes
    
    // nope, just a character; remove it and keep looking
    sv.remove_prefix(1);
    
  } // while sv
  
  return { make_view(sv.end(), sv.end()), nullptr };
} // icarus::ParsingToolkit::findQuotationStart()


// -----------------------------------------------------------------------------
std::string_view icarus::ParsingToolkit::findQuotationEnd
  (std::string_view sv, std::string const& quotEnd) const
{
  while (!sv.empty()) {
    
    std::size_t const pos = sv.find(quotEnd);
    if (pos == std::string_view::npos) break;
    
    if (!isCharacterEscaped(sv.begin(), sv.begin() + pos)) {
      sv.remove_prefix(pos);
      return sv;
    }
    
    sv.remove_prefix(pos + 1);
    
  } // while
  
  return make_view(sv.end(), sv.end());
} // icarus::ParsingToolkit::findQuotationEnd()


// -----------------------------------------------------------------------------
bool icarus::ParsingToolkit::isQuotationUnclosed(std::string_view sv) const {
  
  while (!sv.empty()) {
    
    auto [ qsv, qptr ] = findQuotationStart(sv);
    if (!qptr) return false;
    
    qsv.remove_prefix(qptr->first.length()); // remove the opening quote
    
    qsv = findQuotationEnd(qsv, qptr->second);
    if (qsv.empty()) return true;
    
    qsv.remove_prefix(qptr->second.length());
    sv = qsv;
  } // while
  
  return false;
} // icarus::ParsingToolkit::isQuotationUnclosed()


// -----------------------------------------------------------------------------
auto icarus::ParsingToolkit::splitOn(std::string_view sv, std::string_view sep)
  -> SplitView_t
{
  return
    { make_view(sv.begin(), sep.begin()), sep, make_view(sep.end(), sv.end()) };
} // icarus::ParsingToolkit::splitOn()


// -----------------------------------------------------------------------------
std::string icarus::ParsingToolkit::removeWordEscapes(std::string&& s) const {
  
  // replace in place
  std::string::const_iterator iSrc = s.begin(), send = s.end();
  std::string::iterator iDest = s.begin();
  
  // if the last character is an escape, it's kept
  while (iSrc != send) {
    char const ch = *iSrc++;
    *iDest++ = (isEscape(ch) && (iSrc != send))? *iSrc++: ch;
  } // while
  
  s.erase(iDest, send);
  return std::move(s);
} // icarus::ParsingToolkit::removeWordEscapes()


// -----------------------------------------------------------------------------
std::string icarus::ParsingToolkit::removeWordQuotations(std::string&& s) const
{
  std::string_view sv = make_view(s);
  std::string::iterator iDest = s.begin();
  
  while (!sv.empty()) {
    
    // find the next quotation
    auto const [ fromQ, qptr ] = findQuotationStart(sv);
    
    // copy the material until the next quotation
    iDest = std::copy(sv.begin(), fromQ.begin(), iDest);
    sv = fromQ;
    
    if (!qptr) break; // if there is no quotation, we are done
    
    sv.remove_prefix(qptr->first.length()); // skip the quotation start
    
    // find the end of quotation
    std::string_view const afterQ = findQuotationEnd(sv, qptr->second);
    
    if (afterQ.empty()) { // begin of quotation, but no end: no good
      // leave the "begin of quotation" as is
      iDest = std::copy(fromQ.begin(), fromQ.end(), iDest);
      sv.remove_prefix(sv.length()); // note: quote start was already removed
      break;
    }

    // copy the quoted material
    iDest = std::copy(sv.begin(), afterQ.begin(), iDest);
    sv = afterQ;
    
    sv.remove_prefix(qptr->second.length()); // skip the quotation end
    
  } // while
  
  assert(sv.empty());
  
  s.erase(iDest, s.end());
  return std::move(s);
} // icarus::ParsingToolkit::removeWordQuotations()


// -----------------------------------------------------------------------------
void icarus::ParsingToolkit::adoptParams(Params_t params) {
  
  fParams = std::move(params);
  
  // sort the quotations by length
  auto const byOpeningLength = [](QuotSpec_t const& a, QuotSpec_t const& b)
    {
      std::size_t const al = a.first.length(), bl = b.first.length();
      return (al != bl)? (al > bl): (a < b);
    };
  std::sort(fParams.quotes.begin(), fParams.quotes.end(), byOpeningLength);
  
  // collect the first character of each of the opening quotes
  // (sorted and with no duplicates)
  for (QuotSpec_t const& quotSpec: fParams.quotes)
    fQuoteStarts += quotSpec.first.front();
  std::sort(fQuoteStarts.begin(), fQuoteStarts.end());
  fQuoteStarts.erase
    (std::unique(fQuoteStarts.begin(), fQuoteStarts.end()), fQuoteStarts.end());
  
} // icarus::ParsingToolkit::adoptParams()


// -----------------------------------------------------------------------------
