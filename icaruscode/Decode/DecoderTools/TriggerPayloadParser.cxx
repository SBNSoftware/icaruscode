/**
 * @file   icaruscode/Decode/DecoderTools/TriggerPayloadParser.cxx
 * @brief  Provides a parser for trigger packet content (implementation file).
 * @author Jacob Zettlemoyer,
 *         Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   April 6, 2021
 * @see    icaruscode/Decode/DecoderTools/TriggerPayloadParser.h
 */


// Library header
#include "icaruscode/Decode/DecoderTools/TriggerPayloadParser.h"

// C/C++ standard library
#include <vector>
#include <cctype> // std::isspace()

// -----------------------------------------------------------------------------
// ---  TriggerPayloadParser implementation
// -----------------------------------------------------------------------------
auto daq::TriggerPayloadParser::parse(std::string payload) -> TriggerData_t {
  /*
   * The parser uses the payload as a stream, reading word by word.
   * It is also a rudimental state machine: depending on which (key)word it
   * finds, it will decide what to do next.
   * 
   * This guarantees a lot of flexibility, allowing the order of the keywords
   * to be not mandated, but makes it harder to verify that all required data
   * has been, sooner or later, eventually read.
   * A check concludes the parsing to enforce the validity of the parsed data.
   * 
   * The data is returned in a "proprietary" format that is probably not ideal
   * for the downstream code; it is recommended to copy it in a more suitable
   * format if needed.
   */
  std::istringstream buffer { std::move(payload) };
  std::optional<std::string> token;
  std::vector<std::string> tokens;
  TriggerData_t triggerData;
  while (token = nextToken(buffer)) {
    std::string const& word = token.value();
    // doing exception handling right takes so much time...
    if (word == "Local_TS1")
      triggerData.Local_TS1.emplace(parseTSdata(buffer, word));
    else if (word == "WR_TS1")
      triggerData.WR_TS1.emplace(parseTSdata(buffer, word));
    else if (word == "Gate ID")
      parseGateID(buffer); // TODO store it
    else if (word == "Gate Type")
      parseGateType(buffer); // TODO store it
    else {
#if 1
      // FIXME this is for TEMPORARY support for recent runs (581x)
      // with some buffer overflow problem
      // we are currently not depending on message facility... so I even skip the warning
      // mf::LogWarning("TriggerPayloadParser")
      //   << "Unsupported token: '" << word << "'";
      continue;
#else
      throw cet::exception("TriggerPayloadParser")
        << "Unsupported token: '" << word << "'\n";
#endif // 1?
    }
  }
  validate(triggerData); // throws if incomplete/invalid/inconsistent
  return triggerData;
} // daq::TriggerPayloadParser::parse()


// -----------------------------------------------------------------------------
void daq::TriggerPayloadParser::validate(TriggerData_t const& data) {
  std::vector<std::string> errors;
  if (!data.Local_TS1.has_value())
    errors.push_back("Local_TS1 time stamp not found");
  if (!data.WR_TS1.has_value())
    errors.push_back("WR_TS1 time stamp not found");
  if (errors.empty()) return;
  cet::exception e { "TriggerPayloadParser" };
  e << "Parser found " << errors.size() << " errors:";
  for (std::string const& msg: errors) e << "\n - " << msg;
  throw e << '\n';
} // daq::TriggerPayloadParser::validate()


// -----------------------------------------------------------------------------
std::optional<std::string> daq::TriggerPayloadParser::nextToken
  (std::istream& sstr)
{
  if (wasteSpaces(sstr).fail()) return std::nullopt; // got to EOF (and back)

  // there are only two ways to complete a token:
  // a separator character, or the end of file
  char c;
  std::string token;
  while (sstr.get(c)) { // spaces are merged verbatim into the token
    if (issep(c)) break;
    token += c;
  } // while
  
  // remove trailing spaces from the token
  while (!token.empty() && std::isspace(token.back())) token.pop_back();
  
  // if we reached the end of file and the token is still empty,
  // i.e. if the current token was terminated by EOF rather than a separator,
  // we don't emit the (empty) token: an empty token can be only be made
  // by closing it explicitly with a separator
  return (sstr.eof() && token.empty())
    ? std::nullopt: std::optional(std::move(token));
} // daq::TriggerPayloadParser::nextToken()


// -----------------------------------------------------------------------------
auto daq::TriggerPayloadParser::parseTSdata
  (std::istringstream& sstr, std::string const& TSname) -> TSdata_t
{
  TSdata_t data;
  try {
    data.eventNo = nextRequiredValue<int>(sstr);
  }
  catch (cet::exception const& e) {
    throw cet::exception("TriggerPayloadParser", "", e)
      << "Error trying to read event number for TS1 '" << TSname << "'\n";
  }
  try {
    data.timeStampHigh = nextRequiredValue<int>(sstr);
  }
  catch (cet::exception const& e) {
    throw cet::exception("TriggerPayloadParser", "", e)
      << "Error trying to read high timestamp (seconds) for TS1 '"
      << TSname << "'\n";
  }
  try {
    data.timeStampLow = nextRequiredValue<int>(sstr);
  }
  catch (cet::exception const& e) {
    throw cet::exception("TriggerPayloadParser", "", e)
      << "Error trying to read low timestamp (nanoseconds) for TS1 '"
      << TSname << "'\n";
  }
  return data;
} // daq::TriggerPayloadParser::parseTSdata()


// -----------------------------------------------------------------------------
auto daq::TriggerPayloadParser::parseGateID
  (std::istringstream& sstr) -> GateID_t
{
  return { nextRequiredValue<int>(sstr) };
} // daq::TriggerPayloadParser::parseGateID()


// -----------------------------------------------------------------------------
auto daq::TriggerPayloadParser::parseGateType
  (std::istringstream& sstr) -> GateType_t
{
  GateType_t gateType;
  gateType.A = nextRequiredValue<int>(sstr);
//   gateType.B = nextRequiredValue<int>(sstr);
  return gateType;
} // daq::TriggerPayloadParser::parseGateType()


// -----------------------------------------------------------------------------
/* constexpr */ bool daq::TriggerPayloadParser::issep(char c) // C++20: constexpr
  { return std::find(seps.begin(), seps.end(), c) != seps.end(); }


// -----------------------------------------------------------------------------
std::istream& daq::TriggerPayloadParser::wasteSpaces(std::istream& is) {
  constexpr auto EoF = std::istream::traits_type::eof();
  do {
    auto const ch = is.get();
    if (ch == EoF) return is;
    if (!std::isspace(ch)) return is.unget();
  } while (is);
  throw std::logic_error("daq::TriggerPayloadParser::wasteSpaces() unexpected");
} // daq::TriggerPayloadParser::wasteSpaces()


// -----------------------------------------------------------------------------
