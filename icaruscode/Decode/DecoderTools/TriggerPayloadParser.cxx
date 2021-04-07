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
    else {
      throw cet::exception("TriggerPayloadParser")
        << "Unsupported token: '" << word << "'\n";
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
  std::string token;
  while (sstr >> token) {
    if ((token.length() == 1) && (token.front() == sep)) continue;
    if (!token.empty() && token.back() == sep)
      token.erase(token.length() - 1);
    while (!token.empty() && std::isspace(token.back()))
      token.erase(token.length() - 1);
    if (!token.empty() && token.front() == sep)
      token.erase(0);
    while (!token.empty() && std::isspace(token.front()))
      token.erase(0);
    return std::optional(std::move(token));
  } // while
  return std::nullopt;
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
