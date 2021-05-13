/**
 * @file   icaruscode/Decode/DecoderTools/TriggerPayloadParser.h
 * @brief  Provides a parser for trigger packet content.
 * @author Jacob Zettlemoyer,
 *         Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   April 6, 2021
 * @see    icaruscode/Decode/DecoderTools/TriggerPayloadParser.cxx
 */

#ifndef ICARUSCODE_DECODE_DECODERTOOLS_TRIGGERPAYLOADPARSER_H
#define ICARUSCODE_DECODE_DECODERTOOLS_TRIGGERPAYLOADPARSER_H


// framework libraries
#include "larcorealg/CoreUtils/DebugUtils.h" // lar::debug::demangle()
#include "cetlib_except/exception.h"

// C/C++ standard library
#include <string>
#include <istream>
#include <sstream> // std::istringstream
#include <optional>
#include <utility> // std::move()


// -----------------------------------------------------------------------------
namespace daq { class TriggerPayloadParser; }
/**
 * @brief Class to parse the trigger readout payload.
 * 
 * This class allows parsing a text to extract trigger data information.
 * Example of usage:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * std::string payloadAsText
 *   = reinterpret_cast<char const*>(fragment.dataBeginBytes());
 * std::replace(payloadAsText.begin(), payloadAsText.end(), '\x0D', '\n');
 * daq::TriggerPayloadParser parser;
 * daq::TriggerPayloadParser::TriggerData_t triggerData = parser(payloadAsText);
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * 
 * The example assumes that `fragment` is a `art::Fragment` originating from
 * the trigger readout, and its data ("payload") is text, terminating with a
 * null character.
 * The returned value, of type `TriggerData_t`, contains the parsed information.
 * 
 * The supported format is documented in the `parse()` method.
 */
class daq::TriggerPayloadParser {
  
    public:
  
  using GateID_t = int; ///< Type for gate ID storage. TODO
  
  /// Type for gate type storage. TODO
  struct GateType_t {
    static constexpr int noType = std::numeric_limits<int>::max();
    int A = noType;
    int B = noType;
  }; // GateType_t
  
  /// Data for a single time stamp.
  struct TSdata_t {
    int eventNo = -1;
    int timeStampHigh = 0;
    int timeStampLow = 0;
    double time() const
      { return double(timeStampHigh) + double(timeStampLow) * 1e-9; }
    long long int nanoseconds() const
      { return (timeStampHigh * 1'000'000'000LL) + timeStampLow; }
  }; // struct TSdata_t
  
  /// The complete interpreted trigger data.
  struct TriggerData_t {
    
    std::optional<TSdata_t> Local_TS1; ///< Local time stamp information.
    std::optional<TSdata_t> WR_TS1;    ///< White Rabbit time stamp information.
    
  }; // struct TriggerData_t
  
  
  /**
   * @brief Extracts trigger data from the specified payload.
   * @param payload the string to read data from
   * @return trigger data in a "proprietary" format
   * 
   * The payload is expected to be text, a comma-separated list of words.
   * Some words are keywords that specify the meaning of the following content.
   * Actually, the words are effectively separated by spaces, and commas
   * (at the beginning or end of the word, and isolated) are ignored. Therefore,
   * `"a, b and c"` are four words (`"a"`, `"b"`, `"and"` and `"c"`), and
   * `"d ,e,f, that's it" are also four (`"d"`, `"e,f"`, `"that's"` and `"it"`).
   * 
   * The data is returned in a "proprietary" format that is probably not ideal
   * for the downstream code; it is recommended to copy it in a more suitable
   * format if needed.
   * 
   * Expected format
   * ----------------
   * 
   * A "section" is a sequence of words introduced by a keyword that identifies
   * which section that is and how to interpret the following words.
   * Sections may in principle be nested, if that is how the sections are
   * designed.
   * 
   * The following sections can be generally in any order:
   * * `Local_TS1` keyword for local timestamp information, followed by:
   *     * integer: an event number
   *     * integer: the time stamp, in seconds, truncated
   *     * integer: the rest of the time stamp, in nanoseconds
   * * `WR_TS1` keyword for White Rabbit timestamp information, in the same
   *       format as for the local timestamp `Local_TS1`
   * 
   */
  TriggerData_t parse(std::string payload);
  
  /// An alias of `parse()`.
  TriggerData_t operator() (std::string const& payload)
    { return parse(payload); }
  
  
    private:
  
  /// Separator used to split words (actually, not really).
  static constexpr std::array<char, 3U> seps = { ',', '\r', '\n' };
  
  /// Throws an exception if `data` is not complete.
  static void validate(TriggerData_t const& data);
  
  /// Reads one word from the specified stream of comma-separated words.
  static std::optional<std::string> nextToken(std::istream& sstr);
  
  /// Fills `value` with the `nextToken()`, interpreted as type `T`.
  /// @return whether a value was available
  /// @throw cet::exception if token is available but interpretation fails
  template <typename T>
  static bool nextValue(std::istream& sstr, T& value);
  
  
  /// Returns the `nextToken()`, interpreted as type `T`, if available.
  /// @throw cet::exception if token is available but interpretation fails
  template <typename T>
  static std::optional<T> nextValue(std::istream& sstr);
  
  /// Returns the `nextToken()`, interpreted as type `T`.
  /// @throw cet::exception if token is available but interpretation fails
  /// @throw cet::exception if token is not available
  template <typename T>
  static T nextRequiredValue(std::istream& sstr);
  
  
  /// Parses and returns a `TSdata_t` record (`TSname` is used for messages).
  TSdata_t parseTSdata(std::istringstream& sstr, std::string const& TSname);
  
  
  /// Extracts information about the gate ID from the stream `sstr`.
  static GateID_t parseGateID(std::istringstream& sstr);
  
  /// Extracts information about the gate type from the stream `sstr`.
  static GateType_t parseGateType(std::istringstream& sstr);
  
  /// Extracts all spaces from `is`, returns `is` itself (may end up in EOF).
  static std::istream& wasteSpaces(std::istream& is);
  
  /// Returns whether c` is any of the separators in `seps`.
  static /* constexpr */ bool issep(char c); // C++20: constexpr

}; // class daq::TriggerPayloadParser



// -----------------------------------------------------------------------------
// ---  template implementation
// -----------------------------------------------------------------------------
template <typename T>
bool daq::TriggerPayloadParser::nextValue(std::istream& sstr, T& value) {
  // wish I could use std::string_view... C++23, apparently
  std::optional<std::string> valueStr = nextToken(sstr);
  if (!valueStr) return false;
  std::istringstream tokensstr { valueStr.value() };
  tokensstr >> value;
  if (tokensstr.fail()) {
    throw cet::exception("TriggerPayloadParser")
      << "Failed to convert '" << valueStr.value() << "' to "
      << lar::debug::demangle<T>() << "\n";
  }
  tokensstr.peek();
  if (!tokensstr.eof()) {
    throw cet::exception("TriggerPayloadParser")
      << "Spurious characters (starting with: '" << tokensstr.peek()
      << "') after converting '" << valueStr.value() << "' to "
      << lar::debug::demangle<T>() << "( " << value << ")\n";
  }
  return true;
} // daq::TriggerPayloadParser::nextValue(sstream, T)


// -----------------------------------------------------------------------------
template <typename T>
std::optional<T> nextValue(std::istream& sstr) {
  T value;
  return nextValue(sstr, value)
    ? std::optional{ std::move(value) }: std::nullopt;
} // daq::TriggerPayloadParser::nextValue(sstream)


// -----------------------------------------------------------------------------
template <typename T>
T daq::TriggerPayloadParser::nextRequiredValue(std::istream& sstr) {
  T value;
  if (nextValue(sstr, value)) return value;
  throw cet::exception("TriggerPayloadParser") 
    << "No further data available.\n";
} // daq::TriggerPayloadParser::nextRequiredValue()


// -----------------------------------------------------------------------------


#endif // ICARUSCODE_DECODE_DECODERTOOLS_TRIGGERPAYLOADPARSER_H
