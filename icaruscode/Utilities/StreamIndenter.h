/**
 * @file icaruscode/Utilities/StreamIndenter.h
 * @brief Utility to have simple indentation in a stream.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date March 22, 2022
 * 
 * This is a header-only library.
 */

#ifndef ICARUSCODE_UTILITIES_STREAMINDENTER_H
#define ICARUSCODE_UTILITIES_STREAMINDENTER_H


// C/C++ standard libraries
#include <utility> // std::move(), std::forward()
#include <sstream>
#include <string>


namespace util {
  
  /**
   * @brief Stream modifier that makes it "indented".
   * 
   * The intended use of this class is to be "inserted" into an output stream
   * object in order to gain a simple indentation:
   *  * on the first insertion after this modifier, `firstIndent` will be used
   *    to indent the stream;
   *  * on the following insertions, `indent` will be used instead.
   * 
   * Example:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * std::cout << util::addIndent("> ", "") << "First line: not indented."
   *   << "\nSecond line: indented by '> '."
   *   << "\nThird line: indented by '> '."
   *   ;
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * 
   * Technically, after the first insertion (`std::cout << addIndent(...)`)
   * the returned object is not `std::cout` any more, but a wrapper around it
   * which intercepts all the insertion operations (`operator<<`).
   * 
   * The wrapper is designed to steal the stream object if it is a temporary one
   * (e.g. `std::ofstream{ "test.txt" } << addIndent("> ") << ...`) and to
   * reference to it otherwise. Unless the wrapper object is somehow saved,
   * it is destroyed (together with the stream if it was "stolen") as soon as
   * the statement falls out of scope. Therefore, in this two-lines example:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * std::cout << util::addIndent("> ", "");
   * std::cout << "First line: not indented."
   *   << "\nSecond line: also not indented."
   *   ;
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * the first line does effectively nothing (the wrapper is created, ready to
   * indent, but it is then immediately destroyed) while the second line
   * performs normal `std::cout` output.
   * 
   * @note The "first line" is actually the first insertion into the wrapped
   *       stream. For example,
   *       `std::cout << "A!" << util::addIndent("> ", "$ ") << "First line.";`
   *       will produce `A!$ First line.`.
   * 
   * 
   * @note It is in principle possible to create modifiers with a special
   *       behaviour to interact directly with the wrapper (for example, change
   *       the indentation level or indentation string, reset to the first line
   *       indentation, ...). Such modifiers are not implemented so far.
   */
  struct addIndent {
    std::string firstIndent, indent;
    
    addIndent(std::string indent, std::string firstIndent)
      : firstIndent{ std::move(firstIndent) }, indent{ std::move(indent) }
      {}
    addIndent(std::string const& indent)
      : addIndent{ indent, indent } {}
    
  }; // addIndent
  
  
  //@{
  /**
   * @brief Creates an indented stream wrapper.
   * @tparam Stream the type of stream being wrapped
   * @param out the stream to wrap
   * @param indent string inserted at the beginning of each new line
   * @param firstIndent string inserted at the beginning of the first line
   * @return an indented stream wrapper
   * @see `addIndent()`
   * 
   * The use of the wrapper stream is explained in `addIndent()`.
   * 
   * This function wraps a stream `out`, stealing it if it's a temporary,
   * and returns the wrapping object. This object can be used for indented
   * output to `out` until it is destroyed (if `out` is referenced as opposed
   * to stolen, `out` itself needs to stay valid).
   * 
   * If not specified, `firstIndent` is assigned to be the same as `indent`.
   * 
   * The equivalent two-line example of `addIndent()` becomes:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * auto out = util::makeIndented(std::cout, "> ", "$ ");
   * out << "First line: indented with '$ '."
   *   << "\nSecond line: indented with '> '."
   *   ;
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */
  template <typename Stream>
  auto makeIndented(Stream&& out, std::string indent, std::string firstIndent);
  
  template <typename Stream>
  auto makeIndented(Stream&& out, std::string const& indent);
  //@}
  
  
  /// Helper function for `addIndent` (@see `addIndent`).
  template <typename Stream>
  auto operator<< (Stream&& out, addIndent&& adder);
  
} // namespace util


namespace util::details {

  // ---------------------------------------------------------------------------
  /**
   * @brief Stream wapper.
   */
  template <typename Stream>
  class IndentAdder {
    
    Stream&& fOut;
    std::string const fFirstIndent;
    std::string const fIndent;
    
    std::stringstream fSStr; // don't even try to be thread-safe
    std::string const* fCurrentIndent { &fFirstIndent };
    
      public:
    IndentAdder(Stream&& out, std::string indent, std::string firstIndent)
      : fOut{ std::forward<Stream>(out) }
      , fFirstIndent{ std::move(firstIndent) }, fIndent{ std::move(indent) }
      {}
    IndentAdder(Stream& out, std::string const& indent = "")
      : IndentAdder{ std::forward<Stream>(out), indent, indent } {}
    
    template <typename T>
    IndentAdder& operator<< (T&& value)
      { fSStr << std::forward<T>(value); streamOut(); return *this; }
    
      private:
    
    void streamOut()
      {
        bool newLine = true;
        char ch;
        while(fSStr.get(ch)) {
          if (newLine) { fOut << *fCurrentIndent; fCurrentIndent = &fIndent; }
          fOut << ch;
          newLine = (ch == '\n');
        } // while
        fSStr.clear();
      }
    
  }; // class IndentAdder
  
  
  // ---------------------------------------------------------------------------
  
} // namespace util::details


// -----------------------------------------------------------------------------
template <typename Stream>
auto util::makeIndented
  (Stream&& out, std::string indent, std::string firstIndent)
{
  return details::IndentAdder
    { std::forward<Stream>(out), std::move(indent), std::move(firstIndent) };
}


// -----------------------------------------------------------------------------
template <typename Stream>
auto util::makeIndented(Stream&& out, std::string const& indent)
  { return makeIndented(out, indent, indent); }


// -----------------------------------------------------------------------------
template <typename Stream>
auto util::operator<< (Stream&& out, addIndent&& adder) {
  return details::IndentAdder{
    std::forward<Stream>(out),
    std::move(adder.indent), std::move(adder.firstIndent)
    };
} // operator<< (Stream, addIndent)


// -----------------------------------------------------------------------------


#endif // ICARUSCODE_UTILITIES_STREAMINDENTER_H
