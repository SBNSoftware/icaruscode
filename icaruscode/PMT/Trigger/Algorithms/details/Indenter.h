/**
 * @file   icaruscode/PMT/Trigger/Algorithms/details/Indenter.h
 * @brief  Simple helper to track the indentation level.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   November 12, 2025
 * @see    icaruscode/PMT/Trigger/Algorithms/Indenter.cxx
 * 
 */

#ifndef ICARUSCODE_PMT_TRIGGER_ALGORITHMS_DETAILS_INDENTER_H
#define ICARUSCODE_PMT_TRIGGER_ALGORITHMS_DETAILS_INDENTER_H


// C/C++ standard libraries
#include <iosfwd> // std::ostream
#include <string>
#include <utility> // std::move()


// -----------------------------------------------------------------------------
// forward declarations
namespace icarus::trigger::details {
  class Indenter;
  
  /// Inserts a new line using the specified indenter.
  std::ostream& operator<< (std::ostream& out, Indenter& indenter);
}

// -----------------------------------------------------------------------------
/**
 * @brief Simple helper to apply the indentation.
 * 
 * This object tracks whether it's the first indentation request or not,
 * and applies/returns an indentation string accordingly.
 * There are two indentation strings supported, the one on the first line
 * and the one on all the others.
 * 
 * The first line is continued, i.e. no line break is inserted. For all others,
 * before the indentation a line break (``'\n'``) is inserted (and it is also
 * included in the value returned by `indent()`).
 * 
 * There are several ways to exercise the functionality of this object;
 * the most convenient is arguably via the insertion operator free function:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * template <typename Coll>
 * void dumpVector(std::ostream& out, Coll const& data) {
 *   using std::size, std::begin, std::end;
 * 
 *   icarus::trigger::details::Indenter newline{ "   ", "" };
 *   
 *   std::size_t const n = size(data);
 *   if (n == 0) out << newline << "(empty)";
 *   else {
 *     out << newline << n << " elements:";
 *     std::size_t i = 0;
 *     for (auto const& elem: data) {
 *       out << newline << "[#" << (i++) << "] " << elem;
 *     }
 *   }
 *   // deliberately do not terminate the last line
 * }
 * 
 * int main() {
 *   std::cout << "Data content ";
 *   dumpVector(std::cout, std::valarray{ 1, 2, 3 });
 *   std::cout << std::endl;
 * }
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * will output:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Data content >>> 3 elements:
 *    [#0] 1
 *    [#1] 2
 *    [#2] 3
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * 
 */
class icarus::trigger::details::Indenter {
  std::string fIndent;
  std::string fFirstIndent;
  unsigned int fCount = 0;
  
    public:
  
  /// Constructor: the same string for first and following indentations.
  Indenter(std::string indent)
    : fIndent{ "\n" + indent }, fFirstIndent{ std::move(indent) }
    {}
  
  /// Constructor: two separate indentations for first and following lines.
  Indenter(std::string indent, std::string firstIndent)
    : fIndent{ "\n" + indent }, fFirstIndent{ std::move(firstIndent) } {}
  
  /// Registers the start of a new line, returns whether it's not the first.
  bool useLine() { return bool(fCount++); }
  
  /// String to start a line (e.g. `std::cout << indenter() << ...).
  /// @see nextLine()
  std::string const& operator() () { return nextLine(); }
  
  /// String to start a line (e.g. `std::cout << indenter.nextLine() << ...).
  std::string const& nextLine() { return useLine()? fIndent: fFirstIndent; }
  
  /// Inserts a new line into a stream.
  std::ostream& insertNewLine(std::ostream& out);
  
  /// Returns the internally used indentation string (starts with `'\\n'`).
  std::string const& indent() const { return fIndent; }
  
  /// Returns the internally used first line indentation string.
  std::string const& firstIndent() const { return fFirstIndent; }
  
  /**
   * @brief Returns a new indenter with expanded indentation strings.
   * @param addIndent string to be added to the general indenting
   * @param first first line indentation (replacing, not extending)
   * @return a new indenter object
   * 
   * This function helps to carry an existing indentation on.
   * The line counter of the returned indenter is reset, so it will use
   * `first` string as first indentation string, and then an extended version
   * of the current one form the others.
   * With some care, this allows some fine tuning of the output.
   * 
   * Example of usage:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * void Algorithm::dumpConfig(std::ostream& out, Indenter& indent) const {
   *   
   *   out << indent << "Algorithm configuration:"
   *     << indent << " * input algorithm configuration";
   *   fInputAlgorithm->dumpConfig(out, indent.nested("  ", ": "));
   *   out
   *     << indent << " * output algorithm configuration:"
   *     << indent << "   ";
   *   fOutputAlgorithm->dumpConfig(out, indent.nested("  "));
   * }
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * Here an hypothetical `Algorithm` object holding two sub-algorithms for
   * input and output calls their `dumpConfig()` functions (same signature as
   * `Algorithm::dumpConfig` above) and adds two spaces to the indentation,
   * with the input algorithm using `": "` as first line indentation,
   * and the output algorithm leaving the first line indentation empty.
   * The tricky part here is the first line. Nested indenter will start over
   * and use its own first line indentation, so the input algorithm will have
   * its first line starting from the current output line and adding a colon and
   * a space, plus the text of its first line. The output algorithm instead
   * will start from a new line under the `o` of `output` (no first line
   * indentation) and then its next lines will go two spaces further (regular
   * indentation).
   * 
   * Note that in this example `Algorithm::dumpConfig()` can leave the last
   * line partially filled (i.e. it does not break it before returning).
   * This is standard practice to give control of the following layout to the
   * caller of `dumpConfig()`.
   */
  Indenter nested(std::string addIndent, std::string first = "") const;
  
    protected:
  
  /// Default constructor; for internal use only.
  Indenter() = default;
  
}; // icarus::trigger::details::Indenter


// -----------------------------------------------------------------------------

#endif // ICARUSCODE_PMT_TRIGGER_ALGORITHMS_DETAILS_INDENTER_H
