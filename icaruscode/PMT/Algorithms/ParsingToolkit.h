/**
 * @file   icaruscode/PMT/Algorithms/ParsingToolkit.h
 * @brief  Simple text parsing utilities.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   May 13, 2022
 * @see    icaruscode/PMT/Algorithms/ParsingToolkit.cxx
 */

#ifndef ICARUSCODE_PMT_ALGORITHMS_PARSINGTOOLKIT_H
#define ICARUSCODE_PMT_ALGORITHMS_PARSINGTOOLKIT_H

// C/C++ standard libraries
#include <algorithm> // std::count_if()
#include <istream>
#include <stdexcept> // std::runtime_error
#include <vector>
#include <initializer_list>
#include <string>
#include <string_view>
#include <utility> // std::pair, std::move()
#include <cctype> // std::isblank()
#include <cassert>


// -----------------------------------------------------------------------------
namespace icarus { class ParsingToolkit; }
/**
 * @brief Utilities for text parsing.
 * 
 * This "class" is a glorified namespace with some configuration inside.
 * 
 * 
 * Quotation
 * ----------
 * 
 * A quoted string is the content in between an opening quoting sequence and
 * the matching closing sequence. Each sequence may be any string, including
 * but not limited to a one-character long string. Escaping the first character
 * of an opening or closing quotation string will turn it in common string data
 * carrying no quotation meaning.
 * 
 * 
 * Escaping rules
 * ---------------
 * 
 * Any single character following the escape character is "escaped".
 * The escaped characters lose their standard function and are replaced by a
 * substitute character. For example, escaping the first character of a opening
 * quotation makes that a standard character. An escaped escape character is
 * always replaced by the character itself, without its escape function.
 * 
 * 
 */
struct icarus::ParsingToolkit {
  
  /// Base type for errors in the toolkit.
  class Error;
  
  /// Record of a split token: pre-separator, separator and post-separator.
  struct SplitView_t { std::string_view pre, sep, post; };
  
  /// Specification of quotation: opening and closing.
  using QuotSpec_t = std::pair<std::string, std::string>;
  
  /// All parsing parameters.
  struct Params_t {
    
    
    char escape { '\\' }; ///< Escape character.
    
    std::string comment { "#" }; ///< Word introducing a comment.
    
    char EOL { '\n' };
    
    /// List of matching start and end of quote.
    std::vector<QuotSpec_t> quotes {
      QuotSpec_t{ R"(")", R"(")" },
      QuotSpec_t{ R"(')", R"(')" }
      };
    
  }; // Params_t
  
  // Adapter converting argument of functions like `std::isblank()` properly.
  template <int (*CCTF)(int)>
  struct CCTypeAdapter {
    template <typename Ch>
    constexpr bool operator() (Ch c) const noexcept
      { return CCTF(static_cast<unsigned char>(c)); }
  }; // CCTypeAdapter
  
  /// Adapter for determining if a character is a blank (see `std::isblank()`).
  static constexpr CCTypeAdapter<&std::isblank> isBlank {};
  
  
  // --- BEGIN --- Initialization ----------------------------------------------
  
  static Params_t const DefaultParameters; /// Default parsing parameters.
  
  /// Creates a parser with the default parsing parameters.
  ParsingToolkit() { adoptParams(DefaultParameters); }
  
  /// Creates a parser with the specified parsing parameters.
  ParsingToolkit(Params_t params) { adoptParams(std::move(params)); }
  
  // --- END ----- Initialization ----------------------------------------------
  
  
  // --- BEGIN --- Query -------------------------------------------------------
  
  /// Returns the current parameters of parsing.
  Params_t const& params() const noexcept { return fParams; }
  
  // --- END ----- Query -------------------------------------------------------
  
  
  // --- BEGIN --- Input -------------------------------------------------------
  /// @name Input
  /// @{
  
  /**
   * @brief Returns a single line of text from the input stream.
   * @param in the input stream
   * @return the string read, and the number of lines read
   * @throw Error on fatal parsing errors
   * 
   * This function reads entire lines from `in`, where a line is defined as
   * in `std::getline()`. If the line ends with an unescaped escape character,
   * another line is read and appended (the escape character is dropped).
   * The return value is the merged string with no end-of-line characters,
   * and the number of lines read.
   * If there is no string to be read, it returns an empty string and `0U`.
   * 
   * ### Special behaviour
   * 
   * * If the line ends while a quotation is still open, the next line is also
   *   merged, and the line break is kept; to merge quoted lines without
   *   preserving the line break character, end the quote on the first line,
   *   immediately break the line escaping it, and then next line should
   *   immediately start with opening a quotation.
   * * If the line ends while a quotation is still open, it is a parsing error
   *   to have the line break character escaped (an exception will be thrown)
   *   merged, and the line break is kept.
   * * If the file ends while a quotation is still open, the line is preserved
   *   as such.
   */
  std::pair<std::string, unsigned int> readMultiline(std::istream& in) const;
  
  /// @}
  // --- END ----- Input -------------------------------------------------------
  
  // --- BEGIN --- Tokenization ------------------------------------------------
  /// @name Tokenization
  /// @{
  /**
   * @brief Splits a string into words.
   * @tparam Delim type of delimiter functor
   * @param s the string to be split
   * @param isDelimiter (default: `isblank()`) determines if a character is a
   *                    word delimiter
   * @return a sequence of views, one per word
   * 
   * The splitter algorithm defines a word separator as a sequence of one or
   * more unescaped, unquoted delimiter characters, where a delimiter is a
   * character `ch` for which `isDelimiter(ch)` is `true`.
   * 
   * Note that this function does not change the content of the data, and in
   * particular it does not remove escaping nor quoting (although it interprets
   * both).
   * 
   * A character used as delimiter can appear in a word only if escaped or
   * within quotation. Contiguous non-delimiter elements of a string, including
   * quoted strings, belong to the same word (for example, `a" and "b` is a
   * single word when delimitation is by blank characters).
   * An empty word can be introduced only in quotations (e.g. `""`).
   * 
   * The `Delim` type is a functor so that `isDelimiter(ch)` returns something
   * convertible to `bool`, `true` if the `ch` character should be considered
   * a delimiter. Note that no context is provided for the answer, so the
   * use of each character as delimiter is fixed, and modified only by the
   * hard-coded quotation and escaping rules.
   * 
   * The first characters of quotation starts and the escape characters must not
   * be classified as delimiters, or the algorithm will give wrong results.
   */
  template <typename Delim>
  std::vector<std::string_view> splitWords
    (std::string const& s, Delim isDelimiter) const;
  
  /// Helper version of `splitWords(std::string const&, Delim)`.
  std::vector<std::string_view> splitWords(std::string const& s) const
    { return splitWords(s, isBlank); }
  
  
  /**
   * @brief Finds the first word starting with a comment marker.
   * @tparam Iter type of iterator to the words
   * @param beginWord iterator to the first word to consider
   * @param endWord iterator past the lasy word to consider
   * @return an iterator to the comment word, or `endWord` if not found
   * 
   * The original list is modified, the word starting with a comment marker and
   * all the following ones are removed.
   */
  template <typename Iter>
  Iter findCommentWord(Iter beginWord, Iter endWord) const;
  
  
  /**
   * @brief Removes all the words from the one starting with a comment marker.
   * @param words list of words
   * 
   * The original list is modified, the word starting with a comment marker and
   * all the following ones are removed.
   */
  template <typename WordType>
  void removeCommentLine(std::vector<WordType>& words) const
    { words.erase(findCommentWord(words.begin(), words.end()), words.end()); }
  
  
  /**
   * @brief Finds the start of the next quotation in `sv`.
   * @param sv the buffer to look the quotation start into
   * @return a subview of `sv` starting from the quotation found, empty if none
   */
  std::pair<std::string_view, QuotSpec_t const*> findQuotationStart
    (std::string_view sv) const;
  
  /**
   * @brief Finds the quotation end in `sv`.
   * @param sv the buffer to look the quotation end into
   * @param quotEnd the quotation end to be searched
   * @return a view of `sv` from the quotation end, included, empty if not found
   * 
   * Note that `sv` should not include the quotation start.
   */
  std::string_view findQuotationEnd
    (std::string_view sv, std::string const& quotEnd) const;
  
  /// Returns if the sequence `sv` has unclosed quotation at its end.
  bool isQuotationUnclosed(std::string_view sv) const;
  
  /**
   * @brief Finds the first of the specified keys in `sv`.
   * @tparam BIter type of iterator to the keys
   * @tparam EIter type of key end-iterator
   * @param sv string to be parsed
   * @param beginKey iterator to the first key
   * @param endKey iterator past the last key
   * @return a view of the key found within `sv`, empty if none
   * 
   * The `keys` are required to be sorted, longest first, since they are tested
   * in order and the first match is kept (e.g. if the first key is `=` and the
   * second is `==`, the second key is never matched since the first one matches
   * first).
   * The first character of the key must not be escaped. Escaped characters in
   * the key are not supported.
   * 
   * If no `key` is found, the returned view is zero-length and pointing to the
   * end of `sv`.
   * 
   * The quoting in `sv` is ignored.
   */
  template <typename BIter, typename EIter>
  std::string_view findFirstUnescaped
    (std::string_view sv, BIter beginKey, EIter endKey) const;
  
  // @{
  /**
   * @brief Finds the first of the specified keys in `sv`.
   * @tparam BIter type of iterator to the keys
   * @tparam EIter type of key end-iterator
   * @param sv string to be parsed
   * @param beginKey iterator to the first key
   * @param endKey iterator past the last key
   * @return a view of the key found within `sv`, empty if none
   * 
   * The `keys` are required to be sorted, longest first, since they are tested
   * in order and the first match is kept (e.g. if the first key is `=` and the
   * second is `==`, the second key is never matched since the first one matches
   * first).
   * The first character of the key must not be escaped. Escaped characters in
   * the key are not supported.
   * 
   * If no `key` is found, the returned view is zero-length and pointing to the
   * end of `sv`.
   * 
   * The quoting in `sv` is ignored.
   */
  template <typename Keys>
  std::string_view findFirstUnescaped
    (std::string_view sv, Keys const& keys) const;
  
  template <typename Key>
  std::string_view findFirstUnescaped
    (std::string_view sv, std::initializer_list<Key> keys) const;
  
  //@}
  
  
  /**
   * @brief Finds the first of the specified keys in the unquoted part of `sv`.
   * @tparam BIter type of iterator to the keys
   * @tparam EIter type of key end-iterator
   * @param sv string to be parsed
   * @param beginKey iterator to the first key
   * @param endKey iterator past the last key
   * @return the view pointing to the key in `sv`, or empty to its end if none
   * 
   * The `keys` are required to be sorted, longest first, since they are tested
   * in order and the first match is kept (e.g. if the first key is `=` and the
   * second is `==`, the second key is never matched since the first one matches
   * first).
   * 
   * If no `key` is found, the returned view is zero-length and pointing to the
   * end of `sv`.
   */
  template <typename BIter, typename EIter>
  std::string_view findFirstUnquoted
    (std::string_view sv, BIter beginKey, EIter endKey) const;
  
  // @{
  /**
   * @brief Finds the first of the specified keys in the unquoted part of `sv`.
   * @tparam BIter type of iterator to the keys
   * @tparam EIter type of key end-iterator
   * @param sv string to be parsed
   * @param beginKey iterator to the first key
   * @param endKey iterator past the last key
   * @return the view pointing to the key in `sv`, or empty to its end if none
   * 
   * The `keys` are required to be sorted, longest first, since they are tested
   * in order and the first match is kept (e.g. if the first key is `=` and the
   * second is `==`, the second key is never matched since the first one matches
   * first).
   * 
   * If no `key` is found, the returned view is zero-length and pointing to the
   * end of `sv`.
   */
  template <typename Keys>
  std::string_view findFirstUnquoted
    (std::string_view sv, Keys const& keys) const;
  
  template <typename Key>
  std::string_view findFirstUnquoted
    (std::string_view sv, std::initializer_list<Key> keys) const;
  
  //@}
  
  
  /**
   * @brief Splits the view `sv` in three: before `sep`, `sep` and after `sep`.
   * @param sv view of the string to split
   * @param sep a subview of `sv` to split at
   * @return a `SplitView_t` object with the three parts split, empty if needed
   * 
   * The view `sep` is required to be a subview of `sv`: it's not enough for it
   * to have as content a substring of `sv`. For example, `splitOn("a:1", ":")`
   * will not work, because the string `"a:1"` does not share data in memory
   * with `":"`.
   * 
   * Even if `sep` is empty, it's still required to point with both `begin()`
   * and `end()` within `sv`, and `sv` will be split according to that point.
   */
  static SplitView_t splitOn(std::string_view sv, std::string_view sep);
  
  /// @}
  // --- END ----- Tokenization ------------------------------------------------
  
  
  // --- BEGIN --- Characters --------------------------------------------------
  /// @name Characters
  /// @{
  
  /// Returns whether `ch` is an escape character.
  bool isEscape(char ch) const { return ch == fParams.escape; }
  
  /**
   * @brief Returns whether the character pointed by `itCh` is escaped or not.
   * @tparam BIter iterator type
   * @param begin iterator to the beginning of the string
   * @param itCh iterator to the character to be investigated.
   * @return whether there is an unescaped escape character before `itCh`
   * 
   * Note that `itCh` may be a end iterator (for an empty string, the result
   * is `false`).
   */
  template <typename BIter>
  bool isCharacterEscaped(BIter begin, BIter itCh) const;
  
  /**
   * @brief Finds the next character satisfying the specified criterion.
   * @tparam Sel type of functor determining which character to consider blank
   * @param s view of the string to be parsed
   * @param select functor determining which character(s) to look for
   * @return an iterator to the first character, `s.end()` if none
   * 
   * By default, the selected character is a blank character `ch`, which has
   * `std::isblank(ch)` `true`.
   */
  template <typename Sel>
  std::string_view::const_iterator findNextCharacter
    (std::string_view s, Sel select) const;
  
  /// Helper function for `findNextCharacter(std::string_view, Sel)`.
  std::string_view::const_iterator findNextBlank(std::string_view s) const
    { return findNextCharacter(s, isBlank); }
  
  /**
   * @brief Consumes the blank characters a the beginning of `s`.
   * @tparam CType type of functor determining which type of character to remove
   * @param s view of the string to be parsed
   * @param charType functor determining which characters to remove
   * @return a view of `s` starting after its trailing `charType` characters
   * @see `removeTrailingBlanks()`
   */
  template <typename CType>
  std::string_view removeTrailingCharacters
    (std::string_view s, CType charType) const;
  
  /// @brief Consumes the blank characters a the beginning of `s`.
  /// @see `removeTrailingCharacters()`
  std::string_view removeTrailingBlanks(std::string_view s) const
    { return removeTrailingCharacters(s, isBlank); }
  
  
  /**
   * @brief Returns a copy of `w` with all escape characters removed.
   * @param w the string to change
   * @return a copy of `w` without escaping
   * @see `removeEscapes(Word const&)`
   * 
   * The escaping scheme that is applied is just to remove the escape
   * character (no replacement table supported here).
   * An unescaped escape character at the end of the string will not be removed.
   * 
   * It is recommended that this be done as the last step of the parsing, since
   * it changes the meaning of the parsing elements like quotations, comments
   * etc.
   * 
   * Note that applying `removeEscapes()` more than once will keep removing
   * characters that in the earlier passes were not considered escapes (for
   * example, four escape characters become two in the first pass, one in the
   * second and disappear in the following passes).
   * 
   */
  std::string removeWordEscapes(std::string&& w) const;
  std::string removeWordEscapes(std::string_view w) const
    { return removeWordEscapes(std::string{ w }); }
  std::string removeWordEscapes(const char* w) const
    { return removeWordEscapes(std::string{ w }); }
  // @}
  
  /**
   * @brief Returns a copy of `words` with all escape characters removed.
   * @tparam Words type of list of words
   * @param words the list of words to change
   * @return the list of words without escaping
   * @see `removeEscapes(std::string)`
   * 
   * The escaping is removed from each of the `words` in the list, which are
   * treated as independent.
   * See `removeEscapes(std::string)` for the details.
   */
  template <typename Words>
  std::vector<std::string> removeEscapes(Words const& words) const;
  
  
  //@{
  /**
   * @brief Returns a copy of `w` with no quotation starts and ends.
   * @param w the string to change
   * @return the word without quotations
   * @see `removeQuotations(Words const&)`
   * 
   * Escaping is still honored (if present).
   * 
   * Note that applying `removeQuotations` more than once will keep removing
   * quotation markings that in the earlier passes were not considered such (for
   * example, `a1 << "b1 << 'c1 << " or " << c2' << b2" << a2` will become first
   * `a1 << b1 << 'c1 <<  or  << c2' << b2 << a2`, and eventually
   * `a1 << b1 << c1 <<  or  << c2 << b2 << a2`).
   * 
   */
  std::string removeWordQuotations(std::string&& w) const;
  std::string removeWordQuotations(std::string_view w) const
    { return removeWordQuotations(std::string{ w }); }
  std::string removeWordQuotations(const char* w) const
    { return removeWordQuotations(std::string{ w }); }
  // @}
  
  /**
   * @brief Returns a copy of `words` with no quotation starts and ends.
   * @tparam Words type of list of words
   * @param words the list of words to change
   * @return the list of words without quotations
   * @see `removeQuotations(std::string)`
   * 
   * The substitution is applied on each of the `words` in the list, which are
   * treated as independent.
   * See `removeQuotations(std::string)` for the details.
   */
  template <typename Words>
  std::vector<std::string> removeQuotations(Words const& words) const;
  
  /// @}
  // --- END ----- Characters --------------------------------------------------
  
  
  /// Creates a `std::string_view` from an entire string `s`.
  static std::string_view make_view(std::string const& s)
    { return make_view(s.begin(), s.end()); }
  
  /// Creates a `std::string_view` from two string iterators `b` and `e`.
  template <typename BIter, typename EIter>
  static std::string_view make_view(BIter b, EIter e)
    { return { &*b, static_cast<std::size_t>(std::distance(b, e)) }; }
  
    private:
  Params_t fParams; ///< Parsing parameters.
  
  // --- BEGIN -- Cache --------------------------------------------------------
  
  /// Start characters of all supported quotations.
  std::string fQuoteStarts;
  
  // --- END ---- Cache --------------------------------------------------------
  
  /// Initializes the parameters and caches.
  void adoptParams(Params_t params);
  
}; // icarus::ParsingToolkit


// -----------------------------------------------------------------------------
struct icarus::ParsingToolkit::Error: std::runtime_error {
  Error(std::string msg): std::runtime_error{ std::move(msg) } {}
};


// -----------------------------------------------------------------------------
// ---  template implementation
// -----------------------------------------------------------------------------
template <typename Delim>
std::vector<std::string_view> icarus::ParsingToolkit::splitWords
  (std::string const& s, Delim isDelimiter /* = isBlank */) const
{
  // REQUIREMENT: escape character must not be classified as delimiter
  assert(!isDelimiter(fParams.escape));
  // REQUIREMENT: the first character of no quotation start must be classified
  //              as delimiter
  assert(
    std::count_if(fQuoteStarts.cbegin(), fQuoteStarts.cend(), isDelimiter) == 0
    );
  
  
  // helper class:
  // stores the word as collected so far, updates `sv` and starts new words
  class WordTracker {
    ParsingToolkit const& tk;
    Delim const& isDelimiter;
    std::string_view& sv;
    std::vector<std::string_view> words;
    std::string_view::const_iterator wStart;
      public:
    WordTracker(ParsingToolkit const& tk, Delim const& d, std::string_view& sv)
      : tk{ tk }, isDelimiter{ d }, sv{ consumeDelim(sv) }, wStart{ sv.begin() }
      {}
    void startNew()
      {
        words.push_back(make_view(wStart, sv.begin()));
        wStart = consumeDelim().begin();
      }
    void moveEndTo(std::string_view::const_iterator it)
      { moveEndBy(it - sv.begin()); }
    void moveEndBy(std::size_t n) { sv.remove_prefix(n); }
    std::vector<std::string_view> finish()
      { if (wStart != sv.begin()) startNew(); return std::move(words); }
    std::string_view& consumeDelim(std::string_view& s) const
      { return s = tk.removeTrailingCharacters(s, isDelimiter); }
    std::string_view& consumeDelim() { return consumeDelim(sv); }
  }; // WordTracker
  
  std::string_view sv = make_view(s);
  WordTracker words { *this, isDelimiter, sv }; // shares sv management
  
  // sv.begin() is kept updated to the candidate end of word;
  // the beginning of the current word is always cached as words.wStart
  while (!sv.empty()) {

    // process up to the next quotation
    auto const [ qsv, qptr ] = findQuotationStart(sv);
    
    // parse and split until the quotation start:
    auto const qstart = qsv.begin();
    while(true) {
      
      // find next space;
      // if next space is past the quotation, stop to the quotation instead
      words.moveEndTo
        (findNextCharacter(make_view(sv.begin(), qstart), isDelimiter));
      
      if (sv.begin() == qstart) break;
      
      // not the quote? it's a delimiter! new word found:
      words.startNew();
      
    } // while(true)
    
    // handle the quoted part
    if (qptr) {
      assert(sv.substr(0, qptr->first.length()) == qptr->first);
      
      words.moveEndBy(qptr->first.length());
      
      // find the end of the quote, and swallow it into the current word
      std::string_view const quotEnd = findQuotationEnd(sv, qptr->second);
      words.moveEndTo(quotEnd.begin());
      
      // if we have found a end of quote, swallow it too (otherwise it's over)
      if (!quotEnd.empty()) words.moveEndBy(qptr->second.length());
      
    } // if quotation found
    
  } // while
  
  return words.finish();

} // icarus::ParsingToolkit::splitWords()


// -----------------------------------------------------------------------------
template <typename BIter, typename EIter>
std::string_view icarus::ParsingToolkit::findFirstUnescaped
  (std::string_view sv, BIter beginKey, EIter endKey) const
{
  
  typename std::iterator_traits<BIter>::value_type const* key = nullptr;
  std::size_t keyPos = std::string_view::npos;
  
  for (auto iKey = beginKey; iKey != endKey; ++iKey) {
    // find where this key is (unescaped)
    std::size_t pos = 0;
    while (pos < sv.length()) {
      pos = sv.find(*iKey, pos);
      if (!isCharacterEscaped(sv.begin(), sv.begin() + pos)) break;
      ++pos;
    }
    // is this the first among the keys?
    if (pos >= std::min(keyPos, sv.length())) continue;
    key = &*iKey;
    keyPos = pos;
  } // for keys
  
  // return a substring of sv, not key
  if (key) {
    using std::begin, std::end;
    std::size_t const keyLength = make_view(*key).length();
    return { sv.data() + keyPos, keyLength };
  }
  else return { sv.data() + sv.length(), 0 };
} // icarus::ParsingToolkit::findFirstUnescaped()


// -----------------------------------------------------------------------------
template <typename Keys>
std::string_view icarus::ParsingToolkit::findFirstUnescaped
  (std::string_view sv, Keys const& keys) const
{
  using std::begin, std::end;
  return findFirstUnescaped(sv, begin(keys), end(keys)); 
} // icarus::ParsingToolkit::findFirstUnescaped(Keys)


// -----------------------------------------------------------------------------
template <typename Key>
std::string_view icarus::ParsingToolkit::findFirstUnescaped
  (std::string_view sv, std::initializer_list<Key> keys) const
  { return findFirstUnescaped(sv, keys.begin(), keys.end()); }


// -----------------------------------------------------------------------------
template <typename BIter, typename EIter>
std::string_view icarus::ParsingToolkit::findFirstUnquoted
  (std::string_view sv, BIter beginKey, EIter endKey) const
{
  
  // if a key is found between `b` and `e`, returns `sv` split around the key;
  // otherwise, all `sv` is in post
  auto findKey = [this,beginKey,endKey]
    (std::string_view::const_iterator b, std::string_view::const_iterator e)
    { return findFirstUnescaped(make_view(b, e), beginKey, endKey); };
  
  std::string_view key{ sv.data() + sv.length(), 0 };
  while (!sv.empty()) {
    
    // find the next quotation
    auto const [ fromQ, qptr ] = findQuotationStart(sv);
    
    // search in the unquoted part
    key = findKey(sv.begin(), fromQ.begin());
    if (!key.empty()) break;
    
    // skip the quotation; if there is no quotation, we are done
    if (!qptr) break;
    
    sv = fromQ;
    sv.remove_prefix(qptr->first.length()); // skip the quotation start
    
    // find the end of quotation
    std::string_view const afterQ = findQuotationEnd(sv, qptr->second);
    
    if (afterQ.empty()) { // begin of quotation, but no end: no good
      // so we don't consider this as quotation: search in the "quoted" part
      key = findKey(fromQ.begin(), fromQ.end());
      break;
    } // if

    // skip the quoted material, and the quotation end too
    sv = afterQ;
    sv.remove_prefix(qptr->second.length());
    
  } // while
  
  return key;
  
} // icarus::ParsingToolkit::findFirstUnquoted(Iter)


// -----------------------------------------------------------------------------
template <typename Keys>
std::string_view icarus::ParsingToolkit::findFirstUnquoted
  (std::string_view sv, Keys const& keys) const
{
  using std::begin, std::end;
  return findFirstUnquoted(sv, begin(keys), end(keys)); 
} // icarus::ParsingToolkit::findFirstUnquoted(Keys)


// -----------------------------------------------------------------------------
template <typename Key>
std::string_view icarus::ParsingToolkit::findFirstUnquoted
  (std::string_view sv, std::initializer_list<Key> keys) const
  { return findFirstUnquoted(sv, keys.begin(), keys.end()); }


// -----------------------------------------------------------------------------
template <typename Iter>
Iter icarus::ParsingToolkit::findCommentWord(Iter beginWord, Iter endWord) const
{
  for (auto it = beginWord; it != endWord; ++it) {
    if (std::equal(fParams.comment.begin(), fParams.comment.end(), begin(*it)))
      return it;
  } // for
  return endWord;
} // icarus::ParsingToolkit::findCommentWord()


// -----------------------------------------------------------------------------
template <typename Iter>
bool icarus::ParsingToolkit::isCharacterEscaped(Iter begin, Iter itCh) const
{
  unsigned int nEscapes = 0U;
  while (itCh-- != begin) {
    
    if (!isEscape(*itCh)) break;
    ++nEscapes;
    
  } // while
  
  return (nEscapes & 1) == 1; // odd number of escapes means escaped
  
} // icarus::ParsingToolkit::isCharacterEscaped()


// -----------------------------------------------------------------------------
template <typename Sel>
std::string_view::const_iterator icarus::ParsingToolkit::findNextCharacter
  (std::string_view s, Sel selector) const
{
  auto const sbegin = s.begin(), send = s.end();
  auto it = sbegin;
  while (it != send) {
    it = std::find_if(it, send, selector);
    if (!isCharacterEscaped(sbegin, it)) return it;
    ++it; // skip the escaped character and move on
  } // while
  return send;
} // icarus::ParsingToolkit::findNextCharacter()


// -----------------------------------------------------------------------------
template <typename CType>
std::string_view icarus::ParsingToolkit::removeTrailingCharacters
  (std::string_view s, CType charType) const
{
  // REQUIREMENT: escape character must not be classified as delimiter
  assert(!charType(fParams.escape));
  
  while (!s.empty()) {
    if (!charType(s.front())) break; // escape character triggers this too
    s.remove_prefix(1U);
  } // while
  return s;
} // icarus::ParsingToolkit::removeTrailingCharacters()


// -----------------------------------------------------------------------------
template <typename Words>
std::vector<std::string> icarus::ParsingToolkit::removeEscapes
  (Words const& words) const
{
  using std::size;
  std::vector<std::string> nv;
  nv.reserve(size(words));
  for (auto const& word: words) nv.push_back(removeWordEscapes(word));
  return nv;
} // icarus::ParsingToolkit::removeEscapes()


// -----------------------------------------------------------------------------
template <typename Words>
std::vector<std::string> icarus::ParsingToolkit::removeQuotations
  (Words const& words) const
{
  using std::size;
  std::vector<std::string> nv;
  nv.reserve(size(words));
  for (auto const& word: words) nv.push_back(removeWordQuotations(word));
  return nv;
} // icarus::ParsingToolkit::removeEscapes()


// -----------------------------------------------------------------------------

#endif // ICARUSCODE_PMT_ALGORITHMS_PARSINGTOOLKIT_H
