/**
 * @file ParsingToolkit_test.cc
 * @brief Unit test for utilities in `ParsingToolkit.h`
 * @date May 13, 2022
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @see icaruscode/PMT/Algorithms/ParsingToolkit.h
 * 
 */

// ICARUS libraries
#include "icaruscode/PMT/Algorithms/ParsingToolkit.h"

// Boost libraries
#define BOOST_TEST_MODULE ( ParsingToolkit_test )
#include <boost/test/unit_test.hpp>

// C/C++ standard libraries
#include <tuple> // std::tie()
#include <vector>
#include <string_view>
#include <string>
#include <sstream>
#include <cassert>


// -----------------------------------------------------------------------------
// --- implementation detail tests
// -----------------------------------------------------------------------------
void isCharacterEscaped_test() {
  
  icarus::ParsingToolkit const tk; // default configuration
  
  std::string const s { R"(\a\\a\\\a\\\\a\\\\\a\\\\\\a\\\)" };
  auto const b = s.begin();
  auto it = b;
  
  BOOST_TEST(!tk.isCharacterEscaped(b, it++)); // \                              !
  BOOST_TEST( tk.isCharacterEscaped(b, it++)); // \a                             !
  BOOST_TEST(!tk.isCharacterEscaped(b, it++)); // \a\                            !
  BOOST_TEST( tk.isCharacterEscaped(b, it++)); // \a\\                           !
  BOOST_TEST(!tk.isCharacterEscaped(b, it++)); // \a\\a                          !
  BOOST_TEST(!tk.isCharacterEscaped(b, it++)); // \a\\a\                         !
  BOOST_TEST( tk.isCharacterEscaped(b, it++)); // \a\\a\\                        !
  BOOST_TEST(!tk.isCharacterEscaped(b, it++)); // \a\\a\\\                       !
  BOOST_TEST( tk.isCharacterEscaped(b, it++)); // \a\\a\\\a                      !
  BOOST_TEST(!tk.isCharacterEscaped(b, it++)); // \a\\a\\\a\                     !
  BOOST_TEST( tk.isCharacterEscaped(b, it++)); // \a\\a\\\a\\                    !
  BOOST_TEST(!tk.isCharacterEscaped(b, it++)); // \a\\a\\\a\\\                   !
  BOOST_TEST( tk.isCharacterEscaped(b, it++)); // \a\\a\\\a\\\\                  !
  BOOST_TEST(!tk.isCharacterEscaped(b, it++)); // \a\\a\\\a\\\\a                 !
  BOOST_TEST(!tk.isCharacterEscaped(b, it++)); // \a\\a\\\a\\\\a\                !
  BOOST_TEST( tk.isCharacterEscaped(b, it++)); // \a\\a\\\a\\\\a\\               !
  BOOST_TEST(!tk.isCharacterEscaped(b, it++)); // \a\\a\\\a\\\\a\\\              !
  BOOST_TEST( tk.isCharacterEscaped(b, it++)); // \a\\a\\\a\\\\a\\\\             !
  BOOST_TEST(!tk.isCharacterEscaped(b, it++)); // \a\\a\\\a\\\\a\\\\\            !
  BOOST_TEST( tk.isCharacterEscaped(b, it++)); // \a\\a\\\a\\\\a\\\\\a           !
  BOOST_TEST(!tk.isCharacterEscaped(b, it++)); // \a\\a\\\a\\\\a\\\\\a\          !
  BOOST_TEST( tk.isCharacterEscaped(b, it++)); // \a\\a\\\a\\\\a\\\\\a\\         !
  BOOST_TEST(!tk.isCharacterEscaped(b, it++)); // \a\\a\\\a\\\\a\\\\\a\\\        !
  BOOST_TEST( tk.isCharacterEscaped(b, it++)); // \a\\a\\\a\\\\a\\\\\a\\\\       !
  BOOST_TEST(!tk.isCharacterEscaped(b, it++)); // \a\\a\\\a\\\\a\\\\\a\\\\\      !
  BOOST_TEST( tk.isCharacterEscaped(b, it++)); // \a\\a\\\a\\\\a\\\\\a\\\\\\     !
  BOOST_TEST(!tk.isCharacterEscaped(b, it++)); // \a\\a\\\a\\\\a\\\\\a\\\\\\a    !
  BOOST_TEST(!tk.isCharacterEscaped(b, it++)); // \a\\a\\\a\\\\a\\\\\a\\\\\\a\   !
  BOOST_TEST( tk.isCharacterEscaped(b, it++)); // \a\\a\\\a\\\\a\\\\\a\\\\\\a\\  !
  BOOST_TEST(!tk.isCharacterEscaped(b, it++)); // \a\\a\\\a\\\\a\\\\\a\\\\\\a\\\ !
  assert(it == s.end());
  BOOST_TEST( tk.isCharacterEscaped(b, it));
  
} // isCharacterEscaped_test()


// -----------------------------------------------------------------------------
void removeTrailingBlanks_test() {
  
  using namespace std::string_literals;
  using namespace std::string_view_literals;
  icarus::ParsingToolkit tk; // default configuration
  
  BOOST_TEST(tk.removeTrailingBlanks(""s) == ""sv);
  
  BOOST_TEST(tk.removeTrailingBlanks("a b c "s) == "a b c "sv);
  
  BOOST_TEST(tk.removeTrailingBlanks(" \\  a b c "s) == "\\  a b c "sv);
  
  BOOST_TEST(tk.removeTrailingBlanks(" \t  a b c "s) == "a b c "sv);
  
} // removeTrailingBlanks_test()


// -----------------------------------------------------------------------------
void findNextBlank_test() {
  
  using namespace std::string_literals;
  using namespace std::string_view_literals;
  icarus::ParsingToolkit tk; // default configuration
  
  std::string_view s;
  BOOST_TEST(std::distance(s.cbegin(), tk.findNextBlank(s)) == 0);
  
  s = "a b c "sv;
  BOOST_TEST(std::distance(s.cbegin(), tk.findNextBlank(s)) == 1);
  
  s = "ab\\ cd e "sv;
  BOOST_TEST(std::distance(s.cbegin(), tk.findNextBlank(s)) == 6);
  
  s = "a   "sv;
  BOOST_TEST(std::distance(s.cbegin(), tk.findNextBlank(s)) == 1);
  
} // findNextBlank_test()


// -----------------------------------------------------------------------------
void readMultiline_test() {
  
  using namespace std::string_literals;
  
  std::vector const lines {
      "This line is not empty. The previous and the next two are."s
    , "This is a multiline line of one, "s
    , "two, "s
    , "three lines."s
    , R"(This "is a" "multi\"quoted\" strange)"s
    , R"( line")"s
    , "This is a normal line, but its not terminated end gets to: EOF"s
    };
  
  std::istringstream stream {
      "\n"                 // <--
    + lines[0] + "\n"      // <--
    + "\n"                 // <--
    + "\n"                 // <--
    + lines[1] + "\\\n"    // -.
    + lines[2] + "\\\n"    //  |<--
    + lines[3] + "\n"      // -'
    + lines[4] + "\n"      // -.___
    + lines[5] + "\n"      // -'
    + lines[6]             // <--
    };
  
  icarus::ParsingToolkit tk; // default configuration
  
  assert(stream);
  
  auto [ line, nPieces ] = tk.readMultiline(stream);
  BOOST_TEST(line == ""s);
  BOOST_TEST(nPieces == 1U);
  BOOST_TEST(bool(stream));
  
  std::tie(line, nPieces) = tk.readMultiline(stream);
  BOOST_TEST(line == lines[0]);
  BOOST_TEST(nPieces == 1U);
  BOOST_TEST(bool(stream));
  
  std::tie(line, nPieces) = tk.readMultiline(stream);
  BOOST_TEST(line == "");
  BOOST_TEST(nPieces == 1U);
  BOOST_TEST(bool(stream));
  
  std::tie(line, nPieces) = tk.readMultiline(stream);
  BOOST_TEST(line == "");
  BOOST_TEST(nPieces == 1U);
  BOOST_TEST(bool(stream));
  
  std::tie(line, nPieces) = tk.readMultiline(stream);
  BOOST_TEST(line == lines[1] + lines[2] + lines[3]);
  BOOST_TEST(nPieces == 3U);
  BOOST_TEST(bool(stream));
  
  std::tie(line, nPieces) = tk.readMultiline(stream);
  BOOST_TEST(line == lines[4] + "\n"s + lines[5]);
  BOOST_TEST(nPieces == 2U);
  BOOST_TEST(bool(stream));
  
  std::tie(line, nPieces) = tk.readMultiline(stream);
  BOOST_TEST(line == lines[6]);
  BOOST_TEST(nPieces == 1U);
  BOOST_TEST(bool(stream));
  
  std::tie(line, nPieces) = tk.readMultiline(stream);
  BOOST_TEST(line.empty());
  BOOST_TEST(nPieces == 0U);
  BOOST_TEST(!stream);
  
} // readMultiline_test()


void readMultiline_endOfInput_test() {
  
  using namespace std::string_literals;
  
  icarus::ParsingToolkit tk; // default configuration
  
  std::istringstream stream;
  auto const refillStream = [&stream](std::string s) -> std::istringstream&
    { stream.clear(); stream.str(std::move(s)); return stream; };
  
  refillStream(
    R"()"
    );
  auto [ line, nPieces ] = tk.readMultiline(stream);
  BOOST_TEST(line == "");
  BOOST_TEST(nPieces == 0U);
  
  refillStream(
    R"(The end.)" "\n"
    );
  std::tie(line, nPieces) = tk.readMultiline(stream);
  BOOST_TEST(line == "The end.");
  BOOST_TEST(nPieces == 1U);
  std::tie(line, nPieces) = tk.readMultiline(stream);
  BOOST_TEST(line == "");
  BOOST_TEST(nPieces == 0U);
  
  refillStream(
    R"(Unfinished business)"
    );
  std::tie(line, nPieces) = tk.readMultiline(stream);
  BOOST_TEST(line == "Unfinished business");
  BOOST_TEST(nPieces == 1U);
  
  refillStream(
    R"(Unfinished business\)"
    );
  std::tie(line, nPieces) = tk.readMultiline(stream);
  BOOST_TEST(line == "Unfinished business");
  BOOST_TEST(nPieces == 1U);
  
  refillStream(
    R"(Unfinished "business)"
    );
  std::tie(line, nPieces) = tk.readMultiline(stream);
  BOOST_TEST(line == R"(Unfinished "business)");
  BOOST_TEST(nPieces == 1U);
  
  refillStream(
    R"(Unfinished\)"
    "\n"
    R"( business)"
    );
  std::tie(line, nPieces) = tk.readMultiline(stream);
  BOOST_TEST(line == "Unfinished business");
  BOOST_TEST(nPieces == 2U);
  
  refillStream(
    R"("Unfinished)"
    "\n"
    R"( business")"
    );
  std::tie(line, nPieces) = tk.readMultiline(stream);
  BOOST_TEST(line == "\"Unfinished\n business\""s);
  BOOST_TEST(nPieces == 2U);
  
  refillStream(
    R"("Unfinished)"
    "\n"
    R"( business)"
    );
  std::tie(line, nPieces) = tk.readMultiline(stream);
  BOOST_TEST(line == "\"Unfinished\n business");
  BOOST_TEST(nPieces == 2U);
  
} // readMultiline_endOfInput_test()


// -----------------------------------------------------------------------------
void splitWords_test() {
  
  using namespace std::string_literals;
  
  using Words_t = std::vector<std::string>;
  
  icarus::ParsingToolkit tk; // default configuration
  
  {
    std::string const s { R"(a b c)" };
    auto const& words = tk.splitWords(s);
    Words_t const expected { R"(a)"s, R"(b)"s, R"(c)"s };
    BOOST_CHECK_EQUAL_COLLECTIONS
      (words.begin(), words.end(), expected.begin(), expected.end());
  }
  
  {
    std::string const s { R"( a "bb"  c )" };
    auto const& words = tk.splitWords(s);
    Words_t const expected { R"(a)"s, R"("bb")"s, R"(c)"s };
    BOOST_CHECK_EQUAL_COLLECTIONS
      (words.begin(), words.end(), expected.begin(), expected.end());
  }
  
  {
    std::string const s { R"(a "b c")" };
    auto const& words = tk.splitWords(s);
    Words_t const expected { R"(a)"s, R"("b c")"s };
    BOOST_CHECK_EQUAL_COLLECTIONS
      (words.begin(), words.end(), expected.begin(), expected.end());
  }
  
  {
    std::string const s { R"("")" };
    auto const& words = tk.splitWords(s);
    Words_t const expected { R"("")"s };
    BOOST_CHECK_EQUAL_COLLECTIONS
      (words.begin(), words.end(), expected.begin(), expected.end());
  }
  
  {
    std::string const s { R"(a\ b c)" };
    auto const& words = tk.splitWords(s);
    Words_t const expected { R"(a\ b)"s, R"(c)"s };
    BOOST_CHECK_EQUAL_COLLECTIONS
      (words.begin(), words.end(), expected.begin(), expected.end());
  }
  
  {
    std::string const s { R"(a b" c")" };
    auto const& words = tk.splitWords(s);
    Words_t const expected { R"(a)"s, R"(b" c")"s };
    BOOST_CHECK_EQUAL_COLLECTIONS
      (words.begin(), words.end(), expected.begin(), expected.end());
  }
  
  {
    std::string const s { R"("a "b \"c d"")" };
    auto const& words = tk.splitWords(s);
    Words_t const expected { R"("a "b)"s, R"(\"c)"s, R"(d"")"s };
    BOOST_CHECK_EQUAL_COLLECTIONS
      (words.begin(), words.end(), expected.begin(), expected.end());
  }
  
  
} // splitWords_test()


// -----------------------------------------------------------------------------
void findFirstUnescaped_test() {
  
  using namespace std::string_view_literals;
  
  icarus::ParsingToolkit const tk; // default configuration
  
  {
    std::string_view const sv { R"()" };
    std::string_view const key = tk.findFirstUnescaped(sv, { ":+", ":" });
    BOOST_TEST(key == "");
    BOOST_TEST(std::distance(sv.begin(), key.begin()) == 0);
  }
  
  {
    std::string_view const sv { R"(a :: b)" };
    std::string_view const key = tk.findFirstUnescaped(sv, { ":+", ":" });
    BOOST_TEST(key == ":");
    BOOST_TEST(std::distance(sv.begin(), key.begin()) == 2);
  }
  
  {
    std::string_view const sv { R"(a ::+ b)" };
    std::string_view const key = tk.findFirstUnescaped(sv, { ":+", ":" });
    BOOST_TEST(key == ":");
    BOOST_TEST(std::distance(sv.begin(), key.begin()) == 2);
  }
  
  {
    std::string_view const sv { R"(a \::+ b)" };
    std::string_view const key = tk.findFirstUnescaped(sv, { ":+", ":" });
    BOOST_TEST(key == ":+");
    BOOST_TEST(std::distance(sv.begin(), key.begin()) == 4);
  }
  
  {
    std::string_view const sv { R"(a\ ::+ b)" };
    std::string_view const key = tk.findFirstUnescaped(sv, { ":+", ":" });
    BOOST_TEST(key == ":");
    BOOST_TEST(std::distance(sv.begin(), key.begin()) == 3);
  }
  
  {
    std::string_view const sv { R"(a :+)" };
    std::string_view const key = tk.findFirstUnescaped(sv, { ":+", ":" });
    BOOST_TEST(key == ":+");
    BOOST_TEST(std::distance(sv.begin(), key.begin()) == 2);
  }
  
  {
    std::string_view const sv { R"(a \:+)" };
    std::string_view const key = tk.findFirstUnescaped(sv, { ":+", ":" });
    BOOST_TEST(key == "");
    BOOST_TEST(std::distance(sv.begin(), key.begin()) == 5);
  }
  
  {
    std::string_view const sv { R"(:+ b)" };
    std::string_view const key = tk.findFirstUnescaped(sv, { ":+", ":" });
    BOOST_TEST(key == ":+");
    BOOST_TEST(std::distance(sv.begin(), key.begin()) == 0);
  }
  
  {
    std::string_view const sv { R"(nope)" };
    std::string_view const key = tk.findFirstUnescaped(sv, { ":+", ":" });
    BOOST_TEST(key == "");
    BOOST_TEST(std::distance(sv.begin(), key.begin()) == 4);
  }
  
  
} // findFirstUnescaped_test()


// -----------------------------------------------------------------------------
void findFirstUnquoted_test() {
  
  using namespace std::string_view_literals;
  
  icarus::ParsingToolkit const tk; // default configuration
  
  auto const findAndSplit = [tk](std::string_view sv)
    { return tk.splitOn(sv, tk.findFirstUnquoted(sv, { ":+", ":" })); };
  
  {
    auto const [ pre, sep, post ] = findAndSplit(""sv);
    BOOST_TEST(pre  == R"()"sv);
    BOOST_TEST(sep  == R"()"sv);
    BOOST_TEST(post == R"()"sv);
  }
  
  {
    auto const [ pre, sep, post ] = findAndSplit(R"(a:+b)"sv);
    BOOST_TEST(pre  == R"(a)"sv);
    BOOST_TEST(sep  == R"(:+)"sv);
    BOOST_TEST(post == R"(b)"sv);
  }
  
  {
    auto const [ pre, sep, post ] = findAndSplit(R"(a::b)"sv);
    BOOST_TEST(pre  == R"(a)"sv);
    BOOST_TEST(sep  == R"(:)"sv);
    BOOST_TEST(post == R"(:b)"sv);
  }
  
  {
    auto const [ pre, sep, post ] = findAndSplit(R"(a\:b)"sv);
    BOOST_TEST(pre  == R"(a\:b)"sv);
    BOOST_TEST(sep  == R"()"sv);
    BOOST_TEST(post == R"()"sv);
  }
  
  {
    auto const [ pre, sep, post ] = findAndSplit(R"(a :+ b)"sv);
    BOOST_TEST(pre  == R"(a )"sv);
    BOOST_TEST(sep  == R"(:+)"sv);
    BOOST_TEST(post == R"( b)"sv);
  }
  
  {
    auto const [ pre, sep, post ] = findAndSplit(R"(a\ :+ b)"sv);
    BOOST_TEST(pre  == R"(a\ )"sv);
    BOOST_TEST(sep  == R"(:+)"sv);
    BOOST_TEST(post == R"( b)"sv);
  }
  
  {
    auto const [ pre, sep, post ] = findAndSplit(R"(a\ \::+ b)"sv);
    BOOST_TEST(pre  == R"(a\ \:)"sv);
    BOOST_TEST(sep  == R"(:+)"sv);
    BOOST_TEST(post == R"( b)"sv);
  }
  
  {
    auto const [ pre, sep, post ] = findAndSplit(R"(a\ \:+ b)"sv);
    BOOST_TEST(pre  == R"(a\ \:+ b)"sv);
    BOOST_TEST(sep  == R"()"sv);
    BOOST_TEST(post == R"()"sv);
  }
  
  {
    auto const [ pre, sep, post ] = findAndSplit(R"(a\ :+\ b)"sv);
    BOOST_TEST(pre  == R"(a\ )"sv);
    BOOST_TEST(sep  == R"(:+)"sv);
    BOOST_TEST(post == R"(\ b)"sv);
  }
  
  {
    auto const [ pre, sep, post ] = findAndSplit(R"("a":+"b")"sv);
    BOOST_TEST(pre  == R"("a")"sv);
    BOOST_TEST(sep  == R"(:+)"sv);
    BOOST_TEST(post == R"("b")"sv);
  }
  
  {
    auto const [ pre, sep, post ] = findAndSplit(R"("a:":+"b")"sv);
    BOOST_TEST(pre  == R"("a:")"sv);
    BOOST_TEST(sep  == R"(:+)"sv);
    BOOST_TEST(post == R"("b")"sv);
  }
  
  {
    auto const [ pre, sep, post ] = findAndSplit(R"(a:":+"b)"sv);
    BOOST_TEST(pre  == R"(a)"sv);
    BOOST_TEST(sep  == R"(:)"sv);
    BOOST_TEST(post == R"(":+"b)"sv);
  }
  
  {
    auto const [ pre, sep, post ] = findAndSplit(R"("a:"aa:+b"b)"sv);
    BOOST_TEST(pre  == R"("a:"aa)"sv);
    BOOST_TEST(sep  == R"(:+)"sv);
    BOOST_TEST(post == R"(b"b)"sv);
  }
  
} // findFirstUnquoted_test()


// -----------------------------------------------------------------------------
void removeCommentLine_test() {
  
  using namespace std::string_view_literals;
  
  icarus::ParsingToolkit const tk; // default configuration
  
  std::vector<std::string_view> words, expected;
  BOOST_TEST(
    std::distance(words.begin(), tk.findCommentWord(words.begin(), words.end()))
    == 0
    );
  BOOST_CHECK_EQUAL_COLLECTIONS
    (words.begin(), words.end(), expected.begin(), expected.end());
  
  
  words = { " One"sv, "Two#"sv, "Th#ree"sv, " #Four"sv, "Five"sv, "Six"sv };
  expected = words;
  BOOST_TEST(
    std::distance(words.begin(), tk.findCommentWord(words.begin(), words.end()))
    == 6
    );
  tk.removeCommentLine(words);
  BOOST_CHECK_EQUAL_COLLECTIONS
    (words.begin(), words.end(), expected.begin(), expected.end());
  
  
  words = { "#"sv, "a"sv, "long"sv, "comment"sv };
  expected = {};
  BOOST_TEST(
    std::distance(words.begin(), tk.findCommentWord(words.begin(), words.end()))
    == 0
    );
  tk.removeCommentLine(words);
  BOOST_CHECK_EQUAL_COLLECTIONS
    (words.begin(), words.end(), expected.begin(), expected.end());
  
  
  words = { "#acompactcomment"sv, "!"sv };
  expected = {};
  BOOST_TEST(
    std::distance(words.begin(), tk.findCommentWord(words.begin(), words.end()))
    == 0
    );
  tk.removeCommentLine(words);
  BOOST_CHECK_EQUAL_COLLECTIONS
    (words.begin(), words.end(), expected.begin(), expected.end());
  
  
  words = { " One"sv, "Two#"sv, "Th#ree"sv, " #Four"sv, "#Five"sv, "Six"sv };
  expected = { " One"sv, "Two#"sv, "Th#ree"sv, " #Four"sv };
  BOOST_TEST(
    std::distance(words.begin(), tk.findCommentWord(words.begin(), words.end()))
    == 4
    );
  tk.removeCommentLine(words);
  BOOST_CHECK_EQUAL_COLLECTIONS
    (words.begin(), words.end(), expected.begin(), expected.end());
  
  
  words = { " One"sv, "Two#"sv, "Th#ree"sv, "\\#Four"sv, "#"sv, "Six"sv };
  expected = { " One"sv, "Two#"sv, "Th#ree"sv, "\\#Four"sv };
  BOOST_TEST(
    std::distance(words.begin(), tk.findCommentWord(words.begin(), words.end()))
    == 4
    );
  tk.removeCommentLine(words);
  BOOST_CHECK_EQUAL_COLLECTIONS
    (words.begin(), words.end(), expected.begin(), expected.end());
  
  
} // removeCommentLine_test()


// -----------------------------------------------------------------------------
void findQuotationStart_test() {
  
  using namespace std::string_literals;
  
  std::string_view sv {
    R"(No quote. Still no \"quote\". "a" "b c" 'd' "e 'f' g" "h 'i" "j' k" "unfinished 'm')"
    };
  
  icarus::ParsingToolkit tk; // default configuration
  
  std::string_view qsv;
  icarus::ParsingToolkit::QuotSpec_t const* qptr = nullptr;
  
  std::tie(qsv, qptr) = tk.findQuotationStart(sv);
  BOOST_TEST_REQUIRE(qsv == R"("a" "b c" 'd' "e 'f' g" "h 'i" "j' k" "unfinished 'm')");
  BOOST_TEST_REQUIRE(qptr != nullptr);
  BOOST_TEST_REQUIRE(qptr->first == "\""s);
  
  sv = qsv;
  sv.remove_prefix(3);
  std::tie(qsv, qptr) = tk.findQuotationStart(sv);
  BOOST_TEST_REQUIRE(qsv == R"("b c" 'd' "e 'f' g" "h 'i" "j' k" "unfinished 'm')");
  BOOST_TEST_REQUIRE(qptr != nullptr);
  BOOST_TEST_REQUIRE(qptr->first == "\""s);
  
  sv = qsv;
  sv.remove_prefix(5);
  std::tie(qsv, qptr) = tk.findQuotationStart(sv);
  BOOST_TEST_REQUIRE(qsv == R"('d' "e 'f' g" "h 'i" "j' k" "unfinished 'm')");
  BOOST_TEST_REQUIRE(qptr != nullptr);
  BOOST_TEST_REQUIRE(qptr->first == "'"s);
  
  sv = qsv;
  sv.remove_prefix(3);
  std::tie(qsv, qptr) = tk.findQuotationStart(sv);
  BOOST_TEST_REQUIRE(qsv == R"("e 'f' g" "h 'i" "j' k" "unfinished 'm')");
  BOOST_TEST_REQUIRE(qptr != nullptr);
  BOOST_TEST_REQUIRE(qptr->first == "\""s);
  
  sv = qsv;
  sv.remove_prefix(9);
  std::tie(qsv, qptr) = tk.findQuotationStart(sv);
  BOOST_TEST_REQUIRE(qsv == R"("h 'i" "j' k" "unfinished 'm')");
  BOOST_TEST_REQUIRE(qptr != nullptr);
  BOOST_TEST_REQUIRE(qptr->first == "\""s);
  
  sv = qsv;
  sv.remove_prefix(6);
  std::tie(qsv, qptr) = tk.findQuotationStart(sv);
  BOOST_TEST_REQUIRE(qsv == R"("j' k" "unfinished 'm')");
  BOOST_TEST_REQUIRE(qptr != nullptr);
  BOOST_TEST_REQUIRE(qptr->first == "\""s);
  
  sv = qsv;
  sv.remove_prefix(6);
  std::tie(qsv, qptr) = tk.findQuotationStart(sv);
  BOOST_TEST_REQUIRE(qsv == R"("unfinished 'm')");
  BOOST_TEST_REQUIRE(qptr != nullptr);
  BOOST_TEST_REQUIRE(qptr->first == "\""s);
  
  sv = qsv;
  sv.remove_prefix(1);
  std::tie(qsv, qptr) = tk.findQuotationStart(sv);
  BOOST_TEST_REQUIRE(qsv == R"('m')");
  BOOST_TEST_REQUIRE(qptr != nullptr);
  BOOST_TEST_REQUIRE(qptr->first == "'"s);
  
  sv = qsv;
  sv.remove_prefix(3);
  std::tie(qsv, qptr) = tk.findQuotationStart(sv);
  BOOST_TEST_REQUIRE(qsv.empty());
  BOOST_TEST_REQUIRE(qptr == nullptr);
  
} // findQuotationStart_test()


void findQuotationStart_noquote_test() {
  
  icarus::ParsingToolkit tk; // default configuration
  
  std::string_view sv { "No quotes at all." };
  auto const [ qsv, qptr ] = tk.findQuotationStart(sv);
  BOOST_TEST(qsv.empty());
  BOOST_TEST(qptr == nullptr);
  
  // special feature: in case there is no match, we still point to the end
  BOOST_TEST((qsv.begin() == sv.end()));
  
} // findQuotationStart_noquote_test()


// -----------------------------------------------------------------------------
void findQuotationEnd_test() {
  
  using namespace std::string_literals;
  using namespace std::string_view_literals;
  
  icarus::ParsingToolkit tk; // default configuration
  
  BOOST_TEST(tk.findQuotationEnd(R"('already'!)",    "'"s) == "'already'!"sv);
  BOOST_TEST(tk.findQuotationEnd(R"(already'!)",     "'"s) == "'!"sv);
  BOOST_TEST(tk.findQuotationEnd(R"(\'already\''!)", "'"s) == "'!"sv);
  BOOST_TEST(tk.findQuotationEnd(R"(Nope.)",         "'"s).empty());
  
} // findQuotationEnd_test()


void findQuotationEnd_noquote_test() {
  
  icarus::ParsingToolkit tk; // default configuration
  
  std::string_view const sv { "No end quotes at all." };
  
  // the test string should have no quotes at all:
  BOOST_TEST_REQUIRE(tk.findQuotationStart(sv).second == nullptr);
  
  std::string const& firstEndQuote = tk.params().quotes.front().second; // any
  
  std::string_view const qsv = tk.findQuotationEnd(sv, firstEndQuote);
  BOOST_TEST(qsv.empty());
  
  // special feature: in case there is no match, we still point to the end
  BOOST_TEST((qsv.begin() == sv.end()));
  
} // findQuotationEnd_noquote_test()


// -----------------------------------------------------------------------------
void isQuotationUnclosed_test() {
  
  using namespace std::string_literals;
  
  std::string_view sv {
    R"(No quote. Still no \"quote\". "a" "b c" 'd' "e 'f' g" "h 'i" "j' k" "unfinished 'm'".)"
    };
  
  icarus::ParsingToolkit tk; // default configuration
  
  BOOST_TEST(!tk.isQuotationUnclosed(sv));
  sv.remove_suffix(1);
  BOOST_TEST(!tk.isQuotationUnclosed(sv));
  sv.remove_suffix(1);
  BOOST_TEST( tk.isQuotationUnclosed(sv));
  sv.remove_suffix(1);
  BOOST_TEST( tk.isQuotationUnclosed(sv));
  sv.remove_suffix(1);
  BOOST_TEST( tk.isQuotationUnclosed(sv));
  sv.remove_suffix(1);
  BOOST_TEST( tk.isQuotationUnclosed(sv));
  
} // isQuotationUnclosed_test()


// -----------------------------------------------------------------------------
void splitOn_test() {
  
  std::string_view const sv { "aa:bbb" };
  // key must always be a substring of sv (not checked)
  {
    std::string_view const key = sv.substr(2,1);
    BOOST_TEST_REQUIRE(key == ":");
    auto const [ pre, sep, post ] = icarus::ParsingToolkit::splitOn(sv, key);
    BOOST_TEST(sep == key);
    BOOST_TEST(pre == "aa");
    BOOST_TEST(post == "bbb");
  }
  {
    std::string_view const key = sv.substr(0,2);
    BOOST_TEST_REQUIRE(key == "aa");
    auto const [ pre, sep, post ] = icarus::ParsingToolkit::splitOn(sv, key);
    BOOST_TEST(sep == key);
    BOOST_TEST(pre == "");
    BOOST_TEST(post == ":bbb");
  }
  
  {
    std::string_view const key = sv.substr(4, 1);
    BOOST_TEST_REQUIRE(key == "b");
    auto const [ pre, sep, post ] = icarus::ParsingToolkit::splitOn(sv, key);
    BOOST_TEST(sep == key);
    BOOST_TEST(pre == "aa:b");
    BOOST_TEST(post == "b");
  }
  
  {
    std::string_view const key
      = icarus::ParsingToolkit::make_view(sv.end(), sv.end());
    BOOST_TEST_REQUIRE(key == "");
    auto const [ pre, sep, post ] = icarus::ParsingToolkit::splitOn(sv, key);
    BOOST_TEST(sep == key);
    BOOST_TEST(pre == sv);
    BOOST_TEST(post == "");
  }
  
  {
    std::string_view const key
      = icarus::ParsingToolkit::make_view(sv.begin(), sv.begin());
    BOOST_TEST_REQUIRE(key == "");
    auto const [ pre, sep, post ] = icarus::ParsingToolkit::splitOn(sv, key);
    BOOST_TEST(sep == key);
    BOOST_TEST(pre == "");
    BOOST_TEST(post == sv);
  }
  
  {
    std::string_view const key = sv;
    auto const [ pre, sep, post ] = icarus::ParsingToolkit::splitOn(sv, key);
    BOOST_TEST(sep == key);
    BOOST_TEST(pre == "");
    BOOST_TEST(post == "");
  }
  
} // splitOn_test()


// -----------------------------------------------------------------------------
void removeEscapes_test() {
  
  icarus::ParsingToolkit tk; // default configuration
  
  std::string s { R"()" };
  
  BOOST_TEST(tk.removeWordEscapes(R"()") == R"()");
  
  BOOST_TEST(tk.removeWordEscapes(R"(\a)") == R"(a)");
  BOOST_TEST(tk.removeWordEscapes(R"(\a\h\a)") == R"(aha)");
  BOOST_TEST(tk.removeWordEscapes(R"(\a\\h\a)") == R"(a\ha)");
  BOOST_TEST(tk.removeWordEscapes(R"(aha)") == R"(aha)");
  BOOST_TEST(tk.removeWordEscapes(R"(aha\)") == R"(aha\)");
  
} // removeEscapes_test()


// -----------------------------------------------------------------------------
void removeEscapesDocumentation_test() {
  
  // "\\\\a" -> "\\a" -> "\a" -> "a" -> "a"
  
  icarus::ParsingToolkit tk; // default configuration
  
  std::string s { R"(\\\\a)" };
  
  s = tk.removeWordEscapes(s);
  BOOST_TEST(s == R"(\\a)");
  
  s = tk.removeWordEscapes(s);
  BOOST_TEST(s == R"(\a)");
  
  s = tk.removeWordEscapes(s);
  BOOST_TEST(s == R"(a)");
  
  s = tk.removeWordEscapes(s);
  BOOST_TEST(s == R"(a)");
  
} // removeEscapesDocumentation_test()


// -----------------------------------------------------------------------------
void removeQuotations_test() {
  
  using namespace std::string_literals;
  
  icarus::ParsingToolkit tk; // default configuration
  
  BOOST_TEST(tk.removeWordQuotations(R"()") == R"()");
  
  std::string s = tk.removeWordQuotations(R"("a'b"c'b"a)");
  BOOST_TEST(s == R"(a'bc'b"a)");
  
  s = tk.removeWordQuotations(s);
  BOOST_TEST(s == R"(abcb"a)");
  
  s = tk.removeWordQuotations(s);
  BOOST_TEST(s == R"(abcb"a)");
  
  
  s = tk.removeWordQuotations(R"("a'b""c'b"a)");
  BOOST_TEST(s == R"(a'bc'ba)");
  
  
  s = tk.removeWordQuotations(R"("a'b\"c'b"a)");
  BOOST_TEST(s == R"(a'b\"c'ba)");
  
  s = tk.removeWordQuotations(s);
  BOOST_TEST(s == R"(ab\"cba)");
  
  s = tk.removeWordQuotations(s);
  BOOST_TEST(s == R"(ab\"cba)");
  
} // removeQuotations_test()


void removeQuotationsDocumentation_test() {
  
  using namespace std::string_literals;
  
  icarus::ParsingToolkit tk; // default configuration
  /*
   * example, `a1 << "b1 << 'c1 << " or " << c2' << b2" << a2` will become first
   * `a1 << b1 << 'c1 << " or " << c2' << b2 << a2`, then
   * `a1 << b1 << c1 << " or " << c2 << b2 << a2`, and eventually
   * `a1 << b1 << c1 <<  or  << c2 << b2 << a2`).
   */
  
  std::string s { R"(a1 << "b1 << 'c1 << " or " << c2' << b2" << a2)" };
  
  s = tk.removeWordQuotations(s);
  BOOST_TEST(s == R"(a1 << b1 << 'c1 <<  or  << c2' << b2 << a2)");
  
  s = tk.removeWordQuotations(s);
  BOOST_TEST(s == R"(a1 << b1 << c1 <<  or  << c2 << b2 << a2)");
  
  s = tk.removeWordQuotations(s);
  BOOST_TEST(s == R"(a1 << b1 << c1 <<  or  << c2 << b2 << a2)");
  
} // removeQuotationsDocumentation_test()


// -----------------------------------------------------------------------------
// BEGIN Test cases  -----------------------------------------------------------
// -----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(isCharacterEscaped_testcase) {
  
  isCharacterEscaped_test();
  
} // BOOST_AUTO_TEST_CASE(isCharacterEscaped_testcase)


BOOST_AUTO_TEST_CASE(findNextCharacter_testcase) {
  
  findNextBlank_test();
  
} // BOOST_AUTO_TEST_CASE(findNextCharacter_testcase)


BOOST_AUTO_TEST_CASE(removeTrailingCharacters_testcase) {
  
  removeTrailingBlanks_test();
  
} // BOOST_AUTO_TEST_CASE(removeTrailingCharacters_testcase)


BOOST_AUTO_TEST_CASE(readMultiline_testcase) {
  
  readMultiline_test();
  readMultiline_endOfInput_test();
  
} // BOOST_AUTO_TEST_CASE(readMultiline_testcase)


BOOST_AUTO_TEST_CASE(splitWords_testcase) {
  
  splitWords_test();
  
} // BOOST_AUTO_TEST_CASE(splitWords_testcase)


BOOST_AUTO_TEST_CASE(findFirstUnescaped_testcase) {
  
  findFirstUnescaped_test();
  
} // BOOST_AUTO_TEST_CASE(findFirstUnescaped_testcase)


BOOST_AUTO_TEST_CASE(findFirstUnquoted_testcase) {
  
  findFirstUnquoted_test();
  
} // BOOST_AUTO_TEST_CASE(findFirstUnquoted_testcase)


BOOST_AUTO_TEST_CASE(removeCommentLine_testcase) {
  
  removeCommentLine_test();
  
} // BOOST_AUTO_TEST_CASE(removeCommentLine_testcase)


BOOST_AUTO_TEST_CASE(findQuotationStart_testcase) {
  
  findQuotationStart_test();
  findQuotationStart_noquote_test();
  
} // BOOST_AUTO_TEST_CASE(findQuotationStart_testcase)


BOOST_AUTO_TEST_CASE(findQuotationEnd_testcase) {
  
  findQuotationEnd_test();
  findQuotationEnd_noquote_test();
  
} // BOOST_AUTO_TEST_CASE(findQuotationEnd_testcase)


BOOST_AUTO_TEST_CASE(isQuotationUnclosed_testcase) {
  
  isQuotationUnclosed_test();
  
} // BOOST_AUTO_TEST_CASE(isQuotationUnclosed_testcase)


BOOST_AUTO_TEST_CASE(splitOn_testcase) {
  
  splitOn_test();
  
} // BOOST_AUTO_TEST_CASE(splitOn_testcase)


BOOST_AUTO_TEST_CASE(removeEscapes_testcase) {
  
  removeEscapes_test();
  removeEscapesDocumentation_test();
  
} // BOOST_AUTO_TEST_CASE(removeEscapes_testcase)


BOOST_AUTO_TEST_CASE(removeQuotations_testcase) {
  
  removeQuotations_test();
  removeQuotationsDocumentation_test();
  
} // BOOST_AUTO_TEST_CASE(removeQuotations_testcase)


// -----------------------------------------------------------------------------
// END Test cases  -------------------------------------------------------------
// -----------------------------------------------------------------------------
