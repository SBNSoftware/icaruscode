/**
 * @file   icaruscode/Decode/DecoderTools/triggerPacketParser.cxx
 * @brief  Simple parser for trigger raw data packets.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   March 4, 2022
 */

// ICARUS libraries
#include "icaruscode/Decode/DecoderTools/details/KeyedCSVparser.h"


// C++/Boost libraries
#include "boost/program_options.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <optional>
#include <cstdlib> // std::exit()

/*
 * Notable changes here:
 * 
 * [20220304] (petrillo@slac.stanford.edu) [1.0]
 *     initial version
 * 
 */
static std::string const ProgramVersion = "v1.0";

// -----------------------------------------------------------------------------
boost::program_options::variables_map parseCommandLine(int argc, char** argv) {
  
  namespace po = boost::program_options;
  
  //
  // Declare the supported options.
  //
  po::options_description inputopt("Input/output");
  inputopt.add_options()
    ("input", po::value<std::vector<std::string>>(), "input file")
    ;

  po::options_description genopt("General");
  genopt.add_options()
    ("help,?", "print usage instructions and exit")
    ("version,V", "print usage instructions")
    ;
  
  po::options_description allopt("Options");
  allopt.add(inputopt).add(genopt);

  po::positional_options_description pos;
  pos.add("input", -1); // all positional parameters get option name "input"
  
  //
  // proceed with parsing
  //
  po::variables_map optmap;
  po::store(
    po::command_line_parser(argc, argv)
      .options(allopt).positional(pos).run(),
    optmap
    );
  po::notify(optmap);

  //
  // deal with general options
  //
  std::optional<int> exitWithCode;
  if (optmap.count("verbose")) {
    // C++17 (with supporting Clang/GCC): use std::filesystem::path{ argv[0] }.filename()
    std::cout << argv[0] << " version " << ProgramVersion << std::endl;
    exitWithCode = 0;
  }
  if (optmap.count("help")) {
    std::cout
      <<   "Parses the trigger packet strings from the input file(s) and prints its interpretation."
      << "\nIf no file is specified in the command line (option or positional argument),"
      << "\nstrings are read from the standard input stream."
      << "\n" << allopt
      << std::endl
      ;
    exitWithCode = 0;
  }
  
  if (exitWithCode) std::exit(*exitWithCode);
  return optmap;
  
} // parseCommandLine()


// -----------------------------------------------------------------------------
int processTriggerData(std::string const& triggerString) {
  
  //
  // parse
  //
  
  icarus::details::KeyedCSVparser parser;
  
  parser.addPatterns({
      { "Cryo. (EAST|WEST) Connector . and .", 1U }
    , { "Trigger Type", 1U }
    });
  
  icarus::details::KeyValuesData parsedData;
  try {
    parsedData = parser(triggerString);
  }
  catch (icarus::details::KeyedCSVparser::Error const& e) {
    
    std::cerr << "Error parsing trigger data:\n" << std::string(80, '-')
      << "\n" << triggerString << "\n" << std::string(80, '-')
      << "Error: " << e.what() << std::endl;
    return 1;
    
  }
  
  
  //
  // dump
  //
  
  // for this printout, all keys are treated as lists of integers or,
  // if that conversion fails, vectors of strings;
  // exceptions are listed here:
  
  std::map<std::string, std::vector<std::string>> typeKeys = {
    {
      "hex64", {
          "Cryo1 EAST Connector 0 and 1"
        , "Cryo1 EAST Connector 2 and 3"
        , "Cryo2 WEST Connector 0 and 1"
        , "Cryo2 WEST Connector 2 and 3"
      }
    }
    };
  
  // reversed map
  std::map<std::string, std::string> keyType;
  for (auto const& [ type, keys ]: typeKeys)
    for (std::string const& key: keys) keyType[key] = type;
  
  std::cout
    << "Trigger data (" << triggerString.length() << " char):"
    << "\n" << triggerString
    << "\nParsed as:";
  for (icarus::details::KeyValuesData::Item const& item: parsedData.items()) {
    
    std::cout << "\n '" << item.key() << "':";
    
    std::string type = keyType.count(item.key())? keyType[item.key()]: "auto";
    if (type.empty()) type = "auto";
    
    if (type == "hex64") {
      
      std::vector<std::uint64_t> const& values = item.getVector<std::uint64_t>(
        icarus::details::KeyValuesData::Item::StringConverter<std::uint64_t>(16)
        );
      if (values.empty()) std::cout << " <no number>";
      else if (values.size() == 1) std::cout << " " << std::hex << values[0];
      else {
        auto iValue = values.begin(), vend = values.end();
        std::cout << " (" << values.size() << " numbers)  "
          << std::hex << *iValue;
        while (++iValue != vend) std::cout << " , " << *iValue;
      }
      std::cout << std::dec << "  (" << type << ")";
      continue;
    } // hex
    
    if (type != "str") {
      
      try { // we try `int`
        
        std::vector<int> values = item.getVector<int>();
        if (values.empty()) std::cout << " <no number>";
        else if (values.size() == 1) std::cout << " " << values[0];
        else {
          auto iValue = values.begin(), vend = values.end();
          std::cout << " (" << values.size() << " numbers)  " << *iValue;
          while (++iValue != vend) std::cout << " , " << *iValue;
        }
        std::cout << "  (" << type << ")";
        continue;
        
      }
      catch (icarus::details::KeyValuesData::Error const& e) { // ... nope
        if (type == "int") {
          std::cout << " <error: not int>";
        }
        else if (type != "auto") type = "?";
      }
    }
    
    // as strings, at last
    std::vector<std::string> const& values = item.values();
    if (values.empty()) std::cout << " <no value>";
    else if (values.size() == 1) std::cout << " " << values[0];
    else {
      auto iValue = values.begin(), vend = values.end();
      std::cout << " (" << values.size() << " values)  " << *iValue;
      while (++iValue != vend) std::cout << " , " << *iValue;
    }
    
    std::cout << "  (" << type << ")";
    
  } // for
  std::cout << std::endl;
  
  return 0; // happy
  
} // processTriggerData()


// -----------------------------------------------------------------------------
int main(int argc, char** argv) {
  
  boost::program_options::variables_map const options
    = parseCommandLine(argc, argv);
  
  std::vector<std::string> inputFilePaths;
  if (options.count("input"))
    inputFilePaths = options["input"].as<std::vector<std::string>>();
  if (inputFilePaths.empty()) inputFilePaths.push_back("");
  
  unsigned int nErrors = 0U;
  unsigned int lineCount = 0U;
  for (std::string const& inputFilePath: inputFilePaths) {
    
    // select and open the input file
    std::optional<std::ifstream> realInputFile;
    if (inputFilePath.empty())
      std::clog << "Reading data from standard input." << std::endl;
    else {
      realInputFile.emplace(inputFilePath);
      if (!realInputFile->good()) {
        std::cerr << "FATAL: can't open input file '" << inputFilePath << "'."
          << std::endl;
        return 2;
      }
      std::clog << "Reading data from '" << inputFilePath << "'." << std::endl;
    }
    std::istream& inputFile = realInputFile? *realInputFile: std::cin;
    
    unsigned int fileLineCount = 0U;
    std::string line;
    while (std::getline(inputFile, line)) {
      ++fileLineCount;
      
      if (processTriggerData(line) != 0) ++nErrors;
      
    } // while
    
    lineCount += fileLineCount;
    
    
  } // for inputFilePaths
  
  if (nErrors > 0) {
    std::cerr << "Parsing failed for " << nErrors << "/" << lineCount
      << " lines." << std::endl;
  }
  
  return (nErrors == 0)? 0: 1;
} // main()


// -----------------------------------------------------------------------------
