/**
 * @file   icaruscode/Decode/DecoderTools/TriggerConfigurationExtractor.h
 * @brief  Utility to extract Trigger readout configuration from data.
 * @author Andrea Scarpelli (ascarpell@bnl.gov)
 * @date   March 23 2022
 * @file   icaruscode/Decode/DecoderTools/TriggerConfigurationExtractor.h
 * 
 */

#ifndef ICARUSCODE_DECODE_DECODERTOOLS_TRIGGERCONFIGURATIONEXTRACTOR_H
#define ICARUSCODE_DECODE_DECODERTOOLS_TRIGGERCONFIGURATIONEXTRACTOR_H

// ICARUS libraries
#include "icaruscode/Utilities/ReadArtConfiguration.h" // util::readConfigurationFromArtPrincipal()
#include "icaruscode/Decode/DataProducts/TriggerConfiguration.h"
#include "sbnobj/Common/Trigger/BeamBits.h"

// framework libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib_except/exception.h"

// ROOT libraries
#include "TFile.h"

// C/C++ standard libraries
#include <regex>
#include <string>
#include <optional>
#include <utility> // std::move(), std::pair<>
#include <initializer_list>

// -----------------------------------------------------------------------------
namespace icarus {
  
  class TriggerConfigurationExtractorBase;
  class TriggerConfigurationExtractor;
  
  icarus::TriggerConfiguration extractTriggerReadoutConfiguration
    (std::string const& srcFileName, icarus::TriggerConfigurationExtractor extractor);
  icarus::TriggerConfiguration extractTriggerReadoutConfiguration
    (TFile& srcFile, icarus::TriggerConfigurationExtractor extractor);
  template <typename Principal>
  icarus::TriggerConfiguration extractTriggerReadoutConfiguration
    (Principal const& data, icarus::TriggerConfigurationExtractor extractor);
  
} // namespace icarus


// -----------------------------------------------------------------------------
class icarus::TriggerConfigurationExtractorBase {
  
    public:
  
  using ConfigurationData_t = icarus::TriggerConfiguration;
  
  // --- BEGIN -- Interface ----------------------------------------------------
  /// @name Interface
  /// @{
  
  /// Returns whether `pset` may contain the needed configuration.
  static bool mayHaveConfiguration(fhicl::ParameterSet const& pset)
    { return pset.has_key("configuration_documents"); }
  
  
  /**
   * @brief Extracts all supported Trigger configuration from `config`.
   * @param config a FHiCL parameter set with component configuration
   * @return an object with the supported Trigger configuration
   * 
   * All Trigger-related configuration that is known to this code is extracted and
   * returned.
   * 
   * This function is undefined here: it must be overridden.
   */
  ConfigurationData_t extract(fhicl::ParameterSet const& config) const;
    
  /// @}
  // --- END ---- Interface ----------------------------------------------------
  
  
  // --- BEGIN -- Utility ------------------------------------------------------
  /// @name Utility
  /// @{
  
  /**
   * @brief Returns a parameter set with the content of
   *        `configuration_documents` key from `container`.
   * @param container parameter set including a table with key `configListKey`
   * @param configListKey name of the key in `container` with the configuration
   * @param components keys to be converted (as regular expressions)
   * @return a parameter set
   * 
   * The `configuration_documents` element of `container` is processed: for
   * each of its keys which match at least one of the `components` regular
   * expression patterns (`std::regex_match()`), the associated string value
   * is parsed with FHiCL parser, and the result is set as a FHiCL table in
   * the output parameter set.
   * For example, if the `components` are
   * `{ std::regex{".*pmt.*"}, std::regex{".*trigger.*"} }`, the returned
   * value is a parameter set that may have keys like `icaruspmtee01`,
   * `icaruspmtew02`, `icarustrigger` etc., each one with a FHiCL table as
   * `value.
   */
  static fhicl::ParameterSet convertConfigurationDocuments(
    fhicl::ParameterSet const& container,
    std::string const& configListKey,
    std::initializer_list<std::regex const> components
    );
  
  /// Returns whether `key` matches at least one of the regular expressions
  /// in the [ `rbegin`, `rend` [ range.
  template <typename RBegin, typename REnd>
  static bool matchKey(std::string const& key, RBegin rbegin, REnd rend);
  
  /// @}
  // --- END ---- Utility ------------------------------------------------------
  
}; // icarus::TriggerConfigurationExtractorBase


// -----------------------------------------------------------------------------
/**
 * @brief Class to extract PMT readout board configuration.
 * 
 * This is an example of trigger configuration taken from a pre-Run 1 ICARUS run
 * (`destinations` and `metrics` have been omitted):
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * daq: {
 *   fragment_receiver: {
 *  
 *      ....
 *  
 *      fpga_init_params: {
 *         MajLevelBeamCryo1: {
 *            name: "Majority level inside Beam Gate cryo1"
 *            value: "5"
 *         }
 *         MajLevelBeamCryo2: {
 *            name: "Majority level inside Beam Gate cryo2"
 *            value: "5"
 *         }
 *         MajLevelEnableCryo1: {
 *            name: "Majority level inside enable and outside Beam Gate cryo1"
 *            value: "10"
 *         }
 *         MajLevelEnableCryo2: {
 *            name: "Majority level inside enable and outside Beam Gate cryo2"
 *            value: "10"
 *         }
 *         MajorityTriggerType: {
 *            name: "Majority Trigger Type"
 *            value: "majo1or2"
 *         }
 *         RunType: {
 *            name: "Run type"
 *            value: "Majority"
 *         }
 *         SlidingWindowCryo1: {
 *            name: "Sliding window option cryo1"
 *            value: "Separated Window"
 *         }
 *         SlidingWindowCryo2: {
 *            name: "Sliding window option cryo2"
 *            value: "Separated Window"
 *         }
 *         Veto: {
 *            name: "Veto beam gate [ns]"
 *            value: "4000"
 *         }
 *      }
 *      
 *      ....
 *  
 *      generator: "ICARUSTriggerUDP"
 *  
 *      ....    
 *  
 *      spexi_init_params: {
 *         BNB1DOffset: {
 *            name: "BNB ($1D) Early Warning Offset"
 *            value: "0"
 *         }
 *         BNBBESOffset: {
 *            name: "BNB (gatedBES) Early Warning Offset"
 *            value: "0"
 *         }
 *         BNBBeamWidth: {
 *            name: "BNB Beam Gate Width [ns]"
 *            value: "8000"
 *         }
 *         BNBEnableWidth: {
 *            name: "BNB Enable Gate Width [ns]"
 *            value: "2000000"
 *         }
 *         GateSelection: {
 *            name: "Gate Selection [hex]"
 *            value: "FF"
 *         }
 *         MI12Ethertype: {
 *            name: "MI12 Ethertype"
 *            value: "5752"
 *         }
 *         MI12MACAddressLSB: {
 *            name: "MI12 MAC Address (LSB)"
 *            value: "30E971C8"
 *         }
 *         MI12MACAddressMSB: {
 *            name: "MI12 MAC Address (MSB)"
 *            value: "0800"
 *         }
 *         MI60Ethertype: {
 *            name: "MI60 Ethertype"
 *            value: "5752"
 *         }
 *         MI60MACAddressLSB: {
 *            name: "MI60 MAC Address (LSB)"
 *            value: "30E93B4F"
 *         }
 *         MI60MACAddressMSB: {
 *            name: "MI60 MAC Address (MSB)"
 *            value: "0800"
 *         }
 *         NuMIADOffset: {
 *            name: "NUMI ($AD) Early Warning Offset"
 *            value: "0"
 *         }
 *         NuMIBeamWidth: {
 *            name: "NuMI Beam Gate Width [ns]"
 *            value: "18000"
 *         }
 *         NuMIEnableWidth: {
 *            name: "NuMI Enable Gate Width [ns]"
 *            value: "2000000"
 *         }
 *         NuMIMIBSOffset: {
 *            name: "NuMI (MIBS$74) Early Warning Offset"
 *            value: "0"
 *         }
 *         OffBeamBNBBeamWidth: {
 *            name: "OffBeam-BNB Gate Width [ns]"
 *            value: "8000"
 *         }
 *         OffBeamBNBEnableWidth: {
 *            name: "OffBeam-BNB Enable Gate Width [ns]"
 *            value: "2000000"
 *         }
 *         OffBeamGateRate: {
 *            name: "OffBeam Gate Rate"
 *            value: "00010003"
 *         }
 *         OffBeamNuMIBeamWidth: {
 *            name: "OffBeam-NuMI Gate Width [ns]"
 *            value: "18000"
 *         }
 *         OffBeamNuMIEnableWidth: {
 *            name: "OffBeam-NuMI Enable Gate Width [ns]"
 *            value: "2000000"
 *         }
 *         PPSWidth: {
 *            name: "PPS Pulse Width [ns]"
 *            value: "1000"
 *         }
 *         PreScaleBNBNuMI: {
 *            name: "Prescale BNB-NuMI Number"
 *            value: "00010001"
 *         }
 *         PreScaleOffBeam: {
 *            name: "Prescale OffBeam Number"
 *            value: "00020002"
 *         }
 *         PrescaleZeroBias: {
 *            name: "Prescale ZeroBias Number"
 *            value: "0064"
 *         }
 *         TPCTriggerDelay: {
 *            name: "TPC Trigger Delay/400"
 *            value: "3250"
 *         }
 *         ZeroBiasEnableWidth: {
 *            name: "ZeroBias Enable Gate Width [ns]"
 *            value: "2000000"
 *         }
 *         ZeroBiasFreq: {
 *            name: "ZeroBias Frequency Period [ns]"
 *            value: "1000000000"
 *         }
 *         ZeroBiasWidth: {
 *            name: "ZeroBias Gate Width [ns]"
 *            value: "1600"
 *         }
 *      }
 *      use_wr_time: true
 *      wr_time_offset_ns: 1e9
 *    }
 *  
 *    ....     
 *  
 * }
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * 
 * The extractor identifies the correct configuration by the
 * `daq.fragment_receiver.generator` configuration key (in the above example,
 * it is `ICARUSTriggerUDP`). The key is specified as a regular expression
 * pattern that must match the whole string (as in `std::regex_match()` with
 * standard format).
 * 
 */
class icarus::TriggerConfigurationExtractor
  : public icarus::TriggerConfigurationExtractorBase
{
  
    public:
  
  /// Learns the name of the trigger fragment type.
  TriggerConfigurationExtractor
    (std::string expectedFragmentType = "ICARUSTrigger.*")
    : fExpectedFragmentTypeSpec{ std::move(expectedFragmentType) }
    , fExpectedFragmentType{ fExpectedFragmentTypeSpec }
    {}

  // --- BEGIN -- Interface ----------------------------------------------------
  /// @name Interface
  /// @{
  
  /// Returns whether `pset` may contain the needed configuration.
  static bool mayHaveConfiguration(fhicl::ParameterSet const& pset)
    { return pset.has_key("configuration_documents"); }
  
  /// Returns whether the specified `key` of `pset` is a good configuration.
  static bool isGoodConfiguration
    (fhicl::ParameterSet const& pset, std::string const& key);
  
  
  /**
   * @brief Extracts all supported Trigger configuration from `config`.
   * @param config a FHiCL parameter set with component configuration
   * @return an object with the supported Trigger configuration
   * 
   * All Trigger-related configuration that is known to this code is extracted and
   * returned.
   */
  icarus::TriggerConfiguration extract(fhicl::ParameterSet const& config) const;
    
  /// @}
  // --- END ---- Interface ----------------------------------------------------
  
    private:
  
  /// Specification of the trigger type name pattern.
  std::string const fExpectedFragmentTypeSpec;
  
  /// Regular expression pattern for trigger type name.
  std::regex const fExpectedFragmentType;
  
  /// Regular expressions matching all names of supported Trigger configurations.
  static std::vector<std::regex> const ConfigurationNames;

  /**
   * @brief Extracts trigger readout board configuration from `pset`.
   * @param pset information source (FHiCL configuration)
   * @param boardName name of the board we are given the configuration of
   * @return the requested configuration
   * 
   * 
   */
  icarus::TriggerConfiguration extractTriggerConfiguration
    (fhicl::ParameterSet const& pset)
    const;
      
  /**
   * @brief Returns the specified Trigger readout board configuration.
   * @param pset parameter set including `key`
   * @param key key of the Trigger readout configuration candidate
   * @return the configuration, or an empty object if key does not represent one
   */
  std::optional<fhicl::ParameterSet> readTriggerConfig
    (fhicl::ParameterSet const& pset, std::string const& key) const;

   unsigned long parsePrescaleString( std::string prescaleString, std::size_t source ) const;

   unsigned int parseWindowMode(std::string bitStr) const;
  
}; // icarus::TriggerConfigurationExtractor



// -----------------------------------------------------------------------------
// ---  Template implementation
// -----------------------------------------------------------------------------
namespace icarus::details {
  
  template <typename ConfigMap>
  icarus::TriggerConfiguration extractTriggerReadoutConfigurationImpl
    (ConfigMap const& configMap, icarus::TriggerConfigurationExtractor extractor)
  {
    
    /*
    * Requirements: ConfigMap is a mapping type supporting iteration and whose
    * elements support structured assignment into a pair, the first element
    * being an identifier (process name, parameter set ID...) and the second
    * being a `fhicl::ParameterSet` with the configuration (or something
    * offering the same interface).
    * 
    * The plan is to look in all the FHiCL configuration fragments we can find
    * in the input config, and find all the useful configuration therein.
    * Given that there may be multiple input files, there may also be multiple
    * configurations for the same detector components.
    * In that case, we will extract parameters from each and every one of the
    * configurations, and throw an exception if they are not all consistent.
    * 
    * Consistency is tested only for the extracted parameters, not for the whole
    * FHiCL configuration fragment.
    */
    
    using Key_t = std::tuple_element_t<0U, typename ConfigMap::value_type>;
    
    std::optional<std::pair<Key_t, icarus::TriggerConfiguration>> config;
    
    // look in the global configuration for all parameter sets which contain
    // `configuration_documents` as a (direct) name;
    for (auto const& [ id, pset ]: configMap) {
      if (!extractor.mayHaveConfiguration(pset)) continue;
      
      fhicl::ParameterSet const configDocs
        = extractor.convertConfigurationDocuments
          (pset, "configuration_documents", { std::regex{ "icarustrigger.*" } })
        ;
      
      icarus::TriggerConfiguration candidateConfig = extractor.extract(configDocs);
      if (config) {
        if (config->second == candidateConfig) continue;
        mf::LogError log("extractTriggerReadoutConfiguration");
        log << "Found two candidate configurations differring:"
          "\nFirst [" << config->first << "]:\n" << config->second
          << "\nSecond [" << id << "]:\n" << candidateConfig
          ;
        throw cet::exception("extractTriggerReadoutConfiguration")
          << "extractTriggerReadoutConfiguration() found inconsistent configurations.\n";
      } // if incompatible configurations
      
      config.emplace(std::move(id), std::move(candidateConfig));
    } // for all configuration documents
    
    if (!config) {
      throw cet::exception("extractTriggerReadoutConfiguration")
        << "extractTriggerReadoutConfiguration() could not find a suitable configuration.\n";
    }
    
    return std::move(config->second);
  } // extractPMTreadoutConfigurationImpl(ConfigMap)
  
  
} // namespace icarus::details


// -----------------------------------------------------------------------------
template <typename RBegin, typename REnd>
bool icarus::TriggerConfigurationExtractorBase::matchKey
  (std::string const& key, RBegin rbegin, REnd rend)
{
  for (auto iRegex = rbegin; iRegex != rend; ++iRegex)
    if (std::regex_match(key, *iRegex)) return true;
  return false;
} // icarus::PMTconfigurationExtractorBase::matchKey()


// -----------------------------------------------------------------------------
template <typename Principal>
icarus::TriggerConfiguration icarus::extractTriggerReadoutConfiguration
  (Principal const& principal, icarus::TriggerConfigurationExtractor extractor)
{
  
  return details::extractTriggerReadoutConfigurationImpl
    (util::readConfigurationFromArtPrincipal(principal), std::move(extractor));
  
} // icarus::extracatTriggerReadoutConfiguration(Principal)


// ---------------------------------------------------------------------------


#endif // ICARUSCODE_DECODE_DECODERTOOLS_TRIGGERCONFIGURATIONEXTRACTOR_H
