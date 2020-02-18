///////////////////////////////////////////////////////////////////////
// $Id: SimWireICARUS.cxx,v 1.22 2010/04/23 20:30:53 seligman Exp $
//
// SimWireICARUS class designed to simulate signal on a wire in the TPC
//
// katori@fnal.gov
//
// - Revised to use sim::RawDigit instead of rawdata::RawDigit, and to
// - save the electron clusters associated with each digit.
//
////////////////////////////////////////////////////////////////////////

/**
 * If defined, a hack to make sure DetectorClocksService knows about the new
 * hardware trigger time is enabled.
 * This is violating art/LArSoft recommended practices, and it is not even
 * useful in ICARUS where the
 * @ref DetectorClocksElectronicsStartTime "electronics time start"
 * is _determined_ by the hardware trigger.
 */
#undef ICARUSCODE_SIMWIREICARUS_TRIGGERTIMEHACK


// C/C++ standard library
#include <stdexcept> // std::range_error
#include <vector>
#include <string>
#include <algorithm> // std::fill()
#include <functional>
#include <cassert>
// CLHEP libraries
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Random/JamesRandom.h"
// ROOT libraries
#include "TMath.h"
#include "TComplex.h"
#include "TString.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TFile.h"
#include "TCanvas.h"
// art library and utilities
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art/Utilities/make_tool.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/DelegatedParameter.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/TupleAs.h"
#include "fhiclcpp/types/OptionalTupleAs.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
// art extensions
#include "nurandom/RandomUtils/NuRandomService.h" // `rndm` namespace
// LArSoft libraries
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/CoreUtils/StdUtils.h" // util::begin(), util::end()
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#ifdef ICARUSCODE_SIMWIREICARUS_TRIGGERTIMEHACK
#include "lardataobj/RawData/TriggerData.h"
#include "lardata/DetectorInfoServices/DetectorClocksServiceStandard.h" // FIXME: this is not portable
#endif // ICARUSCODE_SIMWIREICARUS_TRIGGERTIMEHACK
#include "icaruscode/TPC/Utilities/SignalShapingICARUSService_service.h"
#include "icaruscode/TPC/Utilities/ICARUSFFT.h"
#include "lardataalg/Utilities/StatCollector.h" // lar::util::MinMaxCollector<>
#include "larevt/CalibrationDBI/Interface/DetPedestalService.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalProvider.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
#include "tools/IGenNoise.h"
using namespace util;


// TODO move this into a larcoreobj/SimpleTypesAndConstants/geo_types_fhicl.h
#ifndef LARCOREOBJ_SIMPLETYPESANDCONSTANTS_GEO_TYPES_FHICL_H
#define LARCOREOBJ_SIMPLETYPESANDCONSTANTS_GEO_TYPES_FHICL_H


// LArSoft libraries
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

// support libraries
#include "fhiclcpp/types/OptionalTable.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/OptionalSequence.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Atom.h"

// C/C++ standard libraries
#include <vector>
#include <optional>
#include <cstddef> // std::size_t


/// FHiCL objects representing geometry classes as configuration parameters.
namespace geo::fhicl {
  
  // --- BEGIN -- Validated configuration parameters for geometry ID objects ---
  /**
   * @name Validated configuration parameters for geometry ID objects
   * 
   * These data types can be used in a class for validated FHiCL configuration.
   * They are implemented as configuration tables (`fhicl::Table`) of a
   * configuration structure containing one parameter (`fhicl::Atom`) per index
   * in the ID. They do _not_ support default values, but optional parameters
   * may be used as a workaround.
   * An ID described data member can be specified as a table with the same
   * syntax as the standard printout of the IDs, e.g. `{ C:1 T:3 P:2 }`
   * for the plane `C:1 T:3 P:2`.
   * Example:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * struct Config {
   *   
   *   geo::fhicl::PlaneIDsequence Planes {
   *     fhicl::Name("Planes"),
   *     fhicl::Comment("anode planes to process")
   *     };
   *   
   *   geo::fhicl::OptionalPlaneID ReferencePlane {
   *     fhicl::Name("ReferencePlane"),
   *     fhicl::Comment("reference anode plane (first one by default)")
   *     };
   *   
   * }; // struct Config
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * which can be configured as:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * Planes: [
   *   { C:0 T:1 P:0 },
   *   { C:0 T:1 P:1 },
   *   { C:0 T:1 P:2 }
   *   ]
   * ReferencePlane: { C:0 T:1 P:2 }
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * and read as:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * void readParams(art::EDProducer::Table<Config> const& config) {
   *   
   *   std::vector<geo::PlaneID> planes = readIDsequence(config().Planes);
   *   if (planes.empty()) {
   *     throw art::Exception(art::errors::Configuration)
   *       << "At least one plane is needed.\n";
   *   }
   *   
   *   geo::PlaneID refPlane; // invalid by default
   *   if (!config().ReferencePlane(refPlane)) refPlane = planes.front();
   *   
   * } // readParams()
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * 
   * Currently default values are not supported.
   * The workaround is to use the "optional" version of the objects.
   * For parameter sequences the use is a bit more complicate, and an utility
   * is provided for that:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * struct Config {
   *   
   *   geo::fhicl::OptionalPlaneIDsequence Planes {
   *     fhicl::Name("Planes"),
   *     fhicl::Comment("anode planes to process (omit or empty processes all)")
   *     };
   *   
   * }; // struct Config
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * reading as:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * void readParams(art::EDProducer::Table<Config> const& config) {
   *   
   *   std::vector<geo::PlaneID> planes
   *     = readOptionalIDsequence(config().Planes, {});
   *   
   * } // readParams()
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * Note however that the default value will not show in the regular
   * `lar --print-description` output.
   * 
   */
  /// @{
  
  /// Helper class holding the ID validity flag.
  struct ValidIDConfig {
    
    ::fhicl::Atom<bool> isValid {
      ::fhicl::Name("isValid"),
      ::fhicl::Comment("whether the ID is valid"),
      true // default
      };
    
    bool valid() const { return isValid(); }
    
  }; // struct ValidIDConfig
  
  
  // --- BEGIN -- Cryostat ID --------------------------------------------------
  /// Configuration structure for validated `geo::CryostatID` parameter.
  struct CryostatIDConfig: public ValidIDConfig {
    using ID_t = geo::CryostatID; ///< Type read by this configuration.
    
    ::fhicl::Atom<geo::CryostatID::CryostatID_t> C {
      ::fhicl::Name("C"),
      ::fhicl::Comment("cryostat number"),
      [this](){ return isValid(); }
      };
    
    ID_t ID() const { return ID_t{ C() }; }
    operator ID_t() const { return ID(); }
    
  }; // struct CryostatIDConfig
  
  /// Member type of validated `geo::CryostatID` parameter.
  using CryostatID = ::fhicl::Table<CryostatIDConfig>;
  
  /// Member type of optional validated `geo::CryostatID` parameter.
  using OptionalCryostatID = ::fhicl::OptionalTable<CryostatIDConfig>;
  
  /// Member type of sequence of `geo::CryostatID` parameters.
  using CryostatIDsequence = ::fhicl::Sequence<CryostatID>;
  
  /// Member type of optional sequence of `geo::CryostatID` parameters.
  using OptionalCryostatIDsequence = ::fhicl::OptionalSequence<CryostatID>;
  
  // --- END -- Cryostat ID ----------------------------------------------------
  
  
  // --- BEGIN -- TPC ID -------------------------------------------------------
  /// Configuration structure for validated `geo::TPCID` parameter.
  struct TPCIDConfig: public CryostatIDConfig {
    using ID_t = geo::TPCID; ///< Type read by this configuration.
    
    ::fhicl::Atom<geo::TPCID::TPCID_t> T {
      ::fhicl::Name("T"),
      ::fhicl::Comment("TPC number within the cryostat"),
      [this](){ return isValid(); }
      };
    
    ID_t ID() const { return { CryostatIDConfig::ID(), T() }; }
    operator ID_t() const { return ID(); }
  }; // struct TPCIDConfig
  
  /// Member type of validated `geo::TPCID` parameter.
  using TPCID = ::fhicl::Table<TPCIDConfig>;
  
  /// Member type of optional validated `geo::TPCID` parameter.
  using OptionalTPCID = ::fhicl::OptionalTable<TPCIDConfig>;
  
  /// Member type of sequence of `geo::TPCID` parameters.
  using TPCIDsequence = ::fhicl::Sequence<TPCID>;
  
  /// Member type of optional sequence of `geo::TPCID` parameters.
  using OptionalTPCIDsequence = ::fhicl::OptionalSequence<TPCID>;
  
  // --- END -- TPC ID ---------------------------------------------------------
  
  
  // --- BEGIN -- Plane ID -----------------------------------------------------
  /// Configuration structure for validated `geo::PlaneID` parameter.
  struct PlaneIDConfig: public TPCIDConfig {
    using ID_t = geo::PlaneID; ///< Type read by this configuration.
    
    ::fhicl::Atom<geo::PlaneID::PlaneID_t> P {
      ::fhicl::Name("P"),
      ::fhicl::Comment("Plane number within the TPC"),
      [this](){ return isValid(); }
      };
    
    ID_t ID() const { return { TPCIDConfig::ID(), P() }; }
    operator ID_t() const { return ID(); }
  }; // struct PlaneIDConfig
  
  /// Member type of validated `geo::PlaneID` parameter.
  using PlaneID = ::fhicl::Table<PlaneIDConfig>;
  
  /// Member type of optional validated `geo::PlaneID` parameter.
  using OptionalPlaneID = ::fhicl::OptionalTable<PlaneIDConfig>;
  
  /// Member type of sequence of `geo::PlaneID` parameters.
  using PlaneIDsequence = ::fhicl::Sequence<PlaneID>;
  
  /// Member type of optional sequence of `geo::PlaneID` parameters.
  using OptionalPlaneIDsequence = ::fhicl::OptionalSequence<PlaneID>;
  
  // --- END -- Plane ID -------------------------------------------------------
  
  
  // --- BEGIN -- Wire ID -----------------------------------------------------
  /// Configuration structure for validated `geo::PlaneID` parameter.
  struct WireIDConfig: public PlaneIDConfig {
    using ID_t = geo::WireID; ///< Type read by this configuration.
    
    ::fhicl::Atom<geo::WireID::WireID_t> W {
      ::fhicl::Name("W"),
      ::fhicl::Comment("Wire number within the plane"),
      [this](){ return isValid(); }
      };
    
    ID_t ID() const { return { PlaneIDConfig::ID(), W() }; }
    operator ID_t() const { return ID(); }
  }; // struct WireIDConfig
  
  /// Member type of validated `geo::WireID` parameter.
  using WireID = ::fhicl::Table<WireIDConfig>;
  
  /// Member type of optional validated `geo::WireID` parameter.
  using OptionalWireID = ::fhicl::OptionalTable<WireIDConfig>;
  
  /// Member type of sequence of `geo::WireID` parameters.
  using WireIDsequence = ::fhicl::Sequence<WireID>;
  
  /// Member type of optional sequence of `geo::WireID` parameters.
  using OptionalWireIDsequence = ::fhicl::OptionalSequence<WireID>;
  
  // --- END -- Wire ID -------------------------------------------------------
  
  
  // --- BEGIN -- ID sequence parsing ------------------------------------------
  
  //@{
  /// Type of the ID in the ID sequence.
  template <typename IDsequence>
  using IDofSequence = typename IDsequence::value_type::value_type::ID_t;
  //@}
  
  //@{
  /**
   * @brief Returns a vector of IDs extracted from the specified ID sequence.
   * @tparam IDsequence type of FHiCL sequence object
   * @tparam ID type of the element in the returned collection
   *            (default: the ID type of `IDsequence`)
   * @param seq the sequence of ID parameters to convert
   * @return a STL vector of `ID` objects converted from `seq` parameter values
   * 
   * This function returns the value of the specified FHiCL sequence object
   * (`fhicl::Sequence`). It supports both fixed and variable size sequences,
   * but it always returns a STL vector as a result.
   * 
   * Example of usage: the configuration object `Config` and the data member to
   * store the configuration parameter value are defined in a class as:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * struct Config {
   *   
   *   geo::fhicl::TPCIDsequence TPCs
   *     { fhicl::Name("TPCs"), fhicl::Comment("selected TPCs") };
   *   
   * };
   * 
   * std::vector<geo::TPCID> fTPCs;
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * The constructor of that class should have an entry in the initializer list
   * like:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   *   fTPCs(geo::fhicl::readIDsequence(config().TPCs))
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * (note that the argument is just `config().TPCs`, not `config().TPCs()`).
   * 
   * @note The additional template parameter `ID` is provided as an added bonus
   *       to choose which type to convert the configuration parameters into,
   *       and it's not enforced to be a ID type at all.
   */
  template <typename IDsequence, typename ID = IDofSequence<IDsequence>>
  std::vector<IDofSequence<IDsequence>> readIDsequence(IDsequence const& seq);
  //@}
  
  //@{
  /**
   * @brief Returns a vector of IDs extracted from the specified optional ID
   *        sequence.
   * @tparam IDsequence type of FHiCL optional sequence object
   * @tparam ID type of the element in the returned collection
   *            (default: the ID type of `IDsequence`)
   * @param seq the optional sequence of ID parameters to convert
   * @return an optional collection containing a STL vector of `ID` objects
   *         converted from `seq` parameter values, or no value if the parameter
   *         was omitted
   * 
   * This function returns the value of the specified FHiCL optional sequence
   * object (`fhicl::OptionalSequence`). It supports both fixed and variable
   * size optional sequences, but it always returns an optional STL vector as a
   * result.
   * 
   * Example of usage: the configuration object `Config` and the data member to
   * store the configuration parameter value are defined in a class as:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * struct Config {
   *   
   *   geo::fhicl::OptionalTPCIDsequence TPCs
   *     { fhicl::Name("TPCs"), fhicl::Comment("selected TPCs") };
   *   
   * };
   * 
   * std::optional<std::vector<geo::TPCID>> fTPCs;
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * The constructor of that class should have an entry in the initializer list
   * like:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   *   fTPCs(geo::fhicl::readIDsequence(config().TPCs))
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * (note that the argument is just `config().TPCs`, not `config().TPCs()`).
   * If instead a "default value" needs to be provided, the data member is
   * simply:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * std::vector<geo::TPCID> fTPCs;
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * and the value can be assigned via the standard `std::optional` interface:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   *   fTPCs(geo::fhicl::readIDsequence(config().TPCs).value_or(std::vector<geo::TPCID>{}))
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * (in this case the default value is an empty collection of TPC IDs) or using
   * a different overload of `readOptionalIDsequence()`:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   *   fTPCs(geo::fhicl::readIDsequence(config().TPCs, {}))
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * 
   * 
   * @note The additional template parameter `ID` is provided as an added bonus
   *       to choose which type to convert the configuration parameters into,
   *       and it's not enforced to be a ID type at all.
   */
  template <typename IDsequence, typename ID = IDofSequence<IDsequence>>
  std::optional<std::vector<IDofSequence<IDsequence>>>
    readOptionalIDsequence(IDsequence const& seq);
  //@}
  
  //@{
  /**
   * @brief Returns a vector of IDs extracted from the specified optional ID
   *        sequence, or a default value.
   * @tparam IDsequence type of FHiCL optional sequence object
   * @tparam ID type of the element in the returned collection
   *            (default: the ID type of `IDsequence`)
   * @param seq the optional sequence of ID parameters to convert
   * @param defValue value to be returned if the optional parameter was omitted
   * @return a collection containing a STL vector of `ID` objects
   *         converted either from `seq` parameter values or from `defValue`
   * 
   * This function is based on `readOptionalIDsequence(IDsequence const&)`.
   * The operating mode is the same, but if the value is not available from
   * the parameters, a copy of `defValue` is returned, or `defValue` content
   * is moved into the returned value.
   */
  template <typename IDsequence, typename ID = IDofSequence<IDsequence>>
  std::vector<IDofSequence<IDsequence>> readOptionalIDsequence(
    IDsequence const& seq,
    std::vector<IDofSequence<IDsequence>> const& defValue
    );
  
  template <typename IDsequence, typename ID = IDofSequence<IDsequence>>
  std::vector<IDofSequence<IDsequence>> readOptionalIDsequence(
    IDsequence const& seq,
    std::vector<IDofSequence<IDsequence>>&& defValue
    );
  //@}
  
  // --- END -- ID sequence parsing --------------------------------------------
  
  
  /// @}
  // --- END -- Validated configuration parameters for geometry ID objects -----
  
} // namespace geo::fhicl


// -----------------------------------------------------------------------------
// ---  template implementation
// -----------------------------------------------------------------------------
template<
  typename IDsequence,
  typename ID /* = geo::fhicl::IDofSequence<IDsequence> */
  >
auto geo::fhicl::readIDsequence(IDsequence const& seq)
  -> std::vector<IDofSequence<IDsequence>>
{
  using ID_t = ID;
  
  std::vector<ID_t> IDs;
  std::size_t const n = seq.size();
  IDs.reserve(n);
  for (std::size_t i = 0; i < n; ++i)
    IDs.push_back(seq(i)); // seq(i) is TPCIDConfig
  return IDs;
} // geo::fhicl::readIDsequence()


// -----------------------------------------------------------------------------
template <
  typename IDsequence,
  typename ID /* = geo::fhicl::IDofSequence<IDsequence> */
  >
auto geo::fhicl::readOptionalIDsequence(IDsequence const& seq)
  -> std::optional<std::vector<IDofSequence<IDsequence>>>
{
  using values_t = std::vector<ID>;
  
  typename IDsequence::value_type values;
  if (!seq(values)) return std::nullopt;
  
  values_t IDs;
  IDs.reserve(values.size());
  std::copy(values.begin(), values.end(), std::back_inserter(IDs));
  return { std::move(IDs) };
  
} // geo::fhicl::readOptionalIDsequence()


// -----------------------------------------------------------------------------
template <
  typename IDsequence,
  typename ID /* = geo::fhicl::IDofSequence<IDsequence> */
  >
auto geo::fhicl::readOptionalIDsequence
  (IDsequence const& seq, std::vector<IDofSequence<IDsequence>> const& defValue)
  -> std::vector<IDofSequence<IDsequence>>
{
  // making sure `paramValue` is not a r-value; not sure whether it is necessary
  auto paramValue = readOptionalIDsequence(seq);
  return paramValue.value_or(defValue);
} // geo::fhicl::readOptionalIDsequence(IDsequence, std::vector const&)


// -----------------------------------------------------------------------------
template <
  typename IDsequence,
  typename ID /* = geo::fhicl::IDofSequence<IDsequence> */
  >
auto geo::fhicl::readOptionalIDsequence
  (IDsequence const& seq, std::vector<IDofSequence<IDsequence>>&& defValue)
  -> std::vector<IDofSequence<IDsequence>>
{
  return readOptionalIDsequence(seq).value_or(std::move(defValue));
} // geo::fhicl::readOptionalIDsequence(IDsequence, std::vector const&)


// -----------------------------------------------------------------------------


#endif // LARCOREOBJ_SIMPLETYPESANDCONSTANTS_GEO_TYPES_FHICL_H

//
// TODO move this into a larcoreobj/SimpleTypesAndConstants/readout_types_fhicl.h
//
#ifndef LARCOREOBJ_SIMPLETYPESANDCONSTANTS_READOUT_TYPES_FHICL_H
#define LARCOREOBJ_SIMPLETYPESANDCONSTANTS_READOUT_TYPES_FHICL_H


// LArSoft libraries
// #include "larcoreobj/SimpleTypesAndConstants/geo_types_fhicl.h"
#include "larcoreobj/SimpleTypesAndConstants/readout_types.h"

// support libraries
#include "fhiclcpp/types/OptionalTable.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/OptionalSequence.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/Atom.h"


/// FHiCL objects representing readout classes as configuration parameters.
namespace readout::fhicl {
  
  // --- BEGIN -- Validated configuration parameters for readout ID objects ----
  /**
   * @name Validated configuration parameters for readout ID objects
   * 
   * These data types can be used in a class for validated FHiCL configuration.
   * They are implemented as configuration tables (`fhicl::Table`) of a
   * configuration structure containing one parameter (`fhicl::Atom`) per index
   * in the ID. They do _not_ support default values, but optional parameters
   * may be used as a workaround.
   * An ID described data member can be specified as a table with the same
   * syntax as the standard printout of the IDs, e.g. `{ C:1 S:3 R:2 }`
   * for the plane `C:1 S:3 R:2`.
   * Example:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * struct Config {
   *   
   *   geo::fhicl::ROPIDsequence ROPs {
   *     fhicl::Name("ROPs"),
   *     fhicl::Comment("readout planes to process")
   *     };
   *   
   *   geo::fhicl::OptionalROPID ReferenceROP {
   *     fhicl::Name("ReferenceROP"),
   *     fhicl::Comment("reference readout anode plane (first one by default)")
   *     };
   *   
   * }; // struct Config
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * which can be configured as:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * ROPs: [
   *   { C:0 S:1 R:0 },
   *   { C:0 S:1 R:1 },
   *   { C:0 S:1 R:2 }
   *   ]
   * ReferenceROP: { C:0 S:1 R:2 }
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * and read as:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * void readParams(art::EDProducer::Table<Config> const& config) {
   *   
   *   std::vector<geo::ROPID> ROPs = readIDsequence(config().ROPs);
   *   if (ROPs.empty()) {
   *     throw art::Exception(art::errors::Configuration)
   *       << "At least one readout plane is needed.\n";
   *   }
   *   
   *   geo::ROPID refROP; // invalid by default
   *   if (!config().ReferenceROP(refROP)) refROP = planes.front();
   *   
   * } // readParams()
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * 
   * The exception is `readout::fhicl::ChannelID`, which is provided for
   * completeness but is just a `fhicl::Atom` type and reads just as a plain
   * integral number.
   * 
   * Currently default values are not supported.
   * The workaround is to use the "optional" version of the objects.
   * For parameter sequences the use is a bit more complicate, and an utility
   * is provided for that:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * struct Config {
   *   
   *   geo::fhicl::OptionalROPIDsequence ROPs {
   *     fhicl::Name("ROPs"),
   *     fhicl::Comment
   *       ("readout anode planes to process (omit or empty processes all)")
   *     };
   *   
   * }; // struct Config
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * reading as:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * void readParams(art::EDProducer::Table<Config> const& config) {
   *   
   *   std::vector<geo::ROPID> ROPs
   *     = readOptionalIDsequence(config().ROPs, {});
   *   
   * } // readParams()
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * Note however that the default value will not show in the regular
   * `lar --print-description` output.
   */
  /// @{
  
  // --- BEGIN -- Cryostat ID --------------------------------------------------
  /// Configuration structure for validated `readout::CryostatID` parameter.
  struct CryostatIDConfig: public geo::fhicl::ValidIDConfig {
    using ID_t = readout::CryostatID; ///< Type read by this configuration.
    
    ::fhicl::Atom<readout::CryostatID::CryostatID_t> C {
      ::fhicl::Name("C"),
      ::fhicl::Comment("cryostat number"),
      [this](){ return isValid(); }
      };
    
    ID_t ID() const { return ID_t{ C() }; }
    operator ID_t() const { return ID(); }
    
  }; // struct CryostatIDConfig
  
  /// Member type of validated `readout::CryostatID` parameter.
  using CryostatID = ::fhicl::Table<CryostatIDConfig>;
  
  /// Member type of optional validated `readout::CryostatID` parameter.
  using OptionalCryostatID = ::fhicl::OptionalTable<CryostatIDConfig>;
  
  /// Member type of sequence of `readout::CryostatID` parameters.
  using CryostatIDsequence = ::fhicl::Sequence<CryostatID>;
  
  /// Member type of optional sequence of `readout::CryostatID` parameters.
  using OptionalCryostatIDsequence = ::fhicl::OptionalSequence<CryostatID>;
  
  // --- END -- Cryostat ID ----------------------------------------------------
  
  
  // --- BEGIN -- TPC set ID ---------------------------------------------------
  /// Configuration structure for validated `readout::TPCsetID` parameter.
  struct TPCsetIDConfig: public CryostatIDConfig {
    using ID_t = readout::TPCsetID; ///< Type read by this configuration.
    
    ::fhicl::Atom<readout::TPCsetID::TPCsetID_t> S {
      ::fhicl::Name("S"),
      ::fhicl::Comment("TPC set number within the cryostat"),
      [this](){ return isValid(); }
      };
    
    ID_t ID() const { return { CryostatIDConfig::ID(), S() }; }
    operator ID_t() const { return ID(); }
  }; // struct TPCsetIDIDConfig
  
  /// Member type of validated `readout::TPCsetID` parameter.
  using TPCsetID = ::fhicl::Table<TPCsetIDConfig>;
  
  /// Member type of optional validated `readout::TPCsetID` parameter.
  using OptionalTPCsetID = ::fhicl::OptionalTable<TPCsetIDConfig>;
  
  /// Member type of sequence of `readout::TPCsetID` parameters.
  using TPCsetIDsequence = ::fhicl::Sequence<TPCsetID>;
  
  /// Member type of optional sequence of `readout::TPCsetID` parameters.
  using OptionalTPCsetIDsequence = ::fhicl::OptionalSequence<TPCsetID>;
  
  // --- END -- TPC set ID -----------------------------------------------------
  
  
  // --- BEGIN -- Readout plane ID ---------------------------------------------
  /// Configuration structure for validated `readout::ROPID` parameter.
  struct ROPIDConfig: public TPCsetIDConfig {
    using ID_t = readout::ROPID; ///< Type read by this configuration.
    
    ::fhicl::Atom<readout::ROPID::ROPID_t> R {
      ::fhicl::Name("R"),
      ::fhicl::Comment("Readout plane number within the TPC set"),
      [this](){ return isValid(); }
      };
    
    ID_t ID() const { return { TPCsetIDConfig::ID(), R() }; }
    operator ID_t() const { return ID(); }
  }; // struct ROPIDConfig
  
  /// Member type of validated `readout::ROPID` parameter.
  using ROPID = ::fhicl::Table<ROPIDConfig>;
  
  /// Member type of optional validated `readout::ROPID` parameter.
  using OptionalROPID = ::fhicl::OptionalTable<ROPIDConfig>;
  
  /// Member type of sequence of `readout::ROPID` parameters.
  using ROPIDsequence = ::fhicl::Sequence<ROPID>;
  
  /// Member type of optional sequence of `readout::ROPID` parameters.
  using OptionalROPIDsequence = ::fhicl::OptionalSequence<ROPID>;
  
  // --- END -- Readout plane ID -----------------------------------------------
  
  
  // --- BEGIN -- Channel ID ---------------------------------------------------
  
  /// Member type of validated `raw::ChannelID_t` parameter.
  using ChannelID = ::fhicl::Atom<raw::ChannelID_t>;
  
  /// Member type of optional validated `raw::ChannelID_t` parameter.
  using OptionalChannelID = ::fhicl::OptionalAtom<raw::ChannelID_t>;
  
  /// Member type of sequence of `raw::ChannelID_t` parameters.
  using ChannelIDsequence = ::fhicl::Sequence<raw::ChannelID_t>;
  
  /// Member type of optional sequence of `raw::ChannelID_t` parameters.
  using OptionalChannelIDsequence = ::fhicl::OptionalSequence<raw::ChannelID_t>;
  
  // --- END -- Channel ID -----------------------------------------------------
  
  /// @}
  // --- END -- Validated configuration parameters for readout ID objects ------
  
} // namespace readout::fhicl

#endif // LARCOREOBJ_SIMPLETYPESANDCONSTANTS_READOUT_TYPES_FHICL_H


namespace {
  
  template <typename T, typename Src>
  std::vector<T> convertToVectorOf(Src const& src) {
    std::vector<T> dest;
    dest.reserve(src.size());
    std::copy(util::begin(src), util::end(src), dest.begin());
    return dest;
  } // convertToVectorOf(Src)
  
} // local namespace


///Detector simulation of raw signals on wires
namespace detsim {
    
// Base class for creation of raw signals on wires.
class SimWireICARUS : public art::EDProducer
{
public:
    
    /// Module configuration.
    struct Config {
      
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      
      // --- BEGIN -- Source parameters ----------------------------------------
      /// @name Source parameters
      /// @{
      
      /// Parameters for a single test charge.
      struct TestChargeParams {
        
        fhicl::Atom<std::size_t> Index {
          Name("Index"),
          Comment("TDC count to inject the test charge at")
          };
        fhicl::Atom<double> Charge {
          Name("Charge"),
          Comment("charge to be injected")
          };
        
      }; // struct TestChargeParams
      
      
      geo::fhicl::OptionalWireIDsequence TestWires {
        Name("TestWires"),
        Comment(
          "wire IDs to inject test charge into"
          " (e.g. { C:0 T:1 P:1 W:23 } => C:0 T:1 P:1 W:23)"
          ),
        };
      fhicl::Sequence<fhicl::Table<TestChargeParams>> TestCharges {
        Name("TestCharges"),
        Comment("test charges that are injected into all test wires"),
        fhicl::use_if(this, &Config::isTesting)
        };
      
      fhicl::Atom<art::InputTag> DriftEModuleLabel {
        Name("DriftEModuleLabel"),
        Comment(
          "data product tag for input drifted electrons (`sim::SimChannel`)"
          " [forbidden if test wires are specified]"
          ),
        fhicl::use_if(this, &Config::isNotTesting)
        };
      
      bool isTesting() const {
        // FIXME when issue #23652 is resolved
        // return TestWires.hasValue();
        decltype(TestWires)::value_type dummy; return TestWires(dummy);
        }
      
      
      bool isNotTesting() const { return !isTesting(); }
      
      /// @}
      // --- END -- Source parameters ------------------------------------------
      
      
      // --- BEGIN -- Detector region ------------------------------------------
      /// @name Detector region
      /// @{
      
      geo::fhicl::OptionalTPCIDsequence TPCs {
//       fhicl::Sequence<geo::fhicl::TPCID> TPCs {
        Name("TPCs"),
        Comment("only process channels on these TPC's (empty or omitted processes all)")
//        , std::vector<geo::TPCID>{} // default
        };
      
      /// @}
      // --- END -- Source parameters ------------------------------------------
      
      
      // --- BEGIN -- Output format --------------------------------------------
      /// @name Output format
      /// @{
      
      fhicl::Atom<std::string> CompressionType {
        Name("CompressionType"),
        Comment("waveform output compression type: \"none\" or \"Huffman\""),
        "none" // default
        };
      
      fhicl::Atom<bool> SuppressNoSignal {
        Name("SuppressNoSignal"),
        Comment
          ("skip all channels of the boards with only channels with no charge")
        // default
        };
      
      /// @}
      // --- END -- Output format ----------------------------------------------
      
      
      // --- BEGIN -- Readout information --------------------------------------
      /// @name Readout information
      /// @{
      
      fhicl::Atom<int> NumChanPerMB {
        Name("NumChanPerMB"),
        Comment("channels on the same plane in a TPC readout board"),
        32 // default
        };
      
      /// @}
      // --- END -- Readout information ----------------------------------------
      
      
      // --- BEGIN -- Simulation settings --------------------------------------
      /// @name Simulation settings
      /// @{
      
      fhicl::Atom<bool> SimDeadChannels {
        Name("SimDeadChannels"),
        Comment("simulate also channels identified as bad (otherwise skipped)")
        };
      
      fhicl::DelegatedParameter NoiseGenToolVec {
        Name("NoiseGenToolVec"),
        Comment("configuration of noise generator tools, one per plane")
        };
      
      fhicl::Atom<bool> SmearPedestals {
        Name("SmearPedestals"),
        Comment(
          "apply random fluctuations to channel pedestal levels (from database)"
          ),
        true // default
        };
      
      /// @}
      // --- END -- Simulation settings ----------------------------------------
      
      
      // --- BEGIN -- Random generator seeds -----------------------------------
      /// @name Random generator seeds
      /// @{
      
      rndm::SeedAtom Seed {
        Name("Seed"),
        Comment("random engine seed for coherent noise and uncoherent noise")
        };
      
      rndm::SeedAtom SeedPedestal {
        Name("SeedPedestal"),
        Comment("random engine seed for pedestal slow fluctuations")
        };
      
      /// @}
      // --- END -- Random generator seeds -------------------------------------
      
      
      fhicl::Atom<bool> MakeHistograms {
        Name("MakeHistograms"),
        Comment
          ("also produces a few histograms (stored in TFileService output"),
        false // default
        };
      
    }; // struct Config
    
    using Parameters = art::EDProducer::Table<Config>;
    
    
    struct TestChargeParams {
      std::size_t index;
      double charge;
      
      TestChargeParams() = default;
      TestChargeParams(Config::TestChargeParams const& config)
        : index(config.Index()), charge(config.Charge()) {}
      
    }; // struct TestChargeParams
    
    
    explicit SimWireICARUS(Parameters const& config);
    
    // read/write access to event
    void produce (art::Event& evt);
    void beginJob();
    void endJob();
//    void reconfigure(fhicl::ParameterSet const& p);
    
private:
    
    void MakeADCVec(std::vector<short>& adc, const icarusutil::TimeVec& noise, const icarusutil::TimeVec& charge, float ped_mean) const;
    
    art::InputTag const                                  fDriftEModuleLabel; ///< module making the ionization electrons
    std::optional<std::vector<geo::TPCID>>               fTPCs;              ///< Process only these TPCs
    raw::Compress_t                                      fCompression;       ///< compression type to use
    unsigned int                                         fNTimeSamples;      ///< number of ADC readout samples in all readout frames (per event)
    std::map< double, int >                              fShapingTimeOrder;
    
    bool const                                           fSimDeadChannels;   ///< if True, simulate dead channels using the ChannelStatus service.  If false, do not simulate dead channels
    bool const                                           fSuppressNoSignal;  ///< If no signal on wire (simchannel) then suppress the channel
    bool const                                           fSmearPedestals;    ///< If True then we smear the pedestals
    int const                                            fNumChanPerMB;      ///< Number of channels per motherboard
    
    std::vector<std::unique_ptr<icarus_tool::IGenNoise>> fNoiseToolVec;      ///< Tool for generating noise
    std::unique_ptr<icarusutil::ICARUSFFT<double>>       fFFT;               ///< Object to handle thread safe FFT

    bool const                                           fMakeHistograms;
    std::vector<geo::WireID>                             fTestWires; ///< Where to inject test charge.
    std::vector<TestChargeParams>                        fTestParams; ///< When to inject which test charge.
    std::vector<sim::SimChannel>                         fTestSimChannel_v;
    
    
    ///< Range of channels to process: [ `first`, `second` [
    std::pair<raw::ChannelID_t, raw::ChannelID_t>        fChannelRange;
    
    TH1F*                                                fSimCharge;
    TH2F*                                                fSimChargeWire;
    
    // Random engines
    CLHEP::HepRandomEngine&                              fPedestalEngine;
    CLHEP::HepRandomEngine&                              fUncNoiseEngine;
    CLHEP::HepRandomEngine&                              fCorNoiseEngine;

    //define max ADC value - if one wishes this can
    //be made a fcl parameter but not likely to ever change
    const float                                          adcsaturation = 4095;
    
    // little helper class to hold the params of each charge dep
    class ResponseParams {
    public:
        ResponseParams(double charge, size_t time) : m_charge(charge), m_time(time) {}
        double getCharge() { return m_charge; }
        size_t getTime()   { return m_time; }
    private:
        double m_charge;
        size_t m_time;
    };
    
    //services
    const geo::GeometryCore& fGeometry;
    
    
    bool processAllTPCs() const { return !fTPCs.has_value(); }
    
    bool isTesting() const { return !fTestWires.empty(); }

    /// Returns IDs of first and past-the-last channel to process.
    std::pair<raw::ChannelID_t, raw::ChannelID_t> channelRangeToProcess() const;
}; // class SimWireICARUS

DEFINE_ART_MODULE(SimWireICARUS)

//-------------------------------------------------
SimWireICARUS::SimWireICARUS(Parameters const& config)
    : EDProducer(config)
    , fDriftEModuleLabel(config().DriftEModuleLabel())
    , fTPCs             (geo::fhicl::readOptionalIDsequence(config().TPCs))
    , fSimDeadChannels  (config().SimDeadChannels  ())
    , fSuppressNoSignal (config().SuppressNoSignal ())
    , fSmearPedestals   (config().SmearPedestals   ())
    , fNumChanPerMB     (config().NumChanPerMB     ())
    , fMakeHistograms   (config().MakeHistograms   ())
    , fTestWires        (geo::fhicl::readOptionalIDsequence(config().TestWires, {}))
    , fTestParams       (convertToVectorOf<TestChargeParams>(config().TestCharges()))
    , fPedestalEngine   (art::ServiceHandle<rndm::NuRandomService>()->createEngine
                         (*this, "HepJamesRandom", "pedestal", config().SeedPedestal)
                        )
    , fUncNoiseEngine   (art::ServiceHandle<rndm::NuRandomService>()->createEngine
                         (*this, "HepJamesRandom", "noise",    config().Seed)
                        )
    , fCorNoiseEngine   (art::ServiceHandle<rndm::NuRandomService>()->createEngine
                         (*this, "HepJamesRandom", "cornoise", config().Seed)
                        )
    , fGeometry         (*lar::providerFrom<geo::Geometry>())
{
    
    std::vector<fhicl::ParameterSet> noiseToolParamSetVec
      = config().NoiseGenToolVec.get<std::vector<fhicl::ParameterSet>>();
    
    for(auto& noiseToolParams : noiseToolParamSetVec) {
        fNoiseToolVec.push_back(art::make_tool<icarus_tool::IGenNoise>(noiseToolParams));
    }
    //Map the Shaping Times to the entry position for the noise ADC
    //level in fNoiseFactInd and fNoiseFactColl
    fShapingTimeOrder = { {0.6, 0}, {1, 1}, {1.3, 2}, {3.0, 3} };
    //detector properties information
    auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    
    fNTimeSamples = detprop->NumberTimeSamples();

    TString compression(config().CompressionType());
    if (compression.IsNull() || compression.Contains("none", TString::kIgnoreCase)) fCompression = raw::kNone;
    else if (compression.Contains("Huffman", TString::kIgnoreCase))                 fCompression = raw::kHuffman;
    else throw art::Exception(art::errors::Configuration) << "Unsupported compression requested: '" << compression << "'\n";
  
    
    fChannelRange = channelRangeToProcess();
    if (!processAllTPCs()) 
    {
        mf::LogInfo log("SimWireICARUS");
        log << "Only " << fTPCs->size() << " TPC's will be processed:";
      
        for (geo::TPCID const& tpcid: fTPCs.value()) log << " { " << tpcid << " }";
      
        auto const [ firstChannel, endChannel ] = fChannelRange;
      
        log << "\nAll the " << (endChannel - firstChannel) << " channels from "
            << firstChannel << " to " << endChannel
            << " (excluded) will be processed.";
    } // if selected TPCs
    
    // Get instance of FFT machine
    fFFT = std::make_unique<icarusutil::ICARUSFFT<double>>(fNTimeSamples);
    
    //
    // input:
    //
    if(!isTesting()) consumes<std::vector<sim::SimChannel>>(fDriftEModuleLabel);
    
    //
    // output:
    //
    produces<std::vector<raw::RawDigit>>();
    
    
} // SimWireICARUS::SimWireICARUS()
//-------------------------------------------------
void SimWireICARUS::beginJob()
{
    // get access to the TFile service
    art::ServiceHandle<art::TFileService> tfs;
    
    // If in test mode create a test data set
    if(isTesting())
    {
        std::array<double, 3U> xyz;
        xyz.fill(std::numeric_limits<double>::quiet_NaN());
        for (geo::WireID const& wire: fTestWires) {
            
            raw::ChannelID_t const channel = fGeometry.PlaneWireToChannel(wire);
            
            sim::SimChannel sch(channel);
            for (TestChargeParams const& params: fTestParams) {
                sch.AddIonizationElectrons(
                  -1, params.index, params.charge,
                  xyz.data(), std::numeric_limits<double>::max()
                  );
            } // for inject parameters
            
            fTestSimChannel_v.push_back(std::move(sch));
            
        } // for wires in test
        
    } // if testing
    
    fSimCharge     = tfs->make<TH1F>("fSimCharge", "simulated charge", 150, 0, 1500);
    fSimChargeWire = tfs->make<TH2F>("fSimChargeWire", "simulated charge", 5600,0.,5600.,500, 0, 1500);
    
    return;
}
//-------------------------------------------------
void SimWireICARUS::endJob()
{}
void SimWireICARUS::produce(art::Event& evt)
{
    //--------------------------------------------------------------------
    //
    // Get all of the services we will be using
    //
    //--------------------------------------------------------------------
    
    //get pedestal conditions
    const lariov::DetPedestalProvider& pedestalRetrievalAlg = art::ServiceHandle<lariov::DetPedestalService>()->GetPedestalProvider();
    
    //channel status for simulating dead channels
    const lariov::ChannelStatusProvider& ChannelStatusProvider = art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider();
    
#ifdef ICARUSCODE_SIMWIREICARUS_TRIGGERTIMEHACK
    // In case trigger simulation is run in the same job...
    // FIXME:  You should not be calling preProcessEvent
    art::ServiceHandle<detinfo::DetectorClocksServiceStandard>()
      ->preProcessEvent(evt,art::ScheduleContext::invalid());
#endif // ICARUSCODE_SIMWIREICARUS_TRIGGERTIMEHACK

    auto const* ts = lar::providerFrom<detinfo::DetectorClocksService>();
    
    // get the geometry to be able to figure out signal types and chan -> plane mappings
    const raw::ChannelID_t maxChannel = fGeometry.Nchannels();
    
    //Get N_RESPONSES from SignalShapingService, on the fly
    // flag added to use nominal one response per plane or multiple responses
    // per plane and scaling for YZ dependent responses
    // or data driven field responses
    art::ServiceHandle<icarusutil::SignalShapingICARUSService> sss;

    //--------------------------------------------------------------------
    //
    // Get the SimChannels, which we will use to produce RawDigits
    //
    //--------------------------------------------------------------------
    // make a vector of const sim::SimChannel* that has same number
    // of entries as the number of channels in the detector
    // and set the entries for the channels that have signal on them
    // using the chanHandle
    std::vector<const sim::SimChannel*> channels(maxChannel,nullptr);
    if(!isTesting())
    {
        std::vector<const sim::SimChannel*> chanHandle;
        evt.getView(fDriftEModuleLabel,chanHandle);
        
        for(const auto& simChannel : chanHandle) channels.at(simChannel->Channel()) = simChannel;
    }
    else
        for(const auto& testChannel : fTestSimChannel_v) channels.at(testChannel.Channel()) = &testChannel;
    
    // make a unique_ptr of sim::SimDigits that allows ownership of the produced
    // digits to be transferred to the art::Event after the put statement below
    std::unique_ptr< std::vector<raw::RawDigit>> digcol(new std::vector<raw::RawDigit>);
    digcol->reserve(maxChannel);
    //--------------------------------------------------------------------
    //
    // Loop over channels a second time and produce the RawDigits by adding together
    // pedestal, noise, and direct & induced charges
    //
    //--------------------------------------------------------------------
    
    // vectors for working in the following for loop
    std::vector<short>  adcvec(fNTimeSamples, 0);
    icarusutil::TimeVec chargeWork(fNTimeSamples,0.);
    icarusutil::TimeVec zeroCharge(fNTimeSamples,0.);
    icarusutil::TimeVec noisetmp(fNTimeSamples,0.);
    
    // make sure chargeWork is correct size
    if (chargeWork.size() < fNTimeSamples) throw std::range_error("SimWireICARUS: chargeWork vector too small");
    
    //detector properties information
    auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    
    // Let the tools know to update to the next event
    for(const auto& noiseTool : fNoiseToolVec) noiseTool->nextEvent();

    // The original implementation would allow the option to skip channels for which there was no MC signal
    // present. We want to update this so that if there is an MC signal on any wire in a common group (a
    // motherboard) then we keep all of those wires. This so we can implment noise mitigation techniques
    // with the simulation
    //
    // So... first step is to build a map of motherboard and true information
    using MBWithSignalSet = std::set<raw::ChannelID_t>;
    
    MBWithSignalSet mbWithSignalSet;
    
    auto const [ firstChannel, endChannel ] = fChannelRange;
    
    // If we are not suppressing the signal then we need to make sure there is an entry in the set for every motherboard
    if (!fSuppressNoSignal)
    {
        raw::ChannelID_t firstMBIdx(firstChannel / fNumChanPerMB);
        raw::ChannelID_t endMBIdx(endChannel / fNumChanPerMB);
        
        for(raw::ChannelID_t mbIdx = firstMBIdx; mbIdx < endMBIdx; mbIdx++) mbWithSignalSet.insert(mbIdx);
    }
    else
    {
        for(const auto& simChan : channels)
        {
            if (simChan)
            {
                raw::ChannelID_t channel = simChan->Channel();
                
                if (channel >= firstChannel && channel < endChannel) mbWithSignalSet.insert(channel/fNumChanPerMB);
            }
        }
    }
    
    // Ok, now we can simply loop over MB's...
    for(const auto& mb : mbWithSignalSet)
    {
        raw::ChannelID_t baseChannel = fNumChanPerMB * mb;
        
        // And for a given MB we can loop over the channels it contains
        for(raw::ChannelID_t channel = baseChannel; channel < baseChannel + fNumChanPerMB; channel++)
        {
            //clean up working vectors from previous iteration of loop
            adcvec.resize(fNTimeSamples, 0);  //compression may have changed the size of this vector
            noisetmp.resize(fNTimeSamples, 0.);     //just in case
            
            //use channel number to set some useful numbers
            std::vector<geo::WireID> widVec = fGeometry.ChannelToWire(channel);
            size_t                   plane  = widVec[0].Plane;
            
            //Get pedestal with random gaussian variation
            float ped_mean = pedestalRetrievalAlg.PedMean(channel);
            
            if (fSmearPedestals )
            {
                CLHEP::RandGaussQ rGaussPed(fPedestalEngine, 0.0, pedestalRetrievalAlg.PedRms(channel));
                ped_mean += rGaussPed.fire();
            }
            
            //Generate Noise
            double noise_factor(0.);
            auto   tempNoiseVec = sss->GetNoiseFactVec();
            double shapingTime  = sss->GetShapingTime(0);
            
            if (fShapingTimeOrder.find( shapingTime ) != fShapingTimeOrder.end() )
                noise_factor = tempNoiseVec[plane].at( fShapingTimeOrder.find( shapingTime )->second );
            //Throw exception...
            else
            {
                throw cet::exception("SimWireICARUS")
                << "\033[93m"
                << "Shaping Time received from signalservices_microboone.fcl is not one of allowed values"
                << std::endl
                << "Allowed values: 0.6, 1.0, 1.3, 3.0 usec"
                << "\033[00m"
                << std::endl;
            }
            
            // Use the desired noise tool to actually generate the noise on this wire
            fNoiseToolVec[plane]->generateNoise(fUncNoiseEngine,
                                                fCorNoiseEngine,
                                                noisetmp,
                                                noise_factor,
                                                channel);
            
            // Recover the SimChannel (if one) for this channel
            const sim::SimChannel* simChan = channels[channel];
            
            // If there is something on this wire, and it is not dead, then add the signal to the wire
            if(simChan && !(fSimDeadChannels && (ChannelStatusProvider.IsBad(channel) || !ChannelStatusProvider.IsPresent(channel))))
            {
                double gain=sss->GetASICGain(channel) * detprop->SamplingRate() * 1.e-3; // Gain returned is electrons/us, this converts to electrons/tick
                
                std::fill(chargeWork.begin(), chargeWork.end(), 0.);
                
                // loop over the tdcs and grab the number of electrons for each
                for(int tick = 0; tick < int(fNTimeSamples); tick++)
                {
                    int tdc = ts->TPCTick2TDC(tick);
                    
                    // continue if tdc < 0
                    if( tdc < 0 ) continue;
                    
                    double charge = simChan->Charge(tdc);  // Charge returned in number of electrons
                    
                    chargeWork[tick] += charge/gain;  // # electrons / (# electrons/tick)
                } // loop over tdcs
                // now we have the tempWork for the adjacent wire of interest
                // convolve it with the appropriate response function
                fFFT->convolute(chargeWork, sss->GetResponse(channel).getConvKernel(), sss->FieldResponseTOffset(channel));
                
                // "Make" the ADC vector
                MakeADCVec(adcvec, noisetmp, chargeWork, ped_mean);
            }
            // "Make" an ADC vector with zero charge added
            else MakeADCVec(adcvec, noisetmp, zeroCharge, ped_mean);
            
            // add this digit to the collection;
            // adcvec is copied, not moved: in case of compression, adcvec will show
            // less data: e.g. if the uncompressed adcvec has 9600 items, after
            // compression it will have maybe 5000, but the memory of the other 4600
            // is still there, although unused; a copy of adcvec will instead have
            // only 5000 items. All 9600 items of adcvec will be recovered for free
            // and used on the next loop.
            raw::RawDigit rd(channel, fNTimeSamples, adcvec, fCompression);
            
            if(fMakeHistograms && plane==2)
            {
                short area = std::accumulate(adcvec.begin(),adcvec.end(),0,[](const auto& val,const auto& sum){return sum + val - 400;});
                
                if(area>0)
                {
                    fSimCharge->Fill(area);
                    fSimChargeWire->Fill(widVec[0].Wire,area);
                }
            }
            
            rd.SetPedestal(ped_mean);
            digcol->push_back(std::move(rd)); // we do move the raw digit copy, though

        }
    }
    
    evt.put(std::move(digcol));
    
    return;
}
//-------------------------------------------------
void SimWireICARUS::MakeADCVec(std::vector<short>& adcvec, const icarusutil::TimeVec& noisevec,
                               const icarusutil::TimeVec& chargevec, float ped_mean) const
{
    for(unsigned int i = 0; i < fNTimeSamples; ++i)
    {
        float adcval = noisevec[i] + chargevec[i] + ped_mean;

        adcval = std::max(float(0.), std::min(adcval, adcsaturation));

        adcvec[i] = std::round(adcval);
    }// end loop over signal size
    // compress the adc vector using the desired compression scheme,
    // if raw::kNone is selected nothing happens to adcvec
    // This shrinks adcvec, if fCompression is not kNone.
    raw::Compress(adcvec, fCompression);
    
    return;
}
//-------------------------------------------------
std::pair<raw::ChannelID_t, raw::ChannelID_t>
SimWireICARUS::channelRangeToProcess() const {
    
    // return the first and last channel numbers
    // based on whether we are outputting selected TPC's or all of them
    
    raw::ChannelID_t const maxChannel { fGeometry.Nchannels() };
    
    if (processAllTPCs())
        return { raw::ChannelID_t{ 0 }, maxChannel };
    
    //
    // channel selection
    //

    lar::util::MinMaxCollector<raw::ChannelID_t> stats;
    
    for (geo::TPCID const& tpcid: fTPCs.value()) {
        
        for (geo::PlaneGeo const& plane: fGeometry.IteratePlanes(tpcid)) {
            
            raw::ChannelID_t const planeStartChannel
              = fGeometry.PlaneWireToChannel({ plane.ID(), 0U });
            
            stats.add(planeStartChannel);
            
            raw::ChannelID_t const planeEndChannel
              = fGeometry.PlaneWireToChannel({ plane.ID(), plane.Nwires() - 1U }) + 1;
            
            stats.add(planeEndChannel);
            
        } // for planes in TPC
        
    } // for all TPCs
    
    assert(stats.has_data());
    
    return { stats.min(), stats.max() };
    
} // SimWireICARUS::channelRangeToProcess()


//-------------------------------------------------
    
}
