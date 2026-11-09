/**
 * @file   DataProductReader_module.cc
 * @brief  An _art_ producer reading data products and telling which ones.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   June 3, 2026
 *
 */

// LArSoft and ICARUS libraries
#include "icarusalg/Utilities/mfLoggingClass.h"
#include "larcorealg/CoreUtils/DebugUtils.h" // lar::debug::demangle()

// data product headers
#include "sbnobj/Common/CRT/CRTHit.hh"
#include "sbnobj/Common/CRT/CRTPMTMatching.hh"
#include "sbnobj/Common/Trigger/ExtraTriggerInfo.h"
#include "sbnobj/ICARUS/CRT/CRTData.hh"
#include "sbnobj/ICARUS/PMT/Data/WaveformBaseline.h"
#include "sbnobj/ICARUS/PMT/Trigger/Data/OpticalTriggerGate.h"
#include "lardataobj/RawData/TriggerData.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/Simulation/BeamGateInfo.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/Simulation/SimEnergyDepositLite.h"
#include "lardataobj/Simulation/SimPhotons.h"

// framework libraries
#include "art/Framework/Core/ConsumesCollector.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/SharedProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Sequence.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// C/C++ standard libraries
#include <array>
#include <cassert>
#include <memory> // std::unique_ptr
#include <string>
#include <tuple>
#include <utility> // std::move(), std::index_sequence, ...
#include <vector>


//------------------------------------------------------------------------------
// === BEGIN --- Add the data product types to be supported here  ==============
// 
// (and add the appropriate header above and library in CMakeLists.txt`)
//
using SupportedTypes = std::tuple<
  // --- physics ---------------------------------------------------------------
  std::vector<sim::SimEnergyDeposit>,
  std::vector<sim::SimEnergyDepositLite>,
  // --- trigger ---------------------------------------------------------------
  std::vector<raw::Trigger>,
  std::vector<sim::BeamGateInfo>,
  sbn::ExtraTriggerInfo,
    // key icarus_trigger_ReadoutTriggerGate_long_long_unsigned_ints:
  std::vector<icarus::trigger::OpticalTriggerGate::GateData_t>,
  // --- optical ---------------------------------------------------------------
  std::vector<sim::SimPhotons>,
  std::vector<recob::OpHit>,
  std::vector<recob::OpFlash>,
  std::vector<icarus::WaveformBaseline>,
  std::vector<raw::OpDetWaveform>,
  art::Assns<raw::OpDetWaveform, icarus::WaveformBaseline>,
  // --- CRT -------------------------------------------------------------------
  std::vector<sim::AuxDetSimChannel>,
  std::vector<icarus::crt::CRTData>,
  std::vector<sbn::crt::CRTHit>,
  art::Assns<recob::OpFlash, recob::OpHit>,
  // --- matching --------------------------------------------------------------
  std::vector<sbn::crt::CRTPMTMatching>
  // ---------------------------------------------------------------------------
  >;
// === END ----- Add the data product types to be supported here  ==============
//------------------------------------------------------------------------------

// forward declarations
class DataProductHandler;
class DataProductReader;

//------------------------------------------------------------------------------
/**
 * @brief Attempts to retrieve data products by tag and tells which ones.
 *
 * For each specified (supported) data product type and tag, this module prints
 * on screen which data product would be read, reporting the full tag (included
 * instance and process name). If a product is not available (including when
 * it is known but has been dropped), the module will report it as not found.
 * 
 * This module currently does not distinguish between products which were
 * dropped and products that never existed.
 * 
 *
 * Configuration parameters
 * -------------------------
 * 
 * The supported parameters are printed by executing:
 *     
 *     lar --print-description DataProductReader
 *     
 * For each supported type of data product, a sequence is generated with the
 * name of the data product as key. Two examples:
 * 
 * * `sbn_ExtraTriggerInfo` (list of input tags, default: empty): print the
 *     available among the specified tags for data products of type
 *     `sbn::ExtraTriggerInfo`.
 * * `recob::OpFlashs` (list of input tags, default: empty): print the
 *     available among the specified tags for data products of type
 *     `std::vector<recob::OpFlash>`.
 *
 * For each data product type, all tags configured for that type are attempted;
 * if _one of the tags_ in the list is empty (`""`), _all_ data products of that
 * type are listed. If the _list of tags_ is empty, no data product of that type
 * is tested.
 * 
 * The key of the configuration parameter is the full name of its data product,
 * including all template parameters if any. As a special rule, a `std::vector`
 * of data elements will be configured with the name of the type of the
 * contained data, with a letter "s" appended to mock a plural (with the grammar
 * horrific consequences illustrated in the example above for `recob::OpFlash`).
 * 
 * 
 * Output
 * -------
 * 
 * This "producer" prints information on console but it does not produce any
 * data.
 * It is a producer solely to allow it being put at any point of the workflow
 * (`trigger_paths`).
 * 
 * 
 * Implementation notes
 * ---------------------
 * 
 * Due to the static nature of _art_ data model, the type of data product needs
 * to be hard-coded into the module. Therefore, there are only a fixed number
 * of data product types that are supported.
 * 
 * To add support for a data product type, add it in the definition of
 * `SupportedTypes` (the order does not matter).
 * 
 */
class DataProductReader
  : public art::SharedProducer
  , private icarus::ns::util::mfLoggingClass
{
  
    public:
  
  using SupportedTypes_t = ::SupportedTypes;
  static constexpr std::size_t NProductTypes = std::tuple_size_v<SupportedTypes_t>;
  
  struct Config {
    
    using Tags_t = std::array<fhicl::Sequence<art::InputTag>, NProductTypes>;
    
    // All fhicl::Sequence objects are defined here:
    Tags_t tags = tagConfigFor(static_cast<SupportedTypes_t*>(nullptr));

      private:
    template <typename... ProdTypes>
    static Tags_t tagConfigFor(std::tuple<ProdTypes...> const*);
    
  }; // Config
  
  using Parameters = art::SharedProducer::Table<Config>;
  
  
  /// Constructor. No surprise here.
  DataProductReader(Parameters const& params, const art::ProcessingFrame&);
  
  
  virtual void produce(art::Event& event, const art::ProcessingFrame&) override;
  
  
    private:
  using HandlerRegistry_t = std::vector<std::unique_ptr<DataProductHandler>>;
  
  /// All used data product handlers.
  HandlerRegistry_t const fHandlers;
  
  /// Registers all supported data handlers using `allTags` configuration.
  static HandlerRegistry_t makeHandlerRegistry(Config::Tags_t const& allTags);

  /// Implementation detail.
  template <std::size_t... Indices>
  static HandlerRegistry_t makeHandlerRegistryImpl
    (std::index_sequence<Indices...>, Config::Tags_t const& allTags);

  /// Creates a handler for the data product `SupportedTypes_t[Index]`.
  template <typename ProdType>
  static void registerHandler
    (HandlerRegistry_t& registry, std::vector<art::InputTag> const& tags);
  
}; // class DataProductReader


//------------------------------------------------------------------------------
//---  Implementation
//------------------------------------------------------------------------------
namespace {
  
  /// Returns a copy of `s` with all occurrences of `old` replaced by `with`.
  std::string replaceAll
    (std::string s, std::string const& old, std::string const& with)
  {
    
    if (old.empty()) return s;
    
    std::size_t start = 0;
    while (true) {
      
      std::size_t const nextStart = s.find(old, start);
      if (nextStart == std::string::npos) break;
      
      s.replace(nextStart, old.length(), with);
      start = nextStart + with.length();
      assert(s.length() >= start);
    } // while
    return s;
    
  } // replaceAll()
  
  
  /// Human-friendly name for type `T`.
  template <typename T = void, typename = void>
  struct FriendlyClassName {
    template <typename U = T>
    static std::string name() { return lar::debug::demangle<U>(); }
  };
  
  /// Returns a human-friendly name for type `T`.
  template <typename T>
  std::string friendlyClassName() { return FriendlyClassName<T>::name(); }
  
  /// Functor to sanitize a class name into a parameter key
  template <typename T = void>
  struct ClassToKey {
    static std::string makeKey(std::string const& className)
      { return replaceAll(replaceAll(className, "::", "_"), " ", "_"); }
    template <typename U = T>
    static std::string convert()
      { return makeKey(lar::debug::demangle<U>()); }
  };
  
  template <typename T>
  std::string classToKey() { return ClassToKey<T>::convert(); }
  
  /// `FriendlyClassName` for `std::vector`: omit default allocator.
  template <typename T>
  struct FriendlyClassName<std::vector<T>> {
    static std::string name()
      { return "std::vector<" + FriendlyClassName<T>::name() + ">"; }
  };
  
  /// `FriendlyClassName` for `art::Assns`: omit no metadata.
  template <typename L, typename R>
  struct FriendlyClassName<art::Assns<L, R>> {
    static std::string name()
      {
        return "art::Assns<"
          + FriendlyClassName<L>::name() + ", " + FriendlyClassName<R>::name()
          + ">"; 
      }
  };
  
  /// `FriendlyClassName` composition for templates (only type name parameters).
  template <template<typename...> typename Templ, typename... Args>
  struct FriendlyClassName<Templ<Args...>> {
    template <typename First, typename... Others>
    static std::string argumentNames()
      {
        return (
          friendlyClassName<First>()
          + ... + (", " + friendlyClassName<Others>())
        );
      }
    
    static std::string name()
      {
        std::string key = FriendlyClassName<>::name<Templ<Args...>>();
        key.erase(key.find('<'));
        if constexpr(sizeof...(Args) > 0)
          return key + "<" + argumentNames<Args...>() + ">";
        else return key + "<>";
      }
  };
  
  /// Specialization of `ClassToKey` for `std::vector`
  template <typename T, typename... Args>
  struct ClassToKey<std::vector<T, Args...>>: ClassToKey<T> {
    static std::string convert() { return classToKey<T>() + "s"; }
  };
  
  /// Specialization of `ClassToKey` for `art::Assns`: omit no metadata.
  template <typename L, typename R>
  struct ClassToKey<art::Assns<L, R>> {
    static std::string convert()
      { return "art_Assns_" + classToKey<L>() + "_" + classToKey<R>(); }
  };
  
  /// Specialization of `ClassToKey` for templates (only type name parameters).
  template <template<typename...> typename Templ, typename... Args>
  struct ClassToKey<Templ<Args...>> {
    static std::string convert()
      {
        std::string key = ClassToKey<>::convert<Templ<Args...>>();
        key.erase(key.find('<'));
        return (key + ... + ("_" + classToKey<Args>()) );
      }
  };
  
} // local namespace


//------------------------------------------------------------------------------
//---  DataProductHandler
//------------------------------------------------------------------------------
/// Interface of a data product to the module.
class DataProductHandler: protected icarus::ns::util::mfLoggingClass {
  
    protected:
  
  std::vector<art::InputTag> fTags; ///< Tags to be loaded.
  bool fLoadAll = false; ///< Whether to load all data products.
  
    public:
  
  /// Default constructor: no tag registered.
  DataProductHandler(): icarus::ns::util::mfLoggingClass{ "DataProductReader" }
    {}
  
  /// Registers all the requested input tags.
  DataProductHandler(std::vector<art::InputTag> const& tags);
  
  
  /// Returns the class name of the data product.
  virtual std::string dataProductName() const = 0;
  
  /// Processes the consumable.
  virtual void consumes(art::ConsumesCollector& collector) const = 0;
  
  /// Reads the data products and reports.
  virtual void readProducts(art::Event const& event) const = 0;
  
  
  // configuration query
  bool readAllProducts() const { return fLoadAll; }
  std::vector<art::InputTag> const& tags() const { return fTags; }
  std::size_t nTags() const { return fTags.size(); }
  bool empty() const { return fTags.empty() && !readAllProducts(); }
  
}; // DataProductHandler


DataProductHandler::DataProductHandler(std::vector<art::InputTag> const& tags)
  : DataProductHandler()
{
  
  /// Parse which tags to load
  for (art::InputTag const& tag: tags) {
    if (tag.empty()) fLoadAll = true;
    else fTags.push_back(tag);
  }
  
} // DataProductHandler::DataProductHandler()


//------------------------------------------------------------------------------
//---  DataProductHandlerFor<>
//------------------------------------------------------------------------------
template <typename ProdType>
class DataProductHandlerFor: public DataProductHandler {
  
    public:
  
  using Product_t = ProdType; ///< Type being handled.
  using Handle_t = art::Handle<Product_t>; ///< Type of handle to this product.
  
  using DataProductHandler::DataProductHandler; // inherit constructors
  
  
  /// Returns the class name of the data product.
  virtual std::string dataProductName() const override
    { return staticDataProductName(); }
  
  /// Processes the consumable.
  virtual void consumes(art::ConsumesCollector& collector) const override;
  
  /// Reads the data products and reports.
  virtual void readProducts(art::Event const& event) const override;
  
  /// fhicl::Sequence type set up to read tags for this type.
  struct TagConfig: fhicl::Sequence<art::InputTag> { TagConfig(std::string name); };
  
  
  /// Returns the name of the data product type being handle (static function);
  static std::string staticDataProductName()
    { return friendlyClassName<Product_t>(); }
  
  /// Prints into `out` information about `handle` ("NOT FOUND" if invalid).
  template <typename Stream>
  static Stream& presentHandle(Stream& out, Handle_t const& handle);

  /// Prints into `out` information about all `handles`.
  template <typename Stream>
  static Stream& presentAllHandles
    (Stream& out, std::vector<Handle_t> const& handles);

}; // DataProductHandlerFor


//------------------------------------------------------------------------------
template <typename ProdType>
class DataProductHandlerForVectorOf
  : public DataProductHandlerFor<std::vector<ProdType>>
{
  using Base_t = DataProductHandlerFor<std::vector<ProdType>>;
  
    public:
  
  using Base_t::Base_t;
  
}; // DataProductHandlerForVectorOf


//------------------------------------------------------------------------------
template <typename ProdType>
void DataProductHandlerFor<ProdType>::consumes
  (art::ConsumesCollector& collector) const
{
  if (readAllProducts()) collector.consumesMany<Product_t>();
  for (art::InputTag const& tag: tags()) collector.consumes<Product_t>(tag);
}


//------------------------------------------------------------------------------
template <typename ProdType>
void DataProductHandlerFor<ProdType>::readProducts
  (art::Event const& event) const
{
  auto log = mfLogVerbatim();
  log << "Data products of type " << dataProductName() << ":";
  
  if (fLoadAll) {
    log << "\n";
    presentAllHandles(log, event.getMany<Product_t>());
  }
  
  for (art::InputTag const& tag: fTags) {
    log << "\n - '" << tag.encode() << "' -> ";
    presentHandle(log, event.getHandle<Product_t>(tag));
  }
  
} // DataProductHandlerFor<>::readProducts()


//------------------------------------------------------------------------------
template <typename ProdType>
template <typename Stream>
Stream& DataProductHandlerFor<ProdType>::presentHandle
  (Stream& out, Handle_t const& handle)
{
  if (handle) {
    assert(handle.provenance());
    out << "'" << handle.provenance()->inputTag().encode() << "'";
  }
  else out << "NOT FOUND";
  return out;
} // DataProductHandlerFor<>::presentHandle()


//------------------------------------------------------------------------------
template <typename ProdType>
template <typename Stream>
Stream& DataProductHandlerFor<ProdType>::presentAllHandles
  (Stream& out, std::vector<Handle_t> const& handles)
{
  if (handles.empty()) {
    out << " - no data product present";
    return out;
  }
  
  out << " - " << handles.size() << " in total:";
  for (Handle_t const& handle: handles) {
    out << "\n    * ";
    presentHandle(out, handle);
  }
  return out;
  
} // DataProductHandlerFor<>::presentAllHandles()


//------------------------------------------------------------------------------
//---  DataProductReader
//------------------------------------------------------------------------------
template <typename... ProdTypes>
auto DataProductReader::Config::tagConfigFor(std::tuple<ProdTypes...> const*)
  -> Tags_t
{
  return {
    fhicl::Sequence<art::InputTag>{
      fhicl::Name{ classToKey<ProdTypes>() },
      fhicl::Comment{ "tags for data products: " + friendlyClassName<ProdTypes>() },
      std::vector<art::InputTag>{}
    }...
  };
} // DataProductReader::Config::tagConfigFor()


//------------------------------------------------------------------------------
DataProductReader::DataProductReader
  (Parameters const& params, art::ProcessingFrame const&)
  : art::SharedProducer{ params }
  , icarus::ns::util::mfLoggingClass{ "DataProductReader" }
  , fHandlers{ makeHandlerRegistry(params().tags) }
{
  
  async<art::InEvent>();
  
  {
    auto log = mfLogInfo();
    log << "Registered tags for " << fHandlers.size() << " data product types:";
    for (std::unique_ptr<DataProductHandler> const& handler: fHandlers) {
      log << "\n - " << handler->dataProductName();
      if (handler->empty()) continue;
      log << " (";
      if (!handler->tags().empty()) {
        if (handler->readAllProducts()) log << "all available, plus other ";
        log << handler->nTags() << " tags";
      }
      else if (handler->readAllProducts()) log << "all available tags";
      log << ")";
    } // for
  } // local scope
  
} // DataProductReader::DataProductReader()


//------------------------------------------------------------------------------
void DataProductReader::produce(art::Event& event, const art::ProcessingFrame&)
{
  mfLogInfo() << "Reading " << fHandlers.size() << " data product types:";
  for (std::unique_ptr<DataProductHandler> const& handler: fHandlers) {
    
    handler->readProducts(event);
    
  }
  
} // DataProductReader::produce()


//------------------------------------------------------------------------------
auto DataProductReader::makeHandlerRegistry(Config::Tags_t const& allTags)
  -> HandlerRegistry_t
{
  return makeHandlerRegistryImpl
    (std::make_index_sequence<std::tuple_size_v<SupportedTypes_t>>{}, allTags);
}


//------------------------------------------------------------------------------
template <std::size_t... Indices>
auto DataProductReader::makeHandlerRegistryImpl
  (std::index_sequence<Indices...>, Config::Tags_t const& allTags)
  -> HandlerRegistry_t
{
  HandlerRegistry_t registry;
  (
    registerHandler<std::tuple_element_t<Indices, SupportedTypes_t>>
      (registry, allTags[Indices]()),
    ...
  );
  return registry;
}


//------------------------------------------------------------------------------
template <typename ProdType>
void DataProductReader::registerHandler
  (HandlerRegistry_t& registry, std::vector<art::InputTag> const& tags)
{
  if (tags.empty()) return;
  
  registry.push_back(std::make_unique<DataProductHandlerFor<ProdType>>(tags));
  
} // DataProductReader::registerHandler()


//------------------------------------------------------------------------------
DEFINE_ART_MODULE(DataProductReader)


//------------------------------------------------------------------------------

