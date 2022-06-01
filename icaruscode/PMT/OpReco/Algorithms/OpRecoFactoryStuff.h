/**
 * @file   icaruscode/PMT/OpReco/Algorithms/OpRecoFactoryStuff.h
 * @brief  Utility and boilerplate for optical algorithms.
 * @date   May 7, 2022
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * 
 * Some utility and boilerplate to make the use of optical reconstruction
 * algorithms from `larana` more flexible.
 * 
 * This does not attempt to interface those algorithms with _art_ tools, which
 * is definitely possible but just not done here.
 */

#ifndef ICARUSCODE_PMT_OPRECO_ALGORITMS_OPRECOFACTORYSTUFF_H
#define ICARUSCODE_PMT_OPRECO_ALGORITMS_OPRECOFACTORYSTUFF_H

// framework libraries
#include "fhiclcpp/ParameterSet.h"
#include "cetlib_except/exception.h"

// C++ standard libraries
#include <vector>
#include <iterator> // std::empty(), std::begin(), ...
#include <memory> // std::unique_ptr
#include <string>
#include <utility> // std::move()
#include <type_traits> // std::is_base_of_v


// -----------------------------------------------------------------------------
// forward declarations
namespace art { class Event; }


// -----------------------------------------------------------------------------
namespace opdet::factory {
  
  namespace details { struct NoModule_t { explicit NoModule_t() = default; }; }
  
  // ---------------------------------------------------------------------------
  /**
   * @brief Returns a concatenation of strings in `s` separated by `sep`.
   * @tparam S type of string
   * @tparam Coll type of collection of strings
   * @param sep separator string
   * @param s list of strings
   * @return a string of type `S` with the concatenation of all elements of `s`
   * 
   * Example:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * using std::string_literals;
   * std::array const s { "one"s, "2"s, "III"s, "d"s };
   * 
   * std::cout << "Options: '" << join("', '", s) << "'." << std::endl;
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * will print: `Options: 'one', '2', 'III', 'd'.`.
   * 
   * Note that actually `S` can be any copyable type supporting in-place
   * addition (that is concatenation for strings). Indeed, elements in `Coll`
   * do not even need to be of type `S` as long as
   * `S::operator+=(Coll::value_type)` is supported.
   */
  template <typename S, typename Coll>
  S join(S const& sep, Coll const& s);
  
  // some tricks to avoid `S` above to become char* or char[]
  template <typename Coll>
  std::string join(const char* sep, Coll const& s);
  
  // ---------------------------------------------------------------------------
  /// Token to register an algorithm, used in `AlgorithmFactory`.
  template <typename Algo>
  struct Decl {
    using RealAlgo_t = Algo; ///< Type of the algorithm.
    std::string name; ///< Name of the algorithm.
  }; // Decl
  
  template <typename Base>
  class AlgorithmFactory;
  
  // ---------------------------------------------------------------------------
  /**
   * @brief Template type for the third parameter of FWInterfacedIF objects.
   * @tparam Event type of event-level data source
   * @tparam Module type of framework actor calling the algorithm
   *   (e.g. `art::EDProducer` or even `art::ConsumesCollector`)
   */
  template <typename Event = art::Event, typename Module = details::NoModule_t>
  struct FWInterfaceTraits {
    using Event_t = Event;
    using Module_t = Module;
  };
  
  template <typename Base, typename FWTraits = FWInterfaceTraits<>>
  class FWInterfacedIF;
  
  template
    <typename Algo, typename Base, typename FWTraits = FWInterfaceTraits<>>
  class FWInterfacedBase;
  
  template
    <typename Algo, typename Base, typename FWTraits = FWInterfaceTraits<>>
  class FWInterfaced;
  
  
  // ---------------------------------------------------------------------------
  
} // namespace opdet::factory


// ---------------------------------------------------------------------------
/**
 * @brief An algorithm factory class.
 * @tparam Base the base algorithm class
 * 
 * This class merges the utilities above into a coherent interface, minimizing
 * the maintenance cost at the expense of flexibility.
 * 
 * Example of usage:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * opdet::factory::AlgorithmFactory<pmtana::PMTPulseRecoBase> const
 * HitAlgoFactory {
 *     "Name"
 *   , opdet::factory::Decl<pmtana::AlgoThreshold    >{"Threshold"    }
 *   , opdet::factory::Decl<pmtana::AlgoSiPM         >{"SiPM"         }
 *   , opdet::factory::Decl<pmtana::AlgoSlidingWindow>{"SlidingWindow"}
 *   , opdet::factory::Decl<pmtana::AlgoFixedWindow  >{"FixedWindow"  }
 *   , opdet::factory::Decl<pmtana::AlgoCFD          >{"CFD"          }
 *   };
 * 
 * // ...
 * 
 * std::unique_ptr<pmtana::PMTPulseRecoBase> algo = HitAlgoFactory.create(pset);
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * where we assume `pset` is a configuration for the algorithm, and that it
 * contains the name of the algorithm (one of those declared strings) in the key
 * `Name`.
 * 
 */
template <typename Base>
class opdet::factory::AlgorithmFactory {
    public:
  using Algo_t = Base; ///< The algorithm interface this factory produces.
  using Factory_t = AlgorithmFactory<Algo_t>; ///< This type.
  
  
  // --- BEGIN -- Registration of algorithms ---------------------------------
  /// @name Registration of algorithms
  /// @{
  
  /**
   * @brief Derivative of `opdet::factory::Decl` with added type check.
   * @tparam Algo the declared algorithm class
   *
   * This object is used in `AlgorithmFactory` just like
   * `opdet::factory::Decl`, from which it derives. This derivative, whose
   * full identifier is `opdet::factory::AlgorithmFactory<Base>::Decl<Algo>`,
   * is aware of the expected base class of the algorithm and performs a
   * static check on it for added QA.
   */
  template <typename Algo>
  struct Decl: opdet::factory::Decl<Algo> {
    static_assert(std::is_base_of_v<Algo_t, Algo>,
      "The algorithm <Algo> must be a derivate of <Base>.");
  }; // Decl
  
  
  /**
   * @brief Constructor: register the specified algorithms.
   * @tparam Algos the types of algorithms to register
   * @param algorithms the declarators of the algorithms to be registered
   * @see `declare()`
   * 
   * This constructor declares the specified algorithms like with `declare()`.
   * In addition, it may be used for example to create a "constant" instance
   * of a factory.
   * The example in `declare()` becomes:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * using PedAlgoFactory_t
   *   = opdet::factory::AlgorithmFactory<pmtana::PMTPedestalBase>;
   * 
   * PedAlgoFactory_t PedAlgoFactory;
   * PedAlgoFactory.declare(
   *   opdet::factory::Decl<pmtana::PedAlgoEdges>{ "Edges" },
   *   opdet::factory::Decl<pmtana::PedAlgoUB   >{ "UB"    }
   *   );
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * 
   */
  template <typename... Algos>
  AlgorithmFactory(opdet::factory::Decl<Algos>... algorithms)
    : AlgorithmFactory{ "", std::move(algorithms)... } {}

  
  /**
   * @brief Constructor: register the specified algorithms.
   * @tparam Algos the types of algorithms to register
   * @param algoKey name of the configuration key whith the algorithm name
   * @param algorithms the declarators of the algorithms to be registered
   * @see `declare()`, `setAlgorithmConfigurationKey()`
   * 
   * In addition to the declaration of algorithms (see `declare()` or the
   * other constructor), this constructor also sets the configuration key
   * where the algorithm name is expected to be found (@see
   * `setAlgorithmConfigurationKey()`).
   */
  template <typename... Algos>
  AlgorithmFactory
    (std::string algoKey, opdet::factory::Decl<Algos>... algorithms);
  
  
  /**
   * @brief Sets the name of the configuration key with the algorithm name.
   * @param key name of the configuration key with the algorithm name
   * @return this factory
   * @see `create()`
   * 
   * When an algorithm instance is created with `create()`, it is always
   * possible to specify the name of the algorithm explicitly,
   * using `create(std::string const&, fhicl::ParameterSet const&)`.
   * In addition, if it is known that the algorithm name is in a configuration
   * key of that parameter set, the path of that configuration key can be
   * specified here, and the other version of `create()`,
   * `create(fhicl::ParameterSet const&)`, can be used which will discover
   * the required algorithm name from the configuration in argument (which is
   * still the same used to construct the algorithm itself).
   * 
   * Example:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * factory.setAlgorithmConfigurationKey("Name");
   * auto algo = factory.create(pset);
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * is equivalent to:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * auto algo = factory.create(pset.get<std::string>("Name"), pset);
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */
  Factory_t& setAlgorithmConfigurationKey(std::string key)
    { fAlgoNameKey = std::move(key); return *this; }
  
  
  /// Registers an algorithm of type `Algo` associated with a `name`.
  template <typename Algo>
  Factory_t& declare(std::string name);
  
  /**
   * @brief Register a sequence of algorithms
   * @tparam FirstAlgo the type of the first algorithm to register (mandatory)
   * @tparam OtherAlgos the types of more algorithms to register
   * @param first the declarator of the first algorithm to be registered
   * @param others the declarators of more algorithms to be registered
   * @return this factory object
   * 
   * This method allows the declaration of many algorithms at once, one for
   * each declarator in the arguments.
   * For example:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * using PedAlgoFactory_t
   *   = opdet::factory::AlgorithmFactory<pmtana::PMTPedestalBase>;
   * using PedAlgoDecl = PedAlgoFactory_t::Decl;
   * 
   * PedAlgoFactory_t PedAlgoFactory;
   * PedAlgoFactory.declare(
   *   PedAlgoDecl<pmtana::PedAlgoEdges>{ "Edges" },
   *   PedAlgoDecl<pmtana::PedAlgoUB   >{ "UB"    }
   *   );
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * A constructor is likewise available to register algorithms.
   */
  template <typename FirstAlgo, typename... OtherAlgos>
  Factory_t& declare(
    opdet::factory::Decl<FirstAlgo> first,
    opdet::factory::Decl<OtherAlgos>... others
    );
  
  /// @}
  // --- END ---- Registration of algorithms ---------------------------------
  
  
  // --- BEGIN -- Queries and algorithm creation -----------------------------
  /// @name Queries and algorithm creation
  /// @{
  
  /**
   * @brief Returns an instance of algorithm `name` constructed with `pset`.
   * @param name declared name of the algorithm
   * @param pset the configuration of the algorithm
   * @return the newly created algorithm object
   * 
   * The algorithm is constructed with `pset` as only constructor argument.
   */
  std::unique_ptr<Algo_t> create
    (std::string const& name, fhicl::ParameterSet const& pset) const;
  
  /**
   * @brief Creates an instance of the algorithm constructed with `pset`.
   * @param pset the configuration of the algorithm
   * @return the newly created algorithm object
   * 
   * The type of algorithm is discovered from `pset` itself with the mechanism
   * documented in `setAlgorithmConfigurationKey()`.
   * 
   * For example:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * factory.setAlgorithmConfigurationKey("Name");
   * auto algo = factory.create(pset);
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * is equivalent to:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * auto algo = factory.create(pset.get<std::string>("Name"), pset);
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */
  std::unique_ptr<Algo_t> create(fhicl::ParameterSet const& pset) const;
  
  /// Returns a list with the names of supported algorithms.
  std::vector<std::string> names() const;
  
  /// Returns a string with the names of supported algorithms, joint.
  std::string names(std::string const& sep) const
    { return join(sep, names()); }
  
  /// }
  // --- END ---- Queries and algorithm creation -----------------------------
  
    private:
  
  /// Wrapper of an algorithm providing a polymorphic interface to create
  /// objects of the wrapped type.
  struct AlgoMaker;
  
  /**
   * @brief Standard algorithm maker class.
   * @tparam Algo the concrete algorithm class delivered by this object
   * 
   * This class implements `AlgoMaker` by overriding `makeAlgo()` method.
   * Its implementation simply constructs a new object of type `Algo` passing
   * the parameter set in the interface to the constructor of the `Algo` class.
   */
  template <typename Algo>
  struct AlgoMakerFor;
  
  /// Name of the configuration key with algorithm name (empty if not provided).
  std::string fAlgoNameKey;
  
  /// Registry of all maker classes.
  std::vector<std::unique_ptr<AlgoMaker>> fMakers;
  
  /// Returns the maker with the specified `name`, `nullptr` if none.
  AlgoMaker const* getMaker(std::string const& name) const;
  
  /// Adds the `maker` to the list of registered maker classes.
  Factory_t& registerMaker(std::unique_ptr<AlgoMaker>&& maker);
  
}; // class opdet::factory::AlgorithmFactory


// -----------------------------------------------------------------------------
// wrapping for algorithms to allow an interface to the framework
template <typename Base, typename FWTraits /* = FWInterfaceTraits<> */>
class opdet::factory::FWInterfacedIF {
    public:
  using Algo_t = Base;
  using Event_t = typename FWTraits::Event_t;
  using Module_t = typename FWTraits::Module_t;
  
  /// Access to the algorithm.
  Algo_t& algo() const { return *getAlgo(); }
  
  // --- BEGIN -- Framework hooks ----------------------------------------------
  void initialize(Module_t& module) { doInitialize(module); }
  void beginEvent(Event_t const& event) { doBeginEvent(event); }
  void endEvent(Event_t const& event) { doEndEvent(event); }
  // --- END ---- Framework hooks ----------------------------------------------
  
    private:
  
  virtual Algo_t* getAlgo() const = 0;
  
  virtual void doInitialize(Module_t& module) {}
  virtual void doBeginEvent(Event_t const& event) {}
  virtual void doEndEvent(Event_t const& event) {}
  
}; // opdet::factory::FWInterfacedIF


/// Base class for specialization.
template
  <typename Algo, typename Base, typename FWTraits /* = FWInterfaceTraits<> */>
class opdet::factory::FWInterfacedBase
  : public opdet::factory::FWInterfacedIF<Base, FWTraits>
{
  using RealAlgo_t = Algo;
  
    public:
  
  using Algo_t = Base;
  
    protected:
  
  FWInterfacedBase(fhicl::ParameterSet const& pset)
    : fAlgo{ std::make_unique<RealAlgo_t>(pset) } {}
  
  /// Returns the stored algorithm.
  virtual RealAlgo_t* getAlgo() const override { return fAlgo.get(); }
  
    private:
  
  std::unique_ptr<RealAlgo_t> fAlgo; ///< Instance of RealAlgo_t.
  
}; // FWInterfacedBase


template
  <typename Algo, typename Base, typename FWTraits /* = FWInterfaceTraits<> */>
class opdet::factory::FWInterfaced
  : public opdet::factory::FWInterfacedBase<Algo, Base, FWTraits>
{
    public:
  FWInterfaced(fhicl::ParameterSet const& pset)
    : FWInterfacedBase<Algo, Base, FWTraits>{ pset } {}
  
}; // FWInterfaced


// -----------------------------------------------------------------------------
// --- template implementation
// -----------------------------------------------------------------------------
// ---  opdet::factory::AlgorithmFactory::AlgoMaker
// -----------------------------------------------------------------------------
template <typename Base>
struct opdet::factory::AlgorithmFactory<Base>::AlgoMaker {
  
  std::string name; ///< The name associated to this algorithm.
  
  AlgoMaker(std::string name): name{ std::move(name) } {}
  
  /**
   * @brief Algorithm class construction.
   * @param pset the parameter set used to construct the algorithm object
   * @return a pointer to an object derived from `Base`, ready for use
   * 
   * The implementations will do whatever it takes to allocate and construct
   * the `Algo_t`-derived object they cover, based on the content of `pset`,
   * and return it.
   */
  virtual std::unique_ptr<Base> makeAlgo
    (fhicl::ParameterSet const& pset) const = 0;
  
  /// Comparison: sort by name in lexicographic order.
  bool operator< (AlgoMaker const& other) const
    { return name < other.name; }
  
}; // opdet::factory::AlgorithmFactory::AlgoMaker


template <typename Base>
template <typename Algo>
struct opdet::factory::AlgorithmFactory<Base>::AlgoMakerFor
  : opdet::factory::AlgorithmFactory<Base>::AlgoMaker
{
  
  AlgoMakerFor(std::string const& name): AlgoMaker{ name } {}
  
  std::unique_ptr<Base> makeAlgo
    (fhicl::ParameterSet const& pset) const override
    { return std::make_unique<Algo>(pset); }
  
}; // opdet::factory::AlgorithmFactory::AlgoMakerFor


// -----------------------------------------------------------------------------
// ---  opdet::factory::AlgorithmFactory
// -----------------------------------------------------------------------------
template <typename Base>
template <typename... Algos>
opdet::factory::AlgorithmFactory<Base>::AlgorithmFactory
  (std::string algoKey, opdet::factory::Decl<Algos>... algorithms)
  : fAlgoNameKey{ std::move(algoKey) }
{
  if constexpr(sizeof...(Algos) > 0) declare(std::move(algorithms)...);
}


// -----------------------------------------------------------------------------
template <typename Base>
template <typename Algo>
auto opdet::factory::AlgorithmFactory<Base>::declare(std::string name)
  -> Factory_t&
{
  return registerMaker(std::make_unique<AlgoMakerFor<Algo>>(std::move(name)));
} // opdet::factory::AlgorithmFactory::declare(string)


// -----------------------------------------------------------------------------
template <typename Base>
template <typename FirstAlgo, typename... OtherAlgos>
auto opdet::factory::AlgorithmFactory<Base>::declare(
  opdet::factory::Decl<FirstAlgo> first,
  opdet::factory::Decl<OtherAlgos>... others
) -> Factory_t&
{
  declare<FirstAlgo>(std::move(first.name));
  if constexpr(sizeof...(OtherAlgos) > 0) declare(std::move(others)...);
  return *this;
} // opdet::factory::AlgorithmFactory::declare(Decl)


// -----------------------------------------------------------------------------
template <typename Base>
auto opdet::factory::AlgorithmFactory<Base>::create
  (std::string const& name, fhicl::ParameterSet const& pset) const
  -> std::unique_ptr<Algo_t>
{
  AlgoMaker const* maker = getMaker(name);
  if (maker) return maker->makeAlgo(pset);
  
  throw cet::exception{ "AlgorithmFactory" }
    << "Unknown algorithm: '" << name << "'.\nSupported algorithms:\n - '"
    << names("'\n - '") << "'\n";
  
} // opdet::factory::AlgorithmFactory::create()


// -----------------------------------------------------------------------------
template <typename Base>
auto opdet::factory::AlgorithmFactory<Base>::create
  (fhicl::ParameterSet const& pset) const -> std::unique_ptr<Algo_t>
{
  return create(pset.get<std::string>(fAlgoNameKey), pset);
} // opdet::factory::AlgorithmFactory::create()


// -----------------------------------------------------------------------------
template <typename Base>
std::vector<std::string> opdet::factory::AlgorithmFactory<Base>::names() const {
  std::vector<std::string> v;
  for (auto const& maker: fMakers) v.push_back(maker->name);
  return v;
} // opdet::factory::AlgorithmFactory::names()


// -----------------------------------------------------------------------------
template <typename Base>
auto opdet::factory::AlgorithmFactory<Base>::getMaker
  (std::string const& name) const -> AlgoMaker const*
{
  // we could keep the makers ordered and do a binary search;
  // not really worth though
  for (std::unique_ptr<AlgoMaker> const& maker: fMakers)
    if (maker->name == name) return maker.get();
  return nullptr;
} // opdet::factory::AlgorithmFactory::getMaker()


// -----------------------------------------------------------------------------
template <typename Base>
auto opdet::factory::AlgorithmFactory<Base>::registerMaker
  (std::unique_ptr<AlgoMaker>&& maker) -> Factory_t&
{
  fMakers.push_back(std::move(maker)); 
  return *this;
} // opdet::factory::AlgorithmFactory::registerMaker()


// -----------------------------------------------------------------------------
// ---  other utilities
// -----------------------------------------------------------------------------
template <typename S, typename Coll>
S opdet::factory::join(S const& sep, Coll const& s) {
  using std::empty, std::begin, std::end;
  if (empty(s)) return S{};
  auto it = begin(s);
  auto const send = end(s);
  S cat { *it };
  while (++it != send) { cat += sep; cat += *it; }
  return cat;
} // join()

template <typename Coll>
std::string opdet::factory::join(const char* sep, Coll const& s)
  { return join<std::string>(sep, s); }

// -----------------------------------------------------------------------------

#endif // ICARUSCODE_PMT_OPRECO_ALGORITMS_OPRECOFACTORYSTUFF_H
