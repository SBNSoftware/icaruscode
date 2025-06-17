/**
 * @file   icaruscode/Generators/SelectCathodeCrossingGenParticles_module.cc
 * @brief  Implements `sbn::SelectCathodeCrossingGenParticles` filter module.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   June 12, 2025
 *
 */

// ICARUS/LArSoft libraries
#include "icarusalg/Utilities/AtomicPassCounter.h"
#include "icarusalg/Utilities/PassCounter.h"
#include "sbnalg/Utilities/PlaneCrossers.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/geo_vectors_utils.h" // geo::vect::convertTo()
// #include "larcorealg/Geometry/geo_vectors_utils_TVector.h"
#include "larcorealg/CoreUtils/enumerate.h"
#include "larcorealg/CoreUtils/counter.h"
#include "lardataalg/MCDumpers/MCDumperUtils.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
// framework libraries
#include "art/Framework/Core/SharedFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"
// #include "canvas/Utilities/Exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Atom.h"
// #include "cetlib_except/exception.h"

// C/C++ standard libraries
// #include <algorithm> // std::count()
#include <string>
#include <utility> // std::move()
#include <vector>


//------------------------------------------------------------------------------
/**
 * @brief Requires the presence of generated particles crossing a cathode.
 *
 * This module passes only the events with a minimum number of generated
 * particles (from `simb::MCTruth`) predicted to cross any cathode.
 * 
 * All particles in the final state (status `1`) of the specified generated
 * records are considered. The number of particles whose predicted trajectory
 * crosses any of the cathodes of the detector, plus a margin, is compared to
 * the configured requirement.
 * 
 * The counts are currently global: if a particle crosses the cathode of
 * cryostat 1 and another the cathode of cryostat 2, the event passes the
 * requirement of _two_ cathode-crossing particles.
 * 
 * The crossing prediction is based on whether the predicted particle trajectory
 * crosses a rectangularly shaped cathode, with some margin.
 * The trajectory is built as simple inertial extrapolation of the particle at
 * generator time.
 * 
 * Cathode coordinates are not provided by the detector geometry description,
 * so for each TPC a cathode is deduced (see below).
 * On top of the cathode dimensions, each side of each cathode is extended by
 * a configured amount. For example, if a cathode is 9 &times; 3.4 m&sup2;
 * (L &times; H), with `CathodeLengthWiggle` of `0.1` the cathode length (9 m)
 * is extended by 0.9 m on one side _and_ 0.9 m on the other side, for a total
 * width of 10.8 m (that is, 20% longer).
 * 
 * 
 * Cathode size
 * -------------
 * 
 * Cathode coordinates are not provided by the detector geometry description,
 * so _for each TPC_ a cathode is deduced.
 * 
 * The definition of the cathode surface is as follow:
 *  * a cross section of the active volume is determined as the rectangle
 *    with sides along `geo::TPCGeo::LengthDir()` and `geo::TPCGeo::HeightDir()`
 *    and sizes `geo::TPCGeo::ActiveLength()` and `geo::TPCGeo::ActiveHeight()`,
 *    respectively, around the central point of the active volume
 *    (`geo::TPCGeo::GetActiveVolumeCenter()`).
 *  * cathode surface is set translating that rectangle by
 *    `geo::TPCGeo::ActiveHalfWidth()` along the direction `geo::TPCGeo::WidthDir()`,
 *    which is supposed to be opposite to the drift direction.
 *    
 * 
 * 
 *
 * Configuration parameters
 * =========================
 *
 * * `GeneratorTags` (input tags, mandatory): list of tags of generated particle
 *    records to consider.
 * * `MinimumCrossingParticles` (integral, default: `1`): minimum number of
 *    cathode-crossing particles required for the event to pass the filter.
 * * `CathodeLengthWiggle`, `CathodeHeightWiggle` (real, default: `0.0`):
 *    fraction of each dimension the cathode is extended for the purpose of
 *    determining the crossing of the particle (see the documentation above).
 * * `LogCategory` (string, default: `SelectCathodeCrossingGenParticles`):
 *   message facility stream name where to write console messages to.
 *
 * 
 * Input data products
 * ====================
 * 
 * * `std::vector<simb::MCTruth>` (`GeneratorTags`, multiple supported):
 *     the tags of the generated events to be included in the test.
 * 
 * 
 * Service dependencies
 * =====================
 * 
 * * `Geometry`: to learn the position of the cathodes.
 * 
 * 
 * 
 * Output data products
 * =====================
 * 
 * None.
 * 
 */
class SelectCathodeCrossingGenParticles: public art::SharedFilter {
  
    public:
  
  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    fhicl::Sequence<art::InputTag> GeneratorTags{
      Name{ "GeneratorTags" },
      Comment{ "list of generated particle tags to check" },
      std::vector{ art::InputTag{ "generator" } }
      };
    
    fhicl::Atom<unsigned int> MinimumCrossingParticles{
      Name{ "MinimumCrossingParticles" },
      Comment{ "minimum required number of cathode-crossing particles" },
      1U
      };
    
    fhicl::Atom<double> CathodeLengthWiggle{
      Name{ "CathodeLengthWiggle" },
      Comment
        { "extend each cathode length by this fraction (negative contracts)" },
      0.0
      };
    
    fhicl::Atom<double> CathodeHeightWiggle{
      Name{ "CathodeHeightWiggle" },
      Comment
        { "extend each cathode height by this fraction (negative contracts)" },
      0.0
      };
    
    fhicl::Atom<std::string> LogCategory{
      Name{ "LogCategory" },
      Comment{ "name of the message facility stream used by the module" },
      "SelectCathodeCrossingGenParticles"
      };
    
  }; // Config
  
  using Parameters = art::SharedFilter::Table<Config>;
  
  
  SelectCathodeCrossingGenParticles
    (Parameters const& params, art::ProcessingFrame const&);
  
  /// Evaluate the filtering logic.
  virtual bool filter(art::Event& event, art::ProcessingFrame const&) override;
  
  /// Prints a summary.
  virtual void endJob(art::ProcessingFrame const&) override;
  
    private:
  
  /// The algorithm storing plane geometry and providing intersection with lines.
  using PlaneCrosserAlg_t = util::PlaneCrossers<geo::Point_t>;
  
  /// Information identifying the cathode rectangle of one TPC.
  class CathodePlane_t {
    
    PlaneCrosserAlg_t fCrossFinder;
    
    geo::TPCGeo const* fTPC = nullptr;
    
      public:
    
    CathodePlane_t(
      geo::Point_t center, geo::Vector_t length, geo::Vector_t height,
      geo::TPCGeo const& TPC
      )
      : fCrossFinder{ std::move(center), std::move(length), std::move(height) }
      , fTPC{ &TPC }
      {}
    
    geo::Point_t center() const { return fCrossFinder.center(); }
    
    /// Versor of the length direction.
    geo::Vector_t lengthDir() const { return fCrossFinder.Uaxis().unit(); }
    /// Versor of the height direction.
    geo::Vector_t heightDir() const { return fCrossFinder.Vaxis().unit(); }
    
    double length() const { return fCrossFinder.Uaxis().R(); }
    double height() const { return fCrossFinder.Vaxis().R(); }
    
    geo::TPCGeo const& TPC() const { return *fTPC; }
    
    PlaneCrosserAlg_t crossFinder() const
      { return fCrossFinder; }
    
    /// Finds the crossing between the cathode and the specified line.
    /// @see util::PlaneCrossers::findCrossing()
    PlaneCrosserAlg_t::CrossingInfo findCrossing
      (geo::Point_t const& point, geo::Vector_t const& dir) const
      { return fCrossFinder.findCrossing(point, dir); }
    
  }; // CathodePlane_t
  
  
  /// Value returned to represent the possible crossing of a cathode.
  struct CathodeCrossInfo {
    
    /// Which cathode was crossed (if any).
    CathodePlane_t const* cathode = nullptr;
    
    /// Information about the crossing with a line.
    PlaneCrosserAlg_t::CrossingInfo crossing;
    
    /// How many `length()`s the crossing point is from cathode center.
    double cathodeLengthUnits() const { return crossing.u; }
    /// How many `height()`s the crossing point is from cathode center.
    double cathodeHeightUnits() const { return crossing.v; }
    /// How many particle direction units the crossing point is from its position.
    double particleDirUnits() const { return crossing.line; }
    
    /// Returns whether any cathode was crossed.
    bool crossed() const noexcept { return crossing; }
    
    /// Returns whether any cathode was crossed.
    operator bool() const noexcept { return crossed(); }
    
    /// Access to the cathode information (undefined behaviour if `crossed()` is `false`).
    CathodePlane_t const* operator-> () const { return cathode; }
    
    /// Returns the cathode/line intersection point (undefined behaviour if none).
    geo::Point_t intersectionPoint() const
      {
        assert(cathode);
        return cathode->crossFinder().intersectionPoint(crossing);
      }
    
  }; // CathodeCrossInfo
  
  
  /// Counts of particle categories to determine the satisfaction of requirements.
  struct Counters_t {
    
    /// Number of truth records seen.
    unsigned int nTruthRecords = 0;
    
    /// Number of particles from any interaction crossing any cathode.
    icarus::ns::util::PassCounter<> cathodeCrossing;
    
  }; // Counters_t
  
  
  // --- BEGIN -- Configuration parameters -------------------------------------
  
  std::vector<art::InputTag> const fGeneratorTags; ///< Input data products.
  
  /// Minimum number of requested cathode-crossing particles.
  unsigned int const fMinimumCrossingPartcles;
  
  double const fCathodeLengthWiggle; ///< Fractional allowance on cathode length.
  double const fCathodeHeightWiggle; ///< Fractional allowance on cathode height.
  
  std::string const fLogCategory; ///< Name of message facility stream.
  
  // --- END ---- Configuration parameters -------------------------------------
  
  // --- BEGIN -- Caches -------------------------------------------------------
  
  // we don't really need to group by cryostat/TPC, but it comes for ~free
  std::vector<CathodePlane_t> const fCathodes; ///< All cathodes.
  
  // --- END ---- Caches -------------------------------------------------------
  
  /// Count of passed and seen events.
  icarus::ns::util::AtomicPassCounter<> fPassed;
  
  /**
   * @brief Analyzes the `truth` information and updates the `counts` accordingly.
   * @param truth the generated record to extract information from
   * @param[out] counters counters to be updated
   * 
   * This method extracts information but does not take any decision.
   */
  void parseTruth(simb::MCTruth const& truth, Counters_t& counters) const;
  
  
  /**
   * @brief Returns one TPC whose cathode is crossed by a particle.
   * @param pos starting position of the particle
   * @param dir direction of the particle
   * @return information on the crossing, converts to `false` if no crossing
   * 
   * If there are multiple cathodes crossed by the particle, which one is
   * returned is undefined.
   */
  CathodeCrossInfo findCrossingTPCcathode
    (geo::Point_t const& pos, geo::Vector_t const& dir) const;
  
  
  /// Builds and returns information of the cathode of the specified `TPC`.
  static CathodePlane_t buildCathode(geo::TPCGeo const& TPC);
  
  /// Builds and returns all cathodes, one per logical TPC in the geometry.
  static std::vector<CathodePlane_t> buildCathodePlanes
    (geo::GeometryCore const& geom);
  
}; // class SelectCathodeCrossingGenParticles


//------------------------------------------------------------------------------

SelectCathodeCrossingGenParticles::SelectCathodeCrossingGenParticles
  (Parameters const& params, art::ProcessingFrame const&)
  : art::SharedFilter{ params }
  // configuration
  , fGeneratorTags          { params().GeneratorTags() }
  , fMinimumCrossingPartcles{ params().MinimumCrossingParticles() }
  , fCathodeLengthWiggle    { params().CathodeLengthWiggle() * 2.0 } // x2 for
  , fCathodeHeightWiggle    { params().CathodeHeightWiggle() * 2.0 } // half lengths
  , fLogCategory            { params().LogCategory() }
  // caches
  , fCathodes{ buildCathodePlanes(*lar::providerFrom<geo::Geometry>()) }
{
  
  async<art::InEvent>();
  
  //
  // consume declaration
  //
  for (art::InputTag const& tag: fGeneratorTags)
    consumes<std::vector<simb::MCTruth>>(tag);
  
  //
  // dump configuration
  //
  mf::LogInfo log{ fLogCategory };
  log << "Configuration:";
  
  log << "\n * check " << fGeneratorTags.size() << " generated data products:";
  for (art::InputTag const& tag: fGeneratorTags)
    log << " '" << tag.encode() << "'";
  
  log << "\n * minimum number of cathode-crossing particles: "
    << fMinimumCrossingPartcles;
  
  log << "\n * fractional allowance on cathode dimensions: width "
    << (fCathodeLengthWiggle / 2 * 100.0) << "%, height "
    << (fCathodeHeightWiggle / 2 * 100.0) << "%";
  
  //
  // debug dump
  //
  {
    mf::LogTrace log{ fLogCategory };
    log << "Found " << fCathodes.size() << " cathode planes:";
    for (CathodePlane_t const& cathode: fCathodes) {
      log << "\n - " << cathode.TPC().ID()
        << ": ( " << cathode.length() << " x " << cathode.height() << " ) cm at "
        << cathode.center();
    } // for
  }
  
} // SelectCathodeCrossingGenParticles::SelectCathodeCrossingGenParticles()


//------------------------------------------------------------------------------
bool SelectCathodeCrossingGenParticles::filter
  (art::Event& event, art::ProcessingFrame const&)
{
  
  //
  // collect information
  //
  mf::LogDebug{ fLogCategory } << "Evaluation started.";
  
  Counters_t counters;
  
  for (art::InputTag const& truthTag: fGeneratorTags) {
    
    auto const& truths = event.getProduct<std::vector<simb::MCTruth>>(truthTag);
    mf::LogTrace{ fLogCategory }
      << "Now serving: " << truthTag.encode() << " (" << truths.size() << " records)";
    
    for (auto const& [ iTruth, truth ]: util::enumerate(truths)) {
      mf::LogTrace{ fLogCategory }
        << "Truth #" << iTruth << " from '" << truthTag.encode() << ": "
        << truth.NParticles() << " particles";
      
      parseTruth(truth, counters);
      
    } // for truth record
    
  } // for truth tags
  
  
  //
  // draw conclusions
  //
  unsigned int const nCathodeCrossing = counters.cathodeCrossing.passed();
  bool const accepted = (nCathodeCrossing >= fMinimumCrossingPartcles);
  
  mf::LogInfo{ fLogCategory }
    << event.id() << ": " << counters.cathodeCrossing.passed()
    << " / " << counters.cathodeCrossing.total()
    << " generated particles crossed a cathode.";
  
  fPassed.add(accepted);
  
  return accepted;
  
} // SelectCathodeCrossingGenParticles::filter()


//------------------------------------------------------------------------------
void SelectCathodeCrossingGenParticles::endJob(art::ProcessingFrame const&) {
  
  if (fPassed.total() > 0) {
    mf::LogInfo{ fLogCategory }
      << fPassed.passed() << "/" << fPassed.total()
      << " events (" << (fPassed.passed() * 100.0 / fPassed.total())
      << "%) passed the requirements.";
  }
  
} // SelectCathodeCrossingGenParticles::endJob()


//------------------------------------------------------------------------------
void SelectCathodeCrossingGenParticles::parseTruth
  (simb::MCTruth const& truth, Counters_t& counters) const
{
  
  ++(counters.nTruthRecords);
  
  for (auto const iPart: util::counter(truth.NParticles())) {
    simb::MCParticle const& particle = truth.GetParticle(iPart);
    
    bool const isFinal = particle.StatusCode() == 1;
    if (!isFinal) continue;
    
    geo::Point_t const pos = geo::vect::toPoint(particle.EndPosition().Vect());
    geo::Vector_t const dir = geo::vect::toVector(particle.EndMomentum().Vect());
    
    CathodeCrossInfo const crossedTPCcathode = findCrossingTPCcathode(pos, dir);
    
    {
      mf::LogTrace log{ fLogCategory };
      log << "Particle #" << iPart << " (" << sim::ParticleName(particle.PdgCode())
        << " at " << pos << " toward " << dir << ")";
      if (crossedTPCcathode) log << " crosses " << crossedTPCcathode->TPC().ID();
      else                   log << " does not cross any TPC cathode";
    } // block
    
    counters.cathodeCrossing.add(crossedTPCcathode.crossed());
    
  } // for particle
  
} // SelectCathodeCrossingGenParticles::parseTruth()


//------------------------------------------------------------------------------
SelectCathodeCrossingGenParticles::CathodeCrossInfo
SelectCathodeCrossingGenParticles::findCrossingTPCcathode
  (geo::Point_t const& pos, geo::Vector_t const& dir) const
{
  // try all, one after another
  mf::LogTrace{ fLogCategory }
    << "Testing particle from " << pos << " cm toward " << dir;
  
  for (CathodePlane_t const& plane: fCathodes) {
    CathodeCrossInfo const crossing {
        /* .cathode  = */ &plane
      , /* .crossing = */ plane.findCrossing(pos, dir)
      };
    
    if (crossing) {
      mf::LogTrace{ fLogCategory } << " - TPC " << plane.TPC().ID()
        << " met at " << crossing.intersectionPoint()
        << " cm (plane: L=" << crossing.cathodeLengthUnits()
        << " cm, H=" << crossing.cathodeHeightUnits()
        << "; P=" << crossing.particleDirUnits() << " cm)";
    }
    else {
      mf::LogTrace{ fLogCategory } << " - TPC " << plane.TPC().ID()
        << " not crossed";
    }
    
    if (!crossing) continue;
    
    // is crossing point contained in the (expanded) cathode area?
    if (std::abs(crossing.cathodeHeightUnits()) > 1. + fCathodeLengthWiggle) continue;
    if (std::abs(crossing.cathodeHeightUnits()) > 1. + fCathodeHeightWiggle) continue;
    
    // is the point in the future of the particle?
    if (crossing.particleDirUnits() < 0) continue;
    
    return crossing;
  } // for cathode planes
  
  return {};
} // SelectCathodeCrossingGenParticles::findCrossingTPCcathode()


//------------------------------------------------------------------------------
SelectCathodeCrossingGenParticles::CathodePlane_t
SelectCathodeCrossingGenParticles::buildCathode(geo::TPCGeo const& TPC) {
  
  // the width direction is not guaranteed to point from anode to cathode nor
  // the opposite, so we use the drift direction, which should be parallel.
  return {
      // center
      TPC.GetActiveVolumeCenter() - TPC.DriftDir() * TPC.ActiveHalfWidth()
      // length
    , TPC.LengthDir() * TPC.ActiveLength()
      // height
    , TPC.HeightDir() * TPC.ActiveHeight()
    , TPC
  };
  
} // SelectCathodeCrossingGenParticles::buildCathode()


//------------------------------------------------------------------------------
std::vector<SelectCathodeCrossingGenParticles::CathodePlane_t>
SelectCathodeCrossingGenParticles::buildCathodePlanes
  (geo::GeometryCore const& geom)
{
  std::vector<SelectCathodeCrossingGenParticles::CathodePlane_t> cathodePlanes;
  cathodePlanes.reserve(geom.TotalNTPC());
  
  for (auto const& TPC: geom.Iterate<geo::TPCGeo>())
    cathodePlanes.emplace_back(buildCathode(TPC));
  
  return cathodePlanes;
} // SelectCathodeCrossingGenParticles::buildCathodePlanes()


//------------------------------------------------------------------------------
DEFINE_ART_MODULE(SelectCathodeCrossingGenParticles)


//------------------------------------------------------------------------------

