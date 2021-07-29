/**
 * @file   icaruscode/PMT/Trigger/Algorithms/details/EventInfoUtils.cxx
 * @brief  Class hosting selected information about the event (implementation).
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   May 15, 2020
 * @see    icaruscode/PMT/Trigger/Algorithms/details/EventInfoUtils.h
 */


// library header
#include "icaruscode/PMT/Trigger/Algorithms/details/EventInfoUtils.h"

// LArSoft libraries
#include "lardataalg/Utilities/quantities/spacetime.h" // microseconds
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/TPCGeo.h"

// C/C++ standard libraries
#include <utility> // std::move()
#include <cassert>


// -----------------------------------------------------------------------------
namespace {
  
  /// Returns in which active TPC volume `point` falls in (`nullptr` if none).
  geo::TPCGeo const* pointInActiveTPC
    (geo::GeometryCore const& geom, geo::Point_t const& point)
    {
      geo::TPCGeo const* tpc = geom.PositionToTPCptr(point);
      return (tpc && tpc->ActiveBoundingBox().ContainsPosition(point))
        ? tpc: nullptr;
    } // pointInActiveTPC()
  
  
  /// Helper tracking energy depositions.
  struct EnergyAccumulator {
    
    using TimeSpan_t = icarus::trigger::details::EventInfoExtractor::TimeSpan_t;
    using GeV = util::quantities::gigaelectronvolt;
    
    GeV totalEnergy { 0.0 }, inSpillEnergy { 0.0 }, inPreSpillEnergy { 0.0 };
    GeV activeEnergy { 0.0 }, inSpillActiveEnergy { 0.0 },
      inPreSpillActiveEnergy { 0.0 };
    
    EnergyAccumulator(
      TimeSpan_t inSpillTimes,
      TimeSpan_t inPreSpillTimes,
      geo::GeometryCore const& geom
      )
      : fGeom(geom)
      , fInSpillTimes(inSpillTimes)
      , fInPreSpillTimes(inPreSpillTimes)
      {}
    
    void add(
      GeV energy,
      detinfo::timescales::simulation_time time,
      geo::Point_t const& location
      ) {
        
        bool const inSpill
          = (time >= fInSpillTimes.first) && (time <= fInSpillTimes.second);
        bool const inPreSpill =
          (time >= fInPreSpillTimes.first) && (time <= fInPreSpillTimes.second);
        
        totalEnergy += energy;
        if (inSpill) inSpillEnergy += energy;
        if (inPreSpill) inPreSpillEnergy += energy;
        
        if (pointInActiveTPC(location)) {
          activeEnergy += energy;
          if (inSpill) inSpillActiveEnergy += energy;
          if (inPreSpill) inPreSpillActiveEnergy += energy;
        }
      } // add()
    
      private:
    geo::GeometryCore const& fGeom; ///< Geometry service provider.
    
    /// Start and stop time for "in spill" label.
    TimeSpan_t const fInSpillTimes;
    
    /// Start and stop time for "pre-spill" label.
    TimeSpan_t const fInPreSpillTimes;
    
    /// Returns in which active TPC volume `point` falls in (`nullptr` if none).
    geo::TPCGeo const* pointInActiveTPC(geo::Point_t const& point) const
      { return ::pointInActiveTPC(fGeom, point); }

  }; // EnergyAccumulator
  
} // local namespace


// -----------------------------------------------------------------------------
icarus::trigger::details::EventInfoExtractor::EventInfoExtractor(
  std::vector<art::InputTag> truthTags,
  EDepTags_t edepTags,
  TimeSpan_t inSpillTimes,
  TimeSpan_t inPreSpillTimes,
  geo::GeometryCore const& geom,
  detinfo::DetectorPropertiesData const* detProps,
  detinfo::DetectorTimings const* detTimings,
  std::string logCategory
  )
  : fGeneratorTags(std::move(truthTags))
  , fEnergyDepositTags(std::move(edepTags))
  , fLogCategory(std::move(logCategory))
  , fGeom(geom)
  , fDetProps(detProps)
  , fDetTimings(detTimings)
  , fInSpillTimes(std::move(inSpillTimes))
  , fInPreSpillTimes(std::move(inPreSpillTimes))
{
} // icarus::trigger::details::EventInfoExtractor::EventInfoExtractor()


// -----------------------------------------------------------------------------
void icarus::trigger::details::EventInfoExtractor::fillGeneratorInfo
  (EventInfo_t& info, simb::MCTruth const& truth) const
{
  
  if (truth.NeutrinoSet()) fillGeneratorNeutrinoInfo(info, truth);
  
} // icarus::trigger::details::EventInfoExtractor::fillGeneratorInfo()


// -----------------------------------------------------------------------------
void icarus::trigger::details::EventInfoExtractor::fillGeneratorNeutrinoInfo
  (EventInfo_t& info, simb::MCTruth const& truth) const
{
  if (!truth.NeutrinoSet()) return;
  
  simulation_time const interactionTime = getInteractionTime(truth);
  
  if ((info.nVertices() == 0) || (interactionTime < info.InteractionTime()))
    setMainGeneratorNeutrinoInfo(info, truth);
  else
    addGeneratorNeutrinoInfo(info, truth);
  
} // icarus::trigger::details::EventInfoExtractor::fillGeneratorNeutrinoInfo()


// -----------------------------------------------------------------------------
void icarus::trigger::details::EventInfoExtractor::setMainGeneratorNeutrinoInfo
  (EventInfo_t& info, simb::MCTruth const& truth) const
{
  /*
   * Sets the full information of the event, overwriting everything.
   * 
   * Except that the vertex is just inserted into the vertex list, as the
   * first entry.
   */
  
  using GeV = util::quantities::gigaelectronvolt;
  
  //
  // interaction flavor (nu_mu, nu_e)
  // interaction type (CC, NC)
  //

  simb::MCParticle const& nu = truth.GetNeutrino().Nu();
  info.SetNeutrinoPDG(nu.PdgCode());
  info.SetInteractionType(truth.GetNeutrino().InteractionType());
  info.SetInteractionTime(getInteractionTime(truth));

  info.SetNeutrinoEnergy(GeV{ nu.E() });
  info.SetLeptonEnergy(GeV{ truth.GetNeutrino().Lepton().E() });
  //info.SetNucleonEnergy(truth.GetNeutrino().HitNuc().E());

  switch (nu.PdgCode()) {
    case 14:
    case -14:
      info.SetNu_mu(true);
      break;
    case 12:
    case -12:
      info.SetNu_e(true);
      break;
  }

  switch (truth.GetNeutrino().CCNC()) {
    case simb::kCC: info.AddWeakChargedCurrentInteractions(); break;
    case simb::kNC: info.AddWeakNeutralCurrentInteractions(); break;
    default:
      mf::LogWarning(fLogCategory)
        << "Unexpected NC/CC flag (" << truth.GetNeutrino().CCNC() << ")";
  } // switch   
  
  // we do not trust the vertex (`GvX()`) of the neutrino particle,
  // since GenieHelper does not translate the vertex
  // of some of the particles from GENIE to detector frame;
  // trajectory is always translated:
  geo::Point_t const vertex { nu.EndX(), nu.EndY(), nu.EndZ() };
  info.InsertVertex(vertex, 0U);
  
  geo::TPCGeo const* tpc = pointInActiveTPC(vertex);
  if (tpc) info.SetInActiveVolume();
  
} // icarus::trigger::details::EventInfoExtractor::setMainGeneratorNeutrinoInfo()


// -----------------------------------------------------------------------------
void icarus::trigger::details::EventInfoExtractor::addGeneratorNeutrinoInfo
  (EventInfo_t& info, simb::MCTruth const& truth) const
{
  /*
   * The only update to the record so far is the addition of the vertex of the
   * interaction at the end of the current list.
   * No information is overwritten.
   * 
   */
  
  simb::MCParticle const& nu = truth.GetNeutrino().Nu();
  
  // we do not trust the vertex (`GvX()`) of the neutrino particle,
  // since GenieHelper does not translate the vertex
  // of some of the particles from GENIE to detector frame;
  // trajectory is always translated:
  geo::Point_t const vertex { nu.EndX(), nu.EndY(), nu.EndZ() };
  info.AddVertex(vertex);
  
} // icarus::trigger::details::EventInfoExtractor::addGeneratorNeutrinoInfo()


// -----------------------------------------------------------------------------
void icarus::trigger::details::EventInfoExtractor::addEnergyDepositionInfo
  (EventInfo_t& info, std::vector<sim::SimEnergyDeposit> const& energyDeposits)
  const
{
  using MeV = util::quantities::megaelectronvolt;
  
  auto Eacc = EnergyAccumulator(fInSpillTimes, fInPreSpillTimes, fGeom);
  
  for (sim::SimEnergyDeposit const& edep: energyDeposits) {
    
    Eacc.add(
        MeV{ edep.Energy() } // assuming it's stored in MeV
      , detinfo::timescales::simulation_time{ edep.Time() }
      , edep.MidPoint()
      );
    
  } // for all energy deposits in the data product
  
  info.SetDepositedEnergy
    (info.DepositedEnergy() + Eacc.totalEnergy);
  info.SetDepositedEnergyInSpill
    (info.DepositedEnergyInSpill() + Eacc.inSpillEnergy);
  info.SetDepositedEnergyInPreSpill
    (info.DepositedEnergyInPreSpill() + Eacc.inPreSpillEnergy);
  info.SetDepositedEnergyInActiveVolume
    (info.DepositedEnergyInActiveVolume() + Eacc.activeEnergy);
  info.SetDepositedEnergyInSpillInActiveVolume
    (info.DepositedEnergyInSpillInActiveVolume() + Eacc.inSpillActiveEnergy);
  info.SetDepositedEnergyInPreSpillInActiveVolume(
    info.DepositedEnergyInPreSpillInActiveVolume() + Eacc.inPreSpillActiveEnergy
    );
  
} // icarus::trigger::details::EventInfoExtractor::addEnergyDepositionInfo()


// -----------------------------------------------------------------------------
void icarus::trigger::details::EventInfoExtractor::addEnergyDepositionInfo
  (EventInfo_t& info, std::vector<sim::SimChannel> const& channels) const
{
  assert(fDetProps);
  assert(fDetTimings);
  
  using MeV = util::quantities::megaelectronvolt;
  
  auto Eacc = EnergyAccumulator(fInSpillTimes, fInPreSpillTimes, fGeom);
  
  double const driftVel = fDetProps->DriftVelocity(); // cm/us
  
  for (sim::SimChannel const& channel: channels) {
    
    // only channels on any of the first induction planes
    std::vector<geo::WireID> const& wires
      = fGeom.ChannelToWire(channel.Channel());
    if (empty(wires)) continue; // ghost channel or something; move on
    if (wires.front().Plane != 0) continue; // not the first induction plane
    
    geo::PlaneGeo const& plane = fGeom.Plane(wires.front());
    
    for (auto const& [ tdc, IDEs ]: channel.TDCIDEMap()) {
      
      // collection tick: includes also drift time, diffusion and what-not
      detinfo::timescales::electronics_tick const tick { tdc };
      detinfo::timescales::simulation_time const time
         = fDetTimings->toSimulationTime(tick);
      
      for (sim::IDE const& IDE: IDEs) {
        MeV const energy { IDE.energy }; // assuming it's stored in MeV
        geo::Point_t const location { IDE.x, IDE.y, IDE.z };
        
        // tentative estimation of drift length:
        double const d = plane.DistanceFromPlane(location); // cm
        util::quantities::intervals::microseconds const driftTime
          { d / driftVel };
        
        Eacc.add(energy, time - driftTime, location);
      } // all deposits at this time tick
      
    } // all time ticks in this channel
    
  } // for all channels
  
  info.SetDepositedEnergy
    (info.DepositedEnergy() + Eacc.totalEnergy);
  info.SetDepositedEnergyInSpill
    (info.DepositedEnergyInSpill() + Eacc.inSpillEnergy);
  info.SetDepositedEnergyInPreSpill
    (info.DepositedEnergyInPreSpill() + Eacc.inPreSpillEnergy);
  info.SetDepositedEnergyInActiveVolume
    (info.DepositedEnergyInActiveVolume() + Eacc.activeEnergy);
  info.SetDepositedEnergyInSpillInActiveVolume
    (info.DepositedEnergyInSpillInActiveVolume() + Eacc.inSpillActiveEnergy);
  info.SetDepositedEnergyInPreSpillInActiveVolume(
    info.DepositedEnergyInPreSpillInActiveVolume() + Eacc.inPreSpillActiveEnergy
    );
  
} // icarus::trigger::details::EventInfoExtractor::addEnergyDepositionInfo()


// -----------------------------------------------------------------------------
geo::TPCGeo const* icarus::trigger::details::EventInfoExtractor::pointInTPC
  (geo::Point_t const& point) const
{
  return fGeom.PositionToTPCptr(point);
} // icarus::trigger::TriggerEfficiencyPlotsBase::pointInTPC()


//------------------------------------------------------------------------------
geo::TPCGeo const*
icarus::trigger::details::EventInfoExtractor::pointInActiveTPC
  (geo::Point_t const& point) const
{
  return ::pointInActiveTPC(fGeom, point);
} // icarus::trigger::TriggerEfficiencyPlotsBase::pointInActiveTPC()


// -----------------------------------------------------------------------------
auto icarus::trigger::details::EventInfoExtractor::getInteractionTime
  (simb::MCTruth const& truth) const -> simulation_time
{
  assert(truth.NeutrinoSet());
  simb::MCParticle const& nu = truth.GetNeutrino().Nu();
  return simulation_time{ nu.EndT() };
} // icarus::trigger::details::EventInfoExtractor::getInteractionTime()


// -----------------------------------------------------------------------------
// --- icarus::trigger::details::EventInfoExtractorMaker
// -----------------------------------------------------------------------------
icarus::trigger::details::EventInfoExtractorMaker::EventInfoExtractorMaker(
  std::vector<art::InputTag> truthTags,
  EDepTags_t edepTags,
  geo::GeometryCore const& geom,
  detinfo::DetectorPropertiesData const* detProps,
  detinfo::DetectorTimings const* detTimings,
  std::string logCategory
  )
  : fGeneratorTags(std::move(truthTags))
  , fEnergyDepositTags(std::move(edepTags))
  , fLogCategory(std::move(logCategory))
  , fGeom(geom)
  , fDetProps(detProps)
  , fDetTimings(detTimings)
  {}


// -----------------------------------------------------------------------------
icarus::trigger::details::EventInfoExtractor
icarus::trigger::details::EventInfoExtractorMaker::make
  (TimeSpan_t inSpillTimes, TimeSpan_t inPreSpillTimes) const
{
  return EventInfoExtractor{
    fGeneratorTags, fEnergyDepositTags,
    inSpillTimes, inPreSpillTimes,
    fGeom, fDetProps, fDetTimings, fLogCategory
    };
} // icarus::trigger::details::EventInfoExtractorMaker::make()


// -----------------------------------------------------------------------------
