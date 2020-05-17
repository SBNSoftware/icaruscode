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
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/TPCGeo.h"

// C/C++ standard libraries
#include <utility> // std::move()


// -----------------------------------------------------------------------------
icarus::trigger::details::EventInfoExtractor::EventInfoExtractor(
  std::vector<art::InputTag> truthTags,
  std::vector<art::InputTag> edepTags,
  TimeSpan_t inSpillTimes,
  geo::GeometryCore const& geom,
  std::string logCategory
  )
  : fGeneratorTags(std::move(truthTags))
  , fEnergyDepositTags(std::move(edepTags))
  , fLogCategory(std::move(logCategory))
  , fGeom(geom)
  , fInSpillTimes(std::move(inSpillTimes))
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
  
  using GeV = util::quantities::gigaelectronvolt;
  
  //
  // interaction flavor (nu_mu, nu_e)
  // interaction type (CC, NC)
  //

  simb::MCParticle const& nu = truth.GetNeutrino().Nu();
  info.SetNeutrinoPDG(nu.PdgCode());
  info.SetInteractionType(truth.GetNeutrino().InteractionType());

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
  info.AddVertex(vertex);
  
  geo::TPCGeo const* tpc = pointInActiveTPC(vertex);
  if (tpc) info.SetInActiveVolume();
  
} // icarus::trigger::details::EventInfoExtractor::fillGeneratorNeutrinoInfo()


// -----------------------------------------------------------------------------
void icarus::trigger::details::EventInfoExtractor::addEnergyDepositionInfo
  (EventInfo_t& info, std::vector<sim::SimEnergyDeposit> const& energyDeposits)
  const
{
  using MeV = util::quantities::megaelectronvolt;
  using GeV = util::quantities::gigaelectronvolt;
  
  //
  // propagation in the detector
  //
  GeV totalEnergy { 0.0 }, inSpillEnergy { 0.0 };
  GeV activeEnergy { 0.0 }, inSpillActiveEnergy { 0.0 };
  
  for (sim::SimEnergyDeposit const& edep: energyDeposits) {
    
    MeV const e { edep.Energy() }; // assuming it's stored in MeV
    
    detinfo::timescales::simulation_time const t { edep.Time() };
    bool const inSpill
      = (t >= fInSpillTimes.first) && (t <= fInSpillTimes.second);
    
    totalEnergy += e;
    if (inSpill) inSpillEnergy += e;
    
    if (pointInActiveTPC(edep.MidPoint())) {
      activeEnergy += e;
      if (inSpill) inSpillActiveEnergy += e;
    }
    
  } // for all energy deposits in the data product
  
  info.SetDepositedEnergy
    (info.DepositedEnergy() + totalEnergy);
  info.SetDepositedEnergyInSpill
    (info.DepositedEnergyInSpill() + inSpillEnergy);
  info.SetDepositedEnergyInActiveVolume
    (info.DepositedEnergyInActiveVolume() + activeEnergy);
  info.SetDepositedEnergyInSpillInActiveVolume
    (info.DepositedEnergyInSpillInActiveVolume() + inSpillActiveEnergy);
  
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
  geo::TPCGeo const* tpc = pointInTPC(point);
  return
    (tpc && tpc->ActiveBoundingBox().ContainsPosition(point))? tpc: nullptr;
} // icarus::trigger::TriggerEfficiencyPlotsBase::pointInActiveTPC()


// -----------------------------------------------------------------------------
