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
  EdepTags_t edepTags,
  TimeSpan_t inSpillTimes,
  TimeSpan_t inPreSpillTimes,
  geo::GeometryCore const& geom,
  std::string logCategory
  )
  : fGeneratorTags(std::move(truthTags))
  , fEnergyDepositTags(std::move(edepTags))
  , fLogCategory(std::move(logCategory))
  , fGeom(geom)
  , fInSpillTimes(std::move(inSpillTimes))
  , fInPreSpillTimes(std::move(inPreSpillTimes))
{
} // icarus::trigger::details::EventInfoExtractor::EventInfoExtractor()


// -----------------------------------------------------------------------------
void icarus::trigger::details::EventInfoExtractor::fillGeneratorInfo
  (EventInfo_t& info, simb::MCTruth const& truth, std::vector<simb::MCParticle> const& particles) const
{
  
  std::cout << "Cosmic fillGeneratorInfo"<< std::endl;
  if (truth.NeutrinoSet()){ fillGeneratorNeutrinoInfo(info, truth);
  } else fillGeneratorCosmicInfo(info, truth, particles);
  
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
// ------------------- UPDATED for cosmics -------------------------------------
void icarus::trigger::details::EventInfoExtractor::fillGeneratorCosmicInfo
  (EventInfo_t& info, simb::MCTruth const& truth, std::vector<simb::MCParticle> const& particles) const
{
  std::cout << "Cosmic fillGeneratorCosmicInfo"<< std::endl;
  //if (truth.NeutrinoSet()) return;
  
 // simulation_time const interactionTime = getInteractionTime(truth);
  
 // if ((info.nVertices() == 0) || (interactionTime < info.InteractionTime()))
    setMainGeneratorCosmicInfo(info, truth, particles);
  //else
   // addGeneratorCosmicInfo(info, truth);
  
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
  
  //get the outgoing lepton angle w.r.t. the incoming neutrino:
  info.SetLeptonAngle(double{ truth.GetNeutrino().Theta() });

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
  if (tpc){ info.SetInActiveVolume();
  std::cout<<"TPC min X: "<<tpc->MinX()<<std::endl;
  std::cout<<"TPC max X: "<<tpc->MaxX()<<std::endl;
  std::cout<<"TPC min Y: "<<tpc->MinY()<<std::endl;
  std::cout<<"TPC max Y: "<<tpc->MaxY()<<std::endl;
  std::cout<<"TPC min Z: "<<tpc->MinZ()<<std::endl;
  std::cout<<"TPC max Z: "<<tpc->MaxZ()<<std::endl;
  }
  
} // icarus::trigger::details::EventInfoExtractor::setMainGeneratorNeutrinoInfo()

// -----------------------------------------------------------------------------
// ------------------- UPDATED for cosmics -------------------------------------
void icarus::trigger::details::EventInfoExtractor::setMainGeneratorCosmicInfo
  (EventInfo_t& info, simb::MCTruth const& truth, std::vector<simb::MCParticle> const& particles) const
{
 // std::cout << "Cosmic setMainGeneratorCosmicInfo"<< std::endl;
  /*
   * Sets the full information of the event, overwriting everything.
   * 
   * Except that the vertex is just inserted into the vertex list, as the
   * first entry.
   */
  
  
  //
  // interaction flavor (nu_mu, nu_e)
  // interaction type (CC, NC)
  //

  for (simb::MCParticle const& particle: particles){
  //const simb::MCParticle part = truth.GetParticle(i);
  //std::cout << "Cosmic truth info for particle" << i << std::endl;
  //std::cout << part.Momentum().Print() << std::endl;
   //part.Momentum().Print();
  //std::cout << "PDG info for particle" << i << std::endl;
   //std::cout << part.PdgCode() << std::endl;
  //std::cout << "Position info for particle" << i << std::endl;
   //std::cout << part.Position().X() <<std::endl;
  //std::cout << "Position info for particle" << i << std::endl;
   //part.Position().Print();
  const TLorentzVector pos1 = particle.Position();
  const TLorentzVector pos2 = particle.EndPosition();
  bool check_0 = interceptsTPCActiveVolume(pos1, pos2);
  if(check_0){ std::cout<<"Found cosmic in the detector"<<std::endl; getchar();}
  }


//bool a =  interceptsTPCActiveVolume();
//if(interceptsTPCActiveVolume()) std::cout<< "The track intercepts the Active Volume of the TPC" <<std::endl;
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
void icarus::trigger::details::EventInfoExtractor::addGeneratorCosmicInfo
  (EventInfo_t& info, simb::MCTruth const& truth) const
{
  std::cout << "Cosmic addGeneratorCosmicInfo" << std::endl;
  /*
   * The only update to the record so far is the addition of the vertex of the
   * interaction at the end of the current list.
   * No information is overwritten.
   * 
   */
  for (int i=0; i<truth.NParticles(); i++){
  const simb::MCParticle part = truth.GetParticle(i);
 // std::cout << "Cosmic truth info for particle" << i << std::endl;
  //part.Momentum().Print();// << std::endl;
  }
  //simb::MCParticle const& nu = truth.GetNeutrino().Nu();
  
  // we do not trust the vertex (`GvX()`) of the neutrino particle,
  // since GenieHelper does not translate the vertex
  // of some of the particles from GENIE to detector frame;
  // trajectory is always translated:
 // geo::Point_t const vertex { nu.EndX(), nu.EndY(), nu.EndZ() };
 // info.AddVertex(vertex);
  
} // icarus::trigger::details::EventInfoExtractor::addGeneratorNeutrinoInfo()


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
  GeV totalEnergy { 0.0 }, inSpillEnergy { 0.0 }, inPreSpillEnergy { 0.0 };
  GeV activeEnergy { 0.0 }, inSpillActiveEnergy { 0.0 },
    inPreSpillActiveEnergy { 0.0 };
  
  for (sim::SimEnergyDeposit const& edep: energyDeposits) {
    
    MeV const e { edep.Energy() }; // assuming it's stored in MeV
    
    detinfo::timescales::simulation_time const t { edep.Time() };
    bool const inSpill
      = (t >= fInSpillTimes.first) && (t <= fInSpillTimes.second);
    bool const inPreSpill
      = (t >= fInPreSpillTimes.first) && (t <= fInPreSpillTimes.second);
    
    totalEnergy += e;
    if (inSpill) inSpillEnergy += e;
    if (inPreSpill) inPreSpillEnergy += e;
    
    if (pointInActiveTPC(edep.MidPoint())) {
      activeEnergy += e;
      if (inSpill) inSpillActiveEnergy += e;
      if (inPreSpill) inPreSpillActiveEnergy += e;
    }
    
  } // for all energy deposits in the data product
  
  info.SetDepositedEnergy
    (info.DepositedEnergy() + totalEnergy);
  info.SetDepositedEnergyInSpill
    (info.DepositedEnergyInSpill() + inSpillEnergy);
  info.SetDepositedEnergyInPreSpill
    (info.DepositedEnergyInPreSpill() + inPreSpillEnergy);
  info.SetDepositedEnergyInActiveVolume
    (info.DepositedEnergyInActiveVolume() + activeEnergy);
  info.SetDepositedEnergyInSpillInActiveVolume
    (info.DepositedEnergyInSpillInActiveVolume() + inSpillActiveEnergy);
  info.SetDepositedEnergyInPreSpillInActiveVolume
    (info.DepositedEnergyInPreSpillInActiveVolume() + inPreSpillActiveEnergy);
  
} // icarus::trigger::details::EventInfoExtractor::addEnergyDepositionInfo()

// -----------------------------------------------------------------------------
bool icarus::trigger::details::EventInfoExtractor::interceptsTPCActiveVolume
  (const TLorentzVector & starPos, const TLorentzVector & endPos) const
{

double tpcMinX = -368.49; double tpcMaxX = 368.49;
double tpcMinY = -181.86; double tpcMaxY = 134.96;
double tpcMinZ = -894.951; double tpcMaxZ = 894.951;


if (starPos.X() < tpcMinX && endPos.X() < tpcMinX) return false;
if (starPos.X() > tpcMaxX && endPos.X() > tpcMaxX) return false;
if (starPos.Y() < tpcMinY && endPos.Y() < tpcMinY) return false;
if (starPos.Y() > tpcMaxY && endPos.Y() > tpcMaxY) return false;
if (starPos.Z() < tpcMinZ && endPos.Z() < tpcMinZ) return false;
if (starPos.Z() > tpcMaxZ && endPos.Z() > tpcMaxZ) return false;
if (endPos.X() > tpcMinX && endPos.X() < tpcMaxX &&
    endPos.Y() > tpcMinY && endPos.Y() < tpcMaxY &&
    endPos.Z() > tpcMinZ && endPos.Z() < tpcMaxZ) 
    {//Hit = endPos; 
    return true;}


//slope and intercept calculation

double mx = (starPos.Y() - endPos.Y()) / (starPos.X() - endPos.X());
double cx = endPos.Y() - (mx * endPos.X());

double mz = (starPos.Y() - endPos.Y()) / (starPos.Z() - endPos.Z());
double cz = endPos.Y() - (mz * endPos.Z());



double zmu_ybot = (tpcMinY - cz) / mz;
double zmu_ytop = (tpcMaxY - cz) / mz;

double xmu_ybot = (tpcMinY - cx) / mx;
double xmu_ytop = (tpcMaxY - cx) / mx;

if((zmu_ytop <= tpcMaxZ && zmu_ytop >= tpcMinZ) && ( xmu_ytop >= tpcMinX && xmu_ytop <= tpcMaxX)){

std::cout << "Save this cosmic" << std::endl;
return true;

}else if((zmu_ybot <= tpcMaxZ && zmu_ybot >= tpcMinZ) && ( xmu_ybot >= tpcMinX && xmu_ybot <= tpcMaxX)){

std::cout << "Save this cosmic" << std::endl;
return true;

} else { 

std::cout << "Don't save this cosmic" << std::endl;
return false;
}

 // geo::TPCGeo const* tpc = fGeom.PositionToTPC(particle.);

  //geo::Point_t const starPos { particle.Position().X(), particle.Position().Y(), particle.Position().Z() };
  //geo::Point_t const endPos { particle.EndX(), particle.EndY(), particle.EndZ() };
  
return false;

}

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
  EventInfoExtractor::EdepTags_t edepTags,
  geo::GeometryCore const& geom,
  std::string logCategory
  )
  : fGeneratorTags(std::move(truthTags))
  , fEnergyDepositTags(std::move(edepTags))
  , fLogCategory(std::move(logCategory))
  , fGeom(geom)
  {}


// -----------------------------------------------------------------------------
icarus::trigger::details::EventInfoExtractor
icarus::trigger::details::EventInfoExtractorMaker::make
  (TimeSpan_t inSpillTimes, TimeSpan_t inPreSpillTimes) const
{
  return EventInfoExtractor{
    fGeneratorTags, fEnergyDepositTags,
    inSpillTimes, inPreSpillTimes,
    fGeom, fLogCategory
    };
} // icarus::trigger::details::EventInfoExtractorMaker::make()


// -----------------------------------------------------------------------------
