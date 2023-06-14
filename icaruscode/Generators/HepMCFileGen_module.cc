/**
 * @file HepMCFileGen_module.cc
 * @brief Producer generating Monte Carlo truth record in LArSoft format from a text file in HepMC format
 * @date Dec 2019
 * @author Marco Del Tutto
 * @comment adapted from TxtFileGen by Brian Rebel
 */
/**
 * @class evgen::HepMCFileGen
 *
 *  This module assumes that the input file has the hepevt format for
 *  each event to be simulated. In brief each event contains at least two
 *  lines. The first line contains two entries, the event number (which is
 *  ignored in ART/LArSoft) and the number of particles in the event. Each
 *  following line containes 15 entries to describe each particle. The entries
 *  are:
 *
 *  1.  status code (should be set to 1 for any particle to be tracked, others
 *      won't be tracked)
 *  2.  the pdg code for the particle
 *  3.  the entry of the first mother for this particle in the event,
 *      0 means no mother
 *  4.  the entry of the second mother for this particle in the event,
 *      0 means no mother
 *  5. the entry of the first daughter for this particle in the event,
 *      0 means no daughter
 *  6. the entry of the second daughter for this particle in the event,
 *      0 means no daughter
 *  7. x component of the particle momentum
 *  8. y component of the particle momentum
 *  9. z component of the particle momentum
 *  10. energy of the particle
 *  11. mass of the particle
 *  12. x position of the particle initial position
 *  13. y position of the particle initial position
 *  14. z position of the particle initial position
 *  15. time of the particle production
 *
 *  For example, if you want to simulate a single muon with a 5 GeV energy
 *  moving only in the z direction, the entry would be (see onemuon.hepmc):
 *  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  0 1
 *  1 13 0 0 0 0 0. 0. 1.0 5.0011 0.105 1.0 1.0 1.0 0.0
 *  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  Or if you want to simulate a muon neutrino event (see oneneutrino.hepmc): 
 *  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  0 3
 *  0 14 0 0 0 0 0.00350383 0.002469 0.589751 0.589766 0 208.939 63.9671 10.9272 4026.32
 *  1 13 1 0 0 0 -0.168856 -0.0498011 0.44465 0.489765 105.658 208.939 63.9671 10.9272 4026.32
 *  1 2212 1 0 0 0 0.151902 -0.124578 0.0497377 0.959907 938.272 208.939 63.9671 10.9272 4026.32
 *  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  here the first particle is the initial state neutrino (status code 0, meaning initial state, 
 *  not to be prooagated by GEANT). The other particles are the initial state particles.
 *  For the neutrino, write ccnc in place of 1st daugther, and mode (qe, res, ...) 
 *  in place of 2nd daugther.
 *
 *  There are some assumptions that go into using this format that may not
 *  be obvious.  The first is that only particles with status code = 1
 *  are tracked in the LArSoft/Geant4 combination making the mother daughter
 *  relations somewhat irrelevant. That also means that you should let
 *  Geant4 handle any decays.
 *
 *  The units in LArSoft are cm for distances and ns for time.
 *  The use of `TLorentzVector` below does not imply space and time have the same units
 *   (do not use `TLorentzVector::Boost()`).
 */
#include <string>
#include <fstream>
#include <math.h>
#include <map>
#include <iomanip>
#include <algorithm>
#include <sstream>
#include <glob.h>
#include <cstdlib>  // for unsetenv()
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib_except/exception.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "cetlib/search_path.h"
#include "cetlib/getenv.h"
#include "cetlib/split_path.h"
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "TLorentzVector.h"
#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SummaryData/RunData.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "ifdh.h"
#include "ifdh_art/IFDHService/IFDH_service.h" // ifdh_ns::IFDH
#undef USE_IFDH_SERVICE // ifdh for now

namespace evgen {
  class HepMCFileGen;
}
class evgen::HepMCFileGen : public art::EDProducer {
public:
  explicit HepMCFileGen(fhicl::ParameterSet const & p);
  void produce(art::Event & e)        override;
  void beginJob()                                 override;
  void beginRun(art::Run & run)                   override;
  void endSubRun(art::SubRun& sr)     override;
private:
  std::ifstream open_file();
  std::string fInputFilePath; ///< Path to the HEPMC input file, relative to `FW_SEARCH_PATH`.
  std::ifstream* fInputFile;
  
  double         fEventsPerPOT;     ///< Number of events per POT (to be set)
  int            fEventsPerSubRun;  ///< Keeps track of the number of processed events per subrun
  // ifdh_ns::ifdh* fIFDH;             ///< (optional) flux file handling
};
//------------------------------------------------------------------------------
evgen::HepMCFileGen::HepMCFileGen(fhicl::ParameterSet const & p)
  : EDProducer{p}
  , fInputFilePath(p.get<std::string>("InputFilePath"))
  , fInputFile(nullptr)
  , fEventsPerPOT{p.get<double>("EventsPerPOT", -1.)}
  , fEventsPerSubRun(0)
{
  produces< std::vector<simb::MCTruth>   >();
  produces< sumdata::RunData, art::InRun >();
  produces< sumdata::POTSummary, art::InSubRun >();
}
//------------------------------------------------------------------------------

std::ifstream evgen::HepMCFileGen::open_file()
{
  /*
   * The plan:
   *  1. expand the path in FW_SEARCH_PATH (only if relative path)
   *  2. copy it into scratch area (only if starts with `/pnfs`)
   *  3. open the file (original or copy) and return the opened file
   * 
   * Throws a cet::exception if eventually file is not found.
   */
  
  std::string fullFileName = fInputFilePath;
  
  cet::search_path searchPath("FW_SEARCH_PATH");
  if (searchPath.find_file(fInputFilePath, fullFileName)) {
    mf::LogDebug("HepMCFileGen")
      << "Input file '" << fInputFilePath << "' found in FW_SEARCH_PATH:\n"
      << fullFileName
      ;
  }
  
  //
  // prepare the file with IFDH, if path starts with `/pnfs`
  //
  if (fullFileName.compare(0, 6, "/pnfs/") == 0) { 
    fullFileName = art::ServiceHandle<IFDH>()->fetchInput(fullFileName);
    mf::LogDebug("HepMCFileGen")
      << "IFDH fetch: '" << fInputFilePath << "' -> '" << fullFileName << "'";
  }
  
  //
  // attempt to open
  //
  mf::LogDebug("HepMCFileGen")
    << "Reading input file '" << fInputFilePath << "' as:\n" << fullFileName;
  std::ifstream inputFile(fullFileName);
  if (inputFile) return inputFile;
  
  // all attempts failed, give up:
  throw cet::exception("HepMCFileGen")
    << "HEPMC input file '" << fInputFilePath << "' can't be opened.\n";
  
} // evgen::HepMCFileGen::open_file()



//------------------------------------------------------------------------------
void evgen::HepMCFileGen::beginJob()
{
  fInputFile = new std::ifstream(open_file());
  
}
//------------------------------------------------------------------------------
void evgen::HepMCFileGen::beginRun(art::Run& run)
{
    fEventsPerSubRun = 0;
    art::ServiceHandle<geo::Geometry const> geo;
    run.put(std::make_unique<sumdata::RunData>(geo->DetectorName()), art::fullRun());
  }
//------------------------------------------------------------------------------
void evgen::HepMCFileGen::endSubRun(art::SubRun& sr)
  {
    auto p = std::make_unique<sumdata::POTSummary>();
    p->totpot     = fEventsPerSubRun * fEventsPerPOT;
    p->totgoodpot = fEventsPerSubRun * fEventsPerPOT;
    sr.put(std::move(p), art::subRunFragment());
    return;
  }
//------------------------------------------------------------------------------
void evgen::HepMCFileGen::produce(art::Event & e)
{
  // check that the file is still good
  if( !fInputFile->good() || fInputFile->peek() == EOF) {
    open_file();
  }
  if( !fInputFile->good() || fInputFile->peek() == EOF) {
    throw cet::exception("HepMCFileGen") << "input text file "
                                        << " cannot be read in produce().\n";
  }
  std::unique_ptr< std::vector<simb::MCTruth> > truthcol(new std::vector<simb::MCTruth>);
  simb::MCTruth truth;
  // declare the variables for reading in the event record
  int            event          = 0;
  unsigned short nParticles           = 0;
  int            status         = 0;
  int                   pdg            = 0;
  int                   firstMother    = 0;
  int                   secondMother   = 0;
  int                   firstDaughter  = 0;
  int                   secondDaughter = 0;
  double          xMomentum      = 0.;
  double          yMomentum           = 0.;
  double          zMomentum           = 0.;
  double          energy              = 0.;
  double          mass                = 0.;
  double          xPosition           = 0.;
  double          yPosition           = 0.;
  double          zPosition           = 0.;
  double          time                = 0.;
  bool set_neutrino = false;
  // neutrino
  int ccnc = -1, mode = -1, itype = -1, target = -1, nucleon = -1, quark = -1;
  double w = -1, x = -1, y = -1, qsqr = -1;
  // read in line to get event number and number of particles
  std::string oneLine;
  std::getline(*fInputFile, oneLine);
  std::istringstream inputLine;
  inputLine.str(oneLine);
  inputLine >> event >> nParticles;
  // now read in all the lines for the particles
  // in this interaction. only particles with
  // status = 1 get tracked in Geant4. see GENIE GHepStatus
  for(unsigned short i = 0; i < nParticles; ++i){
    std::getline(*fInputFile, oneLine);
    inputLine.clear();
    inputLine.str(oneLine);
    inputLine >> status >> pdg
              >> firstMother >> secondMother >> firstDaughter >> secondDaughter
              >> xMomentum   >> yMomentum    >> zMomentum     >> energy >> mass
              >> xPosition   >> yPosition    >> zPosition     >> time;
    std::cout<<"pdg_all"<<pdg<<"xmomentum_all"<<xMomentum<<"ymomentum_all"<<yMomentum<<"zposition_all"<<zPosition<<std::endl;
    TLorentzVector pos(xPosition, yPosition, zPosition, time);
    TLorentzVector mom(xMomentum, yMomentum, zMomentum, energy);
    simb::MCParticle part(i, pdg, "primary", firstMother, mass, status);
    part.AddTrajectoryPoint(pos, mom);
    //if (abs(pdg) == 18 || abs(pdg) == 12) 
    if (abs(pdg) == 52 )  // Animesh made changes
	{
      set_neutrino = true;
      ccnc = firstDaughter; // for the neutrino we write ccnc in place of 1st daugther
      mode = secondDaughter; // for the neutrino we write mode in place of 2nd daugther
      itype = -1;
      target = nucleon = quark = w = x = y = qsqr = -1;
     std::cout<<"pdg_mother"<<pdg<<"xmomentum_mother"<<xMomentum<<"ymomentum_mother"<<yMomentum<<"zposition_mother"<<zPosition<<std::endl;
    } 
    truth.Add(part);
     std::cout << i << "  Particle added with Pdg " << part.PdgCode() << ", Mother " << part.Mother() << ", track id " << part.TrackId() << ", ene " << part.E() << std::endl;
  }
 
  if (set_neutrino) {
    truth.SetNeutrino(ccnc,
                      mode,
                      itype,
                      target,
                      nucleon,
                      quark,
                      w,
                      x,
                      y,
                      qsqr);
    // set the neutrino information in MCTruth
    truth.SetOrigin(simb::kBeamNeutrino);
    // truth.SetGeneratorInfo(simb::Generator_t::kGENIE,
    //                          __GENIE_RELEASE__,
    //                          {{"tune", fTuneName}});
  }
  //std::cout << " neutrino " << truth.GetNeutrino() << std::endl;
  //std::cout << " lepton pdg " << truth.GetNeutrino().Lepton().PdgCode() << " ene " << truth.GetNeutrino().Lepton().E() << std::endl;
  truthcol->push_back(truth);
  e.put(std::move(truthcol));
  fEventsPerSubRun++;
  return;
}
DEFINE_ART_MODULE(evgen::HepMCFileGen)
