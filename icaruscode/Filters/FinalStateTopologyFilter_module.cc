////////////////////////////////////////////////////////////////////////
//
// FilterTopology class:
// Algorith to filter interactions with specific final state topology
// author: hhausner@fnal.gov
//
//////////////////////////////////////////////////////////////////////////

// ART
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

// ROOT
#include "TH1.h"

// std inculdes
#include <memory>

namespace icarus::filter::topology
{
  class FSTopology
  {
    public:
      // reset - clear the container
      void reset() { primaries.clear(); }
      // count - how many particles of the given pdg code are in the interaction?
      size_t count(const int& pdg) const
      {
        bool inMap = (primaries.count(pdg) == 1);
        return (not inMap) ? 0 : primaries.at(pdg);
      }
      // does the interaction have pdg with mother motherPDG?
      bool particleHasMother(const int& pdg, const int& motherPDG) const
      {
        if ((count(pdg) == 0) || (mothers.count(pdg) == 0)) return false;
        return (mothers.at(pdg).count(motherPDG) == 1);
      }
      // is this interaction inclusive to the selection provided?
      bool inclusive(const std::unordered_map<int, std::pair<size_t, bool>>& request,
                     const std::unordered_map<int, std::vector<int>>& origins = {}) const
      {
        for (auto const& [pdg, countAndExclusivity] : request)
        {
          size_t have = count(pdg);
          if ((have > 0) && (origins.count(pdg) == 1))
          {
            bool hasDesiredMother = false;
            for (auto const& motherPDG : origins.at(pdg))
            {
              if (particleHasMother(pdg, motherPDG))
              {
                hasDesiredMother = true;
                break;
              }
            }
            if (not hasDesiredMother)
              return false;
          }
          if (countAndExclusivity.second)
          {
            if (have != countAndExclusivity.first)
              return false;
          }
          else
          {
            if (have < countAndExclusivity.first) // NB: this may behave strangely if counts = 0 and exclusive = false
              return false;
          }
        }
        return true;
      }
      // is the interaction _exactly_ the topology specified?
      bool match(const std::unordered_map<int, std::pair<size_t, bool>>& request) const
      {
        // must at least be inclusive...
        if (not inclusive(request))
          return false;
        //...and also not have anything extra
        for (const auto& [primaryPDG, primaryCount] : primaries)
        {
          if (request.count(primaryPDG) != 1)
            return false;
          // if primaryPDG is in the request, then b/c inclusive(request) == true
          // we know it must have the quanties we are looking for. no need to double
          // check primaryCount
        }
        // if it has all the things we want, and nothing we don't, it is a match
        return true;
      }
      // constructor - iterate through the particle list and count
      FSTopology(const simb::MCTruth& mcTruth)
      {
        int nParticles = mcTruth.NParticles();
        for (int pIdx = 0; pIdx < nParticles; ++pIdx)
        {
          const simb::MCParticle& pParticle = mcTruth.GetParticle(pIdx);
          if (pParticle.StatusCode() != 1) continue; // keep only the final state particles
          int pPDG = pParticle.PdgCode();
          if ((std::abs(pPDG) == 12) ||
              (std::abs(pPDG) == 14) ||
              (std::abs(pPDG) == 16)  ) continue; // skip neutrinos
          if ( std::abs(pPDG) >= 100000000 ) continue; // skip GENIE junk particles
          increment(pPDG);
          addMother(pPDG, pParticle.Mother());
        }
      }
      // streamer - print the primaries map
      friend std::ostream& operator<<(std::ostream& os, const FSTopology& topology)
      {
        os << "Topology: ";
        for (auto const& [pdg, count] : topology.primaries)
        {
          os << '\n' << "  PDG " << pdg << ", count " << count;
        }
        return os;
      }
    private:
      std::unordered_map<int, size_t> primaries;
      std::unordered_map<int, std::unordered_set<int>> mothers;
      void increment(const int& pdg)
      {
        bool inMap = (primaries.count(pdg) == 1);
        if (not inMap) primaries[pdg] = 0;
        primaries[pdg] += 1;
      }
      void addMother(const int& pdg, const int& pdgMother)
      {
        bool inMapMothers = (mothers.count(pdg) == 1);
        if (not inMapMothers) mothers[pdg] = {};
        mothers[pdg].insert(pdgMother);
      }
  };

  class FinalStateTopologyFilter : public art::EDFilter
  {
    public:
      explicit FinalStateTopologyFilter(fhicl::ParameterSet const&);
      virtual ~FinalStateTopologyFilter();

      void reconfigure(fhicl::ParameterSet const& pset);
      bool filter(art::Event& evt);
    private:
      std::string fGenieModuleLabel;         ///< Where to grab the interactions
      std::unordered_map<int, std::pair<size_t, bool>> fRequest; ///< What are we looking for?
      std::unordered_map<int, std::vector<int>> fOrigins; ///< Are we requesting specific origins?
      bool fInclusive;                       ///< Do we allow interactions to have more than the requested particles?
      bool fVerbose;                         ///< Print debugging statements?
      std::unique_ptr<FSTopology> fTopology; ///< Unique point to FS topology
      std::unique_ptr<TH1D> fSelectedEvents; ///< Count selected events
      std::unique_ptr<TH1D> fTotalEvents;    ///< Count total events
  }; // end FinalStateTopologyFilter class

  //-----------------------------------------------------------------------------------------------
  // FinalStateTopologyFilter constructor
  // Call the reconfigure function and set up counter histograms
  FinalStateTopologyFilter::FinalStateTopologyFilter(fhicl::ParameterSet const& pset) : EDFilter{pset}
  {
    this->reconfigure(pset);

    // set up histograms
    art::ServiceHandle<art::TFileService> tfs;
    fSelectedEvents.reset(tfs->make<TH1D>("fSelectedEvents", "Number of Selected Events", 3, 0, 3));
    fTotalEvents   .reset(tfs->make<TH1D>("fTotalEvents",    "Total Events",              3, 0, 3));
  }

  //-----------------------------------------------------------------------------------------------
  // FinalStateTopologyFilter destructor
  // Don't do nothin'
  FinalStateTopologyFilter::~FinalStateTopologyFilter()
  {
  }

  //-----------------------------------------------------------------------------------------------
  // FinalStateTopologyFilter::reconfigure
  // Read in the fhicl parameters and set up the histograms
  void FinalStateTopologyFilter::reconfigure(fhicl::ParameterSet const& pset)
  {
    fGenieModuleLabel = pset.get<std::string>("GenieModuleLabel");
    fInclusive        = pset.get<bool>       ("IsInclusive", false); // default to false
    fVerbose          = pset.get<bool>       ("IsVerbose",   false); // default to false

    std::vector<int>    PDGCodes         = pset.get<std::vector<int>>   ("PDG");
    std::vector<size_t> PDGCounts        = pset.get<std::vector<size_t>>("PDGCount");
    std::vector<bool>   CountsExclusive  = pset.get<std::vector<bool>>  ("PDGCountExclusivity");

    std::vector<std::vector<int>> PDGOrigins = pset.get<std::vector<std::vector<int>>>("PDGOrigin", {});

    // make sure our vectors are of equal size
    if ((PDGCodes.size() != PDGCounts      .size()) ||
        (PDGCodes.size() != CountsExclusive.size()) )
      throw std::invalid_argument("FinalStateTopologyFilter::reconfigure — size of input vectors must match");
    if ((PDGOrigins.size() > 0) && (PDGCodes.size() != PDGOrigins.size()))
      throw std::invalid_argument("FinalStateTopologyFilter::reconfigure — size of input vectors must match");
    // make request map, check no collisions in PDG codes
    if (fVerbose)
      mf::LogVerbatim("FinalStateTopologyFilter")
        << "Requesting Topology: ";
    for (size_t idx = 0; idx < PDGCodes.size(); ++idx)
    {
      if (fVerbose)
      {
        mf::LogVerbatim("FinalStateTopologyFilter")
         << "  -PDG  " << PDGCodes.at(idx)
         << ", counts " << PDGCounts.at(idx)
         << ", exclusive? " << std::boolalpha << CountsExclusive.at(idx);
      }
      if (not fRequest.insert({PDGCodes.at(idx), {PDGCounts.at(idx), CountsExclusive.at(idx)}}).second)
        throw std::invalid_argument("FinalStateTopologyFilter::reconfigure — cannot repeat PDG codes");
      if (PDGOrigins.size() > 0 && PDGOrigins.at(idx).size() > 0)
      {
        if (fVerbose)
        {
          std::string logStrm;
          for (auto const& origin : PDGOrigins.at(idx))
            logStrm += " " + origin;
          mf::LogVerbatim("FinalStateTopologyFilter")
            << "  Particle " << PDGCodes.at(idx) << " must come from one of" << logStrm;
        }   
        if (not fOrigins.insert({PDGCodes.at(idx), PDGOrigins.at(idx)}).second)
          throw std::invalid_argument("FinalStateTopologyFilter::reconfigure — cannot repeat PDG codes");
      }
    }
  }
  
  //-----------------------------------------------------------------------------------------------
  // FinalStateTopologyFilter::filter
  // filter events based on the supplied topology
  bool FinalStateTopologyFilter::filter(art::Event& evt)
  {
    // get a handle on our input interactions
    art::Handle<std::vector<simb::MCTruth>> mcTruthsHandle;
    evt.getByLabel(fGenieModuleLabel, mcTruthsHandle);

    // increment total counter
    fTotalEvents->Fill(1);

    // loop through the MCTruths and filter only on the beam neutrinos
    if (fVerbose)
      mf::LogVerbatim("FinalStateTopologyFilter")
        << "Have " << mcTruthsHandle->size() << " truths to check";
    for (size_t mcIdx = 0; mcIdx < mcTruthsHandle->size(); mcIdx++)
    {
      art::Ptr<simb::MCTruth> mcTruth(mcTruthsHandle, mcIdx);
      if (fVerbose)
      {
        std::string orgStr;
        switch (mcTruth->Origin())
        {
          case simb::kUnknown:
            orgStr = "kUnknown";
            break;
          case simb::kBeamNeutrino:
            orgStr = "kBeamNeutrino";
            break;
          case simb::kCosmicRay:
            orgStr = "kCosmicRay";
            break;
          case simb::kSuperNovaNeutrino:
            orgStr = "kSuperNovaNeutrino";
            break;
          case simb::kSingleParticle:
            orgStr = "kSingleParticle";
            break;
        }
        mf::LogVerbatim("FinalStateTopologyFilter")
          << "MCTruth(" << mcIdx << ") has origin " << mcTruth->Origin() << " (" << orgStr << ")";
      }
      if (mcTruth->Origin() != simb::kBeamNeutrino)
        continue;
      fTopology = std::make_unique<FSTopology>(*mcTruth); // read in the topology
      if (fVerbose)
        mf::LogVerbatim("FinalStateTopologyFilter") 
          << "Found beam neutrino event at index " << mcIdx << '\n'
          << "  " << *fTopology;
      if (fInclusive && fTopology->inclusive(fRequest))
      {
        fSelectedEvents->Fill(1);
        return true;
      }
      else if (fTopology->match(fRequest))
      {
        fSelectedEvents->Fill(1);
        return true;
      }
      // reset fTopology to save memory
      fTopology->reset();
    }

    return false;
  }

  DEFINE_ART_MODULE(FinalStateTopologyFilter)

} // end icarus::filter::topology namespace
