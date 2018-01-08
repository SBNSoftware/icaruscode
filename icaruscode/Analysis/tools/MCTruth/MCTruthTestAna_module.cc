/////////////////////////////////////////////////////////////////////////////
/// Class:       MCTruthTestAna
/// Module Type: analyzer
/// File:        MCTruthTestAna_module.cc
///
/// Author:         Tracy Usher
/// E-mail address: usher@slac.stanford.edu
///
/// This is a test module which will compare the output of the parallel
/// "backtracker" based on Hit <--> MCParticle assocations to the output
/// of the actual BackTracker.
///
/// It assumes that the event contains both the simchannel information
/// needed by the BackTracker and the Hit <--> MCParticle associations
/// needed by the alternate backtracker.
///
///
/////////////////////////////////////////////////////////////////////////////

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Utilities/make_tool.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib/exception.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/PtrVector.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/SpacePoint.h"

#include "cetlib/cpu_timer.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

// The stuff we really need
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "icaruscode/Analysis/tools/MCTruth/IMCTruthMatching.h"

// ROOT includes. Note: To look up the properties of the ROOT classes,
// use the ROOT web site; e.g.,
// <http://root.cern.ch/root/html532/ClassIndex.html>
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"

// C++ Includes
#include <map>
#include <vector>
#include <tuple>
#include <algorithm>
#include <iostream>
#include <string>
#include <cmath>

#include <iostream>
#include <fstream>

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// class definition

class MCTruthTestAna : public art::EDAnalyzer
{
public:

    // Standard constructor and destructor for an ART module.
    explicit MCTruthTestAna(fhicl::ParameterSet const& pset);
    virtual ~MCTruthTestAna();

    // This method is called once, at the start of the job. In this
    // example, it will define the histograms and n-tuples we'll write.
    void beginJob();
    void endJob();

    // This method is called once, at the start of each run. It's a
    // good place to read databases or files that may have
    // run-dependent information.
    void beginRun(const art::Run& run);

    // This method reads in any parameters from the .fcl files. This
    // method is called 'reconfigure' because it might be called in the
    // middle of a job; e.g., if the user changes parameter values in an
    // interactive event display.
    void reconfigure(fhicl::ParameterSet const& pset);

    // The analysis routine, called once per event.
    void analyze (const art::Event& evt);

private:
    
    // Where the data comes from
    art::InputTag                            fHitProducerLabel;
    art::InputTag                            fTrackProducerLabel;
    
    // For keeping track of the replacement backtracker
    std::unique_ptr<truth::IMCTruthMatching> fMCTruthMatching;
    
    // Some histograms to keep a check on things
    TH1D*                                    fNumIDEHist;
    TH1D*                                    fDeltaIDEHist;
    TH1D*                                    fNegEnergyHist;
    TH1D*                                    fNegNElecHist;
    TH1D*                                    fIDEMisMatchHist;
    TH1D*                                    fBTEfficiencyHist;
    TH1D*                                    fAssocEfficiencyHist;
    TH2D*                                    fEfficiencyCompHist;
    TH1D*                                    fBTPurityHist;
    TH1D*                                    fAssocPurityHist;
    TH2D*                                    fPurityCompHist;

    // Other variables that will be shared between different methods.
    const geo::GeometryCore*                 fGeometry;       // pointer to Geometry service
    const detinfo::DetectorProperties*       fDetectorProperties;
}; // class MCTruthTestAna


//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// class implementation

//-----------------------------------------------------------------------
// Constructor
MCTruthTestAna::MCTruthTestAna(fhicl::ParameterSet const& parameterSet)
    : EDAnalyzer(parameterSet)

{
    fGeometry           = lar::providerFrom<geo::Geometry>();
    fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();

    // Read in the parameters from the .fcl file.
    this->reconfigure(parameterSet);
}

//-----------------------------------------------------------------------
// Destructor
MCTruthTestAna::~MCTruthTestAna() {}

//-----------------------------------------------------------------------
void MCTruthTestAna::beginJob()
{
    // Access ART's TFileService, which will handle creating and writing
    // histograms and n-tuples for us.
    art::ServiceHandle<art::TFileService> tfs;
    
    // The arguments to 'make<whatever>' are the same as those passed
    // to the 'whatever' constructor, provided 'whatever' is a ROOT
    // class that TFileService recognizes.
    
    // Define the histograms. Putting semi-colons around the title
    // causes it to be displayed as the x-axis label if the histogram
    // is drawn.
    fNumIDEHist          = tfs->make<TH1D>("NumIDEs",      ";# ides",                   20,   0.,  20.);
    fDeltaIDEHist        = tfs->make<TH1D>("DeltaIDE",     ";delta",                    20, -10.,  10.);
    fNegEnergyHist       = tfs->make<TH1D>("NegEnergy",    ";energy",                  500,   0.,  10.);
    fNegNElecHist        = tfs->make<TH1D>("NegNumElec",   ";# electrons",             200,   0., 100.);
    fIDEMisMatchHist     = tfs->make<TH1D>("IDEMisMatch",  ";# mismatches",             20,  10.,  10.);
    fBTEfficiencyHist    = tfs->make<TH1D>("BTEfficiency", ";efficiency",              102,   0.,   1.02);
    fAssocEfficiencyHist = tfs->make<TH1D>("AssocEffic",   ";efficiency",              102,   0.,   1.02);
    fEfficiencyCompHist  = tfs->make<TH2D>("Efficiency",   ";BackTrack;Associations",   52,   0.,   1.04, 52, 0., 1.04);
    fBTPurityHist        = tfs->make<TH1D>("BTPurity",     ";Purity",                  102,   0.,   1.02);
    fAssocPurityHist     = tfs->make<TH1D>("AssocPurity",  ";Purity",                  102,   0.,   1.02);
    fPurityCompHist      = tfs->make<TH2D>("Purity",       ";BackTrack;Associations",   52,   0.,   1.04, 52, 0., 1.04);

}

//-----------------------------------------------------------------------
void MCTruthTestAna::beginRun(const art::Run& /*run*/) {}

//-----------------------------------------------------------------------
void MCTruthTestAna::reconfigure(fhicl::ParameterSet const& pset)
{
    // Read parameters from the .fcl file. The names in the arguments
    // to p.get<TYPE> must match names in the .fcl file.
    fHitProducerLabel              = pset.get<art::InputTag>("HitModuleLabel");
    fTrackProducerLabel            = pset.get<art::InputTag>("TrackProducerLabel");

    // Get the tool for MC Truth matching
    fMCTruthMatching = art::make_tool<truth::IMCTruthMatching>(pset.get<fhicl::ParameterSet>("MCTruthMatching"));

    return;
}

//-----------------------------------------------------------------------
void MCTruthTestAna::analyze(const art::Event& event)
{
    // Recover the backtracker
    art::ServiceHandle<cheat::BackTrackerService>       backTracker;
    art::ServiceHandle<cheat::ParticleInventoryService> partInventory;
   
    // "Rebuild" the maps used by the parallel backtracker
    fMCTruthMatching->Rebuild(event);
    
    // Begin the comparisons by simply checking that the number of MCParticles agree between the two
    const sim::ParticleList& particleList = fMCTruthMatching->ParticleList();
    const sim::ParticleList& btPartList   = partInventory->ParticleList();
    
    if (particleList.size() != btPartList.size())
    {
        mf::LogInfo("MCTruthTestAna") << "*** ParticleList size mismatch, BackTracker says: " << btPartList.size() << ", alternate says: " << particleList.size() << std::endl;
    }

    // Flush with the success of passing that test, let's recover the hits and then loop through them
    // and do hit-by-hit comparisons...
    art::Handle<std::vector<recob::Hit> > hitHandle;
    event.getByLabel(fHitProducerLabel, hitHandle);

    if (hitHandle.isValid())
    {
        int nTotalMisMatches(0);
        int nNegTrackIDs(0);
        
        for(const auto& hit : *hitHandle)
        {
            int nIDEMisMatches(0);
            
            // Check the claimed parentage of the current hit
            std::vector<sim::TrackIDE> trackIDEVec = fMCTruthMatching->HitToTrackID(hit);
            std::vector<sim::TrackIDE> btTrkIDEVec = backTracker->HitToTrackIDEs(hit);
            
            int deltaIDs = int(btTrkIDEVec.size()) - int(trackIDEVec.size());
            
            fNumIDEHist->Fill(trackIDEVec.size(),1.);
            fDeltaIDEHist->Fill(deltaIDs, 1.);

            // It can be that the BackTracker is returning IDE's with negative Track ID's and these don't
            // appear in the associations, so bypass them here...
            if (deltaIDs != 0)
            {
                int numNeg = std::count_if(btTrkIDEVec.begin(),btTrkIDEVec.end(),[](const auto& trackIDE){return trackIDE.trackID < 0;});
                
                deltaIDs -= numNeg;
                
                std::sort(btTrkIDEVec.begin(),btTrkIDEVec.end(),[](const auto& left, const auto& right){return left.trackID > right.trackID;});
                
                while(btTrkIDEVec.back().trackID < 0)
                {
                    fNegEnergyHist->Fill(std::min(9.98,double(btTrkIDEVec.back().energy)), 1.);
                    fNegNElecHist->Fill(std::min(99.9,double(btTrkIDEVec.back().numElectrons)), 1.);
                    
                    btTrkIDEVec.pop_back();
                    nNegTrackIDs++;
                }
            }
            
            // Currently assume BackTracker is "correct" so loop over its output and look for match
            for(const auto& btTrkID : btTrkIDEVec)
            {
                std::vector<sim::TrackIDE>::const_iterator trackIDItr = std::find_if(trackIDEVec.begin(),trackIDEVec.end(),[btTrkID](const auto& id){return btTrkID.trackID == id.trackID && btTrkID.energy == id.energy && btTrkID.numElectrons == id.numElectrons;});
                
                if (trackIDItr == trackIDEVec.end()) nIDEMisMatches++;
            }
            
            fIDEMisMatchHist->Fill(nIDEMisMatches, 1.);
            nTotalMisMatches += nIDEMisMatches;
        }
        
        mf::LogInfo("MCTruthTestAna") << "==> Found " << nTotalMisMatches << " between BackTracker and Associations, BT reported " << nNegTrackIDs << " negative track ids \n";
    }
    
    // We're on a roll now!
    // So let's see if we can recover information about tracks
    art::Handle<std::vector<recob::Track>> trackHandle;
    event.getByLabel(fTrackProducerLabel, trackHandle);
    
    if (trackHandle.isValid())
    {
        // We are going to want a vector of art pointers to hits, so build here
        std::vector<art::Ptr<recob::Hit>> hitPtrVector;
        
        art::fill_ptr_vector(hitPtrVector, hitHandle);

        // Recover the associations between tracks and hits
        art::FindManyP<recob::Hit> hitTrackAssns(trackHandle, event, fTrackProducerLabel);
        
        int nBadEffMatches(0);
        int nBadPurMatches(0);

        for(size_t trackIdx = 0; trackIdx < trackHandle->size(); trackIdx++)
        {
            art::Ptr<recob::Track> track(trackHandle,trackIdx);
            
            std::vector<art::Ptr<recob::Hit>> trackHitVec = hitTrackAssns.at(track.key());
            
            std::set<int> trackIDSet = fMCTruthMatching->GetSetOfTrackIDs(trackHitVec);
            std::set<int> btTrkIDSet = backTracker->GetSetOfTrackIds(trackHitVec);
            
            // Check for case were we might have negative track IDs from the BackTracker
            if (btTrkIDSet.size() > trackIDSet.size())
            {
                // Note that the container here is a std::set which will be ordered smallest to largest.
                // In this case we want to remove elements from the front until they are positive
                for(std::set<int>::iterator btTrkIDSetItr = btTrkIDSet.begin(); btTrkIDSetItr != btTrkIDSet.end();)
                {
                    if (!(*btTrkIDSetItr < 0)) break;
                    
                    btTrkIDSetItr = btTrkIDSet.erase(btTrkIDSetItr);
                }
            }
            
            if (btTrkIDSet.size() != trackIDSet.size())
                mf::LogDebug("MCTruthTestAna") << "Mismatch in associated track ids, backtracker: " << btTrkIDSet.size() << ", associations: " << trackIDSet.size() << ", track id: " << track.key() << "\n";

            // Ok, turn the "set" of track ID's into a vector of track ID's
            std::vector<int> btTrkIDVec;
            
            for(auto& trackID : btTrkIDSet) btTrkIDVec.emplace_back(trackID);
            
            // And then use this to recover the vectors of hits associated to each MC track
            std::vector<std::vector<art::Ptr<recob::Hit>>> trkHitVecVec;

            for(const auto& tkID : btTrkIDVec)
            {
                std::vector<art::Ptr<recob::Hit>> hitVec = backTracker->TrackIdToHits_Ps(tkID, hitPtrVector);
                trkHitVecVec.push_back(hitVec);
            }
            // Apply majority logic - we declare the MCParticle with the most hits to be the "winner"
            std::vector<std::vector<art::Ptr<recob::Hit>>>::iterator bestTrkHitVecItr = std::max_element(trkHitVecVec.begin(),trkHitVecVec.end(),[](const auto& a, const auto& b){return a.size() < b.size();});
            
            if (bestTrkHitVecItr == trkHitVecVec.end())
            {
                // I have been watching way too much Star Trek (the original!)
                mf::LogDebug("MCTruthTestAna") << ">>>>>>> ERROR! >>>>>>> ERROR! >>>>>> MUST PURIFY! >>>>>> ERROR!" << "\n";
                continue;
            }
            
            int indexToBestMC = std::distance(trkHitVecVec.begin(),bestTrkHitVecItr);
            int bestMCTrackID = btTrkIDVec.at(indexToBestMC);
            
            std::vector<art::Ptr<recob::Hit>>& bestMCTrackHitVec = *bestTrkHitVecItr;
            
            std::set<int> mcTrackIdxSet = {bestMCTrackID};
            
            double btTrkEffic = backTracker->HitCollectionEfficiency(mcTrackIdxSet, trackHitVec, bestMCTrackHitVec, geo::k3D);
            double trackEffic = fMCTruthMatching->HitCollectionEfficiency(mcTrackIdxSet, trackHitVec, bestMCTrackHitVec, geo::k3D);
            
            if (btTrkEffic != trackEffic)
            {
                mf::LogDebug("MCTruthTestAna") << "Efficiency mismatch, track ID: " << bestMCTrackID << ", track: " << track.key() << ", # hits: " << trackHitVec.size() << ", # MC hits: " << bestMCTrackHitVec.size() << ", btTrkEff: " << btTrkEffic << ", trackEff: " << trackEffic << "\n|";
                nBadEffMatches++;
            }
            
            fBTEfficiencyHist->Fill(std::min(1.01,btTrkEffic), 1.);
            fAssocEfficiencyHist->Fill(std::min(1.01,trackEffic), 1.);
            fEfficiencyCompHist->Fill(std::min(1.01,btTrkEffic), std::min(1.01,trackEffic), 1.);
            
            double btTrkPurity = backTracker->HitCollectionPurity(mcTrackIdxSet, trackHitVec);
            double trackPurity = fMCTruthMatching->HitCollectionPurity(mcTrackIdxSet, trackHitVec);
            
            if (btTrkPurity != trackPurity)
            {
                mf::LogDebug("MCTruthTestAna") << "Purity mismatch, track ID: " << bestMCTrackID << ", track: " << track.key() << ", # hits: " << trackHitVec.size() << ", # MC hits: " << bestMCTrackHitVec.size() << ", btTrkPurity: " << btTrkPurity << ", trackPurity: " << trackPurity << "\n|";
                nBadPurMatches++;
            }
            
            fBTPurityHist->Fill(std::min(1.01,btTrkPurity), 1.);
            fAssocPurityHist->Fill(std::min(1.01,trackPurity), 1.);
            fPurityCompHist->Fill(std::min(1.01,btTrkPurity),std::min(1.01,trackPurity), 1.);
        }
        
        mf::LogDebug("MCTruthTestAna") << "Event with " << trackHandle->size() << " reconstructed tracks, found " << nBadEffMatches << " efficiency mismatches, " << nBadPurMatches << " purity mismatchs \n";
    }

    return;
}

void MCTruthTestAna::endJob()
{
    // Make a call to normalize histograms if so desired

    return;
}

// This macro has to be defined for this module to be invoked from a
// .fcl file; see MCTruthTestAna.fcl for more information.
DEFINE_ART_MODULE(MCTruthTestAna)

