#ifndef CRTBACKTRACKER_H_SEEN

#define CRTBACKTRACKER_H_SEEN

/////////////////////////////////////////////////////////////////
// CRTBackTracker.h
//
// Quick and dirty backtracker for SBND CRT
// written by T Brooks (tbrooks@fnal.gov), November 2018
//
// Ported into icaruscode for use with ICARUS CRT
//  by C Hilgenberg (Chris.Hilgenberg@colostate.edu), April 2020
/////////////////////////////////////////////////////////////////

// framework
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "canvas/Persistency/Common/FindManyP.h"

// Utility libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "icaruscode/CRT/CRTProducts/CRTData.hh"
#include "icaruscode/CRT/CRTProducts/CRTHit.hh"
#include "icaruscode/CRT/CRTProducts/CRTTrack.hh"

// c++
#include <vector>

namespace icarus{
 namespace crt {
    class CRTBackTracker;
 }
}

class icarus::crt::CRTBackTracker {

  public:

    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<art::InputTag> CRTTrueHitLabel {
        Name("CRTTrueHitLabel")
      };
      fhicl::Atom<art::InputTag> CRTDataLabel {
        Name("CRTDataLabel")
      };

      fhicl::Atom<art::InputTag> CRTSimHitLabel {
        Name("CRTSimHitLabel")
      };

      fhicl::Atom<art::InputTag> CRTTrackLabel {
        Name("CRTTrackLabel")
      };

      fhicl::Atom<bool> RollupUnsavedIds {
        Name("RollupUnsavedIds")
      };
    };

    CRTBackTracker(const Config& config);
    CRTBackTracker(const fhicl::ParameterSet& pset) :
    CRTBackTracker(fhicl::Table<Config>(pset, {})()) {}
    CRTBackTracker();
    ~CRTBackTracker();

    void reconfigure(const Config& config);

    // Initialize to speed things up
    void Initialize(const art::Event& event);

    // Check that two CRT data products are the same
    bool DataCompare(const CRTData& data1, const CRTData& data2);

    // Check that two CRT hits are the same
    bool HitCompare(const CRTHit& hit1, const CRTHit& hit2);

    // Check that two CRT tracks are the same
    bool TrackCompare(const CRTTrack& track1, const CRTTrack& track2);

    // Get all the true particle IDs that contributed to the CRT data product
    std::vector<int> AllTrueIds(const art::Event& event, const CRTData& data);

    // Get all the true particle IDs that contributed to the CRT hit
    std::vector<int> AllTrueIds(const art::Event& event, const CRTHit& hit);

    // Get all the true particle IDs that contributed to the CRT track
    std::vector<int> AllTrueIds(const art::Event& event, const CRTTrack& track);

    // Get the true particle ID that contributed the most energy to the CRT data product
    int TrueIdFromTotalEnergy(const art::Event& event, const CRTData& data);

    // Faster function - needs Initialize() to be called first
    int TrueIdFromDataId(const art::Event& event, int data_i);

    // Get the true particle ID that contributed the most energy to the CRT hit
    int TrueIdFromTotalEnergy(const art::Event& event, const CRTHit& hit);

    // Faster function - needs Initialize() to be called first
    int TrueIdFromHitId(const art::Event& event, int hit_i);

    // Get the true particle ID that contributed the most energy to the CRT track
    int TrueIdFromTotalEnergy(const art::Event& event, const CRTTrack& track);

    // Faster function - needs Initialize() to be called first
    int TrueIdFromTrackId(const art::Event& event, int track_i);

  private:

    art::InputTag fCRTTrueHitLabel;
    art::InputTag fCRTDataLabel;
    art::InputTag fCRTSimHitLabel;
    art::InputTag fCRTTrackLabel;

    bool fRollupUnsavedIds;
    std::map<int, std::map<int, double>> fTrueHitTrueIds;
    std::map<int, std::map<int, double>> fDataTrueIds;
    std::map<int, std::map<int, double>> fSimHitTrueIds;
    std::map<int, std::map<int, double>> fTrackTrueIds;

};

#endif
