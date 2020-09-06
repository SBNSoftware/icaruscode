#ifndef CRTTRACKRECOALG_H_SEEN
#define CRTTRACKRECOALG_H_SEEN

//////////////////////////////////////////////////////////////////////////////////
// CRTTrackRecoAlg.h
//
// Functions for CRT track reconstruction
// written by T Brooks (tbrooks@fnal.gov), November 2018
// ported to and modified for use with icaruscode by Chris.Hilgenberg@colostate.edu
///////////////////////////////////////////////////////////////////////////////////

// framework
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "canvas/Persistency/Common/FindManyP.h"

// LArSoft
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

// Utility libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"
#include "cetlib/pow.h" // cet::sum_of_squares()
#include "icaruscode/CRT/CRTProducts/CRTHit.hh"
#include "icaruscode/CRT/CRTProducts/CRTTrack.hh"
#include "icaruscode/CRT/CRTUtils/CRTHitRecoAlg.h"

// c++
#include <iostream>
#include <stdio.h>
#include <sstream>
#include <vector>
#include <map>
#include <utility>
#include <cmath> 
#include <memory>

// ROOT
#include "TVector3.h"
#include "TGeoManager.h"

using std::vector;
using std::pair;
using std::map;

namespace icarus{
namespace crt{

  class CRTTrackRecoAlg {
  public:

    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<double> TimeLimit {
        Name("TimeLimit"),
        Comment("")
      };

      fhicl::Atom<double> AverageHitDistance {
        Name("AverageHitDistance"),
        Comment("Distance to average hits over on same plane")
      };

      fhicl::Atom<double> DistanceLimit {
        Name("DistanceLimit"),
        Comment("Distance to combine CRT hits into track")
      };

    };

    CRTTrackRecoAlg(const Config& config);

    CRTTrackRecoAlg(const fhicl::ParameterSet& pset) :
      CRTTrackRecoAlg(fhicl::Table<Config>(pset, {})()) {}

    CRTTrackRecoAlg(double aveHitDist, double distLim);

    ~CRTTrackRecoAlg();

    void reconfigure(const Config& config);

    vector<vector<art::Ptr<CRTHit>>> CreateCRTTzeros(vector<art::Ptr<CRTHit>>);

    // Function to make creating CRTTracks easier
    CRTTrack FillCrtTrack(CRTHit hit1, CRTHit hit2, bool complete);

    // Function to average hits within a certain distance of each other
    vector<pair<CRTHit, vector<int>>> AverageHits(vector<art::Ptr<CRTHit>> hits, map<art::Ptr<CRTHit>, int> hitIds);
    vector<CRTHit> AverageHits(vector<art::Ptr<CRTHit>> hits);

    // Take a list of hits and find average parameters
    CRTHit DoAverage(vector<art::Ptr<CRTHit>> hits);

    // Create CRTTracks from list of hits
    vector<pair<CRTTrack, vector<int>>> CreateTracks(vector<pair<CRTHit, vector<int>>> hits);
    vector<CRTTrack> CreateTracks(vector<CRTHit> hits);

    // Calculate the tagger crossing point of CRTTrack candidate
    TVector3 CrossPoint(CRTHit hit, TVector3 start, TVector3 diff);

  private:

    geo::GeometryCore const* fGeometryService;

    double fTimeLimit;
    double fAverageHitDistance;
    double fDistanceLimit;

    CRTHitRecoAlg hitAlg;

  };

}//namespace crt
}//namespace icarus

#endif
