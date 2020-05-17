#ifndef ICARUS_CRTTRUEHITRECOALG_H
#define ICARUS_CRTTRUEHITRECOALG_H

// LArSoft includes
#include "lardataobj/Simulation/AuxDetSimChannel.h"

// Utility libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"

// icaruscode includes
#include "icaruscode/CRT/CRTProducts/CRTHit.hh"
#include "icaruscode/CRT/CRTUtils/CRTCommonUtils.h"

// C++ includes
#include <vector>
#include <utility>
#include <map>
#include <set>

// ROOT includes
//#include "Math/GenVector/XYZTVector.h"
//#include "Math/GenVector/LorentzVector.h" 
#include "TLorentzVector.h"

using std::vector;
using std::pair;
using std::map;
using std::set;
//using ROOT::Math::XYZTVector;

namespace icarus {
 namespace crt {
    class CRTTrueHitRecoAlg;
 }
}

struct tagger {
    char type;
    int region;
    std::set<int> layerID;
    map<int,int> stripLayer;
    map<int,sim::AuxDetIDE> stripIDE;
};

class icarus::crt::CRTTrueHitRecoAlg {

 public:

    struct Config {
        using Name = fhicl::Name;
        using Comment = fhicl::Comment;

        fhicl::Atom<bool> UseReadoutWindow {
          Name("UseReadoutWindow"),
          Comment("Only reconstruct hits within readout window")
        };
        fhicl::Atom<double> EDepMin {
          Name("EDepMin"),
          Comment("Lowest consierd deposited energy in a scintillator strip used in hit [MeV]")
        };
        fhicl::Atom<bool> RollupUnusedIds {
          Name("RollupUnusedIds"),
          Comment("merge G4-untracked trackIDs into partent track")
        };
    };

    CRTTrueHitRecoAlg(const Config& config);

    CRTTrueHitRecoAlg(const fhicl::ParameterSet& pset) :
      CRTTrueHitRecoAlg(fhicl::Table<Config>(pset, {})()) {}

    CRTTrueHitRecoAlg();
    ~CRTTrueHitRecoAlg();

    void reconfigure(const Config& config);

    vector<pair<CRTHit,vector<sim::AuxDetIDE>>> CreateCRTHits(vector<art::Ptr<sim::AuxDetSimChannel>> adscList);

    // Function to make filling a CRTHit a bit faster
    CRTHit FillCrtHit(vector<uint8_t> tfeb_id, map<uint8_t,vector<pair<int,float>>> tpesmap, 
                   float peshit, double time0, double time1, int plane,
                   double x, double ex, double y, double ey, double z, double ez, std::string tagger);

 private:

    geo::GeometryCore const* fGeometryService;
    map<int,vector<pair<uint8_t,int>>> fFebMap;

    //config params
    bool   fUseReadoutWindow;
    double fEDepMin;
    bool   fRollupUnusedIds;

};

#endif
