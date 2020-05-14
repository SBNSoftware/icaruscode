////////////////////////////////////////////////////////////////////////
// Class:       CRTTruthMatchAnalysis
// Plugin Type: analyzer (art v3_05_00)
// File:        CRTTruthMatchAnalysis_module.cc
//
// Generated at Fri Apr 24 14:45:31 2020 by Christopher Hilgenberg using cetskelgen
// from cetlib version v3_10_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "icaruscode/CRT/CRTProducts/CRTData.hh"
#include "icaruscode/CRT/CRTProducts/CRTHit.hh"
#include "icaruscode/CRT/CRTUtils/CRTBackTracker.h"

#include <vector>
#include <map>

namespace icarus{ 
 namespace crt {
    class CRTTruthMatchAnalysis;
 }
}

using std::vector;
using std::map;

class icarus::crt::CRTTruthMatchAnalysis : public art::EDAnalyzer {
public:
  explicit CRTTruthMatchAnalysis(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CRTTruthMatchAnalysis(CRTTruthMatchAnalysis const&) = delete;
  CRTTruthMatchAnalysis(CRTTruthMatchAnalysis&&) = delete;
  CRTTruthMatchAnalysis& operator=(CRTTruthMatchAnalysis const&) = delete;
  CRTTruthMatchAnalysis& operator=(CRTTruthMatchAnalysis&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:

  art::InputTag fSimulationLabel;
  art::InputTag fCRTTrueHitLabel;
  art::InputTag fCRTDataLabel;
  art::InputTag fCRTSimHitLabel;
  CRTBackTracker bt;
};


icarus::crt::CRTTruthMatchAnalysis::CRTTruthMatchAnalysis(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    fSimulationLabel(p.get<art::InputTag>("SimulationLabel","largeant")),
    fCRTTrueHitLabel(p.get<art::InputTag>("CRTTrueHitLabel","crttruehit")),
    fCRTDataLabel(p.get<art::InputTag>("CRTDataLabel","crtdaq")),
    fCRTSimHitLabel(p.get<art::InputTag>("CRTSimHitLabel","crtsimhit")),
    bt(p.get<fhicl::ParameterSet>("CRTBackTrack"))
{
}

void icarus::crt::CRTTruthMatchAnalysis::analyze(art::Event const& e)
{

  //std::cout << "initializing BackTracker..." << std::endl;
  bt.Initialize(e);
  //std::cout << "done." << std::endl;

  //std::cout << "loop over AuxDetChannels" << std::endl;
  art::Handle< vector<sim::AuxDetSimChannel> > adscHandle;
  vector< art::Ptr<sim::AuxDetSimChannel> > adscList;
  if( e.getByLabel(fSimulationLabel,adscHandle) )
      art::fill_ptr_vector(adscList,adscHandle);

  map<int,bool> trueRecoedTracks;
  for(auto const& adsc : adscList) {
      for(auto const& ide : adsc->AuxDetIDEs()) {
          trueRecoedTracks[ide.trackID] = false;
      }
  }

  //std::cout << "loop over CRTTrueHits" << std::endl;
  art::Handle< vector<CRTHit> > trueHitHandle;
  vector< art::Ptr<CRTHit> > trueHitList;
  if( e.getByLabel(fCRTTrueHitLabel,trueHitHandle) )
      art::fill_ptr_vector(trueHitList,trueHitHandle);

  for(auto const& hit : trueHitList) {
      //std::cout << "fetching trackIDs from BackTracker..." << std::endl;
      vector<int> btIds = bt.AllTrueIds(e,*hit);
      //std::cout << "  found " << btIds.size() << std::endl;
      for(const int id: btIds) {
          if(trueRecoedTracks.find(id)!=trueRecoedTracks.end()) 
              trueRecoedTracks[id] = true;
          else
              std::cout << "trackID from BackTracker not found in AuxDetSimChannels!" << std::endl;
      }
  }

  //std::cout << "matching trackIDs" << std::endl;
  size_t nmiss=0, ntot=trueRecoedTracks.size();
  for(auto const& trk : trueRecoedTracks) {
      if(trk.second==false) nmiss++;
  }

  std::cout << ntot-nmiss << " (" << 100.0*(ntot-nmiss)/ntot << " %) tracks matched to true hits" << std::endl;

  /*art::Handle< vector<CRTData> > dataHandle;
  //std::vector< art::Ptr<CRTData> > dataList;
  e.getByLabel(fCRTDataLabel,dataHandle);
      //art::fill_ptr_vector(dataList, dataHandle);

  std::cout << "found CRTData vector with " << (*dataHandle).size() << " entries" << std::endl;

  int nsizematch=0, ntracksmatch=0, n0bt=0;
  for(auto const& data : *dataHandle) {

      std::vector<int> btIds = bt.AllTrueIds(e,data);

      std::vector<int> dataIds;
      for(auto const& chdata : data.ChanData()) {
          for(auto const& id : chdata.TrackID())
              dataIds.push_back(id);
      }

      else {
          nsizematch++;
          bool pass=true;
	  for(size_t i=0; i<dataIds.size(); i++) {
                bool found = false;
		for(size_t j=0; j<btIds.size(); j++) {
                    if (dataIds[i] == btIds[j]){
                        found=true;
                        break;
                    }
                }
                if(!found)
                    pass = false;
	  }
          if(pass)
              ntracksmatch++;
      }//else size match
  }//loop over CRTData*/

}

DEFINE_ART_MODULE(icarus::crt::CRTTruthMatchAnalysis)
