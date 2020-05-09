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
#include "icaruscode/CRT/CRTProducts/CRTChannelData.h"
#include "icaruscode/CRT/CRTProducts/CRTData.hh"
//#include "icaruscode/CRT/CRTProducts/CRTHit.h"
#include "icaruscode/CRT/CRTUtils/CRTBackTracker.h"

#include <vector>

namespace icarus{ 
 namespace crt {
    class CRTTruthMatchAnalysis;
 }
}

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

  art::InputTag fCRTDataLabel; 
  CRTBackTracker bt;
};


icarus::crt::CRTTruthMatchAnalysis::CRTTruthMatchAnalysis(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    fCRTDataLabel(p.get<art::InputTag>("CRTDataLabel","crtdaq")),
    bt(p.get<fhicl::ParameterSet>("CRTBackTrack"))
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void icarus::crt::CRTTruthMatchAnalysis::analyze(art::Event const& e)
{

  std::cout << "initializing CRTBackTracker..." << std::endl;
  bt.Initialize(e);
  std::cout << "done." << std::endl;

  art::Handle< std::vector<CRTData> > dataHandle;
  //std::vector< art::Ptr<CRTData> > dataList;
  e.getByLabel(fCRTDataLabel,dataHandle);
      //art::fill_ptr_vector(dataList, dataHandle);

  std::cout << "found CRTData vector with " << (*dataHandle).size() << " entries" << std::endl;

  int nsizematch=0, ntracksmatch=0, n0bt=0;
  for(auto const& data : *dataHandle) {

      std::vector<int> btIds = bt.AllTrueIds(e,data);
      if(btIds.empty()) n0bt++;

      std::vector<int> dataIds;
      for(auto const& chdata : data.ChanData()) {
          for(auto const& id : chdata.TrackID())
              dataIds.push_back(id);
      }

      std::sort(dataIds.begin(),dataIds.end());
      auto last = std::unique(dataIds.begin(),dataIds.end());
      dataIds.erase(last,dataIds.end());

      if(dataIds.size()!=btIds.size()) {
          std::cout << "ids size mismatch data - bt " << dataIds.size()-btIds.size() << std::endl;
          std::cout << "   dataIds size = " << dataIds.size() << ", btIds size = " << btIds.size() << std::endl;
          for(auto const& id : dataIds) std::cout << "    dataId: " << id << std::endl;
          for(auto const& id : btIds) std::cout << "    btId:   " << id << std::endl;
      }//if size mismatch 
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
  }//loop over CRTData

  std::cout << "found " << nsizematch << " size matches" << std::endl;
  std::cout << "found " << ntracksmatch << " track ID set matches" << std::endl;
  std::cout << "found " << n0bt << " empty back tracker results" << std::endl;
}

DEFINE_ART_MODULE(icarus::crt::CRTTruthMatchAnalysis)
