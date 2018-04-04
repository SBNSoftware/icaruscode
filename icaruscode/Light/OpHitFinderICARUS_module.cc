////////////////////////////////////////////////////////////////////////
// Class:       OpHitFinderICARUS
// Plugin Type: producer (art v2_09_06)
// File:        OpHitFinderICARUS_module.cc
//
// Generated at Wed Feb 14 15:51:50 2018 by Andrea Falcone using cetskelgen
// from cetlib version v3_01_03.
////////////////////////////////////////////////////////////////////////

#include <sys/types.h>
#include <sys/stat.h>

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "larana/OpticalDetector/OpHitFinder/PMTPulseRecoBase.h"
#include "larana/OpticalDetector/OpHitFinder/AlgoThreshold.h"
#include "larana/OpticalDetector/OpHitFinder/AlgoSiPM.h"
#include "larana/OpticalDetector/OpHitFinder/AlgoSlidingWindow.h"
#include "larana/OpticalDetector/OpHitFinder/AlgoFixedWindow.h"
#include "larana/OpticalDetector/OpHitFinder/AlgoCFD.h"
#include "larana/OpticalDetector/OpHitFinder/PedAlgoEdges.h"
#include "larana/OpticalDetector/OpHitFinder/PedAlgoRollingMean.h"
#include "larana/OpticalDetector/OpHitFinder/PedAlgoUB.h"
#include "larana/OpticalDetector/OpHitFinder/PulseRecoManager.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larana/OpticalDetector/OpHitFinder/OpHitAlg.h"
#include "lardataobj/Simulation/BeamGateInfo.h"

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Utilities/make_tool.h"
#include "art/Framework/Principal/DataViewImpl.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/Exception.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "icaruscode/Light/OpticalTools/IOpHitFinder.h"

//ROOT includes
#include "TF1.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"

#include <memory>
#include <vector>
#include <cmath>
#include <string>
#include <map>
#include <sstream>
#include <fstream>

namespace ophit{

//class OpHitFinderICARUS;

class OpHitFinderICARUS : public art::EDProducer
{
public:
    explicit OpHitFinderICARUS(fhicl::ParameterSet const & p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    OpHitFinderICARUS(OpHitFinderICARUS const &)               = delete;
    OpHitFinderICARUS(OpHitFinderICARUS &&)                    = delete;
    OpHitFinderICARUS& operator = (OpHitFinderICARUS const &)  = delete;
    OpHitFinderICARUS& operator = (OpHitFinderICARUS &&)       = delete;

    // Required functions.
    void produce(art::Event & e) override;

private:

    size_t fEvNumber;

    std::string fInputModuleName;
    
    std::unique_ptr<light::IOpHitFinder> fOpHitFinder;
};

OpHitFinderICARUS::OpHitFinderICARUS(fhicl::ParameterSet const & p)
{
    produces<std::vector<recob::OpHit>>();

    fInputModuleName = p.get< std::string >("InputModule" );
    
    fOpHitFinder = art::make_tool<light::IOpHitFinder> (p.get<fhicl::ParameterSet>("OpHitFinder"));
}

void OpHitFinderICARUS::produce(art::Event & e)
{
    std::cout << "My module on event #" << e.id().event() << std::endl;

    //art::ServiceHandle<art::TFileService> tfs;
    fEvNumber = e.id().event();

    std::unique_ptr<std::vector<recob::OpHit>> pulseVecPtr(std::make_unique<std::vector<recob::OpHit>>());  

    art::Handle< std::vector< raw::OpDetWaveform > > wfHandle;
    e.getByLabel(fInputModuleName, wfHandle);

    if(!wfHandle.isValid())
    {
      std::cout <<Form("Did not find any G4 photons from a producer: %s", "largeant") << std::endl;
    }

    std::cout << "Dimensione primo " << wfHandle->size() << std::endl; 

//    for(size_t wftime; wftime< wfHandle.size(); wftime++)
    for(auto const& wvf : (*wfHandle))
    {
        light::OpHitVec opHitVec;
        
        fOpHitFinder->FindOpHits(wvf, opHitVec);
        
        for(auto& opHit : opHitVec)
            pulseVecPtr->emplace_back(opHit);
    }
    // Store results into the event
    e.put(std::move(pulseVecPtr));
}

DEFINE_ART_MODULE(OpHitFinderICARUS)

}
