////////////////////////////////////////////////////////////////////////
// Class:       OpHitFinderICARUS
// Plugin Type: producer (art v2_09_06)
// File:        OpHitFinderICARUS_module.cc
//
// Generated at Wed Feb 14 15:51:50 2018 by Andrea Falcone using cetskelgen
// from cetlib version v3_01_03.
////////////////////////////////////////////////////////////////////////

//extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
//}

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

class OpHitFinderICARUS : public art::EDProducer {
public:
  explicit OpHitFinderICARUS(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  OpHitFinderICARUS(OpHitFinderICARUS const &) = delete;
  OpHitFinderICARUS(OpHitFinderICARUS &&) = delete;
  OpHitFinderICARUS & operator = (OpHitFinderICARUS const &) = delete;
  OpHitFinderICARUS & operator = (OpHitFinderICARUS &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

private:

    size_t fEvNumber;
    size_t fChNumber;

    std::string fInputModuleName;
    double fSPEArea; //conversion between phe and Adc*ns 
    double fHitThreshold;
    int fBaselineSample;

    double TimeVector[10000];
    double ADCVector[10000];


};


OpHitFinderICARUS::OpHitFinderICARUS(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{

    produces<std::vector<recob::OpHit>>();

    fInputModuleName = p.get< std::string >("InputModule" );
    fHitThreshold    = p.get< double >("HitThreshold");
    fSPEArea         = p.get< double >("SPEArea");
    fBaselineSample  = p.get< int >("BaselineSample");

}

void OpHitFinderICARUS::produce(art::Event & e)
{
    std::cout << "My module on event #" << e.id().event() << std::endl;

    //art::ServiceHandle<art::TFileService> tfs;
    fEvNumber = e.id().event();

    std::unique_ptr< std::vector< recob::OpHit > > pulseVecPtr(std::make_unique< std::vector< recob::OpHit > > ());  

    art::Handle< std::vector< raw::OpDetWaveform > > wfHandle;
    e.getByLabel(fInputModuleName, wfHandle);

    if(!wfHandle.isValid()){
      std::cout <<Form("Did not find any G4 photons from a producer: %s", "largeant") << std::endl;
    }

    std::cout << "Dimensione primo " << wfHandle->size() << std::endl; 

//    for(size_t wftime; wftime< wfHandle.size(); wftime++)
      for(auto const& wvf : (*wfHandle)){

        fChNumber = wvf.ChannelNumber();
  	std::cout << "Photon channel: " << fChNumber << std::endl;
    	// Load pulses into WaveformVector
    	//std::vector < short_t >  WaveformVector = wvf.Waveform();

	int grsize = wvf.size();

    	//std::vector<Double_t> TimeVector;
	//std::vector<Double_t> ADCVector;

	//TimeVector.size()= grsize;
	//ADCVector.size()= grsize;

    	for (int wtime=0; wtime< grsize; wtime++)
      	{
		TimeVector[wtime] = wtime;
		ADCVector[wtime]  = wvf[wtime];
      	}     

    	float baseline=0;

    	for (int btime =0; btime< fBaselineSample; btime++)
      	{
		baseline = baseline+ ADCVector[btime];
      	}     

    	baseline = baseline/fBaselineSample;

    	std::cout << "Baseline " << baseline << std::endl; 


    	TGraph *gr = new TGraph(grsize,TimeVector,ADCVector);
    
//    	double min = std::min_element(WaveformVector.begin(), WaveformVector.end());
//    	double min_time = distance(WaveformVector.begin(), find(WaveformVector.begin(), WaveformVector.end(), min));    

    	int n_graph = gr->GetN(); // e` 600 ma lo faccio trovare lo stesso
    	double *y_graph = gr->GetY();

    	int min_time = TMath::LocMax(n_graph,y_graph); 
        double min_time_to_put = TMath::LocMax(n_graph,y_graph);    
    	double min = y_graph[min_time];

    	std::cout << "Min " << min << std::endl; 
  
	if (min>baseline+fHitThreshold)
	{
    	TF1 *funz= new TF1("funz", "pol1(0)", min_time-4.0, min_time);

    	gr->Fit("funz","R");

    	double start_moment = funz->GetX(baseline, 0, min_time);

    	std::cout << "Start " << start_moment << std::endl; 

    	TF1 *gauss_start= new TF1("gauss_start", "gaus", 0, min_time);
    	TF1 *gauss_end  = new TF1("gauss_end", "gaus", min_time, n_graph);

	gr->Fit("gauss_start","R");

	double Constant1 = gauss_start->GetParameter(0);
	//double Mean1 = gauss_start->GetParameter(1);
	double Sigma1 = gauss_start->GetParameter(2);

	gr->Fit("gauss_end","R");
   
	double Constant2 = gauss_end->GetParameter(0);
	//double Mean2 = gauss_end->GetParameter(1);
	double Sigma2 = gauss_end->GetParameter(2);

	double Area =  ((Constant1*Sigma1)/2 + (Constant2*Sigma2)/2)*sqrt(2*3.14159);

	unsigned short frame = 1;

	double fasttotal = 3/4;

	double time_abs = sqrt(min_time_to_put);

	double FWHM = 2.35*((Sigma1+Sigma2)/2);

	double phelec= Area/fSPEArea;

	recob::OpHit adcVec(fChNumber, min_time_to_put, time_abs, frame, FWHM, Area, min, phelec, fasttotal);//including hit info
	pulseVecPtr->emplace_back(std::move(adcVec));

	funz->~TF1();
	gauss_start->~TF1();
	gauss_end->~TF1();
	gr->~TGraph();

	}
	else {std::cout << "No OpHit in channel " << fChNumber << std::endl;}
    } 
// Store results into the event
e.put(std::move(pulseVecPtr));
}

DEFINE_ART_MODULE(OpHitFinderICARUS)

}
