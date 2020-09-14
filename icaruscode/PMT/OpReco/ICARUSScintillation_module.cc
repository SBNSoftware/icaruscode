// -*- mode: c++; c-basic-offset: 2; -*-
// Gleb Sinev, Duke, 2016
//
// This module finds periods of time-localized activity
// on each channel of the optical system and creates OpHits.
//
// Split from OpFlashFinder_module.cc
// by Ben Jones, MIT, 2013
//

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
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

#include "art_root_io/TFileService.h"

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/Exception.h"

// ROOT includes
#include <TTree.h>
#include <TFile.h>
// C++ Includes
#include <map>
#include <string>
#include <memory>

namespace opdet {

  class ICARUSScintillation : public art::EDProducer{
  public:

    // Standard constructor and destructor for an ART module.
    explicit ICARUSScintillation(const fhicl::ParameterSet&);
    virtual ~ICARUSScintillation();

    // The producer routine, called once per event.
    void produce(art::Event&);

    void beginJob();
    void endJob();

  private:
    std::map< int, int >  GetChannelMap();
    std::vector< double > GetSPEScales();
    std::vector< double > GetSPEShifts();

    // The parameters we'll read from the .fcl file.
    std::string fInputModule; // Input tag for OpDetWaveform collection
    std::vector< std::string > fInputLabels;
    std::set< unsigned int > fChannelMasks;

    pmtana::PulseRecoManager  fPulseRecoMgr;
    pmtana::PMTPulseRecoBase* fThreshAlg;
    pmtana::PMTPedestalBase*  fPedAlg;


    std::vector<double> _tstart_v,_tmax_v,_tend_v,_tcross_v;
    std::vector<double> _amp_v,_area_v;
    std::vector<double> _ped_mean_v,_ped_sigma_v;
    std::vector<double> _wf;
    double _tstart;
    int _run,_event,_ch;

    TTree* _wftree;
    TTree* _hittree;
    TTree* _geotree;

  };

}

namespace opdet {
  DEFINE_ART_MODULE(ICARUSScintillation)
}

namespace opdet {

  //----------------------------------------------------------------------------
  // Constructor
  ICARUSScintillation::ICARUSScintillation(const fhicl::ParameterSet & pset):
    EDProducer{pset},
    fPulseRecoMgr()
  {
    // Indicate that the Input Module comes from .fcl
    fInputModule   = pset.get< std::string >("InputModule");
    fInputLabels   = pset.get< std::vector< std::string > >("InputLabels");

    for (auto const& ch : pset.get< std::vector< unsigned int > >
                            ("ChannelMasks", std::vector< unsigned int >()))
      fChannelMasks.insert(ch);

    // Initialize the hit finder algorithm
    auto const hit_alg_pset = pset.get< fhicl::ParameterSet >("HitAlgoPset");
    std::string threshAlgName = hit_alg_pset.get< std::string >("Name");
    if      (threshAlgName == "Threshold")
      fThreshAlg = new pmtana::AlgoThreshold(hit_alg_pset);
    else if (threshAlgName == "SiPM")
      fThreshAlg = new pmtana::AlgoSiPM(hit_alg_pset);
    else if (threshAlgName == "SlidingWindow")
      fThreshAlg = new pmtana::AlgoSlidingWindow(hit_alg_pset);
    else if (threshAlgName == "FixedWindow")
      fThreshAlg = new pmtana::AlgoFixedWindow(hit_alg_pset);
    else if (threshAlgName == "CFD" )
      fThreshAlg = new pmtana::AlgoCFD(hit_alg_pset);
    else throw art::Exception(art::errors::UnimplementedFeature)
                    << "Cannot find implementation for "
                    << threshAlgName << " algorithm.\n";

    auto const ped_alg_pset = pset.get< fhicl::ParameterSet >("PedAlgoPset");
    std::string pedAlgName = ped_alg_pset.get< std::string >("Name");
    if      (pedAlgName == "Edges")
      fPedAlg = new pmtana::PedAlgoEdges(ped_alg_pset);
    else if (pedAlgName == "RollingMean")
      fPedAlg = new pmtana::PedAlgoRollingMean(ped_alg_pset);
    else if (pedAlgName == "UB"   )
      fPedAlg = new pmtana::PedAlgoUB(ped_alg_pset);
    else throw art::Exception(art::errors::UnimplementedFeature)
                    << "Cannot find implementation for "
                    << pedAlgName << " algorithm.\n";

    fPulseRecoMgr.AddRecoAlgo(fThreshAlg);
    fPulseRecoMgr.SetDefaultPedAlgo(fPedAlg);

  }

  void ICARUSScintillation::beginJob()
  {

    art::ServiceHandle<art::TFileService> tfs;
    // analyzie tuple prep


    _wftree = tfs->make<TTree>("wftree","wftree");
    _wftree->Branch("run",&_run,"run/I");
    _wftree->Branch("event",&_event,"event/I");
    _wftree->Branch("ch",&_ch,"ch/I");
    _wftree->Branch("wf",&_wf);
    _wftree->Branch("ped_mean",&_ped_mean_v);
    _wftree->Branch("ped_sigma",&_ped_sigma_v);
    _wftree->Branch("tstart",&_tstart,"tstart/D");

    _hittree = tfs->make<TTree>("hittree","hittree");
    _hittree->Branch("run",&_run,"run/I");
    _hittree->Branch("event",&_event,"event/I");
    _hittree->Branch("ch",&_ch,"ch/I");
    _hittree->Branch("tstart",&_tstart_v);
    _hittree->Branch("tmax",&_tmax_v);
    _hittree->Branch("tend",&_tend_v);
    _hittree->Branch("tcross",&_tcross_v);
    _hittree->Branch("amp",&_amp_v);
    _hittree->Branch("area",&_area_v);
    _hittree->Branch("ped_mean",&_ped_mean_v);
    _hittree->Branch("ped_sigma",&_ped_sigma_v);

    _geotree = tfs->make<TTree>("geotree","geotree");
    std::vector<double> pmtX, pmtY, pmtZ;
    std::vector<double> minX, minY, minZ;
    std::vector<double> maxX, maxY, maxZ;
    auto const geop = lar::providerFrom<geo::Geometry>();
    double PMTxyz[3];
    for(size_t opch=0; opch<geop->NOpChannels(); ++opch) {
      geop->OpDetGeoFromOpChannel(opch).GetCenter(PMTxyz);
      pmtX.push_back(PMTxyz[0]);
      pmtY.push_back(PMTxyz[1]);
      pmtZ.push_back(PMTxyz[2]);
    }
    for(auto iter=geop->begin_TPC(); iter!=geop->end_TPC(); ++iter) {
      auto const& tpc = (*iter);
      minX.push_back(tpc.BoundingBox().MinX());
      minY.push_back(tpc.BoundingBox().MinY());
      minZ.push_back(tpc.BoundingBox().MinZ());
      maxX.push_back(tpc.BoundingBox().MaxX());
      maxY.push_back(tpc.BoundingBox().MaxY());
      maxZ.push_back(tpc.BoundingBox().MaxZ());
    }
    _geotree->Branch("pmtX",&pmtX);
    _geotree->Branch("pmtY",&pmtY);
    _geotree->Branch("pmtZ",&pmtZ);
    _geotree->Branch("minX",&minX);
    _geotree->Branch("minY",&minY);
    _geotree->Branch("minZ",&minZ);
    _geotree->Branch("maxX",&maxX);
    _geotree->Branch("maxY",&maxY);
    _geotree->Branch("maxZ",&maxZ);
    _geotree->Fill();

  }

  void ICARUSScintillation::endJob()
  {

  }

  //----------------------------------------------------------------------------
  // Destructor
  ICARUSScintillation::~ICARUSScintillation()
  {

    delete fThreshAlg;
    delete fPedAlg;

  }

  //----------------------------------------------------------------------------
  void ICARUSScintillation::produce(art::Event& evt)
  {

    _event = evt.id().event();
    _run   = evt.id().run();

    // These is the storage pointer we will put in the event
    std::unique_ptr< std::vector< recob::OpHit > >
                      HitPtr(new std::vector< recob::OpHit >);

    //
    // Get the pulses from the event
    //

    if(fInputLabels.empty()) fInputLabels.push_back("");
    for (auto label : fInputLabels) {

      art::Handle< std::vector< raw::OpDetWaveform > > wfHandle;
      evt.getByLabel(fInputModule, label, wfHandle);
      if (!wfHandle.isValid()) continue; // Skip non-existent collections

      for(auto const& waveform : (*wfHandle)) {

	       _ch = static_cast< int >(waveform.ChannelNumber());
	       _tstart = waveform.TimeStamp();

	       _wf.clear();
	       _wf.resize(waveform.size());
	       for(size_t idx=0; idx<_wf.size(); ++idx)
	           _wf[idx] = waveform[idx];

	       fPulseRecoMgr.Reconstruct(waveform);

	       // Record waveforms
	       _ped_mean_v = fPedAlg->Mean();
	       _ped_sigma_v = fPedAlg->Sigma();

         for(size_t idx=0; idx<_wf.size(); ++idx)
             _wf[idx] = waveform[idx] - _ped_mean_v[idx];

	       _wftree->Fill();

         

	       // Record pulses
         /*
         auto const& pulses = fThreshAlg->GetPulses();
	       size_t npulse = pulses.size();


	       _tstart_v.resize(npulse); _tmax_v.resize(npulse); _tend_v.resize(npulse); _tcross_v.resize(npulse);
	       _amp_v.resize(npulse); _area_v.resize(npulse);
	       _ped_mean_v.resize(npulse); _ped_sigma_v.resize(npulse);

	       for(size_t idx=0; idx<npulse; ++idx) {
	          auto const& pulse = pulses[idx];
	          _tstart_v[idx] = pulse.t_start;
	          _tmax_v[idx]   = pulse.t_max;
	          _tend_v[idx]   = pulse.t_end;
	          _tcross_v[idx] = pulse.t_cfdcross;
	          _amp_v[idx]    = pulse.peak;
	          _area_v[idx]   = pulse.area;
	          _ped_mean_v[idx]  = pulse.ped_mean;
	          _ped_sigma_v[idx] = pulse.ped_sigma;
	        }
	       _hittree->Fill();
         */

      } // waveform
    }//input labels



  } // produce
}
