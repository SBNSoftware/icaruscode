////////////////////////////////////////////////////////////////////////
// \file    CVNZlibMakerICARUS_module.cc
// \brief   Analyzer module for creating CVN gzip file objects
// \author  V Hewes - vhewes@fnal.gov
//          Saul Alonso-Monsalve - saul.alonso.monsalve@cern.ch
////////////////////////////////////////////////////////////////////////

// C/C++ includes
#include <iostream>
#include <cstdlib>

#include "boost/filesystem.hpp"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/Exception.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

// Data products
#include "nusimdata/SimulationBase/GTruth.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataalg/DetectorInfo/DetectorClocksData.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "icaruscode/RecoUtils/RecoUtils.h"

// CVN includes
#include "larrecodnn/CVN/func/AssignLabels.h"
#include "larrecodnn/CVN/func/LArTrainingData.h"
#include "larrecodnn/CVN/func/InteractionType.h"
#include "larrecodnn/CVN/func/PixelMap.h"
#include "larrecodnn/CVN/func/CVNImageUtils.h"
#include "icaruscode/TPC/CVN/interfaces/ICVNZlibMakerICARUS.h"

// Compression
#include "zlib.h"
#include "math.h"

#include "TH1.h"
#include "TTree.h"
#include "larcoreobj/SummaryData/POTSummary.h"


namespace lcvn {

  class CVNZlibMakerICARUS : public ICVNZlibMakerICARUS {
  public:

    explicit CVNZlibMakerICARUS(fhicl::ParameterSet const& pset);
    ~CVNZlibMakerICARUS();

    void beginJob();
    void endSubRun(art::SubRun const &sr);
    void analyze(const art::Event& evt);
    void reconfigure(const fhicl::ParameterSet& pset);

  private:

    unsigned int fTopologyHitsCut;
    bool fVerbose;
    bool fUseBackTrackInfo;
    bool fUseNuContainment;
    double fCosEfrac;
    double fNuEfrac;
    double fVolCut;
    bool fSaveTree;
    std::string fGenieGenModuleLabel;
    std::vector<std::string> fPandoraTagSuffixes;
    std::string fPFParticleModuleLabel;
    std::string fHitModuleLabel;
    std::string fT0Label;
    
    // ParticleInventoryService
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;

    // BackTrackerService
    art::ServiceHandle<cheat::BackTrackerService> bt_serv;

    void write_files(LArTrainingNuData td, std::string evtid);
    //bool Is_in_TPC(double stX, double stY, double stZ);
    bool Is_in_Fiducial_Vol(double stX, double stY, double stZ);
    void HitTruth(detinfo::DetectorClocksData const& clockData, art::Ptr<recob::Hit> const& hit, Int_t& truthid);
    bool HitTruthId(detinfo::DetectorClocksData const& clockData, art::Ptr<recob::Hit> const& hit, Int_t& mcid);
    bool TrackIdToMCTruth(Int_t const trkID, art::Ptr<simb::MCTruth>& mctruth);
    double HitEfrmTrkID(detinfo::DetectorClocksData const& clockData, art::Ptr<recob::Hit> const& hit, Int_t truthid);
    double HitTotE(detinfo::DetectorClocksData const& clockData, art::Ptr<recob::Hit> const& hit);
    bool Is_Slice_Nu(art::FindManyP<recob::PFParticle> const& assn, art::Ptr<recob::Slice> const& slce);
    art::Ptr<recob::PFParticle> Get_Nu_like_PFP(art::FindManyP<recob::PFParticle> const& assn, art::Ptr<recob::Slice> const& slce);
    double Get_Slice_Score(art::FindManyP<larpandoraobj::PFParticleMetadata> const& assn, art::Ptr<recob::PFParticle> const& pfp);
    std::vector<int> Get_Best_Slice_ID_PFP_pdg(art::FindManyP<recob::PFParticle> const& pfp_assn, art::FindManyP<larpandoraobj::PFParticleMetadata> const& metadata_assn, std::vector<art::Ptr<recob::Slice>> const& slice_vec);
    int Get_Slice_PFP_ID(art::FindManyP<recob::PFParticle> const& assn, art::Ptr<recob::Slice> const& slce);
    int Get_Hit_Count(unsigned int tpc_no, unsigned int pln_no, std::vector<art::Ptr<recob::Hit>> const& hit_vec);
    double Get_Tot_nu_E(detinfo::DetectorClocksData const& clockData,art::Ptr<simb::MCTruth> const& mcneutrino, std::vector<art::Ptr<recob::Hit>> const& hit_vec);
    
    double HitTotE_frm_given_orgin(detinfo::DetectorClocksData const& clockData, art::Ptr<recob::Hit> const& hit, std::string origin_name);
    double HitTotE_frm_mctruth(detinfo::DetectorClocksData const& clockData, art::Ptr<recob::Hit> const& hit, art::Ptr<simb::MCTruth> const& my_mctruth);
    double TotE_from_mctruth(detinfo::DetectorClocksData const& clockData, std::vector<art::Ptr<recob::Hit>> const& hit_vec, art::Ptr<simb::MCTruth> const& my_mctruth);

    void Clear();

    TH1D* hPOT;
    double fPOT;
    int fRun;
    int fSubRun; 
    
    TTree *fEventTree;
    int frun;
    int fsubrun;
    int fevent;
    int fsliceID;
    double ftotsliceE;
    double ftotsliceNuE;
    double ftotsliceCosE;
    double ftotsliceOthE;
    double fsliceNuEfrac;
    double fsliceCosEfrac;
    double fsliceOthEfrac;	    
    bool fIsSliceNu;
    double fSliceScore;	
    bool fIsbestSlice;
    int fbestpfppdg;
    int fbestsliceID;
    int fpfppdg;
    bool fNuIsContained;
    int fNuID;
    int fNhits_total;
    double fNuSlice_purity;
    double fNuSlice_completeness;
    double fT0;
  };

  //......................................................................
  CVNZlibMakerICARUS::CVNZlibMakerICARUS(fhicl::ParameterSet const& pset):ICVNZlibMakerICARUS(pset)
  {
    reconfigure(pset);
  }

  //......................................................................
  CVNZlibMakerICARUS::~CVNZlibMakerICARUS()
  {

  }

  //......................................................................
  void CVNZlibMakerICARUS::reconfigure(const fhicl::ParameterSet& pset)
  {
    ICVNZlibMakerICARUS::reconfigure(pset);
    fTopologyHitsCut = pset.get<unsigned int>("TopologyHitsCut");
    fVerbose = pset.get<bool>("verbose");
    fUseBackTrackInfo = pset.get<bool>("UseBackTrackInfo");
    fUseNuContainment = pset.get<bool>("UseNuContainment");
    fCosEfrac = pset.get<double>("CosEfrac");
    fNuEfrac = pset.get<double>("NuEfrac");
    fVolCut = pset.get<double>("VolCut");
    fSaveTree = pset.get<bool>("SaveTree");
    fGenieGenModuleLabel = pset.get<std::string>("GenieGenModuleLabel");
    fPandoraTagSuffixes = pset.get<std::vector<std::string>>("PandoraTagSuffixes");
    fPFParticleModuleLabel = pset.get<std::string>("PFParticleModuleLabel");
    fHitModuleLabel = pset.get<std::string>("HitModuleLabel");
    fT0Label = pset.get<std::string>("T0Label");
  }

  //......................................................................
  void CVNZlibMakerICARUS::endSubRun(const art::SubRun & sr)
  {
    std::string fPOTModuleLabel = "generator";
    fRun = sr.run();
    fSubRun = sr.subRun();

    art::Handle< sumdata::POTSummary > potListHandle;
    if(sr.getByLabel(fPOTModuleLabel,potListHandle))
      fPOT = potListHandle->totpot;
    else
      fPOT = 0.;
    if(hPOT) hPOT->Fill(0.5, fPOT);
  }

  //...................................................................... 
  void CVNZlibMakerICARUS::beginJob()
  {
    ICVNZlibMakerICARUS::beginJob();

    art::ServiceHandle<art::TFileService> tfs;
    hPOT = tfs->make<TH1D>("TotalPOT", "Total POT;; POT", 1, 0, 1);
    
  }

  //......................................................................  
  void CVNZlibMakerICARUS::analyze(const art::Event& evt)
  {

    std::cout << "CVNZlibMakerICARUS::analyze" << std::endl;
	  
    std::vector<art::Ptr<lcvn::PixelMap>> pixelmaps;
    art::InputTag itag1(fPixelMapInput, "cvnmap");
    auto pixelmaps_handle = evt.getHandle<std::vector<lcvn::PixelMap>>(itag1);
    if (pixelmaps_handle){
      art::fill_ptr_vector(pixelmaps, pixelmaps_handle);
      std::cout << "if(pixelmaps_handle) fill_ptr_vector" << std::endl;
    }
    if (pixelmaps.size() == 0){
      std::cout << "pixelmaps.size == 0"; return;
    }

    AssignLabels labels;
    TDNuInfo info;
    double event_weight = 1;
    InteractionType interaction = kOther;
    bool Nu_is_contained = false;
       
    std::vector<art::Ptr<simb::MCTruth>> mctruths;
    auto mctruths_handle = evt.getHandle<std::vector<simb::MCTruth>>(fGenieGenModuleLabel);

    if (mctruths_handle){
      art::fill_ptr_vector(mctruths, mctruths_handle);
      std::cout << "got mctruth handle" << std::endl;

      for(auto const& mctruth : mctruths){
        std::cout << "got mctruths" << std::endl;
        
        if(mctruth->Origin() == simb::kBeamNeutrino){
          simb::MCNeutrino true_neutrino = mctruth->GetNeutrino();

          if(fUseNuContainment){
            if(Is_in_Fiducial_Vol(true_neutrino.Nu().Vx(),true_neutrino.Nu().Vy(),true_neutrino.Nu().Vz())){ 
              interaction = labels.GetInteractionType(true_neutrino);
              labels.GetTopology(mctruth, fTopologyHitsCut);
              float nu_energy = true_neutrino.Nu().E();
              float lep_energy = true_neutrino.Lepton().E();
              fNuID = true_neutrino.Nu().TrackId();
              info.SetTruthInfo(nu_energy, lep_energy, 0., event_weight);
              info.SetTopologyInformation(labels.GetPDG(), labels.GetNProtons(), labels.GetNPions(), labels.GetNPizeros(), labels.GetNNeutrons(), labels.GetTopologyType(), labels.GetTopologyTypeAlt());
              Nu_is_contained = true;
              break;
            } // if in fiducial volume
          } // use containment
          else{
            std::cout << "false fUseNuContainment" << std::endl;
            interaction = labels.GetInteractionType(true_neutrino);
            labels.GetTopology(mctruth, fTopologyHitsCut);
            float nu_energy = true_neutrino.Nu().E();
            float lep_energy = true_neutrino.Lepton().E();
            fNuID = true_neutrino.Nu().TrackId();
            info.SetTruthInfo(nu_energy, lep_energy, 0., event_weight);
            info.SetTopologyInformation(labels.GetPDG(), labels.GetNProtons(), labels.GetNPions(), labels.GetNPizeros(), labels.GetNNeutrons(), labels.GetTopologyType(), labels.GetTopologyTypeAlt());
            if(Is_in_Fiducial_Vol(true_neutrino.Nu().Vx(),true_neutrino.Nu().Vy(),true_neutrino.Nu().Vz())) Nu_is_contained = true;
            break;
          } // dont use containment
        } // kBeamNeutrino

        else if(mctruth->Origin() == simb::kCosmicRay){
          interaction = kCosmic;
          break;
        } // kCosmicRay

      } // mctruths
    } // if mctruths_handle

    // collect the input TPC reco tags
    std::vector<std::string> pandora_tag_suffixes = fPandoraTagSuffixes;
    if (pandora_tag_suffixes.size() == 0) pandora_tag_suffixes.push_back("");

    /////////////////////// use slice section ////////////////////////////////////////
    if(fUseSlice){
      std::cout << "***************** Using slice method ****************\n";

      // collect the TPC hits
      std::vector<art::Ptr<recob::Hit>> hits;
      art::Handle<std::vector<recob::Hit>> hits_handle;
      for (unsigned i_tag = 0; i_tag < pandora_tag_suffixes.size(); i_tag++){
        const std::string &pandora_tag_suffix = pandora_tag_suffixes[i_tag];
        if (evt.getByLabel(fHitModuleLabel + pandora_tag_suffix, hits_handle)){
          art::fill_ptr_vector(hits, hits_handle);
        }
      }

      // collect the TPC slices
      std::vector<art::Ptr<recob::Slice>> slices;
      std::vector<std::string> slice_tag_suffixes;
      std::vector<unsigned> slice_tag_indices;
      art::Handle<std::vector<recob::Slice>> slices_handle;
      for (unsigned i_tag = 0; i_tag < pandora_tag_suffixes.size(); i_tag++) {
        const std::string &pandora_tag_suffix = pandora_tag_suffixes[i_tag];
        if (evt.getByLabel(fSliceLabel + pandora_tag_suffix, slices_handle)) {
          art::fill_ptr_vector(slices, slices_handle);
          for (unsigned i = 0; i < slices.size(); i++) {
            slice_tag_suffixes.push_back(pandora_tag_suffix);
            slice_tag_indices.push_back(i_tag);
          }
        }
      }

      // collect the TPC pfps
      std::vector<art::Ptr<recob::PFParticle>> pfps;
      art::Handle<std::vector<recob::PFParticle>> pfps_handle;
      for (unsigned i_tag = 0; i_tag < pandora_tag_suffixes.size(); i_tag++) {
        const std::string &pandora_tag_suffix = pandora_tag_suffixes[i_tag];
        if (evt.getByLabel(fPFParticleModuleLabel + pandora_tag_suffix, pfps_handle)) {
          art::fill_ptr_vector(pfps, pfps_handle);
        }
      }

//      art::FindManyP<recob::Hit> findManyHits(slices_handle, evt, fSliceLabel);
//      art::FindManyP<recob::PFParticle> findManyPFPs(slices_handle, evt, fPFParticleModuleLabel);
//      art::FindManyP<larpandoraobj::PFParticleMetadata> fm_pfpmd(pfps_handle, evt, fPFParticleModuleLabel);
//      art::FindManyP<anab::T0> findManyT0s(pfps_handle, evt, fT0Label);
         
      if(fUseBackTrackInfo){
        std::cout << "***************** Using backtracker information to get neutrino information ****************\n";
            
        auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
           
        for(unsigned int pmID=0; pmID<pixelmaps.size(); pmID++){
          Clear();
          frun = evt.run();
          fsubrun = evt.subRun();
          fevent = evt.id().event();
               
          AssignLabels labels;
          TDNuInfo info;
          double event_weight = 1;
          InteractionType interaction = kOther;	  
               
          for(unsigned int sliceID=0; sliceID<slices.size(); sliceID++){
            art::Ptr<recob::Slice> slice = slices[sliceID];
            const std::string &slice_tag_suff = slice_tag_suffixes[sliceID];
            std::vector<art::Ptr<recob::Slice>> sliceList {slice};
            
            art::FindManyP<recob::PFParticle> findManyPFPs(slices_handle, evt, fPFParticleModuleLabel + slice_tag_suff);
            art::FindManyP<recob::Hit> findManyHits(slices_handle, evt, fSliceLabel + slice_tag_suff);
            art::FindManyP<anab::T0> findManyT0s(pfps_handle, evt, fT0Label + slice_tag_suff);
            art::FindManyP<larpandoraobj::PFParticleMetadata> fm_pfpmd(pfps_handle, evt, fPFParticleModuleLabel + slice_tag_suff);

            if(slice->ID() == pixelmaps[pmID]->fSliceID){

              if(findManyHits.isValid()){
                std::vector<art::Ptr<recob::Hit>> slice_hits = findManyHits.at(slice.key());
                double tot_slice_eng = 0;
                double tot_slice_nu_eng = 0;
                double tot_slice_cos_eng = 0;
                fNhits_total = slice_hits.size();
                std::vector<double> truth_energies(mctruths.size(),0.);

                for(auto const hit : slice_hits){
                //Int_t trkId; // newly deleted
                  tot_slice_eng += HitTotE(clockData,hit);
                  tot_slice_nu_eng += HitTotE_frm_given_orgin(clockData,hit,"nu"); // newly added
                  tot_slice_cos_eng += HitTotE_frm_given_orgin(clockData,hit,"cos"); // newly added
                           
                  // newly added section
                  if(mctruths.size()){
                    for(auto const mctruth : mctruths){
                      if(mctruth->Origin() == simb::kBeamNeutrino){
                        truth_energies[mctruth.key()] += HitTotE_frm_mctruth(clockData,hit,mctruth);
                      }
                    }
                  }

                  // newly deleted section
                  /*if(HitTruthId(clockData,hit,trkId)){
                  art::Ptr<simb::MCTruth> mctruth;
                  if(TrackIdToMCTruth(trkId,mctruth)){
                    if(mctruth->Origin() == simb::kBeamNeutrino){
                      tot_slice_nu_eng += HitEfrmTrkID(clockData,hit,trkId); 
                      truth_energies[mctruth.key()] += HitEfrmTrkID(clockData,hit,trkId);
                    }
                  else if(mctruth->Origin() == simb::kCosmicRay){
                    tot_slice_cos_eng += HitEfrmTrkID(clockData,hit,trkId);
                    }
                  } // found a valid MCTruth object
                  }*/ // hit is matched to a trkID
                
                } // loop over hits in the selected slice
                     
                ftotsliceE = tot_slice_eng;
                if(tot_slice_eng > 0){
                  std::cout << "Total energy : " << tot_slice_eng << "\n";
                  std::cout << "Total cosmic energy : " << tot_slice_cos_eng << "\n";
                  std::cout << "Total nu energy : " << tot_slice_nu_eng << "\n";
                  std::cout << "Cosmic fraction : " << double(tot_slice_cos_eng)/tot_slice_eng << "\n";
                  std::cout << "Neutrino fraction : " << double(tot_slice_nu_eng)/tot_slice_eng << "\n";
                     
                  ftotsliceNuE = tot_slice_nu_eng;
                  ftotsliceCosE = tot_slice_cos_eng;
                  ftotsliceOthE = ftotsliceE - ftotsliceNuE - ftotsliceCosE;
                  fsliceNuEfrac = double(ftotsliceNuE)/ftotsliceE;
                  fsliceCosEfrac = double(ftotsliceCosE)/ftotsliceE;
                  fsliceOthEfrac = double(ftotsliceOthE)/ftotsliceE;

                  ///////////////////////////////////////////////////// Check whether this has pandora T0 ///////////////////////////////
                  std::vector<float> pfp_t0s;

                  if(findManyPFPs.isValid()){
                    std::vector<art::Ptr<recob::PFParticle>> slice_pfps = findManyPFPs.at(slice.key());
                    if(slice_pfps.size()){
                      for(auto const &pfp : slice_pfps){
                        if(findManyT0s.isValid()){
                          std::vector<art::Ptr<anab::T0>> t0s = findManyT0s.at(pfp.key());
                          if(t0s.size()){
                            for(auto const& t0 : t0s){
                              pfp_t0s.push_back(t0->Time());
                            }
                          }
                        }
                      }
                    }
                  }

                  if(pfp_t0s.size()){
                    fT0 = *min_element(pfp_t0s.begin(), pfp_t0s.end());
                  } 

                  ////////////////////////////////////////////////////// ////////////////////////////////////////////////////////////////
                        
                  if((double(tot_slice_cos_eng)/tot_slice_eng) >= fCosEfrac){
                    interaction = kCosmic;
                  } // select cosmic
                  
                  else if((double(tot_slice_nu_eng)/tot_slice_eng) >= fNuEfrac){
                    int index = -1;
                    double min_E = 0.; 
                    std::cout << "Size of the true energy list : " << truth_energies.size() << "\n";
                           
                    for(unsigned int true_energyID=0; true_energyID<truth_energies.size(); true_energyID++){
                      if(truth_energies[true_energyID] > min_E){
                        std::cout << "MC truth energy : " << truth_energies[true_energyID] << "\n";
                        art::Ptr<simb::MCTruth> mctruth = mctruths[true_energyID];
                        simb::MCNeutrino true_neutrino = mctruth->GetNeutrino();
                        std::cout << "Neutrino vertext : " << true_neutrino.Nu().Vx() << "  " << true_neutrino.Nu().Vy() << "  " << true_neutrino.Nu().Vz() << "\n";

                        if(fUseNuContainment){
                          if(Is_in_Fiducial_Vol(true_neutrino.Nu().Vx(),true_neutrino.Nu().Vy(),true_neutrino.Nu().Vz())){ 
                            index = true_energyID;
                            min_E = truth_energies[true_energyID];
                            fNuIsContained = true;
                          }
                        }
                                    
                        else{
                          index = true_energyID;
                          min_E = truth_energies[true_energyID];
                          if(Is_in_Fiducial_Vol(true_neutrino.Nu().Vx(),true_neutrino.Nu().Vy(),true_neutrino.Nu().Vz())) fNuIsContained = true;
                          else fNuIsContained = false;
                        }
                                    
                      }
                    }

                    if(index != -1){
                      std::cout << "index diff -1" << std::endl;
                      std::cout << "GetTopology and SetTruthInfo etc" << std::endl;
                      simb::MCNeutrino true_neutrino = mctruths[index]->GetNeutrino();
                      interaction = labels.GetInteractionType(true_neutrino);
                      labels.GetTopology(mctruths[index], fTopologyHitsCut);
                      float nu_energy = true_neutrino.Nu().E();
                      float lep_energy = true_neutrino.Lepton().E();
                      fNuID = true_neutrino.Nu().TrackId();
                      info.SetTruthInfo(nu_energy, lep_energy, 0., event_weight);
                      info.SetTopologyInformation(labels.GetPDG(), labels.GetNProtons(), labels.GetNPions(), labels.GetNPizeros(), labels.GetNNeutrons(), labels.GetTopologyType(), labels.GetTopologyTypeAlt());
                      //fNuSlice_purity = double(Get_Tot_nu_E(clockData,mctruths[index],slice_hits))/tot_slice_eng; // newly deleted
                      fNuSlice_purity = double(truth_energies[index])/tot_slice_eng; // newly added
                      //if(Get_Tot_nu_E(clockData,mctruths[index],HitList)>0)fNuSlice_completeness = double(Get_Tot_nu_E(clockData,mctruths[index],slice_hits))/Get_Tot_nu_E(clockData,mctruths[index],HitList); // newly deleted
                      
                      if(TotE_from_mctruth(clockData, hits, mctruths[index])) fNuSlice_completeness = double(truth_energies[index])/TotE_from_mctruth(clockData, hits, mctruths[index]); // newly added
                      std::cout << "Nu purity : " << fNuSlice_purity << " Nu completeness : " << fNuSlice_completeness << "\n";
                    }

                  } // select neutrino
                } // total slice energy is greater that 0
              } // valid hit association

              if(findManyPFPs.isValid()){
                std::cout << "findManyPFPs.isValid" << std::endl;
                fIsSliceNu=Is_Slice_Nu(findManyPFPs,slice);

                if(fIsSliceNu){
                  std::cout << "if(fIsSliceNu)" << std::endl;
                  /*art::Handle< std::vector<recob::PFParticle> > PFPListHandle;
                  std::vector< art::Ptr<recob::PFParticle> > PFPList;
                  if( evt.getByLabel(fPFParticleModuleLabel,PFPListHandle))
                  art::fill_ptr_vector(PFPList,PFPListHandle);*/ 
                  //art::FindManyP<larpandoraobj::PFParticleMetadata> fm_pfpmd(PFPListHandle, evt, fPFParticleModuleLabel);

                  if(fm_pfpmd.isValid()){
                    fSliceScore = Get_Slice_Score(fm_pfpmd,Get_Nu_like_PFP(findManyPFPs,slice));
                    if(slices[sliceID]->ID() == Get_Best_Slice_ID_PFP_pdg(findManyPFPs,fm_pfpmd,slices)[0]) fIsbestSlice = true;
                      fbestpfppdg = Get_Best_Slice_ID_PFP_pdg(findManyPFPs,fm_pfpmd,slices)[1];
                      fbestsliceID = Get_Best_Slice_ID_PFP_pdg(findManyPFPs,fm_pfpmd,slices)[0];
                      fpfppdg = Get_Slice_PFP_ID(findManyPFPs, slice);
                  }
                }
              }

              LArTrainingNuData train(interaction, *pixelmaps[pmID], info);
              std::string evtid = "r"+std::to_string(evt.run())+"_s"+std::to_string(evt.subRun())+"_e"+std::to_string(evt.event())+"_sl"+std::to_string(pixelmaps[pmID]->fSliceID)+"_h"+std::to_string(time(0));
              write_files(train, evtid);
              break;
            } // valid pfps
          } // loop over slices
	      } // loop over pixel maps
      } // use backtrack information
         
      else{
        std::cout << "***************** Used truth level information to get neutrino information ****************\n";
        for(unsigned int pmID=0; pmID<pixelmaps.size(); pmID++){
          Clear();
          frun = evt.run();
          fsubrun = evt.subRun();
          fevent = evt.id().event();
               
          if(Nu_is_contained) fNuIsContained = true;
               
          LArTrainingNuData train(interaction, *pixelmaps[pmID], info);
          fsliceID = pixelmaps[pmID]->fSliceID;
           /*art::Handle< std::vector<recob::Slice> > slicesHandle;
           std::vector< art::Ptr<recob::Slice> > slices;
           if( evt.getByLabel(fSliceLabel,slicesHandle))
              art::fill_ptr_vector(slices,slicesHandle);
              art::FindManyP<recob::PFParticle> findManyPFPs(slicesHandle, evt, fPFParticleModuleLabel);*/
               
          for(unsigned int sliceID=0; sliceID<slices.size(); sliceID++){
            art::Ptr<recob::Slice> slice = slices[sliceID];
            const std::string &slice_tag_suff = slice_tag_suffixes[sliceID];
            
            art::FindManyP<recob::PFParticle> findManyPFPs(slices_handle, evt, fPFParticleModuleLabel + slice_tag_suff);
            art::FindManyP<recob::Hit> findManyHits(slices_handle, evt, fSliceLabel + slice_tag_suff);
            art::FindManyP<anab::T0> findManyT0s(pfps_handle, evt, fT0Label + slice_tag_suff);
            art::FindManyP<larpandoraobj::PFParticleMetadata> fm_pfpmd(pfps_handle, evt, fPFParticleModuleLabel + slice_tag_suff);

            if(slice->ID() == fsliceID){
              if(findManyHits.isValid()){
                std::vector<art::Ptr<recob::Hit>> slice_hits = findManyHits.at(slice.key());
                fNhits_total = slice_hits.size();
              }
                     
              if(findManyPFPs.isValid()){
                fIsSliceNu=Is_Slice_Nu(findManyPFPs,slice);
                ////////////////////////////////////////// pandora T0 //////////////////////////////////////////////////////
                std::vector<float> pfp_t0s;
                        
                if(findManyPFPs.isValid()){
                  std::vector<art::Ptr<recob::PFParticle>> slice_pfps = findManyPFPs.at(slice.key());
                  if(slice_pfps.size()){
                    for(auto const &pfp : slice_pfps){
                      if(findManyT0s.isValid()){
                        std::vector<art::Ptr<anab::T0>> t0s = findManyT0s.at(pfp.key());
                        if(t0s.size()){
                          for(auto const& t0 : t0s){
                            pfp_t0s.push_back(t0->Time());
                          }
                        }
                      }
                    }
                  }
                }

                if(pfp_t0s.size()){
                  fT0 = *min_element(pfp_t0s.begin(), pfp_t0s.end());
                } 

                ///////////////////////////////////////////////////////////////////////////////////////////////////////////
                        
                if(fIsSliceNu){
			       /*art::Handle< std::vector<recob::PFParticle> > PFPListHandle;
                               std::vector< art::Ptr<recob::PFParticle> > PFPList;
                               if( evt.getByLabel(fPFParticleModuleLabel,PFPListHandle))
                                   art::fill_ptr_vector(PFPList,PFPListHandle);*/ 

                  if(fm_pfpmd.isValid()){
                    fSliceScore = Get_Slice_Score(fm_pfpmd,Get_Nu_like_PFP(findManyPFPs,slice));
                    if(slice->ID() == Get_Best_Slice_ID_PFP_pdg(findManyPFPs,fm_pfpmd,slices)[0]) fIsbestSlice = true;

                    fbestpfppdg = Get_Best_Slice_ID_PFP_pdg(findManyPFPs,fm_pfpmd,slices)[1];
                    fbestsliceID = Get_Best_Slice_ID_PFP_pdg(findManyPFPs,fm_pfpmd,slices)[0];
                    fpfppdg = Get_Slice_PFP_ID(findManyPFPs, slice);
                  }
                }
              }
              break;
            }
          }

          std::string evtid = "r"+std::to_string(evt.run())+"_s"+std::to_string(evt.subRun())+"_e"+std::to_string(evt.event())+"_sl"+std::to_string(pixelmaps[pmID]->fSliceID)+"_h"+std::to_string(time(0));
          write_files(train, evtid);
        }
      } // use truth info
    } // use slice to make images

    else{
      std::cout << "***************** Using entire event ****************\n";
      LArTrainingNuData train(interaction, *pixelmaps[0], info);
      std::string evtid = "r"+std::to_string(evt.run())+"_s"+std::to_string(evt.subRun())+"_e"+std::to_string(evt.event())+"_h"+std::to_string(time(0));
      write_files(train, evtid);
    } // use event to make images
  }

  //......................................................................
  
  void CVNZlibMakerICARUS::write_files(LArTrainingNuData td, std::string evtid)
  {
    // cropped from 2880 x 500 to 500 x 500 here 
    std::vector<unsigned char> pixel_array(3 * fPlaneLimit * fTDCLimit);
    //fImage.DisableRegionSelection();
    fImage.SetPixelMapSize(td.fPMap.NWire(), td.fPMap.NTdc());
    //fImage.DisableRegionSelection();
    //fImage.EnableRegionSelection();
    //std::cout << "Number of wires : " << td.fPMap.NWire() << "  Number of time ticks : " << td.fPMap.NTdc() << "\n";
    //std::cout << "Number of pixelmap wires : " << fImage.Get_input_pmap_wire_count() << "  Numbr of pixel map tdcs : " << fImage.Get_input_pmap_tdc_count() << "\n";
    fImage.ConvertPixelMapToPixelArray(td.fPMap, pixel_array);

    ulong src_len = 3 * fPlaneLimit * fTDCLimit; // pixelArray length
    ulong dest_len = compressBound(src_len);     // calculate size of the compressed data               
    char* ostream = (char *) malloc(dest_len);  // allocate memory for the compressed data

    int res = compress((Bytef *) ostream, &dest_len, (Bytef *) &pixel_array[0], src_len);
    
    //std::cout << "src_len : " << src_len << " dest_len : " << dest_len << "\n";

    // Buffer error

    if (res == Z_BUF_ERROR)
      std::cout << "Buffer too small!" << std::endl;

    // Memory error
    else if (res ==  Z_MEM_ERROR)
      std::cout << "Not enough memory for compression!" << std::endl;

    // Compression ok 
    else {
      //std::cout << "*************** Creating output files ***************\n";
      //std::cout << "out_dir is : " << out_dir << "\n";
      // Create output files 
      std::string image_file_name = out_dir + "/event_" + evtid + ".gz";
      std::string info_file_name = out_dir + "/event_" +  evtid + ".info";

      std::ofstream image_file (image_file_name, std::ofstream::binary);
      std::ofstream info_file  (info_file_name);

      if(image_file.is_open() && info_file.is_open()) {

        // Write compressed data to file

        image_file.write(ostream, dest_len);

        image_file.close(); // close file

        // Write records to file

        // Category
	
	info_file << td.fInt << std::endl;
        //info_file << td.fInfo << std::endl;
	info_file << td.fInfo;
        info_file << td.fPMap.GetTotHits() << "," << fNhits_total <<  std::endl;   
	info_file << frun << "," << fsubrun << "," << fevent << "," << fsliceID << std::endl;   
	info_file << ftotsliceE << "," << ftotsliceNuE << "," << ftotsliceCosE << "," << ftotsliceOthE << "," << fsliceNuEfrac << "," << fsliceCosEfrac << "," << fsliceOthEfrac << "," << fNuSlice_purity << "," << fNuSlice_completeness <<  std::endl;
	info_file << fIsSliceNu << "," << fSliceScore << "," << fIsbestSlice << "," << fbestpfppdg << "," << fbestsliceID << "," << fpfppdg << std::endl;
	info_file << fNuIsContained << "," << fNuID << "  " << fT0 << std::endl;
	info_file.close(); // close file
      }
      else {

        if (image_file.is_open())
          image_file.close();
        else 
          throw art::Exception(art::errors::FileOpenError)
            << "Unable to open file " << image_file_name << "!" << std::endl;

        if (info_file.is_open())
          info_file.close();
        else
          throw art::Exception(art::errors::FileOpenError)
            << "Unable to open file " << info_file_name << "!" << std::endl;
      }
    }
    
    free(ostream);  // free allocated memory

  } // lcvn::CVNZlibMaker::write_files
  

//////////////////////////////////////////////////////////////////////////////////////////////////////////

bool CVNZlibMakerICARUS::Is_in_Fiducial_Vol(double stX, double stY, double stZ)
{
     double X_bound = 196.5 - fVolCut;
     double Y_bound = 200. - fVolCut;
     double Z_up_bound = 500. - fVolCut;
     double Z_low_bound = fVolCut;
	
     if(TMath::Abs(stX) <= X_bound && TMath::Abs(stY) <= Y_bound && (stZ >=Z_low_bound && stZ <= Z_up_bound)) return true;
     else return false;
}  

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void  CVNZlibMakerICARUS::HitTruth( detinfo::DetectorClocksData const& clockData, art::Ptr<recob::Hit> const& hit, Int_t& truthid){
      std::vector<sim::TrackIDE> trackIDEs = bt_serv->HitToTrackIDEs(clockData, hit);
      if( !trackIDEs.size() ) return;
      Float_t maxe = 0;
      Int_t bestid = 0;
      for(size_t i = 0; i < trackIDEs.size(); ++i){
          if( trackIDEs[i].energy > maxe ) {
              maxe = trackIDEs[i].energy;
              bestid = trackIDEs[i].trackID;
          }
      }
      truthid = bestid;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool  CVNZlibMakerICARUS::HitTruthId(detinfo::DetectorClocksData const& clockData, art::Ptr<recob::Hit> const& hit, Int_t& mcid){
      mcid = std::numeric_limits<int>::lowest();
      HitTruth(clockData,hit,mcid);
      if( mcid > std::numeric_limits<int>::lowest() ) return true;
      else return false;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool  CVNZlibMakerICARUS::TrackIdToMCTruth(Int_t const trkID, art::Ptr<simb::MCTruth>& mctruth){
      bool matchFound = false;
      try {
          mctruth = pi_serv->TrackIdToMCTruth_P(trkID);
          matchFound = true;
      } catch(...) {
      std::cout<<"Exception thrown matching TrackID "<<trkID<<" to MCTruth\n";
      matchFound = false;
      }
      return matchFound;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

double CVNZlibMakerICARUS::HitEfrmTrkID(detinfo::DetectorClocksData const& clockData, art::Ptr<recob::Hit> const& hit, Int_t truthid){
       double hit_E = 0.;
       std::vector<sim::TrackIDE> trackIDEs = bt_serv->HitToTrackIDEs(clockData, hit);
       if( !trackIDEs.size() ) return hit_E;
       for(size_t i = 0; i < trackIDEs.size(); ++i){
           if(trackIDEs[i].trackID == truthid) hit_E += trackIDEs[i].energy;
       }
       return hit_E;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////

double CVNZlibMakerICARUS::HitTotE(detinfo::DetectorClocksData const& clockData, art::Ptr<recob::Hit> const& hit){
       double Tot_E = 0.;
       std::vector<sim::TrackIDE> trackIDEs = bt_serv->HitToTrackIDEs(clockData, hit);
       if( !trackIDEs.size() ) return Tot_E;
       for(size_t i = 0; i < trackIDEs.size(); ++i){
           Tot_E += trackIDEs[i].energy;
       }
       return Tot_E;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool CVNZlibMakerICARUS::Is_Slice_Nu(art::FindManyP<recob::PFParticle> const& assn, art::Ptr<recob::Slice> const& slce){
     if(assn.isValid()){
        std::vector<art::Ptr<recob::PFParticle>> slicePFPs = assn.at(slce.key());
	for(auto const &pfp : slicePFPs){
	    if((pfp->PdgCode() == 12 || pfp->PdgCode() == 14)) return true;
	}
     }
     return false;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double CVNZlibMakerICARUS::Get_Slice_Score(art::FindManyP<larpandoraobj::PFParticleMetadata> const& assn, art::Ptr<recob::PFParticle> const& pfp){
       double best_score = -9999;
       if(assn.isValid() && pfp.isNonnull()){
          const std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> pfpMetaVec = assn.at(pfp->Self());
	  if(pfpMetaVec.size()){
	     for (auto const pfpMeta : pfpMetaVec){
	          larpandoraobj::PFParticleMetadata::PropertiesMap propertiesMap = pfpMeta->GetPropertiesMap();
		  double score = propertiesMap.at("NuScore");
		  if (score > best_score) best_score = score;
	     }
	  }
       }
       return best_score;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
art::Ptr<recob::PFParticle> CVNZlibMakerICARUS::Get_Nu_like_PFP(art::FindManyP<recob::PFParticle> const& assn, art::Ptr<recob::Slice> const& slce){
	art::Ptr<recob::PFParticle> result;
	std::vector<art::Ptr<recob::PFParticle>> slicePFPs = assn.at(slce.key());
	for(auto const &pfp : slicePFPs){
	    if((pfp->PdgCode() == 12 || pfp->PdgCode() == 14)){ 
		result = pfp;
		break;
	    }
	}
	return result;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<int> CVNZlibMakerICARUS::Get_Best_Slice_ID_PFP_pdg(art::FindManyP<recob::PFParticle> const& pfp_assn, art::FindManyP<larpandoraobj::PFParticleMetadata> const& metadata_assn,
		std::vector<art::Ptr<recob::Slice>> const& slice_vec){
    std::vector<int> Sl_ID_pfp_pdf = {-999,-999};
    std::vector<art::Ptr<recob::PFParticle>> nuPfpVec;
    std::vector<art::Ptr<recob::Slice>> nuSliceVec;
    for(auto const &slice : slice_vec){
       std::vector<art::Ptr<recob::PFParticle>> slicePFPs = pfp_assn.at(slice.key());
       for (auto const &pfp : slicePFPs){
            if ((pfp->PdgCode() == 12 || pfp->PdgCode() == 14)){
	         nuPfpVec.push_back(pfp);
                 nuSliceVec.push_back(slice);
		 break;
	    }
       }
    }
    
    if(nuPfpVec.size()){
       float bestScore = -999.;
       art::Ptr<recob::Slice> bestNuSlice;
       art::Ptr<recob::PFParticle> bestNuPfp;
       for(size_t i=0; i<nuPfpVec.size(); i++){
           float score = -999.;
           const std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> pfpMetaVec = metadata_assn.at(nuPfpVec.at(i)->Self());
           for (auto const pfpMeta : pfpMetaVec) {
                larpandoraobj::PFParticleMetadata::PropertiesMap propertiesMap = pfpMeta->GetPropertiesMap();
                score = propertiesMap.at("NuScore");
                if( score > bestScore ) {
                    bestScore = score;
                    bestNuPfp = nuPfpVec.at(i);
                    bestNuSlice = nuSliceVec.at(i);
                }
           }
       }
       
       Sl_ID_pfp_pdf[0] = bestNuSlice->ID();
       Sl_ID_pfp_pdf[1] = bestNuPfp->PdgCode();
    }
    return Sl_ID_pfp_pdf;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int CVNZlibMakerICARUS::Get_Slice_PFP_ID(art::FindManyP<recob::PFParticle> const& assn, art::Ptr<recob::Slice> const& slce){
    std::vector<art::Ptr<recob::PFParticle>> slicePFPs = assn.at(slce.key());
    int result = -9999;
    for(auto const &pfp : slicePFPs){
	    if((pfp->PdgCode() == 12 || pfp->PdgCode() == 14)){ 
		result = pfp->PdgCode();
		break;
	    }
	}
    return result;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int CVNZlibMakerICARUS::Get_Hit_Count(unsigned int tpc_no, unsigned int pln_no, std::vector<art::Ptr<recob::Hit>> const& hit_vec){
    int n_hits = 0;
    if(hit_vec.size()){
       for(auto const hit : hit_vec){
           if(hit->WireID().TPC == tpc_no){
	      if(hit->WireID().Plane == pln_no){
	         n_hits++;
	      }
	   }
       }
    }
    
    return n_hits;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////

double CVNZlibMakerICARUS::Get_Tot_nu_E(detinfo::DetectorClocksData const& clockData,art::Ptr<simb::MCTruth> const& mcneutrino, std::vector<art::Ptr<recob::Hit>> const& hit_vec){
       double tot_nu_E = 0;
       for(auto const hit : hit_vec){
           Int_t trkId;
	   if(HitTruthId(clockData,hit,trkId)){
	      art::Ptr<simb::MCTruth> mctruth;
	      if(TrackIdToMCTruth(trkId,mctruth)){
	         if(mctruth->Origin() == simb::kBeamNeutrino){
		    if(mctruth == mcneutrino){
		       tot_nu_E += HitEfrmTrkID(clockData,hit,trkId);
		    }
		 }
	      }
	   }
       }
       return tot_nu_E;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////

double CVNZlibMakerICARUS::HitTotE_frm_given_orgin(detinfo::DetectorClocksData const& clockData, art::Ptr<recob::Hit> const& hit, std::string origin_name){
       double Tot_E = 0.;
       std::vector<sim::TrackIDE> trackIDEs = bt_serv->HitToTrackIDEs(clockData, hit);
       if( !trackIDEs.size() ) return Tot_E;
       for(size_t i = 0; i < trackIDEs.size(); ++i){
	   art::Ptr<simb::MCTruth> mctruth;
	   if(TrackIdToMCTruth(trackIDEs[i].trackID,mctruth)){
	      if(origin_name == "nu"){
	         if(mctruth->Origin() == simb::kBeamNeutrino) Tot_E += trackIDEs[i].energy;
	      }
	      else if(origin_name == "cos"){
	         if(mctruth->Origin() == simb::kCosmicRay) Tot_E += trackIDEs[i].energy;
	      }
	   }
           
       }
       return Tot_E;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double CVNZlibMakerICARUS::HitTotE_frm_mctruth(detinfo::DetectorClocksData const& clockData, art::Ptr<recob::Hit> const& hit, art::Ptr<simb::MCTruth> const& my_mctruth){
       double Tot_E = 0.;
       std::vector<sim::TrackIDE> trackIDEs = bt_serv->HitToTrackIDEs(clockData, hit);
       if( !trackIDEs.size() ) return Tot_E;
       for(size_t i = 0; i < trackIDEs.size(); ++i){
           art::Ptr<simb::MCTruth> mctruth;
	   if(TrackIdToMCTruth(trackIDEs[i].trackID,mctruth)){
	      if(my_mctruth.isNonnull()){
	         if(mctruth == my_mctruth){
		    Tot_E += trackIDEs[i].energy;
		 }
	      }
	   }
       }
       return Tot_E;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////

double CVNZlibMakerICARUS::TotE_from_mctruth(detinfo::DetectorClocksData const& clockData, std::vector<art::Ptr<recob::Hit>> const& hit_vec, art::Ptr<simb::MCTruth> const& my_mctruth){
       double tot_E = 0;
       if(hit_vec.size() && my_mctruth.isNonnull()){
          for(auto const hit : hit_vec){
	      tot_E += HitTotE_frm_mctruth(clockData,hit,my_mctruth);
          }
       }
       return tot_E;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////

void CVNZlibMakerICARUS::Clear(){
     frun = -9999;
     fsubrun = -9999;
     fevent = -9999;
     fsliceID = -9999;
     ftotsliceE = -9999;
     ftotsliceNuE = -9999;
     ftotsliceCosE = -9999;
     ftotsliceOthE = -9999;
     fsliceNuEfrac = -9999;
     fsliceCosEfrac = -9999;
     fsliceOthEfrac = -9999;	    
     fIsSliceNu = false;
     fSliceScore = -9999;
     fIsbestSlice = false;
     fbestpfppdg = -9999;
     fbestsliceID = -9999;
     fpfppdg = -9999;
     fNuIsContained = false;
     fNuID = -9999;
     fNhits_total = -9999;
     fNuSlice_purity = -9999.;
     fNuSlice_completeness = -9999.;
     fT0 = 0.;
}

} // namespace lcvn

DEFINE_ART_MODULE(lcvn::CVNZlibMakerICARUS)
