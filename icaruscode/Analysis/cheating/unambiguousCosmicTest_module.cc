////////////////////////////////////////////////////////////////////////
// Class:       unambiguousCosmicTest
// Plugin Type: analyzer (Unknown Unknown)
// File:        unambiguousCosmicTest_module.cc
//
// Generated at Thu Feb 20 02:28:49 2025 by Mattia Sotgia using cetskelgen
// from cetlib version 3.18.02.
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

#include "art_root_io/TFileService.h"
#include <TTree.h>

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include "canvas/Persistency/Common/FindManyP.h"

#include "icaruscode/RecoUtils/RecoUtils.h"

namespace icaruscode {
  namespace pandoraCheating {
    class unambiguousCosmicTest;
  }
}


namespace format_utils {
  std::string space4 = "    ";
  std::string space8 = space4 + space4;
  std::string space12 = space4 + space4 + space4;
  std::string space16 = space4 + space4 + space4 + space4;
}

class icaruscode::pandoraCheating::unambiguousCosmicTest : public art::EDAnalyzer {
public:
  explicit unambiguousCosmicTest(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  unambiguousCosmicTest(unambiguousCosmicTest const&) = delete;
  unambiguousCosmicTest(unambiguousCosmicTest&&) = delete;
  unambiguousCosmicTest& operator=(unambiguousCosmicTest const&) = delete;
  unambiguousCosmicTest& operator=(unambiguousCosmicTest&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Output TTree + branches
  TTree* fTree;
  double fNoOfSlice, fTotalUnambiguousCosmic, fEventID;
  double fNoOfSelectedEvents, fNoOfUCSelected;


  // FHiCL settings;
  std::string fTruthLabel, fSliceLabel, fPFParticleLabel, fTrackLabel, fPFParticleMetadataLabel;
  bool fVerbose;
  std::vector<std::string> fPandoraSuffixes;

};


icaruscode::pandoraCheating::unambiguousCosmicTest::unambiguousCosmicTest(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  fTruthLabel(p.get<std::string>("TruthLabel")),
  fSliceLabel(p.get<std::string>("SliceLabel")),
  fPFParticleLabel(p.get<std::string>("PFParticleLabel")),
  fTrackLabel(p.get<std::string>("TrackLabel")),
  fPFParticleMetadataLabel(p.get<std::string>("PFParticleMetadataLabel")),
  fVerbose(p.get<bool>("VerboseOutput")),
  fPandoraSuffixes(p.get<std::vector<std::string>>("PandoraSuffixes"))
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void icaruscode::pandoraCheating::unambiguousCosmicTest::analyze(art::Event const& e)
{
  if (fVerbose)
    std::cout << "--> Event " << e.id().event() << std::endl;

  fEventID = e.id().event();
  fNoOfSlice = 0;
  fTotalUnambiguousCosmic = 0;

  // Check whether the simulation gives a true nuCC event..
  art::ValidHandle<std::vector<simb::MCTruth>> truthHandle =
    e.getValidHandle<std::vector<simb::MCTruth>>(fTruthLabel); 
  std::vector<art::Ptr<simb::MCTruth>> mc_truths;

  if (truthHandle.isValid())
    art::fill_ptr_vector(mc_truths, truthHandle);
  
  int pdg = -1;
  bool iscc = false;
  for (const art::Ptr<simb::MCTruth> &mc: mc_truths) {
    pdg = mc->GetNeutrino().Nu().PdgCode();
    iscc = mc->GetNeutrino().CCNC() == simb::kCC;

    if (std::abs(pdg) == 14 && iscc) {
      if (fVerbose)
        std::cout << format_utils::space4 << "--> True numuCC" << std::endl;
      break;
    } 
  }
  
  if (std::abs(pdg) != 14 || !iscc) {
    if (fVerbose)
      std::cout << format_utils::space4 << "--> No true numuCC in event" << std::endl;
    return;
  }
  
  std::vector<art::Ptr<recob::Slice>> slices;
  std::vector<std::string> slice_tag_suffixes;
  std::vector<unsigned>    slice_tag_id;

  // Collect all slices for TPCs
  for (unsigned tag_i = 0; tag_i<fPandoraSuffixes.size(); tag_i++) {
    const std::string &pandora_tag_suffix = fPandoraSuffixes[tag_i];

    art::ValidHandle<std::vector<recob::Slice>> thisSlice = 
      e.getValidHandle<std::vector<recob::Slice>>(fSliceLabel + pandora_tag_suffix);

    if (thisSlice.isValid()) 
      art::fill_ptr_vector(slices, thisSlice);
    
    for (unsigned i = 0; i < thisSlice->size(); i++)
      slice_tag_suffixes.push_back(pandora_tag_suffix);
    
  }

  if (slices.size() == 0) {
    if (fVerbose)
      std::cout << format_utils::space8 << "--> No slice(s) reco found" << std::endl;
    return;
  }


  // Looping over slices
  for (unsigned sliceID = 0; sliceID < slices.size(); sliceID++) {
    
    art::Ptr<recob::Slice> slice = slices[sliceID];
    const std::string &slice_tag_suff = slice_tag_suffixes[sliceID];


    if (fVerbose)
      std::cout << format_utils::space8 << "--> Slice " << slice.key() 
                << ", in tag=pandora" << slice_tag_suff << std::endl;
    
    
    std::vector<art::Ptr<recob::Slice>> sliceList {slice};

    art::ValidHandle<std::vector<recob::PFParticle>> pfpHandle = 
      e.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleLabel + slice_tag_suff);
    art::FindManyP<recob::PFParticle> 
      slicePfpAssns(sliceList, e, fPFParticleLabel + slice_tag_suff);
    art::FindManyP<larpandoraobj::PFParticleMetadata> 
      pfpMetadataAssns(pfpHandle, e, fPFParticleMetadataLabel + slice_tag_suff);
    art::FindManyP<recob::Track>
      pfpTrackAssns(pfpHandle, e, fTrackLabel + slice_tag_suff);

    // Every slice can be both a nu slice (those should all be nu slices) or
    // a cosmic (IsClearCosmic or Cosmic from second cosmic stage)

    // Get the list of PFPs
    std::vector<art::Ptr<recob::PFParticle>> slicePfps = slicePfpAssns.at(0);

    if (fVerbose) 
      std::cout << format_utils::space12 
                << "This slice has # " << slicePfps.size() << " pfp(s) " 
                << std::endl;

    int unambiguous_cosmics_for_slice = 0;
    bool slice_is_neutrino = false;

    for (const art::Ptr<recob::PFParticle> &pfp: slicePfps) {
      const bool is_primary(pfp->IsPrimary());
      const bool is_neutrino(std::abs(pfp->PdgCode())==14); // Looking for µ neutrinos

      if(!slice_is_neutrino)
        slice_is_neutrino = is_neutrino;

      const larpandoraobj::PFParticleMetadata *pfp_meta = pfpMetadataAssns.at(pfp.key()).at(0).get();
      auto const &properties = pfp_meta->GetPropertiesMap();

      if (properties.count("IsClearCosmic")) {
        unambiguous_cosmics_for_slice++;

        if (fVerbose && is_primary) 
          std::cout << format_utils::space12
                    << "The primary for this slice IsClearCosmic" 
                    << " (this slice will we classified as is_clear_cosmic)" 
                    << std::endl;
        if (fVerbose)
          std::cout << format_utils::space12
                    << "This is a clear cosmic with pdg = " 
                    << pfp->PdgCode() 
                    << std::endl;
      }

      std::vector<art::Ptr<recob::Track>> tracks = pfpTrackAssns.at(pfp.key());

      if (fVerbose)
        std::cout << format_utils::space16 
                  << "Found a total of " << tracks.size() 
                  << " track(s) for this pfp" << std::endl;
                    
      for (const art::Ptr<recob::Track> &track: tracks) {
        if (fVerbose)
          std::cout << format_utils::space16 << format_utils::space4
                    << "-->  starting from (" 
                    << track->Start().X() << ", "
                    << track->Start().Y() << ", "
                    << track->Start().Z() << ") and ending in ("
                    << track->End().X() << ", "
                    << track->End().Y() << ", "
                    << track->End().Z() << ")"
                    <<std::endl;
      }
    } // pfps loop

    if (unambiguous_cosmics_for_slice>0) {
      fTotalUnambiguousCosmic++;
      fNoOfUCSelected += fTotalUnambiguousCosmic;
    }
    
    if (fVerbose) 
      std::cout << format_utils::space16 
                << "This slice " << ((unambiguous_cosmics_for_slice>0) ? "is" : "is NOT")
                << " unambiguous cosmic. " 
                << " Is this Neutrino? " 
                << (slice_is_neutrino ? "yes" : "NO")
                << ". Found # pfp = " 
                << slicePfps.size() 
                << " for this slice "
                << std::endl;
  } // slices loop

  fNoOfSlice = slices.size();
  fNoOfSelectedEvents += slices.size();
  fTree->Fill();

  /////////// OLD IMPLEMENTATION ///////////

  // // loops over two cryostats FUCKING ICARUS
  // for (auto const &pandora_suffix: fPandoraSuffixes) {

  //   if (fVerbose)
  //     std::cout << format_utils::space8 << "--> TPC pandora" << pandora_suffix << std::endl;

  //   art::ValidHandle<std::vector<recob::Slice>> sliceHandle = 
  //     e.getValidHandle<std::vector<recob::Slice>>(fSliceLabel + pandora_suffix);
  //   art::ValidHandle<std::vector<recob::PFParticle>> pfpHandle = 
  //     e.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleLabel + pandora_suffix);
  //   // art::ValidHandle<std::vector<recob::Track>> trackHandle = 
  //   //   e.getValidHandle<std::vector<recob::Track>>(fTrackLabel + pandora_suffix);
    
  //   art::FindManyP<recob::PFParticle> 
  //     slicePfpAssns(sliceHandle, e, fPFParticleLabel + pandora_suffix);
  //   art::FindManyP<larpandoraobj::PFParticleMetadata> 
  //     pfpMetadataAssns(pfpHandle, e, fPFParticleMetadataLabel + pandora_suffix);
  //   art::FindManyP<recob::Track>
  //     pfpTrackAssns(pfpHandle, e, fTrackLabel + pandora_suffix);
    
  //   std::vector<art::Ptr<recob::Slice>> slices;

  //   if (sliceHandle.isValid())
  //     art::fill_ptr_vector(slices, sliceHandle);

  //   if (slices.size() == 0) {
  //     if (fVerbose)
  //       std::cout << format_utils::space12 << "--> No slice(s) reco found" << std::endl;
  //     continue;
  //   }

  //   for (const art::Ptr<recob::Slice> &slice: slices) {
  //     // Every slice can be both a nu slice (those should all be nu slices) or
  //     // a cosmic (IsClearCosmic or Cosmic from second cosmic stage)

  //     if (fVerbose)
  //       std::cout << format_utils::space12 << "--> Slice " << slice.key() << std::endl;
      
  //     // 1. If the slice is neutrino, then check whether it has some PFP
  //     //    identified as clear cosmics
  //     // 2. If not neutrino slice (odd...) check IsClearCosmic anyway...
      
  //     std::vector<art::Ptr<recob::PFParticle>> slicePfps = slicePfpAssns.at(slice.key());

  //     int unambiguous_cosmics_for_slice = 0;
  //     bool slice_is_neutrino = false;
  //     // bool slice_is_primary = false;

  //     for (const art::Ptr<recob::PFParticle> &pfp: slicePfps) {
  //       const bool is_primary(pfp->IsPrimary());
  //       const bool is_neutrino(std::abs(pfp->PdgCode()) == 14);
  //       if(!slice_is_neutrino)
  //         slice_is_neutrino = is_neutrino;
  //       // slice_is_primary = is_primary;

  //       // Also look at 
  //       const larpandoraobj::PFParticleMetadata *pfp_meta = pfpMetadataAssns.at(pfp.key()).at(0).get();
  //       auto const &properties = pfp_meta->GetPropertiesMap();


  //       if (properties.count("IsClearCosmic")) {
  //         unambiguous_cosmics_for_slice++;
          
  //         if (fVerbose && is_primary) 
  //           std::cout << format_utils::space16 
  //                     << "The primary for this slice IsClearCosmic" 
  //                     << " (this slice will we classified as is_clear_cosmic)" 
  //                     << std::endl;
  //         if (fVerbose)
  //           std::cout << format_utils::space16 
  //                     << "This is a clear cosmic with pdg = " 
  //                     << pfp->PdgCode() 
  //                     << std::endl;
  //       }

  //       std::vector<art::Ptr<recob::Track>> tracks = pfpTrackAssns.at(pfp.key());

  //       if (fVerbose)
  //         std::cout << format_utils::space16 
  //                   << "Found a total of " << tracks.size() 
  //                   << " track(s) for this pfp" << std::endl;
                    
  //       for (const art::Ptr<recob::Track> &track: tracks) {
  //         if (fVerbose)
  //           std::cout << format_utils::space16 << format_utils::space4
  //                     << "-->  starting from (" 
  //                     << track->Start().X() << ", "
  //                     << track->Start().Y() << ", "
  //                     << track->Start().Z() << ") and ending in ("
  //                     << track->End().X() << ", "
  //                     << track->End().Y() << ", "
  //                     << track->End().Z() << ")"
  //                     <<std::endl;
  //       }
  //     } // single pfp
      
  //     if (unambiguous_cosmics_for_slice>0)
  //       fTotalUnambiguousCosmic++;
      
  //     if (fVerbose) 
  //       std::cout << format_utils::space16 
  //                 << "This slice " << ((unambiguous_cosmics_for_slice>0) ? "is" : "is NOT")
  //                 << " unambiguous cosmic. " 
  //                 << " Is this Neutrino? " 
  //                 << (slice_is_neutrino ? "yes" : "NO")
  //                 << ". Found # pfp = " 
  //                 << slicePfps.size() 
  //                 << " for this slice "
  //                 << std::endl;
  //   } // single slice
  //   fNoOfSlice = slices.size();
  //   fTree->Fill(); 
  // } // loop pandora_suffixes
}

void icaruscode::pandoraCheating::unambiguousCosmicTest::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("pandora_clearCosmicProblem", "How many clear_cosmic in nu-only run");
  fTree->Branch("eventID", &fEventID);
  fTree->Branch("noOfUnambiguousCosmics", &fTotalUnambiguousCosmic);
  fTree->Branch("noOfSlices", &fNoOfSlice);

  fNoOfSelectedEvents = 0;
  fNoOfUCSelected = 0;
}

void icaruscode::pandoraCheating::unambiguousCosmicTest::endJob()
{
  // Implementation of optional member function here.

  if (fVerbose) 
    std::cout << "Found a total of # cosmic slices = " << fNoOfUCSelected 
              << " over a sampe of # slices = " << fNoOfSelectedEvents
              << ", corresponding to a fraction of " << fNoOfUCSelected/fNoOfSelectedEvents*100 
              << "% " << std::endl;
}

DEFINE_ART_MODULE(icaruscode::pandoraCheating::unambiguousCosmicTest)
