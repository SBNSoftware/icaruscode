////////////////////////////////////////////////////////////////////////
// Class:       DemoRecoAna
// Plugin Type: analyzer
// File:        DemoRecoAna_module.cc
//
// Generated at Thu Nov 11 13:28:47 2021 by Gray Putnam using cetskelgen
// from  version .
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

#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Shower.h"

namespace icarus {
  class DemoRecoAna;
}


class icarus::DemoRecoAna : public art::EDAnalyzer {
public:
  explicit DemoRecoAna(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  DemoRecoAna(DemoRecoAna const&) = delete;
  DemoRecoAna(DemoRecoAna&&) = delete;
  DemoRecoAna& operator=(DemoRecoAna const&) = delete;
  DemoRecoAna& operator=(DemoRecoAna&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:

  art::InputTag fPandoraLabel;
  art::InputTag fTrackLabel;

};


icarus::DemoRecoAna::DemoRecoAna(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    fPandoraLabel(p.get<art::InputTag>("PandoraLabel", "pandoraGausCryoE")),
    fTrackLabel(p.get<art::InputTag>("TrackLabel", "pandoraTrackGausCryoE"))
  // More initializers here.
{
}

void icarus::DemoRecoAna::analyze(art::Event const& e)
{
  std::cout << "Processing Event: " << e.event() << " in run: " << e.run() << " with tag: " << fPandoraLabel << std::endl;

  // Get the list of slices in the event
  std::vector<art::Ptr<recob::Slice>> slices;
  art::ValidHandle<std::vector<recob::Slice>> slice_handle = e.getValidHandle<std::vector<recob::Slice>>(fPandoraLabel);
  art::fill_ptr_vector(slices, slice_handle);

  // Lookup the mapping of Slices -> PFParticle's
  art::FindManyP<recob::PFParticle> fmSliceParticles(slices, e, fPandoraLabel);

  // Iterate over the slices
  for (unsigned i_slc = 0; i_slc < slices.size(); i_slc++) {
    std::cout << "Processing Slice: " << i_slc << std::endl;

    // Get the list of Particles in this slice
    const std::vector<art::Ptr<recob::PFParticle>> &slice_pfps = fmSliceParticles.at(i_slc);

    // Lookup the "Primary" PFParticle
    int primary_pfp_index = -1;
    for (unsigned i_pfp = 0; i_pfp < slice_pfps.size(); i_pfp++) {
      if (slice_pfps[i_pfp]->IsPrimary()) { // Found it!
        primary_pfp_index = i_pfp;
        break;
      }
    }

    // Didn't find a Primary PFParticle: bad slice, ignore it
    if (primary_pfp_index < 0) continue;

    // Let's look at the Primary PFParticle
    const recob::PFParticle &primary = *slice_pfps[primary_pfp_index];

    std::cout << "Primary PFParticle! ID: " << primary.Self() << " PDG: " << primary.PdgCode() << " Parent: " << primary.Parent() << std::endl;

    // The "PdgCode" of the PFParticle tells us what kind of slice this is
    //
    // 13, 11 -- Cosmic (denoted "Clear" or "Obvious" cosmic)
    // 12/14/16 -- Neutrino (Note that the particular Pdg of the neutrino is not reliable
    bool is_cosmic = abs(primary.PdgCode()) == 13 || abs(primary.PdgCode()) == 11;

    if (is_cosmic) {
      std::cout << "Cosmic Slice\n";
    }
    else {
      std::cout << "Neutrino Slice\n";
    }

    // Let's go one level deeper into the reconstruction: the Tracks
    //
    // Lookup all the tracks in the Slice and print their lengths
    art::FindManyP<recob::Track> fmPFPTracks(slice_pfps, e, fTrackLabel);

    for (unsigned i_pfp = 0; i_pfp < slice_pfps.size(); i_pfp ++) {
      // Lookup the Track associated with the PFP
      //
      // There may not be a track associated with the PFP (if it is a
      // "Neutrino" or "Shower" PFP). In this case the vector 
      // will be empty
      const std::vector<art::Ptr<recob::Track>> maybe_track = fmPFPTracks.at(i_pfp);

      if (maybe_track.empty()) continue; // Not a Track

      const recob::Track &track = *maybe_track[0];

      std::cout << "Track ID: " << track.ID() << " len: " << track.Length() << std::endl; 

      // Now that you have the track, you can pull in further information like:
      //   -Calorimetry
      //   -Particle ID
      //   -CRT matching
      //   -Momentum reconstruction
      //
      // Provided your input file has the necessary data products
    }

  }
}

DEFINE_ART_MODULE(icarus::DemoRecoAna)
