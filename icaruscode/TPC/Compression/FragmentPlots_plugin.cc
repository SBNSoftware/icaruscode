//std includes
#include <vector>

//ROOT includes
#include "TFile.h"
#include "TH1.h"

//Framework includes
#include "art/Framework/Core/ResultsProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Results.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Utilities/InputTag.h"

//icarus includes
#include "icaruscode/TPC/Compression/PhysCrateCompressedFragment.cc"
#include "icaruscode/Decode/ChannelMapping/IICARUSChannelMap.h"

//sbndaq includes
#include "sbndaq-artdaq-core/Overlays/ICARUS/PhysCrateFragment.cc"

namespace tcpCompression {
  class FragmentPlots : public art::ResultsProducer {
  public:
    explicit FragmentPlots(fhicl::ParameterSet const& pset);
    ~FragmentPlots() override = default;
    void event(art::Event const& evt) override; 
    void reconfigure(fhicl::ParameterSet const& p);
    void writeResults(art::Results& r) override;
    void clear() override;

  private:
    art::InputTag fFragmentsLabel;
    bool fInitiallyCompressed;
    std::unique_ptr<icarusDB::IICARUSChannelMap> fChannelMap;
    std::unique_ptr<TH1D> adcDiffs;
    std::unique_ptr<TH1D> fragSize_compressed;
    std::unique_ptr<TH1D> fragSize_uncompressed;
    std::unique_ptr<TH3D> randReads_compressed;
    std::unique_ptr<TH3D> randReads_uncompressed;
    std::unique_ptr<TH1D> fullReads_compressed;
    std::unique_ptr<TH1D> fullReads_uncompressed;
  };// end FragmentPlots class

  //------------------------------------------------------------------
  FragmentPlots::FragmentPlots(fhicl::ParameterSet const& pset)
  {
    this->reconfigure(pset);
  }

  //------------------------------------------------------------------
  void FragmentPlots::reconfigure(fhicl::ParameterSet const& pset)
  {
    fFragmentsLabel    = pset.get<art::InputTag>("FragmentsLabel"   , "daq:PHYSCRATEDATA");
    fCheckOldFragments = pset.get<bool>         ("CheckOldFragments", false);
    fDumpADCs          = pset.get<bool>         ("DumpADCs"         , false);

    art::ServiceHandle<art::TFileService> tfs;
    fChannelMap.reset(art::ServiceHandle<icarusDB::IICARUSChannelMap const>{}.get());
    
    std::string fragSize_name = "fragSize";
    std::string fragSize_title = ";Fragment Size (mb);Fragments";
    int nBins = 100;
    double binMin = 0;
    double binMax = 100;
    fragSize_compressed  .reset(tfs->make<TH1D>((fragSize_name + "_compressed").c_str(), fragSize_title, nBins, binMin, binMax));
    fragSize_uncompressed.reset(tfs->make<TH1D>((fragSize_name+"_uncompressed").c_str(), fragSize_title, nBins, binMin, binMax));
  }

  //------------------------------------------------------------------
  void FragmentPlots::event(art::Event const& evt)
  {
  }

  //------------------------------------------------------------------
  void FragmentPlots::writeResults(art::Results& r)
  {
  }

  //------------------------------------------------------------------
  void FragmentPlots::clear()
  {
  }

  DEFINE_ART_RESULTS_PLUGIN(FragmentPlots)
}// end tcpCompression namespace
