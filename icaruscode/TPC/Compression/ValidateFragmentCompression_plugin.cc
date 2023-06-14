//std includes
#include <vector>

//ROOT includes
#include "TFile.h"

//Framework includes
#include "art/Framework/Core/ResultsProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Results.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
//#include "art/Framework/Services/Optional/TFileService.h"
//#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//icarus includes
#include "icaruscode/TPC/Compression/PhysCrateCompressedFragment.cc"

namespace tcpCompression {
  class ValidateCompression : public art::ResultsProducer {
  public:
    explicit ValidateCompression(fhicl::ParameterSet const& pset);
    ~ValidateCompression() override = default;
    void event(art::event& evt) override; 
    void writeResults(art::Results& r) override;
    void reconfigure(fhicl::ParameterSet const& p);
    void FillDataSpectrum(cmf::Spectrum & spec);
    void clear()    override;
    void endJob()   override;

  private:
    art::InputTag fFragmentsLabel;
  };// end ValidateCompression class

  //------------------------------------------------------------------
  ValidateCompression::ValidateCompression(fhicl::ParameterSet const& pset)
  {
    this->reconfigure(pset);
  }

  //------------------------------------------------------------------
  void ValidateCompression::reconfigure(fhicl::ParameterSet const& pset)
  {
    fFragmentsLabel = pset.get<art::InputTag>("FragmentsLabel", "daq:PHYSCRATEDATA");
  }

  //------------------------------------------------------------------
  void ValidateCompression::event(art::event& evt)
  {
    // get the uncompressed fragments from the event
    art::Handle<std::vector<artdaq::Fragment>> originalFragmentsHandle;
    std::vector<art::Ptr<artdaq::Fragment>>    originalFragments;

    for (auto const& frag : originalFragments)
    {
      // put the fragment in an overlay
      icarus::PhysCrateCompressedFragment fragOverlay(*frag);

      // print the compression scheme
      MF_LOG_VERBATIM("ValidateCompression")
        << "Compression scheme is " << fragOverlay.CompressionScheme() << ", data size is " << fragOverlay.DataPayloadSize();
    }
  }

  DEFINE_ART_RESULTS_PLUGIN(ValidateCompression)
}// end tcpCompression namespace
