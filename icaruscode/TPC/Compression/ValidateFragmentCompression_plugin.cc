//std includes
#include <vector>

//ROOT includes
#include "TFile.h"

//Framework includes
#include "art/Framework/Core/ResultsProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Results.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
//#include "art/Framework/Services/Optional/TFileService.h"
//#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Utilities/InputTag.h"

//icarus includes
#include "icaruscode/TPC/Compression/PhysCrateCompressedFragment.cc"

//sbndaq includes
#include "sbndaq-artdaq-core/Overlays/ICARUS/PhysCrateFragment.hh"

namespace tcpCompression {
  class ValidateCompression : public art::ResultsProducer {
  public:
    explicit ValidateCompression(fhicl::ParameterSet const& pset);
    ~ValidateCompression() override = default;
    void event(art::Event const& evt) override; 
    void reconfigure(fhicl::ParameterSet const& p);
    void writeResults(art::Results& r) override;
    void clear() override;

  private:
    art::InputTag fFragmentsLabel;
    bool          fCheckOldFragments;
    bool          fFoundFirstUncompressed = false;
    bool          fFoundFirstCompressed   = false;
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
    fCheckOldFragments = pset.get<bool>("CheckOldFragments", false);
  }

  //------------------------------------------------------------------
  void ValidateCompression::event(art::Event const& evt)
  {
    // get the uncompressed fragments from the event
    //art::Handle<std::vector<artdaq::Fragment>> originalFragmentsHandle;
    //std::vector<art::Ptr<artdaq::Fragment>>    originalFragments;
    //if (e.getByLabel(fFragmentsLabel, originalFragmentsHandle))
    //  art::fill_ptr_vector(originalFragments, originalFragmentsHandle);
    std::vector<artdaq::Fragment> originalFragments = evt.getProduct<std::vector<artdaq::Fragment>>(fFragmentsLabel);

    for (auto const& frag : originalFragments)
    {
      // test that the old overlays work...
      if (fCheckOldFragments)
      {
        icarus::PhysCrateFragment oldOverlay(frag);
        MF_LOG_VERBATIM("ValidateCompression")
          << "Verifying old overlays..." << '\n'
          << " * This fragment has " << oldOverlay.nBoards() << " boards," << '\n'
          << " *                   " << oldOverlay.nChannelsPerBoard() << " channels per board," << '\n'
          << " *               and " << oldOverlay.nSamplesPerChannel() << " samples per channel" << '\n'
          << "-----------------------------------------------";
        oldOverlay.Verify();
        MF_LOG_VERBATIM("ValidateCompression")
          << "-----------------------------------------------";
        continue;
      }

      // put the fragment in a compression compliant overlay
      icarus::PhysCrateCompressedFragment fragOverlay(frag);
      icarus::PhysCrateFragment fragOldOverlay(frag);
      bool isComp = fragOverlay.isCompressed();
      if (not isComp)
      {
        MF_LOG_VERBATIM("ValidateCompression")
          << "Checking first few ADC values for uncompressed fragment" << '\n'
          << " ADC board 0, channel 0, sample 0: " << fragOverlay.adc_val(0, 0, 0) << " | " << fragOldOverlay.adc_val(0, 0, 0) << '\n'
          << " ADC board 0, channel 0, sample 1: " << fragOverlay.adc_val(0, 0, 1) << " | " << fragOldOverlay.adc_val(0, 0, 1) << '\n'
          << " ADC board 0, channel 0, sample 2: " << fragOverlay.adc_val(0, 0, 2) << " | " << fragOldOverlay.adc_val(0, 0, 2) << '\n'
          << " ADC board 0, channel 1, sample 0: " << fragOverlay.adc_val(0, 1, 0) << " | " << fragOldOverlay.adc_val(0, 1, 0) << '\n'
          << " ADC board 0, channel 1, sample 1: " << fragOverlay.adc_val(0, 1, 1) << " | " << fragOldOverlay.adc_val(0, 1, 1) << '\n'
          << " ADC board 0, channel 1, sample 2: " << fragOverlay.adc_val(0, 1, 2) << " | " << fragOldOverlay.adc_val(0, 1, 2) << '\n'
          << " ADC board 0, channel 2, sample 0: " << fragOverlay.adc_val(0, 2, 0) << " | " << fragOldOverlay.adc_val(0, 2, 0) << '\n'
          << " ADC board 0, channel 2, sample 1: " << fragOverlay.adc_val(0, 2, 1) << " | " << fragOldOverlay.adc_val(0, 2, 1) << '\n'
          << " ADC board 0, channel 2, sample 2: " << fragOverlay.adc_val(0, 2, 2) << " | " << fragOldOverlay.adc_val(0, 2, 2) << '\n'
          << " ADC board 1, channel 0, sample 0: " << fragOverlay.adc_val(1, 0, 0) << " | " << fragOldOverlay.adc_val(0, 0, 0) << '\n'
          << " ADC board 1, channel 0, sample 1: " << fragOverlay.adc_val(1, 0, 1) << " | " << fragOldOverlay.adc_val(1, 0, 1) << '\n'
          << " ADC board 1, channel 0, sample 2: " << fragOverlay.adc_val(1, 0, 2) << " | " << fragOldOverlay.adc_val(1, 0, 2) << '\n'
          << " ADC board 1, channel 1, sample 0: " << fragOverlay.adc_val(1, 1, 0) << " | " << fragOldOverlay.adc_val(1, 1, 0) << '\n'
          << " ADC board 1, channel 1, sample 1: " << fragOverlay.adc_val(1, 1, 1) << " | " << fragOldOverlay.adc_val(1, 1, 1) << '\n'
          << " ADC board 1, channel 1, sample 2: " << fragOverlay.adc_val(1, 1, 2) << " | " << fragOldOverlay.adc_val(1, 1, 2) << '\n'
          << " ADC board 1, channel 2, sample 0: " << fragOverlay.adc_val(1, 2, 0) << " | " << fragOldOverlay.adc_val(1, 2, 0) << '\n'
          << " ADC board 1, channel 2, sample 1: " << fragOverlay.adc_val(1, 2, 1) << " | " << fragOldOverlay.adc_val(1, 2, 1) << '\n'
          << " ADC board 1, channel 2, sample 2: " << fragOverlay.adc_val(1, 2, 2) << " | " << fragOldOverlay.adc_val(1, 2, 2) << '\n'
          << " ADC board 2, channel 0, sample 0: " << fragOverlay.adc_val(2, 0, 0) << " | " << fragOldOverlay.adc_val(2, 0, 0) << '\n'
          << " ADC board 2, channel 0, sample 1: " << fragOverlay.adc_val(2, 0, 1) << " | " << fragOldOverlay.adc_val(2, 0, 1) << '\n'
          << " ADC board 2, channel 0, sample 2: " << fragOverlay.adc_val(2, 0, 2) << " | " << fragOldOverlay.adc_val(2, 0, 2) << '\n'
          << " ADC board 2, channel 1, sample 0: " << fragOverlay.adc_val(2, 1, 0) << " | " << fragOldOverlay.adc_val(2, 1, 0) << '\n'
          << " ADC board 2, channel 1, sample 1: " << fragOverlay.adc_val(2, 1, 1) << " | " << fragOldOverlay.adc_val(2, 1, 1) << '\n'
          << " ADC board 2, channel 1, sample 2: " << fragOverlay.adc_val(2, 1, 2) << " | " << fragOldOverlay.adc_val(2, 1, 2) << '\n'
          << " ADC board 2, channel 2, sample 0: " << fragOverlay.adc_val(2, 2, 0) << " | " << fragOldOverlay.adc_val(2, 2, 0) << '\n'
          << " ADC board 2, channel 2, sample 1: " << fragOverlay.adc_val(2, 2, 1) << " | " << fragOldOverlay.adc_val(2, 2, 1) << '\n'
          << " ADC board 2, channel 2, sample 2: " << fragOverlay.adc_val(2, 2, 2) << " | " << fragOldOverlay.adc_val(2, 2, 2) << '\n';
          continue;
      }
      if (isComp)
      {
        MF_LOG_VERBATIM("ValidateCompression")
          << "Checking first few ADC values for compressed fragment" << '\n'
          << " ADC board 0, channel 0, sample 0: " << fragOverlay.adc_val(0, 0, 0) << " | " << fragOldOverlay.adc_val(0, 0, 0) << '\n'
          << " ADC board 0, channel 0, sample 1: " << fragOverlay.adc_val(0, 0, 1) << " | " << fragOldOverlay.adc_val(0, 0, 1) << '\n'
          << " ADC board 0, channel 0, sample 2: " << fragOverlay.adc_val(0, 0, 2) << " | " << fragOldOverlay.adc_val(0, 0, 2) << '\n'
          << " ADC board 0, channel 1, sample 0: " << fragOverlay.adc_val(0, 1, 0) << " | " << fragOldOverlay.adc_val(0, 1, 0) << '\n'
          << " ADC board 0, channel 1, sample 1: " << fragOverlay.adc_val(0, 1, 1) << " | " << fragOldOverlay.adc_val(0, 1, 1) << '\n'
          << " ADC board 0, channel 1, sample 2: " << fragOverlay.adc_val(0, 1, 2) << " | " << fragOldOverlay.adc_val(0, 1, 2) << '\n'
          << " ADC board 0, channel 2, sample 0: " << fragOverlay.adc_val(0, 2, 0) << " | " << fragOldOverlay.adc_val(0, 2, 0) << '\n'
          << " ADC board 0, channel 2, sample 1: " << fragOverlay.adc_val(0, 2, 1) << " | " << fragOldOverlay.adc_val(0, 2, 1) << '\n'
          << " ADC board 0, channel 2, sample 2: " << fragOverlay.adc_val(0, 2, 2) << " | " << fragOldOverlay.adc_val(0, 2, 2) << '\n'
          << " ADC board 1, channel 0, sample 0: " << fragOverlay.adc_val(1, 0, 0) << " | " << fragOldOverlay.adc_val(1, 0, 0) << '\n'
          << " ADC board 1, channel 0, sample 1: " << fragOverlay.adc_val(1, 0, 1) << " | " << fragOldOverlay.adc_val(1, 0, 1) << '\n'
          << " ADC board 1, channel 0, sample 2: " << fragOverlay.adc_val(1, 0, 2) << " | " << fragOldOverlay.adc_val(1, 0, 2) << '\n'
          << " ADC board 1, channel 1, sample 0: " << fragOverlay.adc_val(1, 1, 0) << " | " << fragOldOverlay.adc_val(1, 1, 0) << '\n'
          << " ADC board 1, channel 1, sample 1: " << fragOverlay.adc_val(1, 1, 1) << " | " << fragOldOverlay.adc_val(1, 1, 1) << '\n'
          << " ADC board 1, channel 1, sample 2: " << fragOverlay.adc_val(1, 1, 2) << " | " << fragOldOverlay.adc_val(1, 1, 2) << '\n'
          << " ADC board 1, channel 2, sample 0: " << fragOverlay.adc_val(1, 2, 0) << " | " << fragOldOverlay.adc_val(1, 2, 0) << '\n'
          << " ADC board 1, channel 2, sample 1: " << fragOverlay.adc_val(1, 2, 1) << " | " << fragOldOverlay.adc_val(1, 2, 1) << '\n'
          << " ADC board 1, channel 2, sample 2: " << fragOverlay.adc_val(1, 2, 2) << " | " << fragOldOverlay.adc_val(1, 2, 2) << '\n'
          << " ADC board 2, channel 0, sample 0: " << fragOverlay.adc_val(2, 0, 0) << " | " << fragOldOverlay.adc_val(2, 0, 0) << '\n'
          << " ADC board 2, channel 0, sample 1: " << fragOverlay.adc_val(2, 0, 1) << " | " << fragOldOverlay.adc_val(2, 0, 1) << '\n'
          << " ADC board 2, channel 0, sample 2: " << fragOverlay.adc_val(2, 0, 2) << " | " << fragOldOverlay.adc_val(2, 0, 2) << '\n'
          << " ADC board 2, channel 1, sample 0: " << fragOverlay.adc_val(2, 1, 0) << " | " << fragOldOverlay.adc_val(2, 1, 0) << '\n'
          << " ADC board 2, channel 1, sample 1: " << fragOverlay.adc_val(2, 1, 1) << " | " << fragOldOverlay.adc_val(2, 1, 1) << '\n'
          << " ADC board 2, channel 1, sample 2: " << fragOverlay.adc_val(2, 1, 2) << " | " << fragOldOverlay.adc_val(2, 1, 2) << '\n'
          << " ADC board 2, channel 2, sample 0: " << fragOverlay.adc_val(2, 2, 0) << " | " << fragOldOverlay.adc_val(2, 2, 0) << '\n'
          << " ADC board 2, channel 2, sample 1: " << fragOverlay.adc_val(2, 2, 1) << " | " << fragOldOverlay.adc_val(2, 2, 1) << '\n'
          << " ADC board 2, channel 2, sample 2: " << fragOverlay.adc_val(2, 2, 2) << " | " << fragOldOverlay.adc_val(2, 2, 2) << '\n';
          continue;
      }
      

      //bool originalValid = fragOverlay.Verify();
      //icarus::PhysCrateCompressedFragment newFragOverlay = (isComp) ? fragOverlay.makeUncompressedFragment() :
      //                                                                fragOverlay.  makeCompressedFragment() ;
      //bool newValid = newFragOverlay.Verify();

      //std::string oldValidStr = (originalValid) ? "Original fragment is valid" : "Original fragment is not valid";
      //std::string newValidStr = (newValid)      ? "new fragment is valid"      : "new fragment is not valid";
      //if (originalValid == newValid)
      //{
      //  MF_LOG_VERBATIM("ValidateCompression")
      //    << oldValidStr << ", and " << newValidStr;
      //} else {
      //  MF_LOG_VERBATIM("ValidateCompression")
      //    << oldValidStr << ", but " << newValidStr;
      //}

      //if (not originalValid)
      //{
      //  MF_LOG_VERBATIM("ValidateCompression")
      //    << "Attempting to make valid fragment...";
      //  icarus::PhysCrateCompressedFragment reFragmentOverlay = (isComp) ? newFragOverlay.  makeCompressedFragment() :
      //                                                                     newFragOverlay.makeUncompressedFragment() ;
      //  if (isComp)
      //  {
      //    MF_LOG_VERBATIM("ValidateCompression")
      //      << "Validating recompressed fragment...";
      //    reFragmentOverlay.Verify();
      //  } else {
      //    MF_LOG_VERBATIM("ValidateCompression")
      //      << "Validating re-decompressed fragment...";
      //    reFragmentOverlay.Verify();
      //  }
      //}

    }
  }

  //------------------------------------------------------------------
  void ValidateCompression::writeResults(art::Results& r)
  {
    MF_LOG_VERBATIM("ValidateCompression")
      << "Not writing results (yet)";
  }

  //------------------------------------------------------------------
  void ValidateCompression::clear()
  {
  }

  DEFINE_ART_RESULTS_PLUGIN(ValidateCompression)
}// end tcpCompression namespace
