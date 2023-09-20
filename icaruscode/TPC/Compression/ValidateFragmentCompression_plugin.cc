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
    bool          fDumpADCs;
    bool          fFoundFirstUncompressed = false;
    bool          fFoundFirstCompressed   = false;
    TH1F*         fCompHist;
    TH1F*         fDcmpHist;
    const icarusDB::IICARUSChannelMap* fChannelMap;
  };// end ValidateCompression class

  //------------------------------------------------------------------
  ValidateCompression::ValidateCompression(fhicl::ParameterSet const& pset)
  {
    this->reconfigure(pset);
  }

  //------------------------------------------------------------------
  void ValidateCompression::reconfigure(fhicl::ParameterSet const& pset)
  {
    mf::LogVerbatim("DaqDecoderICARUSTPCwROI") << "NO HARRY YOU ARE NOT CRAZY THIS IS THE RIGHT FILE (DaqDecoderICARUSTPCwROI) TO EDIT";

    fFragmentsLabel    = pset.get<art::InputTag>("FragmentsLabel"   , "daq:PHYSCRATEDATA");
    fCheckOldFragments = pset.get<bool>         ("CheckOldFragments", false);
    fDumpADCs          = pset.get<bool>         ("DumpADCs"         , false);

    art::ServiceHandle<art::TFileService> tfs;
    fChannelMap = art::ServiceHandle<icarusDB::IICARUSChannelMap const>{}.get();

    fCompHist = tfs->make<TH1F>(("CompHist_"+fFragmentsLabel.instance()).c_str(), ";Size (bytes);Number of Fragments", 100, 0, 10000000);
    fDcmpHist = tfs->make<TH1F>(("DcmpHist_"+fFragmentsLabel.instance()).c_str(), ";Size (bytes);Number of Fragments", 100, 0, 10000000);
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

    // if we want to make a histogram for the waveform
    art::ServiceHandle<art::TFileService> tfs;

    // this is useful
    size_t fragNumber = 0;
    for (auto const& frag : originalFragments)
    {
      // test that the old overlays work...
      //if (fCheckOldFragments)
      //{
      //  icarus::PhysCrateFragment oldOverlay(frag);
      //  MF_LOG_VERBATIM("ValidateCompression")
      //    << "Verifying old overlays..." << '\n'
      //    << " * This fragment has " << oldOverlay.nBoards() << " boards," << '\n'
      //    << " *                   " << oldOverlay.nChannelsPerBoard() << " channels per board," << '\n'
      //    << " *               and " << oldOverlay.nSamplesPerChannel() << " samples per channel" << '\n'
      //    << "-----------------------------------------------";
      //  oldOverlay.Verify();
      //  MF_LOG_VERBATIM("ValidateCompression")
      //    << "-----------------------------------------------";
      //  continue;
      //}

      // put the fragment in a compression compliant overlay
      icarus::PhysCrateCompressedFragment fragOverlay(frag);
      icarus::PhysCrateFragment fragOldOverlay(frag);
      bool isComp = fragOverlay.isCompressed();

      if (isComp)
      {
        std::string compStr = (fragOldOverlay.isCompressed()) ? "Yes!" : "no :(";
        MF_LOG_VERBATIM("ValidateCompression")
          << "THe icaruscode overlay is compressed." << '\n'
          << "... Is the sbncode overlay compressed? " << compStr;
      }

      //if (fDumpADCs && not isComp)
      //{
      //  MF_LOG_VERBATIM("ValidateCompression")
      //    << "Checking first few ADC values for uncompressed fragment" << '\n'
      //    << "-----------------------NEW HEADER-----------------------" << '\n'
      //    << *(fragOverlay.DataTileHeader(0)) << '\n'
      //    << *(fragOverlay.DataTileHeader(1)) << '\n'
      //    << *(fragOverlay.DataTileHeader(2)) << '\n'
      //    << "-----------------------OLD HEADER-----------------------" << '\n'
      //    << *(fragOldOverlay.DataTileHeader(0)) << '\n'
      //    << *(fragOldOverlay.DataTileHeader(1)) << '\n'
      //    << *(fragOldOverlay.DataTileHeader(2)) << '\n'
      //    << "--------------------------------------------------------" << '\n'
      //    << "Crate ID (new): " << std::hex << fragOverlay.DataTileHeader()->crate_id << std::dec << '\n'
      //    << "Crate ID (old): " << std::hex << fragOldOverlay.DataTileHeader()->crate_id << std::dec << '\n'
      //    << " ADC board 0, channel 0, sample 0: " << fragOverlay.adc_val(0, 0, 0) << " | " << fragOldOverlay.adc_val(0, 0, 0) << '\n'
      //    << " ADC board 0, channel 0, sample 1: " << fragOverlay.adc_val(0, 0, 1) << " | " << fragOldOverlay.adc_val(0, 0, 1) << '\n'
      //    << " ADC board 0, channel 0, sample 2: " << fragOverlay.adc_val(0, 0, 2) << " | " << fragOldOverlay.adc_val(0, 0, 2) << '\n'
      //    << " ADC board 0, channel 1, sample 0: " << fragOverlay.adc_val(0, 1, 0) << " | " << fragOldOverlay.adc_val(0, 1, 0) << '\n'
      //    << " ADC board 0, channel 1, sample 1: " << fragOverlay.adc_val(0, 1, 1) << " | " << fragOldOverlay.adc_val(0, 1, 1) << '\n'
      //    << " ADC board 0, channel 1, sample 2: " << fragOverlay.adc_val(0, 1, 2) << " | " << fragOldOverlay.adc_val(0, 1, 2) << '\n'
      //    << " ADC board 0, channel 2, sample 0: " << fragOverlay.adc_val(0, 2, 0) << " | " << fragOldOverlay.adc_val(0, 2, 0) << '\n'
      //    << " ADC board 0, channel 2, sample 1: " << fragOverlay.adc_val(0, 2, 1) << " | " << fragOldOverlay.adc_val(0, 2, 1) << '\n'
      //    << " ADC board 0, channel 2, sample 2: " << fragOverlay.adc_val(0, 2, 2) << " | " << fragOldOverlay.adc_val(0, 2, 2) << '\n'
      //    << " ADC board 1, channel 0, sample 0: " << fragOverlay.adc_val(1, 0, 0) << " | " << fragOldOverlay.adc_val(0, 0, 0) << '\n'
      //    << " ADC board 1, channel 0, sample 1: " << fragOverlay.adc_val(1, 0, 1) << " | " << fragOldOverlay.adc_val(1, 0, 1) << '\n'
      //    << " ADC board 1, channel 0, sample 2: " << fragOverlay.adc_val(1, 0, 2) << " | " << fragOldOverlay.adc_val(1, 0, 2) << '\n'
      //    << " ADC board 1, channel 1, sample 0: " << fragOverlay.adc_val(1, 1, 0) << " | " << fragOldOverlay.adc_val(1, 1, 0) << '\n'
      //    << " ADC board 1, channel 1, sample 1: " << fragOverlay.adc_val(1, 1, 1) << " | " << fragOldOverlay.adc_val(1, 1, 1) << '\n'
      //    << " ADC board 1, channel 1, sample 2: " << fragOverlay.adc_val(1, 1, 2) << " | " << fragOldOverlay.adc_val(1, 1, 2) << '\n'
      //    << " ADC board 1, channel 2, sample 0: " << fragOverlay.adc_val(1, 2, 0) << " | " << fragOldOverlay.adc_val(1, 2, 0) << '\n'
      //    << " ADC board 1, channel 2, sample 1: " << fragOverlay.adc_val(1, 2, 1) << " | " << fragOldOverlay.adc_val(1, 2, 1) << '\n'
      //    << " ADC board 1, channel 2, sample 2: " << fragOverlay.adc_val(1, 2, 2) << " | " << fragOldOverlay.adc_val(1, 2, 2) << '\n'
      //    << " ADC board 2, channel 0, sample 0: " << fragOverlay.adc_val(2, 0, 0) << " | " << fragOldOverlay.adc_val(2, 0, 0) << '\n'
      //    << " ADC board 2, channel 0, sample 1: " << fragOverlay.adc_val(2, 0, 1) << " | " << fragOldOverlay.adc_val(2, 0, 1) << '\n'
      //    << " ADC board 2, channel 0, sample 2: " << fragOverlay.adc_val(2, 0, 2) << " | " << fragOldOverlay.adc_val(2, 0, 2) << '\n'
      //    << " ADC board 2, channel 1, sample 0: " << fragOverlay.adc_val(2, 1, 0) << " | " << fragOldOverlay.adc_val(2, 1, 0) << '\n'
      //    << " ADC board 2, channel 1, sample 1: " << fragOverlay.adc_val(2, 1, 1) << " | " << fragOldOverlay.adc_val(2, 1, 1) << '\n'
      //    << " ADC board 2, channel 1, sample 2: " << fragOverlay.adc_val(2, 1, 2) << " | " << fragOldOverlay.adc_val(2, 1, 2) << '\n'
      //    << " ADC board 2, channel 2, sample 0: " << fragOverlay.adc_val(2, 2, 0) << " | " << fragOldOverlay.adc_val(2, 2, 0) << '\n'
      //    << " ADC board 2, channel 2, sample 1: " << fragOverlay.adc_val(2, 2, 1) << " | " << fragOldOverlay.adc_val(2, 2, 1) << '\n'
      //    << " ADC board 2, channel 2, sample 2: " << fragOverlay.adc_val(2, 2, 2) << " | " << fragOldOverlay.adc_val(2, 2, 2) << '\n';
      //    continue;
      //}
      //if (fDumpADCs && isComp && (fragOverlay.CompressionKey(0,0) == 0))
      //{
      //  MF_LOG_VERBATIM("ValidateCompression")
      //    << "Checking first few ADC values for compressed fragment" << '\n'
      //    << "-----------------------NEW HEADER-----------------------" << '\n'
      //    << *(fragOverlay.DataTileHeader(0)) << '\n'
      //    << *(fragOverlay.DataTileHeader(1)) << '\n'
      //    << *(fragOverlay.DataTileHeader(2)) << '\n'
      //    << "-----------------------OLD HEADER-----------------------" << '\n'
      //    << *(fragOldOverlay.DataTileHeader(0)) << '\n'
      //    << *(fragOldOverlay.DataTileHeader(1)) << '\n'
      //    << *(fragOldOverlay.DataTileHeader(2)) << '\n'
      //    << "--------------------------------------------------------" << '\n'
      //    << "Crate ID (new): " << std::hex << fragOverlay.DataTileHeader()->crate_id << std::dec << '\n'
      //    << "Crate ID (old): " << std::hex << fragOldOverlay.DataTileHeader()->crate_id << std::dec << '\n'
      //    << " ADC board 0, channel 0, sample 0: " << fragOverlay.adc_val(0, 0, 0) << " | " << fragOldOverlay.adc_val(0, 0, 0) << '\n'
      //    << " ADC board 0, channel 0, sample 1: " << fragOverlay.adc_val(0, 0, 1) << " | " << fragOldOverlay.adc_val(0, 0, 1) << '\n'
      //    << " ADC board 0, channel 0, sample 2: " << fragOverlay.adc_val(0, 0, 2) << " | " << fragOldOverlay.adc_val(0, 0, 2) << '\n'
      //    << " ADC board 0, channel 1, sample 0: " << fragOverlay.adc_val(0, 1, 0) << " | " << fragOldOverlay.adc_val(0, 1, 0) << '\n'
      //    << " ADC board 0, channel 1, sample 1: " << fragOverlay.adc_val(0, 1, 1) << " | " << fragOldOverlay.adc_val(0, 1, 1) << '\n'
      //    << " ADC board 0, channel 1, sample 2: " << fragOverlay.adc_val(0, 1, 2) << " | " << fragOldOverlay.adc_val(0, 1, 2) << '\n'
      //    << " ADC board 0, channel 2, sample 0: " << fragOverlay.adc_val(0, 2, 0) << " | " << fragOldOverlay.adc_val(0, 2, 0) << '\n'
      //    << " ADC board 0, channel 2, sample 1: " << fragOverlay.adc_val(0, 2, 1) << " | " << fragOldOverlay.adc_val(0, 2, 1) << '\n'
      //    << " ADC board 0, channel 2, sample 2: " << fragOverlay.adc_val(0, 2, 2) << " | " << fragOldOverlay.adc_val(0, 2, 2) << '\n'
      //    << " ADC board 1, channel 0, sample 0: " << fragOverlay.adc_val(1, 0, 0) << " | " << fragOldOverlay.adc_val(1, 0, 0) << '\n'
      //    << " ADC board 1, channel 0, sample 1: " << fragOverlay.adc_val(1, 0, 1) << " | " << fragOldOverlay.adc_val(1, 0, 1) << '\n'
      //    << " ADC board 1, channel 0, sample 2: " << fragOverlay.adc_val(1, 0, 2) << " | " << fragOldOverlay.adc_val(1, 0, 2) << '\n'
      //    << " ADC board 1, channel 1, sample 0: " << fragOverlay.adc_val(1, 1, 0) << " | " << fragOldOverlay.adc_val(1, 1, 0) << '\n'
      //    << " ADC board 1, channel 1, sample 1: " << fragOverlay.adc_val(1, 1, 1) << " | " << fragOldOverlay.adc_val(1, 1, 1) << '\n'
      //    << " ADC board 1, channel 1, sample 2: " << fragOverlay.adc_val(1, 1, 2) << " | " << fragOldOverlay.adc_val(1, 1, 2) << '\n'
      //    << " ADC board 1, channel 2, sample 0: " << fragOverlay.adc_val(1, 2, 0) << " | " << fragOldOverlay.adc_val(1, 2, 0) << '\n'
      //    << " ADC board 1, channel 2, sample 1: " << fragOverlay.adc_val(1, 2, 1) << " | " << fragOldOverlay.adc_val(1, 2, 1) << '\n'
      //    << " ADC board 1, channel 2, sample 2: " << fragOverlay.adc_val(1, 2, 2) << " | " << fragOldOverlay.adc_val(1, 2, 2) << '\n'
      //    << " ADC board 2, channel 0, sample 0: " << fragOverlay.adc_val(2, 0, 0) << " | " << fragOldOverlay.adc_val(2, 0, 0) << '\n'
      //    << " ADC board 2, channel 0, sample 1: " << fragOverlay.adc_val(2, 0, 1) << " | " << fragOldOverlay.adc_val(2, 0, 1) << '\n'
      //    << " ADC board 2, channel 0, sample 2: " << fragOverlay.adc_val(2, 0, 2) << " | " << fragOldOverlay.adc_val(2, 0, 2) << '\n'
      //    << " ADC board 2, channel 1, sample 0: " << fragOverlay.adc_val(2, 1, 0) << " | " << fragOldOverlay.adc_val(2, 1, 0) << '\n'
      //    << " ADC board 2, channel 1, sample 1: " << fragOverlay.adc_val(2, 1, 1) << " | " << fragOldOverlay.adc_val(2, 1, 1) << '\n'
      //    << " ADC board 2, channel 1, sample 2: " << fragOverlay.adc_val(2, 1, 2) << " | " << fragOldOverlay.adc_val(2, 1, 2) << '\n'
      //    << " ADC board 2, channel 2, sample 0: " << fragOverlay.adc_val(2, 2, 0) << " | " << fragOldOverlay.adc_val(2, 2, 0) << '\n'
      //    << " ADC board 2, channel 2, sample 1: " << fragOverlay.adc_val(2, 2, 1) << " | " << fragOldOverlay.adc_val(2, 2, 1) << '\n'
      //    << " ADC board 2, channel 2, sample 2: " << fragOverlay.adc_val(2, 2, 2) << " | " << fragOldOverlay.adc_val(2, 2, 2) << '\n';
      //    for (size_t b = 0; b < fragOverlay.nBoards(); ++b)
      //    {
      //      for (size_t s = 0; s < 10; ++s)
      //      {
      //        MF_LOG_VERBATIM("ValidateCompression")
      //          << "  Board " << b << " Sample " << s << " Key: " << std::bitset<16>(fragOverlay.CompressionKey(b, s));
      //      }
      //    }
      //    continue;
      //}

      // here's where we fill the plots
      bool isActuallyCompresed = (fragOverlay.CompressionKey(0,0) == 0);
      std::string fragCrateName = fChannelMap->getCrateName(frag.fragmentID());
      if (isComp && isActuallyCompresed)
      {
        MF_LOG_VERBATIM("ValidateCompression")
          << "Filling Compressed Hist with " << frag.sizeBytes() << '\n'
          << "|-> Fragment comes from crate " << fragCrateName;
        fCompHist->Fill(frag.sizeBytes());
        
        // look through the fragment and try to find peaks
        // for each channel loop over the keys and look for where there's no compression
        // ignore sample zero because that's always uncompressed
        for (size_t b = 0; b < fragOverlay.nBoards(); ++b)
        {
          for (size_t c = 0; c < fragOverlay.nChannelsPerBoard(); ++c)
          {
            std::string title = "ADC_"+fFragmentsLabel.label()+":"+fFragmentsLabel.instance()+"_Event_"+std::to_string(evt.id().event())+"_Fragment_"+std::to_string(fragNumber)+"_Board_"+std::to_string(b)+"_Channel_"+std::to_string(c);
            TH1F* hist = tfs->make<TH1F>(title.c_str(), ";Sample;ADC Counts", fragOverlay.nSamplesPerChannel(), 0, fragOverlay.nSamplesPerChannel());
            for (size_t s = 0; s < fragOverlay.nSamplesPerChannel(); ++s)
            {
              hist->SetBinContent(s, fragOverlay.adc_val(b, c, s));
              if (fragOverlay.adc_val(b, c, s) > 4000)
              {
                MF_LOG_VERBATIM("ValidateCompression")
                  << "ADV Value for Board " << b << " Channel " << c << " Sample " << s << " is " << fragOverlay.adc_val(b, c, s) << '\n'
                  << "  Cmpression key for Board/Sample is " << std::bitset<16>(fragOverlay.CompressionKey(b, s));
              }
            }
          }
        }
      } else {
        MF_LOG_VERBATIM("ValidateCompression")
          << "Filling Decompressed Hist with " << frag.sizeBytes() << '\n'
          << "|-> Fragment comes from crate " << fragCrateName;
        fDcmpHist->Fill(frag.sizeBytes());
        for (size_t b = 0; b < fragOldOverlay.nBoards(); ++b)
        {
          for (size_t c = 0; c < fragOldOverlay.nChannelsPerBoard(); ++c)
          {
            std::string title = "ADC_"+fFragmentsLabel.label()+":"+fFragmentsLabel.instance()+"_Event_"+std::to_string(evt.id().event())+"_Fragment_"+std::to_string(fragNumber)+"_Board_"+std::to_string(b)+"_Channel_"+std::to_string(c);
            TH1F* hist = tfs->make<TH1F>(title.c_str(), ";Sample;ADC Counts", fragOldOverlay.nSamplesPerChannel(), 0, fragOldOverlay.nSamplesPerChannel());
            for (size_t s = 0; s < fragOldOverlay.nSamplesPerChannel(); ++s)
            {
              hist->SetBinContent(s, fragOldOverlay.adc_val(b, c, s));
            }
          }
        }
      }

      ++fragNumber;
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
