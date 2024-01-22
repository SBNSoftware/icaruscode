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
#include "icaruscode/Decode/ChannelMapping/IICARUSChannelMap.h"

//sbndaq includes
#include "sbndaq-artdaq-core/Overlays/ICARUS/PhysCrateFragment.h"

//lardataobj inlcudes
#include "lardataobj/RawData/RawDigit.h"

namespace tcpCompression {
  class CheckRawDigits : public art::ResultsProducer {
  public:
    explicit CheckRawDigits(fhicl::ParameterSet const& pset);
    ~CheckRawDigits() override = default;
    void event(art::Event const& evt) override; 
    void reconfigure(fhicl::ParameterSet const& p);
    void writeResults(art::Results& r) override;
    void clear() override;

  private:
    std::vector<art::InputTag>         fRawDigitLabelVec;
    std::vector<TH1F*>                 fADCHistVec;
    const icarusDB::IICARUSChannelMap* fChannelMap;
    bool                               fDoFragments;
  };// end CheckRawDigits class

  //------------------------------------------------------------------
  CheckRawDigits::CheckRawDigits(fhicl::ParameterSet const& pset)
  {
    this->reconfigure(pset);
  }

  //------------------------------------------------------------------
  void CheckRawDigits::reconfigure(fhicl::ParameterSet const& pset)
  {
    fRawDigitLabelVec = pset.get<std::vector<art::InputTag>>("RawDigitLabelVec");
    fDoFragments      = pset.get<bool>("DoFragments", false); 

    art::ServiceHandle<art::TFileService> tfs;
    fChannelMap = art::ServiceHandle<icarusDB::IICARUSChannelMap const>{}.get();

    for (size_t idx = 0; idx < fRawDigitLabelVec.size(); ++idx)
    {
      if (not fDoFragments)
      {
        fADCHistVec.emplace_back(tfs->make<TH1F>(("ADC_"+fRawDigitLabelVec[idx].instance()).c_str(), (fRawDigitLabelVec[idx].instance()+";ADC;Samples").c_str(), 500, -300, 200));
      } else {
        fADCHistVec.emplace_back(tfs->make<TH1F>(("ADC_"+fRawDigitLabelVec[idx].instance()).c_str(), (fRawDigitLabelVec[idx].instance()+";ADC;Samples").c_str(), 4000, 1000, 5000));
      }
    }
  }

  //------------------------------------------------------------------
  void CheckRawDigits::event(art::Event const& evt)
  {
    for (size_t idx = 0; idx < fRawDigitLabelVec.size(); ++idx)
    {
      if (not fDoFragments)
      {
        // get the RawDigits
        std::vector<raw::RawDigit> rawDigits = evt.getProduct<std::vector<raw::RawDigit>>(fRawDigitLabelVec[idx]);
 
        // get the Histogram
        TH1F* hist = fADCHistVec[idx];

        // loop over the raw digits and fill the histogram
        for (auto const& rd : rawDigits)
        {
          // loop over the ADCs in the RawDigit
          for (size_t sample = 0; sample < rd.Samples(); ++sample)
          {
            hist->Fill(rd.ADC(sample));
          }
        }
      } else {
        // get the fragments
        std::vector<artdaq::Fragment> fragments = evt.getProduct<std::vector<artdaq::Fragment>>(fRawDigitLabelVec[idx]);

        // get the histogram
        TH1F* hist = fADCHistVec[idx];

        // loop over the raw digits and fill the histogram
        for (auto const& frag : fragments)
        {
          // put on an overlay
          icarus::PhysCrateFragment overlay(frag);
          MF_LOG_VERBATIM("CheckRawDigits")
            << "Is Fragment Compressed? " << overlay.CompressionScheme();

          // for each board/channel/sample get the ADC value
          for (size_t board = 0; board < overlay.nBoards(); ++board)
            for (size_t channel = 0; channel < overlay.nChannelsPerBoard(); ++channel)
              for (size_t sample = 0; sample < overlay.nSamplesPerChannel(); ++sample)
              {
                //MF_LOG_VERBATIM("CheckRawDigits")
                //  << "  getting ADC for board " << board << ", channel " << channel << ", sample " << sample;
                hist->Fill(overlay.adc_val(board, channel, sample));
              }
        }
      }
    }
  }

  //------------------------------------------------------------------
  void CheckRawDigits::writeResults(art::Results& r)
  {
    MF_LOG_VERBATIM("CheckRawDigits")
      << "Not writing results (yet)";
  }

  //------------------------------------------------------------------
  void CheckRawDigits::clear()
  {
  }

  DEFINE_ART_RESULTS_PLUGIN(CheckRawDigits)
}// end tcpCompression namespace
