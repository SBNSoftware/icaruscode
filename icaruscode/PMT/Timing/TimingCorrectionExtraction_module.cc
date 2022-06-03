/**
 * @file   TimingCorrectionExtraction_module.cc
 * @brief  Extract timing correction and adjust waveform starting point.
 * @author Andrea Scarpelli (ascarpell@bnl.gov)
 * @date   June 03, 2022
 */

// ICARUS/SBN libraries
#include "sbnobj/Common/PMT/Data/PMTconfiguration.h" // sbn::PMTconfiguration
#include "icaruscode/Decode/ChannelMapping/IICARUSChannelMap.h"
//#include "icaruscode/PMT/Timing/PMTTimeCorrection.h" // data product holding the corrections

// framework libraries
#include "canvas/Persistency/Provenance/ProcessConfiguration.h"
#include "canvas/Persistency/Provenance/ProcessHistory.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"

#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/FileBlock.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "cetlib_except/exception.h"

// LArSoft libraries
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataalg/DetectorInfo/DetectorTimings.h"
#include "lardataalg/DetectorInfo/DetectorClocks.h"
#include "larcorealg/CoreUtils/enumerate.h"
#include "larcorealg/CoreUtils/counter.h"

// C/C++ standard libraries
#include <memory> // std::unique_ptr<>
#include <optional>
#include <string>
#include <utility> // std::move()
#include <cassert>


// -----------------------------------------------------------------------------
namespace icarus { class TimingCorrectionExtraction; }
/**
 * @brief 
 * 
 * This module reads 
 * 
 * Input
 * ------
 * 
 * 
 * Output
 * -------
 * 
 * 
 * 
 * Configuration parameters
 * -------------------------
 * 
 * 
 * 
 * Multithreading
 * ---------------
 * 
 * 
 */
class icarus::TimingCorrectionExtraction: public art::EDProducer {
  
public:
  
  /// Configuration of the module.
  struct Config {

    fhicl::Sequence<art::InputTag> InputLabels {
        fhicl::Name("InputLabels"), 
        fhicl::Comment("list of the input lables to be used")
    };

    fhicl::Atom<art::InputTag> WaveformsLabel {
        fhicl::Name("WaveformsLabel"),
        fhicl::Comment("")
    };

    fhicl::Atom<bool> RegenerateWaveforms {
        fhicl::Name("RegenerateWaveforms"),
        fhicl::Comment("")
    };

    fhicl::Atom<bool> Verbose {
      fhicl::Name("Verbose"),
      fhicl::Comment("print the times read and the associated channels"),
      false // default
    };
    
    fhicl::Atom<std::string> LogCategory {
      fhicl::Name("LogCategory"),
      fhicl::Comment("category tag used for messages to message facility"),
      "TimingCorrectionExtraction" // default
    };
    
  }; // struct Config

  using Parameters = art::EDProducer::Table<Config>;

  /// Constructor: just reads the configuration.
  explicit TimingCorrectionExtraction(Parameters const& config);
    
  /// process the event
  void produce(art::Event& event ) override;

private:
  
  std::vector<art::InputTag> fInputLabels;

  art::InputTag fWaveformsLabel;

  bool fRegenerateWaveforms = true; 

  bool fVerbose = false; ///< Whether to print the configuration we read.
  
  std::string fLogCategory; ///< Category tag for messages.

  /// Interface to LArSoft configuration for detector timing.
  detinfo::DetectorClocksData const fClocksData;

  /// Pointer to the online channel mapping service.
  icarusDB::IICARUSChannelMap const& fChannelMap;

  // To be ported to data product ? 
  struct PMTTimeCorrection { 

    std::map<std::string, double> timeCorrection;

  };

  std::vector<PMTTimeCorrection> corrections;

  std::map<unsigned int, std::vector<unsigned int>> fCrateFragmentMap {

    {0x1070 , { 0 , 1 , 2  }}, 
    {0x1060 , { 3 , 4 , 5  }}, 
    {0x1050 , { 6 , 7 , 8  }}, 
    {0x1040 , { 9 , 10, 11 }},
    {0x1030 , { 12, 13, 14 }}, 
    {0x1020 , { 15, 16, 17 }}, 
    {0x1010 , { 18, 19, 20 }}, 
    {0x1000 , { 21, 22, 23 }}, 

  }; 

  template<typename T>
    size_t getMaxBin( std::vector<T> vv, size_t startElement, size_t endElement);

  template<typename T>
    size_t getMinBin( std::vector<T> vv, size_t startElement, size_t endElement);

  template<typename T>
    size_t getStartSample( std::vector<T> vv );

  void findTimeCorrection( 
    raw::OpDetWaveform const & wave, 
    std::vector<PMTTimeCorrection> & corrs, 
    std::string const & instance ) ;


}; // icarus::TimingCorrectionExtractor


// -----------------------------------------------------------------------------
icarus::TimingCorrectionExtraction::TimingCorrectionExtraction( Parameters const& config ) 
    : art::EDProducer(config)
    , fInputLabels{ config().InputLabels() }
    , fWaveformsLabel{ config().WaveformsLabel() }
    , fRegenerateWaveforms{ config().RegenerateWaveforms() }
    , fVerbose{ config().Verbose() }
    , fLogCategory{ config().LogCategory() }
    , fClocksData
        { art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob() }
    , fChannelMap{ *(art::ServiceHandle<icarusDB::IICARUSChannelMap const>{}) }
{

    /// Consumes
    for ( auto const & tag : fInputLabels )
        consumes<std::vector<raw::OpDetWaveform>>(tag);
    if( fRegenerateWaveforms )
        consumes<std::vector<raw::OpDetWaveform>>(fWaveformsLabel);

    /// Produces
    //produces<icarus::PMTTimeCorrection>();
    if ( fRegenerateWaveforms )
        produces<std::vector<raw::OpDetWaveform>>();

}


// -----------------------------------------------------------------------------
template<typename T>
  size_t icarus::TimingCorrectionExtraction::getMinBin( 
        std::vector<T> vv, size_t startElement, size_t endElement ){

    auto minel = 
        std::min_element( vv.begin()+startElement, vv.begin()+endElement );
    size_t minsample = std::distance( vv.begin()+startElement, minel );

    return minsample;
}


// -----------------------------------------------------------------------------
template<typename T>
  size_t icarus::TimingCorrectionExtraction::getMaxBin( 
            std::vector<T> vv, size_t startElement, size_t endElement){

    auto maxel = 
        std::max_element( vv.begin()+startElement, vv.begin()+endElement );
    
    size_t maxsample = std::distance( vv.begin()+startElement, maxel );

    return maxsample;
} 


// -----------------------------------------------------------------------------
template<typename T>
  size_t icarus::TimingCorrectionExtraction::getStartSample( std::vector<T> vv ){

    // We are thinking in inverted polarity
      // Todo: restric the reagion over wich the max is searched: 
      // depends on the readout settings 
    size_t minbin = getMinBin( vv, 0, vv.size() );

    //Search only a cropped region of the waveform backward from the min
    size_t maxbin =  minbin-20; //getMaxBin( wave, minbin-20, minbin );

    // Now we crawl betweem maxbin and minbin and we stop when:
      //( maxbin value - minbin value )*0.05 > (maxbin value - bin value )
    size_t startbin = maxbin;
    auto delta = vv[maxbin]-vv[minbin];
    for( size_t bin=maxbin; bin<minbin; bin++ ){
      auto val = vv[maxbin]-vv[bin];
      if( val >= 0.2*delta ){
        startbin = bin - 1;
        break;
      }
    }

    if( startbin < maxbin ){
      startbin=maxbin;
    }

    return startbin;
}


// -----------------------------------------------------------------------------
void icarus::TimingCorrectionExtraction::findTimeCorrection( 
                                raw::OpDetWaveform const & wave, 
                                std::vector<PMTTimeCorrection> & corrs, 
                                std::string const & instance 
) {

    unsigned int channelID = wave.ChannelNumber() & 0xFFF0;

    if( fCrateFragmentMap.find(channelID) == fCrateFragmentMap.end() ){ 

        mf::LogError(fLogCategory) << 
            "Invalid channel number. Verify fCrateFragmentMap" << std::endl;
        throw;
    }

    // This will be the first sample of the falling edge of the special channel signal
    // Which corresponds to the global trigger time. 
    int startSampleSignal = static_cast<int>( getStartSample( wave ) );
    
    double newTriggerTime = wave.TimeStamp() 
                                + startSampleSignal * fClocksData.OpticalClock().TickPeriod();

    // we now access the channels that we need
    for( auto const & fragId : fCrateFragmentMap[channelID] ){
      for( auto const & mapRow : fChannelMap.getChannelIDPairVec(fragId) ){
        
        corrs[std::get<1U>(mapRow)].timeCorrection[instance] = 
                                        fClocksData.TriggerTime() - newTriggerTime;

        if(fVerbose){
            std::cout << std::get<1U>(mapRow)                                         << ", " 
                      << instance                                                     << ", " 
                      << startSampleSignal*fClocksData.OpticalClock().TickPeriod()    << ", " 
                      << newTriggerTime                                               << ", " 
                      << fClocksData.TriggerTime() - newTriggerTime                   << std::endl;
        }

      }
    }

} // icarus::TimingCorrectionExtraction::findTimeCorrection 


// -----------------------------------------------------------------------------
void icarus::TimingCorrectionExtraction::produce( art::Event& event ) {

    corrections.reserve(360); // book enough room on the corrections array

    if( !fInputLabels.empty() ){
        
        for( size_t iTag=0; iTag<fInputLabels.size(); iTag++ ){

            art::InputTag label = fInputLabels[iTag];

            art::Handle<std::vector<raw::OpDetWaveform>> specialWaveformHandle;
            event.getByLabel( label, specialWaveformHandle );

            if( specialWaveformHandle.isValid() && !specialWaveformHandle->empty() ){

                for(  auto const & specialWaveform : *specialWaveformHandle  ){

                     // Extract the time corrections from the special waveforms
                    findTimeCorrection( specialWaveform, corrections, label.instance() );

                }

            } else {
                mf::LogError(fLogCategory) 
                    << " Not found OpDetWaveform data product with label '" 
                    << label.encode() << "'" << std::endl;
                throw;
            }

        }
    } else {
        mf::LogError(fLogCategory) 
            << " InputLabels array should contain more than 1 valid entry " << std::endl;
        throw;
    }

    // Create a copy of the waveforms 
    std::vector<raw::OpDetWaveform> corrWaveforms;

    // Now we read the opdet waveform collection and re-apply the timing corrections 
    // A new collection is hence created 
    art::Handle<std::vector<raw::OpDetWaveform>> waveformHandle;
    event.getByLabel( fWaveformsLabel, waveformHandle );

    if( waveformHandle.isValid() && !waveformHandle->empty() ){

        for ( size_t iWave=0; iWave<waveformHandle->size(); iWave++ ){

            // Not-const copy of the previous object
            raw::OpDetWaveform waveform = (*waveformHandle)[iWave];

            // Recalculate the timestamp
            raw::TimeStamp_t tCorrection = 
                waveform.TimeStamp() - corrections[ waveform.ChannelNumber() ].timeCorrection["trgprim"];
            waveform.SetTimeStamp( tCorrection );

            corrWaveforms.push_back( waveform );
        }

    } else {
        mf::LogError(fLogCategory) 
            << " Not found OpDetWaveform data product with label '" 
            << fWaveformsLabel.encode() << "'" << std::endl;
        throw;
    }

    // The new waveform collection is also saved in the event stream
    event.put(
      std::make_unique<std::vector<raw::OpDetWaveform>>(std::move(corrWaveforms)) 
    );

    corrections.clear(); // Free the memory 

} //produce


// -----------------------------------------------------------------------------
DEFINE_ART_MODULE(icarus::TimingCorrectionExtraction)


// -----------------------------------------------------------------------------