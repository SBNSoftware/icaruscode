/**
 * @file   TimingCorrectionExtraction_module.cc
 * @brief  Extract timing correction and adjust waveform starting point.
 * @author Andrea Scarpelli (ascarpell@bnl.gov)
 * @date   June 03, 2022
 */

// ICARUS/SBN libraries
#include "sbnobj/Common/PMT/Data/PMTconfiguration.h" // sbn::PMTconfiguration
#include "icaruscode/Decode/ChannelMapping/IICARUSChannelMap.h"
#include "icaruscode/Timing/PMTTimingCorrections.h"

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

// Database interface helpers
#include "wda.h"

// C/C++ standard libraries
#include <memory> // std::unique_ptr<>
#include <optional>
#include <string>
#include <utility> // std::move()
#include <cassert>
#include <tuple>


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

    fhicl::Atom<bool> CorrectCablesDelay {
        fhicl::Name("CorrectCablesDelay"),
        fhicl::Comment("Use the calibration database to correct for the cables delays"),
        true //default 
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
    
  /// process the run
  void beginRun(art::Run& run) override;

  /// process the event
  void produce(art::Event& event ) override;

private:

  std::vector<art::InputTag> fInputLabels;

  art::InputTag fWaveformsLabel;

  bool fRegenerateWaveforms = true; 

  bool fCorrectCablesDelay;

  bool fVerbose = false; ///< Whether to print the configuration we read.
  
  std::string fLogCategory; ///< Category tag for messages.

  /// Interface to LArSoft configuration for detector timing.
  detinfo::DetectorClocksData const fClocksData;

  /// Pointer to the online channel mapping service.
  icarusDB::IICARUSChannelMap const& fChannelMap;

  /// Pointer to the online pmt corrections service
  icarusDB::PMTTimingCorrections & fPMTTimingCorrectionsService;

  // To be ported to data product ? 
  struct PMTTimeCorrection { 

    double startTime;
    double triggerReferenceDelay; // trigger reference signal cables delay
    double cablesDelay; // PPS reset cables relay
    double laserDelay; // Singal + electron transit time delay

  };

// -----------------------------------------------------------------------------
icarus::TimingCorrectionExtraction::TimingCorrectionExtraction( Parameters const& config ) 
    : art::EDProducer(config)
    , fInputLabels{ config().InputLabels() }
    , fWaveformsLabel{ config().WaveformsLabel() }
    , fRegenerateWaveforms{ config().RegenerateWaveforms() }
    , fCorrectCablesDelay{ config().CorrectCablesDelay() }  
    , fVerbose{ config().Verbose() }
    , fLogCategory{ config().LogCategory() }
    , fClocksData{ art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob() }
    , fChannelMap{ *(art::ServiceHandle<icarusDB::IICARUSChannelMap const>{}) }
    , fPMTTimingCorrectionsService{ *(art::ServiceHandle<icarusDB::PMTTimingCorrections>{}) }
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

/*
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
                                std::vector<PMTTimeCorrection> & corrections 
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
    
    // we now access the channels that we need
    for( auto const & fragId : fCrateFragmentMap[channelID] ){
      
      for( auto const & mapRow : fChannelMap.getChannelIDPairVec(fragId) ){
        
        unsigned int channel_id = std::get<1U>(mapRow);

        double triggerReferenceDelay = fPMTTimingCorrectionsService.getTriggerCableDelay(channel_id);

        double newTriggerTime = wave.TimeStamp() 
                    + startSampleSignal * fClocksData.OpticalClock().TickPeriod() 
                    + triggerReferenceDelay; // << The correction is saved already with a minus sign

        corrections[channel_id].startTime = fClocksData.TriggerTime() - newTriggerTime;

        if(fVerbose){
            std::cout << channel_id                                                   << ", " 
                      << wave.TimeStamp()                                             << ", "
                      << startSampleSignal                                            << ", "
                      << fClocksData.OpticalClock().TickPeriod()                      << ", "
                      << startSampleSignal*fClocksData.OpticalClock().TickPeriod()    << ", " 
                      << newTriggerTime                                               << ", " 
                      << fClocksData.TriggerTime() - newTriggerTime                   << std::endl;
        }
      }
    }

} // icarus::TimingCorrectionExtraction::findTimeCorrection 
*/

// -----------------------------------------------------------------------------
void icarus::TimingCorrectionExtraction::beginRun( art::Run& run ) {

    fPMTTimingCorrectionsService.readTimeCorrectionDatabase(run);

}

// -----------------------------------------------------------------------------
void icarus::TimingCorrectionExtraction::produce( art::Event& event ) {

    std::vector<PMTTimeCorrection> corrections;
    corrections.resize(360);

    if( !fInputLabels.empty() ){
        
        for( size_t iTag=0; iTag<fInputLabels.size(); iTag++ ){

            art::InputTag label = fInputLabels[iTag];

            art::Handle<std::vector<raw::OpDetWaveform>> specialWaveformHandle;
            event.getByLabel( label, specialWaveformHandle );

            if( specialWaveformHandle.isValid() && !specialWaveformHandle->empty() ){

                for(  auto const & specialWaveform : *specialWaveformHandle  ){

                     // Extract the time corrections from the special waveforms
                    findTimeCorrection( specialWaveform, corrections );

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
    std::vector<raw::OpDetWaveform> correctedWaveforms;

    // Now we read the opdet waveform collection and re-apply the timing corrections 
    // A new collection is hence created 
    art::Handle<std::vector<raw::OpDetWaveform>> waveformHandle;
    event.getByLabel( fWaveformsLabel, waveformHandle );

    if( waveformHandle.isValid() && !waveformHandle->empty() ){

        for ( size_t iWave=0; iWave<waveformHandle->size(); iWave++ ){

            // Not-const copy of the previous object
            raw::OpDetWaveform waveform = (*waveformHandle)[iWave];

            auto channel_id = waveform.ChannelNumber();

            // Recalculate the timestamp
            raw::TimeStamp_t correctT = 
                waveform.TimeStamp() + corrections[channel_id].startTime;

            /*
            if ( fCorrectCablesDelay ){
                correctT += 
                    ( corrections[channel_id].cablesDelay - corrections[channel_id].laserDelay );
            }
            */

            waveform.SetTimeStamp( correctT );

            correctedWaveforms.push_back( waveform );

            //std::cout << (*waveformHandle)[iWave].TimeStamp() << " " << waveform.TimeStamp() << std::endl;

        }

    } else {
        mf::LogError(fLogCategory) 
            << " Not found OpDetWaveform data product with label '" 
            << fWaveformsLabel.encode() << "'" << std::endl;
        throw;
    }

    // The new waveform collection is also saved in the event stream
    event.put(
      std::make_unique<std::vector<raw::OpDetWaveform>>(std::move(correctedWaveforms)) 
    );

} //icarus::TimingCorrectionExtraction::produce


// -----------------------------------------------------------------------------
DEFINE_ART_MODULE(icarus::TimingCorrectionExtraction)


// -----------------------------------------------------------------------------