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
#include "icaruscode/Timing/DataProducts/PMTWaveformTimeCorrection.h"
#include "icaruscode/Timing/PMTWaveformTimeCorrectionExtractor.h"

// framework libraries
#include "canvas/Persistency/Provenance/ProcessConfiguration.h"
#include "canvas/Persistency/Provenance/ProcessHistory.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"

#include "messagefacility/MessageLogger/MessageLogger.h"
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
#include <tuple>

// -----------------------------------------------------------------------------
icarus::timing::PMTWaveformTimeCorrectionExtractor::Error::Error
  (std::string const& msg /* == "" */)
    : cet::exception{
      "PMTWaveformTimeCorrectionExtractor::MultipleCorrectionsForChannel", msg
      }
{}


icarus::timing::PMTWaveformTimeCorrectionExtractor::MultipleCorrectionsForChannel
::MultipleCorrectionsForChannel
  (unsigned int existing, unsigned int additional)
    : Error{
      "Attempt to overwrite correction from channel "
      + std::to_string(existing) + " with one from channel "
      + std::to_string(additional) + "\n"
      }
{}


// -----------------------------------------------------------------------------


icarus::timing::PMTWaveformTimeCorrectionExtractor::PMTWaveformTimeCorrectionExtractor(
            detinfo::DetectorClocksData const detTimingService,
            icarusDB::IICARUSChannelMap const & channelMapService,
            icarusDB::PMTTimingCorrections const& pmtTimingCorrectionsService, 
            bool const & verbose )
: fClocksData( detTimingService )
, fChannelMap( channelMapService )
, fPMTTimingCorrectionsService( pmtTimingCorrectionsService )
, fVerbose( verbose )
{}


// -----------------------------------------------------------------------------


template<typename T>
  size_t icarus::timing::PMTWaveformTimeCorrectionExtractor::getMinBin( 
        std::vector<T> vv, size_t startElement, size_t endElement ){

    auto minel = 
        std::min_element( vv.begin()+startElement, vv.begin()+endElement );
    size_t minsample = std::distance( vv.begin()+startElement, minel );

    return minsample;
}


// -----------------------------------------------------------------------------


template<typename T>
  size_t icarus::timing::PMTWaveformTimeCorrectionExtractor::getMaxBin( 
            std::vector<T> vv, size_t startElement, size_t endElement){

    auto maxel = 
        std::max_element( vv.begin()+startElement, vv.begin()+endElement );
    
    size_t maxsample = std::distance( vv.begin()+startElement, maxel );

    return maxsample;
} 


// -----------------------------------------------------------------------------


template<typename T>
  size_t icarus::timing::PMTWaveformTimeCorrectionExtractor::getStartSample( std::vector<T> vv ){

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


//---------------------------------------------------------------------------------------


void icarus::timing::PMTWaveformTimeCorrectionExtractor::findWaveformTimeCorrections
(   raw::OpDetWaveform const & wave,
    std::string const & cateogry, 
    unsigned int const & waveChannelID, 
    bool const & correctCableDelay,
    std::vector<PMTWaveformTimeCorrection> & corrections  ) const
{

    unsigned int crateSignalID = waveChannelID & 0x00F0;

    auto const itCrateFragment = fCrateFragmentMap.find(crateSignalID);
    if( itCrateFragment == fCrateFragmentMap.end() ){ 

        mf::LogError("icarus::timing::PMTWaveformTimeCorrectionExtractor") << 
            "Invalid special channel number: " << std::hex << waveChannelID << 
            " And category: " << cateogry;

        throw; // FIXME
    }

    // This will be the first sample of the falling edge of the special channel signal
    // Which corresponds to the global trigger time. 
    int startSampleSignal = static_cast<int>( getStartSample( wave ) );
    
    // allocates room for correction for `channel`; intermediate ones are defaulted
    auto correctionFor
      = [&corrections](unsigned int channel) -> PMTWaveformTimeCorrection& 
      {
        if (channel >= corrections.size()) corrections.resize(channel + 1);
        return corrections[channel];
      };
    
    // we now access the channels that we need
    for( auto const & crateFragID : itCrateFragment->second ){
      
        for( auto const & mapRow : fChannelMap.getChannelIDPairVec(crateFragID) ){
        
            unsigned int channelID = std::get<1U>(mapRow);

            double cableDelay = 0; 

            if( correctCableDelay ){
                cableDelay = fPMTTimingCorrectionsService.getTriggerCableDelay(channelID);
            }

            // time in electronics scale when trigger signal arrived to readout;
            // ideally, it would be fClocksData.TriggerTime()
            double newStartTime = wave.TimeStamp() 
                    + startSampleSignal * fClocksData.OpticalClock().TickPeriod() 
                    + cableDelay; // << The correction is saved already with a minus sign

            PMTWaveformTimeCorrection& correction = correctionFor(channelID);
            if (correction.isValid()) {
              throw MultipleCorrectionsForChannel
                { correction.channelID, waveChannelID };
            }
            
            correction.channelID = waveChannelID;
            correction.sample = startSampleSignal * fClocksData.OpticalClock().TickPeriod();
            correction.startTime = fClocksData.TriggerTime() - newStartTime;
            // TODO add check that we are not overwriting a correction

            if(fVerbose){
                std::cout << channelID                                              << ", " 
                    << wave.TimeStamp()                                             << ", "
                    << startSampleSignal                                            << ", "
                    << fClocksData.OpticalClock().TickPeriod()                      << ", "
                    << startSampleSignal*fClocksData.OpticalClock().TickPeriod()    << ", " 
                    << newStartTime                                                 << ", " 
                    << cableDelay                                                   << ", "
                    << fClocksData.TriggerTime() - newStartTime                     << std::endl;
            }
        }
    }
}
