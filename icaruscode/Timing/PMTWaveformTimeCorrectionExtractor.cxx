/**
 * @file   icaruscode/Timing/PMTWaveformTimeCorrectionExtractor.cxx
 * @brief  Extract timing correction and adjust waveform starting point.
 * @author Andrea Scarpelli (ascarpell@bnl.gov)
 * @date   June 03, 2022
 */

// ICARUS/SBN libraries
#include "sbnobj/Common/PMT/Data/PMTconfiguration.h" // sbn::PMTconfiguration
#include "icaruscode/Decode/ChannelMapping/IICARUSChannelMap.h"
#include "icaruscode/Timing/PMTTimingCorrections.h"
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
    : cet::exception{ "PMTWaveformTimeCorrectionExtractor", msg }
{}


icarus::timing::PMTWaveformTimeCorrectionExtractor::MultipleCorrectionsForChannel
::MultipleCorrectionsForChannel
  (unsigned int channel, unsigned int existing, unsigned int additional)
    : Error{
      "Attempt to overwrite the correction for channel "
      + std::to_string(channel) + " (from channel " + std::to_string(existing)
      + ") with another (from channel " + std::to_string(additional) + ")\n"
      }
{}


icarus::timing::PMTWaveformTimeCorrectionExtractor::NotASpecialChannel::NotASpecialChannel
  (unsigned int channel)
    : Error{ makeBaseException(channel) }
    {}
    

auto icarus::timing::PMTWaveformTimeCorrectionExtractor::NotASpecialChannel::makeBaseException
  (unsigned int channel) -> Error
{
    return Error{}
      << "PMT readout channel ID " << channel
      << " (0x" << std::hex << channel << ") is not a special channel.\n";
}


icarus::timing::PMTWaveformTimeCorrectionExtractor::UnknownCrate::UnknownCrate
  (unsigned int channel)
    : Error{ makeBaseException(channel) }
    {}
    

auto icarus::timing::PMTWaveformTimeCorrectionExtractor::UnknownCrate::makeBaseException
  (unsigned int channel) -> Error
{
    return Error{}
      << "PMT readout crate for special channel ID " << channel
      << " (0x" << std::hex << channel << ") not known.\n";
}



// -----------------------------------------------------------------------------


icarus::timing::PMTWaveformTimeCorrectionExtractor::PMTWaveformTimeCorrectionExtractor(
            detinfo::DetectorClocksData const detTimingService,
            icarusDB::IICARUSChannelMap const & channelMapService,
            icarusDB::PMTTimingCorrections const* pmtTimingCorrectionsService, 
            bool verbose )
: fClocksData( detTimingService )
, fChannelMap( channelMapService )
, fPMTTimingCorrectionsService( pmtTimingCorrectionsService )
, fVerbose( verbose )
{}


// -----------------------------------------------------------------------------


template<typename T>
  size_t icarus::timing::PMTWaveformTimeCorrectionExtractor::getMinBin( 
        std::vector<T> const& vv, size_t startElement, size_t endElement ){

    auto minel = 
        std::min_element( vv.begin()+startElement, vv.begin()+endElement );
    size_t minsample = std::distance( vv.begin()+startElement, minel );

    return minsample;
}


// -----------------------------------------------------------------------------


template<typename T>
  size_t icarus::timing::PMTWaveformTimeCorrectionExtractor::getMaxBin( 
            std::vector<T> const& vv, size_t startElement, size_t endElement){

    auto maxel = 
        std::max_element( vv.begin()+startElement, vv.begin()+endElement );
    
    size_t maxsample = std::distance( vv.begin()+startElement, maxel );

    return maxsample;
} 


// -----------------------------------------------------------------------------


template<typename T>
  size_t icarus::timing::PMTWaveformTimeCorrectionExtractor::getStartSample( std::vector<T> const& vv ){
    
    // NOTE: when changing this algorithm, also update the documentation
    // in the section "Signal timing extraction" of the class documentation

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
    bool correctCableDelay,
    std::vector<PMTWaveformTimeCorrection> & corrections  ) const
{
    if (!fPMTTimingCorrectionsService && correctCableDelay) {
      throw Error{ "Requested cable delay correction without providing a correction database!\n" };
    }

    unsigned int const waveChannelID = wave.ChannelNumber();
    if ((waveChannelID & 0xF000) == 0)
      throw NotASpecialChannel{ waveChannelID };
    
    unsigned int crateSignalID = waveChannelID & 0x00F0;

    auto const itCrateFragment = fCrateFragmentMap.find(crateSignalID);
    if( itCrateFragment == fCrateFragmentMap.end() )
        throw UnknownCrate{ waveChannelID };

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
      
        for( auto const & mapRow : fChannelMap.getPMTchannelInfo(crateFragID) ){
        
            unsigned int channelID = mapRow.channelID;

            double cableDelay = 0; 

            if( correctCableDelay ){
                cableDelay = fPMTTimingCorrectionsService->getTriggerCableDelay(channelID);
            }

            // time in electronics scale when trigger signal arrived to readout;
            // ideally, it would be fClocksData.TriggerTime()
            double newStartTime = wave.TimeStamp() 
                    + startSampleSignal * fClocksData.OpticalClock().TickPeriod() 
                    + cableDelay; // << The correction is saved already with a minus sign

            PMTWaveformTimeCorrection& correction = correctionFor(channelID);
            if (correction.isValid()) {
              throw MultipleCorrectionsForChannel
                { channelID, correction.channelID, waveChannelID };
            }
            
            correction.channelID = waveChannelID;
            correction.sample = startSampleSignal * fClocksData.OpticalClock().TickPeriod();
            correction.startTime = fClocksData.TriggerTime() - newStartTime;

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
