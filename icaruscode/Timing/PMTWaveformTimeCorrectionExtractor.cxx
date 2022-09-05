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


icarus::timing::PMTWaveformTimeCorrectionExtractor::PMTWaveformTimeCorrectionExtractor(
            detinfo::DetectorClocksData const detTimingService,
            icarusDB::IICARUSChannelMap const & channelMapService,
            icarusDB::PMTTimingCorrections & pmtTimingCorrectionsService, 
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
    std::vector<PMTWaveformTimeCorrection> & corrections  )
{

    unsigned int crateSignalID = waveChannelID & 0x00F0;

    std::cout << cateogry << " " << crateSignalID << " " << crateSignalID << std::endl;

    if( fCrateFragmentMap.find(crateSignalID) == fCrateFragmentMap.end() ){ 

        std::cout << "!!! " << cateogry << " " << crateSignalID << " " << crateSignalID << std::endl;

        mf::LogError("icarus::timing::PMTWaveformTimeCorrectionExtractor") << 
            "Invalid special channel number: " << std::hex << waveChannelID;

        throw;
    }

    // This will be the first sample of the falling edge of the special channel signal
    // Which corresponds to the global trigger time. 
    int startSampleSignal = static_cast<int>( getStartSample( wave ) );
    
    // we now access the channels that we need
    for( auto const & crateFragID : fCrateFragmentMap[crateSignalID] ){
      
        for( auto const & mapRow : fChannelMap.getChannelIDPairVec(crateFragID) ){
        
            unsigned int channelID = std::get<1U>(mapRow);

            double cableDelay = 0; 

            if( correctCableDelay ){
                cableDelay = fPMTTimingCorrectionsService.getTriggerCableDelay(channelID);
            }

            double newStartTime = wave.TimeStamp() 
                    + startSampleSignal * fClocksData.OpticalClock().TickPeriod() 
                    + cableDelay; // << The correction is saved already with a minus sign

            corrections[channelID].channelID = waveChannelID;
            corrections[channelID].sample = startSampleSignal * fClocksData.OpticalClock().TickPeriod();
            corrections[channelID].startTime = fClocksData.TriggerTime() - newStartTime;

            if(fVerbose){
                std::cout << channelID                                              << ", " 
                    << wave.TimeStamp()                                             << ", "
                    << startSampleSignal                                            << ", "
                    << fClocksData.OpticalClock().TickPeriod()                      << ", "
                    << startSampleSignal*fClocksData.OpticalClock().TickPeriod()    << ", " 
                    << newStartTime                                                 << ", " 
                    << fClocksData.TriggerTime() - newStartTime                     << std::endl;
            }
        }
    }
}
