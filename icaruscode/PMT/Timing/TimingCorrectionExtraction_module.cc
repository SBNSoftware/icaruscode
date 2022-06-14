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

    fhicl::Atom<std::string> DatabaseURL {
        fhicl::Name("DatabaseURL"),
        fhicl::Comment("The entry point url for the calibration database"),
        "https://dbdata0vm.fnal.gov:9443/icarus_con_prod/app/data?" //default
    };

    fhicl::Atom<unsigned int> Timeout {
        fhicl::Name("Timeout"),
        fhicl::Comment("Database query timeout"),
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

  std::string fUrl;

  unsigned int fTimeout;

  bool fCorrectCablesDelay;

  bool fVerbose = false; ///< Whether to print the configuration we read.
  
  std::string fLogCategory; ///< Category tag for messages.

  /// Interface to LArSoft configuration for detector timing.
  detinfo::DetectorClocksData const fClocksData;

  /// Pointer to the online channel mapping service.
  icarusDB::IICARUSChannelMap const& fChannelMap;

  // To be ported to data product ? 
  struct PMTTimeCorrection { 

    double startTime;
    double triggerReferenceDelay; // trigger reference signal cables delay
    double cablesDelay; // PPS reset cables relay
    double laserDelay; // Singal + electron transit time delay

  };

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

  std::map< unsigned int, std::tuple<double, double, double> > fDatabaseTimingCorrections;

  int GetDataset(
    const std::string& name, const uint32_t &run, Dataset& dataset ) const;

  void GetCablesDelay( const uint32_t & run ) ;

  void GetLaserDelay( const uint32_t & run ) ;

  template<typename T>
    size_t getMaxBin( std::vector<T> vv, size_t startElement, size_t endElement);

  template<typename T>
    size_t getMinBin( std::vector<T> vv, size_t startElement, size_t endElement);

  template<typename T>
    size_t getStartSample( std::vector<T> vv );

  void findTimeCorrection( 
    raw::OpDetWaveform const & wave, 
    std::vector<PMTTimeCorrection> & corrs ) ;

}; // icarus::TimingCorrectionExtractor


// -----------------------------------------------------------------------------
icarus::TimingCorrectionExtraction::TimingCorrectionExtraction( Parameters const& config ) 
    : art::EDProducer(config)
    , fInputLabels{ config().InputLabels() }
    , fWaveformsLabel{ config().WaveformsLabel() }
    , fRegenerateWaveforms{ config().RegenerateWaveforms() }
    , fUrl{ config().DatabaseURL() }
    , fTimeout{ config().Timeout() }
    , fCorrectCablesDelay{ config().CorrectCablesDelay() }  
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
int icarus::TimingCorrectionExtraction::GetDataset(
    const std::string& name, const uint32_t &run, Dataset& dataset ) const {

    int error = 0;
    std::string url = fUrl + "f="+name +"&t="+std::to_string(run);
    dataset = getDataWithTimeout( url.c_str(), "", fTimeout, &error );
    if ( error ){
        throw cet::exception(fLogCategory) 
          << "Calibration database access failed. URL (" << url 
          << ") Error code: " << error;
    }
    if ( getHTTPstatus(dataset) !=200 ){
        throw cet::exception(fLogCategory)
            << "Calibration database access failed. URL (" << url
            << "). HTTP error status: " << getHTTPstatus(dataset) 
            << ". HTTP error message: " << getHTTPmessage(dataset); 
    }

    return error;

}


// -----------------------------------------------------------------------------
void icarus::TimingCorrectionExtraction::GetCablesDelay( const uint32_t & run ) { 

    // pmt_cables_delay: delays of the cables relative to trigger 
    // and reset distribution
    const std::string name("pmt_cables_delays_data");
    Dataset dataset;
    int error = GetDataset( name, run, dataset );

    if (error) throw(std::exception());

    for( int row=4; row < getNtuples(dataset); row++ ) {

        Tuple tuple = getTuple(dataset, row);

        if( tuple != NULL ) {

            unsigned int channel_id = getLongValue( tuple, 0, &error );
            if( error ) throw std::runtime_error( "Encountered error while trying to access 'channel_id' on table " + name );
            
            // PPS reset correction
            double reset_distribution_delay = getDoubleValue( tuple, 11, &error );
            if( error ) throw std::runtime_error( "Encountered error '" + std::to_string(error) + "' while trying to access 'reset_distribution_delay' on table " + name );

            // Trigger cable delay
            double trigger_reference_delay = getDoubleValue( tuple, 10, &error );
            if( error ) throw std::runtime_error( "Encountered error while trying to access 'trigger_reference_delay' on table " + name );

            // Phase correction
            double phase_correction = getDoubleValue( tuple, 13, &error );
            if( error ) throw std::runtime_error( "Encountered error while trying to access 'phase_correction' on table " + name );

            std::tuple rowValues = std::make_tuple( trigger_reference_delay/1000., (reset_distribution_delay-phase_correction)/1000., 0 );
            fDatabaseTimingCorrections[channel_id] = rowValues;

            releaseTuple(tuple);
        }

    }

}


void icarus::TimingCorrectionExtraction::GetLaserDelay( const uint32_t & run ) { 

    // pmt_laser_delay: delay from the Electron Transit time inside the PMT 
    // and the PMT signal cable 

    const std::string name("pmt_laser_timing_data");
    Dataset dataset;
    int error = GetDataset( name, run, dataset );

    if (error) throw(std::exception());

    for( int row=4; row < getNtuples(dataset); row++ ) {

        Tuple tuple = getTuple(dataset, row);
        if( tuple != NULL ) {

            unsigned int channel_id = getLongValue( tuple, 0, &error );
            if( error ) throw std::runtime_error( "Encountered error while trying to access channel_id on table " + name );

            double t_signal = getDoubleValue( tuple, 5, &error );
            if( error ) throw std::runtime_error( "Encountered error while trying to access 't_signal' on table " + name );

            std::get<2u>(fDatabaseTimingCorrections[channel_id]) = t_signal/1000.;

            releaseTuple(tuple);
        }

    }


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
    
    double newTriggerTime = wave.TimeStamp() 
                    + startSampleSignal * fClocksData.OpticalClock().TickPeriod();

    // we now access the channels that we need
    for( auto const & fragId : fCrateFragmentMap[channelID] ){
      
      for( auto const & mapRow : fChannelMap.getChannelIDPairVec(fragId) ){
        
        unsigned int channel_id = std::get<1U>(mapRow);

        corrections[channel_id].triggerReferenceDelay = std::get<0U>(fDatabaseTimingCorrections[channel_id]);

        corrections[channel_id].cablesDelay = std::get<1U>(fDatabaseTimingCorrections[channel_id]);

        corrections[channel_id].laserDelay = std::get<2U>(fDatabaseTimingCorrections[channel_id]);

        corrections[channel_id].startTime = fClocksData.TriggerTime() - 
                            ( newTriggerTime + corrections[channel_id].triggerReferenceDelay );

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


// -----------------------------------------------------------------------------
void icarus::TimingCorrectionExtraction::beginRun( art::Run& run ) {

    GetCablesDelay( run.id().run() );
    GetLaserDelay ( run.id().run() );

    if( fVerbose ) {

        std::cout << "Dump information from database " << std::endl;
        std::cout << "channel_id, trigger_reference_delay, reset_distribution_delay-phase_shift, t_signal" << std::endl;

        for( auto const & [key, values] : fDatabaseTimingCorrections ){
            std::cout << key << " " 
                  << std::get<0U>(values) << "," 
                  << std::get<1U>(values) << ", " 
                  << std::get<2U>(values) << ", "
                  << std::endl; 
        }
    }

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

            if ( fCorrectCablesDelay ){
                correctT += 
                    ( corrections[channel_id].cablesDelay - corrections[channel_id].laserDelay );
            }

            waveform.SetTimeStamp( correctT );

            correctedWaveforms.push_back( waveform );
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