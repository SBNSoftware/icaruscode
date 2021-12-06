////////////////////////////////////////////////////////////////////////
//
// WaveformIntegrity class - variant of CalWire that deconvolves in
// Regions Of Interest
//
// baller@fnal.gov
//
//
////////////////////////////////////////////////////////////////////////

// C/C++ standard libraries
#include <string>
#include <vector>
#include <utility> // std::pair<>
#include <memory> // std::unique_ptr<>
#include <iomanip>
#include <fstream>
#include <random>

// ROOT libraries
#include "TH1D.h"
#include "TProfile.h"

// framework libraries
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/ReplicatedProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/cpu_timer.h"

// LArSoft libraries
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/ArtDataHelper/WireCreator.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalService.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalProvider.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"

#include "icaruscode/TPC/SignalProcessing/RecoWire/DeconTools/IROIFinder.h"
#include "icaruscode/TPC/SignalProcessing/RecoWire/DeconTools/IDeconvolution.h"
#include "icaruscode/TPC/SignalProcessing/RecoWire/DeconTools/IBaseline.h"
#include "icarus_signal_processing/WaveformTools.h"

#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
#include "tbb/task_arena.h"
#include "tbb/spin_mutex.h"

///creation of calibrated signals on wires
namespace caldata {
    
tbb::spin_mutex deconvolutionSpinMutex;

class WaveformIntegrity : public art::ReplicatedProducer
{
  public:
    // create calibrated signals on wires. this class runs 
    // an fft to remove the electronics shaping.     
    explicit WaveformIntegrity(fhicl::ParameterSet const &, art::ProcessingFrame const&);
    
    void produce(art::Event& evt, art::ProcessingFrame const&) override; 
 //   void beginJob() override;  
 //   void endJob() override;                 
    void reconfigure(fhicl::ParameterSet const& p);
    
  private:

    // It seems there are pedestal shifts that need correcting
    float fixTheFreakingWaveform(const std::vector<float>&, raw::ChannelID_t, std::vector<float>&) const;
    
    float getTruncatedRMS(const std::vector<float>&) const;
    
    std::vector<art::InputTag>                                 fNewRawDigitLabelVec;           ///< New Raw Digits
    std::vector<art::InputTag>                                 fOldRawDigitLabelVec;           ///< From previous run

    icarus_signal_processing::WaveformTools<float>             fWaveformTool;

    const geo::GeometryCore*                                   fGeometry        = lar::providerFrom<geo::Geometry>();
    const lariov::ChannelStatusProvider*                       fChannelFilter   = lar::providerFrom<lariov::ChannelStatusService>();
    const lariov::DetPedestalProvider*                         fPedRetrievalAlg = lar::providerFrom<lariov::DetPedestalService>();
    
}; // class WaveformIntegrity

DEFINE_ART_MODULE(WaveformIntegrity)
  
//-------------------------------------------------
WaveformIntegrity::WaveformIntegrity(fhicl::ParameterSet const & pset, art::ProcessingFrame const& frame) : art::ReplicatedProducer(pset,frame)
{
    this->reconfigure(pset);
}

//////////////////////////////////////////////////////
void WaveformIntegrity::reconfigure(fhicl::ParameterSet const& pset)
{
    // Recover the parameters
    fNewRawDigitLabelVec = pset.get< std::vector<art::InputTag>> ("NewRawDigitLabelVec", {"daqTPC"});
    fOldRawDigitLabelVec = pset.get< std::vector<art::InputTag>> ("OldRawDigitLabelVec", {"daqTPC"});
    
    return;
}
  
//////////////////////////////////////////////////////
void WaveformIntegrity::produce(art::Event& evt, art::ProcessingFrame const& frame)
{
    // We loop over the collection of RawDigits in our input list
    // This is not done multi threaded as a way to cut down on overall job memory usage...
    if (fNewRawDigitLabelVec.size() != fOldRawDigitLabelVec.size()) 
    {
        std::cout << "WaveformIntegrity not given equal size product label vectors" << std::endl;
        return;
    }

    for(size_t prodIdx = 0; prodIdx < fNewRawDigitLabelVec.size(); prodIdx++)
    {
        art::Handle<std::vector<raw::RawDigit>> newRawDigitHandle;

        evt.getByLabel(fNewRawDigitLabelVec[prodIdx],newRawDigitHandle);

        art::Handle<std::vector<raw::RawDigit>> oldRawDigitHandle;

        evt.getByLabel(fOldRawDigitLabelVec[prodIdx],oldRawDigitHandle);

        if (!newRawDigitHandle->size() || !oldRawDigitHandle->size())
        {
            std::cout << "WaveformIntegrity found a zero length buffer" << std::endl;
            continue;
        }

        for(size_t chanIdx = 0; chanIdx < newRawDigitHandle->size(); chanIdx++)
        {
            // get the reference to the current raw::RawDigit
            art::Ptr<raw::RawDigit> newDigitVec(newRawDigitHandle, chanIdx);
            art::Ptr<raw::RawDigit> oldDigitVec(oldRawDigitHandle, chanIdx);

            raw::ChannelID_t newChannel = newDigitVec->Channel();
            raw::ChannelID_t oldChannel = oldDigitVec->Channel();

            if (newChannel != oldChannel)
            {
                std::cout << "WaveformIntegrity finds channel mismatch, idx: " << chanIdx << ", new channel: " << newChannel << ", old channel: " << oldChannel << std::endl;

                continue;
            }

            size_t dataSize = newDigitVec->Samples();

            std::vector<short> newRawADC(dataSize);
            std::vector<short> oldRawADC(dataSize);
    
            // uncompress the data
            raw::Uncompress(newDigitVec->ADCs(), newRawADC, newDigitVec->Compression());
            raw::Uncompress(oldDigitVec->ADCs(), oldRawADC, oldDigitVec->Compression());

            if (newRawADC != oldRawADC)
            {
                std::vector<short> diffVec;
                std::vector<short> idxVec;

                for(size_t tickIdx = 0; tickIdx < newRawADC.size(); tickIdx++)
                {
                    if (newRawADC[tickIdx] != oldRawADC[tickIdx])
                    {
                        diffVec.push_back(newRawADC[tickIdx] - oldRawADC[tickIdx]);
                        idxVec.push_back(tickIdx);
                    }
                }

                short maxDiff = *std::max_element(diffVec.begin(),diffVec.end());
                short minDiff = *std::min_element(diffVec.begin(),diffVec.end());
            
                std::vector<geo::WireID> wireIDVec = fGeometry->ChannelToWire(newChannel);

                std::cout << "==> Channel: " << newChannel << " - " << wireIDVec[0] << " - has " << diffVec.size() << " max/min: " << maxDiff << "/" << minDiff << std::endl;
//                for(size_t idx = 0; diffVec.size(); idx++) std::cout << idxVec[idx] << "/" << diffVec[idx] << " ";
//                std::cout << std::endl;
            }

        }

    }

    return;
} // produce
    
float WaveformIntegrity::getTruncatedRMS(const std::vector<float>& waveform) const
{
    // do rms calculation - the old fashioned way and over all adc values
    std::vector<float> locWaveform = waveform;
    
    // sort in ascending order so we can truncate the sume
    std::sort(locWaveform.begin(), locWaveform.end(),[](const auto& left, const auto& right){return std::fabs(left) < std::fabs(right);});
    
    float threshold = 10.; //fTruncRMSThreshold;
    
    std::vector<float>::iterator threshItr = std::find_if(locWaveform.begin(),locWaveform.end(),[threshold](const auto& val){return std::fabs(val) > threshold;});
    
    //int minNumBins = std::max(int(fTruncRMSMinFraction * locWaveform.size()),int(std::distance(locWaveform.begin(),threshItr)));
    int minNumBins = std::max(int(0.8 * locWaveform.size()),int(std::distance(locWaveform.begin(),threshItr)));
    
    // Get the truncated sum
    float truncRms = std::inner_product(locWaveform.begin(), locWaveform.begin() + minNumBins, locWaveform.begin(), 0.);
    
    truncRms = std::sqrt(std::max(0.,truncRms / double(minNumBins)));
    
    return truncRms;
}
    
float WaveformIntegrity::fixTheFreakingWaveform(const std::vector<float>& waveform, raw::ChannelID_t channel, std::vector<float>& fixedWaveform) const
{
    // do rms calculation - the old fashioned way and over all adc values
    std::vector<float> locWaveform = waveform;
    
    // sort in ascending order so we can truncate the sume
    std::sort(locWaveform.begin(), locWaveform.end(),[](const auto& left, const auto& right){return std::fabs(left) < std::fabs(right);});
    
    // Get the mean of the waveform we're checking...
    float sumWaveform  = std::accumulate(locWaveform.begin(),locWaveform.begin() + locWaveform.size()/2, 0.);
    float meanWaveform = sumWaveform / float(locWaveform.size()/2);
    
    std::vector<float> locWaveformDiff(locWaveform.size()/2);
    
    std::transform(locWaveform.begin(),locWaveform.begin() + locWaveform.size()/2,locWaveformDiff.begin(), std::bind(std::minus<float>(),std::placeholders::_1,meanWaveform));
    
    float localRMS = std::inner_product(locWaveformDiff.begin(), locWaveformDiff.end(), locWaveformDiff.begin(), 0.);
    
    localRMS = std::sqrt(std::max(float(0.),localRMS / float(locWaveformDiff.size())));
    
    float threshold = 6. * localRMS;
    
    std::vector<float>::iterator threshItr = std::find_if(locWaveform.begin(),locWaveform.end(),[threshold](const auto& val){return std::fabs(val) > threshold;});
    
    //int minNumBins = std::max(int(fTruncRMSMinFraction * locWaveform.size()),int(std::distance(locWaveform.begin(),threshItr)));
    int minNumBins = std::max(int(0.8 * locWaveform.size()),int(std::distance(locWaveform.begin(),threshItr)));
    
    // recalculate the mean
    float aveSum      = std::accumulate(locWaveform.begin(), locWaveform.begin() + minNumBins, 0.);
    float newPedestal = aveSum / minNumBins;
    
    // recalculate the rms
    locWaveformDiff.resize(minNumBins);
    
    std::transform(locWaveform.begin(),locWaveform.begin() + minNumBins,locWaveformDiff.begin(), std::bind(std::minus<float>(),std::placeholders::_1,newPedestal));
    
    localRMS = std::inner_product(locWaveform.begin(), locWaveform.begin() + minNumBins, locWaveform.begin(), 0.);
    localRMS = std::sqrt(std::max(float(0.),localRMS / float(minNumBins)));
    
    // Set the waveform to the new baseline
    std::transform(waveform.begin(), waveform.end(), fixedWaveform.begin(), [newPedestal](const auto& val){return val - newPedestal;});
   
    return localRMS;
}

} // end namespace caldata
