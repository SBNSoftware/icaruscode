////////////////////////////////////////////////////////////////////////
//
// Decon1DROI class - variant of CalWire that deconvolves in
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
#include "fhiclcpp/ParameterSet.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h" 
#include "art/Framework/Principal/Handle.h" 
#include "art/Utilities/make_tool.h"
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art_root_io/TFileService.h" 
#include "canvas/Utilities/Exception.h"

// LArSoft libraries
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RecoBase/Wire.h"
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

class Decon1DROI : public art::EDProducer
{
  public:
    // create calibrated signals on wires. this class runs 
    // an fft to remove the electronics shaping.     
    explicit Decon1DROI(fhicl::ParameterSet const& pset);
    virtual ~Decon1DROI();
    
    void produce(art::Event& evt); 
    void beginJob(); 
    void endJob();                 
    void reconfigure(fhicl::ParameterSet const& p);
    
  private:

    // Define a class to handle processing for individual threads
    class multiThreadDeconvolutionProcessing 
    {
    public:
        multiThreadDeconvolutionProcessing(Decon1DROI const&                        parent,
                                           art::Event&                              event,
                                           art::Handle<std::vector<raw::RawDigit>>& rawDigitHandle, 
                                           std::vector<recob::Wire>&                wireColVec,
                                           art::Assns<raw::RawDigit,recob::Wire>&   wireAssns)
            : fDecon1DROI(parent),
              fEvent(event),
              fRawDigitHandle(rawDigitHandle),
              fWireColVec(wireColVec),
              fWireAssns(wireAssns)
        {}

        void operator()(const tbb::blocked_range<size_t>& range) const
        {
            for (size_t idx = range.begin(); idx < range.end(); idx++)
                fDecon1DROI.processChannel(idx, fEvent, fRawDigitHandle, fWireColVec, fWireAssns);
        }
    private:
        const Decon1DROI&                        fDecon1DROI;
        art::Event&                              fEvent;
        art::Handle<std::vector<raw::RawDigit>>& fRawDigitHandle;
        std::vector<recob::Wire>&                fWireColVec;
        art::Assns<raw::RawDigit,recob::Wire>&   fWireAssns;
    };

    // It seems there are pedestal shifts that need correcting
    float fixTheFreakingWaveform(const std::vector<float>&, raw::ChannelID_t, std::vector<float>&) const;
    
    float getTruncatedRMS(const std::vector<float>&) const;

    // Function to do the work
    void  processChannel(size_t,
                         art::Event&,
                         art::Handle<std::vector<raw::RawDigit>>, 
                         std::vector<recob::Wire>&, 
                         art::Assns<raw::RawDigit,recob::Wire>&) const;
    
    std::string                                                fDigitModuleLabel;           ///< module that made digits
    std::string                                                fSpillName;                  ///< nominal spill is an empty string
                                                                                         ///< it is set by the DigitModuleLabel
                                                                                         ///< ex.:  "daq:preSpill" for prespill data
    unsigned short                                             fNoiseSource;                ///< Used to determine ROI threshold
    int                                                        fSaveWireWF;                 ///< Save recob::wire object waveforms
    size_t                                                     fEventCount;                 ///< count of event processed
    int                                                        fMinAllowedChanStatus;       ///< Don't consider channels with lower status
    
    float                                                      fTruncRMSThreshold;          ///< Calculate RMS up to this threshold...
    float                                                      fTruncRMSMinFraction;        ///< or at least this fraction of time bins
    bool                                                       fOutputHistograms;           ///< Output histograms?
    
    std::vector<std::unique_ptr<icarus_tool::IROIFinder>>      fROIFinderVec;               ///< ROI finders per plane
    std::unique_ptr<icarus_tool::IDeconvolution>               fDeconvolution;
    std::unique_ptr<icarus_tool::IBaseline>                    fBaseline;

    icarus_signal_processing::WaveformTools<float>                        fWaveformTool;

    const geo::GeometryCore*                                   fGeometry        = lar::providerFrom<geo::Geometry>();
    const lariov::ChannelStatusProvider*                       fChannelFilter   = lar::providerFrom<lariov::ChannelStatusService>();
    const lariov::DetPedestalProvider*                         fPedRetrievalAlg = lar::providerFrom<lariov::DetPedestalService>();
    
    // Define here a temporary set of histograms...
    std::vector<TH1F*>     fPedestalOffsetVec;
    std::vector<TH1F*>     fFullRMSVec;
    std::vector<TH1F*>     fTruncRMSVec;
    std::vector<TH1F*>     fNumTruncBinsVec;
    std::vector<TProfile*> fPedByChanVec;
    std::vector<TProfile*> fTruncRMSByChanVec;
    std::vector<TH1F*>     fNumROIsHistVec;
    std::vector<TH1F*>     fROILenHistVec;
    
}; // class Decon1DROI

DEFINE_ART_MODULE(Decon1DROI)
  
//-------------------------------------------------
Decon1DROI::Decon1DROI(fhicl::ParameterSet const& pset) : EDProducer{pset}
{
  this->reconfigure(pset);

  produces< std::vector<recob::Wire> >(fSpillName);
  produces<art::Assns<raw::RawDigit, recob::Wire>>(fSpillName);
}

//-------------------------------------------------
Decon1DROI::~Decon1DROI()
{
}

//////////////////////////////////////////////////////
void Decon1DROI::reconfigure(fhicl::ParameterSet const& pset)
{
    // Recover the parameters
    fDigitModuleLabel     = pset.get< std::string >   ("DigitModuleLabel", "daq");
    fNoiseSource          = pset.get< unsigned short >("NoiseSource",          3);
    fSaveWireWF           = pset.get< int >           ("SaveWireWF"             );
    fMinAllowedChanStatus = pset.get< int >           ("MinAllowedChannelStatus");
    fTruncRMSThreshold    = pset.get< float >         ("TruncRMSThreshold",    6.);
    fTruncRMSMinFraction  = pset.get< float >         ("TruncRMSMinFraction", 0.6);
    fOutputHistograms     = pset.get< bool  >         ("OutputHistograms",   true);
    
    // Recover the vector of fhicl parameters for the ROI tools
    const fhicl::ParameterSet& roiFinderTools = pset.get<fhicl::ParameterSet>("ROIFinderToolVec");
    
    fROIFinderVec.resize(roiFinderTools.get_pset_names().size());

    unsigned short roiPadding(std::numeric_limits<unsigned short>::max());
    
    for(const std::string& roiFinderTool : roiFinderTools.get_pset_names())
    {
        const fhicl::ParameterSet& roiFinderToolParamSet = roiFinderTools.get<fhicl::ParameterSet>(roiFinderTool);
        size_t                     planeIdx              = roiFinderToolParamSet.get<size_t>("Plane");

        roiPadding = std::min(roiPadding,roiFinderToolParamSet.get< std::vector<unsigned short>>("roiLeadTrailPad")[0]);
        
        fROIFinderVec.at(planeIdx) = art::make_tool<icarus_tool::IROIFinder>(roiFinderToolParamSet);
    }
    
    std::sort(fROIFinderVec.begin(),fROIFinderVec.end(),[](const auto& left,const auto& right){return left->plane() < right->plane();});

    fDeconvolution = art::make_tool<icarus_tool::IDeconvolution>(pset.get<fhicl::ParameterSet>("Deconvolution"));
    
    // Recover the baseline tool 
    fhicl::ParameterSet baselineParams = pset.get<fhicl::ParameterSet>("Baseline");

    // Check if we need to set the length for setting the baseline
    if (baselineParams.has_key("MaxROILength")) baselineParams.put_or_replace("MaxROILength",size_t(roiPadding));

    fBaseline  = art::make_tool<icarus_tool::IBaseline> (baselineParams);

    fSpillName.clear();
    
    size_t pos = fDigitModuleLabel.find(":");
    if( pos!=std::string::npos )
    {
        fSpillName = fDigitModuleLabel.substr( pos+1 );
        fDigitModuleLabel = fDigitModuleLabel.substr( 0, pos );
    }
    
    if (fOutputHistograms)
    {
        // Access ART's TFileService, which will handle creating and writing
        // histograms and n-tuples for us.
        art::ServiceHandle<art::TFileService> tfs;
    
        fPedestalOffsetVec.resize(3);
        fFullRMSVec.resize(3);
        fTruncRMSVec.resize(3);
        fNumTruncBinsVec.resize(3);
        fPedByChanVec.resize(3);
        fTruncRMSByChanVec.resize(3);
        fNumROIsHistVec.resize(3);
        fROILenHistVec.resize(3);
    
        for(size_t planeIdx = 0; planeIdx < 3; planeIdx++)
        {
            fPedestalOffsetVec[planeIdx] = tfs->make<TH1F>(    Form("PedPlane_%02zu",planeIdx),            ";Pedestal Offset (ADC);", 100, -5., 5.);
            fFullRMSVec[planeIdx]        = tfs->make<TH1F>(    Form("RMSFPlane_%02zu",planeIdx),           "Full RMS;RMS (ADC);", 400, 0., 40.);
            fTruncRMSVec[planeIdx]       = tfs->make<TH1F>(    Form("RMSTPlane_%02zu",planeIdx),           "Truncated RMS;RMS (ADC);", 100, 0., 10.);
            fNumTruncBinsVec[planeIdx]   = tfs->make<TH1F>(    Form("NTruncBins_%02zu",planeIdx),          ";# bins",     640, 0., 6400.);
            fPedByChanVec[planeIdx]      = tfs->make<TProfile>(Form("PedByWirePlane_%02zu",planeIdx),      ";Wire#", fGeometry->Nwires(planeIdx), 0., fGeometry->Nwires(planeIdx), -5., 5.);
            fTruncRMSByChanVec[planeIdx] = tfs->make<TProfile>(Form("TruncRMSByWirePlane_%02zu",planeIdx), ";Wire#", fGeometry->Nwires(planeIdx), 0., fGeometry->Nwires(planeIdx),  0., 10.);
            fNumROIsHistVec[planeIdx]    = tfs->make<TH1F>(    Form("NROISplane_%02zu",planeIdx),          ";# ROIs;",   100, 0, 100);
            fROILenHistVec[planeIdx]     = tfs->make<TH1F>(    Form("ROISIZEplane_%02zu",planeIdx),        ";ROI size;", 500, 0, 500);
        }
    }
    
    return;
}

//-------------------------------------------------
void Decon1DROI::beginJob()
{
    fEventCount = 0;
} // beginJob

//////////////////////////////////////////////////////
void Decon1DROI::endJob()
{
}
  
//////////////////////////////////////////////////////
void Decon1DROI::produce(art::Event& evt)
{
    // make a collection of Wires
    std::unique_ptr<std::vector<recob::Wire>> wireCol(new std::vector<recob::Wire>);
    // ... and an association set
    std::unique_ptr<art::Assns<raw::RawDigit,recob::Wire>> wireDigitAssn(new art::Assns<raw::RawDigit,recob::Wire>);

    // Read in the digit List object(s). 
    art::Handle< std::vector<raw::RawDigit>> digitVecHandle;
    
    if(fSpillName.size()>0) evt.getByLabel(fDigitModuleLabel, fSpillName, digitVecHandle);
    else                    evt.getByLabel(fDigitModuleLabel, digitVecHandle);

    if (!digitVecHandle->size())
    {
        evt.put(std::move(wireCol), fSpillName);
        evt.put(std::move(wireDigitAssn), fSpillName);
        fEventCount++;
        
        return;
    }

    // Reserve the room for the output
    wireCol->reserve(digitVecHandle->size());

    // ... Launch multiple threads with TBB to do the deconvolution and find ROIs in parallel
    multiThreadDeconvolutionProcessing deconvolutionProcessing(*this, evt, digitVecHandle, *wireCol, *wireDigitAssn);

    tbb::parallel_for(tbb::blocked_range<size_t>(0, digitVecHandle->size()), deconvolutionProcessing);
    
    // Time to stroe everything
    if(wireCol->size() == 0)
      mf::LogWarning("Decon1DROI") << "No wires made for this event.";

    //Make Histogram of recob::wire objects from Signal() vector
    // get access to the TFile service
    if ( fSaveWireWF ){
        art::ServiceHandle<art::TFileService> tfs;
        for (size_t wireN = 0; wireN < wireCol->size(); wireN++){
            std::vector<float> sigTMP = wireCol->at(wireN).Signal();
            TH1D* fWire = tfs->make<TH1D>(Form("Noise_Evt%04zu_N%04zu",fEventCount,wireN), ";Noise (ADC);",
				      sigTMP.size(),-0.5,sigTMP.size()-0.5);
            for (size_t tick = 0; tick < sigTMP.size(); tick++){
                fWire->SetBinContent(tick+1, sigTMP.at(tick) );
            }
        }
    }
    
    evt.put(std::move(wireCol), fSpillName);
    evt.put(std::move(wireDigitAssn), fSpillName);

    fEventCount++;

    return;
} // produce
    
float Decon1DROI::getTruncatedRMS(const std::vector<float>& waveform) const
{
    // do rms calculation - the old fashioned way and over all adc values
    std::vector<float> locWaveform = waveform;
    
    // sort in ascending order so we can truncate the sume
    std::sort(locWaveform.begin(), locWaveform.end(),[](const auto& left, const auto& right){return std::fabs(left) < std::fabs(right);});
    
    float threshold = fTruncRMSThreshold;
    
    std::vector<float>::iterator threshItr = std::find_if(locWaveform.begin(),locWaveform.end(),[threshold](const auto& val){return std::fabs(val) > threshold;});
    
    int minNumBins = std::max(int(fTruncRMSMinFraction * locWaveform.size()),int(std::distance(locWaveform.begin(),threshItr)));
    
    // Get the truncated sum
    float truncRms = std::inner_product(locWaveform.begin(), locWaveform.begin() + minNumBins, locWaveform.begin(), 0.);
    
    truncRms = std::sqrt(std::max(0.,truncRms / double(minNumBins)));
    
    return truncRms;
}
    
float Decon1DROI::fixTheFreakingWaveform(const std::vector<float>& waveform, raw::ChannelID_t channel, std::vector<float>& fixedWaveform) const
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
    
    int minNumBins = std::max(int(fTruncRMSMinFraction * locWaveform.size()),int(std::distance(locWaveform.begin(),threshItr)));
    
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
    
    // Fill histograms
    if (fOutputHistograms)
    {
        std::vector<geo::WireID> wids;
        try
        {
           wids = fGeometry->ChannelToWire(channel);
        }
        catch(...)
        {
            std::cout << "Caught exception looking up channel" << std::endl;
            return localRMS;
        }
    
        // Recover plane and wire in the plane
        size_t plane = wids[0].Plane;
        size_t wire  = wids[0].Wire;
        
//        float fullRMS = std::inner_product(locWaveform.begin(), locWaveform.end(), locWaveform.begin(), 0.);
        
//        fullRMS = std::sqrt(std::max(float(0.),fullRMS / float(locWaveform.size())));
    
        fPedestalOffsetVec[plane]->Fill(newPedestal,1.);
//        fFullRMSVec[plane]->Fill(fullRMS, 1.);
        fTruncRMSVec[plane]->Fill(localRMS, 1.);
        fNumTruncBinsVec[plane]->Fill(minNumBins, 1.);
        fPedByChanVec[plane]->Fill(wire, newPedestal, 1.);
        fTruncRMSByChanVec[plane]->Fill(wire, localRMS, 1.);
    }
    
    return localRMS;
}

void  Decon1DROI::processChannel(size_t                                  idx,
                                 art::Event&                             event,
                                 art::Handle<std::vector<raw::RawDigit>> digitVecHandle, 
                                 std::vector<recob::Wire>&               wireColVec, 
                                 art::Assns<raw::RawDigit,recob::Wire>&  wireAssns) const
{
    // vector that will be moved into the Wire object
    recob::Wire::RegionsOfInterest_t deconVec;
    recob::Wire::RegionsOfInterest_t ROIVec;

    // get the reference to the current raw::RawDigit
    art::Ptr<raw::RawDigit> digitVec(digitVecHandle, idx);

    raw::ChannelID_t channel = digitVec->Channel();
    
    // The following test is meant to be temporary until the "correct" solution is implemented
    if (!fChannelFilter->IsPresent(channel)) return;

    // Testing an idea about rejecting channels
    if (digitVec->GetPedestal() < 0.) return;

    float pedestal = 0.;
    
    // skip bad channels
    if( fChannelFilter->Status(channel) >= fMinAllowedChanStatus)
    {
        size_t dataSize = digitVec->Samples();
        
        // Recover the plane info
        std::vector<geo::WireID> wids; //    = fGeometry->ChannelToWire(channel);
        try
        {
            wids = fGeometry->ChannelToWire(channel);
        }
        catch(...)
        {
            std::cout << "Not able to find channel: " << channel << std::endl;
            return;
        }
        
        const geo::PlaneID&      planeID = wids[0].planeID();

        // vector holding uncompressed adc values
        std::vector<short> rawadc(dataSize);
        
        // uncompress the data
        raw::Uncompress(digitVec->ADCs(), rawadc, digitVec->Compression());
        
        // loop over all adc values and subtract the pedestal
        // When we have a pedestal database, can provide the digit timestamp as the third argument of GetPedestalMean
        try
        {
            pedestal = fPedRetrievalAlg->PedMean(channel);
        }
        catch(...)
        {
            mf::LogDebug("Decon1DROI_module") << "Pedestal lookup fails with channel: " << channel << std::endl;
            return;
        }
        
        
        // Get the pedestal subtracted data, centered in the deconvolution vector
        std::vector<float> rawAdcLessPedVec(dataSize);
        
        std::transform(rawadc.begin(),rawadc.end(),rawAdcLessPedVec.begin(),std::bind(std::minus<short>(),std::placeholders::_1,pedestal));
        
        // It seems there are deviations from the pedestal when using wirecell for noise filtering
        float raw_noise = fixTheFreakingWaveform(rawAdcLessPedVec, channel, rawAdcLessPedVec);
        
        // Recover a measure of the noise on the channel for use in the ROI finder
        //float raw_noise = getTruncatedRMS(rawAdcLessPedVec);
        
        // Try smoothing the input waveform
//        std::vector<float> rawAdcSmoothVec;
//        fWaveformTool->medianSmooth(rawAdcLessPedVec,rawAdcSmoothVec);
        
        // Make a dummy candidate roi vec
        icarus_tool::IROIFinder::CandidateROIVec deconROIVec;
        
        deconROIVec.push_back(icarus_tool::IROIFinder::CandidateROI(0,rawAdcLessPedVec.size()));
        
        // Do the deconvolution on the full waveform
        fDeconvolution->Deconvolve(rawAdcLessPedVec, channel, deconROIVec, deconVec);
        
        // Recover the deconvolved waveform
        const std::vector<float>& deconvolvedWaveform = deconVec.get_ranges().front().data();

        // vector of candidate ROI begin and end bins
        icarus_tool::IROIFinder::CandidateROIVec candRoiVec;
        
        // Now find the candidate ROI's
        fROIFinderVec.at(planeID.Plane)->FindROIs(deconvolvedWaveform, channel, fEventCount, raw_noise, candRoiVec);
        
        icarusutil::TimeVec holder;
        
        // We need to copy the deconvolved (and corrected) waveform ROI's
        for(const auto& candROI : candRoiVec)
        {
            // First up: copy out the relevent ADC bins into the ROI holder
            size_t roiLen = candROI.second - candROI.first;
            
            holder.resize(roiLen);
            
            std::copy(deconvolvedWaveform.begin()+candROI.first, deconvolvedWaveform.begin()+candROI.second, holder.begin());
            
            // Now we do the baseline determination and correct the ROI
            //float base = fBaseline->GetBaseline(holder, channel, roiStart, roiLen);
            float base = fBaseline->GetBaseline(holder, channel, 0, roiLen);
            
            std::transform(holder.begin(),holder.end(),holder.begin(),[base](auto& adcVal){return adcVal - base;});

            // add the range into ROIVec
            ROIVec.add_range(candROI.first, std::move(holder));
        }

        // Make some histograms?
        if (fOutputHistograms)
        {
            fNumROIsHistVec.at(planeID.Plane)->Fill(candRoiVec.size(), 1.);
            
            for(const auto& pair : candRoiVec)
                fROILenHistVec.at(planeID.Plane)->Fill(pair.second-pair.first, 1.);
        
            float fullRMS = std::inner_product(deconvolvedWaveform.begin(), deconvolvedWaveform.end(), deconvolvedWaveform.begin(), 0.);
        
            fullRMS = std::sqrt(std::max(float(0.),fullRMS / float(deconvolvedWaveform.size())));
    
            fFullRMSVec[planeID.Plane]->Fill(fullRMS, 1.);
        }
    } // end if not a bad channel
        
    // Don't save empty wires
    if (ROIVec.empty()) return;

    // First get a lock to make sure we are clear to run
    tbb::spin_mutex::scoped_lock lock(deconvolutionSpinMutex);

    // create the new wire directly in wirecol
    wireColVec.push_back(recob::WireCreator(std::move(ROIVec),*digitVec).move());

    // add an association between the last object in wirecol
    // (that we just inserted) and digitVec
    if (!util::CreateAssn(*this, event, wireColVec, digitVec, wireAssns, fSpillName))
    {
        throw art::Exception(art::errors::ProductRegistrationFailure)
            << "Can't associate wire #" << (wireColVec.size() - 1)
            << " with raw digit #" << digitVec.key();
    } // if failed to add association

    return;
}

} // end namespace caldata
