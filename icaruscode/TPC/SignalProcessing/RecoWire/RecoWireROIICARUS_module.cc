////////////////////////////////////////////////////////////////////////
//
// RecoWireROIICARUS class - variant of CalWire that deconvolves in
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
#include "icaruscode/TPC/Utilities/tools/IWaveformTool.h"

///creation of calibrated signals on wires
namespace caldata {

class RecoWireROIICARUS : public art::EDProducer
{
  public:
    // create calibrated signals on wires. this class runs 
    // an fft to remove the electronics shaping.     
    explicit RecoWireROIICARUS(fhicl::ParameterSet const& pset);
    virtual ~RecoWireROIICARUS();
    
    void produce(art::Event& evt); 
    void beginJob(); 
    void endJob();                 
    void reconfigure(fhicl::ParameterSet const& p);
    
  private:
    // It seems there are pedestal shifts that need correcting
    float fixTheFreakingWaveform(const std::vector<float>&, raw::ChannelID_t, std::vector<float>&) const;
    
    float getTruncatedRMS(const std::vector<float>&) const;
    
    std::string                                             fDigitModuleLabel;           ///< module that made digits
    std::string                                             fSpillName;                  ///< nominal spill is an empty string
                                                                                         ///< it is set by the DigitModuleLabel
                                                                                         ///< ex.:  "daq:preSpill" for prespill data
    unsigned short                                          fNoiseSource;                ///< Used to determine ROI threshold
    int                                                     fSaveWireWF;                 ///< Save recob::wire object waveforms
    size_t                                                  fEventCount;                 ///< count of event processed
    int                                                     fMinAllowedChanStatus;       ///< Don't consider channels with lower status
    
    float                                                   fTruncRMSThreshold;          ///< Calculate RMS up to this threshold...
    float                                                   fTruncRMSMinFraction;        ///< or at least this fraction of time bins
    bool                                                    fOutputHistograms;           ///< Output histograms?
    
    std::vector<std::unique_ptr<icarus_tool::IROIFinder>>   fROIFinderVec;               ///< ROI finders per plane
    std::unique_ptr<icarus_tool::IDeconvolution>            fDeconvolution;
    std::unique_ptr<icarus_tool::IWaveformTool>             fWaveformTool;
    
    const geo::GeometryCore*                                fGeometry = lar::providerFrom<geo::Geometry>();
    
    // Define here a temporary set of histograms...
    std::vector<TH1F*>     fPedestalOffsetVec;
    std::vector<TH1F*>     fFullRMSVec;
    std::vector<TH1F*>     fTruncRMSVec;
    std::vector<TH1F*>     fNumTruncBinsVec;
    std::vector<TProfile*> fPedByChanVec;
    std::vector<TProfile*> fTruncRMSByChanVec;
    std::vector<TH1F*>     fNumROIsHistVec;
    std::vector<TH1F*>     fROILenHistVec;
    
}; // class RecoWireROIICARUS

DEFINE_ART_MODULE(RecoWireROIICARUS)
  
//-------------------------------------------------
RecoWireROIICARUS::RecoWireROIICARUS(fhicl::ParameterSet const& pset) : EDProducer{pset}
{
  this->reconfigure(pset);

  produces< std::vector<recob::Wire> >(fSpillName);
  produces<art::Assns<raw::RawDigit, recob::Wire>>(fSpillName);
}

//-------------------------------------------------
RecoWireROIICARUS::~RecoWireROIICARUS()
{
}

//////////////////////////////////////////////////////
void RecoWireROIICARUS::reconfigure(fhicl::ParameterSet const& pset)
{
    // Recover the vector of fhicl parameters for the ROI tools
    const fhicl::ParameterSet& roiFinderTools = pset.get<fhicl::ParameterSet>("ROIFinderToolVec");
    
    fROIFinderVec.resize(roiFinderTools.get_pset_names().size());
    
    for(const std::string& roiFinderTool : roiFinderTools.get_pset_names())
    {
        const fhicl::ParameterSet& roiFinderToolParamSet = roiFinderTools.get<fhicl::ParameterSet>(roiFinderTool);
        size_t                     planeIdx              = roiFinderToolParamSet.get<size_t>("Plane");
        
        fROIFinderVec.at(planeIdx) = art::make_tool<icarus_tool::IROIFinder>(roiFinderToolParamSet);
    }
    
    std::sort(fROIFinderVec.begin(),fROIFinderVec.end(),[](const auto& left,const auto& right){return left->plane() < right->plane();});

    fDeconvolution = art::make_tool<icarus_tool::IDeconvolution>(pset.get<fhicl::ParameterSet>("Deconvolution"));
    
    // Let's apply some smoothing as an experiment... first let's get the tool we need
    fhicl::ParameterSet waveformToolParams;
    
    waveformToolParams.put<std::string>("tool_type","Waveform");
    
    fWaveformTool = art::make_tool<icarus_tool::IWaveformTool>(waveformToolParams);

    fDigitModuleLabel           = pset.get< std::string >   ("DigitModuleLabel", "daq");
    fNoiseSource                = pset.get< unsigned short >("NoiseSource",          3);
    fSaveWireWF                 = pset.get< int >           ("SaveWireWF"             );
    fMinAllowedChanStatus       = pset.get< int >           ("MinAllowedChannelStatus");
    fTruncRMSThreshold          = pset.get< float >         ("TruncRMSThreshold",    6.);
    fTruncRMSMinFraction        = pset.get< float >         ("TruncRMSMinFraction", 0.6);
    fOutputHistograms           = pset.get< bool  >         ("OutputHistograms",   true);
    
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
            fFullRMSVec[planeIdx]        = tfs->make<TH1F>(    Form("RMSFPlane_%02zu",planeIdx),           "Full RMS;RMS (ADC);", 100, 0., 10.);
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
void RecoWireROIICARUS::beginJob()
{
    fEventCount = 0;
} // beginJob

//////////////////////////////////////////////////////
void RecoWireROIICARUS::endJob()
{
}
  
//////////////////////////////////////////////////////
void RecoWireROIICARUS::produce(art::Event& evt)
{
    //get pedestal conditions
    const lariov::DetPedestalProvider& pedestalRetrievalAlg = art::ServiceHandle<lariov::DetPedestalService>()->GetPedestalProvider();
    
    // make a collection of Wires
    std::unique_ptr<std::vector<recob::Wire> > wirecol(new std::vector<recob::Wire>);
    // ... and an association set
    std::unique_ptr<art::Assns<raw::RawDigit,recob::Wire> > WireDigitAssn(new art::Assns<raw::RawDigit,recob::Wire>);

    // Read in the digit List object(s). 
    art::Handle< std::vector<raw::RawDigit> > digitVecHandle;
    
    if(fSpillName.size()>0) evt.getByLabel(fDigitModuleLabel, fSpillName, digitVecHandle);
    else                    evt.getByLabel(fDigitModuleLabel, digitVecHandle);

    if (!digitVecHandle->size())
    {
        evt.put(std::move(wirecol), fSpillName);
        evt.put(std::move(WireDigitAssn), fSpillName);
        return;
    }
    
    raw::ChannelID_t channel = raw::InvalidChannelID; // channel number
    
    const lariov::ChannelStatusProvider& chanFilt = art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider();
    
    // loop over all wires
    wirecol->reserve(digitVecHandle->size());
    for(size_t rdIter = 0; rdIter < digitVecHandle->size(); ++rdIter)
    {
        // vector that will be moved into the Wire object
        recob::Wire::RegionsOfInterest_t ROIVec;

        // get the reference to the current raw::RawDigit
        art::Ptr<raw::RawDigit> digitVec(digitVecHandle, rdIter);
        channel = digitVec->Channel();
      
        // The following test is meant to be temporary until the "correct" solution is implemented
        if (!chanFilt.IsPresent(channel)) continue;

        // Testing an idea about rejecting channels
        if (digitVec->GetPedestal() < 0.) continue;

        float pedestal = 0.;
        
        // skip bad channels
        if( chanFilt.Status(channel) >= fMinAllowedChanStatus)
        {
            size_t dataSize = digitVec->Samples();
            
            // Recover the plane info
            std::vector<geo::WireID> wids    = fGeometry->ChannelToWire(channel);
            const geo::PlaneID&      planeID = wids[0].planeID();

            // vector holding uncompressed adc values
            std::vector<short> rawadc(dataSize);
            
            // uncompress the data
            raw::Uncompress(digitVec->ADCs(), rawadc, digitVec->Compression());
            
            // loop over all adc values and subtract the pedestal
            // When we have a pedestal database, can provide the digit timestamp as the third argument of GetPedestalMean
            pedestal = pedestalRetrievalAlg.PedMean(channel);
            
            // Get the pedestal subtracted data, centered in the deconvolution vector
            std::vector<float> rawAdcLessPedVec(dataSize);
            
            std::transform(rawadc.begin(),rawadc.end(),rawAdcLessPedVec.begin(),std::bind(std::minus<short>(),std::placeholders::_1,pedestal));
            
            // It seems there are deviations from the pedestal when using wirecell for noise filtering
            float raw_noise = fixTheFreakingWaveform(rawAdcLessPedVec, channel, rawAdcLessPedVec);
            
            // Recover a measure of the noise on the channel for use in the ROI finder
            //float raw_noise = getTruncatedRMS(rawAdcLessPedVec);
            
            // Try smoothing the input waveform
//            std::vector<float> rawAdcSmoothVec;
//            fWaveformTool->medianSmooth(rawAdcLessPedVec,rawAdcSmoothVec);
            
            // vector of candidate ROI begin and end bins
            icarus_tool::IROIFinder::CandidateROIVec candRoiVec;

            // Now find the candidate ROI's
            fROIFinderVec.at(planeID.Plane)->FindROIs(rawAdcLessPedVec, channel, fEventCount, raw_noise, candRoiVec);
            
            // Do the deconvolution
            fDeconvolution->Deconvolve(rawAdcLessPedVec, channel, candRoiVec, ROIVec);
            
            // Make some histograms?
            if (fOutputHistograms)
            {
                // First up, determine what kind of wire we have
                std::vector<geo::WireID> wids    = fGeometry->ChannelToWire(channel);
                const geo::PlaneID&      planeID = wids[0].planeID();
                
                fNumROIsHistVec.at(planeID.Plane)->Fill(candRoiVec.size(), 1.);
                
                for(const auto& pair : candRoiVec)
                    fROILenHistVec.at(planeID.Plane)->Fill(pair.second-pair.first, 1.);
            }
        } // end if not a bad channel

        // create the new wire directly in wirecol
        wirecol->push_back(recob::WireCreator(std::move(ROIVec),*digitVec).move());

        // add an association between the last object in wirecol
        // (that we just inserted) and digitVec
        if (!util::CreateAssn(*this, evt, *wirecol, digitVec, *WireDigitAssn, fSpillName))
        {
            throw art::Exception(art::errors::ProductRegistrationFailure)
                << "Can't associate wire #" << (wirecol->size() - 1)
                << " with raw digit #" << digitVec.key();
        } // if failed to add association
        //  DumpWire(wirecol->back()); // for debugging
    }

    if(wirecol->size() == 0)
      mf::LogWarning("RecoWireROIICARUS") << "No wires made for this event.";

    //Make Histogram of recob::wire objects from Signal() vector
    // get access to the TFile service
    if ( fSaveWireWF ){
        art::ServiceHandle<art::TFileService> tfs;
        for (size_t wireN = 0; wireN < wirecol->size(); wireN++){
            std::vector<float> sigTMP = wirecol->at(wireN).Signal();
            TH1D* fWire = tfs->make<TH1D>(Form("Noise_Evt%04zu_N%04zu",fEventCount,wireN), ";Noise (ADC);",
				      sigTMP.size(),-0.5,sigTMP.size()-0.5);
            for (size_t tick = 0; tick < sigTMP.size(); tick++){
                fWire->SetBinContent(tick+1, sigTMP.at(tick) );
            }
        }
    }
    
    evt.put(std::move(wirecol), fSpillName);
    evt.put(std::move(WireDigitAssn), fSpillName);

    fEventCount++;

    return;
} // produce
    
float RecoWireROIICARUS::getTruncatedRMS(const std::vector<float>& waveform) const
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
    
float RecoWireROIICARUS::fixTheFreakingWaveform(const std::vector<float>& waveform, raw::ChannelID_t channel, std::vector<float>& fixedWaveform) const
{
    // Get the truncated mean and rms
    float fullRMS;
    float truncRMS;
    float truncMean;
    float nSig(2.0);  // make tight constraint
    int   nTrunc;
    
    fWaveformTool->getTruncatedMeanRMS(waveform, nSig, truncMean, fullRMS, truncRMS, nTrunc);
    
    // Set the waveform to the new baseline
    std::transform(waveform.begin(), waveform.end(), fixedWaveform.begin(), std::bind(std::minus<float>(),std::placeholders::_1,truncMean));
    
    // Fill histograms
    if (fOutputHistograms)
    {
        std::vector<geo::WireID> wids = fGeometry->ChannelToWire(channel);
    
        // Recover plane and wire in the plane
        size_t plane = wids[0].Plane;
        size_t wire  = wids[0].Wire;
    
        fPedestalOffsetVec[plane]->Fill(truncMean,1.);
        fFullRMSVec[plane]->Fill(fullRMS, 1.);
        fTruncRMSVec[plane]->Fill(truncRMS, 1.);
        fNumTruncBinsVec[plane]->Fill(nTrunc, 1.);
        fPedByChanVec[plane]->Fill(wire, truncMean, 1.);
        fTruncRMSByChanVec[plane]->Fill(wire, truncRMS, 1.);
    }
    
    return truncRMS;
}

} // end namespace caldata
