////////////////////////////////////////////////////////////////////////
//
// RecoWireROI class - variant of RecoWire that deconvolves in 
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

// ROOT libraries
#include "TComplex.h"
#include "TH1D.h"

// framework libraries
#include "fhiclcpp/ParameterSet.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h" 
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art_root_io/TFileService.h" 
#include "canvas/Utilities/Exception.h"

#include "art/Framework/Services/Registry/ServiceMacros.h"

// LArSoft libraries
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardata/ArtDataHelper/WireCreator.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "icaruscode/TPC/Utilities/SignalShapingICARUSService_service.h"
#include "icarus_signal_processing/ICARUSFFT.h"

#include "larevt/CalibrationDBI/Interface/DetPedestalService.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalProvider.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"


///creation of calibrated signals on wires
namespace recowire {

class RecoWireROI : public art::EDProducer
{
public:
    // create calibrated signals on wires. this class runs 
    // an fft to remove the electronics shaping.     
    explicit RecoWireROI(fhicl::ParameterSet const& pset); 
    virtual ~RecoWireROI();
    
    void produce(art::Event& evt); 
    void beginJob(); 
    void endJob();                 
    void reconfigure(fhicl::ParameterSet const& p);
    
private:
    
    std::string                 fDigitModuleLabel;     ///< module that made digits
    std::string                 fSpillName;            ///< nominal spill is an empty string
                                                       ///< it is set by the DigitModuleLabel
                                                       ///< ex.:  "daq:preSpill" for prespill data
    unsigned short              fNoiseSource;          ///< Used to determine ROI threshold
    unsigned short              fNumBinsHalf;          ///< Determines # bins in ROI running sum
    std::vector<unsigned short> fThreshold;            ///< abs(threshold) ADC counts for ROI
    std::vector<int>            fNumSigma;             ///< "# sigma" rms noise for ROI threshold
    std::vector<unsigned short> fPreROIPad;            ///< ROI padding
    std::vector<unsigned short> fPostROIPad;           ///< ROI padding
    bool                        fDoBaselineSub;        ///< Do baseline subtraction after deconvolution?
    bool                        fuPlaneRamp;           ///< set true for correct U plane wire response
    int                         fSaveWireWF;           ///< Save recob::wire object waveforms
    size_t                      fEventCount;           ///< count of event processed
    int                         fMinAllowedChanStatus; ///< Don't consider channels with lower status
    bool                        fDodQdxCalib;          ///< Do we apply wire-by-wire calibration?
    std::string                 fdQdxCalibFileName;    ///< Text file for constants to do wire-by-wire calibration
    std::map<unsigned int, float> fdQdxCalib;          ///< Map to do wire-by-wire calibration, key is channel number, content is correction factor

    float fMinROIAverageTickThreshold; // try to remove bad ROIs

    void doDecon(icarusutil::TimeVec&                          holder,
                 raw::ChannelID_t                              channel,
                 unsigned int                                  thePlane,
                 const std::vector<std::pair<size_t, size_t>>& rois,
                 const std::vector<std::pair<size_t, size_t>>& holderInfo,
                 recob::Wire::RegionsOfInterest_t&             ROIVec);
    
    float SubtractBaseline(std::vector<float>& holder,
                           float               basePre,
                           float               basePost,
                           size_t              roiStart,
                           size_t              roiLen,
                           size_t              dataSize);

    float SubtractBaseline(const std::vector<float>& holder);
    
    const geo::GeometryCore&                        fGeometry;
    icarusutil::SignalShapingICARUSService&         fSignalServices;
    const lariov::ChannelStatusProvider&            fChanFilt;
    std::unique_ptr<icarus_signal_processing::ICARUSFFT<double>>  fFFT;                  ///< Object to handle thread safe FFT
}; // class RecoWireROI

DEFINE_ART_MODULE(RecoWireROI)
  
//-------------------------------------------------
RecoWireROI::RecoWireROI(fhicl::ParameterSet const& pset) : EDProducer{pset},
    fGeometry(*lar::providerFrom<geo::Geometry>()),
    fSignalServices(*art::ServiceHandle<icarusutil::SignalShapingICARUSService>()),
    fChanFilt(art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider())
{
    this->reconfigure(pset);

    produces< std::vector<recob::Wire>  >(fSpillName);
    produces<art::Assns<raw::RawDigit, recob::Wire>>(fSpillName);
}

//-------------------------------------------------
RecoWireROI::~RecoWireROI()
{
}

//////////////////////////////////////////////////////
void RecoWireROI::reconfigure(fhicl::ParameterSet const& p)
{
    std::vector<unsigned short> uin;    std::vector<unsigned short> vin;
    std::vector<unsigned short> zin;
    
    fDigitModuleLabel     = p.get< std::string >                   ("DigitModuleLabel", "daq");
    fNoiseSource          = p.get< unsigned short >                ("NoiseSource",          3);
    fNumBinsHalf          = p.get< unsigned short >                ("NumBinsHalf",          3);
    fThreshold            = p.get< std::vector<unsigned short> >   ("Threshold"              );
    fNumSigma             = p.get< std::vector<int> >              ("NumSigma"               );
    uin                   = p.get< std::vector<unsigned short> >   ("uPlaneROIPad"           );
    vin                   = p.get< std::vector<unsigned short> >   ("vPlaneROIPad"           );
    zin                   = p.get< std::vector<unsigned short> >   ("zPlaneROIPad"           );
    fDoBaselineSub        = p.get< bool >                          ("DoBaselineSub"          );
    fuPlaneRamp           = p.get< bool >                          ("uPlaneRamp"             );
    fSaveWireWF           = p.get< int >                           ("SaveWireWF"             );
    fMinAllowedChanStatus = p.get< int >                           ("MinAllowedChannelStatus");
    fMinROIAverageTickThreshold = p.get<float>("MinROIAverageTickThreshold",-0.5);
    
    if(uin.size() != 2 || vin.size() != 2 || zin.size() != 2) {
      throw art::Exception(art::errors::Configuration)
        << "u/v/z plane ROI pad size != 2";
    }

    fPreROIPad.resize(3);
    fPostROIPad.resize(3);
    
    // put the ROI pad sizes into more convenient vectors
    fPreROIPad[0]  = uin[0];
    fPostROIPad[0] = uin[1];
    fPreROIPad[1]  = vin[0];
    fPostROIPad[1] = vin[1];
    fPreROIPad[2]  = zin[0];
    fPostROIPad[2] = zin[1];
    
    fSpillName.clear();
    
    size_t pos = fDigitModuleLabel.find(":");
    if( pos!=std::string::npos ) {
      fSpillName = fDigitModuleLabel.substr( pos+1 );
      fDigitModuleLabel = fDigitModuleLabel.substr( 0, pos );
    }

    //wire-by-wire calibration
    fDodQdxCalib        = p.get< bool >                          ("DodQdxCalib", false);
    if (fDodQdxCalib){
      fdQdxCalibFileName = p.get< std::string >                  ("dQdxCalibFileName");
      std::string fullname;
      cet::search_path sp("FW_SEARCH_PATH");
      sp.find_file(fdQdxCalibFileName, fullname);

      if (fullname.empty()) {
	std::cout << "Input file " << fdQdxCalibFileName << " not found" << std::endl;
	throw cet::exception("File not found");
      }
      else
	std::cout << "Applying wire-by-wire calibration using file " << fdQdxCalibFileName << std::endl;

      std::ifstream inFile(fullname, std::ios::in);
      std::string line;
      
      while (std::getline(inFile,line)) {
	unsigned int channel;
	float        constant;
	std::stringstream linestream(line);
	linestream >> channel >> constant;
	fdQdxCalib[channel] = constant;
	if (channel%1000==0) std::cout<<"Channel "<<channel<<" correction factor "<<fdQdxCalib[channel]<<std::endl;
      }
    }

    auto const* detprop      = lar::providerFrom<detinfo::DetectorPropertiesService>();
	
    // Now set up our plans for doing the convolution
    fFFT = std::make_unique<icarus_signal_processing::ICARUSFFT<double>>(detprop->NumberTimeSamples());
}

//-------------------------------------------------
void RecoWireROI::beginJob()
{
    fEventCount = 0;
} // beginJob

//////////////////////////////////////////////////////
void RecoWireROI::endJob()
{
}
  
//////////////////////////////////////////////////////
void RecoWireROI::produce(art::Event& evt)
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
    else evt.getByLabel(fDigitModuleLabel, digitVecHandle);

    if (!digitVecHandle->size())
    {
        evt.put(std::move(wirecol), fSpillName);
        evt.put(std::move(WireDigitAssn), fSpillName);
        return;
    }
    
    raw::ChannelID_t channel = raw::InvalidChannelID; // channel number

    auto const* detprop      = lar::providerFrom<detinfo::DetectorPropertiesService>();
    double      samplingRate = detprop->SamplingRate()/1000.;
    double      deconNorm    = fSignalServices.GetDeconNorm();

    // We'll need to set the transform size once we get the waveform and know its size
    size_t transformSize = 0;
    
    // loop over all wires
    wirecol->reserve(digitVecHandle->size());
    for(size_t rdIter = 0; rdIter < digitVecHandle->size(); ++rdIter)
    {
        // vector that will be moved into the Wire object
        recob::Wire::RegionsOfInterest_t ROIVec;
      
        // vector of ROI begin and end bins
        std::vector<std::pair<size_t, size_t>> rois;
      
        // get the reference to the current raw::RawDigit
        art::Ptr<raw::RawDigit> digitVec(digitVecHandle, rdIter);
        channel = digitVec->Channel();

        // The following test is meant to be temporary until the "correct" solution is implemented
        if (!fChanFilt.IsPresent(channel)) continue;

        // Testing an idea about rejecting channels
        if (digitVec->GetPedestal() < 0.) continue;
        
        // skip bad channels
        if( fChanFilt.Status(channel) >= fMinAllowedChanStatus)
        {
            size_t dataSize = digitVec->Samples();
            
            if (!dataSize)
            {
                std::cout << "***>> zero length data buffer, channel: " << channel << std::endl;
                continue;
            }
            
            // vector holding uncompressed adc values
            std::vector<short> rawadc(dataSize);
            
            std::vector<geo::WireID> wids     = fGeometry.ChannelToWire(channel);
            size_t                   thePlane = wids[0].Plane;
            
            // Set up the deconvolution and the vector to deconvolve
          
            // Set up the deconvolution and the vector to deconvolve
            // This is called only once per event, but under the hood nothing happens
            //   unless the FFT vector length changes (which it shouldn't for a run)
            if (!transformSize)
            {
                fSignalServices.SetDecon(dataSize, channel);
                transformSize = dataSize;
            }
            
            icarusutil::TimeVec rawAdcLessPedVec;
            
            rawAdcLessPedVec.resize(transformSize,0.);
            
            // uncompress the data
            raw::Uncompress(digitVec->ADCs(), rawadc, digitVec->Compression());
            
            // loop over all adc values and subtract the pedestal
            // When we have a pedestal database, can provide the digit timestamp as the third argument of GetPedestalMean
            size_t numBins   = 2 * fNumBinsHalf + 1;
            float  pdstl     = pedestalRetrievalAlg.PedMean(channel);
            float  rms_noise = digitVec->GetSigma();
            float  raw_noise = fSignalServices.GetRawNoise(channel);
            
            if      (fNoiseSource == 1) raw_noise = rms_noise;
            else if (fNoiseSource != 2) raw_noise = std::max(raw_noise,rms_noise);
            
            size_t binOffset(0);
            
            if (transformSize > dataSize) binOffset = (transformSize - dataSize) / 2;
            
            size_t startBin(binOffset);
            size_t stopBin(binOffset+numBins);
            
            //float startThreshold = sqrt(float(numBins)) * (fNumSigma[thePlane] * raw_noise + fThreshold[thePlane]);
            float startThreshold = float(numBins) * (fNumSigma[thePlane] * raw_noise + fThreshold[thePlane]);
            float stopThreshold  = startThreshold;
            
            // Get the pedestal subtracted data, centered in the deconvolution vector
            std::transform(rawadc.begin(),rawadc.end(),rawAdcLessPedVec.begin()+startBin,[pdstl](const short& adc){return std::round(float(adc) - pdstl);});
            std::fill(rawAdcLessPedVec.begin(),rawAdcLessPedVec.begin()+startBin,0.);        //rawAdcLessPedVec.at(startBin));
            std::fill(rawAdcLessPedVec.begin()+startBin+dataSize,rawAdcLessPedVec.end(),0.); //rawAdcLessPedVec.at(startBin+dataSize-1));
            
            // Try a loose cut to see if there is a potential for activity on this channel
            icarusutil::TimeVec::iterator overThreshItr = std::find_if(rawAdcLessPedVec.begin(),rawAdcLessPedVec.end(),[raw_noise](const auto& val){return val > 2.5 * raw_noise;});
            
            if (overThreshItr == rawAdcLessPedVec.end()) continue;

            float runningSum = std::accumulate(rawAdcLessPedVec.begin()+startBin,rawAdcLessPedVec.begin()+stopBin, 0.);
            
            size_t roiStartBin(0);
            bool   roiCandStart(false);

            // search for ROIs - follow prescription from Bruce B using a running sum to make faster
            // Note that we start in the middle of the running sum... if we find an ROI padding will extend
            // past this to take care of ends of the waveform
            for(size_t bin = fNumBinsHalf + 1; bin < dataSize - fNumBinsHalf; bin++)
            {
                // handle the running sum
                // Case, we are at start of waveform
                runningSum -= rawAdcLessPedVec[startBin++];
                
                // Case, we are at end of waveform
                runningSum += rawAdcLessPedVec[stopBin++];
                
                // We have already started a candidate ROI
                if (roiCandStart)
                {
                    if (fabs(runningSum) < stopThreshold)
                    {
                        if (bin - roiStartBin > 2) rois.push_back(std::make_pair(roiStartBin, bin));
                        
                        roiCandStart = false;
                    }
                }
                // Not yet started a candidate ROI
                else
                {
                    if (fabs(runningSum) > startThreshold)
                    {
                        roiStartBin  = bin;
                        roiCandStart = true;
                    }
                }
            } // bin
            // add the last ROI if existed
            if (roiCandStart)
            {
                //unsigned int roiLen = dataSize -1 - roiStart;
                // if(roiLen > fMinWid)
                rois.push_back(std::make_pair(roiStartBin, dataSize-1));
            }

            // skip deconvolution if there are no ROIs
            if(rois.size() == 0) continue;

            // pad the ROIs
            for(size_t ii = 0; ii < rois.size(); ++ii)
            {
                // low ROI end
                rois[ii].first  = std::max(int(rois[ii].first - fPreROIPad[thePlane]),0);
                // high ROI end
                rois[ii].second = std::min(rois[ii].second + fPostROIPad[thePlane],dataSize - 1);
            }

            // merge the ROIs?
            if(rois.size() > 1)
            {
                // temporary vector for merged ROIs
                std::vector<std::pair<size_t,size_t>> trois;
          
                // Loop through candidate roi's
                size_t startRoi = rois.front().first;
                size_t stopRoi  = rois.front().second;
                
                for(size_t idx = 1; idx < rois.size(); idx++)
                {
                    if (rois[idx].first < stopRoi)
                    {
                        stopRoi = rois[idx].second;
                    }
                    else
                    {
                        trois.push_back(std::pair<size_t,size_t>(startRoi,stopRoi));
                        
                        startRoi = rois[idx].first;
                        stopRoi  = rois[idx].second;
                    }
                }
                
                // Make sure to get the last one
                trois.push_back(std::pair<size_t,size_t>(startRoi,stopRoi));
	  
                rois = trois;
            }
            
            // Strategy is to run deconvolution on the entire channel and then pick out the ROI's we found above
            // Deconvolute the raw signal using the channel's nominal response
            fFFT->deconvolute(rawAdcLessPedVec, fSignalServices.GetResponse(channel).getDeconvKernel(), fSignalServices.FieldResponseTOffset(channel));
            
            std::vector<float> holder;
            
            for(size_t roiIdx = 0; roiIdx < rois.size(); roiIdx++)
            {
                const auto roi = rois[roiIdx];
                
                // First up: copy out the relevent ADC bins into the ROI holder
                size_t roiLen = roi.second - roi.first;
                
                holder.resize(roiLen);
                
                std::copy(rawAdcLessPedVec.begin()+binOffset+roi.first, rawAdcLessPedVec.begin()+binOffset+roi.second, holder.begin());
                float dnorm=samplingRate*deconNorm;
               // std::cout << " dnorm " << dnorm << std::endl;
                std::transform(holder.begin(),holder.end(),holder.begin(),[dnorm](float& deconVal){return deconVal/dnorm;});

                // Now we do the baseline determination (and I'm left wondering if there is a better way using the entire waveform?)
                bool  baseSet(false);
                float base(0.);
                if(fDoBaselineSub && fPreROIPad[thePlane] > 0 )
                {
                    //1. Check Baseline match?
                    // If not, include next ROI(if none, go to the end of signal)
                    // If yes, proceed
                    size_t binsToAve(20);
                    float  basePre  = std::accumulate(holder.begin(),holder.begin()+binsToAve,0.) / float(binsToAve);
                    float  basePost = std::accumulate(holder.end()-binsToAve,holder.end(),0.) / float(binsToAve);
                    
                    // emulate method for refining baseline from original version of RecoWireROI
                    float deconNoise = 1.26491 * fSignalServices.GetDeconNoise(channel);    // 4./sqrt(10) * noise

                    // If the estimated baseline from the front of the roi does not agree well with that from the end
                    // of the roi then we'll extend the roi hoping for good agreement
                    if (!(fabs(basePre - basePost) < deconNoise))
                    {
                        int   nTries(0);

                        // get start of roi and find the maximum we can extend to
                        icarusutil::TimeVec::iterator rawAdcRoiStartItr = rawAdcLessPedVec.begin() + binOffset + roi.first;
                        icarusutil::TimeVec::iterator rawAdcMaxItr      = rawAdcLessPedVec.end()   - binOffset;

                        // if this is not the last roi then limit max range to start of next roi
                        if (roiIdx < rois.size() - 1)
                            rawAdcMaxItr = rawAdcLessPedVec.begin() + binOffset + rois[roiIdx+1].first;

                        // Basically, allow roi to be extended until we get good agreement unless it seems pointless
                        while (!(fabs(basePre - basePost) < deconNoise) && nTries++ < 3)
                        {
                            size_t nBinsToAdd(100);
                            int    nBinsLeft      = std::distance(rawAdcRoiStartItr+roiLen,rawAdcMaxItr) > 0
                                                  ? std::distance(rawAdcRoiStartItr+roiLen,rawAdcMaxItr) : 0;
                            size_t roiLenAddition = std::min(nBinsToAdd, size_t(nBinsLeft));
                        
                            if (roiLenAddition > 0)
                            {
                                std::vector<float> additionVec(roiLenAddition);
                        
                                std::copy(rawAdcRoiStartItr + roiLen, rawAdcRoiStartItr + roiLen + roiLenAddition, additionVec.begin());
                            
                                holder.resize(holder.size() + roiLenAddition);
                                std::transform(additionVec.begin(),additionVec.end(),holder.begin() + roiLen,[deconNorm](float& deconVal){return deconVal/deconNorm;});
                            
                                basePost = std::accumulate(holder.end()-binsToAve,holder.end(),0.) / float(binsToAve);
                        
                                roiLen = holder.size();
                            }
                            else break;
                        }
                    }
                   
                    baseSet = true;
                    base    = SubtractBaseline(holder, basePre,basePost,roi.first,roiLen,dataSize);
                } // fDoBaselineSub ...
                
                if (baseSet) std::transform(holder.begin(),holder.end(),holder.begin(),[base](float& adcVal){return adcVal - base;});

                // apply wire-by-wire calibration
                if (fDodQdxCalib)
                {
                    if(fdQdxCalib.find(channel) != fdQdxCalib.end())
                    {
                        float constant = fdQdxCalib[channel];
                        //std::cout<<channel<<" "<<constant<<std::endl;
                        for (size_t iholder = 0; iholder < holder.size(); ++iholder) holder[iholder] *= constant;
                    }
                }

                //wes 23.12.2016 --- sum up the roi, and if it's very negative get rid of it
                float average_val = std::accumulate(holder.begin(),holder.end(),0.0) / holder.size();
                float min = *std::min_element(holder.begin(),holder.end());
                float max = *std::max_element(holder.begin(),holder.end());
                if(average_val>fMinROIAverageTickThreshold || std::abs(min)<std::abs(max)){
                    // add the range into ROIVec
                    ROIVec.add_range(roi.first, std::move(holder));
                }
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
      mf::LogWarning("RecoWireROI") << "No wires made for this event.";

    //Make Histogram of recob::wire objects from Signal() vector
    // get access to the TFile service
    if ( fSaveWireWF )
    {
        art::ServiceHandle<art::TFileService> tfs;
        for (size_t wireN = 0; wireN < wirecol->size(); wireN++)
        {
            std::vector<float> sigTMP = wirecol->at(wireN).Signal();
            TH1D* fWire = tfs->make<TH1D>(Form("Noise_Evt%04zu_N%04zu",fEventCount,wireN), ";Noise (ADC);", sigTMP.size(),-0.5,sigTMP.size()-0.5);
            for (size_t tick = 0; tick < sigTMP.size(); tick++) fWire->SetBinContent(tick+1, sigTMP.at(tick) );
        }
    }
    
    evt.put(std::move(wirecol),             fSpillName);
    evt.put(std::move(WireDigitAssn),       fSpillName);

    fEventCount++;

    return;
} // produce


float RecoWireROI::SubtractBaseline(std::vector<float>& holder,
                                   float               basePre,
                                   float               basePost,
                                   size_t              roiStart,
                                   size_t              roiLen,
                                   size_t              dataSize)
{
    float base=0;

    //can not trust the early part
    if (roiStart < 20 && roiStart + roiLen < dataSize - 20){
        base = basePost;
    // can not trust the later part
    }else if (roiStart >= 20 && roiStart + roiLen >= dataSize - 20){
        base = basePre;
    // can trust both
    }else if (roiStart >= 20 && roiStart + roiLen < dataSize - 20){
        if (fabs(basePre-basePost)<3){
            base = (basePre+basePost)/2.;
        }else{
            if (basePre < basePost){
                base = basePre;
            }else{
                base = basePost;
            }
        }
    // can not use both
    }else{
        float min = 0,max=0;
        for (unsigned int bin = 0; bin < roiLen; bin++){
            if (holder[bin] > max) max = holder[bin];
            if (holder[bin] < min) min = holder[bin];
        }
        int nbin = max - min;
        if (nbin!=0){
            TH1F *h1 = new TH1F("h1","h1",nbin,min,max);
            for (unsigned int bin = 0; bin < roiLen; bin++){
                h1->Fill(holder[bin]);
            }
            float ped = h1->GetMaximum();
            float ave=0,ncount = 0;
            for (unsigned int bin = 0; bin < roiLen; bin++){
                if (fabs(holder[bin]-ped)<2){
                    ave +=holder[bin];
                    ncount ++;
                }
            }
            if (ncount==0) ncount=1;
            ave = ave/ncount;
            h1->Delete();
            base = ave;
        }
    }
    
    return base;
}

  
void RecoWireROI::doDecon(icarusutil::TimeVec&                         holder,
                          raw::ChannelID_t                             channel,
                          unsigned int                                 thePlane,
                          const std::vector<std::pair<size_t,size_t>>& rois,
                          const std::vector<std::pair<size_t,size_t>>& holderInfo,
                          recob::Wire::RegionsOfInterest_t&            ROIVec)
{
    // Deconvolute the raw signal using the channel's nominal response
    fFFT->deconvolute(holder, fSignalServices.GetResponse(channel).getDeconvKernel(), fSignalServices.FieldResponseTOffset(channel));

    // transfer the ROIs and start bins into the vector that will be
    // put into the event
    for(size_t jr = 0; jr < holderInfo.size(); ++jr)
    {
        std::vector<float> sigTemp;
        size_t bBegin = holderInfo[jr].first;
        size_t theROI = holderInfo[jr].second;
        size_t roiLen = rois[theROI].second - rois[theROI].first;
        size_t bEnd = bBegin + roiLen;
        float basePre = 0., basePost = 0.;
        // Baseline subtraction if requested and the ROIs are padded.
        // Can't baseline subtract signals when the ROI start is close to 0 either
        if(fDoBaselineSub && fPreROIPad[thePlane] > 0 &&
           rois[theROI].first > fPreROIPad[thePlane])
        {
            // find the baseline from the first few bins in the leading Pad region
            unsigned short bbins = fPreROIPad[thePlane];
            size_t bin;
            if(bbins > 5) bbins = 5;
            for(bin = 0; bin < bbins; ++bin) {
                basePre  += holder[bBegin + bin];
                basePost += holder[bEnd - bin];
            }
            basePre /= (float)bbins;
            basePost /= (float)bbins;
            float slp = (basePost - basePre) / (float)(roiLen - bbins);
            float base;
            for(size_t jj = bBegin; jj < bEnd; ++jj) {
                base = basePre + slp * (jj - bBegin);
                sigTemp.push_back(holder[jj] - base);
            } // jj
        } // fDoBaselineSub ...
        else {
            for(size_t jj = bBegin; jj < bEnd; ++jj) sigTemp.push_back(holder[jj]);
        } // !fDoBaselineSub ...
      
        // add the range into ROIVec
        ROIVec.add_range(rois[theROI].first, std::move(sigTemp));
    } // jr
} // doDecon


} // end namespace recowire
