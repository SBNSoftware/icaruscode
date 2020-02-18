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
#include <iostream>
#include <iomanip>
#include <fstream>
#include <random>

// Intel Threading Building Blocks
#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"

// ROOT libraries
#include "TH1D.h"
#include "TProfile.h"

// framework libraries
#include "fhiclcpp/ParameterSet.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/ReplicatedProducer.h"
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
#include "lardata/Utilities/LArFFTWPlan.h"
#include "icaruscode/TPC/Utilities/SignalShaperServiceICARUS.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalService.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalProvider.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"

#include "icaruscode/TPC/SignalProcessing/RecoWire/DeconTools/IROIFinder.h"
#include "icaruscode/TPC/SignalProcessing/RecoWire/DeconTools/IBaseline.h"
#include "icaruscode/TPC/Utilities/tools/IWaveformTool.h"

using std::cout;
using std::endl;
using std::ofstream;

///creation of calibrated signals on wires
namespace caldata {

class mtDecon1DROI : public art::ReplicatedProducer
{
  public:
    // create calibrated signals on wires. this class runs 
    // an fft to remove the electronics shaping.     
    explicit mtDecon1DROI(fhicl::ParameterSet const& pset, art::ProcessingFrame const& frame); 
    virtual ~mtDecon1DROI();
    
    void produce(art::Event& evt, art::ProcessingFrame const& frame); 
    void beginJob(art::ProcessingFrame const& frame);
    void endJob(art::ProcessingFrame const& frame);
    void reconfigure(fhicl::ParameterSet const& p);
    void Deconvolute(size_t rdIter, 
                    art::Handle< std::vector<raw::RawDigit> >&digitVecHandle,std::vector<recob::Wire>& wirecol,
                    art::Assns<raw::RawDigit,recob::Wire>& WireDigitAssn, art::Event& evt, std::vector<int>& dgvindx,
                    void* fplan, void* rplan)const;
    
  private:
    // It seems there are pedestal shifts that need correcting
    float fixTheFreakingWaveform(const std::vector<float>&, raw::ChannelID_t, std::vector<float>&) const;
    
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
    
    bool                                                    fDodQdxCalib;                ///< Do we apply wire-by-wire calibration?
    std::string                                             fdQdxCalibFileName;          ///< Text file for constants to do wire-by-wire calibration
    std::map<unsigned int, float>                           fdQdxCalib;                  ///< Map to do wire-by-wire calibration, key is channel

    art::ServiceHandle<util::SignalShaperServiceICARUS>     fSignalShaper;
    const geo::GeometryCore*                                fGeometry = lar::providerFrom<geo::Geometry>();

    std::vector<std::unique_ptr<icarus_tool::IROIFinder>>   fROIFinderVec;               ///< ROI finders per plane
    std::unique_ptr<icarus_tool::IWaveformTool>             fWaveformTool;
    std::unique_ptr<icarus_tool::IBaseline>                 fBaseline;

    int                                                     fLogLevel;                   ///< Log level: 0=none, 1=init
    ofstream _fileout1;
    
}; // class mtDecon1DROI

DEFINE_ART_MODULE(mtDecon1DROI)

//////////////////////////////////////////////////////

class lartbb_Deconvolute {
  public:
    lartbb_Deconvolute(mtDecon1DROI const & prod,
      art::Handle< std::vector<raw::RawDigit> >&dgvHandle, std::vector<recob::Wire>& wcol,
      art::Assns<raw::RawDigit,recob::Wire>& wdassn, art::Event& evnt, std::vector<int>& dgvindx,
      void* fplan, void* rplan)
      : prod(prod),
        digitVecHandle(dgvHandle),
        wirecol(wcol),
        WireDigitAssn(wdassn),
        evt(evnt),
	dgvindx(dgvindx),
	fplan(fplan),
	rplan(rplan){}
    void operator()(const tbb::blocked_range<size_t>& range) const{
      //std::cout << " !!!!!!!!!! range.begin(): " << range.begin() << " and range.end(): " << range.end() << std::endl;
      for (size_t i = range.begin(); i < range.end(); ++i)
        prod.Deconvolute(i, digitVecHandle, wirecol, WireDigitAssn,  evt, dgvindx, fplan, rplan);
    }
  private:
    mtDecon1DROI const & prod;
    art::Handle< std::vector<raw::RawDigit> > &digitVecHandle;
    std::vector<recob::Wire> &wirecol;
    art::Assns<raw::RawDigit,recob::Wire> &WireDigitAssn;
    art::Event& evt;
    std::vector<int> &dgvindx;
    void* fplan;
    void* rplan;
};

//-------------------------------------------------
mtDecon1DROI::mtDecon1DROI(fhicl::ParameterSet const& pset, art::ProcessingFrame const& frame) : art::ReplicatedProducer(pset, frame)
{
  this->reconfigure(pset);

  produces< std::vector<recob::Wire> >(fSpillName);
  produces<art::Assns<raw::RawDigit, recob::Wire>>(fSpillName);
}

//-------------------------------------------------
mtDecon1DROI::~mtDecon1DROI()
{
}

//////////////////////////////////////////////////////
void mtDecon1DROI::reconfigure(fhicl::ParameterSet const& pset)
{
    // Get signal shaping service.
    fSignalShaper = art::ServiceHandle<util::SignalShaperServiceICARUS>();

    // Let's apply some smoothing as an experiment... first let's get the tool we need
    fhicl::ParameterSet waveformToolParams;
    waveformToolParams.put<std::string>("tool_type","Waveform");    
    fWaveformTool = art::make_tool<icarus_tool::IWaveformTool>(waveformToolParams);

    // Recover the parameters for wire-by-wire calibration
    fDodQdxCalib   = pset.get< bool >("DodQdxCalib", false);
    
    if (fDodQdxCalib)
    {
        fdQdxCalibFileName = pset.get< std::string >("dQdxCalibFileName");
        std::string fullname;
        cet::search_path sp("FW_SEARCH_PATH");
        sp.find_file(fdQdxCalibFileName, fullname);
        
        if (fullname.empty())
        {
            std::cout << "Input file " << fdQdxCalibFileName << " not found" << std::endl;
            throw cet::exception("File not found");
        }
        else
            std::cout << "Applying wire-by-wire calibration using file " << fdQdxCalibFileName << std::endl;
        
        std::ifstream inFile(fullname, std::ios::in);
        std::string line;
        
        while (std::getline(inFile,line))
        {
            unsigned int channel;
            float        constant;
            std::stringstream linestream(line);
            linestream >> channel >> constant;
            fdQdxCalib[channel] = constant;
            if (channel%1000==0) std::cout<<"Channel "<<channel<<" correction factor "<<fdQdxCalib[channel]<<std::endl;
        }
    }

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

    // Recover the baseline tool
    fBaseline  = art::make_tool<icarus_tool::IBaseline> (pset.get<fhicl::ParameterSet>("Baseline"));

    fDigitModuleLabel           = pset.get< std::string >   ("DigitModuleLabel", "daq");
    fNoiseSource                = pset.get< unsigned short >("NoiseSource",          3);
    fSaveWireWF                 = pset.get< int >           ("SaveWireWF"             );
    fMinAllowedChanStatus       = pset.get< int >           ("MinAllowedChannelStatus");
    fTruncRMSThreshold          = pset.get< float >         ("TruncRMSThreshold",    6.);
    fTruncRMSMinFraction        = pset.get< float >         ("TruncRMSMinFraction", 0.6);
    fOutputHistograms           = pset.get< bool  >         ("OutputHistograms",   true);
    
    fLogLevel = 0;
    pset.get_if_present<int>("LogLevel", fLogLevel);

    fSpillName.clear();
    
    size_t pos = fDigitModuleLabel.find(":");
    if( pos!=std::string::npos )
    {
        fSpillName = fDigitModuleLabel.substr( pos+1 );
        fDigitModuleLabel = fDigitModuleLabel.substr( 0, pos );
    }

    return;
}

//-------------------------------------------------
void mtDecon1DROI::beginJob(art::ProcessingFrame const&)
{
    fEventCount = 0;
    if ( fLogLevel >= 1 )_fileout1.open("dumpsignals.txt");
} // beginJob

//////////////////////////////////////////////////////
void mtDecon1DROI::endJob(art::ProcessingFrame const&)
{
    if ( fLogLevel >= 1 )_fileout1.close();
}
  
//////////////////////////////////////////////////////
void mtDecon1DROI::produce(art::Event& evt, art::ProcessingFrame const&)
{
    if ( fLogLevel >= 1 ){
      _fileout1 << "Event "
    	        << " " << evt.id().run()
    	        << " " << evt.id().subRun()
    	        << " " << evt.id().event() << endl;
    }
    
    // make a collection of Wires and an association set
    std::unique_ptr<std::vector<recob::Wire> > wirecol(new std::vector<recob::Wire>);
    std::unique_ptr<art::Assns<raw::RawDigit,recob::Wire> > WireDigitAssn(new art::Assns<raw::RawDigit,recob::Wire>);

    // Read in the digit List object(s). 
    art::Handle< std::vector<raw::RawDigit> > digitVecHandle;
    if(fSpillName.size()>0) evt.getByLabel(fDigitModuleLabel, fSpillName, digitVecHandle);
    else                    evt.getByLabel(fDigitModuleLabel, digitVecHandle);

    if (!digitVecHandle->size()){
        evt.put(std::move(wirecol), fSpillName);
        evt.put(std::move(WireDigitAssn), fSpillName);
        return;
    }

    // mwang added
    art::Ptr<raw::RawDigit> digitVec0(digitVecHandle, 0);
    unsigned int dataSize0 = digitVec0->Samples();
    int transformSize = fSignalShaper->GetFFTSize();
    if (int(dataSize0) != transformSize){
      throw art::Exception(art::errors::Configuration)
        << "transform size "<< transformSize << " does not match data size " << dataSize0;
    }
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // ---- Loop over all wires ----  
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    wirecol->resize(digitVecHandle->size());
    std::vector<int> dgvindx(digitVecHandle->size(),-1);

    // ... Create the fftw "plan" to execute on
    util::LArFFTWPlan fftp(transformSize,"ES");

    //for(size_t rdIter = 0; rdIter < digitVecHandle->size(); ++rdIter)
    //{
    //  Deconvolute(rdIter, digitVecHandle, *wirecol, *WireDigitAssn,  evt, dgvindx, fftp.fPlan, fftp.rPlan);
    //} // end of loop over wires

    // ... Launch multiple threads with TBB to do the deconvolution and find ROIs in parallel
    auto func = lartbb_Deconvolute(*this,digitVecHandle,*wirecol,*WireDigitAssn,evt,dgvindx,fftp.fPlan,fftp.rPlan);
    tbb::parallel_for(tbb::blocked_range<size_t>(0, digitVecHandle->size()), func);    

    // ... If wire collection is not empty, compress the wirecol and dgvindx vectors by removing all slots with no entries
    if(wirecol->size() == 0) {
      mf::LogWarning("mtDecon1DROI") << "No wires made for this event.";
    } else {
      wirecol->erase(std::remove_if(wirecol->begin(),wirecol->end(),[](const recob::Wire & w){return w.SignalROI().n_ranges()==0;}),
                     wirecol->end());
      dgvindx.erase(std::remove_if(dgvindx.begin(),dgvindx.end(),[](const int & dvidx){return dvidx==-1;}),dgvindx.end());
      if ( fLogLevel >= 1 ){
        int widx=0;
        for(std::vector<recob::Wire>::iterator it=wirecol->begin(); it != wirecol->end(); ++it){
          int isempty=1;
          if(it->SignalROI().n_ranges()>0)isempty=0;
          _fileout1 << "  Wire " << widx << ": channel= " << it->Channel() << ", view= " << it->View()
               << " has number of ranges= " << it->SignalROI().n_ranges() << " isEmpty=" << isempty << endl;
      
          int ridx=0;
          for(const auto& range : it->SignalROI().get_ranges()){
            const std::vector<float>& signal = range.data();
            _fileout1 << "  --> Range # " << ridx << " with " << signal.size() << " bins:" << endl;
            int bidx=0;
            for(std::vector<float>::const_iterator its=signal.begin(); its != signal.end(); ++its){
              _fileout1 << "      Bin " << bidx << ": " << *its << endl;
              bidx++;
            }
            ridx++;
          }
          widx++;
        }
      }
    }

    for ( size_t i = 0; i < wirecol->size(); ++i ) {

      int rdIter=dgvindx[i];
      art::Ptr<raw::RawDigit> digitVec(digitVecHandle, rdIter);

      //if ( ! util::CreateAssn(*this, evt, *wirecol, digitVec, *WireDigitAssn, fSpillName, i) ) {
      if ( ! util::CreateAssn(evt, *wirecol, digitVec, *WireDigitAssn, fSpillName, i) ) {
        throw art::Exception(art::errors::ProductRegistrationFailure)
          << "Can't associate wire #" << i
          << " with raw digit #" << digitVec.key();
      }
    }

    evt.put(std::move(wirecol), fSpillName);
    evt.put(std::move(WireDigitAssn), fSpillName);

    fEventCount++;

    return;
} // produce
    
float mtDecon1DROI::fixTheFreakingWaveform(const std::vector<float>& waveform, raw::ChannelID_t channel, std::vector<float>& fixedWaveform) const
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
    
    return localRMS;
}

void mtDecon1DROI::Deconvolute(size_t rdIter, 
			       art::Handle< std::vector<raw::RawDigit> >&digitVecHandle,std::vector<recob::Wire>& wirecol,
                               art::Assns<raw::RawDigit,recob::Wire>& WireDigitAssn, art::Event& evt, std::vector<int>& dgvindx,
                               void* fplan, void* rplan)const{

  // get the reference to the current raw::RawDigit
  art::Ptr<raw::RawDigit> digitVec(digitVecHandle, rdIter);

  raw::ChannelID_t channel = raw::InvalidChannelID; // channel number
  channel = digitVec->Channel();

  const lariov::ChannelStatusProvider& chanFilt = art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider();
  if (!chanFilt.IsPresent(channel)) return; // temporary until "correct" solution implemented
  if (digitVec->GetPedestal() < 0.) return; // Testing an idea about rejecting channels

  float pedestal = 0.;
  recob::Wire::RegionsOfInterest_t ROIVec;
  
  // skip bad channels
  if( chanFilt.Status(channel) >= fMinAllowedChanStatus)
  {
    // uncompress the data
    int dataSize = digitVec->Samples();
    std::vector<short> rawadc(dataSize);
    raw::Uncompress(digitVec->ADCs(), rawadc, digitVec->Compression());
  
    // Get the pedestal subtracted data, centered in the deconvolution vector
    // get pedestal conditions
    const lariov::DetPedestalProvider& pedestalRetrievalAlg = art::ServiceHandle<lariov::DetPedestalService>()->GetPedestalProvider();
    pedestal = pedestalRetrievalAlg.PedMean(channel);
    std::vector<float> rawAdcLessPedVec(dataSize);
    std::transform(rawadc.begin(),rawadc.end(),rawAdcLessPedVec.begin(),std::bind(std::minus<short>(),std::placeholders::_1,pedestal));
  
    // It seems there are deviations from the pedestal when using wirecell for noise filtering
    float raw_noise = fixTheFreakingWaveform(rawAdcLessPedVec, channel, rawAdcLessPedVec);
  
    int transformSize = fSignalShaper->GetFFTSize();

    // Do the deconvolution on the full waveform
    std::vector<float> deconvolvedWaveform(transformSize);
  
    int binOffset    = transformSize > dataSize ? (transformSize - dataSize) / 2 : 0;
    float  deconNorm	   = fSignalShaper->GetDeconNorm();
    float  normFactor	   = 1. / deconNorm; // This is what we had previously: (samplingRate * deconNorm);
    bool   applyNormFactor = std::abs(normFactor - 1.) > std::numeric_limits<float>::epsilon() ? true : false;
  
    // Copy the input (assumed pedestal subtracted) waveforms into our zero padded deconvolution buffer
    std::copy(rawAdcLessPedVec.begin(),rawAdcLessPedVec.end(),deconvolvedWaveform.begin()+binOffset);
  
    // Strategy is to run deconvolution on the entire channel and then pick out the ROI's we found above
    util::LArFFTW fft(transformSize, fplan, rplan, 0);
    fft.Convolute(deconvolvedWaveform, fSignalShaper->GetDeconvKernel(channel));

    //negative number
    int time_offset = fSignalShaper->FieldResponseTOffset(channel);
    std::vector<float> temp;
    if (time_offset <=0){
      temp.assign(deconvolvedWaveform.end()+time_offset,deconvolvedWaveform.end());
      deconvolvedWaveform.erase(deconvolvedWaveform.end()+time_offset,deconvolvedWaveform.end());
      deconvolvedWaveform.insert(deconvolvedWaveform.begin(),temp.begin(),temp.end());
    }else{
      temp.assign(deconvolvedWaveform.begin(),deconvolvedWaveform.begin()+time_offset);
      deconvolvedWaveform.erase(deconvolvedWaveform.begin(),deconvolvedWaveform.begin()+time_offset);
      deconvolvedWaveform.insert(deconvolvedWaveform.end(),temp.begin(),temp.end());
    }

    if (applyNormFactor) std::transform(deconvolvedWaveform.begin(),deconvolvedWaveform.end(),deconvolvedWaveform.begin(), std::bind(std::multiplies<float>(),std::placeholders::_1,normFactor));
  
    // Get the truncated mean and rms
    float truncMean;
    int   nTrunc;
  
    fWaveformTool->getTruncatedMean(deconvolvedWaveform, truncMean, nTrunc);
  
    std::transform(deconvolvedWaveform.begin(),deconvolvedWaveform.end(),deconvolvedWaveform.begin(), std::bind(std::minus<float>(),std::placeholders::_1,truncMean));

    // apply wire-by-wire calibration
    if (fDodQdxCalib){
      if(fdQdxCalib.find(channel) != fdQdxCalib.end()){
  	float constant = fdQdxCalib.at(channel);
  	for (size_t iholder = 0; iholder < deconvolvedWaveform.size(); ++iholder){
  	  deconvolvedWaveform[iholder] *= constant;
  	}
      }
    }

    // Now find the candidate ROI's
    std::vector<geo::WireID> wids    = fGeometry->ChannelToWire(channel);
    const geo::PlaneID&      planeID = wids[0].planeID();
    icarus_tool::IROIFinder::CandidateROIVec candRoiVec;
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
  	double base = fBaseline->GetBaseline(holder, channel, 0, roiLen);
  	std::transform(holder.begin(),holder.end(),holder.begin(),[base](double& adcVal){return adcVal - base;});

  	// add the range into ROIVec
  	ROIVec.add_range(candROI.first, std::move(holder));
    }
  } // end if not a bad channel

  // Don't save empty wires
  if (!ROIVec.empty()){
    // save them
    wirecol[rdIter]=recob::WireCreator(std::move(ROIVec),*digitVec).move();
    dgvindx[rdIter]=rdIter;
  }
}

//////////////////////////////////////////////////////

} // end namespace caldata
