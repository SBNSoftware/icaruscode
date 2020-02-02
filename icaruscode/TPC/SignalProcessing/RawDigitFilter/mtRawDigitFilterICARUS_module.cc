#include <cmath>
#include <algorithm>
#include <vector>

#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"

#include "TComplex.h"

#include "art/Framework/Core/ReplicatedProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalService.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalProvider.h"
#include "lardata/Utilities/LArFFTWPlan.h"
#include "lardata/Utilities/LArFFTW.h"

#include "icaruscode/TPC/SignalProcessing/RawDigitFilter/Algorithms/RawDigitNoiseFilterDefs.h"
#include "icaruscode/TPC/SignalProcessing/RawDigitFilter/Algorithms/RawDigitBinAverageAlg.h"
#include "icaruscode/TPC/SignalProcessing/RawDigitFilter/Algorithms/RawDigitCharacterizationAlg.h"
#include "icaruscode/TPC/SignalProcessing/RawDigitFilter/Algorithms/RawDigitCorrelatedCorrectionAlg.h"
#include "icaruscode/TPC/SignalProcessing/RawDigitFilter/Algorithms/IRawDigitFilter.h"
#include "icaruscode/TPC/SignalProcessing/RawDigitFilter/Algorithms/ChannelGroups.h"
#include "icaruscode/Utilities/tools/IFilter.h"

#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"

using std::vector;
using std::cout;
using std::endl;

//  raw::ChannelID_t channel;
struct WireChar {
  float truncMean;
  float truncRms;
  short mean;
  short median;
  short mode;
  float skewness;
  float fullRms;
  short minMax;
  float neighborRatio;
  float pedCor;
  unsigned int tcka;
  unsigned int tckb;
  unsigned int plane;
  unsigned int wire;
  unsigned int wireIdx;
  raw::ChannelID_t channel;
  int irawdig;
};

struct GroupWireDigIndx {
  int group;
  int windx;
  int irawdig;
  int qgroup;
  unsigned int qgindx;
};

class RawDigitFilterICARUS : public art::ReplicatedProducer
{
public:

    // Copnstructors, destructor.
    explicit RawDigitFilterICARUS(fhicl::ParameterSet const & pset, art::ProcessingFrame const& frame);
    virtual ~RawDigitFilterICARUS();

    // Overrides.
    virtual void configure(fhicl::ParameterSet const & pset);
    virtual void produce(art::Event & e, art::ProcessingFrame const& frame);
    virtual void beginJob(art::ProcessingFrame const& frame);
    virtual void endJob(art::ProcessingFrame const& frame);
    void WaveformChar(unsigned int i, unsigned int& fDataSize, unsigned int& fftsize, void* fplan, void* rplan,
                      vector<GroupWireDigIndx>& igwvec,
                      std::vector<const raw::RawDigit*>& rawDigitVec,
                      vector<vector<caldata::RawDigitVector>>& rawadcgvec,
                      vector<vector<WireChar>>& wgcvec,
                      vector<vector<vector <int>>>& wgqvec,
                      std::unique_ptr<std::vector<raw::RawDigit> >& filteredRawDigit)const;
    void RemoveCorrelatedNoise(unsigned int igrp, unsigned int& fftSize, unsigned int& halfFFTSize, void* fplan, void* rplan,
                               vector<vector<caldata::RawDigitVector>>& rawadcgvec,
                               vector<vector<WireChar>>& wgcvec,
                               vector<vector<vector <int>>>& wgqvec,
                               std::unique_ptr<std::vector<raw::RawDigit> >& filteredRawDigit)const;

private:

    template <typename T> void findPeaks(typename std::vector<T>::iterator startItr,
                                         typename std::vector<T>::iterator stopItr,
                                         std::vector<std::tuple<size_t,size_t,size_t>>& peakTupleVec,
                                         T                                 threshold,
                                         size_t                            firstTick) const;

    void saveRawDigits(std::unique_ptr<std::vector<raw::RawDigit> >&, raw::ChannelID_t&, caldata::RawDigitVector&, float, float);

    // Fcl parameters.
    std::string          fDigitModuleLabel;      ///< The full collection of hits
    bool                 fTruncateTicks;         ///< If true then drop ticks off ends of wires
    bool                 fTruncateChannels;      ///< If true then we drop channels with "no signal"
    bool                 fDoFFTCorrection;       ///< Run the FFT noise correction
    bool                 fSmoothCorrelatedNoise; ///< Should we smooth the noise?
    bool                 fDoCorrelatedNoise;     ///< Process the noise
    bool                 fApplyCorSmoothing;     ///< Attempt to smooth the correlated noise correction?
    bool                 fApplyFFTCorrection;    ///< Use an FFT to get the correlated noise correction
    bool                 fApplyTopHatFilter;     ///< Apply the top hat filter
    unsigned int         fWindowSize;            ///< # ticks to keep in window
    unsigned int         fNumTicksToDropFront;   ///< # ticks to drop from front of waveform
    float                fTruncMeanFraction;     ///< Fraction for truncated mean
    std::vector<size_t>  fNumWiresToGroup;       ///< If smoothing, the number of wires to look at
    std::vector<short>   fMinMaxSelectionCut;    ///< Plane by plane cuts for spread cut
    std::vector<float>   fNRmsChannelReject;     ///< # rms to reject channel as no signal
    std::vector<float>   fRmsRejectionCutHi;     ///< Maximum rms for input channels, reject if larger
    std::vector<float>   fRmsRejectionCutLow;    ///< Minimum rms to consider channel "alive"
    std::vector<float>   fNumRmsToSmoothVec;     ///< # "sigma" to smooth correlated correction vec
    std::vector<double>  fFFTMinPowerThreshold;  ///< Threshold for trimming FFT power spectrum

    // Statistics.
    int fNumEvent;        ///< Number of events seen.

    // Correction algorithms
    caldata::RawDigitBinAverageAlg               fBinAverageAlg;
    caldata::RawDigitCharacterizationAlg       fCharacterizationAlg;
    caldata::RawDigitCorrelatedCorrectionAlg   fCorCorrectAlg;

    std::unique_ptr<caldata::IRawDigitFilter>    fRawDigitFilterTool;

    std::map<size_t,std::vector<std::complex<double>>>      fFilterVec;
    std::map<size_t,std::unique_ptr<icarus_tool::IFilter>> fFilterToolMap;

    // Useful services, keep copies for now (we can update during begin run periods)
    geo::GeometryCore const*           fGeometry;             ///< pointer to Geometry service
    detinfo::DetectorProperties const* fDetectorProperties;   ///< Detector properties service
    const lariov::DetPedestalProvider& fPedestalRetrievalAlg; ///< Keep track of an instance to the pedestal retrieval alg

    // mwang added
    caldata::ChannelGroups fChannelGroups;
};

DEFINE_ART_MODULE(RawDigitFilterICARUS)

//----------------------------------------------------------------------------
class lartbb_WaveformChar {
  public:
    lartbb_WaveformChar(RawDigitFilterICARUS const & prod,
      unsigned int & fdatasize,
      unsigned int & fftsize,
      void* fplan, void* rplan,
      vector<GroupWireDigIndx>& igwv,
      std::vector<const raw::RawDigit*>& rawdigitvec,
      vector<vector<caldata::RawDigitVector>>& rawadcgv,
      vector<vector<WireChar>>& wgcv,
      vector<vector<vector <int>>>& wgqv,
      std::unique_ptr<std::vector<raw::RawDigit> >& filteredrawdigit)
      : prod(prod),
        fDataSize(fdatasize),
        fftSize(fftsize),
        fplan(fplan),
        rplan(rplan),
        igwvec(igwv),
        rawDigitVec(rawdigitvec),
        rawadcgvec(rawadcgv),
        wgcvec(wgcv),
        wgqvec(wgqv),
        filteredRawDigit(filteredrawdigit){}
    void operator()(const tbb::blocked_range<size_t>& range) const{
      //std::cout << " !!!!!!!!!! range.begin(): " << range.begin() << " and range.end(): " << range.end() << std::endl;
      for (size_t i = range.begin(); i < range.end(); ++i)
        prod.WaveformChar(i, fDataSize, fftSize, fplan, rplan, igwvec, rawDigitVec, rawadcgvec, wgcvec, wgqvec, filteredRawDigit);
    }
  private:
    RawDigitFilterICARUS const & prod;
    unsigned int & fDataSize;
    unsigned int & fftSize;
    void* fplan;
    void* rplan;
    vector<GroupWireDigIndx>& igwvec;
    std::vector<const raw::RawDigit*>& rawDigitVec;
    vector<vector<caldata::RawDigitVector>>& rawadcgvec;
    vector<vector<WireChar>>& wgcvec;
    vector<vector<vector <int>>>& wgqvec;
    std::unique_ptr<std::vector<raw::RawDigit> >& filteredRawDigit;
};

//----------------------------------------------------------------------------
class lartbb_RemoveCorrelatedNoise {
  public:
    lartbb_RemoveCorrelatedNoise(RawDigitFilterICARUS const & prod,
      unsigned int & fftsize,
      unsigned int & halffftsize,
      void* fplan, void* rplan,
      vector<vector<caldata::RawDigitVector>>& rawadcgv,
      vector<vector<WireChar>>& wgcv,
      vector<vector<vector <int>>>& wgqv,
      std::unique_ptr<std::vector<raw::RawDigit> >& filteredrawdigit)
      : prod(prod),
        fftSize(fftsize),
        halfFFTSize(halffftsize),
        fplan(fplan),
        rplan(rplan),
        rawadcgvec(rawadcgv),
        wgcvec(wgcv),
        wgqvec(wgqv),
        filteredRawDigit(filteredrawdigit){}
    void operator()(const tbb::blocked_range<size_t>& range) const{
      for (size_t i = range.begin(); i < range.end(); ++i)
        prod.RemoveCorrelatedNoise(i, fftSize, halfFFTSize, fplan, rplan, rawadcgvec, wgcvec, wgqvec, filteredRawDigit);
    }
  private:
    RawDigitFilterICARUS const & prod;
    unsigned int & fftSize;
    unsigned int & halfFFTSize;
    void* fplan;
    void* rplan;
    vector<vector<caldata::RawDigitVector>>& rawadcgvec;
    vector<vector<WireChar>>& wgcvec;
    vector<vector<vector <int>>>& wgqvec;
    std::unique_ptr<std::vector<raw::RawDigit> >& filteredRawDigit;
};

//----------------------------------------------------------------------------
RawDigitFilterICARUS::RawDigitFilterICARUS(fhicl::ParameterSet const & pset, art::ProcessingFrame const& frame) :
                      art::ReplicatedProducer(pset, frame),
                      fNumEvent(0),
                      fBinAverageAlg(pset),
                      fCharacterizationAlg(pset.get<fhicl::ParameterSet>("CharacterizationAlg")),
                      fCorCorrectAlg(pset.get<fhicl::ParameterSet>("CorrelatedCorrectionAlg")),
                      fPedestalRetrievalAlg(*lar::providerFrom<lariov::DetPedestalService>()),
                      fChannelGroups(pset)
{

    fGeometry = lar::providerFrom<geo::Geometry>();
    fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();

    configure(pset);
    produces<std::vector<raw::RawDigit> >();

    // Report.
    mf::LogInfo("RawDigitFilterICARUS") << "RawDigitFilterICARUS configured\n";
}

//----------------------------------------------------------------------------
RawDigitFilterICARUS::~RawDigitFilterICARUS()
{}

//----------------------------------------------------------------------------
void RawDigitFilterICARUS::configure(fhicl::ParameterSet const & pset)
{
    fDigitModuleLabel      = pset.get<std::string>        ("DigitModuleLabel",                                        "daq");
    fWindowSize            = pset.get<size_t>             ("WindowSize",                                               6400);
    fTruncateTicks         = pset.get<bool>               ("TruncateTicks",                                           false);
    fNumTicksToDropFront   = pset.get<size_t>             ("NumTicksToDropFront",                                      2400);
    fTruncateChannels      = pset.get<bool>               ("TruncateChannels",                                        false);
    fDoFFTCorrection       = pset.get<bool>               ("FFTNoise",                                                 true);
    fNRmsChannelReject     = pset.get<std::vector<float>> ("NRMSChannelReject",     std::vector<float>() = {3.,  3.,  3.  });
    fSmoothCorrelatedNoise = pset.get<bool>               ("SmoothCorrelatedNoise",                                    true);
    // .. waveform characterization
    fRmsRejectionCutHi     = pset.get<std::vector<float>> ("RMSRejectionCutHi",     std::vector<float>() = {25.0,25.0,25.0});
    fRmsRejectionCutLow    = pset.get<std::vector<float>> ("RMSRejectionCutLow",    std::vector<float>() = {0.70,0.70,0.70});
    fMinMaxSelectionCut    = pset.get<std::vector<short>> ("MinMaxSelectionCut",    std::vector<short>() = {13, 13, 11});
    // .. correlated noise correction
    fNumWiresToGroup       = pset.get<std::vector<size_t>>("NumWiresToGroup",          std::vector<size_t>() = {48, 48, 96});
    fDoCorrelatedNoise     = pset.get<bool>               ("ProcessNoise",                                             true);
    fApplyCorSmoothing     = pset.get<bool>               ("ApplyCorSmoothing",                                        true);
    fTruncMeanFraction     = pset.get<float>              ("TruncMeanFraction",                                        0.15);
    fNumRmsToSmoothVec     = pset.get<std::vector<float>> ("NumRmsToSmooth",          std::vector<float>() = {3.6, 3.6, 4.});
    fApplyFFTCorrection    = pset.get<bool>               ("ApplyFFTCorrection",                                       true);
    fFFTMinPowerThreshold  = pset.get<std::vector<double>>("FFTPowerThreshold",     std::vector<double>() = {100.,75.,500.});
    fApplyTopHatFilter     = pset.get<bool>               ("ApplyTopHatFilter",                                        true);

    fRawDigitFilterTool = art::make_tool<caldata::IRawDigitFilter>(pset.get<fhicl::ParameterSet>("RawDigitFilterTool"));

    // Implement the tools for handling the responses
    const fhicl::ParameterSet& filterTools = pset.get<fhicl::ParameterSet>("FilterTools");

    for(const std::string& filterTool : filterTools.get_pset_names())
    {
        const fhicl::ParameterSet& filterToolParamSet = filterTools.get<fhicl::ParameterSet>(filterTool);
        size_t                     planeIdx           = filterToolParamSet.get<size_t>("Plane");

        fFilterToolMap.insert(std::pair<size_t,std::unique_ptr<icarus_tool::IFilter>>(planeIdx,art::make_tool<icarus_tool::IFilter>(filterToolParamSet)));
    }
}

//----------------------------------------------------------------------------
void RawDigitFilterICARUS::beginJob(art::ProcessingFrame const&)
{
    // Access ART's TFileService, which will handle creating and writing
    // histograms and n-tuples for us.
    art::ServiceHandle<art::TFileService> tfs;

    fCharacterizationAlg.initializeHists(tfs);
    fCorCorrectAlg.initializeHists(tfs);

    art::TFileDirectory dir = tfs->mkdir(Form("RawDigitFilter"));

    fRawDigitFilterTool->initializeHistograms(dir);

    return;
}

//----------------------------------------------------------------------------
void RawDigitFilterICARUS::produce(art::Event & event, art::ProcessingFrame const&)
{
  ++fNumEvent;

  // Read in the digit List object(s).
  art::Handle< std::vector<raw::RawDigit> > digitVecHandle;
  event.getByLabel(fDigitModuleLabel, digitVecHandle);

  // Agreed convention is to ALWAYS output to the event store so get a pointer to our collection
  std::unique_ptr<std::vector<raw::RawDigit> > filteredRawDigit(new std::vector<raw::RawDigit>);
  filteredRawDigit->resize(digitVecHandle->size());
  //std::cout << " ~~~~~ Initial size of filteredRawDigit: " << filteredRawDigit->size() << std::endl;

  // ... Require a valid handle
  if (digitVecHandle.isValid() && digitVecHandle->size()>0 ){
    unsigned int maxChannels	= fGeometry->Nchannels();
    //unsigned int maxTimeSamples = fDetectorProperties->NumberTimeSamples();

    // .. Let's first sort the rawDigitVec
    std::vector<const raw::RawDigit*> rawDigitVec;
    for(size_t idx = 0; idx < digitVecHandle->size(); idx++) rawDigitVec.push_back(&digitVecHandle->at(idx));
    std::sort(rawDigitVec.begin(),rawDigitVec.end(),
    [](const raw::RawDigit* left, const raw::RawDigit* right) {return left->Channel() < right->Channel();});

    // .. Use the handle to get a particular (0th) element of collection.
    art::Ptr<raw::RawDigit> digitVec0(digitVecHandle, 0);
    unsigned int fDataSize = digitVec0->Samples(); //size of raw data vectors
    unsigned int fftSize;
    if (fTruncateTicks) {
      fftSize = fWindowSize;
    } else {
      fftSize = fDataSize;
    }
    vector<short> rawadc;
    rawadc.resize(fftSize);

    vector<vector<caldata::RawDigitVector>> rawadcgvec;
    vector<vector<WireChar>> wgcvec;
    vector<vector<vector <int>>> wgqvec;
    vector<GroupWireDigIndx>igwvec;

    vector<caldata::RawDigitVector> rawadcvec;
    vector<WireChar> wcvec;
    vector<vector<int>>wqvec;
    wqvec.push_back({});
    wqvec.push_back({});
    int irawdig=-1;

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // ... Do a first loop over all the rawDigits to set up the grouped vectors
    //     for:
    //     wgcvec: wire charactestics
    //     wgqvec: wire quality
    //     rawdcgvec: uncompressed raw adcs
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for(const auto& rawDigit : rawDigitVec){

      irawdig++;

      raw::ChannelID_t channel = rawDigit->Channel();
      bool goodChan(true);
      std::vector<geo::WireID> wids;
      try {
        wids = fGeometry->ChannelToWire(channel);
      }
      catch(...){
        goodChan = false;
      }
      if (channel >= maxChannels || !goodChan) continue;

      unsigned int plane = wids[0].Plane;
      unsigned int wire  = wids[0].Wire;

      // .. Verify that dataSize looks fine
      unsigned int dataSize = rawDigit->Samples();
      if (dataSize < 1){
        std::cout << "****>> Found zero length raw digit buffer, channel: "
	          << channel << ", plane: " << plane << ", wire: " << wire << std::endl;
        continue;
      }else if (dataSize!=fDataSize) {
        std::cout << "****>> DataSize has changed from " << fDataSize << " to " << dataSize
	          << " for channel: " << channel << ", plane: " << plane << ", wire: " << wire << std::endl;
        continue;
      }
      rawadcvec.push_back(rawadc);

      unsigned int wireIdx  = wire % fNumWiresToGroup[plane];
      //std::cout << "Channel = " << channel << " wire = " << wire
      //          << " plane = " << plane << " wireIdx = " << wireIdx << std::endl;

      WireChar wc;
      wc.wire = wire;
      wc.plane = plane;
      wc.channel = channel;
      wc.wireIdx = wireIdx;
      wc.irawdig = irawdig;
      wcvec.push_back(wc);

      GroupWireDigIndx igw;
      igw.windx=int(wcvec.size()-1);
      igw.irawdig=irawdig;
      igw.group=int(wgcvec.size());

      igw.qgroup=-1;
      size_t group = fChannelGroups.channelGroup(plane,wire);
      if (group==0) {
        wqvec[0].push_back(wcvec.size()-1);
        igw.qgroup=0;
        igw.qgindx=wqvec[0].size()-1;
      } else if (group==1){
        wqvec[1].push_back(wcvec.size()-1);
        igw.qgroup=1;
        igw.qgindx=wqvec[1].size()-1;
      }

      igwvec.push_back(igw);

      // Are we at the correct boundary for dealing with the noise?
      if (!((wireIdx + 1) % fNumWiresToGroup[plane])){
    	//std::cout << "!!!! Reached boundary to group wires" << std::endl;
        rawadcgvec.push_back(rawadcvec);
        wgcvec.push_back(wcvec);
        wgqvec.push_back(wqvec);
        rawadcvec.clear();
        wcvec.clear();
        wqvec[0].clear();
        wqvec[1].clear();
      }
    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // ... Now that we have set up the data structures above, we can peform
    //     the FFT correction and determine the waveform parameters for each
    //     individual wire.
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // .. First set up the filters
    unsigned int halfFFTSize(fftSize/2 + 1);
    for(unsigned int plne = 0; plne < 3; plne++){
      fFilterToolMap.at(plne)->setResponse(fftSize,1.,1.);
      const std::vector<TComplex>& filter = fFilterToolMap.at(plne)->getResponseVec();
      fFilterVec[plne] = std::vector<std::complex<double>>();
      fFilterVec.at(plne).reserve(halfFFTSize);
      for(auto& rootComplex : filter) fFilterVec.at(plne).emplace_back(rootComplex.Re(),rootComplex.Im());
    }

    // .. Now set up the fftw plan
    util::LArFFTWPlan lfftwp(fftSize,"ES");

    //int nwavedump = 0;

    //for (std::size_t i=0; i<igwvec.size(); i++){
    //  WaveformChar(i, fDataSize, igwvec, rawDigitVec, rawadcgvec, wgcvec, filteredRawDigit);
    //}
    // ... Launch multiple threads with TBB to do the waveform characterization and fft correction in parallel
    auto func = lartbb_WaveformChar(*this, fDataSize, fftSize, lfftwp.fPlan, lfftwp.rPlan, igwvec, rawDigitVec,
                                    rawadcgvec, wgcvec, wgqvec, filteredRawDigit);
    tbb::parallel_for(tbb::blocked_range<size_t>(0, igwvec.size()), func);

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // ... Next, we can do the correlated noise correction for each wire group
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (fDoCorrelatedNoise && fSmoothCorrelatedNoise){

      // .. Loop over each group of wires
      //for (size_t igrp = 0; igrp < wgcvec.size(); igrp++) {
      //  RemoveCorrelatedNoise(igrp, fftSize, halfFFTSize, lfftwp.fPlan, lfftwp.rPlan, rawadcgvec, wgcvec, wgqvec, filteredRawDigit);
      //} // loop over igrp
      auto func = lartbb_RemoveCorrelatedNoise(*this, fftSize, halfFFTSize, lfftwp.fPlan, lfftwp.rPlan,
                                               rawadcgvec, wgcvec, wgqvec, filteredRawDigit);
      tbb::parallel_for(tbb::blocked_range<size_t>(0, wgcvec.size()), func);
    } // if do and smooth correlated noise

    filteredRawDigit->erase(std::remove_if(filteredRawDigit->begin(),filteredRawDigit->end(),
                            [](const raw::RawDigit & frd){return frd.ADCs().size()==0;}),
			    filteredRawDigit->end());
    
    rawadcgvec.clear();
    wgcvec.clear();
    wgqvec.clear();
    igwvec.clear();

  }

  //std::cout << " ~~~~~ Final size of filteredRawDigit: " << filteredRawDigit->size() << std::endl;

  // Add tracks and associations to event.
  event.put(std::move(filteredRawDigit));
}

//----------------------------------------------------------------------------
void RawDigitFilterICARUS::RemoveCorrelatedNoise(unsigned int igrp, unsigned int& fftSize, unsigned int& halfFFTSize, void* fplan, void* rplan,
                                                 vector<vector<caldata::RawDigitVector>>& rawadcgvec,
                                                 vector<vector<WireChar>>& wgcvec,
                                                 vector<vector<vector <int>>>& wgqvec,
                                                 std::unique_ptr<std::vector<raw::RawDigit> >& filteredRawDigit)const{

  for (size_t iq = 0; iq < 2; iq++) {

    if (wgqvec[igrp][iq].size() == 0) continue;

    // .. Remove the unclassified wires whose indices have been set to negative
    wgqvec[igrp][iq].erase(std::remove_if(wgqvec[igrp][iq].begin(),wgqvec[igrp][iq].end(),[](const int& i){return i<0;}),wgqvec[igrp][iq].end());

    // .. Don't try to do correction if too few wires unless they have gaps
    size_t nwq = wgqvec[igrp][iq].size();
    if (nwq <= 2) continue;

    std::vector<float> corValVec;
    corValVec.resize(fftSize, 0.);

    // ----------------------------------------------------
    // .. Build the vector of corrections for each time bin
    // ----------------------------------------------------
    for(size_t itck = 0; itck < fftSize; itck++){
      std::vector<float> adcValuesVec;
      // .. Loop over each entry in wgqvec[igrp][iq]
      for (size_t i = 0; i < nwq; i++) {
  	// .. get index into the wgvec[igrp] array
  	size_t iwdx = wgqvec[igrp][iq][i];
  	// .. Check that we should be doing something in this range
  	//    Note that if the wire is not to be considered then the "start" bin will be after the last bin
  	if (itck < wgcvec[igrp][iwdx].tcka || itck >= wgcvec[igrp][iwdx].tckb) continue;
  	// .. Accumulate
  	adcValuesVec.push_back(float(rawadcgvec[igrp][iwdx][itck]) - wgcvec[igrp][iwdx].truncMean);
      }
      // ... Get the median for this time tick across all wires in the group
      float medval(-10000);
      if (!adcValuesVec.empty()) {
  	  std::sort(adcValuesVec.begin(),adcValuesVec.end());
  	  size_t medidx = adcValuesVec.size() / 2;
  	  medval = adcValuesVec[medidx];
  	  if (adcValuesVec.size() > 1 && medidx % 2) medval = (medval + adcValuesVec[medidx+1]) / 2;
      }
      corValVec[itck] = std::max(medval,float(-10000.));
    } // loop over itck

    // .. get the plane number for first wire in this set, for use below
    size_t iwdx0 = wgqvec[igrp][iq][0];
    unsigned int plane = wgcvec[igrp][iwdx0].plane;

    // --------------------------------------
    // ... Try to eliminate any real outliers
    // --------------------------------------
    if (fApplyCorSmoothing) {
      std::vector<float> localCorValVec = corValVec;
      std::sort(localCorValVec.begin(),localCorValVec.end());

      int   nTruncVal  = (1. - fTruncMeanFraction) * localCorValVec.size();
      float corValSum  = std::accumulate(localCorValVec.begin(),localCorValVec.begin() + nTruncVal,0.);
      float meanCorVal = corValSum / float(nTruncVal);

      std::vector<float> diffVec(nTruncVal);
      std::transform(localCorValVec.begin(),localCorValVec.begin() + nTruncVal, diffVec.begin(),
  		     std::bind(std::minus<float>(),std::placeholders::_1,meanCorVal));

      float rmsValSq   = std::inner_product(diffVec.begin(),diffVec.end(),diffVec.begin(),0.);
      float rmsVal     = std::sqrt(rmsValSq / float(nTruncVal));

      // .. Now set up to run through and do a "simple" interpolation over outliers
      std::vector<float>::iterator lastGoodItr = corValVec.begin();
      bool wasOutlier(false);
      for(std::vector<float>::iterator corValItr = lastGoodItr+1; corValItr != corValVec.end(); corValItr++){
  	if (fabs(*corValItr - meanCorVal) < fNumRmsToSmoothVec.at(plane)*rmsVal){
  	      if (wasOutlier){
  		float lastVal  = *lastGoodItr;
  		float curVal   = *corValItr;
  		float numTicks = std::distance(lastGoodItr,corValItr);
  		float slope    = (curVal - lastVal) / numTicks;
  		while(lastGoodItr != corValItr){
  		  *lastGoodItr++ = (numTicks - std::distance(lastGoodItr,corValItr)) * slope + lastVal;
  		}
  	      }
  	      wasOutlier  = false;
  	      lastGoodItr = corValItr;
  	} else {
  	  wasOutlier = true;
  	}
      }
    }  // fApplyCorSmoothing

    // ... Get the FFT correction
    if (fApplyFFTCorrection) {
      std::vector<std::complex<double>> fftOutputVec(halfFFTSize);
      util::LArFFTW lfftw(fftSize, fplan, rplan, 0);
      lfftw.DoFFT(corValVec, fftOutputVec);

      std::vector<double> powerVec(halfFFTSize);
      std::transform(fftOutputVec.begin(), fftOutputVec.begin() + halfFFTSize, powerVec.begin(), [](const auto& val){return std::abs(val);});

      // Want the first derivative
      std::vector<double> firstDerivVec(powerVec.size(), 0.);
    
      //fWaveformTool->firstDerivative(powerVec, firstDerivVec);
      for(size_t idx = 1; idx < firstDerivVec.size() - 1; idx++)
          firstDerivVec.at(idx) = 0.5 * (powerVec.at(idx + 1) - powerVec.at(idx - 1));

      // Find the peaks
      std::vector<std::tuple<size_t,size_t,size_t>> peakTupleVec;
    
      findPeaks(firstDerivVec.begin(),firstDerivVec.end(),peakTupleVec,fFFTMinPowerThreshold[plane],0);
    
      if (!peakTupleVec.empty())
      {
          for(const auto& peakTuple : peakTupleVec)
          {
              size_t startTick = std::get<0>(peakTuple);
              size_t stopTick  = std::get<2>(peakTuple);
        
              if (stopTick > startTick)
              {
        	  std::complex<double> slope = (fftOutputVec[stopTick] - fftOutputVec[startTick]) / double(stopTick - startTick);
            
        	  for(size_t tick = startTick; tick < stopTick; tick++)
        	  {
        	      std::complex<double> interpVal = fftOutputVec[startTick] + double(tick - startTick) * slope;
        	
        	      fftOutputVec[tick]		   = interpVal;
        	      //fftOutputVec[fftDataSize - tick - 1] = interpVal;
        	  }
              }
          }
      
          std::vector<double> tmpVec(corValVec.size());
      
          lfftw.DoInvFFT(fftOutputVec, tmpVec);
      
          std::transform(corValVec.begin(),corValVec.end(),tmpVec.begin(),corValVec.begin(),std::minus<double>());
      }
    } // fApplyFFTCorrection

    // -------------------------------------------
    // ... Now go through and apply the correction
    // -------------------------------------------
    for(size_t itck = 0; itck < fftSize; itck++){
      for (size_t i = 0; i < nwq; i++) {
  	float corVal;
  	size_t iwdx = wgqvec[igrp][iq][i];
  	// .. If the "start" bin is after the "stop" bin then we are meant to skip this wire in the averaging process
  	//    Or if the sample index is in a chirping section then no correction is applied.
  	//    Both cases are handled by looking at the sampleIdx
  	if (itck < wgcvec[igrp][iwdx].tcka || itck >= wgcvec[igrp][iwdx].tckb) {
  	  corVal=0.;
  	} else {
  	  corVal = corValVec[itck];
  	}
  	// .. Probably doesn't matter, but try to get slightly more accuracy by doing float math and rounding
  	float newAdcValueFloat = float(rawadcgvec[igrp][iwdx][itck]) - corVal - wgcvec[igrp][iwdx].pedCor;
  	rawadcgvec[igrp][iwdx][itck] = std::round(newAdcValueFloat);
      }
    }
  } // loop over iq

  // ----------------------------------------------------
  // ... One more pass through to store the good channels
  // ----------------------------------------------------
  for (size_t iwdx = 0; iwdx < wgcvec[igrp].size(); iwdx++) {

    unsigned int plane = wgcvec[igrp][iwdx].plane;

    // Try baseline correction?
    if (fApplyTopHatFilter && plane != 2 && wgcvec[igrp][iwdx].skewness > 0.) {
  	fRawDigitFilterTool->FilterWaveform(rawadcgvec[igrp][iwdx], iwdx, plane);
    }

    // recalculate rms for the output
    float rmsVal   = 0.;
    float pedestal = wgcvec[igrp][iwdx].truncMean;
    float pedCor   = wgcvec[igrp][iwdx].pedCor;
    float deltaPed = pedestal - pedCor;

    caldata::RawDigitVector& rawDataVec = rawadcgvec[igrp][iwdx];
    fCharacterizationAlg.getTruncatedRMS(rawDataVec, deltaPed, rmsVal);

    // The ultra high noise channels are simply zapped
    raw::ChannelID_t channel = wgcvec[igrp][iwdx].channel;
    if (rmsVal < fRmsRejectionCutHi[plane]) { // && ImAGoodWire(plane,baseWireIdx + locWireIdx))
  	int irdg = wgcvec[igrp][iwdx].irawdig;
  	//saveRawDigits(filteredRawDigit, channelWireVec[locWireIdx], rawDataVec, pedestal, rmsVal);
  	filteredRawDigit->at(irdg) = raw::RawDigit(channel, rawDataVec.size(), rawDataVec, raw::kNone);
  	filteredRawDigit->at(irdg).SetPedestal(pedestal, rmsVal);
    } else {
  	mf::LogInfo("RawDigitFilterICARUS") <<  "--> Rejecting channel for large rms, channel: "
  	<< channel << ", rmsVal: " << rmsVal << ", truncMean: " << pedestal
  	<< ", pedestal: " << pedCor << std::endl;
    }
  } // loop over igrp

}

//----------------------------------------------------------------------------
void RawDigitFilterICARUS::WaveformChar(unsigned int i, unsigned int& fDataSize, unsigned int& fftSize, void* fplan, void* rplan,
                                        vector<GroupWireDigIndx>& igwvec,
                                        std::vector<const raw::RawDigit*>& rawDigitVec,
                                        vector<vector<caldata::RawDigitVector>>& rawadcgvec,
                                        vector<vector<WireChar>>& wgcvec,
                                        vector<vector<vector <int>>>& wgqvec,
                                        std::unique_ptr<std::vector<raw::RawDigit> >& filteredRawDigit)const{
  int igrp = igwvec[i].group;
  int iwdx = igwvec[i].windx;
  int irdg = igwvec[i].irawdig;
  const raw::RawDigit* rawDigit = rawDigitVec.at(irdg);

  // .. Uncompress the RawDigit
  caldata::RawDigitVector& rawADC = rawadcgvec[igrp][iwdx];
  caldata::RawDigitVector tempVec(fDataSize);
  if (fTruncateTicks){
    raw::Uncompress(rawDigit->ADCs(), tempVec, rawDigit->Compression());
    std::copy(tempVec.begin() + fNumTicksToDropFront, tempVec.begin() + fNumTicksToDropFront + fWindowSize, rawADC.begin());
  } else {
    raw::Uncompress(rawDigit->ADCs(), rawADC, rawDigit->Compression());
  }

  // .. Do the FFT correction

  raw::ChannelID_t channel = rawDigit->Channel();
  unsigned int plane = wgcvec[igrp][iwdx].plane;

  if (fDoFFTCorrection){
    // .. Subtract the pedestal
    float pedestal = fPedestalRetrievalAlg.PedMean(channel);
    std::vector<float> holder(fftSize);
    std::transform(rawADC.begin(),rawADC.end(),holder.begin(),[pedestal](const auto& val){return float(float(val) - pedestal);});

    // .. Do the correction
    util::LArFFTW lfftw(fftSize, fplan, rplan, 0);
    lfftw.Convolute(holder, fFilterVec.at(plane));

    // .. Restore the pedestal
    std::transform(holder.begin(), holder.end(), rawADC.begin(), [pedestal](const float& adc){return std::round(adc + pedestal);});
  }

  fCharacterizationAlg.getWaveformParams(rawADC,
                                         channel,
                                         plane,
                                         wgcvec[igrp][iwdx].wire,
                                         wgcvec[igrp][iwdx].truncMean,
                                         wgcvec[igrp][iwdx].truncRms,
                                         wgcvec[igrp][iwdx].mean,
                                         wgcvec[igrp][iwdx].median,
                                         wgcvec[igrp][iwdx].mode,
                                         wgcvec[igrp][iwdx].skewness,
                                         wgcvec[igrp][iwdx].fullRms,
                                         wgcvec[igrp][iwdx].minMax,
                                         wgcvec[igrp][iwdx].neighborRatio,
                                         wgcvec[igrp][iwdx].pedCor);

  // This allows the module to be used simply to truncate waveforms with no noise processing
  if (!fDoCorrelatedNoise)
  {
    // Is this channel "quiet" and should be rejected?
    // Note that the "max - min" range is to be compared to twice the rms cut
    if (fTruncateChannels && wgcvec[igrp][iwdx].minMax < 2. * fNRmsChannelReject[plane] * wgcvec[igrp][iwdx].truncRms) return;

    caldata::RawDigitVector pedCorrectedVec;
    pedCorrectedVec.resize(rawADC.size(),0);
    std::transform(rawADC.begin(),rawADC.end(),pedCorrectedVec.begin(),std::bind(std::minus<short>(),std::placeholders::_1,wgcvec[igrp][iwdx].pedCor));

    //saveRawDigits(filteredRawDigit, channel, pedCorrectedVec, truncMeanWireVec[wireIdx], truncRmsWireVec[wireIdx]);
    filteredRawDigit->at(irdg) = raw::RawDigit(channel, pedCorrectedVec.size(), pedCorrectedVec, raw::kNone);
    filteredRawDigit->at(irdg).SetPedestal(wgcvec[igrp][iwdx].truncMean,wgcvec[igrp][iwdx].truncRms);
    return;
  }

  // If we are not performing noise corrections then we are done with this wire
  // Store it and move on
  if (!fSmoothCorrelatedNoise)
  {
    // Filter out the very high noise wires
    if (wgcvec[igrp][iwdx].truncRms < fRmsRejectionCutHi[plane]) {
      //saveRawDigits(filteredRawDigit, channel, rawadc, truncMeanWireVec[wireIdx], truncRmsWireVec[wireIdx]);
      filteredRawDigit->at(irdg) = raw::RawDigit(channel, rawADC.size(), rawADC, raw::kNone);
      filteredRawDigit->at(irdg).SetPedestal(wgcvec[igrp][iwdx].truncMean,wgcvec[igrp][iwdx].truncRms);
    } else {
      // Eventually we'll interface to some sort of channel status communication mechanism.
      // For now use the log file
      mf::LogInfo("RawDigitFilterICARUS") <<  "--> Rejecting channel for large rms, channel: " << channel
      << ", rmsVal: " << wgcvec[igrp][iwdx].truncRms << ", truncMean: " << wgcvec[igrp][iwdx].truncMean
      << ", pedestal: " << wgcvec[igrp][iwdx].pedCor << std::endl;
    }

    return;
  }

  // .. Classify the waveform
  if (wgcvec[igrp][iwdx].minMax > fMinMaxSelectionCut[plane] && wgcvec[igrp][iwdx].truncRms < fRmsRejectionCutHi[plane]){
    wgcvec[igrp][iwdx].tcka = 0;
    wgcvec[igrp][iwdx].tckb = rawADC.size();
    // .. Look for chirping wire sections. Confine this to only the V plane
    if (plane == 1){
      // .. Do wire shape corrections to look for chirping wires & other oddities to avoid. Recover our objects...
      short threshold(6);
      short mean = wgcvec[igrp][iwdx].mean;

      // .. If going from quiescent to on again, then the min/max will be large
      if (wgcvec[igrp][iwdx].skewness > 0. && wgcvec[igrp][iwdx].neighborRatio < 0.7 && wgcvec[igrp][iwdx].minMax > 50){
          raw::RawDigit::ADCvector_t::iterator stopChirpItr = std::find_if(rawADC.begin(),rawADC.end(),
	  			   [mean,threshold](const short& elem){return abs(elem - mean) > threshold;});
          size_t threshIndex = std::distance(rawADC.begin(),stopChirpItr);
          if (threshIndex > 60) wgcvec[igrp][iwdx].tcka = threshIndex;
      } else if (wgcvec[igrp][iwdx].minMax > 20 && wgcvec[igrp][iwdx].neighborRatio < 0.7){ // .. Check in the reverse direction?
          threshold = 3;
          raw::RawDigit::ADCvector_t::reverse_iterator startChirpItr = std::find_if(rawADC.rbegin(),rawADC.rend(),
	  				   [mean,threshold](const short& elem){return abs(elem - mean) > threshold;});
          size_t threshIndex = std::distance(rawADC.rbegin(),startChirpItr);
          if (threshIndex > 60) wgcvec[igrp][iwdx].tckb = rawADC.size() - threshIndex;
      }
    }
  } else {
    // .. If unable to classify, then set the wire index in wgqvec to negative to skip coherent noise correction
    int iqgrp = igwvec[i].qgroup;
    unsigned int iqdx = igwvec[i].qgindx;
    if ( iqgrp == 0 || iqgrp == 1 ) {
      wgqvec[igrp][iqgrp][iqdx]=-1;
    }
    // .. and apply the pedestal correction
    std::transform(rawADC.begin(),rawADC.end(),rawADC.begin(),std::bind(std::minus<short>(),std::placeholders::_1,wgcvec[igrp][iwdx].pedCor));
  }

  return;
}

//----------------------------------------------------------------------------
template <typename T> void RawDigitFilterICARUS::findPeaks(typename std::vector<T>::iterator startItr,
                                                           typename std::vector<T>::iterator stopItr,
                                                           std::vector<std::tuple<size_t,size_t,size_t>>& peakTupleVec,
                                                           T threshold,
                                                           size_t firstTick) const
{
    // Need a minimum distance or else nothing to do
    if (std::distance(startItr,stopItr) > 4)
    {
        // This is a divide and conquer algorithm, start by finding the maximum element.
        typename std::vector<T>::iterator firstItr = std::max_element(startItr,stopItr,[](float left, float right){return std::fabs(left) < std::fabs(right);});

        // Are we over threshold?
        if (std::fabs(*firstItr) > threshold)
        {
            // What am I thinking?
            // First task is to find the "other" lobe max point
            // Set one to the "first", the other to the "second"
            // Search backward from first to find start point, forward from second to find end point
            // Set mid point between first and second as "peak"?
            typename std::vector<T>::iterator secondItr = firstItr;
        
            // Assume if max bin is positive then second lobe is later
            if (*firstItr > 0)
            {
                typename std::vector<T>::iterator tempItr = secondItr;
            
                while(tempItr != stopItr)
                {
                    if (*++tempItr < -threshold)
                    {
                        if (*tempItr < *secondItr) secondItr = tempItr;
                    }
                    else if (secondItr != firstItr) break;
                }
            }
            // Otherwise it goes the other way
            else
            {
                typename std::vector<T>::iterator tempItr = secondItr;
            
                while(tempItr != startItr)
                {
                    if (*--tempItr > threshold)
                    {
                        if (*tempItr > *secondItr) secondItr = tempItr;
                    }
                    else if (secondItr != firstItr) break;
                }
            
                std::swap(firstItr,secondItr);
            }
        
            // It might that no real pulse was found
            if (firstItr != secondItr)
            {
                // Get the "peak" position
                size_t peakBin = std::distance(startItr,firstItr) + std::distance(firstItr,secondItr) / 2;
        
                // Advance (forward or backward) the first and second iterators to get back to zero crossing
                while(firstItr  != startItr) if (*--firstItr  < 0.) break;
                while(secondItr != stopItr)  if (*++secondItr > 0.) break;
        
                size_t firstBin = std::distance(startItr,firstItr);
                size_t lastBin  = std::distance(startItr,secondItr);
        
                // Find leading peaks
                findPeaks(startItr, firstItr, peakTupleVec, threshold, firstTick);
        
                // Save this peak
                peakTupleVec.push_back(std::tuple<size_t,size_t,size_t>(firstBin+firstTick,peakBin+firstTick,lastBin+firstTick));
        
                // Find downstream peaks
                findPeaks(secondItr, stopItr, peakTupleVec, threshold, firstTick + std::distance(startItr,secondItr));
            }
        }
    }

    return;
}

//----------------------------------------------------------------------------
void RawDigitFilterICARUS::saveRawDigits(std::unique_ptr<std::vector<raw::RawDigit> >& filteredRawDigit,
                                         raw::ChannelID_t&                             channel,
                                         caldata::RawDigitVector&                      rawDigitVec,
                                         float                                         pedestal,
                                         float                                         rms)
{
    //filteredRawDigit->emplace_back(raw::RawDigit(channel, rawDigitVec.size(), rawDigitVec, raw::kNone));
    filteredRawDigit->emplace_back(channel, rawDigitVec.size(), rawDigitVec, raw::kNone);
    filteredRawDigit->back().SetPedestal(pedestal,rms);

    return;
}

//----------------------------------------------------------------------------
void RawDigitFilterICARUS::endJob(art::ProcessingFrame const&)
{
    mf::LogInfo("RawDigitFilterICARUS") << "Looked at " << fNumEvent << " events" << std::endl;
}
