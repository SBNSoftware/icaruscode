
#include "icaruscode/Analysis/tools/IWireHistogramTool.h"

#include "fhiclcpp/ParameterSet.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Utilities/make_tool.h"
#include "art_root_io/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalService.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalProvider.h"
#include "larreco/HitFinder/HitFinderTools/IWaveformTool.h"

#include "icaruscode/TPC/Utilities/SignalShapingICARUSService_service.h"
#include "icarussigproc/WaveformTools.h"

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TF1.h"
#include "TVirtualFFT.h"

#include <cmath>
#include <algorithm>

namespace BasicWireAnalysis
{
    ////////////////////////////////////////////////////////////////////////
    //
    // Class:       BasicWireAnalysis
    // Module Type: producer
    // File:        BasicWireAnalysis.h
    //
    //              The intent of this module is to provide methods for
    //              "analyzing" hits on waveforms
    //
    // Configuration parameters:
    //
    // TruncMeanFraction     - the fraction of waveform bins to discard when
    //
    // Created by Tracy Usher (usher@slac.stanford.edu) on February 19, 2016
    //
    ////////////////////////////////////////////////////////////////////////

class BasicWireAnalysis : virtual public IWireHistogramTool
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    explicit BasicWireAnalysis(fhicl::ParameterSet const & pset);
    
    /**
     *  @brief  Destructor
     */
    ~BasicWireAnalysis();
    
    // provide for initialization
    void configure(fhicl::ParameterSet const & pset) override;

    /**
     *  @brief Interface for initializing the histograms to be filled
     *
     *  @param TFileService   handle to the TFile service
     *  @param string         subdirectory to store the hists in
     */
    void initializeHists(art::ServiceHandle<art::TFileService>&, const std::string&) override;
    
    /**
     *  @brief Interface for method to executve at the end of run processing
     *
     *  @param int            number of events to use for normalization
     */
    void endJob(int numEvents) override;
    
    /**
     *  @brief Interface for filling histograms
     */
    void fillHistograms(const IWireHistogramTool::WirePtrVec&, const IWireHistogramTool::SimChannelMap&, int) const override;
    
private:
    
    // Define a structure to contain hits
    using HitCandidate_t = struct HitCandidate
    {
        size_t startTick;
        size_t stopTick;
        size_t maxTick;
        size_t minTick;
        float  maxDerivative;
        float  minDerivative;
        float  hitCenter;
        float  hitSigma;
        float  hitHeight;
    };
    
    using HitCandidateVec      = std::vector<HitCandidate_t>;
//    using MergeHitCandidateVec = std::vector<HitCandidateVec>;
    
    using Waveform = std::vector<float>;
    
    // Internal functions
    void findHitCandidates(Waveform::const_iterator,
                           Waveform::const_iterator,
                           size_t,
                           size_t,
                           HitCandidateVec&) const;
    void findHitCandidates(Waveform::const_iterator,
                           Waveform::const_iterator,
                           Waveform::const_iterator,
                           Waveform::const_iterator,
                           size_t,
                           size_t,
                           HitCandidateVec&) const;

    // Finding the nearest maximum/minimum from current point
    Waveform::const_iterator findNearestMax(Waveform::const_iterator, Waveform::const_iterator) const;
    Waveform::const_iterator findNearestMin(Waveform::const_iterator, Waveform::const_iterator) const;
    
    // handle finding the "start" and "stop" of a candidate hit
    Waveform::const_iterator findStartTick(Waveform::const_iterator, Waveform::const_iterator) const;
    Waveform::const_iterator findStopTick(Waveform::const_iterator, Waveform::const_iterator)  const;
    
    // some control variables
    std::vector<int>     fMinDeltaTicks;       //< minimum ticks from max to min to consider
    std::vector<int>     fMaxDeltaTicks;       //< maximum ticks from max to min to consider
    std::vector<float>   fMinDeltaPeaks;       //< minimum maximum to minimum peak distance
    float                fMinHitHeight;        //< Drop candidate hits with height less than this
    size_t               fNumInterveningTicks; //< Number ticks between candidate hits to merge
    int                  fStructuringElement;  //< Window size for morphologcial filter
    
    // Member variables from the fhicl file
    std::unique_ptr<reco_tool::IWaveformTool> fWaveformTool;

    // Pointers to the histograms we'll create.
    std::vector<TH1D*>                  fTruncMeanHist;
    std::vector<TH1D*>                  fTruncRmsHist;
    std::vector<TH1D*>                  fFullRmsHist;
    std::vector<TH1D*>                  fNumTruncHist;
    
    std::vector<TH1D*>                  fDeltaTicksHist;
    std::vector<TH1D*>                  fRangeHist;

    art::TFileDirectory*                fHistDirectory;

    // Useful services, keep copies for now (we can update during begin run periods)
    const geo::GeometryCore&                fGeometry;             ///< pointer to Geometry service
    icarusutil::SignalShapingICARUSService& fSignalServices;       ///< The signal shaping service
    const detinfo::DetectorProperties*      fDetectorProperties;   ///< Detector properties service
    const lariov::DetPedestalProvider&      fPedestalRetrievalAlg; ///< Keep track of an instance to the pedestal retrieval alg
};
    
//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// pset - Fcl parameters.
///
BasicWireAnalysis::BasicWireAnalysis(fhicl::ParameterSet const & pset) :
    fGeometry(*lar::providerFrom<geo::Geometry>()),
    fSignalServices(*art::ServiceHandle<icarusutil::SignalShapingICARUSService>()),
    fPedestalRetrievalAlg(*lar::providerFrom<lariov::DetPedestalService>())
{
    fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
    
    configure(pset);
    
    // Report.
    mf::LogInfo("BasicWireAnalysis") << "BasicWireAnalysis configured\n";
}

//----------------------------------------------------------------------------
/// Destructor.
BasicWireAnalysis::~BasicWireAnalysis()
{}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void BasicWireAnalysis::configure(fhicl::ParameterSet const & pset)
{
    fMinDeltaTicks       = pset.get<std::vector<int>   >("MinDeltaTicks",       std::vector<int>()   = {    0,     0,     0});
    fMaxDeltaTicks       = pset.get<std::vector<int>   >("MaxDeltaTicks",       std::vector<int>()   = {   30,    30,    30});
    fMinDeltaPeaks       = pset.get<std::vector<float> >("MinDeltaPeaks",       std::vector<float>() = {0.025, 0.025, 0.025});
    fMinHitHeight        = pset.get<            float  >("MinHitHeight",        2.0);
    fNumInterveningTicks = pset.get<            size_t >("NumInterveningTicks", 6 );
    fStructuringElement  = pset.get<            int    >("StructuringElement",  20);
    
    // Recover an instance of the waveform tool
    fhicl::ParameterSet waveformToolParams;
    
    waveformToolParams.put<std::string>("tool_type","WaveformTools");
    
    fWaveformTool = art::make_tool<reco_tool::IWaveformTool>(waveformToolParams);
}

//----------------------------------------------------------------------------
/// Begin job method.
void BasicWireAnalysis::initializeHists(art::ServiceHandle<art::TFileService>& tfs, const std::string& dirName)
{
    fHistDirectory = tfs.get();
    
    // Make a directory for these histograms
    art::TFileDirectory dir = tfs->mkdir(dirName.c_str());
    
    // Define the histograms. Putting semi-colons around the title
    // causes it to be displayed as the x-axis label if the histogram
    // is drawn.
    
    fTruncMeanHist.resize(3);
    fTruncRmsHist.resize(3);
    fFullRmsHist.resize(3);
    fNumTruncHist.resize(3);
    fDeltaTicksHist.resize(3);
    fRangeHist.resize(3);

    for(size_t plane = 0; plane < fGeometry.Nplanes(); plane++)
    {
        std::string histName = "TruncMean_" + std::to_string(plane);
        
        fTruncMeanHist[plane] = dir.make<TH1D>(histName.c_str(), ";ADC", 200, -50., 50.);
        
        histName = "TruncRMS_" + std::to_string(plane);
        
        fTruncRmsHist[plane] = dir.make<TH1D>(histName.c_str(), ";ADC", 100, 0., 20.);
        
        histName = "FullRMS_" + std::to_string(plane);
        
        fFullRmsHist[plane] = dir.make<TH1D>(histName.c_str(), ";ADC", 100, 0., 20.);
        
        histName = "NTrunc_" + std::to_string(plane);
        
        fNumTruncHist[plane] = dir.make<TH1D>(histName.c_str(), ";Truncated Fraction", 100, 0., 1.);
        
        histName = "DeltaTicks_" + std::to_string(plane);

        fDeltaTicksHist[plane] = dir.make<TH1D>(histName.c_str(), ";Delta", 500, 0., 500.);
        
        histName = "Range_" + std::to_string(plane);
        
        fRangeHist[plane] = dir.make<TH1D>(histName.c_str(), ";Range", 200, 0., 200.);
    }

    return;
}
    
void BasicWireAnalysis::fillHistograms(const IWireHistogramTool::WirePtrVec&     wirePtrVec,
                                       const IWireHistogramTool::SimChannelMap&  channelMap,
                                       int                                       eventNum) const
{
    // Sadly, the RawDigits come to us in an unsorted condition which is not optimal for
    // what we want to do here. So we make a vector of pointers to the input raw digits and sort them
    std::vector<const recob::Wire*> wireVec;
    
    // Ugliness to fill the pointer vector...
    for(size_t idx = 0; idx < wirePtrVec.size(); idx++) wireVec.push_back(wirePtrVec.at(idx).get());
    
    // Sort (use a lambda to sort by channel id)
    std::sort(wireVec.begin(),wireVec.end(),[](const recob::Wire* left, const recob::Wire* right) {return left->Channel() < right->Channel();});

    // Commence looping over raw digits
    for(const auto& wire : wireVec)
    {
        raw::ChannelID_t channel = wire->Channel();

        // Try to limit to the wire number (since we are already segregated by plane)
        std::vector<geo::WireID> wids    = fGeometry.ChannelToWire(channel);
        size_t                   cryo    = wids[0].Cryostat;
        size_t                   tpc     = wids[0].TPC;
        size_t                   plane   = wids[0].Plane;
        size_t                   wireNum = wids[0].Wire;
        
        // Make a directory for these histograms
        art::TFileDirectory dir = fHistDirectory->mkdir(Form("WavePlane_%1zu/c%1zu/c%1zut%1zuwire_%05zu",plane,size_t(eventNum),cryo,tpc,wireNum));

        // If MC, does this channel have signal?
        bool hasSignal = channelMap.find(channel) != channelMap.end();
        
        if (!hasSignal) continue;
        
        const recob::Wire::RegionsOfInterest_t& signalROI = wire->SignalROI();
        
        for(const auto& range : signalROI.get_ranges())
        {
            const Waveform& waveform = range.data();
            
            // Get mean rms and stuff
            float truncMean(0.);
            float truncRMS(0.);
            float fullRMS(0.);
            int   nTrunc(0);
            
            fWaveformTool->getTruncatedMeanRMS(waveform, truncMean, fullRMS, truncRMS, nTrunc);
            
            fTruncMeanHist.at(plane)->Fill(truncMean, 1.);
            fTruncRmsHist.at(plane)->Fill(truncRMS, 1.);
            fFullRmsHist.at(plane)->Fill(fullRMS, 1.);
            fNumTruncHist.at(plane)->Fill(float(nTrunc)/float(waveform.size()),1.);

            // ROI start time
            raw::TDCtick_t roiStartTick = range.begin_index();
            
            HitCandidateVec hitCandidateVec;

            // In this case we want to find hit candidates based on the derivative of of the input waveform
            // We get this from our waveform algs too...
            Waveform rawDerivativeVec;
            Waveform derivativeVec;

            fWaveformTool->firstDerivative(waveform, rawDerivativeVec);
            fWaveformTool->triangleSmooth(rawDerivativeVec,derivativeVec);
            
            // We keep track of the waveform and derivative:
            TProfile* waveHist =
              dir.make<TProfile>(Form("WWfm_%03zu_ctw%01zu-%01zu-%01zu-%05zu-%05zu",size_t(eventNum),cryo,tpc,plane,wireNum,size_t(roiStartTick)), "Waveform", waveform.size(), 0, waveform.size(), -500., 500.);
            TProfile* derivHist =
              dir.make<TProfile>(Form("WDer_%03zu_ctw%01zu-%01zu-%01zu-%05zu-%05zu",size_t(eventNum),cryo,tpc,plane,wireNum,size_t(roiStartTick)), "Derivative", waveform.size(), 0, waveform.size(), -500., 500.);
            TProfile* candHitHist =
              dir.make<TProfile>(Form("WPeak_%03zu_ctw%01zu-%01zu-%01zu-%05zu-%05zu",size_t(eventNum),cryo,tpc,plane,wireNum,size_t(roiStartTick)), "Peaks", waveform.size(), 0, waveform.size(), -500., 500.);
            TProfile* edgeHitHist =
              dir.make<TProfile>(Form("WEdge_%03zu_ctw%01zu-%01zu-%01zu-%05zu-%05zu",size_t(eventNum),cryo,tpc,plane,wireNum,size_t(roiStartTick)), "Edges", waveform.size(), 0, waveform.size(), -500., 500.);

            for(size_t idx = 0; idx < waveform.size(); idx++)
            {
                waveHist->Fill(idx, waveform.at(idx), 1.);
                derivHist->Fill(idx, derivativeVec.at(idx), 1.);
            }

            // Now find the hits
            findHitCandidates(derivativeVec.begin(),derivativeVec.end(),0,plane,hitCandidateVec);
        
            // Reset the hit height from the input waveform
            for(auto& hitCandidate : hitCandidateVec)
            {
                size_t centerIdx = hitCandidate.hitCenter;
                
                candHitHist->Fill(hitCandidate.maxTick, hitCandidate.maxDerivative, 1.);
                candHitHist->Fill(hitCandidate.minTick, hitCandidate.minDerivative, 1.);
                edgeHitHist->Fill(hitCandidate.startTick, 1.);
                edgeHitHist->Fill(hitCandidate.stopTick, 1.);

                hitCandidate.hitHeight = waveform.at(centerIdx);
            }
            
            // We make lots of vectors... erosion, dilation, average and difference
            Waveform erosionVec;
            Waveform dilationVec;
            Waveform averageVec;
            Waveform differenceVec;
            Waveform closingVec;
            Waveform openingVec;
            
            // Define histograms for this particular channel?
            reco_tool::HistogramMap histogramMap;
            
            histogramMap[reco_tool::WAVEFORM] = waveHist;
            histogramMap[reco_tool::EROSION] =
              dir.make<TProfile>(Form("WEro_%03zu_ctw%01zu-%01zu-%01zu-%05zu-%05zu",size_t(eventNum),cryo,tpc,plane,wireNum,size_t(roiStartTick)), "Erosion", waveform.size(), 0, waveform.size(), -500., 500.);
            histogramMap[reco_tool::DILATION] =
              dir.make<TProfile>(Form("WDil_%03zu_ctw%01zu-%01zu-%01zu-%05zu-%05zu",size_t(eventNum),cryo,tpc,plane,wireNum,size_t(roiStartTick)), "Dilation", waveform.size(), 0, waveform.size(), -500., 500.);
            histogramMap[reco_tool::AVERAGE] =
              dir.make<TProfile>(Form("WAve_%03zu_ctw%01zu-%01zu-%01zu-%05zu-%05zu",size_t(eventNum),cryo,tpc,plane,wireNum,size_t(roiStartTick)), "Average", waveform.size(), 0, waveform.size(), -500., 500.);
            histogramMap[reco_tool::DIFFERENCE] =
              dir.make<TProfile>(Form("WDif_%03zu_ctw%01zu-%01zu-%01zu-%05zu-%05zu",size_t(eventNum),cryo,tpc,plane,wireNum,size_t(roiStartTick)), "Difference", waveform.size(), 0, waveform.size(), -500., 500.);
            histogramMap[reco_tool::CLOSING] =
              dir.make<TProfile>(Form("WClo_%03zu_ctw%01zu-%01zu-%01zu-%05zu-%05zu",size_t(eventNum),cryo,tpc,plane,wireNum,size_t(roiStartTick)), "Closing", waveform.size(), 0, waveform.size(), -500., 500.);
            histogramMap[reco_tool::OPENING] =
              dir.make<TProfile>(Form("WOpe_%03zu_ctw%01zu-%01zu-%01zu-%05zu-%05zu",size_t(eventNum),cryo,tpc,plane,wireNum,size_t(roiStartTick)), "Opening", waveform.size(), 0, waveform.size(), -500., 500.);
            histogramMap[reco_tool::DOPENCLOSING] =
                dir.make<TProfile>(Form("WDOC_%03zu_ctw%01zu-%01zu-%01zu-%05zu-%05zu",size_t(eventNum),cryo,tpc,plane,wireNum,size_t(roiStartTick)), "Difference", waveform.size(), 0, waveform.size(), -500., 500.);

            // Compute the morphological filter vectors
            fWaveformTool->getErosionDilationAverageDifference(waveform, fStructuringElement, histogramMap, erosionVec, dilationVec, averageVec, differenceVec);
            
            fWaveformTool->getOpeningAndClosing(erosionVec, dilationVec, fStructuringElement, histogramMap, closingVec, openingVec);

            // Initialial hit finding
            HitCandidateVec hitCandidateMorphVec;

            // Now find the hits
//            findHitCandidates(erosionVec.begin(),erosionVec.end(),differenceVec.begin(),differenceVec.end(),0,plane,hitCandidateMorphVec);
            findHitCandidates(openingVec.begin(),openingVec.end(),differenceVec.begin(),differenceVec.end(),0,plane,hitCandidateMorphVec);

            TProfile* candMorphHist =
            dir.make<TProfile>(Form("MPeak_%03zu_ctw%01zu-%01zu-%01zu-%05zu-%05zu",size_t(eventNum),cryo,tpc,plane,wireNum,size_t(roiStartTick)), "Peaks", waveform.size(), 0, waveform.size(), -500., 500.);
            TProfile* edgeMorphHist =
            dir.make<TProfile>(Form("MEdge_%03zu_ctw%01zu-%01zu-%01zu-%05zu-%05zu",size_t(eventNum),cryo,tpc,plane,wireNum,size_t(roiStartTick)), "Edges", waveform.size(), 0, waveform.size(), -500., 500.);
            
            // Reset the hit height from the input waveform
            for(auto& hitCandidate : hitCandidateMorphVec)
            {
                size_t centerIdx = hitCandidate.hitCenter;
                
                candMorphHist->Fill(hitCandidate.maxTick, hitCandidate.maxDerivative, 1.);
                candMorphHist->Fill(hitCandidate.minTick, hitCandidate.minDerivative, 1.);
                candMorphHist->Fill(hitCandidate.hitCenter, waveform.at(hitCandidate.hitCenter), 1.);
                edgeMorphHist->Fill(hitCandidate.startTick, 1.);
                edgeMorphHist->Fill(hitCandidate.stopTick, 1.);
                
                hitCandidate.hitHeight = waveform.at(centerIdx);
            }
        }
    }
    
    return;
}
    
void BasicWireAnalysis::findHitCandidates(Waveform::const_iterator startItr,
                                          Waveform::const_iterator stopItr,
                                          size_t                   roiStartTick,
                                          size_t                   planeIdx,
                                          HitCandidateVec&         hitCandidateVec) const
{
    // Search for candidate hits...
    // The idea will be to find the largest deviation in the input derivative waveform as the starting point. Depending
    // on if a maximum or minimum, we search forward or backward to find the minimum or maximum that our extremum
    // corresponds to.
    std::pair<Waveform::const_iterator, Waveform::const_iterator> minMaxPair = std::minmax_element(startItr, stopItr);
    
    Waveform::const_iterator maxItr = minMaxPair.second;
    Waveform::const_iterator minItr = minMaxPair.first;
    
    // Use the larger of the two as the starting point and recover the nearest max or min
    if (std::fabs(*maxItr) > std::fabs(*minItr)) minItr = findNearestMin(maxItr, stopItr);
    else                                         maxItr = findNearestMax(minItr,startItr);
    
    int   deltaTicks = std::distance(maxItr,minItr);
    float range      = *maxItr - *minItr;
    
    fDeltaTicksHist.at(planeIdx)->Fill(deltaTicks, 1.);
    fRangeHist.at(planeIdx)->Fill(range, 1.);
    
    //    std::cout << "** max at tick: " << std::distance(startItr,maxItr) << ", val: " << *maxItr << ", min at tick: " << std::distance(startItr,minItr) << ", val: " << *minItr << ", delta: " << deltaTicks << ", range: " << range << std::endl;
    
    // At some point small rolling oscillations on the waveform need to be ignored...
    if (deltaTicks >= fMinDeltaTicks.at(planeIdx) && range > fMinDeltaPeaks.at(planeIdx))
    {
        // Need to back up to find zero crossing, this will be the starting point of our
        // candidate hit but also the endpoint of the pre sub-waveform we'll search next
        Waveform::const_iterator newEndItr = findStartTick(maxItr, startItr);
        
        int startTick = std::distance(startItr,newEndItr);
        
        // Now need to go forward to again get close to zero, this will then be the end point
        // of our candidate hit and the starting point for the post sub-waveform to search
        Waveform::const_iterator newStartItr = findStopTick(minItr, stopItr);
        
        int stopTick = std::distance(startItr,newStartItr);
        
        // Find hits in the section of the waveform leading up to this candidate hit
        if (startTick > 2) findHitCandidates(startItr,newEndItr,roiStartTick,planeIdx,hitCandidateVec);
        
        // Create a new hit candidate and store away
        HitCandidate_t hitCandidate;
        
        Waveform::const_iterator peakItr = std::min_element(maxItr,minItr,[](const auto& left, const auto& right){return std::fabs(left) < std::fabs(right);});
        
        // Check balance
        if      (2 * std::distance(peakItr,minItr) < std::distance(maxItr,peakItr)) peakItr--;
        else if (2 * std::distance(maxItr,peakItr) < std::distance(peakItr,minItr)) peakItr++;
        
        hitCandidate.startTick     = roiStartTick + startTick;
        hitCandidate.stopTick      = roiStartTick + stopTick;
        hitCandidate.maxTick       = roiStartTick + std::distance(startItr,maxItr);
        hitCandidate.minTick       = roiStartTick + std::distance(startItr,minItr);
        hitCandidate.maxDerivative = maxItr != stopItr ? *maxItr : 0.;
        hitCandidate.minDerivative = minItr != stopItr ? *minItr : 0.;
        hitCandidate.hitCenter     = roiStartTick + std::distance(startItr,peakItr) + 0.5;
        hitCandidate.hitSigma      = 0.5 * float(hitCandidate.minTick - hitCandidate.maxTick);
        hitCandidate.hitHeight     = hitCandidate.hitSigma * (hitCandidate.maxDerivative - hitCandidate.minDerivative) / 1.2130;
        
        hitCandidateVec.push_back(hitCandidate);
        
        // Finally, search the section of the waveform following this candidate for more hits
        if (std::distance(newStartItr,stopItr) > 2) findHitCandidates(newStartItr,stopItr,roiStartTick + stopTick,planeIdx,hitCandidateVec);
    }
    
    return;
}
    
void BasicWireAnalysis::findHitCandidates(Waveform::const_iterator eroStartItr,
                                          Waveform::const_iterator eroStopItr,
                                          Waveform::const_iterator diffStartItr,
                                          Waveform::const_iterator diffStopItr,
                                          size_t                   roiStartTick,
                                          size_t                   planeIdx,
                                          HitCandidateVec&         hitCandidateVec) const
{
    // Search for candidate hits...
    // First task is to recover the maximum and minimum difference and reject waveform if no chance for an actual hit
    std::pair<Waveform::const_iterator,Waveform::const_iterator> diffMinMaxPair = std::minmax_element(diffStartItr,diffStopItr);
    
    float deltaDiff = *diffMinMaxPair.second - *diffMinMaxPair.first;
    
    if (deltaDiff < 7) return;
    
    // The idea will be to find the largest deviation in the input derivative waveform as the starting point. Depending
    // on if a maximum or minimum, we search forward or backward to find the minimum or maximum that our extremum
    // corresponds to.
    Waveform::const_iterator peakLeftItr  = std::max_element(eroStartItr,eroStopItr);
    Waveform::const_iterator peakRightItr = peakLeftItr;
    Waveform::const_iterator diffLeftItr  = diffStartItr;

    std::advance(diffLeftItr, std::distance(eroStartItr,peakLeftItr));
    
    Waveform::const_iterator diffRightItr  = diffLeftItr;
    
    // Look for the case of a large waveform, handle differently
    if (*diffLeftItr < *peakLeftItr)
    {
        // Find the position where the erosion vector is clearly less than the difference
        // should this be a ratio test?
        while(std::distance(eroStartItr,peakLeftItr) > 0 && *diffLeftItr  - 1. < *peakLeftItr)  {diffLeftItr--;  peakLeftItr--;}
        while(std::distance(peakRightItr,eroStopItr) > 0 && *diffRightItr - 1. < *peakRightItr) {diffRightItr++; peakRightItr++;}
    }
    
    // Find maximum in difference to left of peak
    Waveform::const_iterator leftMaxItr = diffLeftItr;
    float                    leftMaxVal = *leftMaxItr;
    
    while(std::distance(diffStartItr,leftMaxItr) > 0 && *(leftMaxItr - 1) >= leftMaxVal) leftMaxVal = *(--leftMaxItr);
    
    // Find the rise point in the difference distribution
    Waveform::const_iterator leftRiseItr = leftMaxItr;
    float                    leftRiseVal = *leftRiseItr;
    
    while(std::distance(diffStartItr,leftRiseItr) >= 0 && *(leftRiseItr - 1) < leftRiseVal) leftRiseVal = *(--leftRiseItr);

    // Switch gears and look to the right of the peak
    Waveform::const_iterator rightMaxItr = diffRightItr;
    float                    rightMaxVal = *rightMaxItr;

    while(std::distance(rightMaxItr,diffStopItr) > 0 && *(rightMaxItr + 1) >= rightMaxVal) rightMaxVal = *(++rightMaxItr);
    
    // Find the rise point in the difference distribution
    Waveform::const_iterator rightRiseItr = rightMaxItr;
    float                    rightRiseVal = *rightRiseItr;
    
    while(std::distance(rightRiseItr,diffStopItr) > 0 && *(rightRiseItr + 1) < rightRiseVal) rightRiseVal = *(++rightRiseItr);
    
    // Find hits in the section of the waveform leading up to this candidate hit
    if (std::distance(diffStartItr,leftRiseItr) > 4)
    {
        Waveform::const_iterator newEroStopItr = eroStartItr;
        
        std::advance(newEroStopItr,std::distance(diffStartItr,leftRiseItr));
        
        findHitCandidates(eroStartItr,newEroStopItr,diffStartItr,leftRiseItr,roiStartTick,planeIdx,hitCandidateVec);
    }

    // Fill the data structure
    HitCandidate_t hitCandidate;
    
    hitCandidate.startTick     = roiStartTick + std::distance(diffStartItr,leftRiseItr);
    hitCandidate.stopTick      = roiStartTick + std::distance(diffStartItr,rightRiseItr);
    hitCandidate.maxTick       = roiStartTick + std::distance(diffStartItr,leftMaxItr);
    hitCandidate.minTick       = roiStartTick + std::distance(diffStartItr,rightMaxItr);
    hitCandidate.maxDerivative = *leftMaxItr;
    hitCandidate.minDerivative = *rightMaxItr;
    hitCandidate.hitCenter     = roiStartTick + (std::distance(eroStartItr,peakLeftItr) + std::distance(eroStartItr,peakRightItr) + 0.5)/2;
    hitCandidate.hitSigma      = 0.5 * float(hitCandidate.minTick - hitCandidate.maxTick);
    hitCandidate.hitHeight     = hitCandidate.hitSigma * (hitCandidate.maxDerivative - hitCandidate.minDerivative) / 1.2130;
    
    hitCandidateVec.push_back(hitCandidate);
    
    // Finally, search the section of the waveform following this candidate for more hits
    if (std::distance(rightRiseItr,diffStopItr) > 4)
    {
        Waveform::const_iterator newEroStartItr = eroStartItr;
        int                      newStartTick   = std::distance(diffStartItr,rightRiseItr);
        
        std::advance(newEroStartItr, newStartTick);
        
        findHitCandidates(newEroStartItr,eroStopItr,rightRiseItr,diffStopItr,roiStartTick + newStartTick,planeIdx,hitCandidateVec);
    }

    return;
}

// Useful for normalizing histograms
void BasicWireAnalysis::endJob(int numEvents)
{
    // A task to complete is to fit the average power displays with aim to develop a "good" filter function and
    // get the signal to noise ratio
    
    return;
}
    
BasicWireAnalysis::Waveform::const_iterator BasicWireAnalysis::findNearestMin(Waveform::const_iterator maxItr, Waveform::const_iterator stopItr) const
{
    // reset the min iterator and search forward to find the nearest minimum
    Waveform::const_iterator lastItr = maxItr;
    
    // The strategy is simple, loop forward over ticks until we find the point where the waveform is increasing again
    while((lastItr + 1) != stopItr)
    {
        if (*(lastItr + 1) > *lastItr) break;
        
        lastItr++;
    }
    
    // The minimum will be the last iterator value...
    return lastItr;
}

BasicWireAnalysis::Waveform::const_iterator BasicWireAnalysis::findNearestMax(Waveform::const_iterator minItr, Waveform::const_iterator startItr) const
{
    // Set the internal loop variable...
    Waveform::const_iterator lastItr = minItr;
    
    // One extra condition to watch for here, make sure we can actually back up!
    if (std::distance(startItr,minItr) > 0)
    {
        // Similar to searching for a maximum, we loop backward over ticks looking for the waveform to start decreasing
        while((lastItr - 1) != startItr)
        {
            if (*(lastItr - 1) < *lastItr) break;
            
            lastItr--;
        }
    }
    
    return lastItr;
}

BasicWireAnalysis::Waveform::const_iterator BasicWireAnalysis::findStartTick(Waveform::const_iterator maxItr, Waveform::const_iterator startItr) const
{
    Waveform::const_iterator lastItr = maxItr;
    
    // If we can't back up then there is nothing to do
    if (std::distance(startItr,lastItr) > 0)
    {
        // In theory, we are starting at a maximum and want to find the "start" of the candidate peak
        // Ideally we would look to search backward to the point where the (derivative) waveform crosses zero again.
        // However, the complication is that we need to watch for the case where two peaks are merged together and
        // we might run through another peak before crossing zero...
        // So... loop until we hit the startItr...
        Waveform::const_iterator loopItr = lastItr - 1;
        
        while(loopItr != startItr)
        {
            // Ideal world case, we cross zero... but we might encounter a minimum... or an inflection point
            if (*loopItr < 0. || !(*loopItr < *lastItr)) break;
            
            lastItr = loopItr--;
        }
    }
    else lastItr = startItr;
    
    return lastItr;
}

BasicWireAnalysis::Waveform::const_iterator BasicWireAnalysis::findStopTick(Waveform::const_iterator minItr, Waveform::const_iterator stopItr)   const
{
    Waveform::const_iterator lastItr = minItr;
    
    // If we can't go forward then there is really nothing to do
    if (std::distance(minItr,stopItr) > 1)
    {
        // Pretty much the same strategy as for finding the start tick...
        Waveform::const_iterator loopItr = lastItr + 1;
        
        while(loopItr != stopItr)
        {
            // Ideal case that we have crossed zero coming from a minimum... but watch for a maximum as well
            if (*loopItr > 0. || !(*loopItr > *lastItr)) break;
            
            lastItr = loopItr++;
        }
    }
    
    return lastItr;
}
    
DEFINE_ART_CLASS_TOOL(BasicWireAnalysis)
}
