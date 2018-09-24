
#include "icaruscode/Analysis/tools/IHitEfficiencyHistogramTool.h"

#include "fhiclcpp/ParameterSet.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "lardataobj/RecoBase/Hit.h"

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TProfile2D.h"

#include <cmath>
#include <algorithm>

namespace HitEfficiencyAnalysis
{
    ////////////////////////////////////////////////////////////////////////
    //
    // Class:       HitEfficiencyAnalysis
    // Module Type: producer
    // File:        HitEfficiencyAnalysis.h
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
    
// The following typedefs will, obviously, be useful
using HitPtrVec       = std::vector<art::Ptr<recob::Hit>>;
using ViewHitMap      = std::map<size_t,HitPtrVec>;
using TrackViewHitMap = std::map<int,ViewHitMap>;

class HitEfficiencyAnalysis : virtual public IHitEfficiencyHistogramTool
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    explicit HitEfficiencyAnalysis(fhicl::ParameterSet const & pset);
    
    /**
     *  @brief  Destructor
     */
    ~HitEfficiencyAnalysis();
    
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
    void fillHistograms(const std::vector<recob::Hit>&, const std::vector<sim::SimChannel>&)  const override;
    
private:
    
    // Fcl parameters.
    std::string                 fLocalDirName;     ///< Fraction for truncated mean
    std::vector<unsigned short> fOffsetVec;        ///< Allow offsets for each plane
    
    // Pointers to the histograms we'll create.
    std::vector<TH1F*>     fTotalElectronsVec;
    std::vector<TH1F*>     fMaxElectronsVec;
    std::vector<TH1F*>     fHitElectronsVec;
    std::vector<TH1F*>     fHitSumADCVec;
    std::vector<TH1F*>     fHitPulseHeightVec;
    std::vector<TH1F*>     fHitPulseWidthVec;
    std::vector<TH1F*>     fSimNumTDCVec;
    std::vector<TH1F*>     fHitNumTDCVec;
    std::vector<TH1F*>     fNMatchedHitVec;
    std::vector<TH1F*>     fDeltaMidTDCVec;
    std::vector<TProfile*> fHitEfficVec;
    std::vector<TProfile*> fHitEfficPHVec;
    std::vector<TH2F*>     fHitVsSimChgVec;
    
    std::vector<TH1F*>     fNSimChannelHitsVec;
    std::vector<TH1F*>     fNRecobHitVec;
    std::vector<TH1F*>     fHitEfficiencyVec;
    
    // Useful services, keep copies for now (we can update during begin run periods)
    const geo::GeometryCore*           fGeometry;             ///< pointer to Geometry service
    const detinfo::DetectorProperties* fDetectorProperties;   ///< Detector properties service
    const detinfo::DetectorClocks*     fClockService;         ///< Detector clocks service
};
    
//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// pset - Fcl parameters.
///
HitEfficiencyAnalysis::HitEfficiencyAnalysis(fhicl::ParameterSet const & pset)
{
    fGeometry           = lar::providerFrom<geo::Geometry>();
    fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
    fClockService       = lar::providerFrom<detinfo::DetectorClocksService>();
    
    configure(pset);
    
    // Report.
    mf::LogInfo("HitEfficiencyAnalysis") << "HitEfficiencyAnalysis configured\n";
}

//----------------------------------------------------------------------------
/// Destructor.
HitEfficiencyAnalysis::~HitEfficiencyAnalysis()
{}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void HitEfficiencyAnalysis::configure(fhicl::ParameterSet const & pset)
{
    fLocalDirName = pset.get<std::string>("LocalDirName", std::string("wow"));
    fOffsetVec    = pset.get<std::vector<unsigned short>>("OffsetVec", std::vector<unsigned short>()={0,0,0});
}

//----------------------------------------------------------------------------
/// Begin job method.
void HitEfficiencyAnalysis::initializeHists(art::ServiceHandle<art::TFileService>& tfs, const std::string& dirName)
{
    // Make a directory for these histograms
    art::TFileDirectory dir = tfs->mkdir(dirName.c_str());
    
    fTotalElectronsVec.resize(fGeometry->Nplanes());
    
    fTotalElectronsVec[0] = dir.make<TH1F>("TotalElecs0", ";# electrons", 250, 0., 100000.);
    fTotalElectronsVec[1] = dir.make<TH1F>("TotalElecs1", ";# electrons", 250, 0., 100000.);
    fTotalElectronsVec[2] = dir.make<TH1F>("TotalElecs2", ";# electrons", 250, 0., 100000.);
    
    fMaxElectronsVec.resize(fGeometry->Nplanes());
    
    fMaxElectronsVec[0]   = dir.make<TH1F>("MaxElecs0", ";# electrons", 250, 0., 20000.);
    fMaxElectronsVec[1]   = dir.make<TH1F>("MaxElecs1", ";# electrons", 250, 0., 20000.);
    fMaxElectronsVec[2]   = dir.make<TH1F>("MaxElecs2", ";# electrons", 250, 0., 20000.);

    fHitElectronsVec.resize(fGeometry->Nplanes());
    
    fHitElectronsVec[0]   = dir.make<TH1F>("HitElecs0", ";# electrons", 250, 0., 100000.);
    fHitElectronsVec[1]   = dir.make<TH1F>("HitElecs1", ";# electrons", 250, 0., 100000.);
    fHitElectronsVec[2]   = dir.make<TH1F>("HitElecs2", ";# electrons", 250, 0., 100000.);
    
    fHitSumADCVec.resize(fGeometry->Nplanes());
    
    fHitSumADCVec[0]      = dir.make<TH1F>("SumADC0",   "Sum ADC",    500,  0., 5000.);
    fHitSumADCVec[1]      = dir.make<TH1F>("SumADC1",   "Sum ADC",    500,  0., 5000.);
    fHitSumADCVec[2]      = dir.make<TH1F>("SumADC2",   "Sum ADC",    500,  0., 5000.);
    
    fHitPulseHeightVec.resize(fGeometry->Nplanes());
    
    fHitPulseHeightVec[0] = dir.make<TH1F>("PulseHeight0", "PH (ADC)", 150,  0.,  150.);
    fHitPulseHeightVec[1] = dir.make<TH1F>("PulseHeight1", "PH (ADC)", 150,  0.,  150.);
    fHitPulseHeightVec[2] = dir.make<TH1F>("PulseHeight2", "PH (ADC)", 150,  0.,  150.);

    fHitPulseWidthVec.resize(fGeometry->Nplanes());
    
    fHitPulseWidthVec[0]  = dir.make<TH1F>("PulseWidth0",  ";RMS",     40,  0.,  20.);
    fHitPulseWidthVec[1]  = dir.make<TH1F>("PulseWidth1",  ";RMS",     40,  0.,  20.);
    fHitPulseWidthVec[2]  = dir.make<TH1F>("PulseWidth2",  ";RMS",     40,  0.,  20.);

    fSimNumTDCVec.resize(fGeometry->Nplanes());
    
    fSimNumTDCVec[0]      = dir.make<TH1F>("SimNumTDC0", ";TDC ticks", 100, 0., 100.);
    fSimNumTDCVec[1]      = dir.make<TH1F>("SimNumTDC1", ";TDC ticks", 100, 0., 100.);
    fSimNumTDCVec[2]      = dir.make<TH1F>("SimNumTDC2", ";TDC ticks", 100, 0., 100.);

    fHitNumTDCVec.resize(fGeometry->Nplanes());
    
    fHitNumTDCVec[0]      = dir.make<TH1F>("HitNumTDC0", ";TDC ticks", 100, 0., 100.);
    fHitNumTDCVec[1]      = dir.make<TH1F>("HitNumTDC1", ";TDC ticks", 100, 0., 100.);
    fHitNumTDCVec[2]      = dir.make<TH1F>("HitNumTDC2", ";TDC ticks", 100, 0., 100.);
    
    fNMatchedHitVec.resize(fGeometry->Nplanes());
    
    fNMatchedHitVec[0]    = dir.make<TH1F>("NMatched0", ";# hits", 20, 0., 20.);
    fNMatchedHitVec[1]    = dir.make<TH1F>("NMatched1", ";# hits", 20, 0., 20.);
    fNMatchedHitVec[2]    = dir.make<TH1F>("NMatched2", ";# hits", 20, 0., 20.);

    fDeltaMidTDCVec.resize(fGeometry->Nplanes());
    
    fDeltaMidTDCVec[0]    = dir.make<TH1F>("DeltaMid0", ";# hits", 50, -25., 25.);
    fDeltaMidTDCVec[1]    = dir.make<TH1F>("DeltaMid1", ";# hits", 50, -25., 25.);
    fDeltaMidTDCVec[2]    = dir.make<TH1F>("DeltaMid2", ";# hits", 50, -25., 25.);

    fHitEfficVec.resize(fGeometry->Nplanes());
    
    fHitEfficVec[0]       = dir.make<TProfile>("HitEffic0", "Hit Efficiency;# electrons", 200, 0., 100000., 0., 1.);
    fHitEfficVec[1]       = dir.make<TProfile>("HitEffic1", "Hit Efficiency;# electrons", 200, 0., 100000., 0., 1.);
    fHitEfficVec[2]       = dir.make<TProfile>("HitEffic2", "Hit Efficiency;# electrons", 200, 0., 100000., 0., 1.);

    fHitEfficPHVec.resize(fGeometry->Nplanes());
    
    fHitEfficPHVec[0]     = dir.make<TProfile>("HitEfficPH0", "Hit Efficiency;# electrons", 200, 0., 20000., 0., 1.);
    fHitEfficPHVec[1]     = dir.make<TProfile>("HitEfficPH1", "Hit Efficiency;# electrons", 200, 0., 20000., 0., 1.);
    fHitEfficPHVec[2]     = dir.make<TProfile>("HitEfficPH2", "Hit Efficiency;# electrons", 200, 0., 20000., 0., 1.);

    fHitVsSimChgVec.resize(fGeometry->Nplanes());
    
    fHitVsSimChgVec[0]    = dir.make<TH2F>("HitVSimQ0", "Sim;Hit", 250, 0., 5000., 250, 0., 100000.);
    fHitVsSimChgVec[1]    = dir.make<TH2F>("HitVSimQ1", "Sim;Hit", 250, 0., 5000., 250, 0., 100000.);
    fHitVsSimChgVec[2]    = dir.make<TH2F>("HitVSimQ2", "Sim;Hit", 250, 0., 5000., 250, 0., 100000.);
    
    fNSimChannelHitsVec.resize(fGeometry->Nplanes());
    
    fNSimChannelHitsVec[0] = dir.make<TH1F>("NSimChan0", ";# hits", 300, 0., 1200.);
    fNSimChannelHitsVec[1] = dir.make<TH1F>("NSimChan1", ";# hits", 500, 0., 2000.);
    fNSimChannelHitsVec[2] = dir.make<TH1F>("NSimChan2", ";# hits", 500, 0., 2000.);

    fNRecobHitVec.resize(fGeometry->Nplanes());
    
    fNRecobHitVec[0] = dir.make<TH1F>("NRecobHit0", ";# hits", 300, 0., 1200.);
    fNRecobHitVec[1] = dir.make<TH1F>("NRecobHit1", ";# hits", 500, 0., 2000.);
    fNRecobHitVec[2] = dir.make<TH1F>("NRecobHit2", ";# hits", 500, 0., 2000.);

    fHitEfficiencyVec.resize(fGeometry->Nplanes());
    
    fHitEfficiencyVec[0] = dir.make<TH1F>("PlnEffic0", ";# hits", 101, 0., 1.01);
    fHitEfficiencyVec[1] = dir.make<TH1F>("PlnEffic1", ";# hits", 101, 0., 1.01);
    fHitEfficiencyVec[2] = dir.make<TH1F>("PlnEffic2", ";# hits", 101, 0., 1.01);

    return;
}
    
void HitEfficiencyAnalysis::fillHistograms(const std::vector<recob::Hit>& hitVec, const std::vector<sim::SimChannel>& simChannelVec) const
{
    // If there is no sim channel informaton then exit
    if (simChannelVec.empty()) return;
    
    // what needs to be done?
    // First we should map out all hits by channel so we can easily look up from sim channels
    // Then go through the sim channels and match hits
    using ChanToHitVecMap = std::map<raw::ChannelID_t,std::vector<const recob::Hit*>>;
    ChanToHitVecMap channelToHitVec;
    
    for(const auto& hit : hitVec) channelToHitVec[hit.Channel()].push_back(&hit);
    
    // Now go through the sim channels
    // There are several things going on here... for each channel we have particles (track id's) depositing energy in a range to ticks
    // So... for each channel we want to build a structure that relates particles to tdc ranges and deposited energy (or electrons)
    // Here is a complicated structure:
    
    using TDCToIDEMap             = std::map<unsigned short, float>;
    using ChanToTDCToIDEMap       = std::map<raw::ChannelID_t, TDCToIDEMap>;
    using PartToChanToTDCToIDEMap = std::map<int, ChanToTDCToIDEMap>;
    
    PartToChanToTDCToIDEMap partToChanToTDCToIDEMap;
    
    // Build out the above data structure
    for(const auto& simChannel : simChannelVec)
    {
        for(const auto& tdcide : simChannel.TDCIDEMap())
        {
            for(const auto& ide : tdcide.second) partToChanToTDCToIDEMap[ide.trackID][simChannel.Channel()][tdcide.first] = ide.numElectrons;
        }
    }

    // Find the longest which is meant to be the primary
    PartToChanToTDCToIDEMap::iterator longestItr = partToChanToTDCToIDEMap.begin();
    
    for(PartToChanToTDCToIDEMap::iterator chanItr = partToChanToTDCToIDEMap.begin(); chanItr != partToChanToTDCToIDEMap.end(); chanItr++)
    {
        if (chanItr->second.size() > longestItr->second.size()) longestItr = chanItr;
    }
    
    std::vector<int> nSimChannelHitVec = {0,0,0};
    std::vector<int> nRecobHitVec      = {0,0,0};
    
    // Go through the longest iterator and match to hits
    for(const auto& chanToTDCToIDEMap : longestItr->second)
    {
        ChanToHitVecMap::iterator hitIter     = channelToHitVec.find(chanToTDCToIDEMap.first);
        TDCToIDEMap               tdcToIDEMap = chanToTDCToIDEMap.second;
        float                     totalElectrons(0.);
        float                     maxElectrons(0.);
        int                       nMatchedHits(0);
        
        // The below try-catch block may no longer be necessary
        // Decode the channel and make sure we have a valid one
        std::vector<geo::WireID> wids = fGeometry->ChannelToWire(chanToTDCToIDEMap.first);
        
        // Recover plane and wire in the plane
        unsigned int plane = wids[0].Plane;
//        unsigned int wire  = wids[0].Wire;
        
        for(const auto& ideVal : tdcToIDEMap)
        {
            totalElectrons += ideVal.second;
            
            maxElectrons = std::max(maxElectrons,ideVal.second);
        }
        
        totalElectrons = std::min(totalElectrons, float(99900.));
        
        fTotalElectronsVec.at(plane)->Fill(totalElectrons, 1.);
        fMaxElectronsVec.at(plane)->Fill(maxElectrons, 1.);
        
        nSimChannelHitVec.at(plane)++;

        if (hitIter != channelToHitVec.end())
        {
            unsigned short startTDC = tdcToIDEMap.begin()->first;
            unsigned short stopTDC  = tdcToIDEMap.rbegin()->first;
            unsigned short midTDC   = (startTDC + stopTDC) / 2;
            
            fSimNumTDCVec.at(plane)->Fill(stopTDC - startTDC, 1.);

            // Set up to extract the "best" parameters in the event of more than one hit for this pulse train
            float          nElectronsTotalBest(0.);
            float          hitChargeBest(0.);
            float          hitPulseHeightBest(0.);
            float          hitWidthBest(0.);
            unsigned short hitStopTDCBest(0);
            unsigned short hitStartTDCBest(0);
            unsigned short midHitTDCBest(0);

            // Loop through the hits for this channel and look for matches
            // In the event of more than one hit associated to the sim channel range, keep only
            // the best match (assuming the nearby hits are "extra")
            // Note that assumption breaks down for long pulse trains but worry about that later
            for(const auto& hit : hitIter->second)
            {
                unsigned short hitStartTick = hit->PeakTime() - 1. * hit->RMS();
                unsigned short hitStopTick  = hit->PeakTime() + 1. * hit->RMS();
                unsigned short hitStartTDC  = fClockService->TPCTick2TDC(hitStartTick) - fOffsetVec.at(plane);
                unsigned short hitStopTDC   = fClockService->TPCTick2TDC(hitStopTick)  - fOffsetVec.at(plane);
                unsigned short midHitTDC    = (hitStopTDC + hitStartTDC) / 2;
                
                // If hit is out of range then skip, it is not related to this particle
                if (hitStartTDC > stopTDC || hitStopTDC < startTDC) continue;
                
                float hitHeight = std::min(hit->PeakAmplitude(), float(149.5));
                
                // Use the hit with the largest pulse height as the "best"
                if (hitHeight < hitPulseHeightBest) continue;

                nElectronsTotalBest = 0.;
                
                hitPulseHeightBest = hitHeight;
                hitChargeBest      = std::min(hit->SummedADC(),float(4999.));
                hitWidthBest       = std::min(hit->RMS(), float(19.8));
                hitStartTDCBest    = hitStartTDC;
                hitStopTDCBest     = hitStopTDC;
                midHitTDCBest      = midHitTDC;
                
                nMatchedHits++;
                
                // Get the number of electrons
                for(unsigned short tick = hitStartTDC; tick <= hitStopTDC; tick++)
                {
                    TDCToIDEMap::iterator ideIterator = tdcToIDEMap.find(tick);
                    
                    if (ideIterator != tdcToIDEMap.end()) nElectronsTotalBest += ideIterator->second;
                }
            }

            if (nMatchedHits > 0)
            {
                fHitSumADCVec.at(plane)->Fill(hitChargeBest, 1.);
                fHitVsSimChgVec.at(plane)->Fill(hitChargeBest, totalElectrons, 1.);
                fHitPulseHeightVec.at(plane)->Fill(hitPulseHeightBest, 1.);
                fHitPulseWidthVec.at(plane)->Fill(hitWidthBest, 1.);
                fHitElectronsVec.at(plane)->Fill(nElectronsTotalBest, 1.);
                fHitNumTDCVec.at(plane)->Fill(hitStopTDCBest - hitStartTDCBest, 1.);
                fDeltaMidTDCVec.at(plane)->Fill(midHitTDCBest - midTDC, 1.);
                
                nRecobHitVec.at(plane)++;
            }
        }
        
        fNMatchedHitVec.at(plane)->Fill(nMatchedHits, 1.);
        fHitEfficVec.at(plane)->Fill(totalElectrons, std::min(nMatchedHits,1),1.);
        fHitEfficPHVec.at(plane)->Fill(maxElectrons, std::min(nMatchedHits,1),1.);
    }
    
    for(size_t idx = 0; idx < fGeometry->Nplanes();idx++)
    {
        if (nSimChannelHitVec.at(idx) > 10)
        {
            float hitEfficiency = float(nRecobHitVec.at(idx)) / float(nSimChannelHitVec.at(idx));
            
            fNSimChannelHitsVec.at(idx)->Fill(std::min(nSimChannelHitVec.at(idx),1999),1.);
            fNRecobHitVec.at(idx)->Fill(std::min(nRecobHitVec.at(idx),1999), 1.);
            fHitEfficiencyVec.at(idx)->Fill(hitEfficiency, 1.);
        }
    }

    return;
}
    
// Useful for normalizing histograms
void HitEfficiencyAnalysis::endJob(int numEvents)
{
    return;
}
    
DEFINE_ART_CLASS_TOOL(HitEfficiencyAnalysis)
}
