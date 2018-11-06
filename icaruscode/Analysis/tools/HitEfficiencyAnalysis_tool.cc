
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
#include "TTree.h"

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
     *  @brief Interface for initializing the tuple variables
     *
     *  @param TTree          pointer to a TTree object to which to add variables
     */
    void initializeTuple(TTree*) override;

    /**
     *  @brief Interface for method to executve at the end of run processing
     *
     *  @param int            number of events to use for normalization
     */
    void endJob(int numEvents) override;
    
    /**
     *  @brief Interface for filling histograms
     */
    void fillHistograms(const std::vector<recob::Hit>&, const std::vector<simb::MCParticle>&, const std::vector<sim::SimChannel>&)  const override;
    
private:
    
    // Clear mutable variables
    void clear() const;
    
    // Fcl parameters.
    std::string                 fLocalDirName;     ///< Fraction for truncated mean
    std::vector<unsigned short> fOffsetVec;        ///< Allow offsets for each plane
    std::vector<float>          fSigmaVec;         ///< Window size for matching to SimChannels
    
    // Pointers to the histograms we'll create.
    std::vector<TH1F*>          fTotalElectronsHistVec;
    std::vector<TH1F*>          fMaxElectronsHistVec;
    std::vector<TH1F*>          fHitElectronsVec;
    std::vector<TH1F*>          fHitSumADCVec;
    std::vector<TH1F*>          fHitPulseHeightVec;
    std::vector<TH1F*>          fHitPulseWidthVec;
    std::vector<TH1F*>          fSimNumTDCVec;
    std::vector<TH1F*>          fHitNumTDCVec;
    std::vector<TH1F*>          fNMatchedHitVec;
    std::vector<TH1F*>          fDeltaMidTDCVec;
    std::vector<TProfile*>      fHitEfficVec;
    std::vector<TProfile*>      fHitEfficPHVec;
    std::vector<TH2F*>          fHitVsSimChgVec;
    
    std::vector<TH1F*>          fNSimChannelHitsVec;
    std::vector<TH1F*>          fNRecobHitVec;
    std::vector<TH1F*>          fHitEfficiencyVec;
    
    // TTree variables
    mutable TTree*             fTree;
    
    mutable std::vector<int>   fTPCVec;
    mutable std::vector<int>   fCryoVec;
    mutable std::vector<int>   fPlaneVec;
    mutable std::vector<int>   fWireVec;
    
    mutable std::vector<float> fTotalElectronsVec;
    mutable std::vector<float> fMaxElectronsVec;
    mutable std::vector<int>   fStartTickVec;
    mutable std::vector<int>   fStopTickVec;
    mutable int                fNMatchedHits;
    
    mutable std::vector<float> fHitPeakTimeVec;
    mutable std::vector<float> fHitPeakAmpVec;
    mutable std::vector<float> fHitPeakRMSVec;
    mutable std::vector<float> fHitBaselinevec;
    mutable std::vector<float> fHitSummedADCVec;
    mutable std::vector<float> fHitIntegralVec;
    mutable std::vector<int>   fHitStartTickVec;
    mutable std::vector<int>   fHitStopTickVec;
    
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
HitEfficiencyAnalysis::HitEfficiencyAnalysis(fhicl::ParameterSet const & pset) : fTree(nullptr)
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
    fSigmaVec     = pset.get<std::vector<float>>("SigmaVec", std::vector<float>()={1.,1.,1.});
}

//----------------------------------------------------------------------------
/// Begin job method.
void HitEfficiencyAnalysis::initializeHists(art::ServiceHandle<art::TFileService>& tfs, const std::string& dirName)
{
    // Make a directory for these histograms
    art::TFileDirectory dir = tfs->mkdir(dirName.c_str());
    
    fTotalElectronsHistVec.resize(fGeometry->Nplanes());
    fMaxElectronsHistVec.resize(fGeometry->Nplanes());
    fHitElectronsVec.resize(fGeometry->Nplanes());
    fHitSumADCVec.resize(fGeometry->Nplanes());
    fHitPulseHeightVec.resize(fGeometry->Nplanes());
    fHitPulseWidthVec.resize(fGeometry->Nplanes());
    fSimNumTDCVec.resize(fGeometry->Nplanes());
    fHitNumTDCVec.resize(fGeometry->Nplanes());
    fNMatchedHitVec.resize(fGeometry->Nplanes());
    fDeltaMidTDCVec.resize(fGeometry->Nplanes());
    fHitVsSimChgVec.resize(fGeometry->Nplanes());
    fNSimChannelHitsVec.resize(fGeometry->Nplanes());
    fNRecobHitVec.resize(fGeometry->Nplanes());
    fHitEfficiencyVec.resize(fGeometry->Nplanes());

    fHitEfficVec.resize(fGeometry->Nplanes());
    fHitEfficPHVec.resize(fGeometry->Nplanes());

    for(size_t plane = 0; plane < fGeometry->Nplanes(); plane++)
    {
        fTotalElectronsHistVec.at(plane)  = dir.make<TH1F>(("TotalElecs"  + std::to_string(plane)).c_str(), ";# electrons", 250,   0.,  100000.);
        fMaxElectronsHistVec.at(plane)    = dir.make<TH1F>(("MaxElecs"    + std::to_string(plane)).c_str(), ";# electrons", 250,   0.,  20000.);
        fHitElectronsVec.at(plane)        = dir.make<TH1F>(("HitElecs"    + std::to_string(plane)).c_str(), ";# electrons", 250,   0.,  100000.);
        fHitSumADCVec.at(plane)           = dir.make<TH1F>(("SumADC"      + std::to_string(plane)).c_str(), "Sum ADC",      500,   0.,  5000.);
        fHitPulseHeightVec.at(plane)      = dir.make<TH1F>(("PulseHeight" + std::to_string(plane)).c_str(), "PH (ADC)",     150,   0.,  150.);
        fHitPulseWidthVec.at(plane)       = dir.make<TH1F>(("PulseWidth"  + std::to_string(plane)).c_str(), ";RMS",          40,   0.,  20.);
        fSimNumTDCVec.at(plane)           = dir.make<TH1F>(("SimNumTDC"   + std::to_string(plane)).c_str(), ";TDC ticks",   100,   0.,  100.);
        fHitNumTDCVec.at(plane)           = dir.make<TH1F>(("HitNumTDC"   + std::to_string(plane)).c_str(), ";TDC ticks",   100,   0.,  100.);
        fNMatchedHitVec.at(plane)         = dir.make<TH1F>(("NMatched"    + std::to_string(plane)).c_str(), ";# hits",       20,   0.,  20.);
        fDeltaMidTDCVec.at(plane)         = dir.make<TH1F>(("DeltaMid"    + std::to_string(plane)).c_str(), ";# hits",       50, -25.,  25.);
        fNSimChannelHitsVec.at(plane)     = dir.make<TH1F>(("NSimChan"    + std::to_string(plane)).c_str(), ";# hits",      300,   0.,  1200.);
        fNRecobHitVec.at(plane)           = dir.make<TH1F>(("NRecobHit"   + std::to_string(plane)).c_str(), ";# hits",      300,   0.,  1200.);
        fHitEfficiencyVec.at(plane)       = dir.make<TH1F>(("PlnEffic"    + std::to_string(plane)).c_str(), ";# hits",      101,   0.,  1.01);
    
        fHitVsSimChgVec.at(plane)         = dir.make<TH2F>(("HitVSimQ" + std::to_string(plane)).c_str(), "Sim;Hit", 250, 0., 5000., 250, 0., 100000.);
    
        fHitEfficVec.at(plane)            = dir.make<TProfile>(("HitEffic"   + std::to_string(plane)).c_str(), "Hit Efficiency;# electrons", 200, 0., 100000., 0., 1.);
        fHitEfficPHVec.at(plane)          = dir.make<TProfile>(("HitEfficPH" + std::to_string(plane)).c_str(), "Hit Efficiency;# electrons", 200, 0.,  20000., 0., 1.);
    }
    
    return;
}
    
void HitEfficiencyAnalysis::initializeTuple(TTree* tree)
{
    fTree = tree;
    
    fTree->Branch("CryostataVec",      "std::vector<int>",   &fCryoVec);
    fTree->Branch("TPCVec",            "std::vector<int>",   &fTPCVec);
    fTree->Branch("PlaneVec",          "std::vector<int>",   &fPlaneVec);
    fTree->Branch("WireVec",           "std::vector<int>",   &fWireVec);

    fTree->Branch("TotalElectronsVec", "std::vector<float>", &fTotalElectronsVec);
    fTree->Branch("MaxElectronsVec",   "std::vector<float>", &fMaxElectronsVec);
    fTree->Branch("StartTick",         "std::vector<int>",   &fStartTickVec);
    fTree->Branch("StopTick",          "std::vector<int>",   &fStopTickVec);
    fTree->Branch("NMatchedHits",      &fNMatchedHits,       "NMatchedHits/I");
    
    fTree->Branch("HitPeakTimeVec",    "std::vector<float>", &fHitPeakTimeVec);
    fTree->Branch("HitPeakAmpVec",     "std::vector<float>", &fHitPeakAmpVec);
    fTree->Branch("HitPeakRMSVec",     "std::vector<float>", &fHitPeakRMSVec);
    fTree->Branch("HitBaselineVec",    "std::vector<float>", &fHitBaselinevec);
    fTree->Branch("HitSummedADCVec",   "std::vector<float>", &fHitSummedADCVec);
    fTree->Branch("HitIntegralVec",    "std::vector<float>", &fHitIntegralVec);
    fTree->Branch("HitStartTickVec",   "std::vector<float>", &fHitStartTickVec);
    fTree->Branch("HitStopTickVec",    "std::vector<float>", &fHitStopTickVec);
                  
    clear();

    return;
}
    
void HitEfficiencyAnalysis::clear() const
{
    fTPCVec.clear();
    fCryoVec.clear();
    fPlaneVec.clear();
    fWireVec.clear();

    fTotalElectronsVec.clear();
    fMaxElectronsVec.clear();
    fStartTickVec.clear();
    fStopTickVec.clear();
    fNMatchedHits = 0;

    fHitPeakTimeVec.clear();
    fHitPeakAmpVec.clear();
    fHitPeakRMSVec.clear();
    fHitBaselinevec.clear();
    fHitSummedADCVec.clear();
    fHitIntegralVec.clear();
    fHitStartTickVec.clear();
    fHitStopTickVec.clear();

    return;
}
    
void HitEfficiencyAnalysis::fillHistograms(const std::vector<recob::Hit>& hitVec, const std::vector<simb::MCParticle>& mcParticleVec, const std::vector<sim::SimChannel>& simChannelVec) const
{
    // If there is no sim channel informaton then exit
    if (simChannelVec.empty()) return;
    
    // what needs to be done?
    // First we should map out all hits by channel so we can easily look up from sim channels
    // Then go through the sim channels and match hits
    using ChanToHitVecMap = std::map<raw::ChannelID_t,std::vector<const recob::Hit*>>;
    ChanToHitVecMap channelToHitVec;
    
    for(const auto& hit : hitVec) channelToHitVec[hit.Channel()].push_back(&hit);

    // It is useful to create a mapping between trackID and MCParticle
    using TrackIDToMCParticleMap = std::map<int, const simb::MCParticle*>;

    TrackIDToMCParticleMap trackIDToMCParticleMap;

    for(const auto& mcParticle : mcParticleVec) trackIDToMCParticleMap[mcParticle.TrackId()] = &mcParticle;
    
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
    
    std::vector<int> nSimChannelHitVec = {0,0,0};
    std::vector<int> nRecobHitVec      = {0,0,0};

    for(const auto& partToChanInfo : partToChanToTDCToIDEMap)
    {
        TrackIDToMCParticleMap::const_iterator trackIDToMCPartItr = trackIDToMCParticleMap.find(partToChanInfo.first);

        if (trackIDToMCPartItr == trackIDToMCParticleMap.end()) continue;

        int         trackPDGCode = trackIDToMCPartItr->second->PdgCode();
        std::string processName  = trackIDToMCPartItr->second->Process();

        // Looking for primary muons (e.g. CR Tracks)
        if (fabs(trackPDGCode) != 13 || processName != "primary") continue;

        for(const auto& chanToTDCToIDEMap : partToChanInfo.second)
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
//            unsigned int wire  = wids[0].Wire;
            
            for(const auto& ideVal : tdcToIDEMap)
            {
                totalElectrons += ideVal.second;
                
                maxElectrons = std::max(maxElectrons,ideVal.second);
            }
            
            totalElectrons = std::min(totalElectrons, float(99900.));
            
            fTotalElectronsHistVec.at(plane)->Fill(totalElectrons, 1.);
            fMaxElectronsHistVec.at(plane)->Fill(maxElectrons, 1.);
            
            nSimChannelHitVec.at(plane)++;

            unsigned short startTDC = tdcToIDEMap.begin()->first;
            unsigned short stopTDC  = tdcToIDEMap.rbegin()->first;
            
            // Convert to ticks to get in same units as hits
            unsigned short startTick = fClockService->TPCTDC2Tick(startTDC) + fOffsetVec.at(plane);
            unsigned short stopTick  = fClockService->TPCTDC2Tick(stopTDC)  + fOffsetVec.at(plane);
            unsigned short midTick   = (startTick + stopTick) / 2;

            fSimNumTDCVec.at(plane)->Fill(stopTick - startTick, 1.);
    
            // Set up to extract the "best" parameters in the event of more than one hit for this pulse train
            float          nElectronsTotalBest(0.);
            float          hitSummedADCBest(0.);
            float          hitIntegralBest(0.);
            float          hitPeakTimeBest(0.);
            float          hitPeakAmpBest(-100.);
            float          hitRMSBest(0.);
            float          hitBaselineBest(0.);
            unsigned short hitStopTickBest(0);
            unsigned short hitStartTickBest(0);
            unsigned short midHitTickBest(0);
            
            const recob::Hit* rejectedHit = 0;
            const recob::Hit* bestHit     = 0;
            
            if (hitIter != channelToHitVec.end())
            {

                // Loop through the hits for this channel and look for matches
                // In the event of more than one hit associated to the sim channel range, keep only
                // the best match (assuming the nearby hits are "extra")
                // Note that assumption breaks down for long pulse trains but worry about that later
                for(const auto& hit : hitIter->second)
                {
                    unsigned short hitStartTick = hit->PeakTime() - fSigmaVec.at(plane) * hit->RMS();
                    unsigned short hitStopTick  = hit->PeakTime() + fSigmaVec.at(plane) * hit->RMS();
                    unsigned short midHitTick   = (hitStopTick + hitStartTick) / 2;
                    
                    // If hit is out of range then skip, it is not related to this particle
                    if (hitStartTick > stopTick || hitStopTick < startTick)
                    {
                        if (plane == 1) rejectedHit = hit;
                        continue;
                    }
                    
                    float hitHeight = hit->PeakAmplitude();
                    
                    // Use the hit with the largest pulse height as the "best"
                    if (hitHeight < hitPeakAmpBest) continue;
                    
                    hitPeakAmpBest   = hitHeight;
                    bestHit          = hit;
                    hitStartTickBest = hitStartTick;
                    hitStopTickBest  = hitStopTick;
                    midHitTickBest   = midHitTick;
                }
                
                // Find a match?
                if (bestHit)
                {
                    nElectronsTotalBest = 0.;
                    hitPeakTimeBest     = bestHit->PeakTime();
                    hitIntegralBest     = bestHit->Integral();
                    hitSummedADCBest    = bestHit->SummedADC();
                    hitRMSBest          = bestHit->RMS();
                    hitBaselineBest     = 0.;  // To do...
                    
                    nMatchedHits++;
                    
                    // Get the number of electrons
                    for(unsigned short tick = hitStartTickBest; tick <= hitStopTickBest; tick++)
                    {
                        unsigned short hitTDC = fClockService->TPCTick2TDC(tick - fOffsetVec.at(plane));
                        
                        TDCToIDEMap::iterator ideIterator = tdcToIDEMap.find(hitTDC);
                        
                        if (ideIterator != tdcToIDEMap.end()) nElectronsTotalBest += ideIterator->second;
                    }
                }

                if (nMatchedHits > 0)
                {
                    fHitSumADCVec.at(plane)->Fill(hitSummedADCBest, 1.);
                    fHitVsSimChgVec.at(plane)->Fill(std::min(hitSummedADCBest,float(4999.)), totalElectrons, 1.);
                    fHitPulseHeightVec.at(plane)->Fill(std::min(hitPeakAmpBest,float(149.5)), 1.);
                    fHitPulseWidthVec.at(plane)->Fill(std::min(hitRMSBest,float(19.8)), 1.);
                    fHitElectronsVec.at(plane)->Fill(nElectronsTotalBest, 1.);
                    fHitNumTDCVec.at(plane)->Fill(hitStopTickBest - hitStartTickBest, 1.);
                    fDeltaMidTDCVec.at(plane)->Fill(midHitTickBest - midTick, 1.);
                    
                    nRecobHitVec.at(plane)++;
                }
                else if (rejectedHit)
                {
                    unsigned short hitStartTick = rejectedHit->PeakTime() - fSigmaVec.at(plane) * rejectedHit->RMS();
                    unsigned short hitStopTick  = rejectedHit->PeakTime() + fSigmaVec.at(plane) * rejectedHit->RMS();

                    std::cout << "**> TPC: " << rejectedHit->WireID().TPC << ", Plane " << rejectedHit->WireID().Plane << ", wire: " << rejectedHit->WireID().Wire << ", hit start/stop tick: " << hitStartTick << "/" << hitStopTick << ", start/stop ticks: " << startTick << "/" << stopTick << std::endl;
                    std::cout << "    TPC/Plane/Wire: " << wids[0].TPC << "/" << plane << "/" << wids[0].Wire << ", Track # hits: " << partToChanInfo.second.size() << ", # hits: " << hitIter->second.size() << ", # electrons: " << totalElectrons << ", pulse Height: " << rejectedHit->PeakAmplitude() << ", charge: " << rejectedHit->Integral() << ", " << rejectedHit->SummedADC() << std::endl;
                }
                else
                {
                    std::cout << "==> No match, TPC/Plane/Wire: " << "/" << wids[0].TPC << "/" << wids[0].Plane << "/" << wids[0].Wire << ", # electrons: " << totalElectrons << ", startTick: " << startTick << ", stopTick: " << stopTick << std::endl;
                }
            }
        
            fNMatchedHitVec.at(plane)->Fill(nMatchedHits, 1.);
            fHitEfficVec.at(plane)->Fill(totalElectrons, std::min(nMatchedHits,1), 1.);
            fHitEfficPHVec.at(plane)->Fill(maxElectrons, std::min(nMatchedHits,1), 1.);
            
            // Store tuple variables
            fTPCVec.push_back(wids[0].TPC);
            fCryoVec.push_back(wids[0].Cryostat);
            fPlaneVec.push_back(wids[0].Plane);
            fWireVec.push_back(wids[0].Wire);
            
            fTotalElectronsVec.push_back(totalElectrons);
            fMaxElectronsVec.push_back(maxElectrons);
            fStartTickVec.push_back(startTick);
            fStopTickVec.push_back(stopTick);
            fNMatchedHits = nMatchedHits;

            fHitPeakTimeVec.push_back(hitPeakTimeBest);
            fHitPeakAmpVec.push_back(hitPeakAmpBest);
            fHitPeakRMSVec.push_back(hitRMSBest);
            fHitBaselinevec.push_back(hitBaselineBest);
            fHitSummedADCVec.push_back(hitSummedADCBest);
            fHitIntegralVec.push_back(hitIntegralBest);
            fHitStartTickVec.push_back(hitStartTickBest);
            fHitStopTickVec.push_back(hitStopTickBest);

        }
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
