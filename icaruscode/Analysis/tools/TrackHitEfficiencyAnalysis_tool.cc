#include "icaruscode/Analysis/tools/IHitEfficiencyHistogramTool.h"

#include "fhiclcpp/ParameterSet.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"

#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include "sbnobj/ICARUS/TPC/ChannelROI.h"

// Eigen
#include <Eigen/Dense>

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TTree.h"

#include <cmath>
#include <algorithm>
#include <tuple>

namespace TrackHitEfficiencyAnalysis
{
    ////////////////////////////////////////////////////////////////////////
    //
    // Class:       TrackHitEfficiencyAnalysis
    // Module Type: producer
    // File:        TrackHitEfficiencyAnalysis.cc
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

class TrackHitEfficiencyAnalysis : virtual public IHitEfficiencyHistogramTool
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    explicit TrackHitEfficiencyAnalysis(fhicl::ParameterSet const & pset);
    
    /**
     *  @brief  Destructor
     */
    ~TrackHitEfficiencyAnalysis();
    
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
    void fillHistograms(const art::Event&)  const override;
    
private:
    
    // Clear mutable variables
    void clear() const;
    
    // Fcl parameters.
    std::vector<art::InputTag>  fRawDigitProducerLabelVec;
    std::vector<art::InputTag>  fWireProducerLabelVec;
    std::vector<art::InputTag>  fHitProducerLabelVec;
    art::InputTag               fMCParticleProducerLabel;
    art::InputTag               fSimChannelProducerLabel;
    art::InputTag               fBadChannelProducerLabel;
    bool                        fUseBadChannelDB;
    std::string                 fLocalDirName;           ///< Fraction for truncated mean
    std::vector<int>            fOffsetVec;              ///< Allow offsets for each plane
    std::vector<float>          fSigmaVec;               ///< Window size for matching to SimChannels
    int                         fMinAllowedChanStatus;   ///< Don't consider channels with lower status
    float                       fSimChannelMinEnergy;

    // TTree variables
    mutable TTree*             fTree;
    
    mutable std::vector<short>   fTPCVec;
    mutable std::vector<short>   fCryoVec;
    mutable std::vector<short>   fPlaneVec;
    mutable std::vector<short>   fWireVec;
    
    mutable std::vector<float>   fTotalElectronsVec;
    mutable std::vector<float>   fMaxElectronsVec;
    mutable std::vector<short>   fmaxElectronsTickVec;
    mutable std::vector<short>   fStartTickVec;
    mutable std::vector<short>   fStopTickVec;
    mutable std::vector<short>   fIDELenVec;
    mutable std::vector<float>   fPartDirX;
    mutable std::vector<float>   fPartDirY;
    mutable std::vector<float>   fPartDirZ;
    mutable std::vector<float>   fCosThetaXZVec;

    mutable std::vector<short>   fRawDigitPulseHeightVec;
    mutable std::vector<short>   fRawDigitMaxTickVec;
    mutable std::vector<short>   fRawDigitMinTickVec;

    mutable std::vector<short>   fNMatchedWires;
    mutable std::vector<short>   fNMatchedHits;

    mutable std::vector<short>   fROIMaxValVec;
    mutable std::vector<short>   fROIMaxTickVec;
    mutable std::vector<short>   fROILenVec;
    mutable std::vector<short>   fROIDeltaTVec;
    
    mutable std::vector<float>   fHitPeakTimeVec;
    mutable std::vector<float>   fHitPeakAmpVec;
    mutable std::vector<float>   fHitPeakRMSVec;
    mutable std::vector<float>   fHitBaselinevec;
    mutable std::vector<float>   fHitSummedADCVec;
    mutable std::vector<float>   fHitIntegralVec;
    mutable std::vector<short>   fHitStartTickVec;
    mutable std::vector<short>   fHitStopTickVec;
    mutable std::vector<short>   fHitDeltaTVec;
    mutable std::vector<short>   fSnippetLengthVec;
    mutable std::vector<short>   fHitMultiplicityVec;
    mutable std::vector<short>   fHitLocalIndexVec;
    mutable std::vector<float>   fHitGoodnessVec;
    mutable std::vector<short>   fNumDegreesVec;

    // Useful services, keep copies for now (we can update during begin run periods)
    const geo::GeometryCore*           fGeometry;             ///< pointer to Geometry service
};
    
//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// pset - Fcl parameters.
///
TrackHitEfficiencyAnalysis::TrackHitEfficiencyAnalysis(fhicl::ParameterSet const & pset) : fTree(nullptr)
{
    fGeometry           = lar::providerFrom<geo::Geometry>();
    
    configure(pset);
    
    // Report.
    mf::LogInfo("TrackHitEfficiencyAnalysis") << "TrackHitEfficiencyAnalysis configured\n";
}

//----------------------------------------------------------------------------
/// Destructor.
TrackHitEfficiencyAnalysis::~TrackHitEfficiencyAnalysis()
{}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void TrackHitEfficiencyAnalysis::configure(fhicl::ParameterSet const & pset)
{
    fRawDigitProducerLabelVec = pset.get< std::vector<art::InputTag>>("RawDigitLabelVec",   std::vector<art::InputTag>() = {"rawdigitfilter"});
    fWireProducerLabelVec     = pset.get< std::vector<art::InputTag>>("WireModuleLabelVec", std::vector<art::InputTag>() = {"decon1droi"});
    fHitProducerLabelVec      = pset.get< std::vector<art::InputTag>>("HitModuleLabelVec",  std::vector<art::InputTag>() = {"gauss"});
    fMCParticleProducerLabel  = pset.get< art::InputTag             >("MCParticleLabel",    "largeant");
    fSimChannelProducerLabel  = pset.get< art::InputTag             >("SimChannelLabel",    "largeant");
    fBadChannelProducerLabel  = pset.get< art::InputTag             >("BadChannelLabel",    "simnfspl1:badchannels");
    fUseBadChannelDB          = pset.get< bool                      >("UseBadChannelDB",    true);
    fLocalDirName             = pset.get<std::string                >("LocalDirName",       std::string("wow"));
    fOffsetVec                = pset.get<std::vector<int>           >("OffsetVec",          std::vector<int>()={0,0,0});
    fSigmaVec                 = pset.get<std::vector<float>         >("SigmaVec",           std::vector<float>()={1.,1.,1.});
    fMinAllowedChanStatus     = pset.get< int                       >("MinAllowedChannelStatus");
    fSimChannelMinEnergy      = pset.get<float                      >("SimChannelMinEnergy", std::numeric_limits<float>::epsilon());
}

//----------------------------------------------------------------------------
/// Begin job method.
void TrackHitEfficiencyAnalysis::initializeHists(art::ServiceHandle<art::TFileService>& tfs, const std::string& dirName)
{
    // Make a directory for these histograms
//    art::TFileDirectory dir = tfs->mkdir(dirName.c_str());
    
    return;
}
    
void TrackHitEfficiencyAnalysis::initializeTuple(TTree* tree)
{
    fTree = tree;
    
    fTree->Branch("CryostataVec",        "std::vector<short>",   &fCryoVec);
    fTree->Branch("TPCVec",              "std::vector<short>",   &fTPCVec);
    fTree->Branch("PlaneVec",            "std::vector<short>",   &fPlaneVec);
    fTree->Branch("WireVec",             "std::vector<short>",   &fWireVec);
    
    fTree->Branch("TotalElectronsVec",   "std::vector<float>",   &fTotalElectronsVec);
    fTree->Branch("MaxElectronsVec",     "std::vector<float>",   &fMaxElectronsVec);
    fTree->Branch("StartTick",           "std::vector<short>",   &fStartTickVec);
    fTree->Branch("StopTick",            "std::vector<short>",   &fStopTickVec);
    fTree->Branch("maxElectronsTick",    "std::vector<short>",   &fmaxElectronsTickVec);
    fTree->Branch("IDELength",           "std::vector<short",    &fIDELenVec);
    fTree->Branch("PartDirX",            "std::vector<float>",   &fPartDirX);
    fTree->Branch("PartDirY",            "std::vector<float>",   &fPartDirY);
    fTree->Branch("PartDirZ",            "std::vector<float>",   &fPartDirZ);
    fTree->Branch("CosThetaXZ",          "std::vector<float>",   &fCosThetaXZVec);
    
    fTree->Branch("NMatchedWires",       "std::vector<short>",   &fNMatchedWires);
    fTree->Branch("NMatchedHits",        "std::vector<short>",   &fNMatchedHits);

    fTree->Branch("RawDigitPulseHeight", "std::vector<short>",   &fRawDigitPulseHeightVec);
    fTree->Branch("RawDigitMaxTick",     "std::vector<short>",   &fRawDigitMaxTickVec);
    fTree->Branch("RawDigitMinTick",     "std::vector<short>",   &fRawDigitMinTickVec);

    fTree->Branch("ROIMaxValue",         "std::vector<short>",   &fROIMaxValVec);
    fTree->Branch("ROITickValue",        "std::vector<short>",   &fROIMaxTickVec);
    fTree->Branch("ROILength",           "std::vector<short>",   &fROILenVec);
    fTree->Branch("ROIDeltaT",           "std::vector<short>",   &fROIDeltaTVec);

    fTree->Branch("HitPeakTimeVec",      "std::vector<float>",   &fHitPeakTimeVec);
    fTree->Branch("HitPeakAmpVec",       "std::vector<float>",   &fHitPeakAmpVec);
    fTree->Branch("HitPeakRMSVec",       "std::vector<float>",   &fHitPeakRMSVec);
    fTree->Branch("HitBaselineVec",      "std::vector<float>",   &fHitBaselinevec);
    fTree->Branch("HitSummedADCVec",     "std::vector<float>",   &fHitSummedADCVec);
    fTree->Branch("HitIntegralVec",      "std::vector<float>",   &fHitIntegralVec);
    fTree->Branch("HitStartTickVec",     "std::vector<short>",   &fHitStartTickVec);
    fTree->Branch("HitStopTickVec",      "std::vector<short>",   &fHitStopTickVec);
    fTree->Branch("HitDeltaTVec",        "std::vector<short>",   &fHitDeltaTVec);
    fTree->Branch("SnippetLengthkVec",   "std::vector<short>",   &fSnippetLengthVec);
    fTree->Branch("HitMultiplicity",     "std::vector<short>",   &fHitMultiplicityVec);
    fTree->Branch("HitLocalIndex",       "std::vector<short>",   &fHitLocalIndexVec);
    fTree->Branch("HitGoodness",         "std::vector<float>",   &fHitGoodnessVec);
    fTree->Branch("HitNumDegrees",       "std::vector<short>",   &fNumDegreesVec);
  
    clear();

    return;
}
    
void TrackHitEfficiencyAnalysis::clear() const
{
    fTPCVec.clear();
    fCryoVec.clear();
    fPlaneVec.clear();
    fWireVec.clear();

    fTotalElectronsVec.clear();
    fMaxElectronsVec.clear();
    fStartTickVec.clear();
    fStopTickVec.clear();
    fmaxElectronsTickVec.clear();
    fIDELenVec.clear();
    fPartDirX.clear();
    fPartDirY.clear();
    fPartDirZ.clear();
    fCosThetaXZVec.clear();

    fRawDigitPulseHeightVec.clear();
    fRawDigitMaxTickVec.clear();
    fRawDigitMinTickVec.clear();
    
    fNMatchedWires.clear();
    fNMatchedHits.clear();

    fROIMaxValVec.clear();
    fROIMaxTickVec.clear();
    fROILenVec.clear();
    fROIDeltaTVec.clear();

    fHitPeakTimeVec.clear();
    fHitPeakAmpVec.clear();
    fHitPeakRMSVec.clear();
    fHitBaselinevec.clear();
    fHitSummedADCVec.clear();
    fHitIntegralVec.clear();
    fHitStartTickVec.clear();
    fHitStopTickVec.clear();
    fHitDeltaTVec.clear();
    fSnippetLengthVec.clear();
    fHitMultiplicityVec.clear();
    fHitLocalIndexVec.clear();
    fHitGoodnessVec.clear();
    fNumDegreesVec.clear();
    
    return;
}

void TrackHitEfficiencyAnalysis::fillHistograms(const art::Event& event) const
{
   // std::cout << " filling histos " << std::endl;
    // Basic assumption is that the producer label vecs for RawDigits and Wire data are
    // all the same length and in the same order. Here we just check for length
    if (fRawDigitProducerLabelVec.size() != fWireProducerLabelVec.size()) return;
    
    // Always clear the tuple
    clear();
    
    art::Handle< std::vector<sim::SimChannel>> simChannelHandle;
    event.getByLabel(fSimChannelProducerLabel, simChannelHandle);
    
    art::Handle< std::vector<simb::MCParticle>> mcParticleHandle;
    event.getByLabel(fMCParticleProducerLabel, mcParticleHandle);

    // If there is no sim channel informaton then exit
    if (!simChannelHandle.isValid() || simChannelHandle->empty() || !mcParticleHandle.isValid()) return;
    
    // what needs to be done?
    // First we define a straightforward channel to Wire map so we can look up a given
    // channel's Wire data as we loop over SimChannels.
    using ChanToWireMap = std::unordered_map<raw::ChannelID_t,const recob::ChannelROI*>;
    
    ChanToWireMap channelToWireMap;
    
    // We will use the presence of a RawDigit as an indicator of a good channel... So
    // we want a mapping between channel and RawDigit
    using ChanToRawDigitMap = std::unordered_map<raw::ChannelID_t,const raw::RawDigit*>;
    
    ChanToRawDigitMap chanToRawDigitMap;

    // Now start a loop over the individual TPCs to build out the structures for RawDigits and Wires
    for(size_t tpcID = 0; tpcID < fRawDigitProducerLabelVec.size(); tpcID++)
    {
        art::Handle< std::vector<raw::RawDigit> > rawDigitHandle;
        event.getByLabel(fRawDigitProducerLabelVec[tpcID], rawDigitHandle);

        art::Handle< std::vector<recob::ChannelROI> > wireHandle;
        event.getByLabel(fWireProducerLabelVec[tpcID], wireHandle);

        if (!rawDigitHandle.isValid() || !wireHandle.isValid()) return;
        
        for(const auto& wire : *wireHandle) channelToWireMap[wire.Channel()] = &wire;
        
        for(const auto& rawDigit : *rawDigitHandle) chanToRawDigitMap[rawDigit.Channel()] = &rawDigit;
    }
    
    // Now we create a data structure to relate hits to their channel ID
    using ChanToHitVecMap = std::unordered_map<raw::ChannelID_t,std::vector<const recob::Hit*>>;
    
    ChanToHitVecMap channelToHitVec;
    
    // And now fill it
    for(const auto& hitLabel : fHitProducerLabelVec)
    {
        art::Handle< std::vector<recob::Hit> > hitHandle;
        event.getByLabel(hitLabel, hitHandle);

        for(const auto& hit : *hitHandle) channelToHitVec[hit.Channel()].emplace_back(&hit);
    }
    
    // It is useful to create a mapping between trackID and MCParticle
    using TrackIDToMCParticleMap = std::unordered_map<int, const simb::MCParticle*>;
    
    TrackIDToMCParticleMap trackIDToMCParticleMap;
    
    for(const auto& mcParticle : *mcParticleHandle) trackIDToMCParticleMap[mcParticle.TrackId()] = &mcParticle;
    
    std::vector<int> nSimChannelHitVec  = {0,0,0};
    std::vector<int> nRecobHitVec       = {0,0,0};
    std::vector<int> nFakeHitVec        = {0,0,0};
    std::vector<int> nSimulatedWiresVec = {0,0,0};
    
    unsigned int lastwire=-1;
    
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event);
    
//    const lariov::ChannelStatusProvider& chanFilt = art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider();
    
    // There are several things going on here... for each channel we have particles (track id's) depositing energy in a range to ticks
    // So... for each channel we want to build a structure that relates particles to tdc ranges and deposited energy (or electrons)
    // Here is a complicated structure:
    using TickTDCIDEVec           = std::vector<const sim::TDCIDE*>;
    using ChanToTDCIDEMap         = std::unordered_map<raw::ChannelID_t,TickTDCIDEVec>;
    using PartToChanToTDCToIDEMap = std::unordered_map<int, ChanToTDCIDEMap>;
    using ChannelToSimChannelMap  = std::unordered_map<raw::ChannelID_t,const sim::SimChannel*>;
    
    PartToChanToTDCToIDEMap partToChanToTDCToIDEMap;
    ChannelToSimChannelMap  channelToSimChannelMap;

    using StartStopChargeTuple    = std::tuple<short,short,float,int>;  // start, stop, total charge
    using StartStopChargeTupleMap = std::map<int,StartStopChargeTuple>; // reference by track id
    
    // Build out the above data structure
    for(const auto& simChannel : *simChannelHandle)
    {
        raw::ChannelID_t channel = simChannel.Channel();

        // If no deposited charge then skip the channel
        if (simChannel.TDCIDEMap().size() == 0) continue;

        sim::SimChannel::TDCIDEs_t::const_iterator tdcIdeItr = simChannel.TDCIDEMap().begin();

//        std::cout << "****************** Channel " << channel << " ***********************" << std::endl;

        // Loop through all of the TDC values in the "map"
        while(true)
        {
            // If we hit the end and haven't processed anything then exit
            if (tdcIdeItr == simChannel.TDCIDEMap().end()) break;

            sim::SimChannel::TDCIDEs_t::const_iterator startItr(tdcIdeItr);
            sim::SimChannel::TDCIDEs_t::const_iterator stopItr(tdcIdeItr);
            sim::SimChannel::TDCIDEs_t::const_iterator maxElectronsItr(tdcIdeItr);

            float totalElectrons(0.);
            float maxElectrons(0.);

            StartStopChargeTupleMap startStopChargeTupleMap;

//            std::cout << "  --> First TDC value for map: " << simChannel.TDCIDEMap().begin()->first << std::endl;

            // Loop through continuous groups of TDC values
            while(tdcIdeItr != simChannel.TDCIDEMap().end())
            {
                int   bestTrackID(0);

                short curTDC = tdcIdeItr->first;

                // If there is a gap to the last TDC then we process those values 
                // (we'll then re-enter this loop looking for the next set)
                if (curTDC - stopItr->first > 7) break;

                float electronsThisTDC(0.);
                float maxElectronsByTrackID(0.);

                for(const auto& ide : tdcIdeItr->second) 
                {
                    electronsThisTDC += ide.numElectrons;

                    if (maxElectronsByTrackID < ide.numElectrons)
                    {
                        bestTrackID           = ide.trackID;
                        maxElectronsByTrackID = ide.numElectrons;
                    }

                    StartStopChargeTupleMap::iterator mapItr = startStopChargeTupleMap.find(ide.trackID);

                    if (mapItr == startStopChargeTupleMap.end()) startStopChargeTupleMap[ide.trackID] = StartStopChargeTuple(tdcIdeItr->first,tdcIdeItr->first,ide.numElectrons,ide.origTrackID);
                    else
                    {
                        StartStopChargeTuple& startStopChargeTuple = mapItr->second;

                        if (tdcIdeItr->first < std::get<0>(startStopChargeTuple)) std::get<0>(startStopChargeTuple) = tdcIdeItr->first;
                        if (tdcIdeItr->first > std::get<1>(startStopChargeTuple)) std::get<1>(startStopChargeTuple) = tdcIdeItr->first;

                        std::get<2>(startStopChargeTuple) += ide.numElectrons;
                    }
                }

                totalElectrons += electronsThisTDC;

                if (electronsThisTDC > maxElectrons)
                {
                    maxElectrons    = electronsThisTDC;
                    maxElectronsItr = tdcIdeItr;
                }

                partToChanToTDCToIDEMap[bestTrackID][channel].emplace_back(&(*tdcIdeItr));

                stopItr = tdcIdeItr++;
            }

            int   nMatchedWires(0);
            int   nMatchedHits(0);
            int   bestTrackID(0);
//            int   bestOrigTrackID(0);
//            short bestStartTDC(0);
//            short bestStopTDC(0);

            float trackElectrons(0.);

            // Look for the MC Track making the largest overall contribution to deposited charge
            for(const auto& mapItr : startStopChargeTupleMap)
            {
                if (std::get<2>(mapItr.second) > trackElectrons)
                {
                    bestTrackID     = mapItr.first;
//                    bestOrigTrackID = std::get<3>(mapItr.second);
//                    bestStartTDC    = std::get<0>(mapItr.second);
//                    bestStopTDC     = std::get<1>(mapItr.second);
                    trackElectrons  = std::get<2>(mapItr.second);
                }
            }

//            std::cout << "  --> start/stop TDC: " << startItr->first << "/" << stopItr->first << ", bestTrackID: " << bestTrackID << ", start/stop: " << bestStartTDC << "/" << bestStopTDC << ", largestCharge: " << maxElectrons << ", total e: " << totalElectrons << std::endl;

            // Recover the MCParticle associated to the best track ID
            TrackIDToMCParticleMap::const_iterator trackIDToMCPartItr = trackIDToMCParticleMap.find(std::abs(bestTrackID));
            //TrackIDToMCParticleMap::const_iterator trackIDToMCPartItr = trackIDToMCParticleMap.find(bestOrigTrackID);

            // Require we have an MCParticle
            if (trackIDToMCPartItr != trackIDToMCParticleMap.end())
            {
                int         trackPDGCode = trackIDToMCPartItr->second->PdgCode();
//                int         trackID      = trackIDToMCPartItr->first;
//                std::string processName  = trackIDToMCPartItr->second->Process();

//                std::cout << "#### trackID: " << bestTrackID << ", orig: " << bestOrigTrackID << ", trackPDGCode: " << trackPDGCode << ", process: " << processName << std::endl;

                // Looking for primary muons (e.g. CR Tracks)
//                if (fabs(trackPDGCode) != 13 || processName != "primary") continue;
                if (std::abs(trackPDGCode) != 13) continue;

                // Recover initial/final particle direction information
                Eigen::Vector3f partStartDir(trackIDToMCPartItr->second->Px(),trackIDToMCPartItr->second->Py(),trackIDToMCPartItr->second->Pz());
                Eigen::Vector3f partEndDir(trackIDToMCPartItr->second->EndPx(),trackIDToMCPartItr->second->EndPy(),trackIDToMCPartItr->second->EndPz());

                partStartDir /= trackIDToMCPartItr->second->P();
                partEndDir   /= trackIDToMCPartItr->second->EndMomentum()[3];

                Eigen::Vector3f partAveDir = 0.5 * (partStartDir + partEndDir);

                // Make sure normalized
                partAveDir.normalize();
                
                float cosThetaXZ = partAveDir[0];

                // The below try-catch block may no longer be necessary
                // Decode the channel and make sure we have a valid one
                std::vector<geo::WireID> wids = fGeometry->ChannelToWire(channel);

                // Recover plane and wire in the plane
                unsigned int plane = wids[0].Plane;
                unsigned int wire  = wids[0].Wire;

                nSimChannelHitVec[plane]++;

                if(wire!=lastwire) nSimulatedWiresVec[plane]++;

                lastwire=wire;

                // Convert to ticks to get in same units as hits
//                unsigned short startTick        = clockData.TPCTDC2Tick(startItr->first)        + fOffsetVec[plane];
//                unsigned short stopTick         = clockData.TPCTDC2Tick(stopItr->first)         + fOffsetVec[plane];
//                unsigned short maxElectronsTick = clockData.TPCTDC2Tick(maxElectronsItr->first) + fOffsetVec[plane];
                unsigned short startTick        = clockData.TPCTDC2Tick(startItr->first);
                unsigned short stopTick         = clockData.TPCTDC2Tick(stopItr->first);
                unsigned short maxElectronsTick = clockData.TPCTDC2Tick(maxElectronsItr->first);

                if (startTick < 6)        startTick    = 6;
                if (stopTick < 6)         stopTick     = 6;
                if (maxElectronsTick < 6) maxElectrons = 6;

                startTick        += fOffsetVec[plane];
                stopTick         += fOffsetVec[plane];
                maxElectronsTick += fOffsetVec[plane];

//                std::cout << "  --> startTick/stopTick/maxElectronsTick: " << startTick << "/" << stopTick << "/" << maxElectronsTick 
//                          << " which is from TDC values: " << startItr->first << "/" << stopItr->first << "/" << maxElectronsItr->first << std::endl;

                // Apparently it can happen that we have a start tick that exceeds the length of the input waveform. 
                // We should also make cuts at edges of readout so we don't have distorted waveforms
                // When this happens skip
                if (startTick > 10 && stopTick < 4085) 
                {
                    // Set up to extract the "best" parameters in the event of more than one hit for this pulse train
                    short          roiMaxValue(0);
                    short          roiMaxValueTick(0);
                    short          roiLen(0);
                    short          roiDeltaT(-2048);
                    float          hitSummedADCBest(0.);
                    float          hitIntegralBest(0.);
                    float          hitPeakTimeBest(0.);
                    float          hitPeakAmpBest(-100.);
                    float          hitRMSBest(0.);
                    int            hitMultiplicityBest(0);
                    int            hitLocalIndexBest(0);
                    float          hitGoodnessBest(0.);
                    int            hitNumDegreesBest(0);
                    float          hitBaselineBest(0.);
                    float          hitSnippetLenBest(0.);
                    unsigned short hitStopTickBest(0);
                    unsigned short hitStartTickBest(0);
                    int            hitDeltaTBest(4096);
                    unsigned short rawDigitPulseHeight(0);
                    int            rawDigitMaxTick(0);
                    int            rawDigitMinTick(0);

                    // Start by getting an estimate of the pulse height from the RawDigits
                    ChanToRawDigitMap::const_iterator rawDigitItr = chanToRawDigitMap.find(channel); 

                    if (rawDigitItr != chanToRawDigitMap.end())
                    {
                        const raw::RawDigit*      rawDigit = rawDigitItr->second;
                        const std::vector<short>& waveform = rawDigit->ADCs();

                        std::vector<short>::const_iterator firstItr  = waveform.begin() + std::max(0,int(startTick)-15);
                        std::vector<short>::const_iterator lastItr   = waveform.begin() + std::min(4095,int(stopTick)+10);
                        std::vector<short>::const_iterator maxADCItr = std::max_element(firstItr,lastItr);
                        std::vector<short>::const_iterator minADCItr = std::min_element(maxADCItr,lastItr);

                        if (plane == 2)
                            rawDigitPulseHeight = *maxADCItr;
                        else
                        {
                            maxADCItr = std::max_element(firstItr,lastItr-10);
                            minADCItr = std::min_element(maxADCItr,lastItr);
                            rawDigitPulseHeight = *maxADCItr - *minADCItr;
                        }
                        rawDigitMaxTick = std::distance(waveform.begin(),maxADCItr);
                        rawDigitMinTick = std::distance(waveform.begin(),minADCItr);
                    }

                    // Now recover the Wire associated to this channel
                    ChanToWireMap::const_iterator wireItr = channelToWireMap.find(channel);

                    if (wireItr != channelToWireMap.end())
                    {
                        const recob::ChannelROI::RegionsOfInterest_t& signalROI = wireItr->second->SignalROI();
                        const lar::sparse_vector<short>::datarange_t* wireRangePtr(NULL);

                        // Here we need to match the range of the ROI's on the given Wire with the tick range from the SimChannel
                        for(const auto& range : signalROI.get_ranges())
                        {
                            // #################################################
                            // ### Getting a vector of signals for this wire ###
                            // #################################################
                            //std::vector<float> signal(wire->Signal());
                            raw::TDCtick_t roiFirstBinTick = range.begin_index();
                            raw::TDCtick_t roiLastBinTick  = range.end_index();

                            if (roiFirstBinTick > stopTick)  break;
                            if (roiLastBinTick  < startTick) continue;

                            // Require the simulated charge deposit is full contained in the ROI
                            if (startTick > roiFirstBinTick && stopTick < roiLastBinTick) wireRangePtr = &range;
                            break;
                        }

                        // Check that we have found the wire range
                        // Note that if we have not matched an ROI then we can't have a hit either so skip search for that...
                        if (wireRangePtr)
                        {
                            // Get the ROI parameters we want to save
                            std::vector<short>::const_iterator maxItr = std::max_element(wireRangePtr->data().begin(),wireRangePtr->data().end());

                            roiMaxValue     = *maxItr;
                            roiMaxValueTick = std::distance(wireRangePtr->data().begin(),maxItr) + wireRangePtr->begin_index();
                            roiLen          = wireRangePtr->data().size();
                            roiDeltaT       = roiMaxValueTick - maxElectronsTick;

                            const recob::Hit* rejectedHit = 0;
                            const recob::Hit* bestHit     = 0;

                            nMatchedWires++;
                            // The next mission is to recover the hits associated to this Wire
                            // The easiest way to do this is to simply look up all the hits on this channel and then match
                            ChanToHitVecMap::iterator hitIter = channelToHitVec.find(channel);

//                            std::cout << "  --> roiLen, roiDeltaT: " << roiLen << ", " << roiDeltaT << std::endl;

                            if (hitIter != channelToHitVec.end())
                            {
                                // Loop through the hits for this channel and look for matches
                                // In the event of more than one hit associated to the sim channel range, keep only
                                // the best match (assuming the nearby hits are "extra")
                                // Note that assumption breaks down for long pulse trains but worry about that later
                                for(const auto& hit : hitIter->second)
                                {
                                    int            hitPeakTime  = hit->PeakTime();
                                    unsigned short hitStartTick = hitPeakTime - fSigmaVec[plane] * hit->RMS();
                                    unsigned short hitStopTick  = hitPeakTime + fSigmaVec[plane] * hit->RMS();
                                    // If hit is out of range then skip, it is not related to this particle
                                    if (hitStartTick > stopTick || hitStopTick < startTick || std::abs(hitPeakTime - maxElectronsTick) > std::abs(hitDeltaTBest))
                                    {
                                        nFakeHitVec[plane]++;
                                        rejectedHit = hit;
                                        continue;
                                    }
                                    float hitHeight  = hit->PeakAmplitude();
                                    hitPeakAmpBest   = hitHeight;
                                    bestHit          = hit;
                                    hitStartTickBest = hitStartTick;
                                    hitStopTickBest  = hitStopTick;
                                    hitDeltaTBest    = hitPeakTime - maxElectronsTick;
                                }

                                // Find a match?
                                if (bestHit)
                                {
                                    hitPeakTimeBest     = bestHit->PeakTime();
                                    hitIntegralBest     = bestHit->Integral();
                                    hitSummedADCBest    = bestHit->SummedADC();
                                    hitRMSBest          = bestHit->RMS();
                                    hitMultiplicityBest = bestHit->Multiplicity();
                                    hitLocalIndexBest   = bestHit->LocalIndex();
                                    hitGoodnessBest     = bestHit->GoodnessOfFit();
                                    hitNumDegreesBest   = bestHit->DegreesOfFreedom();
                                    hitSnippetLenBest   = bestHit->EndTick() - bestHit->StartTick();
                                    hitBaselineBest     = 0.;  // To do...
                                    nMatchedHits++;
                                }

                                if (nMatchedHits > 0)
                                    nRecobHitVec[plane]++;
                                else if (rejectedHit)
                                { 
                                    unsigned short hitStartTick = rejectedHit->PeakTime() - fSigmaVec[plane] * rejectedHit->RMS();
                                    unsigned short hitStopTick  = rejectedHit->PeakTime() + fSigmaVec[plane] * rejectedHit->RMS();

                                    mf::LogDebug("TrackHitEfficiencyAnalysis") << "**> TPC: " << rejectedHit->WireID().TPC << ", Plane " << rejectedHit->WireID().Plane << ", wire: " << rejectedHit->WireID().Wire << ", hit startstop            tick: " << hitStartTick << "/" << hitStopTick << ", start/stop ticks: " << startTick << "/" << stopTick << std::endl;
                                    mf::LogDebug("TrackHitEfficiencyAnalysis") << "    TPC/Plane/Wire: " << wids[0].TPC << "/" << plane << "/" << wids[0].Wire  << ", # hits: "<< hitIter->second.size() << ", # electrons: " << totalElectrons << ", pulse Height: " << rejectedHit->PeakAmplitude() << ", charge: " << rejectedHit->Integral()      << ", " <<rejectedHit->SummedADC() << std::endl;
                                }
                                else
                                {
                                    mf::LogDebug("TrackHitEfficiencyAnalysis") << "==> No match, TPC/Plane/Wire: " << "/" << wids[0].TPC << "/" << wids[0].Plane << "/" << wids[0].Wire << ", # electrons: " << totalElectrons << ",           startTick: " << startTick << ", stopTick: " << stopTick << std::endl;
                                }
                            }
                        }
                    }

                    // Store tuple variables
                    fTPCVec.emplace_back(wids[0].TPC);
                    fCryoVec.emplace_back(wids[0].Cryostat);
                    fPlaneVec.emplace_back(wids[0].Plane);
                    fWireVec.emplace_back(wids[0].Wire);
                    fTotalElectronsVec.emplace_back(totalElectrons);
                    fMaxElectronsVec.emplace_back(maxElectrons);
                    fStartTickVec.emplace_back(startTick);
                    fStopTickVec.emplace_back(stopTick);
                    fmaxElectronsTickVec.emplace_back(maxElectronsTick);
                    fIDELenVec.emplace_back(stopTick-startTick);
                    fPartDirX.emplace_back(partStartDir[0]);
                    fPartDirY.emplace_back(partStartDir[1]);
                    fPartDirZ.emplace_back(partStartDir[2]);
                    fCosThetaXZVec.emplace_back(cosThetaXZ);
                    fRawDigitPulseHeightVec.emplace_back(rawDigitPulseHeight);
                    fRawDigitMaxTickVec.emplace_back(rawDigitMaxTick);
                    fRawDigitMinTickVec.emplace_back(rawDigitMinTick);
                    fNMatchedWires.emplace_back(nMatchedWires);
                    fNMatchedHits.emplace_back(nMatchedHits);
                    fROIMaxValVec.emplace_back(roiMaxValue);
                    fROIMaxTickVec.emplace_back(roiMaxValueTick);
                    fROILenVec.emplace_back(roiLen);
                    fROIDeltaTVec.emplace_back(roiDeltaT);
                    fHitPeakTimeVec.emplace_back(hitPeakTimeBest);
                    fHitPeakAmpVec.emplace_back(hitPeakAmpBest);
                    fHitPeakRMSVec.emplace_back(hitRMSBest);
                    fHitBaselinevec.emplace_back(hitBaselineBest);
                    fHitSummedADCVec.emplace_back(hitSummedADCBest);
                    fHitIntegralVec.emplace_back(hitIntegralBest);
                    fHitStartTickVec.emplace_back(hitStartTickBest);
                    fHitStopTickVec.emplace_back(hitStopTickBest);
                    fHitDeltaTVec.emplace_back(hitDeltaTBest);
                    fSnippetLengthVec.emplace_back(hitSnippetLenBest);
                    fHitMultiplicityVec.emplace_back(hitMultiplicityBest);
                    fHitLocalIndexVec.emplace_back(hitLocalIndexBest);
                    fHitGoodnessVec.emplace_back(hitGoodnessBest);
                    fNumDegreesVec.emplace_back(hitNumDegreesBest);
                }
            }
        } // Looping over TDC bins

        channelToSimChannelMap[channel] = &simChannel;
    } // end of loop over SimChannels

    return;
}
    
// Useful for normalizing histograms
void TrackHitEfficiencyAnalysis::endJob(int numEvents)
{
    return;
}
    
DEFINE_ART_CLASS_TOOL(TrackHitEfficiencyAnalysis)
}
