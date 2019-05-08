////////////////////////////////////////////////////////////////////////
//
// Class:       HitMerger
// Module Type: producer
// File:        HitMerger_module.cc
//
// This module merges hit collections from multiple producers to create a
// single output hit collection (including associations)
//
// Configuration parameters:
//
// HitProducerLabels        - the producers of the recob::Hit objects
//
// Created by Tracy Usher (usher@slac.stanford.edu) on July 17, 2018
//
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <algorithm>
#include <vector>

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Persistency/Common/PtrMaker.h"

#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/ArtDataHelper/HitCreator.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"

class HitMerger : public art::EDProducer
{
public:
    
    // Copnstructors, destructor.
    explicit HitMerger(fhicl::ParameterSet const & pset);
    virtual ~HitMerger();
    
    // Overrides.
    virtual void reconfigure(fhicl::ParameterSet const & pset);
    virtual void produce(art::Event & e);
    virtual void beginJob();
    virtual void endJob();
    
private:
    using RecobHitToPtrMap = std::unordered_map<const recob::Hit*, art::Ptr<recob::Hit>>;
    
    /**
     *  @brief Create recob::Wire to recob::Hit associations
     */
    void makeWireAssns(const art::Event&, art::Assns<recob::Wire, recob::Hit>&, RecobHitToPtrMap&) const;
    
    /**
     *  @brief Create raw::RawDigit to recob::Hit associations
     */
    
    void makeRawDigitAssns(const art::Event&, art::Assns<raw::RawDigit, recob::Hit>&, RecobHitToPtrMap&) const;
    // define vector for hits to make sure of uniform use
    using HitPtrVector = std::vector<art::Ptr<recob::Hit>>;
    
    // Fcl parameters.
    std::vector<art::InputTag>  fHitProducerLabels;         ///< The full collection of hits
};

DEFINE_ART_MODULE(HitMerger)

//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// pset - Fcl parameters.
///
HitMerger::HitMerger(fhicl::ParameterSet const & pset) : EDProducer{pset},
fNumEvent(0)
{
    reconfigure(pset);
    
    produces< std::vector<recob::Hit>>();
    produces< art::Assns<recob::Wire,   recob::Hit>>();
    produces< art::Assns<raw::RawDigit, recob::Hit>>();
    
    // Report.
    mf::LogInfo("HitMerger") << "HitMerger configured\n";
}

//----------------------------------------------------------------------------
/// Destructor.
HitMerger::~HitMerger()
{}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void HitMerger::reconfigure(fhicl::ParameterSet const & pset)
{
    fHitProducerLabels = pset.get<std::vector<art::InputTag>>("HitProducerLabels");
}

//----------------------------------------------------------------------------
/// Begin job method.
void HitMerger::beginJob()
{
//    auto const* detp = lar::providerFrom<detinfo::DetectorPropertiesService>();
//    auto const* geo  = lar::providerFrom<geo::Geometry>();
//    auto const* ts   = lar::providerFrom<detinfo::DetectorClocksService>();
}

//----------------------------------------------------------------------------
/// Produce method.
///
/// Arguments:
///
/// evt - Art event.
///
/// This is the primary method. The goal is to merge input hit collections and
/// output a single hit collection on the backside.
///
void HitMerger::produce(art::Event & evt)
{
    // container for our new hit collection
    std::unique_ptr<std::vector<recob::Hit>> outputHitPtrVec(new std::vector<recob::Hit>);
    
    /// Associations with wires.
    std::unique_ptr<art::Assns<recob::Wire, recob::Hit>> wireAssns(new art::Assns<recob::Wire, recob::Hit>);
    
    /// Associations with raw digits.
    std::unique_ptr<art::Assns<raw::RawDigit, recob::Hit>> rawDigitAssns(new art::Assns<raw::RawDigit, recob::Hit>);

    std::vector<recob::Hit> hitPtrVec;
    RecobHitToPtrMap        recobHitToPtrMap;
    
    // Outside loop over the input hit producers
    for(const auto& inputTag : fHitProducerLabels)
    {
        // Start by looking up the original hits
        art::Handle< std::vector<recob::Hit> > hitHandle;
        evt.getByLabel(inputTag, hitHandle);
        
        for(size_t hitIdx = 0; hitIdx < hitHandle->size(); hitIdx++)
        {
            // Create and save the new recob::Hit with the correct WireID
            hitPtrVec.emplace_back(recob::HitCreator(*hit2D->getHit(), hit3D.getWireIDs()[idx]).copy());
            
            // Recover a pointer to it...
            recob::Hit* newHit = &hitPtrVec.back();
            
            // Create a mapping from this hit to an art Ptr representing it
            recobHitToPtrMap[newHit] = ptrMaker(hitPtrVec.size()-1);
        }
    }
    
    // Set up to make the associations (if desired)
    makeWireAssns(evt, *wireAssns, recobHitToPtrMap);
    
    makeRawDigitAssns(evt, *rawDigitAssns, recobHitToPtrMap);
    
    // Move everything into the event
    evt.put(std::move(outputHitPtrVec));
    evt.put(std::move(wireAssns));
    evt.put(std::move(rawDigitAssns));
        
    // put the hit collection and associations into the event
    hcol.put_into(evt);
    
    return;
}
    
void HitMerger::makeWireAssns(const art::Event& evt, art::Assns<recob::Wire, recob::Hit>& wireAssns, RecobHitToPtrMap& recobHitPtrMap) const
{
    // Let's make sure the input associations container is empty
    wireAssns = art::Assns<recob::Wire, recob::Hit>();
    
    // First task is to recover all of the previous wire <--> hit associations and map them by channel number
    // Create the temporary container
    std::unordered_map<raw::ChannelID_t, art::Ptr<recob::Wire>> channelToWireMap;
    
    // Go through the list of input sources and fill out the map
    for(const auto& inputTag : m_hitFinderTagVec)
    {
        art::ValidHandle<std::vector<recob::Hit>> hitHandle = evt.getValidHandle<std::vector<recob::Hit>>(inputTag);
        
        art::FindOneP<recob::Wire> hitToWireAssns(hitHandle, evt, inputTag);
        
        if (hitToWireAssns.isValid())
        {
            for(size_t wireIdx = 0; wireIdx < hitToWireAssns.size(); wireIdx++)
            {
                art::Ptr<recob::Wire> wire = hitToWireAssns.at(wireIdx);
                
                channelToWireMap[wire->Channel()] = wire;
            }
        }
    }
    
    // Now fill the container
    for(const auto& hitPtrPair : recobHitPtrMap)
    {
        raw::ChannelID_t channel = hitPtrPair.first->Channel();
        
        std::unordered_map<raw::ChannelID_t, art::Ptr<recob::Wire>>::iterator chanWireItr = channelToWireMap.find(channel);
        
        if (!(chanWireItr != channelToWireMap.end()))
        {
            std::cout << "******>> Did not find channel to wire match! Skipping..." << std::endl;
            continue;
        }
        
        wireAssns.addSingle(chanWireItr->second, hitPtrPair.second);
    }
    
    return;
}
    
void HitMerger::makeRawDigitAssns(const art::Event& evt, art::Assns<raw::RawDigit, recob::Hit>& rawDigitAssns, RecobHitToPtrMap& recobHitPtrMap) const
{
    // Let's make sure the input associations container is empty
    rawDigitAssns = art::Assns<raw::RawDigit, recob::Hit>();
    
    // First task is to recover all of the previous wire <--> hit associations and map them by channel number
    // Create the temporary container
    std::unordered_map<raw::ChannelID_t, art::Ptr<raw::RawDigit>> channelToRawDigitMap;
    
    // Go through the list of input sources and fill out the map
    for(const auto& inputTag : m_hitFinderTagVec)
    {
        art::ValidHandle<std::vector<recob::Hit>> hitHandle = evt.getValidHandle<std::vector<recob::Hit>>(inputTag);
        
        art::FindOneP<raw::RawDigit> hitToRawDigitAssns(hitHandle, evt, inputTag);
        
        if (hitToRawDigitAssns.isValid())
        {
            for(size_t rawDigitIdx = 0; rawDigitIdx < hitToRawDigitAssns.size(); rawDigitIdx++)
            {
                art::Ptr<raw::RawDigit> rawDigit = hitToRawDigitAssns.at(rawDigitIdx);
                
                channelToRawDigitMap[rawDigit->Channel()] = rawDigit;
            }
        }
    }
    
    // Now fill the container
    for(const auto& hitPtrPair : recobHitPtrMap)
    {
        raw::ChannelID_t channel = hitPtrPair.first->Channel();
        
        std::unordered_map<raw::ChannelID_t, art::Ptr<raw::RawDigit>>::iterator chanRawDigitItr = channelToRawDigitMap.find(channel);
        
        if (!(chanRawDigitItr != channelToRawDigitMap.end()))
        {
            std::cout << "******>> Did not find channel to RawDigit match! Skipping..." << std::endl;
            continue;
        }
        
        rawDigitAssns.addSingle(chanRawDigitItr->second, hitPtrPair.second);
    }
    
    return;
}
    
//----------------------------------------------------------------------------
/// End job method.
void HitMerger::endJob()
{
    mf::LogInfo("HitMerger")
    << "HitMerger statistics:\n"
    << "  Number of events = " << fNumEvent;
}
