////////////////////////////////////////////////////////////////////////
//
// Class:       HitSelector
// Module Type: producer
// File:        HitSelector_module.cc
//
// This module tries to remove spurious hits which might have been
// created during the deconvolution process
//
// Configuration parameters:
//
// HitProducerLabel        - the producer of the recob::Hit objects
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

#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/ArtDataHelper/HitCreator.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"

class HitSelector : public art::EDProducer
{
public:
    
    // Copnstructors, destructor.
    explicit HitSelector(fhicl::ParameterSet const & pset);
    virtual ~HitSelector();
    
    // Overrides.
    virtual void reconfigure(fhicl::ParameterSet const & pset);
    virtual void produce(art::Event & e);
    virtual void beginJob();
    virtual void endJob();
    
private:
    // define vector for hits to make sure of uniform use
    using HitPtrVector = std::vector<art::Ptr<recob::Hit>>;
    
    // Fcl parameters.
    art::InputTag fHitProducerLabel;         ///< The full collection of hits
    
    // Statistics.
    int           fNumEvent;        ///< Number of events seen.
};

DEFINE_ART_MODULE(HitSelector)

//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// pset - Fcl parameters.
///
HitSelector::HitSelector(fhicl::ParameterSet const & pset) :
fNumEvent(0)
{
    reconfigure(pset);
    
    // let HitCollectionCreator declare that we are going to produce
    // hits and associations with wires and raw digits
    // (with no particular product label)
    recob::HitCollectionCreator::declare_products(*this);
    
    // Report.
    mf::LogInfo("HitSelector") << "HitSelector configured\n";
}

//----------------------------------------------------------------------------
/// Destructor.
HitSelector::~HitSelector()
{}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void HitSelector::reconfigure(fhicl::ParameterSet const & pset)
{
    fHitProducerLabel = pset.get<art::InputTag>("HitProducerLabel");
}

//----------------------------------------------------------------------------
/// Begin job method.
void HitSelector::beginJob()
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
/// This is the primary method. The goal is to produce a list of recob::Hit
/// objects which are a "clean" subset of all hits and which are believed to
/// be due to a neutrino interaction. It does this by considering input CosmicTag
/// objects, relating them to PFParticles/Tracks and removing the hits
/// associated to those objects which are believed to be Cosmic Rays.
///
void HitSelector::produce(art::Event & evt)
{
    ++fNumEvent;
    
    // Start by looking up the original hits
    art::Handle< std::vector<recob::Hit> > hitHandle;
    evt.getByLabel(fHitProducerLabel, hitHandle);
    
    // also get the associated wires and raw digits;
    // we assume they have been created by the same module as the hits
    art::FindOneP<raw::RawDigit> hitToRawDigitAssns(hitHandle, evt, fHitProducerLabel);
    art::FindOneP<recob::Wire>   hitToWireAssns(    hitHandle, evt, fHitProducerLabel);
    
    // this object contains the hit collection
    // and its associations to wires and raw digits:
    recob::HitCollectionCreator hcol(*this,
                                     evt,
                                     hitToWireAssns.isValid(),
                                     hitToRawDigitAssns.isValid()
                                     );

    // If there are no hits then there should be no output
    if (hitHandle.isValid())
    {
        HitPtrVector hitPtrVec;
        HitPtrVector hitSnippetVec;
        
        art::fill_ptr_vector(hitPtrVec, hitHandle);

        int lastHitIndex = 0;
        
        // Loop the hits and make some plots
        for(const auto& curHitPtr : hitPtrVec)
        {
            // Have we collected all of the hits on the same snippet?
            if (!hitSnippetVec.empty() && lastHitIndex >= curHitPtr->LocalIndex())
            {
                // Only worried about multi hit snippets
                if (hitSnippetVec.size() > 1)
                {
                    // Sort in order of largest to smallest pulse height
                    HitPtrVector::iterator maxPulseHeightItr = std::max_element(hitSnippetVec.begin(),hitSnippetVec.end(),[](const auto& left, const auto& right){return left->PeakAmplitude() < right->PeakAmplitude();});
                    
                    float maxPulseHeight = (*maxPulseHeightItr)->PeakAmplitude();
                    
                    if (maxPulseHeight > 6.)
                    {
//                        int maxHitIdx = std::distance(hitSnippetVec.begin(),maxPulseHeightItr);
                        
                        for(size_t idx = 0; idx < hitSnippetVec.size(); idx++)
                        {
                            art::Ptr<recob::Hit> hitPtr = hitSnippetVec.at(idx);
                            
                            float pulseHeight = hitPtr->PeakAmplitude();
                            float pulseWid    = hitPtr->RMS();
                            int   numDOF      = hitPtr->DegreesOfFreedom();
                            
                            if (numDOF == 1 || (pulseHeight > 12. && pulseWid > 3.8))
                            {
                                art::Ptr<recob::Wire>   wire      = hitToWireAssns.at(hitPtr.key());
                                art::Ptr<raw::RawDigit> rawdigits = hitToRawDigitAssns.at(hitPtr.key());
                                
                                hcol.emplace_back(*hitPtr, wire, rawdigits);
                            }
                        }
                    }
                }
                else
                {
                    art::Ptr<recob::Hit> hitPtr = hitSnippetVec.front();
                    
                    float  pulseHeight = hitPtr->PeakAmplitude();
                    float  pulseWid    = hitPtr->RMS();
                    
                    if (pulseHeight > 10. && pulseWid > 3.6)
                    {
                        art::Ptr<recob::Wire>   wire      = hitToWireAssns.at(hitPtr.key());
                        art::Ptr<raw::RawDigit> rawdigits = hitToRawDigitAssns.at(hitPtr.key());
                        
                        hcol.emplace_back(*hitPtr, wire, rawdigits);
                   }
                }
                
                hitSnippetVec.clear();
            }
 
            lastHitIndex = curHitPtr->LocalIndex();
            hitSnippetVec.push_back(curHitPtr);
        }
    }
    
    std::cout << "*** Read " << hitHandle->size() << " hits, writing " << hcol.size() << " selected hits" << std::endl;
    
    // put the hit collection and associations into the event
    hcol.put_into(evt);
    
    return;
}

//----------------------------------------------------------------------------
/// End job method.
void HitSelector::endJob()
{
    mf::LogInfo("HitSelector")
    << "HitSelector statistics:\n"
    << "  Number of events = " << fNumEvent;
}
