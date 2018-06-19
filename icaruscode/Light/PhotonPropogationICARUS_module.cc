////////////////////////////////////////////////////////////////////////
//
// Class:       PhotonPropogationICARUS
// Module Type: producer
// File:        PhotonPropogationICARUS_module.cc
//
//              The intent of this module is to modify the SimPhotons
//              (e.g. to account for timing offsets)
//
// Configuration parameters:
//
// SimPhotonModuleLabel      - the input source of the SimPhoton collection
//
//
// Created by Andrea Falcone (andrea.falcone@uta.edu) on June 19, 2018
//
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <algorithm>
#include <vector>

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "nutools/RandomUtils/NuRandomService.h"

#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larsim/PhotonPropagation/PhotonVisibilityService.h"

#include "lardataobj/Simulation/SimPhotons.h"

class PhotonPropogationICARUS : public art::EDProducer
{
public:

    // Copnstructors, destructor.
    explicit PhotonPropogationICARUS(fhicl::ParameterSet const & pset);
    virtual ~PhotonPropogationICARUS();

    // Overrides.
    virtual void configure(fhicl::ParameterSet const & pset);
    virtual void produce(art::Event & e);
    virtual void beginJob();
    virtual void endJob();

private:

    // Fcl parameters.
    art::InputTag fSimPhotonModuleLabel;      ///< The full collection of SimPhotons

    // Statistics.
    int fNumEvent;        ///< Number of events seen.
    
    // Useful services, keep copies for now (we can update during begin run periods)
    geo::GeometryCore const*           fGeometry;             ///< pointer to Geometry service
    detinfo::DetectorProperties const* fDetectorProperties;   ///< Detector properties service
};

DEFINE_ART_MODULE(PhotonPropogationICARUS)

//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// pset - Fcl parameters.
///
PhotonPropogationICARUS::PhotonPropogationICARUS(fhicl::ParameterSet const & pset) :
                      fNumEvent(0)
{
    
    fGeometry = lar::providerFrom<geo::Geometry>();
    fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
    
    configure(pset);
    
    produces<std::vector<sim::SimPhotons>>();
    
    // Not sure if you need this?
//    art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this, "HepJamesRandom", "icarusphoton", pset, "SeedPhoton");

    // Report.
    mf::LogInfo("PhotonPropogationICARUS") << "PhotonPropogationICARUS configured\n";
}

//----------------------------------------------------------------------------
/// Destructor.
PhotonPropogationICARUS::~PhotonPropogationICARUS()
{}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void PhotonPropogationICARUS::configure(fhicl::ParameterSet const & pset)
{
    fSimPhotonModuleLabel = pset.get<art::InputTag>("SimPhotonsModuleLabel");
}

//----------------------------------------------------------------------------
/// Begin job method.
void PhotonPropogationICARUS::beginJob()
{
    // Access ART's TFileService, which will handle creating and writing
    // histograms and n-tuples for us.
    art::ServiceHandle<art::TFileService> tfs;
    
//    art::TFileDirectory dir = tfs->mkdir(Form("PhotonPropogation"));

    return;
}

//----------------------------------------------------------------------------
/// Produce method.
///
/// Arguments:
///
/// evt - Art event.
///
/// This is the primary method.
///
void PhotonPropogationICARUS::produce(art::Event & event)
{
    ++fNumEvent;
    
    // Agreed convention is to ALWAYS output to the event store so get a pointer to our collection
    std::unique_ptr<std::vector<sim::SimPhotons>> simPhotonsVec(new std::vector<sim::SimPhotons>);
    
    // Read in the digit List object(s).
    art::Handle< std::vector<sim::SimPhotons> > simPhotonsVecHandle;
    event.getByLabel(fSimPhotonModuleLabel, simPhotonsVecHandle);
    
    // Require a valid handle
    if (simPhotonsVecHandle.isValid())
    {
        // Recover useful service...
        art::ServiceHandle<phot::PhotonVisibilityService> pvs;
        
        // Loop through the input photons (this might need to be more complicated...)
        for(const auto& simPhoton : *simPhotonsVecHandle)
        {
            // Make a copy?
            sim::SimPhotons localPhoton = simPhoton;
            
            // Do something to it?
            
            // Add to new output collection?
            simPhotonsVec->emplace_back(localPhoton);
        }
    }
    
    // Add tracks and associations to event.
    event.put(std::move(simPhotonsVec));
}

//----------------------------------------------------------------------------
/// End job method.
void PhotonPropogationICARUS::endJob()
{
    mf::LogInfo("PhotonPropogationICARUS") << "Looked at " << fNumEvent << " events" << std::endl;
}
