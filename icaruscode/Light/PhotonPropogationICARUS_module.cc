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


// LArSoft libraries
#include "larsim/PhotonPropagation/PhotonVisibilityService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/OpDetGeo.h"
#include "larcorealg/Geometry/geo_vectors_utils.h" // geo::vect::toPoint()
#include "lardataobj/Simulation/SimPhotons.h"
#include "nutools/RandomUtils/NuRandomService.h"

// framework libraries
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Utilities/InputTag.h"
// #include "art/Utilities/make_tool.h"
// #include "canvas/Persistency/Common/Ptr.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"

// CLHEP/ROOT libraries
// #include "CLHEP/Random/RandFlat.h"
// #include "CLHEP/Random/RandPoissonQ.h"
#include "CLHEP/Random/RandLandau.h"

// C++ libraries
// #include <cmath>
// #include <algorithm>
#include <vector>
#include <memory> // std::make_unique()


class PhotonPropogationICARUS : public art::EDProducer
{
public:

    // Copnstructors, destructor.
    explicit PhotonPropogationICARUS(fhicl::ParameterSet const & pset);

    // Overrides.
    virtual void configure(fhicl::ParameterSet const & pset);
    virtual void produce(art::Event & e) override;
    virtual void beginJob() override;
    virtual void endJob() override;

private:

    // Fcl parameters.
    art::InputTag fSimPhotonModuleLabel;      ///< The full collection of SimPhotons

    // Statistics.
    unsigned int fNumEvent = 0;        ///< Number of events seen.
    
    // Useful services, keep copies for now (we can update during begin run periods)
    geo::GeometryCore const* fGeometry = nullptr; ///< Pointer to Geometry service.
//    detinfo::DetectorProperties const* fDetectorProperties = nullptr; ///< Detector properties service.
    CLHEP::HepRandomEngine&  fPhotonEngine;
    
    /// We don't keep more than this number of photons per `sim::SimPhoton`.
    static constexpr unsigned int MaxPhotons = 10000000U;
    
}; // class PhotonPropagationICARUS

DEFINE_ART_MODULE(PhotonPropogationICARUS)

//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// pset - Fcl parameters.
///
PhotonPropogationICARUS::PhotonPropogationICARUS(fhicl::ParameterSet const & pset)
  : fGeometry(lar::providerFrom<geo::Geometry>())
  , fPhotonEngine(art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this, "HepJamesRandom", "icarusphoton", pset, "SeedPhoton"))
//  , fDetectorProperties(lar::providerFrom<detinfo::DetectorPropertiesService>())
{
    configure(pset);
    
    produces<std::vector<sim::SimPhotons>>();

    // Report.
    mf::LogDebug("PhotonPropogationICARUS") << "PhotonPropogationICARUS configured";
} // PhotonPropagationICARUS::PhotonPropagationICARUS()

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
//    art::ServiceHandle<art::TFileService> tfs;
    
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
    
    /*
     * In LArSoft, detected photons are organised optical channel by channel,
     * each channel with its own sim::SimPhotons.
     * Each channel contains individual information for each of the photons
     * it has detected.
     * Here we go through all detected photons from all the channels,
     * and for each we assign a delay to the start time, parametrized by the
     * distance between its generation point and the PMT.
     *
     */
    // Agreed convention is to ALWAYS output to the event store so get a pointer to our collection
    auto simPhotonVec = std::make_unique<std::vector<sim::SimPhotons>>();
    
    // Read in the digit List object(s).
    auto const& srcSimPhotons
      = *(event.getValidHandle< std::vector<sim::SimPhotons> >(fSimPhotonModuleLabel));
       
    if (srcSimPhotons.empty()) {
        mf::LogWarning("PhotonPropagationICARUS") << "No photons! Nice event you have here.";
        event.put(std::move(simPhotonVec));
        return;
    }

    // get hold of all needed services
//    auto const& pvs = *(art::ServiceHandle<phot::PhotonVisibilityService>());
    CLHEP::RandLandau landauGen(fPhotonEngine);

    // Loop through the input photons (this might need to be more complicated...)
    for(const auto& simPhoton : srcSimPhotons)
    {
        auto const channel = simPhoton.OpChannel();
        if (simPhoton.empty()) continue;
        
        auto const& PMTcenter = fGeometry->OpDetGeoFromOpChannel(channel).GetCenter();
       
        // TODO restore the "MF_LOG_TRACE" line for normal operations
        // (MF_LOG_TRACE will print only in debug mode, with debug qualifiers and proper messagefacility settings)
      //  MF_LOG_TRACE("LightPropagationICARUS")
//        mf::LogVerbatim("LightPropagationICARUS")
//          << "Processing photon channel #" << channel << ", detector center: " << PMTcenter;
        // 
        // fix the arrival time
        //
        sim::SimPhotons localPhoton = simPhoton; // modify a copy of the original photon
        unsigned int photonNo = 0;
        for (auto& onePhoton: localPhoton) 
        {
             //
             // sanity check: if there are too many photons we are in trouble
             // (like in: not enough computing resources)
             //
             if (photonNo++ >= MaxPhotons) {
               mf::LogError("LightPropagationICARUS")
                 << "Too many photons to handle! only " << (photonNo - 1) << " saved.";
               break;
             }
               
             geo::Point_t const position = geo::vect::toPoint(onePhoton.InitialPosition);
             double const dis = (position - PMTcenter).R();

             double mean  = 0.18*dis; // TODO
             double sigma = 0.75*dis; // TODO
             
             //double time_plus = -1;

	     //while (time_plus<dis/21.74){time_plus=mean + landauGen.fire() * sigma; }// TODO

	    const double minPropTime = dis / 21.74; // d / (c/n) in [ns]
	    //const double maxPropTime = 2000 / 21.74; // dimension / (c/n) in [ns]
     	    double time_plus;
            do {
            time_plus = mean + landauGen.fire() * sigma; // TODO
		//time_plus = landauGen.fire(mean,sigma); // TODO
            } while (time_plus < minPropTime);
             
             // TODO restore the "MF_LOG_TRACE" line for normal operations
             MF_LOG_TRACE("LightPropagationICARUS");
//             mf::LogVerbatim("LightPropagationICARUS")
//               << "Photon #" << photonNo
//               << " (at " << position << ", " << dis << " cm far from PMT) given offset "
//               << time_plus;
             //double time_plus = dis/21.74; 
             //TRandom r;
		
             onePhoton.Time += time_plus;
		//onePhoton.Time = time_plus;
             
        } // for photons in simPhoton
        // move the photon with the new arrival time into the new output collection
        simPhotonVec->push_back(std::move(localPhoton));
    } // for all simPhoton channels
    
    // Add tracks and associations to event.
    event.put(std::move(simPhotonVec));
    
} // LightPropagationICARUS::produce()

//----------------------------------------------------------------------------
/// End job method.
void PhotonPropogationICARUS::endJob()
{
    mf::LogInfo("PhotonPropogationICARUS") << "Looked at " << fNumEvent << " events" << std::endl;
}
