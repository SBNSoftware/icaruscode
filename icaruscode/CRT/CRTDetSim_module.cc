///////////////////////////////////////////////////////////////////////////////
/// Class: CRTDetSim
/// Module Type: producer
/// \file CRTDetSim_module.cc
///
/// Based on LArIAT TOFSimDigits.cc (Author: Lucas Mendes Santos)
/// with modifications for SBND (Author: mastbaum@uchicago.edu)
/// then modified for ICARUS
///
/// Author: Chris.Hilgenberg@colostate.edu
///
/// Extracts position, time, and energy deposited, aint with
///   trackID associated with deposit from AuxDetSimChannel objects
///   produced by G4 stage. The trackID is used for truth matching.
/// Simulated CRT data products are intended to match the known
///   front-end electronics behavior and basic DAQ format anticipated.
///   Each hit strip is assigned a channel and associated with
//      1) T0 - hit time relative to t0 (beam signal)
//      2) T1 - hit time relative to PPS from GPS
//      3) ADC from hit
///   Effects included are
///     * time
///         - light propegation delay in scintillator and WLS fiber
///         - smearing in the arrival time at the SiPM due to scattering
///         - amplitude dependant smearing due to the discriminator
///     * position
///         - AuxDetSimChannel provides true entry and exit point in scintillator strip
///         - take arithmetic mean as true "hit" position
///         - "hit" position defines transverse distance (hit to fiber) and
///           intitudinal distance (length of WLS fiber between fiber entry and SiPM)
///     * energy deposited
///         - take true energy deposited in scintillator strip
///         - convert true energy to photons yielded (taken from measurements on
///               normally incident MIP muons
///         - using known attenuation lengths for bulk scinillator and WLS fiber,
///               apply attenuation correction using transverse and intitudinal
///               propegation distances respectively and a simple exponential model
///         - for the attenuated light arriving at the SiPM, correct for counting
///               statistics sampling from Poisson
///         - convert N photons seen by SiPM into ADCs using known pedestal and
///               gain values for each channel (assumed all the same for now)
///               ADC_i = gain_i * nphotons + ped_i
///         - smear ADC using guasian with ADC_i as mean and
///               width = width_pedestal+sqrt(nphotons)
///     * front-end electonics deadtime
///         - sort vector of data products by time (T0)
///         - for C and D type modules, intermodule coincidence is applied
///            - one strip from each of 2 layers must be above threshold
///         - for M modules, only one channel is required to be above threshold
///         - the earliest channel that was part of the trigger provides the T0
///            defining the readout window (set in FHiCL)
///         - for each channel in the readout window with an entry above threshold,
///            the data (charge and time) is added the "track and hold list"
///         - channels in track and hold do not accept new data during the readout window
///         - after the readout window has passed, no channels can receive data until the
///           (FHiCL congifurable) deadtime has passed. Channels are then reset and
///           the front-end board is again able to trigger
///
/////////////////////////////////////////////////////////////////////////////////////////

//art includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "nurandom/RandomUtils/NuRandomService.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Persistency/Common/PtrMaker.h"

//larsoft includes
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataalg/DetectorInfo/ElecClock.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/AuxDetGeo.h"
#include "larcore/Geometry/AuxDetGeometry.h"
#include "larcorealg/CoreUtils/NumericUtils.h"

//ROOT includes
#include "TFile.h"
#include "TNtuple.h"

//C++ includes
#include <cmath>
#include <memory>
#include <string>
#include <utility>
#include <map>
#include <vector>

//CRT includes
#include "icaruscode/CRT/CRTProducts/CRTData.hh"
#include "icaruscode/CRT/CRTUtils/CRTCommonUtils.h"
#include "icaruscode/CRT/CRTUtils/CRTDetSimAlg.h"

using std::string;

namespace icarus {
 namespace crt {
    class CRTDetSim;
 }
}

class icarus::crt::CRTDetSim : public art::EDProducer {

 public:
    explicit CRTDetSim(fhicl::ParameterSet const & p);

    CRTDetSim(CRTDetSim const &) = delete;
    CRTDetSim(CRTDetSim &&) = delete;
    CRTDetSim& operator = (CRTDetSim const &) = delete;
    CRTDetSim& operator = (CRTDetSim &&) = delete;
    void reconfigure(fhicl::ParameterSet const & p);

    void produce(art::Event & e) override;
    string fG4ModuleLabel;

 private:

    CLHEP::HepRandomEngine& fRandEngine;
    CRTDetSimAlg detAlg;

}; //class CRTDetSim

//getting parameter values from FHiCL
void icarus::crt::CRTDetSim::reconfigure(fhicl::ParameterSet const & p) {
    fG4ModuleLabel = p.get<std::string>("G4ModuleLabel");
}//CRTDetSim::reconfigure()

// constructor
icarus::crt::CRTDetSim::CRTDetSim(fhicl::ParameterSet const & p) : EDProducer{p},
    fRandEngine(art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this, "HepJamesRandom", "crt", p, "Seed")),
    detAlg(p.get<fhicl::ParameterSet>("DetSimAlg"),fRandEngine)
{
    this->reconfigure(p);

    produces< std::vector<icarus::crt::CRTData> >();
    produces<std::vector<sim::AuxDetIDE> >();
    produces< art::Assns<icarus::crt::CRTData, sim::AuxDetIDE> >();
}

//module producer
void icarus::crt::CRTDetSim::produce(art::Event& event) {

    //pointer to vector of CRT data products to be pushed to event
    std::unique_ptr<std::vector<CRTData> > triggeredCRTHits (
        new std::vector<CRTData>);
    //pointer to vector of AuxDetIDE products to be pushed to event
    std::unique_ptr<std::vector<sim::AuxDetIDE> > ides (
        new std::vector<sim::AuxDetIDE>);
    //pointer associations between CRT data products and AuxDetIDEs to be pushed to event
    std::unique_ptr< art::Assns<CRTData,sim::AuxDetIDE> > dataAssn (
        new art::Assns<CRTData,sim::AuxDetIDE> );
    art::PtrMaker<CRTData> makeDataPtr(event);
    art::PtrMaker<sim::AuxDetIDE> makeIDEPtr(event);

    //clear detAlg member data
    detAlg.ClearTaggers();

    // Handle for (truth) AuxDetSimChannels
    art::Handle<std::vector<sim::AuxDetSimChannel> > adChanHandle;;
    std::vector< art::Ptr<sim::AuxDetSimChannel> > adChanList;
    if (event.getByLabel(fG4ModuleLabel, adChanHandle) )
        art::fill_ptr_vector(adChanList, adChanHandle);

    mf::LogInfo("CRTDetSimProducer")
        <<"Number of AuxDetChannels = " << adChanList.size();

    int nide=0;
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(event);
    for(auto const& adsc : adChanList) {

        auto const& auxDetIDEs = adsc->AuxDetIDEs();
        if(auxDetIDEs.size()>0) {
            nide++;
            detAlg.FillTaggers(clockData, adsc->AuxDetID(), adsc->AuxDetSensitiveID(), auxDetIDEs);
        }

    }//loop over AuxDetSimChannels

    //generate CRTData products, associates from filled detAlg member data
    int nData = 0;
    if(nide>0) {
        std::vector<std::pair<CRTData,std::vector<sim::AuxDetIDE>>> data = detAlg.CreateData();

        for(auto const& dataPair : data){

          triggeredCRTHits->push_back(dataPair.first);
          art::Ptr<CRTData> dataPtr = makeDataPtr(triggeredCRTHits->size()-1);
          nData++;

          for(auto const& ide : dataPair.second){
            ides->push_back(ide);
            art::Ptr<sim::AuxDetIDE> idePtr = makeIDEPtr(ides->size()-1);
            dataAssn->addSingle(dataPtr, idePtr);
          }
        }

        std::sort(triggeredCRTHits->begin(),triggeredCRTHits->end(),
            [](const CRTData& d1, const CRTData& d2) {
                return d1.fTs0 > d2.fTs0;
            });
    }

    if(nData==0)
        mf::LogWarning("CRTDetSimProducer")
          << "0 CRTData produced (expected for most neutrino events, never for cosmics)";

    //push products to event
    event.put(std::move(triggeredCRTHits));
    event.put(std::move(ides));
    event.put(std::move(dataAssn));

    mf::LogInfo("CRTDetSimProducer")
      << "Number of CRT data produced = "<< nData << '\n'
      << "  from " << nide << " AuxDetIDEs";
}

DEFINE_ART_MODULE(icarus::crt::CRTDetSim)
