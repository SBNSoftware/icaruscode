// WireAna_module.cc
//
// The aim of this ana module is to do some basic analysis of RawDigit
// waveforms to better understand issues...
//

#ifndef WireAna_module
#define WireAna_module

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RawData/raw.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalService.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalProvider.h"

//#include "cetlib/search_path.h"
#include "cetlib/cpu_timer.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/View.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib_except/exception.h"

#include "icaruscode/Analysis/tools/IWireHistogramTool.h"

// C++ Includes
#include <map>
#include <vector>
#include <tuple>
#include <algorithm>
#include <iostream>
#include <string>
#include <cmath>

#include <iostream>
#include <fstream>

namespace WireAna
{
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// class definition

class WireAna : public art::EDAnalyzer
{
public:

    // Standard constructor and destructor for an ART module.
    explicit WireAna(fhicl::ParameterSet const& pset);
    virtual ~WireAna();

    // This method is called once, at the start of the job. In this
    // example, it will define the histograms and n-tuples we'll write.
    void beginJob() override;
    void endJob()   override;

    // This method is called once, at the start of each run. It's a
    // good place to read databases or files that may have
    // run-dependent information.
    void beginRun(const art::Run& run) override;

    // This method reads in any parameters from the .fcl files. This
    // method is called 'reconfigure' because it might be called in the
    // middle of a job; e.g., if the user changes parameter values in an
    // interactive event display.
    void reconfigure(fhicl::ParameterSet const& pset);

    // The analysis routine, called once per event.
    void analyze (const art::Event& evt) override;

private:

    // The parameters we'll read from the .fcl file.
    art::InputTag fWireProducerLabel;
    art::InputTag fSimChannelProducerLabel;

    // The variables that will go into the n-tuple.
    int fEvent;
    int fRun;
    int fSubRun;
    int fNumEvents;
    
    // Keep track of the hit histogramming tools here
    std::vector<std::unique_ptr<IWireHistogramTool>> fWireHistogramToolVec;

    std::vector<std::vector<double>> fChannelPedVec;

    // Other variables that will be shared between different methods.
    const geo::GeometryCore*           fGeometry;       // pointer to Geometry service
    const detinfo::DetectorProperties* fDetectorProperties;
    const lariov::DetPedestalProvider& fPedestalRetrievalAlg; ///< Keep track of an instance to the pedestal retrieval alg
}; // class WireAna


//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// class implementation

//-----------------------------------------------------------------------
// Constructor
WireAna::WireAna(fhicl::ParameterSet const& parameterSet)
    : EDAnalyzer(parameterSet),
      fPedestalRetrievalAlg(*lar::providerFrom<lariov::DetPedestalService>())

{
    fGeometry = lar::providerFrom<geo::Geometry>();
    fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();

    // Read in the parameters from the .fcl file.
    this->reconfigure(parameterSet);
}

//-----------------------------------------------------------------------
// Destructor
WireAna::~WireAna()
{}

//-----------------------------------------------------------------------
void WireAna::beginJob()
{
    // Get the detector length, to determine the maximum bin edge of one
    // of the histograms.
    //double detectorLength = fGeometry->DetLength();

    // Access ART's TFileService, which will handle creating and writing
    // histograms and n-tuples for us.
    art::ServiceHandle<art::TFileService> tfs;

    // The arguments to 'make<whatever>' are the same as those passed
    // to the 'whatever' constructor, provided 'whatever' is a ROOT
    // class that TFileService recognizes.
    for (auto& wireHistTool : fWireHistogramToolVec) wireHistTool->initializeHists(tfs, "WireAna");

    // zero out the event counter
    fNumEvents = 0;
}

//-----------------------------------------------------------------------
void WireAna::beginRun(const art::Run& /*run*/)
{
    // How to convert from number of electrons to GeV.  The ultimate
    // source of this conversion factor is
    // ${LARSIM_DIR}/include/SimpleTypesAndConstants/PhysicalConstants.h.
//    art::ServiceHandle<sim::LArG4Parameters> larParameters;
//    fElectronsToGeV = 1./larParameters->GeVToElectrons();
}

//-----------------------------------------------------------------------
void WireAna::reconfigure(fhicl::ParameterSet const& p)
{
    // Read parameters from the .fcl file. The names in the arguments
    // to p.get<TYPE> must match names in the .fcl file.
    fWireProducerLabel       = p.get< std::string >("WireModuleLabel",       "recowire");
    fSimChannelProducerLabel = p.get< std::string >("SimChannelModuleLabel", "daq"     );

    // Implement the tools for handling the responses
    const std::vector<fhicl::ParameterSet>& wireHistogramToolVec = p.get<std::vector<fhicl::ParameterSet>>("WireHistogramToolList");
    
    for(auto& wireHistogramTool : wireHistogramToolVec)
        fWireHistogramToolVec.push_back(art::make_tool<IWireHistogramTool>(wireHistogramTool));

    return;
}

//-----------------------------------------------------------------------
void WireAna::analyze(const art::Event& event)
{
    // Start by fetching some basic event information for our n-tuple.
    fEvent  = event.id().event();
    fRun    = event.run();
    fSubRun = event.subRun();

    fNumEvents++;
    
    // Make a pass through all hits to make contrasting plots
    art::Handle< std::vector<recob::Wire> >  wireHandle;
    event.getByLabel(fWireProducerLabel, wireHandle);
    
    // Recover sim channels (if they exist) so we know when a
    // channel has signal (or not)
    art::Handle<std::vector<sim::SimChannel>>  simChannelHandle;
    event.getByLabel(fSimChannelProducerLabel, simChannelHandle);
    
    // Recover list of simChannels mapped by channel to make
    // look up easier below
    IWireHistogramTool::SimChannelMap channelMap;
    
    if (simChannelHandle.isValid())
    {
        for(const auto& simChannel : *simChannelHandle) channelMap[simChannel.Channel()] = &simChannel;
//        {
//            raw::ChannelID_t       channel       = simChannel.Channel();
//            const sim::SimChannel* simChannelPtr = &simChannel;
//
//            channelMap.at(channel) = simChannelPtr;
//        }
    }

    if (wireHandle.isValid())
    {
        IWireHistogramTool::WirePtrVec wireVec;
        art::fill_ptr_vector(wireVec, wireHandle);

        for(auto& wireHistTool : fWireHistogramToolVec) wireHistTool->fillHistograms(wireVec,channelMap,fNumEvents);
    }

    return;
}

void WireAna::endJob()
{
    // Make a call to normalize histograms if so desired
    for(auto& wireHistTool : fWireHistogramToolVec) wireHistTool->endJob(fNumEvents);

    return;
}

// This macro has to be defined for this module to be invoked from a
// .fcl file; see WireAna.fcl for more information.
DEFINE_ART_MODULE(WireAna)

} // namespace WireAna

#endif // WireAna_module
