// RawDigitAna_module.cc
//
// The aim of this ana module is to do some basic analysis of RawDigit
// waveforms to better understand issues...
//

#ifndef RawDigitAna_module
#define RawDigitAna_module

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/RawData/RawDigit.h"
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
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib_except/exception.h"

#include "icaruscode/Analysis/tools/IRawDigitHistogramTool.h"

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

namespace RawDigitAna
{
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// class definition

class RawDigitAna : public art::EDAnalyzer
{
public:

    // Standard constructor and destructor for an ART module.
    explicit RawDigitAna(fhicl::ParameterSet const& pset);
    virtual ~RawDigitAna();

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
    std::vector<art::InputTag> fRawDigitProducerLabelVec;
    art::InputTag              fSimChannelProducerLabel;

    // The variables that will go into the n-tuple.
    int fEvent;
    int fRun;
    int fSubRun;
    int fNumEvents;
    
    // Keep track of the hit histogramming tools here
    std::vector<std::unique_ptr<IRawDigitHistogramTool>> fRawDigitHistogramToolVec;

    std::vector<std::vector<double>> fChannelPedVec;

    // Other variables that will be shared between different methods.
    const geo::GeometryCore*           fGeometry;       // pointer to Geometry service
    const detinfo::DetectorProperties* fDetectorProperties;
    const lariov::DetPedestalProvider& fPedestalRetrievalAlg; ///< Keep track of an instance to the pedestal retrieval alg
}; // class RawDigitAna


//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// class implementation

//-----------------------------------------------------------------------
// Constructor
RawDigitAna::RawDigitAna(fhicl::ParameterSet const& parameterSet)
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
RawDigitAna::~RawDigitAna()
{}

//-----------------------------------------------------------------------
void RawDigitAna::beginJob()
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
    for (auto& rawDigitHistTool : fRawDigitHistogramToolVec) rawDigitHistTool->initializeHists(tfs, "RawDigitAna");

    // zero out the event counter
    fNumEvents = 0;
}

//-----------------------------------------------------------------------
void RawDigitAna::beginRun(const art::Run& /*run*/)
{
    // How to convert from number of electrons to GeV.  The ultimate
    // source of this conversion factor is
    // ${LARSIM_DIR}/include/SimpleTypesAndConstants/PhysicalConstants.h.
//    art::ServiceHandle<sim::LArG4Parameters> larParameters;
//    fElectronsToGeV = 1./larParameters->GeVToElectrons();
}

//-----------------------------------------------------------------------
void RawDigitAna::reconfigure(fhicl::ParameterSet const& p)
{
    // Read parameters from the .fcl file. The names in the arguments
    // to p.get<TYPE> must match names in the .fcl file.
    fRawDigitProducerLabelVec = p.get< std::vector<art::InputTag> >("RawDigitModuleLabel",   std::vector<art::InputTag>() = {"rawdigitfilter"});
    fSimChannelProducerLabel  = p.get< std::string                >("SimChannelModuleLabel", "largeant"      );

    // Implement the tools for handling the responses
    const std::vector<fhicl::ParameterSet>& rawDigitHistogramToolVec = p.get<std::vector<fhicl::ParameterSet>>("RawDigitHistogramToolList");
    
    for(auto& rawDigitHistogramTool : rawDigitHistogramToolVec)
        fRawDigitHistogramToolVec.push_back(art::make_tool<IRawDigitHistogramTool>(rawDigitHistogramTool));

    return;
}

//-----------------------------------------------------------------------
void RawDigitAna::analyze(const art::Event& event)
{
    // Start by fetching some basic event information for our n-tuple.
    fEvent  = event.id().event();
    fRun    = event.run();
    fSubRun = event.subRun();

    fNumEvents++;
    
    // Loop over RawDigits
    for(const auto& rawDigitLabel : fRawDigitProducerLabelVec)
    {
        // Make a pass through all hits to make contrasting plots
        art::Handle< std::vector<raw::RawDigit> > rawDigitHandle;
        event.getByLabel(rawDigitLabel, rawDigitHandle);
        
        // Recover sim channels (if they exist) so we know when a
        // channel has signal (or not)
        art::Handle<std::vector<sim::SimChannel>>  simChannelHandle;
        event.getByLabel(fSimChannelProducerLabel, simChannelHandle);
        
        // Recover list of simChannels mapped by channel to make
        // look up easier below
        IRawDigitHistogramTool::SimChannelMap channelMap;
        
        if (simChannelHandle.isValid())
        {
            for(const auto& simChannel : *simChannelHandle) channelMap[simChannel.Channel()] = &simChannel;
        }
        
        if (rawDigitHandle.isValid())
        {
            IRawDigitHistogramTool::RawDigitPtrVec allRawDigitVec;
            art::fill_ptr_vector(allRawDigitVec, rawDigitHandle);
            
            for(auto& rawDigitHistTool : fRawDigitHistogramToolVec) rawDigitHistTool->fillHistograms(allRawDigitVec,channelMap);
        }
    }

    return;
}

void RawDigitAna::endJob()
{
    // Make a call to normalize histograms if so desired
    for(auto& rawDigitHistTool : fRawDigitHistogramToolVec) rawDigitHistTool->endJob(fNumEvents);

    return;
}

// This macro has to be defined for this module to be invoked from a
// .fcl file; see RawDigitAna.fcl for more information.
DEFINE_ART_MODULE(RawDigitAna)

} // namespace RawDigitAna

#endif // RawDigitAna_module
