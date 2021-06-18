#ifndef IRawDigitHistogramTool_H
#define IRawDigitHistogramTool_H
////////////////////////////////////////////////////////////////////////
//
// Class:       IRawDigitHistogramTool
// Module Type: tool
// File:        IRawDigitHistogramTool.h
//
//              This provides an interface for tools which do histogramming
//              of various quantities associated to recob::Hit objects
//
// Created by Tracy Usher (usher@slac.stanford.edu) on October 16, 2017
//
////////////////////////////////////////////////////////////////////////

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/Simulation/SimChannel.h"
namespace detinfo {
  class DetectorClocksData;
  class DetectorPropertiesData;
}

class IRawDigitHistogramTool
{
public:
    /**
     *  @brief  Virtual Destructor
     */
    virtual ~IRawDigitHistogramTool() noexcept = default;
    
    /**
     *  @brief Interface for configuring the particular algorithm tool
     *
     *  @param ParameterSet  The input set of parameters for configuration
     */
    virtual void configure(const fhicl::ParameterSet&) = 0;

    /**
     *  @brief Interface for initializing the histograms to be filled
     *
     *  @param TFileService   handle to the TFile service
     *  @param string         subdirectory to store the hists in
     */
    virtual void initializeHists(detinfo::DetectorClocksData const& clockData,
                                 detinfo::DetectorPropertiesData const& detProp,
                                 art::ServiceHandle<art::TFileService>&, const std::string&) = 0;

    /**
     *  @brief Interface for method to executve at the end of run processing
     *
     *  @param int            number of events to use for normalization
     */
    virtual void endJob(int numEvents) = 0;
    
    /**
     *  @brief Interface for filling histograms
     */
    using RawDigitPtrVec = std::vector<art::Ptr<raw::RawDigit>>;
    using SimChannelMap  = std::map<raw::ChannelID_t, const sim::SimChannel*>;
    
    virtual void fillHistograms(const detinfo::DetectorClocksData&,
                                const detinfo::DetectorPropertiesData&,
                                const RawDigitPtrVec&, const SimChannelMap&)  const = 0;
};

#endif
