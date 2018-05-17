#ifndef IWireHistogramTool_H
#define IWireHistogramTool_H
////////////////////////////////////////////////////////////////////////
//
// Class:       IWireHistogramTool
// Module Type: tool
// File:        IWireHistogramTool.h
//
//              This provides an interface for tools which do histogramming
//              of various quantities associated to recob::Hit objects
//
// Created by Tracy Usher (usher@slac.stanford.edu) on October 16, 2017
//
////////////////////////////////////////////////////////////////////////

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/Simulation/SimChannel.h"

class IWireHistogramTool
{
public:
    /**
     *  @brief  Virtual Destructor
     */
    virtual ~IWireHistogramTool() noexcept = default;
    
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
    virtual void initializeHists(art::ServiceHandle<art::TFileService>&, const std::string&) = 0;

    /**
     *  @brief Interface for method to executve at the end of run processing
     *
     *  @param int            number of events to use for normalization
     */
    virtual void endJob(int numEvents) = 0;
    
    /**
     *  @brief Interface for filling histograms
     */
    using WirePtrVec     = std::vector<art::Ptr<recob::Wire>>;
    using SimChannelMap  = std::map<raw::ChannelID_t, const sim::SimChannel*>;
    
    virtual void fillHistograms(const WirePtrVec&, const SimChannelMap&, int)  const = 0;
};

#endif
