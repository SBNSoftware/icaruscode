#ifndef ITrackHistogramTOOL_H
#define ITrackHistogramTOOL_H
////////////////////////////////////////////////////////////////////////
//
// Class:       ITrackHistogramTool
// Module Type: tool
// File:        ITrackHistogramTool.h
//
//              This provides an interface for tools which do histogramming
//              of various quantities associated to recob::Track objects
//
// Created by Tracy Usher (usher@slac.stanford.edu) on October 16, 2017
//
////////////////////////////////////////////////////////////////////////

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "lardataobj/RecoBase/Track.h"

class ITrackHistogramTool
{
public:
    /**
     *  @brief  Virtual Destructor
     */
    virtual ~ITrackHistogramTool() noexcept = default;
    
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
    virtual void fillHistograms(const art::Event&)  const = 0;
};

#endif
