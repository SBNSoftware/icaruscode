/**
 *  @file   IDecoder.h
 *
 *  @brief  This provides an art tool interface definition for tools which "decode" artdaq 
 *          fragments into LArSoft data objects
 *
 *  @author usher@slac.stanford.edu
 *
 */
#ifndef IDecoder_h
#define IDecoder_h

// Framework Includes
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Event.h"

// Algorithm includes
#include "artdaq-core/Data/Fragment.hh"

//------------------------------------------------------------------------------------------------------------------------------------------

namespace art
{
    class ProducesCollector;
    class ConsumesCollector;
    class Run;
}

namespace daq
{
/**
 *  @brief  IDecoder interface class definiton
 */
class IDecoder
{
public:
    /**
     *  @brief  Virtual Destructor
     */
    virtual ~IDecoder() noexcept = default;

    /**
     *  @brief Declare to the framework what you expect to read.
     */
    virtual void consumes(art::ConsumesCollector&) {}

    /**
     *  @brief The space point building should output the hit collection
     *         for those hits which combine to form space points - a nice noise filter!
     */
    virtual void produces(art::ProducesCollector&) = 0;

    /**
     *  @brief Interface for configuring the particular algorithm tool
     *
     *  @param ParameterSet  The input set of parameters for configuration
     */
    virtual void configure(const fhicl::ParameterSet&) = 0;

    /**
     *  @brief Preparation to process a new run.
     *
     *  To be called on every _art_ run transition.
     */
    virtual void setupRun(art::Run const& run) {}

    /**
     *  @brief Returns the tag of the input fragment, if known (empty otherwise).
     *
     *  The steering module can optionally use this tag for choosing or
     *  overriding the input fragment.
     */
    virtual std::optional<art::InputTag> preferredInput() const
        { return std::nullopt; }

    /**
     *  @brief Preparation to process a new event.
     *
     *  To be called on every _art_ event transition.
     */
    virtual void setupEvent(art::Event const& event) {}

    /**
     *  @brief Initialize any data products the tool will output
     *
     */
    virtual void initializeDataProducts() = 0;

    /**
     *  @brief Given a set of recob hits, run DBscan to form 3D clusters
     *
     *  @param fragment            The artdaq fragment to process
     */
    virtual void process_fragment(const artdaq::Fragment& fragment) = 0;

    /**
     *  @brief Output the data products to the event store
     * 
     *  @param event The event store objects
     */
    virtual void outputDataProducts(art::Event& event) = 0;
};

} // namespace lar_cluster3d
#endif
