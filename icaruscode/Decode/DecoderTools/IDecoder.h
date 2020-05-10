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
#include "lardataobj/RawData/RawDigit.h"

//------------------------------------------------------------------------------------------------------------------------------------------

namespace art
{
    class ProducesCollector;
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
     *  @brief Initialize any data products the tool will output
     *
     */
    virtual void initializeDataProducts() = 0;

    /**
     *  @brief Given a set of recob hits, run DBscan to form 3D clusters
     *
     *  @param fragment            The artdaq fragment to process
     *  @param rawDigitColllection The output RawDigits
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
