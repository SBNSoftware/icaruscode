/**
 *  @file   IFakeParticle.h
 *
 *  @brief  This provides an art tool interface definition for tools which can create
 *          "fake" particles to overlay onto input daq fragments during decoding
 *
 *  @author usher@slac.stanford.edu
 *
 */
#ifndef IFakeParticle_h
#define IFakeParticle_h

// Framework Includes
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Event.h"

// Algorithm includes
#include "artdaq-core/Data/Fragment.hh"
namespace detinfo {
  class DetectorClocksData;
}

//------------------------------------------------------------------------------------------------------------------------------------------

namespace daq
{
/**
 *  @brief  IFakeParticle interface class definiton
 */
class IFakeParticle
{
public:
    /**
     *  @brief  Virtual Destructor
     */
    virtual ~IFakeParticle() noexcept = default;

    /**
     *  @brief Interface for configuring the particular algorithm tool
     *
     *  @param ParameterSet  The input set of parameters for configuration
     */
    virtual void configure(const fhicl::ParameterSet&) = 0;

    /**
     *  @brief Creates a fake particle and overlays on the input fragment
     *
     *  @param waveforms  The waveform container to place fake particle on
     */
    using VectorShort  = std::vector<short>;
    using VectorFloat  = std::vector<float>;
    using ArrayShort   = std::vector<VectorShort>;
    using ArrayFloat   = std::vector<VectorFloat>;
    virtual void overlayFakeParticle(detinfo::DetectorClocksData const& clockData,
                                     ArrayFloat& waveforms) = 0;
};

} // namespace daq
#endif
