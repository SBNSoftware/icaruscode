/**
 *  @file   IDecoderFilter.h
 *
 *  @brief  This provides an art tool interface definition for tools which "decode" artdaq 
 *          fragments into LArSoft data objects
 *
 *  @author usher@slac.stanford.edu
 *
 */
#ifndef IDecoderFilter_h
#define IDecoderFilter_h

// Framework Includes
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Event.h"

// Algorithm includes
#include "artdaq-core/Data/Fragment.hh"

#include "icarussigproc/ICARUSSigProcDefs.h"

//------------------------------------------------------------------------------------------------------------------------------------------

namespace daq
{
/**
 *  @brief  IDecoderFilter interface class definiton
 */
class IDecoderFilter
{
public:
    /**
     *  @brief  Virtual Destructor
     */
    virtual ~IDecoderFilter() noexcept = default;

    /**
     *  @brief Interface for configuring the particular algorithm tool
     *
     *  @param ParameterSet  The input set of parameters for configuration
     */
    virtual void configure(const fhicl::ParameterSet&) = 0;

    /**
     *  @brief Given a set of recob hits, run DBscan to form 3D clusters
     *
     *  @param fragment            The artdaq fragment to process
     */
   virtual void process_fragment(const artdaq::Fragment& fragment) = 0; 

    /**
     *  @brief Recover the selection values
     */
    virtual const icarussigproc::ArrayBool  getSelectionVals()     const = 0;

    /**
     *  @brief Recover the ROI values
     */
    virtual const icarussigproc::ArrayBool  getROIVals()           const = 0;

    /**
     *  @brief Recover the original raw waveforms
     */
    virtual const icarussigproc::ArrayFloat getRawWaveforms() const = 0;

    /**
     *  @brief Recover the pedestal corrected waveforms
     */
    virtual const icarussigproc::ArrayFloat getPedCorWaveforms() const = 0;

    /**
     *  @brief Recover the "intrinsic" RMS
     */
    virtual const icarussigproc::ArrayFloat getIntrinsicRMS()      const = 0;

    /**
     *  @brief Recover the correction median values
     */
    virtual const icarussigproc::ArrayFloat getCorrectedMedians()  const = 0;

    /**
     *  @brief Recover the waveforms less coherent noise
     */
    virtual const icarussigproc::ArrayFloat getWaveLessCoherent()  const = 0;

    /**
     *  @brief Recover the morphological filter waveforms
     */
    virtual const icarussigproc::ArrayFloat getMorphedWaveforms()  const = 0;

    /**
     *  @brief Recover the pedestals for each channel
     */
    virtual const icarussigproc::VectorFloat getPedestalVals() const = 0;

    /**
     *  @brief Recover the full RMS before coherent noise
     */
    virtual const icarussigproc::VectorFloat getFullRMSVals()  const = 0;
 
    /**
     *  @brief Recover the truncated RMS noise 
     */
    virtual const icarussigproc::VectorFloat getTruncRMSVals() const = 0;

    /**
     *  @brief Recover the number of bins after truncation
     */
    virtual const icarussigproc::VectorInt   getNumTruncBins() const = 0;
 
};

} // namespace lar_cluster3d
#endif
