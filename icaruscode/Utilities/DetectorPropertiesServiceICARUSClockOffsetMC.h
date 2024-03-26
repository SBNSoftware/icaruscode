////////////////////////////////////////////////////////////////////////
// DetectorPropertiesServiceICARUSClockOffsetMC.h
//
// Service interface for DetectorProperties functions
// -- Copied from Standard but pointing to the ICARUS version with clock offsets from the trigger added.
//
// Bruce Howard heavily copied from Jon Paley jpaley@fnal.gov
//
////////////////////////////////////////////////////////////////////////
#ifndef DETECTORPROPERTIESSERVICESICARUSCLOCKOFFSETMC_H
#define DETECTORPROPERTIESSERVICESICARUSCLOCKOFFSETMC_H

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "icaruscode/Utilities/DetectorPropertiesICARUSClockOffsetMC.h"

#include "art/Framework/Principal/Run.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"
#include "art/Persistency/Provenance/ScheduleContext.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/Name.h"

/// General LArSoft Utilities
namespace detinfo {

  /**
   * Configuration parameters
   * -------------------------
   *
   * This service passes the whole configuration down to its service provider,
   * but it also reacts to:
   * - *InheritNumberTimeSamples* (boolean; default: false): if true, the
   *   configuration database in the ROOT input file is queried and if a
   *   configuration for this service is found, it's used instead of the
   *   one from the current FHiCL configuration
   *
   */

  class DetectorPropertiesServiceICARUSClockOffsetMC : public DetectorPropertiesService {

  public:
    // the following is currently not used for validation,
    // but only for documentation
    struct ServiceConfiguration_t {

      // service-specific configuration
      fhicl::Atom<bool> InheritNumberTimeSamples{
        fhicl::Name("InheritNumberTimeSamples"),
        fhicl::Comment(""),
        false /* default value */
      };

      // provider configuration
      detinfo::DetectorPropertiesICARUSClockOffsetMC::Configuration_t ProviderConfiguration;

    }; // ServiceConfiguration_t

    // this enables art to print the configuration help:
    using Parameters = art::ServiceTable<ServiceConfiguration_t>;

    DetectorPropertiesServiceICARUSClockOffsetMC(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);

  private:
    DetectorPropertiesData getDataForJob(DetectorClocksData const& clockData) const override
    {
      return fProp.DataFor(clockData);
    }

    DetectorPropertiesData getDataFor(art::Event const&,
                                      DetectorClocksData const& clockData) const override
    {
      return fProp.DataFor(clockData);
    }

    void postOpenFile(const std::string& filename);

    DetectorPropertiesICARUSClockOffsetMC fProp;
    fhicl::ParameterSet fPS; ///< Original parameter set.

    bool fInheritNumberTimeSamples; ///< Flag saying whether to inherit NumberTimeSamples

    bool isDetectorPropertiesServiceICARUSClockOffsetMC(const fhicl::ParameterSet& ps) const;

  }; // class DetectorPropertiesServiceICARUSClockOffsetMC
} // namespace detinfo

DECLARE_ART_SERVICE_INTERFACE_IMPL(detinfo::DetectorPropertiesServiceICARUSClockOffsetMC,
                                   detinfo::DetectorPropertiesService,
                                   SHARED)
#endif // DETECTORPROPERTIESSERVICESICARUSCLOCKOFFSETMC_H
