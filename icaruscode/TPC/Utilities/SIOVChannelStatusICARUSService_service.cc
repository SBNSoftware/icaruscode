#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
#include "art/Persistency/Provenance/ScheduleContext.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "larcore/CoreUtils/EnsureOnlyOneSchedule.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Providers/SIOVChannelStatusProvider.h"

namespace lariov {

  /**
     \class SIOVChannelStatusICARUSService
     art service implementation of ChannelStatusService.  Implements
     a channel status retrieval service for database scheme in which
     all elements in a database folder share a common interval of validity
  */
  class SIOVChannelStatusICARUSService : public ChannelStatusService,
                                   private lar::EnsureOnlyOneSchedule<SIOVChannelStatusICARUSService> {

  public:
    SIOVChannelStatusICARUSService(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);

    void PreProcessEvent(const art::Event& evt, art::ScheduleContext);

  private:
    const ChannelStatusProvider& DoGetProvider() const override { return fProvider; }

    const ChannelStatusProvider* DoGetProviderPtr() const override { return &fProvider; }

    SIOVChannelStatusProvider fProvider;
  };
} //end namespace lariov

DECLARE_ART_SERVICE_INTERFACE_IMPL(lariov::SIOVChannelStatusICARUSService,
                                   lariov::ChannelStatusService,
                                   SHARED)

namespace lariov {

  SIOVChannelStatusICARUSService::SIOVChannelStatusICARUSService(fhicl::ParameterSet const& pset,
                                                     art::ActivityRegistry& reg)
    : fProvider(pset.get<fhicl::ParameterSet>("ChannelStatusProvider"))
  {

    //register callback to update local database cache before each event is processed
    reg.sPreProcessEvent.watch(this, &SIOVChannelStatusICARUSService::PreProcessEvent);
  }

  void SIOVChannelStatusICARUSService::PreProcessEvent(const art::Event& evt, art::ScheduleContext)
  {
    std::uint64_t timeStamp = std::uint64_t(evt.time().timeHigh()) * 1000000000 + std::uint64_t(evt.time().timeLow());

    mf::LogInfo("SIOVChannelStatusICARUSService") << "==> PreProcessEvent using timestamp (ns): " << timeStamp << ", timeHigh: " << evt.time().timeHigh() << ", timeLow: " << evt.time().timeLow();

    //First grab an update from the database
    fProvider.UpdateTimeStamp(timeStamp);
  }

} //end namespace lariov

DEFINE_ART_SERVICE_INTERFACE_IMPL(lariov::SIOVChannelStatusICARUSService, lariov::ChannelStatusService)
