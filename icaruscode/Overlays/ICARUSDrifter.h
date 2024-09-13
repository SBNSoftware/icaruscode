#ifndef ICARUS_DRIFTER_H
#define ICARUS_DRIFTER_H

#include "WireCellGen/Drifter.h"

#include "larwirecell/Interfaces/IArtEventVisitor.h"
#include "larevt/CalibrationDBI/Providers/DBFolder.h"

#include <map>
#include <string>

namespace wcls {
    class ICARUSDrifter : public WireCell::Gen::Drifter,
                          public IArtEventVisitor
    {
    public:
        ICARUSDrifter();
        virtual ~ICARUSDrifter();

        /// IArtEventVisitor.
        //
        // Note: we don't actually poke at the event but use this
        // entry to refresh info from services in case they change
        // event-to-event.
        virtual void visit(art::Event & event);

        /// IConfigurable.
        //
        // Defer default to parent.  By default this class does not
        // overriding. 
        
    // Defer to parent but override if asked to.
        virtual void configure(const WireCell::Configuration& config);

    private:
        double GetLifetime(uint64_t run, uint64_t itpc);

	std::string fDBFileName;
	std::string fDBTag;
	bool fVerbose;
	lariov::DBFolder *fDB;

        int fTPC;
        bool fELifetimeCorrection;

        std::map<std::pair<int, int>, double> fLifetimes; // cache
    };
}  // wcls
#endif 
