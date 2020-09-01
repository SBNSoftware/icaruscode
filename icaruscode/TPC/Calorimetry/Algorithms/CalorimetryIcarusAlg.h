////////////////////////////////////////////////////////////////////////
// \file CalorimetryIcarusAlg.h
//
// \brief Functions to calculate dE/dx. Based on code in Calorimetry.cxx
//
// \author andrzej.szelc@yale.edu
//
////////////////////////////////////////////////////////////////////////
#ifndef UTIL_CalorimetryIcarusAlg_H
#define UTIL_CalorimetryIcarusAlg_H

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "larcore/Geometry/Geometry.h"
#include <vector>

namespace detinfo {
  class DetectorClocksData;
  class DetectorPropertiesData;
}

namespace recob {
  class Hit;
}


///General LArSoft Utilities
namespace calo{
    class CalorimetryIcarusAlg {
    public:

	struct Config {
		using Name = fhicl::Name;
		using Comment = fhicl::Comment;

		fhicl::Sequence< double > CalAmpConstants {
			Name("CalAmpConstants"),
			Comment("ADC to electrons constants for each plane.")
		};

		fhicl::Sequence< double > CalAreaConstants {
			Name("CalAreaConstants"),
			Comment("Area to electrons constants for each plane.")
		};

		fhicl::Atom< bool > CaloUseModBox {
			Name("CaloUseModBox"),
			Comment("Use modified box model if true, birks otherwise")
		};

                fhicl::Atom< int > CaloLifeTimeForm {
                        Name("CaloLifeTimeForm"),
                        Comment("0 = exponential, 1 = exponential + constant")
		};

		fhicl::Atom< bool > CaloDoLifeTimeCorrection {
			Name("CaloDoLifeTimeCorrection"),
			Comment("Apply lifetime correction if true")
		};

    };

	CalorimetryIcarusAlg(const fhicl::ParameterSet& pset) :
		CalorimetryIcarusAlg(fhicl::Table<Config>(pset, {})())
	{}

    CalorimetryIcarusAlg(const Config& config);

    ~CalorimetryIcarusAlg();

    void   reconfigure(const Config& config);
    void   reconfigure(const fhicl::ParameterSet& pset)
      { reconfigure(fhicl::Table<Config>(pset, {})()); }

    double dEdx_AMP(detinfo::DetectorClocksData const& clockData,
                    detinfo::DetectorPropertiesData const& detProp,
                    art::Ptr< recob::Hit >  hit, double pitch, double T0=0) const;
    double dEdx_AMP(detinfo::DetectorClocksData const& clockData,
                    detinfo::DetectorPropertiesData const& detProp,
                    recob::Hit const&  hit, double pitch, double T0=0) const;
    double dEdx_AMP(detinfo::DetectorClocksData const& clockData,
                    detinfo::DetectorPropertiesData const& detProp,
                    double dQ, double time, double pitch, unsigned int plane, double T0=0) const;
    double dEdx_AMP(detinfo::DetectorClocksData const& clockData,
                    detinfo::DetectorPropertiesData const& detProp,
                    double dQdx,double time, unsigned int plane, double T0=0) const;

    double dEdx_AREA(detinfo::DetectorClocksData const& clockData,
                     detinfo::DetectorPropertiesData const& detProp,
                     art::Ptr< recob::Hit >  hit, double pitch, double T0=0) const;
    double dEdx_AREA(detinfo::DetectorClocksData const& clockData,
                     detinfo::DetectorPropertiesData const& detProp,
                     recob::Hit const&  hit, double pitch, double T0=0) const;
    double dEdx_AREA(detinfo::DetectorClocksData const& clockData,
                     detinfo::DetectorPropertiesData const& detProp,
                     double dQ,double time, double pitch, unsigned int plane, double T0=0) const;
    double dEdx_AREA(detinfo::DetectorClocksData const& clockData,
                     detinfo::DetectorPropertiesData const& detProp,
                     double dQdx,double time, unsigned int plane, double T0=0) const;

    double dEdx_SUMADC(detinfo::DetectorClocksData const& clockData,
                       detinfo::DetectorPropertiesData const& detProp,
                       art::Ptr< recob::Hit >  hit, double pitch, double T0=0) const;
    double dEdx_SUMADC(detinfo::DetectorClocksData const& clockData,
                       detinfo::DetectorPropertiesData const& detProp,
                       recob::Hit const&  hit, double pitch, double T0=0) const;
    double ElectronsFromADCPeak(double adc, unsigned short plane) const
    { return adc / fCalAmpConstants[plane]; }

    double ElectronsFromADCArea(double area, unsigned short plane) const
    { return area / fCalAreaConstants[plane]; }

    double LifetimeCorrection(detinfo::DetectorClocksData const& clockData,
                              detinfo::DetectorPropertiesData const& detProp,
                              double time, double T0=0) const;

  private:

    art::ServiceHandle<geo::Geometry const> geom;

    double dEdx_from_dQdx_e(detinfo::DetectorClocksData const& clockData,
                            detinfo::DetectorPropertiesData const& detProp,
                            double dQdx_e,double time, double T0=0) const;

    std::vector< double > fCalAmpConstants;
    std::vector< double > fCalAreaConstants;
    bool fUseModBox;
    int  fLifeTimeForm;
    bool fDoLifeTimeCorrection;

    }; // class CalorimetryIcarusAlg
} //namespace calo
#endif // UTIL_CalorimetryIcarusAlg_H
