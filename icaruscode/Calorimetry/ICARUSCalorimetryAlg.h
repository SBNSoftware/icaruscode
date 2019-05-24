////////////////////////////////////////////////////////////////////////
// \file ICARUSCalorimetryAlg.h
//
// \brief Functions to calculate dE/dx. Based on code in Calorimetry.cxx
//
// \author afava@fnal.gov
//
////////////////////////////////////////////////////////////////////////
#ifndef UTIL_ICARUSCALORIMETRYALG_H
#define UTIL_ICARUSCALORIMETRYALG_H

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"

#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include <vector>

namespace recob { 
  class Hit; 
}


///General LArSoft Utilities
namespace calo{
    class ICARUSCalorimetryAlg {
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
                        Comment("0 = exponential, 1 = exponential + constant, 2 = manual")
		};

		fhicl::Atom< bool > CaloDoLifeTimeCorrection {
			Name("CaloDoLifeTimeCorrection"),
			Comment("Apply lifetime correction if true")
		};
            
                fhicl::Atom< double > CaloManualLifeTime {
                        Name("CaloManualLifeTime"),
                        Comment("force correction to a given lifetime setting")
                };

    };

	ICARUSCalorimetryAlg(const fhicl::ParameterSet& pset) :
		ICARUSCalorimetryAlg(fhicl::Table<Config>(pset, {})())
	{}

    ICARUSCalorimetryAlg(const Config& config);
    
    ~ICARUSCalorimetryAlg();
      
    void   reconfigure(const Config& config);
    void   reconfigure(const fhicl::ParameterSet& pset)
      { reconfigure(fhicl::Table<Config>(pset, {})()); }
    
    double dEdx_AMP(art::Ptr< recob::Hit >  hit, double pitch, double T0=0) const;
    double dEdx_AMP(recob::Hit const&  hit, double pitch, double T0=0) const;
    double dEdx_AMP(double dQ, double time, double pitch, unsigned int plane, double T0=0) const;
    double dEdx_AMP(double dQdx,double time, unsigned int plane, double T0=0) const;
    
    double dEdx_AREA(art::Ptr< recob::Hit >  hit, double pitch, double T0=0) const;
    double dEdx_AREA(recob::Hit const&  hit, double pitch, double T0=0) const;
    double dEdx_AREA(double dQ,double time, double pitch, unsigned int plane, double T0=0) const;
    double dEdx_AREA(double dQdx,double time, unsigned int plane, double T0=0) const;
      
    double ElectronsFromADCPeak(double adc, unsigned short plane) const
    { return adc / fCalAmpConstants[plane]; }
      
    double ElectronsFromADCArea(double area, unsigned short plane) const
    { return area / fCalAreaConstants[plane]; }
    
    double LifetimeCorrection(double time, double T0=0) const;
    
  private:

    art::ServiceHandle<geo::Geometry> geom; 
    const detinfo::DetectorProperties* detprop;

    double dEdx_from_dQdx_e(double dQdx_e,double time, double T0=0) const;
   
    std::vector< double > fCalAmpConstants;
    std::vector< double > fCalAreaConstants;
    bool fUseModBox;
    int  fLifeTimeForm;
    bool fDoLifeTimeCorrection;
    double fManualLifeTime;
    
    }; // class ICARUSCalorimetryAlg
} //namespace calo
#endif // UTIL_ICARUSCALORIMETRYALG_H
