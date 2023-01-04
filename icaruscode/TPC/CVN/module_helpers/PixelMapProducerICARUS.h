#ifndef PixelMapProducerICARUS_H
#define PixelMapProducerICARUS_H

#include <iostream>
#include <ostream>
#include <list>
#include <algorithm>
#include <numeric>
#include <cstdlib>

#include "larrecodnn/CVN/module_helpers/PixelMapProducer.h"

namespace cvn
{
  template <class T, class U> class PixelMapProducerICARUS : public PixelMapProducer<T,U>
  {
    public:
	PixelMapProducerICARUS(unsigned int nWire, unsigned int nTdc, double tRes, double threshold = 0.):PixelMapProducer<T,U>::PixelMapProducer(nWire, nTdc, tRes, threshold){std::cout << "============ Calling the function PixelMapProducerICARUS::PixelMapProducerICARUS() ==============\n";}
        PixelMapProducerICARUS():PixelMapProducer<T,U>::PixelMapProducer(){std::cout << "============ Calling the function PixelMapProducerICARUS::PixelMapProducerICARUS() ==============\n";}
	PixelMapProducerICARUS(const fhicl::ParameterSet& pset):PixelMapProducer<T,U>::PixelMapProducer(pset),fverbose(pset.get<bool>("verbose")),fChangeWireNo(pset.get<bool>("ChangeWireNo")),fReadoutSize(pset.get<double>("ReadoutSize")),fShiftT(pset.get<float>("ShiftT")),fInductionWires(pset.get<int>("InductionWires")),fFlipInductionView(pset.get<bool>("FlipInductionView")),fUseT(pset.get<bool>("UseT")) {std::cout << "============ Calling the function PixelMapProducerICARUS::PixelMapProducerICARUS() ==============\n";}
	Boundary DefineBoundary(detinfo::DetectorPropertiesData const& detProp,const std::vector< const T* >& cluster) override;
	void ConvertLocaltoGlobal(geo::WireID wireid, unsigned int &globalWire, unsigned int &globalPlane) const override; 
	void ConvertLocaltoGlobalTDC(geo::WireID wireid, double localTDC, unsigned int &globalWire, unsigned int &globalPlane, double &globalTDC) const override;
	PixelMap CreateMapGivenBoundary(detinfo::DetectorPropertiesData const& detProp,const std::vector< const T* >& cluster,const Boundary& bound) override;
	PixelMap CreateMap(detinfo::DetectorPropertiesData const& detProp,const std::vector< art::Ptr< T > >& cluster);
	PixelMap CreateMap(detinfo::DetectorPropertiesData const& detProp,const std::vector< const T* >& cluster);
	void Set_fT0_value(float value);
    protected:
        bool fverbose;
        bool fChangeWireNo; 
        double fReadoutSize; // in time ticks
	float fShiftT; // size of the back/front porch
	int fInductionWires; // number of wires in the first induction plane
	bool fFlipInductionView; // should we flip the induction view
	bool fUseT; // 
	float fT0; // T0 coming from the PFP particles (in time ticks)
  };
  
  typedef PixelMapProducerICARUS<recob::Hit, HitHelper> PixelMapHitProducerICARUS;
  typedef PixelMapProducerICARUS<recob::Wire, WireHelper> PixelMapWireProducerICARUS;
  typedef PixelMapProducerICARUS<sim::SimChannel, SimChannelHelper> PixelMapSimProducerICARUS;
}

#endif
