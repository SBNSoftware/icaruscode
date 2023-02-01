#include <ostream>

#include "TH2D.h"
#include "TVector2.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "larrecodnn/CVN/func/AssignLabels.h"
#include "larrecodnn/CVN/interfaces/PixelMapProducer.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

#include "icaruscode/TPC/CVN/interfaces/PixelMapProducerICARUS.h"

namespace lcvn
{
 template <class T, class U> void PixelMapProducerICARUS<T, U>::ConvertLocaltoGlobal(geo::WireID wireid, 
                                      unsigned int &globalWire, unsigned int &globalPlane) const
  {
    if(fverbose) std::cout << "============ Calling the function PixelMapProducerICARUS::ConvertLocaltoGlobal() ==============\n";
    globalWire = wireid.Wire;
    globalPlane = wireid.Plane;
  }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////		  
  template <class T, class U> Boundary PixelMapProducerICARUS<T, U>::DefineBoundary(detinfo::DetectorPropertiesData const& detProp,
                                            const std::vector< const T*>& cluster)
  {
    if(fverbose) std::cout << "============ Calling the function PixelMapProducerICARUS::DefineBoundary() ==============\n";
    std::vector<double> tmin_0;
    std::vector<double> tmin_1;
    std::vector<double> tmin_2;

    std::vector<int> wire_0, bwire_0;
    std::vector<int> wire_1, bwire_1;
    std::vector<int> wire_2, bwire_2;
    
    std::vector<double> tsum = {0., 0., 0.};
    std::vector<double> tsize = {0., 0., 0.};
    
    for(size_t iHit = 0; iHit < cluster.size(); ++iHit)
    {
      U wraphit(*(cluster[iHit]), this->fThreshold);
      Waveform wf = wraphit.GetWaveform();
      geo::WireID wireid = wraphit.GetID();

      unsigned int tempWire  = wireid.Wire;
      unsigned int tempPlane = wireid.Plane;

      if(!this->fMultipleDrifts)
         ConvertLocaltoGlobal(wireid, tempWire, tempPlane);

      for(auto &pulse: wf){
        double min_tick = (double)INT_MAX;
        for(auto &i: pulse){
          double temptdc = i.first;
          
          if(this->fMultipleDrifts)
            ConvertLocaltoGlobalTDC(wireid, i.first, tempWire, tempPlane, temptdc);
          
          if(temptdc < min_tick) min_tick = temptdc;

          tsum[tempPlane] += temptdc;
          tsize[tempPlane] += 1.;
        }
        
        if(!(pulse.empty())){ 
          if(tempPlane == 0){
            tmin_0.push_back(min_tick);
            wire_0.push_back(tempWire);
          }
          if(tempPlane == 1){ 
            tmin_1.push_back(min_tick);
            //wire_0.push_back(tempWire);
	    wire_1.push_back(tempWire);
          }
          if(tempPlane == 2){ 
            tmin_2.push_back(min_tick);
            //wire_0.push_back(tempWire);
	    wire_2.push_back(tempWire);
          }
        }
      } // end loop over pulses on single wire
    } // end loop over struck wires
    
    double tmean_0 = tsum[0]/tsize[0];
    double tmean_1 = tsum[1]/tsize[1];
    double tmean_2 = tsum[2]/tsize[2]; 
    
    for(int i = 0; i < (int)wire_0.size(); i++){
      if(std::abs(tmin_0[i]-tmean_0) < (double)this->fTRes)
        bwire_0.push_back(wire_0[i]);
    }
    for(int i = 0; i < (int)wire_1.size(); i++){
      if(std::abs(tmin_1[i]-tmean_1) < (double)this->fTRes)
        bwire_1.push_back(wire_1[i]);
    }
    for(int i = 0; i < (int)wire_2.size(); i++){
      //std::cout << wire_2[i] << "  " << tmin_2[i] << "\n";
      if(std::abs(tmin_2[i]-tmean_2) < (double)this->fTRes){
         bwire_2.push_back(wire_2[i]);
	 //std::cout << wire_2[i] << "  " << tmin_2[i] << "\n";
      }
    }
      

    std::cout << "Boundary wire vector sizes: " << bwire_0.size() << ", " << bwire_1.size() << ", " << bwire_2.size() << std::endl;

    int minwire_0 = 0;
    int minwire_1 = 0;
    int minwire_2 = 0;
    auto minwireelement_0 = std::min_element(bwire_0.begin(), bwire_0.end());
    auto minwireelement_1 = std::min_element(bwire_1.begin(), bwire_1.end());
    auto minwireelement_2 = std::min_element(bwire_2.begin(), bwire_2.end());
    
    if(fChangeWireNo){
       if(bwire_0.size() > 0) { minwire_0 = *minwireelement_0-1; std::cout<<"minwire 0: "<<(*minwireelement_0-1) <<std::endl;}
       if(bwire_1.size() > 0) { minwire_1 = *minwireelement_1-1; std::cout<<"minwire 1: "<<(*minwireelement_1-1) <<std::endl;}
       if(bwire_2.size() > 0) { minwire_2 = *minwireelement_2-1; std::cout<<"minwire 2: "<<(*minwireelement_2-1) <<std::endl;}
    }
    
    else{
       if(bwire_0.size() > 0) { minwire_0 = *minwireelement_0; std::cout<<"minwire 0: "<<(*minwireelement_0) <<std::endl;}
       if(bwire_1.size() > 0) { minwire_1 = *minwireelement_1; std::cout<<"minwire 1: "<<(*minwireelement_1) <<std::endl;}
       if(bwire_2.size() > 0) { minwire_2 = *minwireelement_2; std::cout<<"minwire 2: "<<(*minwireelement_2) <<std::endl;}
    }
    
    this->fTotHits = bwire_0.size() + bwire_1.size() + bwire_2.size();

    Boundary bound(this->fNWire,this->fTRes,minwire_0,minwire_1,minwire_2,tmean_0,tmean_1,tmean_2);

    return bound;
  }
		  
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  template <class T, class U> void PixelMapProducerICARUS<T, U>::ConvertLocaltoGlobalTDC(geo::WireID wireid, double localTDC, 
                                       unsigned int &globalWire, unsigned int &globalPlane, 
                                       double &globalTDC) const
  {
    if(fverbose) std::cout << "============ Calling the function PixelMapProducerICARUS::ConvertLocaltoGlobalTDC() ==============\n";
    /*globalWire = wireid.Wire;
    globalPlane = wireid.Plane;
    if(wireid.TPC == 0) globalTDC = localTDC + fShiftT;
    else if(wireid.TPC == 1) globalTDC = fReadoutSize + (fReadoutSize - localTDC) - fShiftT;*/
	
    /////////////////////// calculating T0 correction ///////////////////	    
	        
    float new_shift = 0.;
    
    if(fUseT){
       /*std::cout << "***************** WE ARE USING PFP PARTICLE T0 *********************\n";
       std::cout << "T0 Value given to pixel map producer : " << fT0 << "\n";*/
       if(fT0 < 0.){
          new_shift = float(abs(fT0))/500 + fShiftT;
       }
       else if(fT0 > 0.){
            if(float(abs(fT0)/500 < fShiftT)) new_shift = fShiftT - float(abs(fT0))/500; 
       } 
       
       else new_shift = fShiftT;
    }
    
    else{
      //std::cout << "***************** WE ARE NOT USING PFP PARTICLE T0 *********************\n";
      new_shift = fShiftT;
    }
    
    ///////////////////////// end of calculating T0 correction /////////    
    
    if(wireid.TPC == 0){
       globalWire = wireid.Wire;
       globalPlane = wireid.Plane;
       //globalTDC = localTDC + fShiftT;
       globalTDC = localTDC + new_shift;
       if(fFlipInductionView){
          if(globalPlane == 1)globalWire = fInductionWires-globalWire;
       }
    }
    
    else{
      if(wireid.Plane == 2){
         globalWire = wireid.Wire;
	 globalPlane = wireid.Plane;
	 //globalTDC = fReadoutSize + (fReadoutSize - localTDC) - fShiftT;
	 globalTDC = fReadoutSize + (fReadoutSize - localTDC) - new_shift;
      }
      
      else if(wireid.Plane == 1){
         globalWire = wireid.Wire;
	 globalPlane = 0;
	 //globalTDC = fReadoutSize + (fReadoutSize - localTDC) - fShiftT;
	 globalTDC = fReadoutSize + (fReadoutSize - localTDC) - new_shift;
      }
      
      else{
         globalWire = wireid.Wire;
	 globalPlane = 1;
	 //globalTDC = fReadoutSize + (fReadoutSize - localTDC) - fShiftT;
	 globalTDC = fReadoutSize + (fReadoutSize - localTDC) - new_shift;
	 if(fFlipInductionView){
	    globalWire = fInductionWires-globalWire;
	 }
      }
    }
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  template <class T, class U> PixelMap PixelMapProducerICARUS<T, U>::CreateMapGivenBoundary(detinfo::DetectorPropertiesData const& detProp,const std::vector< const T* >& cluster,
		  const Boundary& bound)
  {
    if(fverbose) std::cout << "============ Calling the function PixelMapProducerICARUS::CreateMapGivenBoundary() ==============\n";
    PixelMap pm(this->fNWire, this->fNTdc, bound);

    for(size_t iHit = 0; iHit < cluster.size(); ++iHit)
    {

      U wraphit(*(cluster[iHit]), this->fThreshold);
      Waveform wf = wraphit.GetWaveform();
      geo::WireID wireid = wraphit.GetID();

      unsigned int tempWire  = wireid.Wire;
      unsigned int tempPlane = wireid.Plane;

      if(!this->fMultipleDrifts)
         ConvertLocaltoGlobal(wireid, tempWire, tempPlane);

      for(auto &pulse: wf){
        // Leigh: Simple modification to unwrap the collection view wire plane
        for(auto &i: pulse){
          const double pe = i.second;
          double temptdc = i.first;
          if(this->fMultipleDrifts)
            ConvertLocaltoGlobalTDC(wireid, i.first, tempWire, tempPlane, temptdc); 
       
          const unsigned int wire = tempWire;
          const unsigned int wirePlane = tempPlane;
          const double tdc = temptdc;
        
          pm.Add(wire, tdc, wirePlane, pe);
        }
      }

    }
    pm.SetTotHits(this->fTotHits);
    return pm;
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  template <class T, class U> PixelMap PixelMapProducerICARUS<T, U>::CreateMap(detinfo::DetectorPropertiesData const& detProp,const std::vector< art::Ptr< T > >& cluster)
  {
    if(fverbose) std::cout << "============ Calling the function PixelMapProducerICARUS::CreateMap() ==============\n";
    std::vector<const T*> newCluster;
    for(const art::Ptr<T> hit : cluster){
      newCluster.push_back(hit.get());
    }
    return CreateMap(detProp, newCluster);
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  template <class T, class U> PixelMap PixelMapProducerICARUS<T, U>::CreateMap(detinfo::DetectorPropertiesData const& detProp,const std::vector< const T* >& cluster)
  {
    if(fverbose) std::cout << "============ Calling the function PixelMapProducerICARUS::CreateMap() ==============\n";
    Boundary bound = DefineBoundary(detProp, cluster);
    return CreateMapGivenBoundary(detProp, cluster, bound);
  }
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  template <class T, class U> void PixelMapProducerICARUS<T, U>::Set_fT0_value(float value)
  {
    fT0 = value;
  }
  
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  template class PixelMapProducerICARUS<recob::Hit, HitHelper>;
  template class PixelMapProducerICARUS<recob::Wire, WireHelper>;
  template class PixelMapProducerICARUS<sim::SimChannel, SimChannelHelper>;

} // namespace lcvn

