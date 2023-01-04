#include "icaruscode/TPC/CVN/module_helpers/ICVNMapperICARUS.h"
#include "lardataobj/RecoBase/Slice.h"

#include "canvas/Persistency/Common/FindManyP.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/T0.h"

#include <algorithm> // std::min_element
#include <iterator>  // std::begin, std::end

namespace cvn
{
  ////////////////////////////////////////////////////////////////////////////////
	
  template <class T, class U> ICVNMapperICARUS<T, U>::~ICVNMapperICARUS()
  {

  }
  
  //////////////////////////////////////////////////////////////////////////////// 
  
  template <class T, class U> void ICVNMapperICARUS<T, U>::produce(art::Event& evt)
  {
   if(fverbose) std::cout << "============ Calling the function ICVNMapperICARUS::produce() ==============\n";
   
   if(fUseSlice){
      if(fverbose) std::cout << "============ Calling the function ICVNMapperICARUS::produce() is using slices ==============\n";
      art::Handle< std::vector<recob::Slice> > SliceListHandle;
      std::vector< art::Ptr<recob::Slice> > SliceList;
      if( evt.getByLabel(fSliceLabel,SliceListHandle))
          art::fill_ptr_vector(SliceList,SliceListHandle);
      
      art::Handle< std::vector<recob::PFParticle> > PFPListHandle;
      std::vector< art::Ptr<recob::PFParticle> > PFPList;
      if( evt.getByLabel(fPFParticleModuleLabel,PFPListHandle))
          art::fill_ptr_vector(PFPList,PFPListHandle);
      
      art::FindManyP<U> findManyHits(SliceListHandle, evt, fSliceLabel);
      
      art::FindManyP<recob::PFParticle> findManyPFPs(SliceListHandle, evt, fPFParticleModuleLabel);
      art::FindManyP<anab::T0> findManyT0s(PFPListHandle, evt, fT0Label);
      
      std::unique_ptr< std::vector<PixelMap> > pmCol(new std::vector<PixelMap>);
       
      auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt);
      
      for(auto const& slice : SliceList){
	  std::vector<float> pfp_T0_vec;
	  if(findManyPFPs.isValid()){
	     std::vector<art::Ptr<recob::PFParticle>> slicePFPs = findManyPFPs.at(slice.key());
	     if(slicePFPs.size()){
	        for(auto const &pfp : slicePFPs){
		    if(findManyT0s.isValid()){
		       std::vector<art::Ptr<anab::T0>> T0_vec = findManyT0s.at(pfp.key());
		       if(T0_vec.size()){
		          for(auto const& T0 : T0_vec){
			      pfp_T0_vec.push_back(T0->Time());
			  }
		       }
	            }
	        }
	     }
	  }    
	  
	  float min_T0 = 0.;
	  if(pfp_T0_vec.size()){
	     min_T0 = *min_element(pfp_T0_vec.begin(), pfp_T0_vec.end());
	  } 		      
	      
          if(findManyHits.isValid()){
	     std::vector<art::Ptr<U>> slicehits = findManyHits.at(slice.key());
	     this->fProducer.Set_fT0_value(min_T0);
	     //std::cout << "============== T0 provided to make the pixel map : " << min_T0 << "  ================\n";
	     PixelMap pm = this->fProducer.CreateMap(detProp, slicehits);
             auto nhits = this->fProducer.NROI();
             pm.SetTotHits(nhits);
	     pm.fSliceID = slice->ID();
	     //pm.fT0 = min_T0;
             
	     if(nhits > this->fMinClusterHits && pmCol->size()<fMapVecSize){
		std::cout << "********* " << evt.run() << "  " << evt.subRun() << "  " << evt.id().event() << "  " << slice->ID() << "  **************\n";
	        pmCol->push_back(pm);
	     }
	  
          }
      }
      
      evt.put(std::move(pmCol), this->fClusterPMLabel);
      
   }
   
   else{
       if(fverbose) std::cout << "============ Calling the function ICVNMapperICARUS::produce() is using full event ==============\n";
       ICVNMapper<T,U>::produce(evt);
   }
   
   //std::cout << "=============================== REACHED the END of the produce() function ==========================\n";
   
  }
  
  ////////////////////////////////////////////////////////////////////////////////
  
  template <class T, class U> void ICVNMapperICARUS<T, U>::beginJob()
  {
    if(fverbose) std::cout << "============ Calling the function ICVNMapperICARUS::beginJob() ==============\n";
  }
  
  ////////////////////////////////////////////////////////////////////////////////
  
  template <class T, class U> void ICVNMapperICARUS<T, U>::endJob()
  {
    if(fverbose) std::cout << "============ Calling the function ICVNMapperICARUS::endJob() ==============\n";
  }
  
  ////////////////////////////////////////////////////////////////////////////////
}
