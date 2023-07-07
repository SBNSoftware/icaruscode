#include "icaruscode/TPC/CVN/interfaces/ICVNMapperICARUS.h"
#include "lardataobj/RecoBase/Slice.h"

#include "canvas/Persistency/Common/FindManyP.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/T0.h"

#include <algorithm> // std::min_element
#include <iterator>  // std::begin, std::end

namespace lcvn
{
  ////////////////////////////////////////////////////////////////////////////////
	
  template <class T, class U> ICVNMapperICARUS<T, U>::~ICVNMapperICARUS()
  {

  }
  
  //////////////////////////////////////////////////////////////////////////////// 
  
  template <class T, class U> void ICVNMapperICARUS<T, U>::produce(art::Event& evt)
  {
    if(fverbose) std::cout << "============ Calling the function ICVNMapperICARUS::produce() ==============\n";
   
    // collect the input TPC reco tags
    std::vector<std::string> pandora_tag_suffixes = fPandoraTagSuffixes;
    if (pandora_tag_suffixes.size() == 0) {
      pandora_tag_suffixes.push_back("");
    }
    if(fUseSlice){
      if(fverbose) std::cout << "============ Calling the function ICVNMapperICARUS::produce() is using slices ==============\n";

      // collect the TPC slices
      std::vector<art::Ptr<recob::Slice>> slices;
      std::vector<std::string> slice_tag_suffixes;
      std::vector<unsigned> slice_tag_indices;
      for (unsigned i_tag = 0; i_tag < pandora_tag_suffixes.size(); i_tag++) {
        const std::string &pandora_tag_suffix = pandora_tag_suffixes[i_tag];
        art::Handle<std::vector<recob::Slice>> slices_handle;
        evt.getByLabel(fSliceLabel + pandora_tag_suffix, slices_handle);
        if(slices_handle.isValid()){
          art::fill_ptr_vector(slices, slices_handle);
          for (unsigned i = 0; i < slices_handle->size(); i++) {
            slice_tag_suffixes.push_back(pandora_tag_suffix); 
            slice_tag_indices.push_back(i_tag);
          }
        }
      }

      for (unsigned i = 0; i < slice_tag_suffixes.size(); i++) {
        std::cout << slice_tag_suffixes[i] << std::endl; 
        std::cout << slice_tag_indices[i] << std::endl; 
      }

      std::unique_ptr< std::vector<PixelMap> > pmCol(new std::vector<PixelMap>);
       
      auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt);
      
      for(unsigned int sliceID=0; sliceID<slices.size(); sliceID++){
        art::Ptr<recob::Slice> slice = slices[sliceID];
        const std::string &slice_tag_suff = slice_tag_suffixes[sliceID];

        // get handle of slices
        art::Handle<std::vector<recob::Slice>> slices_handle;
        evt.getByLabel(fSliceLabel + slice_tag_suff, slices_handle);
        // get handle of pfps
        art::Handle<std::vector<recob::PFParticle>> pfps_handle;
        evt.getByLabel(fPFParticleModuleLabel + slice_tag_suff, pfps_handle);

        std::vector<float> pfp_T0_vec;
        art::FindManyP<recob::PFParticle> findManyPFPs(slices_handle, evt, fPFParticleModuleLabel + slice_tag_suff);
        art::FindManyP<recob::Hit> findManyHits(slices_handle, evt, fSliceLabel + slice_tag_suff);
        art::FindManyP<anab::T0> findManyT0s(pfps_handle, evt, fT0Label + slice_tag_suff);

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
          PixelMap pm = this->fProducer.CreateMap(detProp, slicehits);
          auto nhits = this->fProducer.NROI();
          pm.SetTotHits(nhits);
          pm.fSliceID = slice->ID();

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
  
} // namespace lcvn
