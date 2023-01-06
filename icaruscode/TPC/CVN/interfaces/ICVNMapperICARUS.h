#ifndef ICVNMAPPERICARUS_H
#define ICVNMAPPERICARUS_H

#include <iostream>
#include <ostream>
#include <list>
#include <algorithm>
#include <numeric>

#include "larrecodnn/CVN/interfaces/ICVNMapper.h"

namespace cvn
{
  template <class T, class U> class ICVNMapperICARUS : public ICVNMapper<T,U>
  {
    public:
    explicit ICVNMapperICARUS(fhicl::ParameterSet const&pset):ICVNMapper<T,U>::ICVNMapper(pset),fverbose(pset.get<bool>("verbose")),fUseSlice(pset.get<bool>("UseSlice")),fSliceLabel(pset.get<std::string>("SliceLabel")),fPFParticleModuleLabel(pset.get<std::string>("PFParticleModuleLabel")),fT0Label(pset.get<std::string>("T0Label")),fMapVecSize(pset.get<int>("MapVecSize")) {std::cout << "============ Calling the function ICVNMapperICARUS::ICVNMapperICARUS() ==============\n";} 
    ~ICVNMapperICARUS();
    void produce(art::Event& evt);
    void beginJob();
    void endJob();
    protected:
     bool fverbose;
     bool fUseSlice; 
     std::string fSliceLabel;
     std::string fPFParticleModuleLabel;
     std::string fT0Label;
     unsigned int fMapVecSize;
  };
}

#endif
