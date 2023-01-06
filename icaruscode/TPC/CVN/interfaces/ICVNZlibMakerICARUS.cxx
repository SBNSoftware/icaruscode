#include "icaruscode/TPC/CVN/interfaces/ICVNZlibMakerICARUS.h"
#include "lardataobj/RecoBase/Slice.h"

#include "canvas/Persistency/Common/FindManyP.h"

namespace cvn
{
////////////////////////////////////////////////////////////////////////	
	
ICVNZlibMakerICARUS::ICVNZlibMakerICARUS(fhicl::ParameterSet const& pset)
    : ICVNZlibMaker(pset)
  {
    std::cout << "============ Calling the function ICVNZlibMakerICARUS::ICVNZlibMakerICARUS() ==============\n";
    reconfigure(pset); 
  }
  
///////////////////////////////////////////////////////////////////////
  
ICVNZlibMakerICARUS::~ICVNZlibMakerICARUS()
{
}

///////////////////////////////////////////////////////////////////////

void ICVNZlibMakerICARUS::reconfigure(const fhicl::ParameterSet& pset)
{
  std::cout << "============ Calling the function ICVNZlibMakerICARUS::reconfigure() ==============\n";
  ICVNZlibMaker::reconfigure(pset);
  fverbose = pset.get<bool>("verbose");
  fUseSlice = pset.get<bool>("UseSlice");
  fSliceLabel = pset.get<std::string>("SliceLabel");
}

//////////////////////////////////////////////////////////////////////

void ICVNZlibMakerICARUS::beginJob()
{
     std::cout << "============ Calling the function ICVNZlibMakerICARUS::beginJob() ==============\n";
     ICVNZlibMaker::beginJob();
}

/////////////////////////////////////////////////////////////////////

}

