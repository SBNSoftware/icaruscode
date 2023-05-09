#ifndef ICVNZLIBMAKERICARUS_H
#define ICVNZLIBMAKERICARUS_H

#include <iostream>
#include <ostream>
#include <list>
#include <algorithm>
#include <numeric>
#include <cstdlib>

//#include "zlib.h"
//#include "math.h"

#include "larrecodnn/CVN/interfaces/ICVNZlibMaker.h"

namespace lcvn
{
  class ICVNZlibMakerICARUS : public ICVNZlibMaker
  {
    public:
      explicit ICVNZlibMakerICARUS(fhicl::ParameterSet const& pset);
      ~ICVNZlibMakerICARUS();
      void beginJob();
      void analyze(const art::Event& evt){}
      void reconfigure(const fhicl::ParameterSet& pset);
    protected:
      bool fVerbose;
      bool fUseSlice;
      std::string fSliceLabel; 
  };
} // namespace lcvn

#endif
