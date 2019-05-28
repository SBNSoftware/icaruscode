#ifndef FLASHALGOBASE_CXX
#define FLASHALGOBASE_CXX

#include "FlashAlgoBase.h"

namespace pmtana{

  FlashAlgoBase::FlashAlgoBase(const std::string name)
  {
    _name = name;
    Reset();
  }

  FlashAlgoBase::~FlashAlgoBase() {}

  void FlashAlgoBase::Reset() {}
  
}
#endif

