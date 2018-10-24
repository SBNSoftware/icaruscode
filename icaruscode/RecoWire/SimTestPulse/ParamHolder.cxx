#include "ParamHolder.h"

namespace alternative {

  ParamHolder* ParamHolder::_me = nullptr;

  void ParamHolder::Register(alternative::TruthHit&& hit)
  { 
    _hit_v.emplace_back(std::move(hit)); 
    _hit_v.back().signal_id = _hit_v.size() - 1;
  }
  
  const std::vector<alternative::TruthHit>& ParamHolder::TruthHitArray() const
  { return _hit_v; }

  void ParamHolder::Clear()
  { _hit_v.clear(); }
}
