#ifndef __SIMTESTPULSE_PARAMHOLDER_H__
#define __SIMTESTPULSE_PARAMHOLDER_H__

#include <array>
#include <vector>

namespace alternative {

  struct TruthHit {
    int signal_id;
    std::array<int,3> channel_list;
    unsigned int tdc;
    unsigned int tick;
    double num_electrons;
  };

  class ParamHolder {
  private:
    ParamHolder() {}
  public:
    ~ParamHolder() {}

    static ParamHolder& get()
    { if(!_me) _me = new ParamHolder;
      return *_me;
    }

    static void destroy()
    { if(_me) delete _me; }

    void Register(alternative::TruthHit&& hit);

    const std::vector<alternative::TruthHit>& TruthHitArray() const;

    void Clear();

  private:
    static ParamHolder* _me;
    std::vector<alternative::TruthHit> _hit_v;
  };
}
#endif
