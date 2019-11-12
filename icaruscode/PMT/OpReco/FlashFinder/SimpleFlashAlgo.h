#ifndef SIMPLEFLASHALGO_H
#define SIMPLEFLASHALGO_H

#include "FlashAlgoBase.h"
#include "FlashAlgoFactory.h"
#include <map>
namespace pmtana
{

  class SimpleFlashAlgo : public FlashAlgoBase {

  public:

    SimpleFlashAlgo(const std::string name);

    void Configure(const Config_t &p);
    
    virtual ~SimpleFlashAlgo();

    LiteOpFlashArray_t RecoFlash(const LiteOpHitArray_t ophits);

    bool Veto(double t) const;

    const std::vector<double>& PESumArray() const { return _pesum_v; }

    const double TimeRes() const { return _time_res; }

  private:

    double TotalCharge(const std::vector<double>& PEs);

    // minimum PE to account for a hit
    double _min_pe_hit;

    // minimum PE to make a flash
    double _min_pe_flash;

    // minimum PE to make a flash candidate
    double _min_pe_coinc;

    // minimum PE to use an OpHit
    double _min_mult_coinc;

    // integral period
    double _integral_time;

    // veto time
    double _veto_time;

    // time resolution of pe sum
    double _time_res;

    // time pre-sample
    double _pre_sample;

    // pw aum array
    std::vector<double> _pesum_v;

    // calibration: PEs to be subtracted from each opdet
    std::vector<double> _pe_baseline_v;

    // veto window start
    std::map<double,double> _flash_veto_range_m;

    // debug mode flag
    bool _debug;

    // list of opchannel to use
    std::vector<int> _opch_to_index_v;
    std::vector<int> _index_to_opch_v;
    
  };

  /**
     \class pmtana::SimpleFlashAlgoFactory
     \brief A concrete factory class for pmtana::SimpleFlashAlgo
  */
  class SimpleFlashAlgoFactory : public FlashAlgoFactoryBase {
  public:
    /// ctor
    SimpleFlashAlgoFactory() { FlashAlgoFactory::get().add_factory("SimpleFlashAlgo",this); }
    /// dtor
    ~SimpleFlashAlgoFactory() {}
    /// creation method
    FlashAlgoBase* create(const std::string instance_name) { return new SimpleFlashAlgo(instance_name); }
  };

}
#endif

/** @} */ // end of doxygen group
