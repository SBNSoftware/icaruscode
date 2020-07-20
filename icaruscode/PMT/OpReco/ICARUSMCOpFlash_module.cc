////////////////////////////////////////////////////////////////////////
// Class:       ICARUSMCOpFlash
// Plugin Type: producer (art v3_01_02)
// File:        ICARUSMCOpFlash_module.cc
//
// Generated at Sun Mar  3 17:53:36 2019 by Kazuhiro Terao using cetskelgen
// from cetlib version v3_05_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/OpFlash.h"

#include <memory>

class ICARUSMCOpFlash;


class ICARUSMCOpFlash : public art::EDProducer {
public:
  explicit ICARUSMCOpFlash(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ICARUSMCOpFlash(ICARUSMCOpFlash const&) = delete;
  ICARUSMCOpFlash(ICARUSMCOpFlash&&) = delete;
  ICARUSMCOpFlash& operator=(ICARUSMCOpFlash const&) = delete;
  ICARUSMCOpFlash& operator=(ICARUSMCOpFlash&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  void GetFlashLocation(std::vector<double> pePerOpChannel,
			double& Ycenter,
			double& Zcenter,
			double& Ywidth,
			double& Zwidth) const;

  // Declare member data here.
  double _merge_period;
  double _pe_threshold;
  bool   _store_empty_flash;
  std::string _hit_label;
  std::string _mct_label;
  std::vector<bool> _enabled_opch_v;
};


ICARUSMCOpFlash::ICARUSMCOpFlash(fhicl::ParameterSet const& p)
  : EDProducer{p}  // ,
  // More initializers here.
{
  _pe_threshold = p.get<double>("PEThresholdHit",0.5);
  _store_empty_flash = p.get<bool>("StoreEmptyFlash",false);
  _mct_label = p.get<std::string>("MCTruthProducer");
  _hit_label = p.get<std::string>("OpHitProducer");
  _merge_period = p.get<double>("MergePeriod");
  _enabled_opch_v.clear();
  std::vector<size_t> enabled_ch_v;
  enabled_ch_v = p.get<std::vector<size_t> >("OpChannel",enabled_ch_v);
  if(enabled_ch_v.empty()) {
    auto opch_range = p.get<std::vector<size_t> >("OpChannelRange");
    if(opch_range.size()!=2) {
      std::cerr << "OpChannelRange must be a vector of length two!" << std::endl;
      throw std::exception();
    }else if(opch_range[0] > opch_range[1]) {
      std::cerr << "OpChannelRange 0th element (" << opch_range[0]
		<< ") must be larger than the 1st element (" << opch_range[1]
		<< ")" << std::endl;
      throw std::exception();
    }
    for(size_t opch=opch_range[0]; opch<=opch_range[1]; ++opch)
      enabled_ch_v.push_back(opch);
  }
  for(auto const& opch : enabled_ch_v) {
    if(_enabled_opch_v.size() <= opch) _enabled_opch_v.resize(opch+1,false);
    _enabled_opch_v[opch] = true;
  }

  produces<std::vector<recob::OpFlash> >();
}

void ICARUSMCOpFlash::produce(art::Event& e)
{
  auto const ts = lar::providerFrom<detinfo::DetectorClocksService>();
  auto const geop = lar::providerFrom<geo::Geometry>();
  auto opf_v = std::unique_ptr<std::vector<recob::OpFlash> >(new std::vector<recob::OpFlash>());

  art::Handle< std::vector< simb::MCTruth > > mct_h;
  e.getByLabel(_mct_label,mct_h);
  if(!mct_h.isValid()) {
    std::cerr << "std::vector<simb::MCTruth> could not be found with a label: " << _mct_label << std::endl;
    throw std::exception();
  }

  art::Handle< std::vector< recob::OpHit > > oph_h;
  e.getByLabel(_hit_label,oph_h);
  if(!oph_h.isValid()) {
    std::cerr << "std::vector<recob::OpHit> could not be found with a label: " << _hit_label << std::endl;
    throw std::exception();
  }

  std::set<double> flash_time_s;
  for(auto const& mct : *mct_h) {
    for(int i=0; i<mct.NParticles(); ++i) {
      auto const& part = mct.GetParticle(i);
      if( (int)(part.StatusCode()) != 1 ) continue;
      double flash_time = ts->G4ToElecTime(part.T()) - ts->TriggerTime();
      flash_time_s.insert(flash_time);
    }
  }
  std::vector<double> flash_time_v;
  flash_time_v.reserve(flash_time_s.size());
  double last_time = -1.1e20;
  for(auto const& time : flash_time_s) {
    if(time > (last_time + _merge_period)) flash_time_v.push_back(time);
    last_time = time;
  }

  for(auto const& flash_time : flash_time_v) {
    std::vector<double> pe_v(geop->NOpChannels(),0.);
    double pe_total=-1.;
    double pe_sum=0.;
    double pe_sum1=0.;
    double pe_sum2=0.;
    for(auto const& oph : *oph_h) {
      auto opch = oph.OpChannel();
      if(oph.PE()<_pe_threshold) continue;
      pe_sum += oph.PE();
      if((int)(_enabled_opch_v.size()) <= opch || !_enabled_opch_v[opch]) continue;
      pe_sum1 += oph.PE();
      if(oph.PeakTime() < flash_time || oph.PeakTime() > (flash_time + _merge_period)) {
	/*
	std::cout << "Ch " << opch << " Time: " << oph.PeakTime() 
		  << " PE " << oph.PE()
		  << " ... dt " << oph.PeakTime() - flash_time <<std::endl;
	*/
	continue;
      }
      pe_sum2 += oph.PE();
      pe_v[opch] += oph.PE();
      pe_total = (pe_total < 0. ? oph.PE() : pe_total + oph.PE());
    }
    //std::cout<<"PE sum " << pe_sum << " " << pe_sum1 << " " << pe_sum2 << std::endl;
    if(pe_total < 0. && !_store_empty_flash) continue;
    double y,z,ywidth,zwidth;
    GetFlashLocation(pe_v,y,z,ywidth,zwidth);
    recob::OpFlash f(flash_time, _merge_period, flash_time + ts->TriggerTime(), 0,
		     pe_v, false, 0, 0, 0, 0, 0, 0);
    opf_v->emplace_back(std::move(f));
  }
  e.put(std::move(opf_v));
}

void ICARUSMCOpFlash::GetFlashLocation(std::vector<double> pePerOpChannel,
				       double& Ycenter,
				       double& Zcenter,
				       double& Ywidth,
				       double& Zwidth) const
{

  // Reset variables
  auto const geop = lar::providerFrom<geo::Geometry>();
  Ycenter = Zcenter = -1.e20;
  Ywidth  = Zwidth  = -1.e20;
  double totalPE = 0.;
  double sumy = 0., sumz = 0., sumy2 = 0., sumz2 = 0.;

  for (unsigned int opch = 0; opch < pePerOpChannel.size(); opch++) {

    // Get physical detector location for this opChannel 
    double PMTxyz[3];
    geop->OpDetGeoFromOpChannel(opch).GetCenter(PMTxyz);

    // Add up the position, weighting with PEs 
    sumy    += pePerOpChannel[opch]*PMTxyz[1];
    sumy2   += pePerOpChannel[opch]*PMTxyz[1]*PMTxyz[1];
    sumz    += pePerOpChannel[opch]*PMTxyz[2];
    sumz2   += pePerOpChannel[opch]*PMTxyz[2]*PMTxyz[2];

    totalPE += pePerOpChannel[opch];
  }

  Ycenter = sumy/totalPE;
  Zcenter = sumz/totalPE;

  // This is just sqrt(<x^2> - <x>^2)
  if ( (sumy2*totalPE - sumy*sumy) > 0. )
    Ywidth = std::sqrt(sumy2*totalPE - sumy*sumy)/totalPE;

  if ( (sumz2*totalPE - sumz*sumz) > 0. )
    Zwidth = std::sqrt(sumz2*totalPE - sumz*sumz)/totalPE;
}

DEFINE_ART_MODULE(ICARUSMCOpFlash)
