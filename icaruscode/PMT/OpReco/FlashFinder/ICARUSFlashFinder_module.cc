/**
 * @file   icaruscode/PMT/OpReco/FlashFinder/ICARUSFlashFinder_module.cc
 * @brief  ICARUS customization of flash reconstruction in LArSoft
 * @author modified from Kazuhiro Terao's original module
 */


#include "icaruscode/PMT/OpReco/Algorithms/OpHitTimeSelector.h"

#include "art/Framework/Core/SharedProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Utilities/Exception.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/DelegatedParameter.h"

#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/OpFlash.h"

#include <iostream> // std::cerr
#include <memory>
#include <string>
#include <utility> // std::move()
#include <vector>
#include "FlashFinderManager.h"
#include "FlashFinderFMWKInterface.h" // ::pmtana::OpDetCenterFromOpChannel()


class ICARUSFlashFinder : public art::SharedProducer {
public:
  
  struct Config {
    
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    fhicl::Atom<art::InputTag> OpHitProducer {
      Name{ "OpHitProducer" },
      Comment{ "input tag for the optical hits to build the flashes from" }
      };
    
    fhicl::Atom<std::string> FlashFinderAlgo {
      Name{ "FlashFinderAlgo" },
      Comment{ "name of the flash finding algorithm" }
      };
    
    fhicl::Atom<std::string> TimeType {
      Name{ "TimeType" },
      Comment{ "Type of hit time to be used for flash building" },
      "Start"
      };
    
    fhicl::DelegatedParameter AlgoConfig {
      Name{ "AlgoConfig" },
      Comment{ "configuration of the flash finding algorithm" }
      };
    
  }; // Config
  
  using Parameters = art::SharedProducer::Table<Config>;
  
  explicit ICARUSFlashFinder(Parameters const & p, art::ProcessingFrame const&);

  void produce(art::Event & e, art::ProcessingFrame const&) override;


private:

  // Declare member data here.
  ::pmtana::FlashFinderManager _mgr;
  art::InputTag _hit_producer;
  
  /// Extracts a configured time from `recob::OpHit`.
  recob::OpHitTimeSelector const fHitTime;

  void GetFlashLocation(std::vector<double>, double&, double&, double&, double&);

};


ICARUSFlashFinder::ICARUSFlashFinder(Parameters const & p, art::ProcessingFrame const&)
  : art::SharedProducer{p}
  , fHitTime{ recob::opHitTimeType(p().TimeType()) }
// Initialize member data here.
{
  async<art::InEvent>();
  
  _hit_producer   = p().OpHitProducer();
  
  auto const flash_algo  = p().FlashFinderAlgo();
  auto const flash_pset = p().AlgoConfig.get<pmtana::Config_t>();
  auto algo_ptr = ::pmtana::FlashAlgoFactory::get().create(flash_algo,flash_algo);
  algo_ptr->Configure(flash_pset);
  _mgr.SetFlashAlgo(algo_ptr);

  produces< std::vector<recob::OpFlash>   >();
  produces< art::Assns <recob::OpHit, recob::OpFlash> >();
}

void ICARUSFlashFinder::produce(art::Event & e, art::ProcessingFrame const&)
{

  // produce OpFlash data-product to be filled within module
  auto opflashes = std::make_unique<std::vector<recob::OpFlash>>();
  auto flash2hit_assn_v = std::make_unique<art::Assns<recob::OpHit, recob::OpFlash>>();
  // load OpHits previously created
  art::Handle<std::vector<recob::OpHit> > ophit_h;
  e.getByLabel(_hit_producer,ophit_h);

  // make sure hits look good
  if(!ophit_h.isValid()) {
    std::cerr<<"\033[93m[ERROR]\033[00m ... could not locate OpHit!"<<std::endl;
    throw art::Exception{ art::errors::ProductNotFound }
      << "Couldn't read OpHit data product from '" << _hit_producer.encode() << "'";
  }

  ::pmtana::LiteOpHitArray_t ophits;
  double trigger_time=1.1e20;
  for(auto const& oph : *ophit_h) {
    ::pmtana::LiteOpHit_t loph;
    if(trigger_time > 1.e20) trigger_time = oph.PeakTimeAbs() - oph.PeakTime();
    loph.peak_time = fHitTime(oph);

    loph.pe = oph.PE();
    loph.channel = oph.OpChannel();
    ophits.emplace_back(std::move(loph));
  }
  
  auto const flash_v = _mgr.RecoFlash(ophits);

  art::PtrMaker<recob::OpFlash> makeFlashPtr{ e };
  for(const auto& lflash :  flash_v) {

    double Ycenter, Zcenter, Ywidth, Zwidth;
    GetFlashLocation(lflash.channel_pe, Ycenter, Zcenter, Ywidth, Zwidth);

    recob::OpFlash flash(lflash.time, lflash.time_err,
			 trigger_time + lflash.time,
			 (trigger_time + lflash.time) / 1600.,
			 lflash.channel_pe,
                         0, 0, 1, // this are just default values
                         Ycenter, Ywidth, Zcenter, Zwidth);
    opflashes->emplace_back(std::move(flash));

    art::Ptr<recob::OpFlash> const flashPtr = makeFlashPtr(opflashes->size() - 1);
    for(auto const& hitidx : lflash.asshit_idx) {
      const art::Ptr<recob::OpHit> hit_ptr(ophit_h, hitidx);
      flash2hit_assn_v->addSingle(hit_ptr, flashPtr);
    }
  }
  
  e.put(std::move(opflashes));
  e.put(std::move(flash2hit_assn_v));
}

void ICARUSFlashFinder::GetFlashLocation(std::vector<double> pePerOpChannel, 
                                     double& Ycenter, 
                                     double& Zcenter, 
                                     double& Ywidth, 
                                     double& Zwidth)
{

  // Reset variables
  Ycenter = Zcenter = 0.;
  Ywidth  = Zwidth  = -999.;
  double totalPE = 0.;
  double sumy = 0., sumz = 0., sumy2 = 0., sumz2 = 0.;

  for (unsigned int opch = 0; opch < pePerOpChannel.size(); opch++) {
  /*
    if (opch > 31 && opch < 200){
      //  std::cout << "Ignoring channel " << opch << " as it's not a real channel" << std::endl;                                                                                             
      continue;
    }
  */
    // Get physical detector location for this opChannel
    double PMTxyz[3];
    ::pmtana::OpDetCenterFromOpChannel(opch, PMTxyz);

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

DEFINE_ART_MODULE(ICARUSFlashFinder)
