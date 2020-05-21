/**
 * @file RecoWireAssn_module.cc
 *
 * @brief A module to associate recob::Wire and raw::RawDigits
 *
 * @author Wenqiang Gu
 * Contact: wenqiang.gu.1@gmail.com
 *
 */
#ifndef RecoWireAssn_H
#define RecoWireAssn_H

// Framework includes
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "cetlib/search_path.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/Assns.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Core/EDProducer.h" // include the proper bit of the framework
#include "art/Framework/Core/ModuleMacros.h" 

// LArSoft includes
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardata/ArtDataHelper/WireCreator.h"
#include "lardata/Utilities/AssociationUtil.h"

namespace recowire {

  class RecoWireAssn : public art::EDProducer {

  public:
    
    // create association of recob::Wire and raw::RawDigits
    // recob::Wire sometimes is from the other modules (eg larwirecell)    
    explicit RecoWireAssn(fhicl::ParameterSet const& pset); 
    virtual ~RecoWireAssn();
    
    void produce(art::Event& evt); 
    void beginJob(); 
    void endJob();                 
    void reconfigure(fhicl::ParameterSet const& p);
 
  private:
    
    std::string  fDigitLabel;  ///< digits label
    std::string  fWireLabel;  ///< wires label


  protected: 
    
  }; // class RecoWireAssn
}

namespace recowire{

  //-------------------------------------------------
  RecoWireAssn::RecoWireAssn(fhicl::ParameterSet const& pset) : EDProducer{pset}
  {
    this->reconfigure(pset);
    
    // produces< std::vector<recob::Wire> >();
    produces<art::Assns<raw::RawDigit, recob::Wire>>();
    
  }
  
  //-------------------------------------------------
  RecoWireAssn::~RecoWireAssn()
  {
  }

  //////////////////////////////////////////////////////
  void RecoWireAssn::reconfigure(fhicl::ParameterSet const& p)
  {
    fDigitLabel = p.get< std::string >("DigitLabel", "daq");
    fWireLabel  = p.get< std::string >("WireLabel", "gauss");
  }

  //-------------------------------------------------
  void RecoWireAssn::beginJob()
  {  
    MF_LOG_DEBUG("RecoWireAssn") << "RecoWireAssn_module: beginJob()";
  }

  //////////////////////////////////////////////////////
  void RecoWireAssn::endJob()
  {  
  }
  
  //////////////////////////////////////////////////////
  void RecoWireAssn::produce(art::Event& evt)
  {
    // make a collection of Wires
    // std::unique_ptr<std::vector<recob::Wire> > wirecol(new std::vector<recob::Wire>);

    // create an association set
    std::unique_ptr<art::Assns<raw::RawDigit,recob::Wire> > WireDigitAssn
      (new art::Assns<raw::RawDigit,recob::Wire>);

    // read in RawDigits
    art::Handle< std::vector<raw::RawDigit> > digits_handle;
    if (! evt.getByLabel(fDigitLabel, digits_handle)) {
        std::cout << "WARNING: no label " << fDigitLabel << " for RawDigit" << std::endl;
        return;
    }
    // std::vector< art::Ptr<raw::RawDigit> >  digits;
    // art::fill_ptr_vector(digits, digits_handle);
    // for (auto const& digit: digits) {
    //     int chanId = digit->Channel();
    //     std::cout << chanId << std::endl;
    // }

    // create a map between channel ID and the location in RawDigit handle
    std::map<raw::ChannelID_t, unsigned int> rdchan_loc;
    for(unsigned int rdIter = 0; rdIter < digits_handle->size(); ++rdIter){
      
      art::Ptr<raw::RawDigit> digit(digits_handle, rdIter);
      auto chanId = digit->Channel();
      rdchan_loc[chanId] = rdIter;
    }


    // read in recob::Wire
    art::Handle< std::vector<recob::Wire> > wires_handle;
    if (! evt.getByLabel(fWireLabel, wires_handle)) {
        std::cout << "WARNING: no label " << fWireLabel << " for recob::Wire" << std::endl;
        return;
    }
    std::vector< art::Ptr<recob::Wire> >  wires;
    art::fill_ptr_vector(wires, wires_handle);

    // create association
    for (auto const& wire: wires) {
        auto chanId = wire->Channel();
        // std::cout << chanId << std::endl;
        art::Ptr<raw::RawDigit> digit(digits_handle, rdchan_loc.at(chanId));
        util::CreateAssn(*this, evt, wire, digit, *WireDigitAssn);
    }
    
    if (!digits_handle->size() || !wires_handle->size())  return;
    mf::LogInfo("RecoWireAssn") << "RecoWireAssn:: digitVecHandle size is " << digits_handle->size()
                                << "wireVecHandle size is " << wires_handle->size();

    evt.put(std::move(WireDigitAssn));
    
    return;
  }
  
} // end namespace recowire


namespace recowire{

  DEFINE_ART_MODULE(RecoWireAssn)
  
} // end namespace recowire


#endif // RecoWireAssn_H

