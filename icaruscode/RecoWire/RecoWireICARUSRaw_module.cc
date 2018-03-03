////////////////////////////////////////////////////////////////////////
// $Id: RecoWireICARUSRaw.cxx,v 1.36 2010/09/15  bpage Exp $
//
// RecoWireICARUSRaw class
//
// brebel@fnal.gov
//
////////////////////////////////////////////////////////////////////////
#ifndef RecoWireICARUSRaw_H
#define RecoWireICARUSRaw_H

// ROOT includes
#include <TFile.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TF1.h>
#include <TComplex.h>

// Framework includes
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"
#include "cetlib/search_path.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/Assns.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Core/EDProducer.h" // include the proper bit of the framework
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"


// LArSoft includes
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardata/ArtDataHelper/WireCreator.h"
#include "lardata/Utilities/LArFFT.h"
#include "lardata/Utilities/AssociationUtil.h"

///creation of calibrated signals on wires
namespace recowireraw {

  class RecoWireICARUSRaw : public art::EDProducer {

  public:
    
    // create calibrated signals on wires. this class runs 
    // an fft to remove the electronics shaping.     
    explicit RecoWireICARUSRaw(fhicl::ParameterSet const& pset);
    virtual ~RecoWireICARUSRaw();
    
    void produce(art::Event& evt); 
    void beginJob(); 
    void endJob();                 
    void reconfigure(fhicl::ParameterSet const& p);
    void OfflineIntegration(std::vector<float>&);
    void DoubleRebinning(std::vector<float>&);

 
  private:
    
    int          fDataSize;          ///< size of raw data on one wire
    std::string  fDigitModuleLabel;  ///< module that made digits
    TH1F* fWireRMS;


  protected: 
    
  }; // class RecoWireICARUSRaw
}

namespace recowireraw{

  //-------------------------------------------------
  RecoWireICARUSRaw::RecoWireICARUSRaw(fhicl::ParameterSet const& pset)
  {
    this->reconfigure(pset);
    
    produces< std::vector<recob::Wire> >();
    produces<art::Assns<raw::RawDigit, recob::Wire>>();
      
    
    
  }
  
  //-------------------------------------------------
  RecoWireICARUSRaw::~RecoWireICARUSRaw()
  {
  }

  //////////////////////////////////////////////////////
  void RecoWireICARUSRaw::reconfigure(fhicl::ParameterSet const& p)
  {
    fDigitModuleLabel = p.get< std::string >("DigitModuleLabel", "daq");
  }

  //-------------------------------------------------
  void RecoWireICARUSRaw::beginJob()
  {
      // get access to the TFile service
      art::ServiceHandle<art::TFileService> tfs;
      
      fWireRMS    = tfs->make<TH1F>("fWireRMS", "RMS(ADC#)", 1000, 0, 10);
  }

  //////////////////////////////////////////////////////
  void RecoWireICARUSRaw::endJob()
  {  
  }
  
  //////////////////////////////////////////////////////
  void RecoWireICARUSRaw::produce(art::Event& evt)
  {
    // get the geometry
    art::ServiceHandle<geo::Geometry> geom;
      
    // make a collection of Wires
    std::unique_ptr<std::vector<recob::Wire> > wirecol(new std::vector<recob::Wire>);
    // ... and an association set
    std::unique_ptr<art::Assns<raw::RawDigit,recob::Wire> > WireDigitAssn
      (new art::Assns<raw::RawDigit,recob::Wire>);
    
      std::cout << " before wirecol size " << wirecol->size() << std::endl;
      
    // Read in the digit List object(s). 
    art::Handle< std::vector<raw::RawDigit> > digitVecHandle;
    evt.getByLabel(fDigitModuleLabel, digitVecHandle);

    if (!digitVecHandle->size())  return;
    mf::LogInfo("RecoWireICARUSRaw") << "RecoWireICARUSRaw:: digitVecHandle size is " << digitVecHandle->size();

    // Use the handle to get a particular (0th) element of collection.
    art::Ptr<raw::RawDigit> digitVec0(digitVecHandle, 0);
        
    unsigned int dataSize = digitVec0->Samples(); //size of raw data vectors
    
    raw::ChannelID_t channel(raw::InvalidChannelID); // channel number
    
    unsigned int bin(0);     // time bin loop variable
    
    std::vector<float> holder;                // holds signal data
    std::vector<short> rawadc(dataSize);  // vector holding uncompressed adc values

    wirecol->reserve(digitVecHandle->size()); 
    // loop over all wires
      
    for(unsigned int rdIter = 0; rdIter < digitVecHandle->size(); ++rdIter) { // ++ move
      holder.clear();
      
      art::Ptr<raw::RawDigit> digitVec(digitVecHandle, rdIter);
     channel = digitVec->Channel();
        std::vector<geo::WireID> wids = geom->ChannelToWire(channel);
        geo::WireID wid = wids[0];
//     std::cout << " cryo " << wid.Cryostat << " TPC " << wid.TPC << " view " << wid.Plane << " wire " << wid.Wire << std::endl;
      holder.resize(dataSize);
      
      // uncompress the data
      raw::Uncompress(digitVec->ADCs(), rawadc, digitVec->Compression());
      
      for(bin = 0; bin < dataSize; ++bin) 
	    holder[bin]=(rawadc[bin]-digitVec->GetPedestal());
        
       // std::cout << " pedestal " << digitVec->GetPedestal() << std::endl;
   
        if(wid.Plane==1) {
            float max_raw=0;
            float max_reb=0;
          for(int j=0;j<4096;j++) {
           // std::cout << " before integration sample " << j << " wave " << holder[j] << std::endl;
              if(holder[j]>max_raw) max_raw=holder[j];
          }
             OfflineIntegration(holder);
             DoubleRebinning(holder);
           for(int j=0;j<4096;j++) {
           // std::cout << " after integration sample " << j << " wave " << holder[j] << std::endl;
            if(holder[j]>max_reb) max_reb=holder[j];
        }
        
           
        }
        float media=0;
        float rms2=0;
        for(int js=0;js<4096;js++)
            media+=(holder[js]/4096.);
        for(int js=0;js<4096;js++)
            rms2+=(holder[js]-media)*(holder[js]-media)/4096.;
        float rms=sqrt(rms2);
        fWireRMS->Fill(rms);
        
      wirecol->push_back(recob::WireCreator(holder,*digitVec).move());
        
        // add an association between the last object in wirecol
        // (that we just inserted) and digitVec
        if (!util::CreateAssn(*this, evt, *wirecol, digitVec, *WireDigitAssn)) {
            throw art::Exception(art::errors::ProductRegistrationFailure)
            << "Can't associate wire #" << (wirecol->size() - 1)
            << " with raw digit #" << digitVec.key();
        } // if failed to add association
        //else std::cout << " creating association channel " << channel << std::endl;

        
    	}

      //std::cout << " after wirecol size " << wirecol->size() << std::endl;

      
    if(wirecol->size() == 0)
      mf::LogWarning("RecoWireICARUSRaw") << "No wires made for this event.";
    
    evt.put(std::move(wirecol));
    evt.put(std::move(WireDigitAssn));
    
    return;
  }

    void RecoWireICARUSRaw::OfflineIntegration(std::vector<float>& data)
    {
        //std::cout << " offline integration " << std::endl;
       
        double signal[4096];
               double integral[4096];
        double check[4096];
      
       // double tau0=100;
        for(int j=0;j<4096;j++)
            signal[j]=data[j];
        
        double tau0=100;
        
        integral[0]=signal[0];
        for(int j=1;j<4096;j++)
           integral[j]=integral[j-1]*exp(-0.4/tau0)+signal[j];
        
        double outl[4096];
        double outb[4096];
        double outs[4096];
        double tl=1000;
        double ts=0;
        
        ts=0;
        tl=27;
        
        
        
        
        outl[0]=integral[0];
        double sum_tot=0;
        for(int j=1;j<4096;j++) {
            sum_tot+=integral[j];
            outl[j]=outl[j-1]*exp(-0.4/tl)+integral[j];
        }
        
        
        outb[4095]=integral[4095];
        
        for(int j=4094;j>=0;j--) {
            outb[j]=outb[j+1]*exp(-0.4/tl)+integral[j];
        }
        
        outs[0]=integral[0];
        sum_tot=0;
        for(int j=1;j<4096;j++) {
            sum_tot+=integral[j];
            if(ts>0)
                outs[j]=outs[j-1]*exp(-0.4/ts)+integral[j];
            else
                outs[j]=integral[j];
            
        }
        double nl=(1-exp(-0.4/tl));
        double ns=(1-exp(-0.4/ts));
        
        double intl=0;
        double ints=0;
        double int0=0;
        double intc=0;
        for(int j=0;j<4096;j++) 
            intl+=outl[j];
        for(int j=0;j<4096;j++) 
            ints+=outs[j];
        for(int j=0;j<4096;j++) 
            int0+=integral[j];
        for(int j=0;j<4096;j++) 
            intc+=check[j];
        
        
        for(int j=0;j<4096;j++) { 
            double media=0.5*(outl[j]+outb[j]);
            data[j]=round(outs[j]*(ns*1.)-media*(nl*1.));
            data[j]=round(integral[j]);
        }
        
        
        //cout << " end simulating amplifier " << endl;
    }
    void RecoWireICARUSRaw::DoubleRebinning(std::vector<float>& data)
    {
        int win=12;
        double wout=100;

        
        double output[4096];

        for(int js=0;js<4096;js++)
        {
            int il=js-wout/2;
            if(il<0)
                il=0;
            int fl=js-(win/2);
            if(fl<0)
                fl=0;
            int ir=js+win/2;
            if(ir>4095)
                ir=4095;
            int fr=js+wout/2;
            if(fr>4095)
                fr=4095;
            
            double lreb=0;
            // double insl=js/wout;
            std::vector<int> vout;
            for(int jj=il;jj<=fr;jj++) {
                lreb+=data[jj];
            }
            
            //int middle=(wout)/2;
            lreb/=(wout+1);
            
            double sreb=0;
            for(int jj=fl;jj<ir;jj++)
                sreb+=data[jj];
            sreb/=(win);
            
            if(win==0)
                sreb=data[js];
            double diff;
            diff=sreb-lreb;
            
            //  std::cout << " before output " << std::endl;
            output[js]=diff;

        }
       for(int j=0;j<3996;j++)
           data[j]=output[j];
       for(int j=3996;j<4096;j++)
            data[j]=0;
      
    }
    
    
} // end namespace recowireraw


namespace recowireraw{

  DEFINE_ART_MODULE(RecoWireICARUSRaw)
  
} // end namespace recowireraw


#endif // RecoWireICARUSRawH

