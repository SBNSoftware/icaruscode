////////////////////////////////////////////////////////////////////////
// $Id: RecoWireICARUS.cxx,v 1.36 2010/09/15  bpage Exp $
//
// RecoWireICARUS class
//
// brebel@fnal.gov
//
////////////////////////////////////////////////////////////////////////
#ifndef RecoWireICARUS_H
#define RecoWireICARUS_H

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

// LArSoft includes
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardata/ArtDataHelper/WireCreator.h"
#include "lardata/Utilities/LArFFT.h"
#include "lardata/Utilities/AssociationUtil.h"

///creation of calibrated signals on wires
namespace recowire {

  class RecoWireICARUS : public art::EDProducer {

  public:
    
    // create calibrated signals on wires. this class runs 
    // an fft to remove the electronics shaping.     
    explicit RecoWireICARUS(fhicl::ParameterSet const& pset); 
    virtual ~RecoWireICARUS();
    
    void produce(art::Event& evt); 
    void beginJob(); 
    void endJob();                 
    void reconfigure(fhicl::ParameterSet const& p);
 
  private:
    
    std::string  fResponseFile;      ///< response file containing transformed 
                                     ///< shape histograms and decay constants
    int          fDataSize;          ///< size of raw data on one wire
    int          fExpEndBins;        ///< number of end bins to consider for tail fit    
    int          fPostsample;        ///< number of postsample bins
    std::string  fDigitModuleLabel;  ///< module that made digits

    std::vector<std::vector<TComplex> > fKernelR;      ///< holds transformed induction 
                                                       ///< response function
    std::vector<std::vector<TComplex> > fKernelS;      ///< holds transformed induction 
                                                       ///< response function
    std::vector<double>                 fDecayConstsR; ///< vector holding RC decay 
                                                       ///< constants
    std::vector<double>                 fDecayConstsS; ///< vector holding RC decay 
                                                       ///< constants
    std::vector<int>                    fKernMapR;     ///< map telling which channels  
                                                       ///< have which response functions
    std::vector<int>                    fKernMapS;     ///< map telling which channels 
                                                       ///< have which response functions
  protected: 
    
  }; // class RecoWireICARUS
}

namespace recowire{

  //-------------------------------------------------
  RecoWireICARUS::RecoWireICARUS(fhicl::ParameterSet const& pset)
  {
    this->reconfigure(pset);
    
    produces< std::vector<recob::Wire> >();
    produces<art::Assns<raw::RawDigit, recob::Wire>>();
    
  }
  
  //-------------------------------------------------
  RecoWireICARUS::~RecoWireICARUS()
  {
  }

  //////////////////////////////////////////////////////
  void RecoWireICARUS::reconfigure(fhicl::ParameterSet const& p)
  {
    fDigitModuleLabel = p.get< std::string >("DigitModuleLabel", "daq");
    cet::search_path sp("FW_SEARCH_PATH");
    sp.find_file(p.get<std::string>("ResponseFile"), fResponseFile);
    fExpEndBins       = p.get< int >        ("ExponentialEndBins");
    fPostsample       = p.get< int >        ("PostsampleBins");
  }

  //-------------------------------------------------
  void RecoWireICARUS::beginJob()
  {  
    
    LOG_DEBUG("RecoWireICARUS") << "RecoWireICARUS_plugin: Opening  Electronics Response File: " 
			 << fResponseFile.c_str();
    
    TFile f(fResponseFile.c_str());
    if( f.IsZombie() )
      mf::LogWarning("RecoWireICARUS") << "Cannot open response file " 
				<< fResponseFile.c_str();
    
    TH2D *respRe       = dynamic_cast<TH2D*>(f.Get("real/RespRe")   );
    TH2D *respIm       = dynamic_cast<TH2D*>(f.Get("real/RespIm")   );
    TH1D *decayHist    = dynamic_cast<TH1D*>(f.Get("real/decayHist"));
    unsigned int wires = decayHist->GetNbinsX();
    unsigned int bins  = respRe->GetYaxis()->GetNbins();
    unsigned int bin   = 0;
    unsigned int wire  = 0;
    fDecayConstsR.resize(wires);
    fKernMapR.resize(wires);
    fKernelR.resize(respRe->GetXaxis()->GetNbins());
    const TArrayD *edges = respRe->GetXaxis()->GetXbins();
    for(int i = 0; i < respRe->GetXaxis()->GetNbins(); ++i)  {
      fKernelR[i].resize(bins);
      for(bin = 0; bin < bins; ++bin) {  
        
	const TComplex a(respRe->GetBinContent(i+1,bin+1), 
			 respIm->GetBinContent(i+1,bin+1));
	fKernelR[i][bin]=a;
      }
      for(; wire < (*edges)[i+1]; ++wire) {
	fKernMapR[wire]=i;
        fDecayConstsR[wire]=decayHist->GetBinContent(wire+1); 
      }
    }
    respRe       = dynamic_cast<TH2D*>(f.Get("sim/RespRe")   );
    respIm       = dynamic_cast<TH2D*>(f.Get("sim/RespIm")   );
    decayHist    = dynamic_cast<TH1D*>(f.Get("sim/decayHist"));
    wires = decayHist->GetNbinsX();
    bins  = respRe->GetYaxis()->GetNbins();
    fDecayConstsS.resize(wires);
    fKernMapS.resize(wires);
    fKernelS.resize(respRe->GetXaxis()->GetNbins()); 
    const TArrayD *edges1 = respRe->GetXaxis()->GetXbins();
    wire =0;
    for(int i = 0; i < respRe->GetXaxis()->GetNbins(); ++i)  {
      fKernelS[i].resize(bins);
      for(bin = 0; bin < bins; ++bin) {  
	const TComplex b(respRe->GetBinContent(i+1,bin+1), 
			 respIm->GetBinContent(i+1,bin+1));
	fKernelS[i][bin]=b;
      }
      for(; wire < (*edges1)[i+1]; ++wire) {
	fKernMapS[wire]=i;
        fDecayConstsS[wire]=decayHist->GetBinContent(wire+1); 
      }
    }
   
    f.Close();
  }

  //////////////////////////////////////////////////////
  void RecoWireICARUS::endJob()
  {  
  }
  
  //////////////////////////////////////////////////////
  void RecoWireICARUS::produce(art::Event& evt)
  {
    
      
    // get the geometry
    art::ServiceHandle<geo::Geometry> geom;

    std::vector<double> decayConsts;  
    std::vector<int> kernMap;
    std::vector<std::vector<TComplex> > kernel; 
    //Put correct response functions and decay constants in place
    if(evt.isRealData()) {
      decayConsts=fDecayConstsR;
      kernMap=fKernMapR;
      kernel=fKernelR;
    }
    else {
      decayConsts=fDecayConstsS;
      kernMap=fKernMapS;
      kernel=fKernelS;
    }

    // get the FFT service to have access to the FFT size
    art::ServiceHandle<util::LArFFT> fFFT;

    // make a collection of Wires
    std::unique_ptr<std::vector<recob::Wire> > wirecol(new std::vector<recob::Wire>);
    // ... and an association set
    std::unique_ptr<art::Assns<raw::RawDigit,recob::Wire> > WireDigitAssn
      (new art::Assns<raw::RawDigit,recob::Wire>);
    
    // Read in the digit List object(s). 
    art::Handle< std::vector<raw::RawDigit> > digitVecHandle;
    evt.getByLabel(fDigitModuleLabel, digitVecHandle);

    if (!digitVecHandle->size())  return;
    mf::LogInfo("RecoWireICARUS") << "RecoWireICARUS:: digitVecHandle size is " << digitVecHandle->size();

    // Use the handle to get a particular (0th) element of collection.
    art::Ptr<raw::RawDigit> digitVec0(digitVecHandle, 0);
        
    unsigned int dataSize = digitVec0->Samples(); //size of raw data vectors
    
    int transformSize = fFFT->FFTSize();
    raw::ChannelID_t channel(raw::InvalidChannelID); // channel number
    unsigned int bin(0);     // time bin loop variable
    
    double decayConst = 0.;  // exponential decay constant of electronics shaping
    double fitAmplitude    = 0.;  //This is the seed value for the amplitude in the exponential tail fit 
    std::vector<float> holder;                // holds signal data
    std::vector<float> shortADC;                // holds signal data
    std::vector<float> longADC;                // holds signal data
    std::vector<short> rawadc(transformSize);  // vector holding uncompressed adc values
    std::vector<TComplex> freqHolder(transformSize+1); // temporary frequency data
    wirecol->reserve(digitVecHandle->size()); 
    // loop over all wires    
    for(unsigned int rdIter = 0; rdIter < digitVecHandle->size(); ++rdIter){ // ++ move
      holder.clear();
      
      art::Ptr<raw::RawDigit> digitVec(digitVecHandle, rdIter);
      channel = digitVec->Channel();
      //std::cout<<channel<<std::endl;
      holder.resize(transformSize);
      shortADC.resize(transformSize);
      longADC.resize(transformSize);
      
      // uncompress the data
      raw::Uncompress(digitVec->ADCs(), rawadc, digitVec->Compression());
      
      for(bin = 0; bin < dataSize; ++bin) 
	holder[bin]=(rawadc[bin]-digitVec->GetPedestal());
      // fExpEndBins only nonzero for detectors needing exponential tail fitting
      geo::SigType_t sigtype = geom->SignalType(channel);
      size_t k;
      if(sigtype == geo::kInduction)
	k = 0;
      else if(sigtype == geo::kCollection)
	k = 1;
      else
	throw cet::exception("RecoWireICARUS") << "Bad signal type = " << sigtype << "\n";
      if (k >= kernel.size())
	throw cet::exception("RecoWireICARUS") << "kernel size < " << k << "!\n";

      if(k==0)
	{
        for(bin = 1; bin < dataSize; ++bin)
        {
            shortADC[bin]=0;
        }
       for(bin = 0; bin < dataSize; ++bin)
        {
            if (holder[bin]>9) {
                for (int ijk=0; ijk<30; ijk++) {
                    if ((bin-ijk)>=0 && (bin+ijk)<4096) {
                        shortADC[bin-ijk]=1;
                        shortADC[bin+ijk]=1;
                    }
                }
            }
        }
        for(bin = 0; bin < dataSize; ++bin)
        {
            if(shortADC[bin]==0)holder[bin]=0;
        }

        for(bin = 1; bin < dataSize; ++bin)
        {
           holder[bin]+= holder[bin-1];
        }
        for(bin = 0; bin < dataSize; ++bin)
        {
            if(shortADC[bin]==0)holder[bin]=0;
        }
        
        /*
        for(bin = 1; bin < dataSize; ++bin)
      {
	    holder[bin]+= holder[bin-1];
        longADC[bin]=0;shortADC[bin]=0;
      }
        
	  for(bin = 128; bin < dataSize; ++bin)
	    {
	      longADC[bin]=0;
	      for(int ijk=0;ijk<128;ijk++)
		{
		  longADC[bin]+=(holder[bin-ijk])/128.;
		}
	    }
        for(bin = 8; bin < dataSize; ++bin)
        {
            shortADC[bin]=0;
            for(int ijk=0;ijk<8;ijk++)
            {
                shortADC[bin]+=(holder[bin-ijk])/8.;
            }
        }
	  
	  //for(bin = 0; bin < dataSize; ++bin)
	    //holder[bin]= shortADC[bin]-longADC[bin];
        for(bin = 0; bin < dataSize; ++bin)holder[bin]=holder[bin]*0.1;
	  
	  //fFFT->Convolute(holder,kernel[k]);
 */
	}

      if(k==2)
	{
	  if(fExpEndBins && std::abs(decayConsts[channel]) > 0.0){
	
	    TH1D expTailData("expTailData","Tail data for fit",
			     fExpEndBins,dataSize-fExpEndBins,dataSize);
	    TF1  expFit("expFit","[0]*exp([1]*x)");
	    
	    for(bin = 0; bin < (unsigned int)fExpEndBins; ++bin) 
	      expTailData.Fill(dataSize-fExpEndBins+bin,holder[dataSize-fExpEndBins+bin]);
	    decayConst = decayConsts[channel];
	    fitAmplitude = holder[dataSize-fExpEndBins]/exp(decayConst*(dataSize-fExpEndBins));
	    expFit.FixParameter(1,decayConst);
	    expFit.SetParameter(0,fitAmplitude);
	    expTailData.Fit(&expFit,"QWN","",dataSize-fExpEndBins,dataSize);
	    expFit.SetRange(dataSize,transformSize);
	    for(bin = 0; bin < dataSize; ++bin)
	      holder[dataSize+bin]= expFit.Eval(bin+dataSize);
	  }
	  // This is actually deconvolution, by way of convolution with the inverted 
	  // kernel.  This code assumes the response function has already been
	  // been transformed and inverted.  This way a complex multiplication, rather
	  // than a complex division is performed saving 2 multiplications and 
	  // 2 divsions
      
	  // the example below is for MicroBooNE, experiments should 
	  // adapt as appropriate
      
	  // Figure out which kernel to use (0=induction, 1=collection).
      
	  fFFT->Convolute(holder,kernel[k]);
    
	  holder.resize(dataSize,1e-5);
	  //This restores the DC component to signal removed by the deconvolution.
	  if(fPostsample) {
	    double average=0.0;
	    for(bin=0; bin < (unsigned int)fPostsample; ++bin) 
	      average+=holder[holder.size()-1-bin]/(double)fPostsample;
	    for(bin = 0; bin < holder.size(); ++bin) holder[bin]-=average;
	  }
	}
      wirecol->push_back(recob::WireCreator(holder,*digitVec).move());
      // add an association between the last object in wirecol
      // (that we just inserted) and digitVec
      if (!util::CreateAssn(*this, evt, *wirecol, digitVec, *WireDigitAssn)) {
        throw art::Exception(art::errors::ProductRegistrationFailure)
          << "Can't associate wire #" << (wirecol->size() - 1)
          << " with raw digit #" << digitVec.key();
      } // if failed to add association
    } // for raw digits

    if(wirecol->size() == 0)
      mf::LogWarning("RecoWireICARUS") << "No wires made for this event.";
    
    evt.put(std::move(wirecol));
    evt.put(std::move(WireDigitAssn));
    
    return;
  }
  
} // end namespace recowire


namespace recowire{

  DEFINE_ART_MODULE(RecoWireICARUS)
  
} // end namespace recowire


#endif // RecoWireICARUS_H

