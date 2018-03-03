#ifndef ICARUSHitFINDER_H
#define ICARUSHitFINDER_H

////////////////////////////////////////////////////////////////////
//HIT FINDER THAT RUNS ON RAW SIGNALS INSTEAD OF DECONVOLUTED ONES.
//DEVELOPED INITIALLY FOR ICARUS-T600 OLD ELECTRONICS AT LNGS
//filippo.varanini@pd.infn.it .
////////////////////////////////////////////////////////////////////

// C/C++ standard libraries
#include <string>
#include <vector>
#include <fstream>
#include <set>

//Framework
#include "fhiclcpp/ParameterSet.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h" 
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"


//LArSoft
#include "larcore/Geometry/Geometry.h"
#include "larevt/Filters/ChannelFilter.h"
#include "lardata/Utilities/LArFFT.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/ArtDataHelper/HitCreator.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"

//LArSoft From FFT
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"

#include "larreco/HitFinder/HitFinderTools/ICandidateHitFinder.h"
#include "larreco/HitFinder/HitFinderTools/IPeakFitter.h"

//ROOT from CalData
#include "TComplex.h"
#include "TFile.h"
#include "TH2D.h"
#include "TF1.h"

//ROOT From Gauss
#include "TH1D.h"
#include "TDecompSVD.h"
#include "TMath.h"

namespace hit {

  class ICARUSHitFinder : public art::EDProducer {

    public:

      explicit ICARUSHitFinder(fhicl::ParameterSet const& pset);
      virtual ~ICARUSHitFinder();

      void produce(art::Event& evt); 
      void beginJob(); 
      void endJob(); 
      void reconfigure(fhicl::ParameterSet const& p);
     
      void expandHit(reco_tool::ICandidateHitFinder::HitCandidate& h, std::vector<float> holder, std::vector<reco_tool::ICandidateHitFinder::HitCandidate> how );
      void computeBestLocalMean(reco_tool::ICandidateHitFinder::HitCandidate& h, std::vector<float> holder, std::vector<reco_tool::ICandidateHitFinder::HitCandidate> how, float& localmean);
      
      static Double_t fitf(Double_t *x, Double_t *par);


    private:
      unsigned int  fDataSize;                  //SIZE OF RAW DATA ON ONE WIRE.
      art::InputTag fDigitModuleLabel;          //MODULE THAT MADE DIGITS.
      std::string   fSpillName;                 //NOMINAL SPILL IS AN EMPTY STRING.

      //FFT COPIED VARIABLES.
      std::string         fCalDataModuleLabel;
      std::string         fHitLabelName;
      
      int              fThetaAngle;
       bool                fUncompressWithPed;   //OPTION TO UNCOMPRESS WITH PEDESTAL.
      
      std::unique_ptr<reco_tool::ICandidateHitFinder> fHitFinderTool;  ///< For finding candidate hits
      std::unique_ptr<reco_tool::IPeakFitter>         fPeakFitterTool; ///< Perform fit to candidate peaks
      
      
      //histograms
      TH1F* fFirstChi2;
      TH1F* fChi2;
      TH1F* fHeightC;
      TH1F* fWidthC;
      TH1F* fNoiseC;
      TH1F* fAreaC;
      TH1F* fnhwC;
      TH1F* fIntegralC;
      TH1F* fAreaInt;


      TH1F* fHeightI2;
      TH1F* fWidthI2;
      TH1F* fNoiseI2;
      TH1F* fHeightI1;
      TH1F* fWidthI1;
      TH1F* fNoiseI1;
     
  }; // class ICARUSHitFinder


  //-------------------------------------------------
  ICARUSHitFinder::ICARUSHitFinder(fhicl::ParameterSet const& pset)
  {
    this->reconfigure(pset);

    //LET HITCOLLECTIONCREATOR DECLARE THAT WE ARE GOING TO PRODUCE
    //HITS AND ASSOCIATIONS TO RAW DIGITS BUT NOT ASSOCIATIONS TO WIRES
    //(WITH NO PARTICULAR PRODUCT LABEL).
    recob::HitCollectionCreator::declare_products(*this, 
        /*instance_name*/"");
  }

  //-------------------------------------------------
  ICARUSHitFinder::~ICARUSHitFinder()
  {
  }

  void ICARUSHitFinder::reconfigure(fhicl::ParameterSet const& p)
  {
       fUncompressWithPed  = p.get< bool         >("UncompressWithPed", true);
    fDigitModuleLabel   = p.get< art::InputTag >("DigitModuleLabel", "daq");
    fCalDataModuleLabel = p.get< std::string  >("CalDataModuleLabel");
 
      fThetaAngle=p.get< double  >("ThetaAngle");
      
      fHitFinderTool  = art::make_tool<reco_tool::ICandidateHitFinder>(p.get<fhicl::ParameterSet>("CandidateHits"));
      fPeakFitterTool = art::make_tool<reco_tool::IPeakFitter>(p.get<fhicl::ParameterSet>("PeakFitter"));
      
    mf::LogInfo("ICARUSHitFinder_module") << "fDigitModuleLabel: " << fDigitModuleLabel << std::endl;
  }

  //-------------------------------------------------
  void ICARUSHitFinder::beginJob()
  {
      // get access to the TFile service
      art::ServiceHandle<art::TFileService> tfs;
      
      
      // ======================================
      // === Hit Information for Histograms ===
      fFirstChi2	= tfs->make<TH1F>("fFirstChi2", "#chi^{2}", 10000, 0, 5000);
      fChi2	        = tfs->make<TH1F>("fChi2", "#chi^{2}", 10000, 0, 5000);
      fHeightC	= tfs->make<TH1F>("fHeightC", "height(ADC#)", 100, 0, 100);
      fWidthC	        = tfs->make<TH1F>("fWidthC", "width(samples)", 100, 0, 100);
      fNoiseC	        = tfs->make<TH1F>("fNoiseC", "Noise Area(ADC#)", 100, 0, 100);
      fAreaC            = tfs->make<TH1F>("fAreaC", "", 150, 0, 1500);
      fAreaC->GetXaxis()->SetTitle("Area(ADC#*ticks)");
      fAreaC->GetYaxis()->SetTitle("events");
      fIntegralC            = tfs->make<TH1F>("fIntegralC", "", 500, 0, 3000);
      fIntegralC->GetXaxis()->SetTitle("Area(ADC#*ticks)");
      fIntegralC->GetYaxis()->SetTitle("events");
      fAreaInt = tfs->make<TH1F>("fAreaInt", "Area/Int", 100,0.,2.);

      fnhwC            = tfs->make<TH1F>("fnhwC", "", 10, 0, 10);
      fnhwC->GetXaxis()->SetTitle("n hits per wire");
      fnhwC->GetYaxis()->SetTitle("wires");
      fHeightI2	= tfs->make<TH1F>("fHeightI2", "height(ADC#)", 100, 0, 100);
      fWidthI2	        = tfs->make<TH1F>("fWidthI2", "width(samples)", 100, 0, 100);
      fNoiseI2	        = tfs->make<TH1F>("fNoiseI2", "Noise Area(ADC#)", 100, 0, 100);
      fHeightI1	= tfs->make<TH1F>("fHeightI1", "height(ADC#)", 100, 0, 100);
      fWidthI1	        = tfs->make<TH1F>("fWidthI1", "width(samples)", 100, 0, 100);
      fNoiseI1	        = tfs->make<TH1F>("fNoiseI1", "Noise Area(ADC#)", 100, 0, 100);
      std::cout << " ICARUSHitfinder begin " << std::endl;
  }

  void ICARUSHitFinder::endJob()
  {
      std::cout << " ICARUSHitFinder endjob " << std::endl;
   
  }

  //-------------------------------------------------
  void ICARUSHitFinder::produce(art::Event& evt)
  {      //0
      
      std::cout << " ICARUSHitFinder produce " << std::endl;
      
    //GET THE GEOMETRY.
    art::ServiceHandle<geo::Geometry> geom;
      
      // ###############################################
      // ### Making a ptr vector to put on the event ###
      // ###############################################
      // this contains the hit collection
      // and its associations to wires and raw digits
      
      // Handle the filtered hits collection...
      recob::HitCollectionCreator  hcol(*this, evt);

      //    if (fAllHitsInstanceName != "") filteredHitCol = &hcol;
      
      std::cout << " before reading wirelist " << std::endl;
      
      // ##########################################
      // ### Reading in the Wire List object(s) ###
      // ##########################################
      art::Handle< std::vector<recob::Wire> > wireVecHandle;
      evt.getByLabel(fCalDataModuleLabel,wireVecHandle);
      
      std::cout << " after reading wirelist " << fCalDataModuleLabel << std::endl;

      
      // #################################################################
      // ### Reading in the RawDigit associated with these wires, too  ###
      // #################################################################
      art::FindOneP<raw::RawDigit> RawDigits
      (wireVecHandle, evt, fCalDataModuleLabel);
      
      std::cout << " after reading rawdigits " << std::endl;

      
      // Channel Number
      raw::ChannelID_t channel = raw::InvalidChannelID;
      
      
      std::vector<float> holder;      //HOLDS SIGNAL DATA.
      std::vector<short> rawadc;      //UNCOMPRESSED ADC VALUES.
      
      std::vector<float> startTimes;  //STORES TIME OF WINDOW START.
      std::vector<float> maxTimes;    //STORES TIME OF LOCAL MAXIMUM.
      std::vector<float> endTimes;    //STORES TIME OF WINDOW END.
      std::vector<float> peakHeight;  //STORES ADC COUNT AT THE MAXIMUM.
      std::vector<float> hitrms;    	 //STORES CHARGE WEIGHTED RMS OF TIME ACROSS THE HIT.
      std::vector<double> charge;     //STORES THE TOTAL CHARGE ASSOCIATED WITH THE HIT.
      
      // geo::SigType_t sigType = geo::kInduction;
      std::stringstream numConv;
      
      unsigned int hwC(0),lwC(5728); //lowest and highest wire with physical deposition in Collection
      unsigned int hwI2(0),lwI2(5728); //same for Induction
      unsigned int hwI1(0),lwI1(2112);
      
      unsigned int nw1hitC(0),nw1hitI2(0),nw1hitI1(0); //number of wires (between lowest and highest) with a single hit (rough definition of hitfinding efficiency)
      int nhWire[5600];
      for(int jw=0;jw<5600;jw++)
          nhWire[jw]=0;
      float wCharge[5600];
      for(int jw=0;jw<5600;jw++)
      wCharge[jw]=0;
      float wInt[5600];
      for(int jw=0;jw<5600;jw++)
      wInt[jw]=0;
      
      
      unsigned int nhitsC(0),nhitsI1(0),nhitsI2(0); //total number of reconstructed hits in a view
      
      unsigned int nnhitsC(0),nnhitsI1(0),nnhitsI2(0); //total number of reconstructed hits in a view

      
      unsigned int minWireC,maxWireC;
      if(fThetaAngle==45) {minWireC=2539; maxWireC=3142;}
      if(fThetaAngle==0) {minWireC=2539; maxWireC=4700;}
      if(fThetaAngle==20) {minWireC=2539; maxWireC=4000;}
      if(fThetaAngle==40) {minWireC=2539; maxWireC=3190;}
      if(fThetaAngle==60) {minWireC=2539; maxWireC=2905;}
      if(fThetaAngle==70) {minWireC=2539; maxWireC=2805;}
      if(fThetaAngle==80) {minWireC=2539; maxWireC=2740;}
      

      //45 deg
      /*int minWireC=2539; //empirical
      int maxWireC=3142;
      int minWireI2=2539; //empirical
      int maxWireI2=3102;
      int minWireI1=0; //empirical
      int maxWireI1=497;*/
    
 
      unsigned int minWireI2=2539; //empirical
      unsigned int maxWireI2=4700;
      unsigned int minDrift=850;
      unsigned int maxDrift=1500;
      
      //### Looping over the wires ###
      //##############################
      for(size_t wireIter = 0; wireIter < wireVecHandle->size(); wireIter++)
      {
          // ####################################
          // ### Getting this particular wire ###
          // ####################################
          art::Ptr<recob::Wire>   wire(wireVecHandle, wireIter);
          art::Ptr<raw::RawDigit> rawdigits = RawDigits.at(wireIter);
          
          // --- Setting Channel Number and Signal type ---
          channel = wire->Channel();
          
          std::vector<float> signal(wire->Signal());

          
          // get the WireID for this hit
          std::vector<geo::WireID> wids = geom->ChannelToWire(channel);
          // for now, just take the first option returned from ChannelToWire
          geo::WireID wid  = wids[0];
          // We need to know the plane to look up parameters
          geo::PlaneID::PlaneID_t plane = wid.Plane;
          size_t cryostat=wid.Cryostat;
          size_t tpc=wid.TPC;
size_t iWire=wid.Wire;


      holder.clear();

      //GET THE REFERENCE TO THE CURRENT raw::RawDigit.
      channel   = rawdigits->Channel();
      fDataSize = rawdigits->Samples();

        std::vector<geo::WireID> widVec = geom->ChannelToWire(channel);
        size_t                   iwire  = widVec[0].Wire;
     //   size_t plane = widVec[0].Plane;

        
      rawadc.resize(fDataSize);
      holder.resize(fDataSize);

      //UNCOMPRESS THE DATA.
      if (fUncompressWithPed) {
        int pedestal = (int)rawdigits->GetPedestal();
        raw::Uncompress(rawdigits->ADCs(), rawadc, pedestal, rawdigits->Compression());
      }
      else{
        raw::Uncompress(rawdigits->ADCs(), rawadc, rawdigits->Compression());
      }

      //GET THE LIST OF BAD CHANNELS.
      lariov::ChannelStatusProvider const& channelStatus
        = art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider();

      lariov::ChannelStatusProvider::ChannelSet_t const BadChannels
        = channelStatus.BadChannels();

      for(unsigned int bin = 0; bin < fDataSize; ++bin){ 
        //holder[bin]=(rawadc[bin]-rawdigits->GetPedestal());
          holder[bin]=signal[bin];
          if(plane == 0) holder[bin]=-holder[bin];
          //if(plane == 1) holder[bin]=-holder[bin];
         // if(plane==2&&iwire==2600&&abs(holder[bin])>0.1) std::cout << " wire 2600 bin " << bin << " signal " << holder[bin] << std::endl;
      }

        if(plane==0&&iwire<lwI1) lwI1=iwire;
        if(plane==0&&iwire>hwI1) hwI1=iwire;
        if(plane==1&&iwire<lwI2) lwI2=iwire;
        if(plane==1&&iwire>hwI2) hwI2=iwire;
        if(plane==2&&iwire<lwC) lwC=iwire;
        if(plane==2&&iwire>hwC) hwC=iwire;
        
        
        
     // sigType = geom->SignalType(channel);

      peakHeight.clear();
      endTimes.clear();
      startTimes.clear();
      maxTimes.clear();
      charge.clear();
      hitrms.clear();
        
             bool channelSwitch = false;

      for(auto it = BadChannels.begin(); it != BadChannels.end(); it++)
      {
        if(channel==*it)
        {
          channelSwitch = true;
          break;
        }
      }

      if(channelSwitch==false)
      { //1
          // Hit finding parameters
      int    hitIndex(0);                      //INDEX OF CURRENT HIT IN SEQUENCE.
          double amplitude(0);       //FIT PARAMETERS.
     // double start(0), end(0);
     // double amplitudeErr(0), positionErr(0);  //FIT ERRORS.
     // double goodnessOfFit(0),
         double  chargeErr(0);   //CHI2/NDF and error on charge.
     // double hrms(0);
      //double totSig;
     

    
      reco_tool::ICandidateHitFinder::HitCandidateVec      hitCandidateVec;

      reco_tool::ICandidateHitFinder::MergeHitCandidateVec mergedCandidateHitVec;
          
    std::vector<geo::WireID> wids = geom->ChannelToWire(channel);
      
          
          fHitFinderTool->findHitCandidates(holder, plane,double(iwire),hitCandidateVec);
          for(auto& hitCand : hitCandidateVec) {
            expandHit(hitCand,holder,hitCandidateVec);
           // std::cout << " after expand start " << hitCand.startTick << " stop " << hitCand.stopTick << std::endl;
          }
          //hitCand->localmean=computeBestLocalMean(hitCand,holder,hitCandidateVec);
          fHitFinderTool->MergeHitCandidates(holder, hitCandidateVec, mergedCandidateHitVec);
         // std::cout << " line 408 " << std::endl;
          
          
      //numHits = hits.size();
          int nghC=0;
          int nghI2=0;
          int nghI1=0;
          
      //  if(cryostat==0&&tpc==0&&plane==2)
        //     std::cout << "  Wire " << iwire << " numhits " << hitCandidateVec.size() <<" mergedhits " << mergedCandidateHitVec.size() << std::endl;
          for(auto& mergedCands : mergedCandidateHitVec)
      {
        //  std::cout << " mergedcands size " << mergedCands.size() <<  std::endl;
          
//int startT= mergedCands.front().startTick;
//int endT  = mergedCands.back().stopTick;
       //   std::cout << " after startT " << std::endl;

       
        
        // if(cryostat!=0||tpc!=0) continue;
        
      
        
          // ##################################################
          // ### Calling the function for fitting ICARUS ###
          // ##################################################
          double                                chi2PerNDF(0.);
          int                                   NDF(1);
          reco_tool::IPeakFitter::PeakParamsVec peakParamsVec;
          
        //  std::cout << " before fit " << std::endl;
          
          fPeakFitterTool->findPeakParameters(signal, mergedCands, peakParamsVec, chi2PerNDF, NDF);
         // std::cout << " after fit " << std::endl;

          if (!(chi2PerNDF < std::numeric_limits<double>::infinity()))
          {
              chi2PerNDF = 200.;
              NDF        = 2;
          }
          fFirstChi2->Fill(chi2PerNDF);
          int jhit=0;
          for(const auto& peakParams : peakParamsVec)
          {
              // Extract values for this hit
              float peakAmp   = peakParams.peakAmplitude;
              float peakMean  = peakParams.peakCenter;
              float peakWidth = peakParams.peakSigma;
              float peakLeft   = peakParams.peakTauLeft;
              float peakRight  = peakParams.peakTauRight;
              float peakBaseline = peakParams.peakBaseline;

         
           //   std::cout << " fitted parameters " << peakAmp << " " << peakMean << " " << peakWidth << std::endl;
              
              // Place one bit of protection here
              if (std::isnan(peakAmp))
              {
                //  std::cout << "**** hit peak amplitude is a nan! Channel: " << channel << ", start tick: " << startT << std::endl;
                  continue;
              }
              
              // Extract errors
              float peakAmpErr   = peakParams.peakAmplitudeError;
             float peakMeanErr  = peakParams.peakCenterError;
            //  float peakWidthErr = peakParams.peakSigmaError;
          
              //float fitCharge=chargeFunc(peakMean, peakAmp, peakWidth, fAreaNormsVec[plane],startT,endT);
              //float fitChargeErr = std::sqrt(TMath::Pi()) * (peakAmpErr*peakWidthErr + peakWidthErr*peakAmpErr);
              int start=mergedCands[jhit].startTick;
              int end=mergedCands[jhit].stopTick;

              TF1 Func("ICARUSfunc",fitf,end-start,start,end);
             // for(int jf=0;jf<5;jf++)
                  Func.SetParameter(0,peakBaseline);
              Func.SetParameter(1,peakAmp);
              Func.SetParameter(2,peakMean);
              Func.SetParameter(3,peakLeft);
              Func.SetParameter(4,peakRight);

            

              float fitCharge=Func.Integral(start,end)/2.5;
              float totSig=std::accumulate(holder.begin() + (int) start, holder.begin() + (int) end, 0.);
         //     float fitChargeErr=0;
              
             // std::cout << " fitCharge " << fitCharge << std::endl;
             // std::cout << " totSig " << totSig << std::endl;

//std::cout << " before hit creator " << std::endl;
        recob::HitCreator hit(
            *wire,                                                                     //RAW DIGIT REFERENCE.
            wid,                                                                           //WIRE ID.
            start,                                                                         //START TICK.
            end,                                                                           //END TICK. 
                (peakLeft+peakRight)/2.,                                                                          //RMS.
            peakMean,                                                                      //PEAK_TIME.
            peakMeanErr,                                                                   //SIGMA_PEAK_TIME.
            peakAmp,                                                                     //PEAK_AMPLITUDE.
            peakAmpErr,                                                                  //SIGMA_PEAK_AMPLITUDE.
            fitCharge,                                                                        //HIT_INTEGRAL.
            chargeErr,                                                                     //HIT_SIGMA_INTEGRAL.
            totSig, //SUMMED CHARGE.
            mergedCands.size(),                                                                             //MULTIPLICITY.
            jhit,                                                                            //LOCAL_INDEX.
            chi2PerNDF,                                                                 //WIRE ID.
            NDF                                                               //DEGREES OF FREEDOM.
            );
          
         // std::cout << " before emplace back " << std::endl;
   //       filteredHitVec.push_back(hit.copy());
          
        hcol.emplace_back(hit.move(), wire, rawdigits);
          
              bool driftWindow=(peakMean)>=minDrift&&(peakMean)<=maxDrift;
              bool wireWindowC=iWire>=minWireC&&iWire<=maxWireC;
              bool wireWindowI2=iWire>=minWireI2&&iWire<=maxWireI2;
              bool outDriftWindow=(peakMean)<=minDrift||(peakMean)>=maxDrift;
              bool outWireWindowC=iWire<=minWireC||iWire>=maxWireC;
              bool outWireWindowI2=iWire<=minWireI2||iWire>=maxWireI2;
    
       
          if(plane==0&&driftWindow)  {
           //   std::cout << " wire " << hits[i].iWire << " ngh " << ngh << std::endl;
              nghI1++;
              if(nghI1==1) nw1hitI1++;
      
             //  std::cout << " wire " << hits[i].iWire << " nw1h " << nw1 << std::endl;
              fHeightI1->Fill(amplitude);
              fWidthI1->Fill(end-start);
          }
          if(plane==1&&wireWindowI2) {
           //   std::cout << " wire " << hits[i].iWire << " ngh " << ngh << std::endl;
              nghI2++;
              if(nghI2==1) nw1hitI2++;
             // std::cout << " filling height histo wire" << hits[i].iWire << " drift " << hits[i].hitCenter << " amplitude " << amplitude << std::endl;
              fHeightI2->Fill(amplitude);
              fWidthI2->Fill(end-start);
          }
          if(plane==2&&wireWindowC) {
              nghC++;
              if(nghC==1) nw1hitC++;
              nhWire[iWire]++;
              fHeightC->Fill(peakAmp);
              fWidthC->Fill(peakWidth);
             // fAreaC->Fill(totSig);
              wCharge[iWire]+=fitCharge;
              wInt[iWire]+=totSig;

           //   if(intSig<0.)
             // std::cout << "  intsig " << intSig << " localmean " << hits[i].localmean << " wire " << hits[i].iWire << std::endl;
          }

          if(plane==0) nhitsI1++;
          if(plane==1) nhitsI2++;
          if(plane==2) nhitsC++;
          
          if(plane==0&&tpc==0&&cryostat==0&&outDriftWindow) {nnhitsI1++; fNoiseI1->Fill(amplitude); }
          if(plane==1&&tpc==0&&cryostat==0&&outWireWindowI2) { nnhitsI2++; fNoiseI2->Fill(amplitude);
            //if(cryostat==0&&tpc==0)
              //    std::cout << " noise hit Wire " << hits[i].iWire << " tick " << hits[i].hitCenter << std::endl;
          }
          if(plane==2&&tpc==0&&cryostat==0&&outWireWindowC) { nnhitsC++; fNoiseC->Fill(amplitude); }
        
          
        ++hitIndex;
      } //end loop on found hits
         // if(plane==2&&cryostat==0&&tpc==0&&wCharge[iwire]>0.)
          // std::cout << " filling  wire  " << iwire << " area " << wCharge[iwire] << std::endl;
          if(plane==2&&cryostat==0&&tpc==0)
           if(wCharge[iwire]>0)
           fAreaC->Fill(wCharge[iwire]);
          if(plane==2&&cryostat==0&&tpc==0)
          if(wInt[iwire]>0)
          fIntegralC->Fill(wInt[iwire]);
          if(wInt[iwire]>0&&cryostat==0&&tpc==0)
              fAreaInt->Fill(wCharge[iwire]/wInt[iwire]);
          //if(wCharge[iwire]>0&&cryostat==0&&tpc==0&&wCharge[iwire]/wInt[iwire]<0.9)
            //  std::cout << " wire " << iwire << " ratio " << wCharge[iwire]/wInt[iwire] << std::endl;
      }
    } //end channel condition

    } //end loop on channels
      
            double PhysWC=maxWireC-minWireC+1;
             double NoiseWC=5119-480-(maxWireC-minWireC+1);
      double PhysWI2=maxWireI2-minWireI2+1;
      double NoiseWI2=5119-480-(maxWireI2-minWireI2+1);
      double PhysWI1=1056;
      double NoiseWI1=1056;
    
      std::cout << " Physical collection wires " << PhysWC << " single hit wires " << nw1hitC << " efficiency " << float(nw1hitC)/PhysWC << std::endl;
      std::cout << " Average collection noise hits per wire " << nnhitsC << " " << NoiseWC << std::endl;
      std::cout << " Physical ind2 wires " << PhysWI2 << " single hit wires " << nw1hitI2 << " noise " << nnhitsI2 << " efficiency " << float(nw1hitI2)/PhysWI2 << std::endl;
      std::cout << " Average ind2 noise hits per wire " << nnhitsI2<< " " << NoiseWI2 << std::endl;
      std::cout << " Physical ind1 wires " << PhysWI1 << " single hit wires " << nw1hitI1 << " efficiency " << float(nw1hitI1)/PhysWI1 << std::endl;
      std::cout << " Average ind1 noise hits per wire " << nnhitsI1/NoiseWI1 << std::endl;

      for(unsigned int jw=minWireC;jw<maxWireC;jw++)
          fnhwC->Fill(nhWire[jw]);
      
      
    hcol.put_into(evt);
      std::cout << " end ICARUSHitfinder " << std::endl;
     
      
  } //end produce

    
    void ICARUSHitFinder::expandHit(reco_tool::ICandidateHitFinder::HitCandidate& h, std::vector<float> holder, std::vector<reco_tool::ICandidateHitFinder::HitCandidate> how)
    {
        // Given a hit or hit candidate <hit> expand its limits to the closest minima
        int nsamp=50;
        int cut=1;
        int upordown;
        int found=0;
        
        unsigned int first=h.startTick;
        unsigned int last =h.stopTick;
        
        std::vector<reco_tool::ICandidateHitFinder::HitCandidate>::iterator hiter;
        std::vector<reco_tool::ICandidateHitFinder::HitCandidate> hlist;
        
        // fill list of existing hits on this wire
        for(unsigned int j=0;j<how.size();j++)
        {
            reco_tool::ICandidateHitFinder::HitCandidate h2=how[j];
            if(h2.hitCenter!=h.hitCenter)
            hlist.push_back(h2);
        }
        
        
        // look for first sample
        while(!found)
        {
            if(first==0)
                break;
            
            for(hiter=hlist.begin();hiter!=hlist.end();hiter++)
            {
                reco_tool::ICandidateHitFinder::HitCandidate h2=*hiter;
                if(first==h2.stopTick)
                  found=1;
            }
            
            if(found==1) break;
            
            upordown=0;
            for(int l=0;l<nsamp/2;l++)
            {
                if(first-nsamp/2+l<4095) {
                    if(holder[first-nsamp/2+l+1]-holder[first-nsamp/2+l]>0) upordown++;
                    else if(holder[first-nsamp/2+l+1]-holder[first-nsamp/2+l]<0) upordown--;
            } }
            //std::cout << " checking " << first << " upordown " << upordown << std::endl;
            if(upordown>cut)
                first--;
            else
                found=1;
        }
        
        // look for last sample
        found=0;
        while(!found)
        {
            if(last==4095)
            {
                found=1;
                break;
            }
            
            for(hiter=hlist.begin();hiter!=hlist.end();hiter++)
            {
                reco_tool::ICandidateHitFinder::HitCandidate h2=*hiter;
                if(first==h2.startTick)
                found=1;
            }
            
            if(found==1) break;
            
            upordown=0;
            for(int l=0;l<nsamp/2;l++)
            {
                if(last+nsamp/2-l>0) {
                    if(holder[last+nsamp/2-l]-holder[last+nsamp/2-l-1]>0) upordown++;
                    else if(holder[last+nsamp/2-l]-holder[last+nsamp/2-l-1]<0) upordown--;
                } }
            
            if(upordown<-cut)
                last++;
            else
                found=1;
        }
        
        h.startTick=first;
        h.stopTick=last;
    }
    void ICARUSHitFinder::computeBestLocalMean(reco_tool::ICandidateHitFinder::HitCandidate& h, std::vector<float> holder, std::vector<reco_tool::ICandidateHitFinder::HitCandidate> how, float& localmean)
    {
        const int bigw=130;   //size of the window where to look for the minimum localmean value
        const int meanw=70;   //size of the window where the mean is calculated
        const int outofbounds=9999;
        
        float samples1[bigw];   //list to contain samples bellow the startTick
        float samples2[bigw];   //list to contain samples above the stopTick
        
        std::vector<reco_tool::ICandidateHitFinder::HitCandidate>::iterator hiter;
        std::vector<reco_tool::ICandidateHitFinder::HitCandidate> hlist;

        float min1;
        float min2;
        int startTick=h.startTick;
        int stopTick=h.stopTick;
        float count1=0;
        float count2=0;
        unsigned int shift1=0;
        unsigned int shift2=0;
        int foundborder1=0;
        int foundborder2=0;
        
        // fill list of existing hits on this wire
        for(unsigned int j=0;j<how.size();j++)
        {
            reco_tool::ICandidateHitFinder::HitCandidate h2=how[j];
            if(h2.hitCenter!=h.hitCenter)
            hlist.push_back(h2);
        }
        
        // fill the arrays of samples to be examined
        for(unsigned int i=0;i<bigw;i++)
        {
            // remove samples from other hits in the wire
            for(hiter=hlist.begin();hiter!=hlist.end();hiter++)
            {
                reco_tool::ICandidateHitFinder::HitCandidate h2=*hiter;
                if(startTick-i-shift1 >= h2.startTick && startTick-i-shift1 <= h2.stopTick)
                shift1+=h2.stopTick-h2.startTick+1;
                else if(stopTick+i+shift2 >= h2.startTick && stopTick+i+shift2 <= h2.stopTick)
                shift2+=h2.stopTick-h2.startTick+1;
            }
            
            // fill the lists
            if(startTick-i-shift1>=0)
            samples1[i]=holder[h.startTick-i-shift1];
            else
            samples1[i]=outofbounds;
            
            if(stopTick+i+shift2<=4095)
            samples2[i]=holder[h.stopTick+i+shift2];
            else
            samples2[i]=outofbounds;
        }
        
        // initialize counters
        for(int j=0;j<meanw;j++)
        {
            if(samples1[j]==outofbounds) foundborder1=1;
            if(samples2[j]==outofbounds) foundborder2=1;
            
            if(!foundborder1) count1+=samples1[j];
            if(!foundborder2) count2+=samples2[j];
        }
        min1=count1;
        min2=count2;
        
        //look for the local minima from the filled lists
        for(int j=0;j<bigw-meanw;j++)
        {
            if(count1<min1) min1=count1;
            if(samples1[j+meanw]==outofbounds) foundborder1=1;
            if(!foundborder1) count1+=samples1[j+meanw]-samples1[j];
            
            if(count2<min2) min2=count2;
            if(samples2[j+meanw]==outofbounds) foundborder2=1;
            if(!foundborder2) count2+=samples2[j+meanw]-samples2[j];
            
            if(foundborder1 && foundborder2) break;
        }
        if(count1<min1) min1=count1;
        if(count2<min2) min2=count2;
        
       // std::cout << " localmean count1 " << count1 << " count2 " << count2 << std::endl;
        //for the moment take the highest value
        if((foundborder1 && !foundborder2)  || (shift1 && !shift2))
        localmean=((float) min2)/meanw;
        else if((!foundborder1 && foundborder2) || (!shift1 && shift2))
        localmean=((float) min1)/meanw;
        else
        localmean=(min1>min2)?((float) min1)/meanw :((float) min2)/meanw;
        
        //for the moment take the average value
        //  h->localmean=(min1+min2)/(2.*meanw);
        
        
    }
Double_t ICARUSHitFinder::fitf(Double_t *x, Double_t *par)
    {
        // Double_t arg = 0;
        
        Double_t fitval = par[0]+par[1]*TMath::Exp(-(x[0]-par[2])/par[3])/(1+TMath::Exp(-(x[0]-par[3])/par[4]));
        return fitval;
    }
  
  DEFINE_ART_MODULE(ICARUSHitFinder)

} // end namespace hit


#endif //ICARUSHitFinder_H
