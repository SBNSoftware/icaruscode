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
//#include "icaruscode/HitFinder/PeakFitterICARUS.h"

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
      static Double_t fitlong(Double_t *x, Double_t *par);

      using ICARUSPeakFitParams_t = struct ICARUSPeakFitParams
      {
          float peakCenter;
          float peakCenterError;
          float peakSigma;
          float peakSigmaError;
          float peakAmplitude;
          float peakAmplitudeError;
          float peakTauLeft;
          float peakTauLeftError;
          float peakTauRight;
          float peakTauRightError;
          float peakFitWidth;
          float peakFitWidthError;
          float peakSlope;
          float peakSlopeError;
          float peakBaseline;
          float peakBaselineError;
      };
      using ICARUSPeakParamsVec = std::vector<ICARUSPeakFitParams_t>;
      void findMultiPeakParameters(const std::vector<float>&,
                                   const reco_tool::ICandidateHitFinder::HitCandidateVec&,
                                   ICARUSPeakParamsVec&,
                                   double&,
                                   int&, int) const;
      void findLongPeakParameters(const std::vector<float>&,
                                  const reco_tool::ICandidateHitFinder::HitCandidateVec&,
                                  ICARUSPeakParamsVec&,
                                  double&,
                                  int&, int) const;
      double ComputeChiSquare(TF1 func, TH1 *histo) const;
      double ComputeNullChiSquare(std::vector<float>) const;


      void setWire(int i) {
          iWire=i;
       //   std::cout << " setting iwire " << iWire << std::endl;
      } ;
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
   //  PeakFitterICARUS*         fPeakFitterTool; ///< Perform fit to candidate peaks
      
      size_t              fMaxMultiHit;              ///<maximum hits for multi fit
double                fChi2NDF;                  ///maximum Chisquared / NDF allowed for a hit to be saved
      std::vector<int>    fLongMaxHitsVec;           ///<Maximum number hits on a really long pulse train
      std::vector<int>    fLongPulseWidthVec;        ///<Sets width of hits used to describe long pulses
      
      //histograms
      TH1F* fFirstChi2;
      TH1F* fNullChi2;

      TH1F* fChi2;
      TH1F* fHeightC;
      TH1F* fWidthC;
      TH1F* fNoiseC;
      TH1F* fAreaC;
      TH1F* fnhwC;
      TH1F* fIntegralC;
      TH1F* fAreaInt;

      TH1F* fWire2771;

      TH1F* fHeightI2;
      TH1F* fWidthI2;
      TH1F* fNoiseI2;
      TH1F* fHeightI1;
      TH1F* fWidthI1;
      TH1F* fNoiseI1;
      
      // Member variables from the fhicl file
      double                   fMinWidth;     ///< minimum initial width for ICARUS fit
      double                   fMaxWidthMult; ///< multiplier for max width for ICARUS fit
      
      int iWire;
      
      const geo::GeometryCore* fGeometry = lar::providerFrom<geo::Geometry>();
     
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
      fMaxMultiHit       = p.get< int             >("MaxMultiHit");
      fChi2NDF           = p.get< double          >("Chi2NDF");
    fLongPulseWidthVec = p.get< std::vector<int>>("LongPulseWidth", std::vector<int>() = {16,16,16});
          fLongMaxHitsVec    = p.get< std::vector<int>>("LongMaxHits",    std::vector<int>() = {25,25,25});
      fThetaAngle=p.get< double  >("ThetaAngle");
      fMinWidth=p.get< double  >("MinWidth");
      fMaxWidthMult=p.get< double  >("MaxWidthMult");

      
      fHitFinderTool  = art::make_tool<reco_tool::ICandidateHitFinder>(p.get<fhicl::ParameterSet>("CandidateHits"));
     // fPeakFitterTool = art::make_tool<reco_tool::PeakFitterICARUS>(p.get<fhicl::ParameterSet>("PeakFitter"));
      
    mf::LogInfo("ICARUSHitFinder_module") << "fDigitModuleLabel: " << fDigitModuleLabel << std::endl;
  }

  //-------------------------------------------------
  void ICARUSHitFinder::beginJob()
  {
      // get access to the TFile service
      art::ServiceHandle<art::TFileService> tfs;
      
      
      // ======================================
      // === Hit Information for Histograms ===
      fFirstChi2	= tfs->make<TH1F>("fFirstChi2", "#chi^{2}", 10000, 0, 1000);
      fNullChi2    = tfs->make<TH1F>("fNullChi2", "#chi^{2}", 100, 0, 10);

      fChi2	        = tfs->make<TH1F>("fChi2", "#chi^{2}", 10000, 0, 1000);
      fHeightC	= tfs->make<TH1F>("fHeightC", "height(ADC#)", 100, 0, 100);
      fWidthC	        = tfs->make<TH1F>("fWidthC", "width(samples)", 200, 0, 200);
      fNoiseC	        = tfs->make<TH1F>("fNoiseC", "Noise Area(ADC#)", 100, 0, 100);
      fAreaC            = tfs->make<TH1F>("fAreaC", "", 150, 0, 1500);
      fAreaC->GetXaxis()->SetTitle("Area(ADC#*ticks)");
      fAreaC->GetYaxis()->SetTitle("events");
      
      fWire2771            = tfs->make<TH1F>("fWire2771", "", 4096, -0.5, 4095.5);
      fWire2771->GetXaxis()->SetTitle("tick");
      fWire2771->GetYaxis()->SetTitle("ADC#");
      
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
   //   std::cout << " ICARUSHitfinder begin " << std::endl;
  }

  void ICARUSHitFinder::endJob()
  {
   //   std::cout << " ICARUSHitFinder endjob " << std::endl;
   
  }

  //-------------------------------------------------
  void ICARUSHitFinder::produce(art::Event& evt)
  {      //0
      //return;
      std::ofstream output("areaFit.out");

  //    std::cout << " ICARUSHitFinder produce " << std::endl;
      
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
      
      // ##########################################
      // ### Reading in the Wire List object(s) ###
      // ##########################################
      art::Handle< std::vector<recob::Wire> > wireVecHandle;
      evt.getByLabel(fCalDataModuleLabel,wireVecHandle);
      
      // #################################################################
      // ### Reading in the RawDigit associated with these wires, too  ###
      // #################################################################
      art::FindOneP<raw::RawDigit> RawDigits
      (wireVecHandle, evt, fCalDataModuleLabel);
      
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
      if(fThetaAngle==0) {minWireC=2535; maxWireC=4486;}
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
     
          double chi2null=ComputeNullChiSquare(holder);
//          std::cout << " wire " << iWire << " chi2null " << chi2null << std::endl;
          fNullChi2->Fill(chi2null);


    
      reco_tool::ICandidateHitFinder::HitCandidateVec      hitCandidateVec;

      reco_tool::ICandidateHitFinder::MergeHitCandidateVec mergedCandidateHitVec;
          
    std::vector<geo::WireID> wids = geom->ChannelToWire(channel);
          
          fHitFinderTool->findHitCandidates(holder, 0,channel,0,hitCandidateVec);
          for(auto& hitCand : hitCandidateVec) {
            expandHit(hitCand,holder,hitCandidateVec);
          }
          
          
          //hitCand->localmean=computeBestLocalMean(hitCand,holder,hitCandidateVec);
          fHitFinderTool->MergeHitCandidates(holder, hitCandidateVec, mergedCandidateHitVec);          
          
      //numHits = hits.size();
          int nghC=0;
          int nghI2=0;
          int nghI1=0;
          
      //  if(cryostat==0&&tpc==0&&plane==2)
        //     std::cout << "  Wire " << iwire << " numhits " << hitCandidateVec.size() <<" mergedhits " << mergedCandidateHitVec.size() << std::endl;
          for(auto& mergedCands : mergedCandidateHitVec)
      {
          int startT= mergedCands.front().startTick;
          int endT  = mergedCands.back().stopTick;
          // ### Putting in a protection in case things went wrong ###
          // ### In the end, this primarily catches the case where ###
          // ### a fake pulse is at the start of the ROI           ###
          if (endT - startT < 5) continue;
          fWidthC->Fill(endT-startT);
          // #######################################################
          // ### Clearing the parameter vector for the new pulse ###
          // #######################################################
          
          // === Setting the number of Gaussians to try ===
          int nGausForFit = mergedCands.size();
        
          // ##################################################
          // ### Calling the function for fitting ICARUS ###
          // ##################################################
          double                                chi2PerNDF(0.);
          double                                chi2Long(0.);

          int                                   NDF(1);
          const int npk=mergedCands.size();
          ICARUSPeakParamsVec peakParamsVec(npk);
          peakParamsVec.clear();
          int islong=0;
          if (mergedCands.size() <= fMaxMultiHit)
          {
              //std::cout << "setting iwire " << iwire << std::endl;
             // fPeakFitterTool->setWire(iwire);
            //  std::cout << " fitting iwire " << iwire << std::endl;
             // std::cout << " cryostat " << cryostat << " tpc " << tpc << " plane " << plane << " wire " << iwire << std::endl;
        findMultiPeakParameters(signal, mergedCands, peakParamsVec, chi2PerNDF, NDF, iwire);

          if (!(chi2PerNDF < std::numeric_limits<double>::infinity()))
          {
              chi2PerNDF = 200.;
              NDF        = 2;
          }
          fFirstChi2->Fill(chi2PerNDF);
            //  std::cout << " wire " << iwire << " first chi2NDF " << chi2PerNDF << std::endl;
          }

          //std::cout << " before longpulse chi2 " << chi2PerNDF << " threshold " << fChi2NDF << std::endl;
          if (chi2PerNDF < 10.)
              fChi2->Fill(chi2PerNDF);
       //   if(chi2PerNDF<0.1)      // change from 10 to reduce output
         //     std::cout << " wire " << iwire << " SMALL chi2NDF " << chi2PerNDF << " thr chi2NDF " << fChi2NDF << std::endl;
          ICARUSPeakParamsVec peakParamsLong(npk);
          peakParamsLong.clear();
          if (chi2PerNDF > fChi2NDF)
          {
              islong=1;
              findLongPeakParameters(signal, mergedCands, peakParamsLong, chi2Long, NDF, iwire);
          //    if(chi2Long<0.3) std::cout << " small chi2long " << chi2Long << std::endl;
              if(chi2Long<chi2PerNDF&&chi2Long>0.1) {
                  fChi2->Fill(chi2Long);
                  peakParamsVec=peakParamsLong;
              }
              else { fChi2->Fill(chi2PerNDF);
          //    if(chi2PerNDF<0.1)      // change from 10 to reduce output
           //       std::cout << " wire " << iwire << " SMALLSMALL chi2NDF " << chi2PerNDF << " thr chi2NDF " << fChi2NDF << std::endl;
          }
          }
          
          int jhit=0;
          
          
          for(const auto& peakParams : peakParamsVec)
          {
              // Extract values for this hit
              float peakAmp   = peakParams.peakAmplitude;
              float peakMean  = peakParams.peakCenter;
              //float peakWidth = peakParams.peakSigma;
              float peakLeft   = peakParams.peakTauLeft;
              float peakRight  = peakParams.peakTauRight;
              float peakBaseline = peakParams.peakBaseline;
              float peakSlope=0;
              float peakFitWidth=0;
              
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
              jhit++;

              TF1 Func("ICARUSfunc",fitf,start,end,1+5*mergedCands.size());
              for(unsigned int jf=0;jf<mergedCands.size();jf++) {
              Func.SetParameter(1+5*jf,peakBaseline);
              Func.SetParameter(2+5*jf,peakAmp);
              Func.SetParameter(3+5*jf,peakMean);
              Func.SetParameter(4+5*jf,peakRight);
              Func.SetParameter(5+5*jf,peakLeft);
              }
              TF1 FuncLong("ICARUSfuncLong",fitlong,start,end,1+7*mergedCands.size());
              for(unsigned int jf=0;jf<mergedCands.size();jf++) {
                  Func.SetParameter(1+7*jf,peakBaseline);
                  Func.SetParameter(2+7*jf,peakAmp);
                  Func.SetParameter(3+7*jf,peakMean);
                  Func.SetParameter(4+7*jf,peakRight);
                  Func.SetParameter(5+7*jf,peakLeft);
                  Func.SetParameter(6+7*jf,peakFitWidth);
                  Func.SetParameter(7+7*jf,peakSlope);
              }

              float fitCharge=0;
              if(!islong) {
              try
              { fitCharge=Func.Integral(start,end);}
              catch(...)
              {mf::LogWarning("ICARUSHitFinder") << "Icarus numerical integration failed";
                  fitCharge=std::accumulate(holder.begin() + (int) start, holder.begin() + (int) end, 0.);
              }
              }
              else {
                  try
                  { fitCharge=FuncLong.Integral(start,end);}
                  catch(...)
                  {mf::LogWarning("ICARUSHitFinder") << "Icarus numerical integration failed";
                      fitCharge=std::accumulate(holder.begin() + (int) start, holder.begin() + (int) end, 0.);
                  }
              }
              if(isnan(fitCharge)&&!islong) fitCharge=std::accumulate(holder.begin() + (int) start, holder.begin() + (int) end, 0.);
              if(isnan(fitCharge)&&islong) fitCharge=std::accumulate(holder.begin() + (int) start, holder.begin() + (int) end, 0.);
              //Func.Integral(start,end);
              float totSig=std::accumulate(holder.begin() + (int) start, holder.begin() + (int) end, 0.);
         //     float fitChargeErr=0;
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
            nGausForFit,                                                                             //MULTIPLICITY.
            jhit,                                                                            //LOCAL_INDEX.
            chi2PerNDF,                                                                 //WIRE ID.
            NDF                                                               //DEGREES OF FREEDOM.
            );
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
             // fWidthC->Fill(peakWidth);
             // fAreaC->Fill(totSig);
              wCharge[iWire]+=fitCharge;
              wInt[iWire]+=totSig;

          }

          if(plane==0) nhitsI1++;
          if(plane==1) nhitsI2++;
          if(plane==2) nhitsC++;
          
          if(plane==0&&tpc==0&&cryostat==0&&outDriftWindow) {nnhitsI1++; fNoiseI1->Fill(amplitude); }
          if(plane==1&&tpc==0&&cryostat==0&&outWireWindowI2) { nnhitsI2++; fNoiseI2->Fill(amplitude);
            //if(cryostat==0&&tpc==0)
              //    std::cout << " noise hit Wire " << hits[i].iWire << " tick " << hits[i].hitCenter << std::endl;
          }
          if(plane==2&&tpc==0&&cryostat==0&&outWireWindowC) {
        nnhitsC++; fNoiseC->Fill(amplitude); }
        
        ++hitIndex;
      } //end loop on found hits
        
          
          if(plane==2&&cryostat==0&&tpc==0)
           if(wCharge[iwire]>0)
           fAreaC->Fill(wCharge[iwire]);
          if(plane==2&&cryostat==0&&tpc==0)
              if(wInt[iwire]>0)
          fIntegralC->Fill(wInt[iwire]);
          if(plane==2&&cryostat==0&&tpc==0)
              if(wCharge[iwire]>0)
                  output << iwire << " " <<wInt[iwire] << std::endl;
          if(wInt[iwire]>0&&cryostat==0&&tpc==0)
              fAreaInt->Fill(wCharge[iwire]/wInt[iwire]);
          //if(wCharge[iwire]>0&&cryostat==0&&tpc==0&&wCharge[iwire]/wInt[iwire]<0.9)
            //  std::cout << " wire " << iwire << " ratio " << wCharge[iwire]/wInt[iwire] << std::endl;
      }
    } //end channel condition

    } //end loop on channels
      
//      std::cout <<  " nhitsI1 " << nhitsI1 <<" nhitsI2 " << nhitsI2 <<" nhitsC " << nhitsC << std::endl;


      for(unsigned int jw=minWireC;jw<maxWireC;jw++)
          fnhwC->Fill(nhWire[jw]);
      
      
    hcol.put_into(evt);
      //std::cout << " end ICARUSHitfinder " << std::endl;
     
      
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
                if(last==h2.startTick)
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
            //if(startTick-i-shift1>=0)
            samples1[i]=holder[h.startTick-i-shift1];
            //else
            //samples1[i]=outofbounds;
            
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
    void ICARUSHitFinder::findMultiPeakParameters(const std::vector<float>&                   roiSignalVec,
                                                   const reco_tool::ICandidateHitFinder::HitCandidateVec& hitCandidateVec,
                                                   ICARUSPeakParamsVec&                              peakParamsVec,
                                                   double&                                     chi2PerNDF,
                                                   int&                                        NDF, int iWire) const
    {
        ICARUSPeakParamsVec                              peakParamsVec0;
        TH1F* fHistogram=new TH1F("","",roiSignalVec.size(),0.,roiSignalVec.size());;
    
        std::string wireName = "PeakFitterHitSignal_" + std::to_string(iWire);
        fHistogram->SetName(wireName.c_str());
        
        if (hitCandidateVec.empty()) return;
        
        // in case of a fit failure, set the chi-square to infinity
        chi2PerNDF = std::numeric_limits<double>::infinity();
        
        int startTime = hitCandidateVec.front().startTick;
        int endTime   = hitCandidateVec.back().stopTick;
        int roiSize   = endTime - startTime;
        
        //std::cout << " roisize " << roiSize << std::endl;
        
        // Check to see if we need a bigger histogram for fitting
        if (roiSize > fHistogram->GetNbinsX())
        {
            std::string histName = "PeakFitterHitSignal_" + std::to_string(iWire);
            fHistogram = new TH1F(histName.c_str(),"",roiSize,0.,roiSize);
            //fHistogram->Sumw2();
        }
        
        fHistogram->Reset();
        for(int idx = 0; idx < roiSize; idx++)
            fHistogram->SetBinContent(idx+1,roiSignalVec.at(startTime+idx));
        for(int idx = 0; idx < roiSize; idx++)
            fHistogram->SetBinError(idx+1,2.4);
      //  for(int idx = 0; idx < roiSize; idx++)
        //   std::cout << " bin " << idx << " Error " << fHistogram->GetBinError(idx+1) << std::endl;
        
        
        // Build the string to describe the fit formula
        std::string equation = "gaus(0)";
        
        
        // Now define the complete function to fit
        // Now define the complete function to fit
        TF1 Func("ICARUSfunc",fitf,0,roiSize,1+5*hitCandidateVec.size());
        
        
        // ### Setting the parameters for the ICARUS Fit ###
        int parIdx(0);
        
        for(auto& candidateHit : hitCandidateVec)
        {
            double peakMean   = candidateHit.hitCenter - float(startTime);
            double peakWidth  = candidateHit.hitSigma;
            // std::cout << " peakWidth " << peakWidth << std::endl;
            // std::cout << " hitcenter " << candidateHit.hitCenter << " starttime " << startTime <<std::endl;
            
            
            double amplitude  = candidateHit.hitHeight;
            // double meanLowLim = std::max(peakMean - fPeakRange * peakWidth,              0.);
            // double meanHiLim  = std::min(peakMean + fPeakRange * peakWidth, double(roiSize));

            Func.SetParameter(0,hitCandidateVec.size());
            Func.SetParLimits(0,hitCandidateVec.size(),hitCandidateVec.size());
            
            Func.SetParameter(1+parIdx,0);
            Func.SetParameter(2+parIdx, amplitude);
            Func.SetParameter(3+parIdx, peakMean);
            Func.SetParameter(4+parIdx,peakWidth);
            Func.SetParameter(5+parIdx,peakWidth);
            
            Func.SetParLimits(1+parIdx, -5, 5);
            Func.SetParLimits(2+parIdx, 0.1 * amplitude,  10. * amplitude);
            Func.SetParLimits(3+parIdx, peakMean-peakWidth,peakMean+peakWidth);
            Func.SetParLimits(4+parIdx, std::max(fMinWidth, 0.01 * peakWidth), fMaxWidthMult * peakWidth);
            Func.SetParLimits(5+parIdx, std::max(fMinWidth, 0.01 * peakWidth), fMaxWidthMult * peakWidth);
        
            parIdx += 5;
            
        }
        
        int fitResult(-1);
        // if(hitCandidateVec.size()>2) return;
        try
        {  fitResult = fHistogram->Fit(&Func,"QRWB","", 0., roiSize);
        }
        catch(...)
        {mf::LogWarning("GausHitFinder") << "Fitter failed finding a hit";}
        
        
       // if(fitResult==0)
       //     std::cout << " icarus fit converges " << iWire << std::endl;
        if(fitResult!=0)
            std::cout << " icarus fit cannot converge " << iWire << std::endl;
        // ##################################################
        // ### Getting the fitted parameters from the fit ###
        // ##################################################
        NDF        = roiSize-5*hitCandidateVec.size();
        chi2PerNDF = (Func.GetChisquare() / NDF);
        
        double chi2mio=ComputeChiSquare(Func,fHistogram);
//        std::cout << " chi2mio " << chi2mio << std::endl;
        chi2PerNDF=chi2mio;
        //  for(int idx = 0; idx < roiSize; idx++)
        //   std::cout << " bin " << idx << " Error " << fHistogram->GetBinError(idx+1) << std::endl;
        
//         std::cout << " chi2 " << Func.GetChisquare() << std::endl;
//         std::cout << " ndf " << NDF << std::endl;
        
        //      std::cout << " chi2ndf " << chi2PerNDF<< std::endl;
        parIdx = 0;
        peakParamsVec0.clear();

        for(size_t idx = 0; idx < hitCandidateVec.size(); idx++)
        {
            ICARUSPeakFitParams_t peakParams;
            
            peakParams.peakAmplitude      = Func.GetParameter(2+parIdx);
            peakParams.peakAmplitudeError = Func.GetParError(2+parIdx);
            peakParams.peakCenter         = Func.GetParameter(3+parIdx) + float(startTime);
            peakParams.peakCenterError    = Func.GetParError(3+parIdx);
            //std::cout << " rising time " << Func.GetParameter(3) << " falling time " <<Func.GetParameter(4) << std::endl;
            peakParams.peakTauRight        = Func.GetParameter(4+parIdx);
            peakParams.peakTauRightError        = Func.GetParError(4+parIdx);
            peakParams.peakTauLeft        = Func.GetParameter(5+parIdx);
            peakParams.peakTauLeftError        = Func.GetParError(5+parIdx);
            peakParams.peakBaseline        = Func.GetParameter(1+parIdx);
            peakParams.peakBaselineError        = Func.GetParError(1+parIdx);
            peakParams.peakFitWidth        =0;
            peakParams.peakFitWidthError        = 0;
            peakParams.peakSlope        = 0;
            peakParams.peakSlopeError        = 0;
            peakParamsVec.emplace_back(peakParams);
            parIdx += 5;
            //std::cout << " first center " << Func.GetParameter(2) + float(startTime) << std::endl;
            // std::cout << " second center " << Func.GetParameter(7) + float(startTime) << std::endl;
        }
        
        //Gaus.Delete();
        //Func.Delete();
        bool writeWaveform=false;
        if(writeWaveform) {
       TFile *f = new TFile("fitICARUS.root","UPDATE");
    fHistogram->Write();
        f->Close();
        f->Delete();
        }

        fHistogram->Delete();
        return;
    }

    void ICARUSHitFinder::findLongPeakParameters(const std::vector<float>&                   roiSignalVec,
                                                  const reco_tool::ICandidateHitFinder::HitCandidateVec& hitCandidateVec,
                                                  ICARUSPeakParamsVec&                              peakParamsVec,
                                                  double&                                     chi2PerNDF,
                                                  int&                                        NDF, int iWire) const
    {
        TH1F* fHistogram=new TH1F("","",roiSignalVec.size(),0.,roiSignalVec.size());;
        std::string wireName = "PeakFitterHitSignal_" + std::to_string(iWire);
        fHistogram->SetName(wireName.c_str());
        
        if (hitCandidateVec.empty()) return;
        
        // in case of a fit failure, set the chi-square to infinity
        chi2PerNDF = std::numeric_limits<double>::infinity();
        
        int startTime = hitCandidateVec.front().startTick;
        int endTime   = hitCandidateVec.back().stopTick;
        int roiSize   = endTime - startTime;
        
        //std::cout << " roisize " << roiSize << std::endl;
        
        // Check to see if we need a bigger histogram for fitting
        if (roiSize > fHistogram->GetNbinsX())
        {
            std::string histName = "PeakFitterHitSignal_" + std::to_string(iWire);
            fHistogram = new TH1F(histName.c_str(),"",roiSize,0.,roiSize);
            fHistogram->Sumw2();
        }
        
        fHistogram->Reset();
        for(int idx = 0; idx < roiSize; idx++)
            fHistogram->SetBinContent(idx+1,roiSignalVec.at(startTime+idx));
        
        // Build the string to describe the fit formula
        std::string equation = "gaus(0)";
        
        
        
        // Now define the complete function to fit
        TF1 Func("ICARUSfunc",fitlong,0,roiSize,1+7*hitCandidateVec.size());
        
        // ### Setting the parameters for the ICARUS Fit ###
        int parIdx(0);
        for(auto& candidateHit : hitCandidateVec)
        {
            double peakMean   = candidateHit.hitCenter - float(startTime);
            double peakWidth  = candidateHit.hitSigma;
            double amplitude  = candidateHit.hitHeight;
            
            Func.SetParameter(0,hitCandidateVec.size());
            Func.SetParLimits(0,hitCandidateVec.size(),hitCandidateVec.size());
            
            Func.SetParameter(1+parIdx,0);
            Func.SetParameter(2+parIdx, amplitude);
            Func.SetParameter(3+parIdx, peakMean);
            Func.SetParameter(4+parIdx,peakWidth);
            Func.SetParameter(5+parIdx,peakWidth);
            Func.SetParameter(6+parIdx,2*peakWidth);
            Func.SetParameter(7+parIdx,0);
            
            
            Func.SetParLimits(1+parIdx, -5, 5);
            Func.SetParLimits(2+parIdx, 0.1 * amplitude,  10. * amplitude);
            Func.SetParLimits(3+parIdx, peakMean-peakWidth,peakMean+peakWidth);
            Func.SetParLimits(4+parIdx, std::max(fMinWidth, 0.01 * peakWidth), fMaxWidthMult * peakWidth);
            Func.SetParLimits(5+parIdx, std::max(fMinWidth, 0.01 * peakWidth), 4 * peakWidth);
            Func.SetParLimits(6+parIdx, 0,4*peakWidth);
            Func.SetParLimits(7+parIdx, -1,1);
            
            parIdx += 7;
            
        }
        int fitResult(-1);
        try
        {  fitResult = fHistogram->Fit(&Func,"QRWB","", 0., roiSize);
        }
        catch(...)
        {mf::LogWarning("GausHitFinder") << "Fitter failed finding a hit";}
        
//        if(fitResult!=0)
//            std::cout << " long fit cannot converge " << iWire << std::endl;
        // ##################################################
        // ### Getting the fitted parameters from the fit ###
        // ##################################################
        NDF        = roiSize-7*hitCandidateVec.size();
        chi2PerNDF = (Func.GetChisquare() / NDF);
        
        parIdx = 0;
        peakParamsVec.clear();
        for(size_t idx = 0; idx < hitCandidateVec.size(); idx++)
        {
            ICARUSPeakFitParams_t peakParams;
            
            peakParams.peakAmplitude      = Func.GetParameter(2+parIdx);
            peakParams.peakAmplitudeError = Func.GetParError(2+parIdx);
            peakParams.peakCenter         = Func.GetParameter(3+parIdx) + float(startTime);
            peakParams.peakCenterError    = Func.GetParError(3+parIdx);
            
            peakParams.peakTauRight        = Func.GetParameter(4+parIdx);
            peakParams.peakTauRightError        = Func.GetParError(4+parIdx);
            peakParams.peakTauLeft        = Func.GetParameter(5+parIdx);
            peakParams.peakTauLeftError        = Func.GetParError(5+parIdx);
            peakParams.peakFitWidth        = Func.GetParameter(6+parIdx);
            peakParams.peakFitWidthError        = Func.GetParError(6+parIdx);
            peakParams.peakSlope        = Func.GetParameter(7+parIdx);
            peakParams.peakSlopeError        = Func.GetParError(7+parIdx);
            peakParams.peakBaseline        = Func.GetParameter(1+parIdx);
            peakParams.peakBaselineError        = Func.GetParError(1+parIdx);
         //   std::cout << " before adding peakparams size  " << peakParamsVec.size() << std::endl;
            peakParamsVec.emplace_back(peakParams);
           // std::cout << " after adding peakparams size  " << peakParamsVec.size() << std::endl;
            
            parIdx += 7;
            
        }
        
        //Func.Delete();
        bool writeWaveform=false;
        if(writeWaveform) {
            TFile *f = new TFile("fitICARUS.root","UPDATE");
            fHistogram->Write();
            f->Close();
            f->Delete();
        }
        
        fHistogram->Delete();
        return;
    }
    

Double_t ICARUSHitFinder::fitf(Double_t *x, Double_t *par)
        {
            int npeaks=(int)(par[0]);
            Double_t fitval=0;
            for(int jp=0;jp<npeaks;jp++)
                fitval += par[5*jp+1]+par[5*jp+2]*TMath::Exp(-(x[0]-par[5*jp+3])/par[5*jp+4])/(1+TMath::Exp(-(x[0]-par[5*jp+3])/par[5*jp+5]));
            return fitval;
        }
Double_t ICARUSHitFinder::fitlong(Double_t *x, Double_t *par)
    {
        int npeaks=(int)(par[0]);
        Double_t fitval=0;
        for(int jp=0;jp<npeaks;jp++) {
            for (int js=0; js<floor(par[7*jp+6]); js++) {
                fitval += (1.+js*par[7*jp+7])*(par[7*jp+1]+par[7*jp+2]*TMath::Exp(-(x[0]-par[7*jp+3])/par[7*jp+4])/(1+TMath::Exp(-(x[0]-par[7*jp+3])/par[7*jp+5])))/(par[7*jp+6]);
            }
        }
        return fitval;
    }
    
    double ICARUSHitFinder::ComputeChiSquare(TF1 func, TH1 *histo) const
    {
        double chi=0;
        int nb=histo->GetNbinsX();
        double wb=histo->GetBinWidth(0);
        
        int jp;
        for( jp=1;jp<nb;jp++) {
            if(histo->GetBinContent(jp)==0) break;
            double xb=jp*wb;
            double fv=(func)(xb);
            double hv=histo->GetBinContent(jp);
            double dv=hv-fv;
            double cv=dv/(histo->GetBinError(jp));
            chi+=cv*cv;
            //std::cout << " chi " << chi << std::endl;
            
        }
//        std::cout << " chi2mio" << chi << std::endl;
        //std::cout << " ndf " << ndf << std::endl;
        return chi/(jp-5);
    }
    double ICARUSHitFinder::ComputeNullChiSquare(std::vector<float> holder) const
    {
        double chi=0;
        int nb=33;
        
        int jp;
        for( jp=1;jp<nb;jp++) {
            
            double hv=holder[jp];
            double dv=hv;
            double cv=dv/2.4;
            chi+=cv*cv;
           // std::cout << " chi " << chi << std::endl;
            
        }
       // std::cout << " chi2null" << chi << std::endl;
        //std::cout << " ndf " << ndf << std::endl;
        return chi/(jp);
    }
    
  
  DEFINE_ART_MODULE(ICARUSHitFinder)

} // end namespace hit


#endif //ICARUSHitFinder_H
