
////////////////////////////////////////////////////////////////////////
//
// ICARUSPurityDQM class
//
// Christian Farnese
//
////////////////////////////////////////////////////////////////////////

//#ifndef ICARUSPURITYDQM_H
//#define ICARUSPURITYDQM_H

#include <iomanip>
#include <TH1F.h>
#include <TProfile.h>
#include <vector>
#include <string>
#include <array>
#include <fstream>

//Framework includes
#include "art/Framework/Core/ModuleMacros.h" 
#include "canvas/Persistency/Common/FindManyP.h"
#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
//#include "art/Framework/Services/Registry/ServiceHandle.h" 
//#include "art/Framework/Services/Optional/TFileService.h" 
#include "art_root_io/TFileService.h"
//#include "art/Framework/Services/Optional/TFileDirectory.h" 
//#include "messagefacility/MessageLogger/MessageLogger.h" 

//LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nug4/ParticleNavigation/ParticleList.h"
#include "nug4/ParticleNavigation/EmEveIdCalculator.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "larsim/MCCheater/BackTracker.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"

//Database Connection Files
//#include "../../MetricManagerShim/MetricManager.hh"
//#include "../../MetricConfig/ConfigureRedis.hh"

//purity info class
#include "sbnobj/Common/Analysis/TPCPurityInfo.hh"

#include "art/Framework/Core/EDProducer.h"
#include <TMath.h>
#include <TH1F.h>
#include "TH2D.h"
#include "TProfile2D.h"
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include "TF1.h"
#include "TCanvas.h"
#include "TNtuple.h"
#include "TFitResultPtr.h"

//Redis connection 
//#include "sbndaq-redis-plugin/Utilities.h"

class TH1F;
class TH2F;


///Cluster finding and building 
namespace icarus {
   
  class ICARUSPurityDQM : public art::EDProducer {
    
  public:
    
    explicit ICARUSPurityDQM(fhicl::ParameterSet const& pset);
    virtual ~ICARUSPurityDQM();
    
    /// read access to event
    void produce(art::Event& evt);
    void beginJob();
    void endJob();
    int Nothere(std::vector<int>* a, int b);
    int Notheref(std::vector<float>* a, float b);
    Double_t FoundMeanLog(std::vector<float>* a,float b);

  private:

    TH1F* puritytpc0;
    TH1F* puritytpc1;
    TH1F* puritytpc2;
    TH1F* puritytpc3;
    



    //TH1F* purityvalues;
    TH1F* purityvalues2;
    TH1F* h_hit_height;
       TH1F* h_hit_area;
     TH2D* h_hit_height_area;

    TH1F* h_basediff;
    TH1F* h_base1;
    TH1F* h_base2;
    TH1F* h_basebase;
       TH1F* h_ratio;
       TH1F* h_ratio_after;
        TH1F* h_ratio_after_2;
    //TH2D* h_ratio_after_3;
       TH2D* h_ratio_3;


    TH1F* h_rms;
    TH1F* h_errors;
    TH1F* h_hittime;
    TH1F* h_hittime_2;
    TH1F* h_hittime_3;

    //TH1F* fRun; 
    //TH2D* fRunSub;

    std::vector<art::InputTag>  fDigitModuleLabel;
    //short fPrintLevel;

    float fValoretaufcl; 
    float fmaxfcl;
    float fminfcl;
    float fdisfcl;
    int fcryofcl;
    int fquantevoltefcl;

    int fplanefcl;
    int fthresholdfcl;
    int fdumphitsfcl;
    int fgruppifcl;
int fdwclusfcl;
int fdsclusfcl;
    bool fPersistPurityInfo;
    
    TNtuple* purityTuple;
    bool fFillAnaTuple;

      TTree* fpurTree;
      int run_tree;
      //int subrun_tree;
      int event_tree;
      int tpc_tree;
      int qhits_tree;
      float dwire_tree;
      float awire_tree;
      float dtime_tree;
      float earea_tree;
      float marea_tree;
      float slope_tree;
      float errslope_tree;
      float chi_tree;
      long time_tree;

      
  }; // class ICARUSPurityDQM
  
}

//#endif 

namespace icarus{

  //--------------------------------------------------------------------
  ICARUSPurityDQM::ICARUSPurityDQM(fhicl::ParameterSet const& pset)
    : EDProducer(pset)
    , fDigitModuleLabel     (pset.get< std::vector<art::InputTag> > ("RawModuleLabel"))
    , fValoretaufcl         (pset.get< float >       ("ValoreTauFCL"))
    , fmaxfcl               (pset.get< float >       ("FracMaxFCL",0.9))
    , fminfcl               (pset.get< float >       ("FraxMinFCL",0.05))
    , fdisfcl               (pset.get< float >       ("DisFCL",3.0))
    , fcryofcl              (pset.get< int >         ("CryostatFCL"))
    , fquantevoltefcl       (pset.get< int >         ("QuanteFCL",5))
    , fplanefcl             (pset.get< int >         ("PlaneFCL"))
    , fthresholdfcl         (pset.get< int >         ("ThresholdFCL"))
    , fdumphitsfcl          (pset.get< int >         ("DumpHitFCL",0))
    , fgruppifcl            (pset.get< int >         ("GruppiFCL",8))
    , fdwclusfcl            (pset.get< int >         ("DeltaWClusFCL",3))
    , fdsclusfcl            (pset.get< int >         ("DeltaSClusFCL",100))
    , fPersistPurityInfo    (pset.get< bool  >       ("PersistPurityInfo",false))
    , fFillAnaTuple         (pset.get< bool  >       ("FillAnaTuple",false))
  {

    //declare what we produce .. allow it to not be persistable to the event
    if(fPersistPurityInfo)
      produces< std::vector<anab::TPCPurityInfo> >("",art::Persistable::Yes);
    else
      produces< std::vector<anab::TPCPurityInfo> >("",art::Persistable::No);      
 
  }
  
  //------------------------------------------------------------------
  ICARUSPurityDQM::~ICARUSPurityDQM()
  {
  
  }
  
  void ICARUSPurityDQM::beginJob()
  {
    
    // get access to the TFile service
    art::ServiceHandle<art::TFileService> tfs;
  
    //purityvalues = tfs->make<TH1F>("purityvalues","purityvalues",20000,-10,10);
    h_basediff = tfs->make<TH1F>("h_basediff","h_basediff",10000,-20,20);
    h_basebase = tfs->make<TH1F>("h_basebase","h_basebase",10000,-20,20);
    h_base1 = tfs->make<TH1F>("h_base1","h_base1",10000,-20,20);
    h_base2 = tfs->make<TH1F>("h_base2","h_base2",10000,-20,20);
    
    h_rms = tfs->make<TH1F>("h_rms","h_rms",2000,0,20);
    //fRun=tfs->make<TH1F>("fRun","Events per run", 4000,0.5 ,4000.5);
    //fRunSub=tfs->make<TH2D>("fRunSub","Events per run", 4000,0.5 ,4000.5,50,0.5,50.5);
    h_errors = tfs->make<TH1F>("h_errors","",1000,0,1000);


   h_hittime = tfs->make<TH1F>("h_hittime","",2500,-0.5,2499.5);
   h_hittime_2 = tfs->make<TH1F>("h_hittime_2","",2500,-0.5,2499.5);
   h_hittime_3 = tfs->make<TH1F>("h_hittime_3","",2500,-0.5,2499.5);

    purityvalues2 = tfs->make<TH1F>("purityvalues2","purityvalues2",20000,-10,10);
    puritytpc0 = tfs->make<TH1F>("puritytpc0","puritytpc0",20000,-10,10);
    puritytpc1 = tfs->make<TH1F>("puritytpc1","puritytpc1",20000,-10,10);
    puritytpc2 = tfs->make<TH1F>("puritytpc2","puritytpc2",20000,-10,10);
    puritytpc3 = tfs->make<TH1F>("puritytpc3","puritytpc3",20000,-10,10);
    h_hit_height = tfs->make<TH1F>("h_hit_height","h_hit_height",2000,0.5,2000.5);
    h_hit_area = tfs->make<TH1F>("h_hit_area","h_hit_area",10000,0,10000);
    h_hit_height_area=tfs->make<TH2D>("h_hit_height_area","", 300,0.5,300.5,2000,0.,2000.);

    h_ratio = tfs->make<TH1F>("h_ratio","h_ratio",20000,0,5);
    h_ratio_after = tfs->make<TH1F>("h_ratio_after","h_ratio_after",20000,0,5);
    h_ratio_after_2 = tfs->make<TH1F>("h_ratio_after_2","",2000,0,2);
    //h_ratio_after_3 = tfs->make<TH2D>("h_ratio_after_3","",1300,0,2600,2000,0,2);
    h_ratio_3 = tfs->make<TH2D>("h_ratio_3","",300,2000,2600,2000,0,2);



      fpurTree = tfs->make<TTree>("puritytree","tree for the purity studies");
      fpurTree->Branch("run",&run_tree,"run/I");
      //fpurTree->Branch("subrun",&subrun_tree,"subrun/I");
      fpurTree->Branch("event",&event_tree,"event/I");
      fpurTree->Branch("time",&time_tree,"time/L");
      fpurTree->Branch("tpc",&tpc_tree,"tpc/I");
      fpurTree->Branch("qhits",&qhits_tree,"qhits/I");
      fpurTree->Branch("deltawire",&dwire_tree,"deltawire/F");
      fpurTree->Branch("averagewire",&awire_tree,"averagewire/F");
      fpurTree->Branch("deltatime",&dtime_tree,"deltatime/F");
      fpurTree->Branch("errarea",&earea_tree,"errarea/F");
      fpurTree->Branch("meanarea",&marea_tree,"meanarea/F");
      fpurTree->Branch("slope",&slope_tree,"slope/F");
      fpurTree->Branch("errorslope",&errslope_tree,"errorslope/F");
      fpurTree->Branch("chisquare",&chi_tree,"chisquare/F");

      
    if(fFillAnaTuple)
      purityTuple = tfs->make<TNtuple>("purityTuple","Purity Tuple","run:ev:tpc:att");

  }
  
  void ICARUSPurityDQM::endJob()
  {
    std::ofstream goodpurofinal("valore_indicativo.out",std::ios::app);
    goodpurofinal << purityvalues2->GetMean() << std::endl;
    //if(fPrintLevel == -1) outFile.close();
  }
  
  
  int ICARUSPurityDQM::Nothere(std::vector<int>* a, int b){
    int Cisono=3;
    for(int i=0; i<(int)a->size();i++)
      {
	if((*a)[i]==b)
	  {
	    Cisono=8;
	  }
      }
    return Cisono;//se "ci sono" vale 8 allora l'intero b e' contenuto nel vettore a mentre se "ci sono" vale 3 a non contiene b.
  }
  
  
  int ICARUSPurityDQM::Notheref(std::vector<float>* a, float b){
    int Cisono=3;
    for(int i=0; i<(int)a->size();i++)
      {
	if((*a)[i]==b)
	  {
	    Cisono=8;
	  }
      }
        
    return Cisono;//se "ci sono" vale 8 allora l'intero b e' contenuto nel vettore a mentre se "ci sono" vale 3 a non contiene b.
  }
    
  Double_t ICARUSPurityDQM::FoundMeanLog(std::vector<float>* a,float b){


  float taglio=0;
  if(a->size()>0){
  int punto_taglio=a->size()*(b)+0.5;
  //cout << punto_taglio << " e " << a->size() << endl;
  if(punto_taglio==0)punto_taglio=1;
std::vector<float>* usedhere=new std::vector<float>;
 for(int j=0;j<(int)a->size();j++)usedhere->push_back((*a)[j]);

  //cout << (*usedhere)[punto_taglio-1] << " PRIMA " << endl;
std::sort(usedhere->begin(),usedhere->end(), [](float &c, float &d){ return c<d; });
  //cout << (*usedhere)[punto_taglio-1] << " DOPO " << endl;
  taglio=(*usedhere)[punto_taglio-1];
  delete usedhere;
  }
  return taglio;

/*        
    int punto_taglio=a->size()*(1-b)+0.5;
    ///////std::cout << punto_taglio << " e " << a->size() << std::endl;
    //if(punto_taglio==0)punto_taglio=1;
    std::vector<float>* usedhere=new std::vector<float>;
    for(int jj=0;jj<punto_taglio+1;jj++)
      {
	//cout << jj << endl;
	float maximum=0;
	for(int j=0;j<(int)a->size();j++)
	  {
	    if((*a)[j]>maximum && Notheref(usedhere,(*a)[j])==3)
	      {
		maximum=(*a)[j];
	      }
	  }
	//cout << maximum << endl;
	usedhere->push_back(maximum);
      }
    //cout << (*usedhere)[punto_taglio] << endl;
    return (*usedhere)[punto_taglio-1];
*/
  }
      
  void ICARUSPurityDQM::produce(art::Event& evt)
  {
    
//      //std::cout << " Inizia Purity ICARUS Ana - upgraded by C.FARNESE, WES and OLIVIA " << std::endl;
      // code stolen from TrackAna_module.cc
      art::ServiceHandle<geo::Geometry>      geom;
      unsigned int  fDataSize;
      std::vector<short> rawadc;      //UNCOMPRESSED ADC VALUES.
      // get all hits in the event
      //InputTag cluster_tag { "fuzzycluster" }; //CH comment trovato con eventdump code
      
      //to get run and event info, you use this "eventAuxillary()" object.
      //art::Timestamp ts = evt.time();
      //std::cout << "Processing for Purity " << " Run " << evt.run() << ", " << "Event " << evt.event() << " and Time " << evt.time().timeHigh() << std::endl;
      
      //fRun->Fill(evt.run());
      //fRunSub->Fill(evt.run(),evt.subRun());
      art::Handle< std::vector<raw::RawDigit> > digitVecHandle;
      std::vector<const raw::RawDigit*> rawDigitVec;

      //setup output vector
      std::unique_ptr< std::vector<anab::TPCPurityInfo> > outputPtrVector(new std::vector<anab::TPCPurityInfo>() );
    
    anab::TPCPurityInfo purity_info;
    purity_info.Run = evt.run();
    purity_info.Subrun = evt.subRun(); //evt.time().timeHigh()-1598580000;//evt.subRun();
    purity_info.Event = evt.event();
    std::ofstream purh("dump_purity_hits.out",std::ios::app);

    for(const auto& digitlabel : fDigitModuleLabel)
      {
	evt.getByLabel(digitlabel, digitVecHandle);
	std::vector<const raw::RawDigit*> rawDigitVec;
	
	if (digitVecHandle.isValid())
	  {
	    // Sadly, the RawDigits come to us in an unsorted condition which is not optimal for
	    // what we want to do here. So we make a vector of pointers to the input raw digits and sort them
	    // Ugliness to fill the pointer vector...
	    for(size_t idx = 0; idx < digitVecHandle->size(); idx++) 
	      rawDigitVec.push_back(&digitVecHandle->at(idx));
	    // Sort (use a lambda to sort by channel id)
	    std::sort(rawDigitVec.begin(),
		      rawDigitVec.end(),
		      [](const raw::RawDigit* left, const raw::RawDigit* right) {return left->Channel() < right->Channel();});
	  }
	
	
	std::vector<int> *www0=new std::vector<int>;
	std::vector<float> *sss0=new std::vector<float>;
	std::vector<float> *hhh0=new std::vector<float>;
	std::vector<float> *ehh0=new std::vector<float>;
	std::vector<int> *www2=new std::vector<int>;
	std::vector<float> *sss2=new std::vector<float>;
	std::vector<float> *hhh2=new std::vector<float>;
	std::vector<float> *ehh2=new std::vector<float>;
	std::vector<int> *ccc0=new std::vector<int>;
	std::vector<int> *ccc2=new std::vector<int>;
	std::vector<float> *aaa0=new std::vector<float>;
	std::vector<float> *aaa2=new std::vector<float>;
	
	
        short int used[4096];	
	for(const auto& rawDigit : rawDigitVec)
	  {
	    raw::ChannelID_t channel = rawDigit->Channel();
	    ////std::cout << channel << std::endl;
	    std::vector<geo::WireID> wids = geom->ChannelToWire(channel);
	    // for now, just take the first option returned from ChannelToWire
	    geo::WireID wid  = wids[0];
	    // We need to know the plane to look up parameters
	    geo::PlaneID::PlaneID_t plane = wid.Plane;
	    //size_t 
	    
	    int cryostat=wid.Cryostat;
	    size_t tpc=wid.TPC;
	    size_t iWire=wid.Wire;
	    ////std::cout << channel << " INFO CHANNEL " << plane << " " << tpc << " " << cryostat << " " << iWire << std::endl;
	    fDataSize = rawDigit->Samples();
	    rawadc.resize(fDataSize);
	    
	    for(int ijk=0;ijk<4096;ijk++)used[ijk]=0;
	    //UNCOMPRESS THE DATA.
	    int pedestal = (int)rawDigit->GetPedestal();
	    float pedestal2 = rawDigit->GetPedestal();
	    float sigma_pedestal = rawDigit->GetSigma();
	    raw::Uncompress(rawDigit->ADCs(), rawadc, pedestal, rawDigit->Compression());
	    ////std::cout << pedestal2 << " " << fDataSize << " " << rawadc.size() << " " << sigma_pedestal << std::endl;
	    float massimo=0;
	    float quale_sample_massimo;
	    //    if (plane==0) {
	    if ((int)plane==fplanefcl && cryostat==fcryofcl){
	      for(int volte=0;volte<fquantevoltefcl;volte++){	
            massimo=0;
           quale_sample_massimo=-1;

              TH1F *h111 = new TH1F("h111","delta aree",20000,0,0);
              for (unsigned int ijk=0; ijk<(fDataSize); ijk++)
                {
		  if ((rawDigit->ADC(ijk)-pedestal2)>massimo && ijk>150 && ijk<(fDataSize-150) && used[ijk]==0)
		    {
		      massimo=(rawDigit->ADC(ijk)-pedestal2);
		      quale_sample_massimo=ijk;
		    }
                }
              float base_massimo_before=0;
              float base_massimo_after=0;
              for (unsigned int ijk=quale_sample_massimo-135; ijk<quale_sample_massimo-35; ijk++)
		{
                  base_massimo_before+=(rawDigit->ADC(ijk)-pedestal2)*0.01;
		}
              for (unsigned int ijk=quale_sample_massimo+35; ijk<quale_sample_massimo+135; ijk++)
		{
                  base_massimo_after+=(rawDigit->ADC(ijk)-pedestal2)*0.01;
		}
              float basebase=(base_massimo_after+base_massimo_before)*0.5;
              h_basediff->Fill(fabs(base_massimo_after-base_massimo_before)); h_base1->Fill(base_massimo_before);
	      h_base2->Fill(base_massimo_after); h_basebase->Fill(basebase);
              float areaarea=0;
              for (unsigned int ijk=0; ijk<(fDataSize); ijk++)
		{
                  if (ijk>(quale_sample_massimo-30) && ijk<(quale_sample_massimo+30)) {
		    if(fabs(base_massimo_after-base_massimo_before)<=sigma_pedestal)areaarea+=(rawDigit->ADC(ijk)-pedestal2-basebase);
                    if(fabs(base_massimo_after-base_massimo_before)>sigma_pedestal && base_massimo_after<base_massimo_before)areaarea+=(rawDigit->ADC(ijk)-pedestal2-base_massimo_after);
                    if(fabs(base_massimo_after-base_massimo_before)>sigma_pedestal && base_massimo_after>base_massimo_before)areaarea+=(rawDigit->ADC(ijk)-pedestal2-base_massimo_before);

                  }
                  else{
		    if(used[ijk]==0)h111->Fill(rawDigit->ADC(ijk)-pedestal2-basebase);
                  }                  
		}
              sigma_pedestal=h111->GetRMS();
              h111->Delete();
              h_rms->Fill(sigma_pedestal);
              if(fdumphitsfcl==1)purh<<evt.event()<< " " <<iWire<<" "<<tpc<<" "<<cryostat<<" "<<basebase<<" "<<base_massimo_after<<" " <<base_massimo_before<<" "<<pedestal2<<" "<<sigma_pedestal<< " "  << quale_sample_massimo << " " << massimo << std::endl;

              if (massimo>(fthresholdfcl*sigma_pedestal)/* && fabs(base_massimo_after-base_massimo_before)<sigma_pedestal*/)
                {
                     if(fdumphitsfcl==1)purh<<iWire<<" "<<tpc<<" "<<cryostat<< " CANDIDATE HIT "  << quale_sample_massimo << " " << massimo << std::endl;

                     if(tpc==0)www0->push_back(iWire+64);            
                     if(tpc==0)sss0->push_back(quale_sample_massimo);
                     if(tpc==0)hhh0->push_back(massimo);
                     if(tpc==0)ehh0->push_back(sigma_pedestal);
                     if(tpc==0)ccc0->push_back(-1);
                     if(tpc==1)www0->push_back(iWire+64+2536);
                     if(tpc==1)sss0->push_back(quale_sample_massimo);
                     if(tpc==1)hhh0->push_back(massimo);
                     if(tpc==1)ehh0->push_back(sigma_pedestal);
                     if(tpc==1)ccc0->push_back(-1);
                     if(tpc==2)www2->push_back(iWire+64);
                     if(tpc==2)sss2->push_back(quale_sample_massimo);
                     if(tpc==2)hhh2->push_back(massimo);
                     if(tpc==2)ehh2->push_back(sigma_pedestal);
                     if(tpc==2)ccc2->push_back(-1);
                     if(tpc==3)www2->push_back(iWire+64+2536);
                     if(tpc==3)sss2->push_back(quale_sample_massimo);
                     if(tpc==3)hhh2->push_back(massimo);
                     if(tpc==3)ehh2->push_back(sigma_pedestal);
                     if(tpc==3)ccc2->push_back(-1);

		     if(tpc==0)aaa0->push_back(areaarea);
		     if(tpc==1)aaa0->push_back(areaarea);
		     if(tpc==2)aaa2->push_back(areaarea);
		     if(tpc==3)aaa2->push_back(areaarea);
                     for(int ijk=0;ijk<330;ijk++)
          {
	int ent_value=quale_sample_massimo+ijk-165;
	used[ent_value]=1;
                }

                }
                }
            }
	  }
////std::cout << "SIZE " << www0->size() << " " << www1->size() << " " << www2->size() << " " << www3->size() << std::endl;	
	for (unsigned int ijk=0; ijk<www0->size(); ijk++)
	  {
	    for (unsigned int ijk2=0; ijk2<www0->size(); ijk2++)
	      {
		if (fabs((*www0)[ijk]-(*www0)[ijk2])<fdwclusfcl && ijk!=ijk2 && fabs((*sss0)[ijk]-(*sss0)[ijk2])<fdsclusfcl)
		  {
		    if ((*ccc0)[ijk]<0 && (*ccc0)[ijk2]<0)
		      {
			(*ccc0)[ijk]=ijk+1;
			(*ccc0)[ijk2]=ijk+1;
		      }
		    else if ((*ccc0)[ijk]>0 && (*ccc0)[ijk2]<0)
		      {
			(*ccc0)[ijk2]=(*ccc0)[ijk];
		      }
		    else if ((*ccc0)[ijk]<0 && (*ccc0)[ijk2]>0)
		      {
			(*ccc0)[ijk]=(*ccc0)[ijk2];
		      }
		    else if ((*ccc0)[ijk]>0 && (*ccc0)[ijk2]>0)
		      {
			for (unsigned int ijk3=0; ijk3<www0->size(); ijk3++)
			  {
			    if(((*ccc0)[ijk3])==((*ccc0)[ijk2]))(*ccc0)[ijk3]=(*ccc0)[ijk];
			  }
		      }
		  }
	      }
	  }
	
	
	for (unsigned int ijk=0; ijk<www2->size(); ijk++)
	  {
            ////std::cout << " WWW2 " << ijk << std::endl;
	    for (unsigned int ijk2=0; ijk2<www2->size(); ijk2++)
	      {
		if (fabs((*www2)[ijk]-(*www2)[ijk2])<fdwclusfcl && ijk!=ijk2 && fabs((*sss2)[ijk]-(*sss2)[ijk2])<fdsclusfcl)
		  {
		    if ((*ccc2)[ijk]<0 && (*ccc2)[ijk2]<0)
		      {
			(*ccc2)[ijk]=ijk+1;
			(*ccc2)[ijk2]=ijk+1;
		      }
		    else if ((*ccc2)[ijk]>0 && (*ccc2)[ijk2]<0)
		      {
			(*ccc2)[ijk2]=(*ccc2)[ijk];
		      }
		    else if ((*ccc2)[ijk]<0 && (*ccc2)[ijk2]>0)
		      {
			(*ccc2)[ijk]=(*ccc2)[ijk2];
		      }
		    else if ((*ccc2)[ijk]>0 && (*ccc2)[ijk2]>0)
		      {
			for (unsigned int ijk3=0; ijk3<www2->size(); ijk3++)
			  {
			    if(((*ccc2)[ijk3])==((*ccc2)[ijk2]))(*ccc2)[ijk3]=(*ccc2)[ijk];
			  }
		      }
		  }
	      }
	  }
	
	
	Int_t clusters_creation[4][30000];
	//Int_t clusters_avewire[4][1000];
	Int_t clusters_swire[4][30000];
	Int_t clusters_lwire[4][30000];
	Int_t clusters_ssample[4][30000];
	Int_t clusters_lsample[4][30000];
	
	Int_t clusters_nn[1000];
	Int_t clusters_vi[1000];
	Int_t clusters_qq[1000];
	//Int_t clusters_avewire[4][1000];
	Int_t clusters_dw[1000];
	Int_t clusters_ds[1000];
        Double_t clusters_mw[1000];
        Double_t clusters_ms[1000];
        Double_t clusters_mintime[1000];
	
	
	for (unsigned int ijk=0; ijk<4; ijk++) {
          for (unsigned int ijk2=0; ijk2<30000; ijk2++) {
	    clusters_creation[ijk][ijk2]=0;
	    clusters_swire[ijk][ijk2]=100000;
	    clusters_lwire[ijk][ijk2]=0;
	    clusters_ssample[ijk][ijk2]=100000;
	    clusters_lsample[ijk][ijk2]=0;
	    if (ijk2<1000 && ijk>2) {
	      clusters_nn[ijk2]=-1;
	      clusters_vi[ijk2]=-1;
	      clusters_dw[ijk2]=-1;
	      clusters_ds[ijk2]=-1;
	      clusters_qq[ijk2]=-1;
              clusters_mw[ijk2]=-1;
              clusters_ms[ijk2]=-1;
              clusters_mintime[ijk2]=-1;
	    }
          }
	}
	
	for (unsigned int ijk=0; ijk<www0->size(); ijk++)
	  {
	    if ((*ccc0)[ijk]>0) {
              int numero=(*ccc0)[ijk];
              clusters_creation[0][numero]+=1;
              if((*www0)[ijk]<clusters_swire[0][numero])clusters_swire[0][numero]=(*www0)[ijk];
              if((*www0)[ijk]>clusters_lwire[0][numero])clusters_lwire[0][numero]=(*www0)[ijk];
              if((*sss0)[ijk]<clusters_ssample[0][numero])clusters_ssample[0][numero]=(*sss0)[ijk];
              if((*sss0)[ijk]>clusters_lsample[0][numero])clusters_lsample[0][numero]=(*sss0)[ijk];
	    }
	  }
	for (unsigned int ijk=0; ijk<www2->size(); ijk++)
	  {
	    if ((*ccc2)[ijk]>0) {
              int numero=(*ccc2)[ijk];
              clusters_creation[2][numero]+=1;
              if((*www2)[ijk]<clusters_swire[2][numero])clusters_swire[2][numero]=(*www2)[ijk];
              if((*www2)[ijk]>clusters_lwire[2][numero])clusters_lwire[2][numero]=(*www2)[ijk];
              if((*sss2)[ijk]<clusters_ssample[2][numero])clusters_ssample[2][numero]=(*sss2)[ijk];
              if((*sss2)[ijk]>clusters_lsample[2][numero])clusters_lsample[2][numero]=(*sss2)[ijk];
	    }
	  }
	
	int quanti_clusters=0;
	for (unsigned int ijk=0; ijk<4; ijk++) {
          for (unsigned int ijk2=0; ijk2<30000; ijk2++) {
	    if(clusters_creation[ijk][ijk2]>50)
              {
		////std::cout<<clusters_creation[ijk][ijk2] << " CLUSTER MINE " << clusters_swire[ijk][ijk2] << " " << clusters_lwire[ijk][ijk2] << " " << clusters_ssample[ijk][ijk2] << " " << clusters_lsample[ijk][ijk2] << std::endl;
		clusters_qq[quanti_clusters]=clusters_creation[ijk][ijk2];
		clusters_vi[quanti_clusters]=ijk;
		clusters_nn[quanti_clusters]=ijk2;
		clusters_dw[quanti_clusters]=clusters_lwire[ijk][ijk2]-clusters_swire[ijk][ijk2];
		clusters_ds[quanti_clusters]=clusters_lsample[ijk][ijk2]-clusters_ssample[ijk][ijk2];
                clusters_mw[quanti_clusters]=(clusters_lwire[ijk][ijk2]+clusters_swire[ijk][ijk2])*0.5;
                clusters_ms[quanti_clusters]=(clusters_lsample[ijk][ijk2]+clusters_ssample[ijk][ijk2])*0.5;
                clusters_mintime[quanti_clusters]=(clusters_ssample[ijk][ijk2]);
		quanti_clusters+=1;
              }
          }
	}
	
	
	for(int icl = 0; icl < quanti_clusters; ++icl){
	  if (clusters_qq[icl]>0) {
	    //std::cout << " CLUSTER INFO " << icl << " " << clusters_qq[icl] << " " << clusters_vi[icl] << " " << clusters_dw[icl] << " " << clusters_ds[icl] << " " << clusters_nn[icl] << std::endl;
	    int tpc_number=clusters_vi[icl];//qui andrebbe messo il numero di TPC
	    
	    if (clusters_ds[icl]>2250 && clusters_dw[icl]>100)
	      {//if analisi
		
		std::vector<float> *whc=new std::vector<float>;
		std::vector<float> *shc=new std::vector<float>;
		std::vector<float> *ahc=new std::vector<float>;
		std::vector<float> *fahc=new std::vector<float>;
		
		if (clusters_vi[icl]==0) {
		  for (unsigned int ijk=0; ijk<www0->size(); ijk++) {
		    if ((*ccc0)[ijk]==clusters_nn[icl]) {
                      whc->push_back((*www0)[ijk]);
                      shc->push_back((*sss0)[ijk]);
                      ///ahc->push_back((*hhh0)[ijk]);
                      ahc->push_back((*aaa0)[ijk]);
                      fahc->push_back((*hhh0)[ijk]);
		    }
		  }
		}
		if (clusters_vi[icl]==2) {
		  for (unsigned int ijk=0; ijk<www2->size(); ijk++) {
		    if ((*ccc2)[ijk]==clusters_nn[icl]) {
                      whc->push_back((*www2)[ijk]);
                      shc->push_back((*sss2)[ijk]);
                      ahc->push_back((*aaa2)[ijk]);
                      //ahc->push_back((*hhh2)[ijk]);
                      fahc->push_back((*hhh2)[ijk]);
		    }
		  }
		}
		
		
		////std::cout << " CLUSTER INFO " << icl << " " << clusters_qq[icl] << " " << whc->size() << std::endl;
		
		
		if(whc->size()>100)//prima 0
		  {
		    float pendenza=0;float intercetta=0;int found_ok=0;
		    std::vector<int> *escluse=new std::vector<int>;
/*
for(int j=0;j<(int)whc->size();j++)
{
if(((*shc)[j]-clusters_mintime[icl]<200) || ((*shc)[j]-clusters_mintime[icl])>2150)escluse->push_back(j);
}
*/
		    for(int j=0;j<(int)whc->size();j++)
		      {
			if(found_ok<1)
			  {
			    Double_t wires[10000];
			    Double_t samples[10000];
			    Double_t ex[10000];
			    Double_t quale[10000];
			    Double_t ey[10000];
			    int quanti=0;
			    for(int k=0;k<(int)whc->size();k++)
			      {
				if(Nothere(escluse,k)==3)
				  {
				    wires[quanti]=(*whc)[k]*3;
				    samples[quanti]=(*shc)[k]*0.628;
				    ex[quanti]=0;
				    ey[quanti]=0;
				    quale[quanti]=k;
				    quanti+=1;
				  }
			      }
			    //std::cout << quanti << " FIRST FIT " << std::endl;
			    TGraphErrors *gr3 = new TGraphErrors(quanti,wires,samples,ex,ey);
			    TFitResultPtr fp1 = gr3->Fit("pol1", "MQ");//std::cout << int(fp1) <<std::endl;

			    //gr3->Fit("pol1","Q");
			    TF1 *fitfunc = gr3->GetFunction("pol1");
			    pendenza=fitfunc->GetParameter(1);
			    intercetta=fitfunc->GetParameter(0);
			    float distance_maximal=fdisfcl;
			    int quella_a_distance_maximal=0;
			    int found_max=0;
			    for(int jj=0;jj<quanti;jj++)
			      {
				if((abs((pendenza)*(wires[jj])-samples[jj]+intercetta)/sqrt((pendenza)*pendenza+1))>distance_maximal)
				  {
				    found_max=1;
				    quella_a_distance_maximal=quale[jj];
				    distance_maximal=(abs(pendenza*wires[jj]-samples[jj]+intercetta)/sqrt(pendenza*pendenza+1));
				  }
			      }
			    if(found_max==1)escluse->push_back(quella_a_distance_maximal);
			    if(found_max==0)found_ok=1;
			  }
		      }
		    ////std::cout << escluse->size() << " escluse " << whc->size() << " " << found_ok << std::endl;
		    delete escluse;
		    std::vector<float> *hittime=new std::vector<float>;
		    std::vector<float> *hitwire=new std::vector<float>;
		    std::vector<float> *hitarea=new std::vector<float>;
		    std::vector<float> *hittimegood=new std::vector<float>;
		    std::vector<float> *hitareagood=new std::vector<float>;
		    std::vector<float> *hitwiregood=new std::vector<float>;
		    ///////std::cout <<  found_ok << std::endl;
		    ///////std::cout <<  pendenza << std::endl;
		    ///////std::cout <<  intercetta << std::endl;
		    if(found_ok==1)
		      {
			for(int kkk=0;kkk<(int)whc->size();kkk++)
			  {
			    if((*ahc)[kkk]>0)// && (((*shc)[kkk]-clusters_mintime[icl]<200) || ((*shc)[kkk]-clusters_mintime[icl])<2150))
			      {
				float distance=(abs(pendenza*((*whc)[kkk]*3)-(*shc)[kkk]*0.628+intercetta))/sqrt(pendenza*pendenza+1);
				if(distance<=fdisfcl)
				  {
				    //cout << log((*ahc)[kkk]/(0.4*peach)) << endl;
				    hittime->push_back((*shc)[kkk]*0.4);
				    hitwire->push_back((*whc)[kkk]);
				    hitarea->push_back((*ahc)[kkk]);
                                    h_hit_height->Fill((*fahc)[kkk]);
                                    h_hit_area->Fill((*ahc)[kkk]);
                                    h_hit_height_area->Fill((*fahc)[kkk],(*ahc)[kkk]);

				  }
			      }
			  }
			//float result_rms=0.14;
			////std::cout << result_rms << endl;
			Double_t area[10000];
			Double_t nologarea[10000];
			Double_t tempo[10000];
			Double_t ex[10000];
			//Double_t quale[10000];
			Double_t ey[10000];
			Double_t ek[10000];
			Double_t ez[10000];
			////std::cout <<  hitarea->size() << " dimensione hitarea" << std::endl;

			////std::cout<<""<<std::endl;
			////std::cout<<"HERE line 802"<<std::endl;
			////std::cout<<""<<std::endl;
                        //std::ofstream purh("dump_purity_hits.out",std::ios::app);

			if(hitarea->size()>100)//prima 30
			  {
			    h_ratio->Fill(((float)whc->size())/((float)clusters_dw[icl]));
                            h_ratio_3->Fill(clusters_ds[icl],((float)whc->size())/((float)clusters_dw[icl]));

                            ////std::cout << "RATIO INFO 1 " << (float)hitarea->size() << " " << (float)clusters_dw[icl] << " " << (((float)hitarea->size())/((float)clusters_dw[icl])) << std::endl;	
			    float minimo=100000;
			    float massimo=0;
                            float wire_minimo=100000;
                            float wire_massimo=0;

                            float delta_sample_selected=0;
                            float delta_wire_selected=0;
                            float wire_del_massimo=-1;
                            float wire_del_minimo=-1;
                            float sample_massimo=-1;
                            float sample_minimo=-1;
                            int quante_hit_nel_range_tempo=0;
			    purh<< evt.run() << " " << evt.subRun() << "  " << evt.event() << "  -1 " <<  " -1 " << " -1 " << " -1 " << std::endl;
                            for(int kk=0;kk<(int)hitarea->size();kk++)
                              {
				if(fdumphitsfcl==1)purh<< evt.run() << " SELECTED " << evt.subRun() << "  " << evt.event() << "  " << tpc_number <<  " " << (*hitwire)[kk] << " " << (*hittime)[kk] << " " << (*hitarea)[kk] << std::endl;

quante_hit_nel_range_tempo+=1;
                                if((*hittime)[kk]>massimo)
                              {
                               massimo=(*hittime)[kk];
                               wire_del_massimo=(*hitwire)[kk];
                              }


                                if((*hittime)[kk]<minimo)
                              {
                               minimo=(*hittime)[kk];
                               wire_del_minimo=(*hitwire)[kk];
                              }


                                if((*hitwire)[kk]>wire_massimo)
                              {
                               wire_massimo=(*hitwire)[kk];
                              }


                                if((*hitwire)[kk]<wire_minimo)
                              {
                               wire_minimo=(*hitwire)[kk];
                              }

                              }

                                delta_sample_selected=(massimo-minimo)/0.4;
                                delta_wire_selected=wire_del_massimo-wire_del_minimo;
                                sample_massimo=massimo/0.4;
                                sample_minimo=minimo/0.4;
                            h_ratio_after->Fill(((float)quante_hit_nel_range_tempo)/fabs(delta_wire_selected));
                            h_ratio_after_2->Fill(((float)quante_hit_nel_range_tempo)/fabs(wire_massimo-wire_minimo+1));
                            //h_ratio_after_3->Fill(delta_sample_selected,((float)quante_hit_nel_range_tempo)/fabs(wire_massimo-wire_minimo+1));

                            ////std::cout << "RATIO INFO 2 " << (float)hitarea->size() << " " << fabs(delta_wire_selected) << " " << (((float)hitarea->size())/fabs(delta_wire_selected)) << std::endl;

			    ////std::cout << hitarea->size() << std::endl;
			    //int gruppi=hitarea->size()/50;
			    int gruppi=fgruppifcl;//originale 8
			    ////std::cout << gruppi << std::endl;
			    
			    float steptime=(massimo-minimo)/(gruppi+1);
			    ////std::cout << steptime << " steptime " << minimo << " " << massimo << std::endl;
			    float starting_value_tau=fValoretaufcl;
			    
			    ////std::cout << starting_value_tau << " VALORE INDICATIVO TAU " << std::endl;
			    //if(tpc_number==2 || tpc_number==5)starting_value_tau=6500;
			    //if(tpc_number==10 || tpc_number==13)starting_value_tau=5700;
			    for(int stp=0;stp<=gruppi;stp++)
			      {
				std::vector<float>* hitpertaglio=new std::vector<float>;
				////std::cout << 500+stp*steptime << " time " << 500+(stp+1)*(steptime) << std::endl;
				///////////std::cout << minimo+stp*steptime << " " << minimo+(stp+1)*(steptime) << std::endl;
				for(int kk=0;kk<(int)hitarea->size();kk++)
				  {
				    if((*hittime)[kk]>=(minimo+stp*steptime) && (*hittime)[kk]<=(minimo+(stp+1)*(steptime))) 
				      hitpertaglio->push_back((*hitarea)[kk]*exp((*hittime)[kk]/starting_value_tau));
				  }
				///////////std::cout << hitpertaglio->size() << std::endl;
				float tagliomax=FoundMeanLog(hitpertaglio,fmaxfcl);//0.9
				float tagliomin=FoundMeanLog(hitpertaglio,fminfcl);//0.05
				//float tagliomin=0;
				//float tagliomax=1000000;
				delete hitpertaglio;
				////std::cout << tagliomax << " t " << std::endl;
				for(int kk=0;kk<(int)hitarea->size();kk++)
				  {
				    ////std::cout << (*hittime)[kk] << " " << (*hitwire)[kk] << " " << (*hitarea)[kk] << " " << (minimo+stp*steptime) << " " << (minimo+(stp+1)*steptime) << " " << (*hitarea)[kk]*exp((*hittime)[kk]/starting_value_tau) << std::endl;
				    if((*hittime)[kk]>(minimo+stp*steptime) && 
				       (*hittime)[kk]<(minimo+(stp+1)*steptime) &&
				       (*hitarea)[kk]*exp((*hittime)[kk]/starting_value_tau)<tagliomax && 
				       (*hitarea)[kk]*exp((*hittime)[kk]/starting_value_tau)>tagliomin)
				      {
					////std::cout << ((*hitarea)[kk]*exp((*hittime)[kk]/1400)) << " GOOD " << (*hitarea)[kk] << " " << (*hittime)[kk] << std::endl;
					hitareagood->push_back((*hitarea)[kk]);
					hittimegood->push_back((*hittime)[kk]);
					hitwiregood->push_back((*hitwire)[kk]);
				      }
				  }
			      }
			    ////std::cout << hitareagood->size() << " hitareagood" << std::endl;    
if(delta_sample_selected>1900)
{
                    for(int k=0;k<(int)whc->size();k++)
                      {
                        h_hittime->Fill((*shc)[k]-clusters_mintime[icl]);
                      }
                    for(int k=0;k<(int)hittime->size();k++)
                      {
                        h_hittime_2->Fill((*hittime)[k]/0.4-clusters_mintime[icl]);
                      }

                    for(int k=0;k<(int)hittimegood->size();k++)
                      {
                        h_hittime_3->Fill((*hittimegood)[k]/0.4-clusters_mintime[icl]);
                      }
}

			    for(int k=0;k<(int)hitareagood->size();k++)
			      {
                                if(fdumphitsfcl==1)purh<< evt.run() << " AFTER CLEAN " << evt.subRun() << "  " << evt.event() << "  " << tpc_number <<  " " << (*hitwiregood)[k] << " " << (*hittimegood)[k] << " " << (*hitareagood)[k] << std::endl;
				//if((*hittimegood)[k]-600*0.4<=1000)//correzione 15/08
				    tempo[k]=(*hittimegood)[k];
				    area[k]=log((*hitareagood)[k]);
				    ////std::cout << (*hitareagood)[k] << " " << area[k] << std::endl;
				    nologarea[k]=((*hitareagood)[k]);
				    ex[k]=0;
				    ez[k]=60;
				    ey[k]=0.23;
			      }
			    ////std::cout<<""<<std::endl;
			    ////std::cout<<"HERE line 872"<<std::endl;
			    ////std::cout<<""<<std::endl;
			    //std::cout<<hitareagood->size() <<" SECOND FIT "<<std::endl;
			    TGraphErrors *gr31 = new TGraphErrors(hitareagood->size(),tempo,area,ex,ey);
			    //TGraphErrors *gr4 = new TGraphErrors(hitareagood->size(),tempo,nologarea,ex,ey);
			    TFitResultPtr fp31 = gr31->Fit("pol1", "MQ");//std::cout << int(fp31) <<std::endl;
			    //gr31->Fit("pol1","Q");
			    TF1 *fit = gr31->GetFunction("pol1");
			    float slope_purity=fit->GetParameter(1);
			    //float error_slope_purity=fit->GetParError(1);
			    float intercetta_purezza=fit->GetParameter(0);
			    
			    TH1F *h111 = new TH1F("h111","delta aree",200,-10,10);
			    //float sum_per_rms_test=0;
                            int quanti_in_h111=0;
			    for(int k=0;k<(int)hitareagood->size();k++)
			      {
				h111->Fill(area[k]-slope_purity*tempo[k]-intercetta_purezza);
				//sum_per_rms_test+=(area[k]-slope_purity*tempo[k]-intercetta_purezza)*(area[k]-slope_purity*tempo[k]-intercetta_purezza);
                                if((area[k]-slope_purity*tempo[k]-intercetta_purezza)>-10 && (area[k]-slope_purity*tempo[k]-intercetta_purezza)<10)quanti_in_h111+=1;
			      }
                       
                        float error=100.;
                        //std::cout << quanti_in_h111 << " quanti h111 " << error << std::endl;
                         if(quanti_in_h111>1){
			   TFitResultPtr fph111 = h111->Fit("gaus", "MQ");//std::cout << int(fph111) <<std::endl;
			 //h111->Fit("gaus","Q");
                        TF1 *fitg = h111->GetFunction("gaus");
                        error=fitg->GetParameter(2);
                        }
                        ////std::cout << " error " << error << std::endl;
                        //float error_2=sqrt(sum_per_rms_test/(hitareagood->size()-2));
                        ////std::cout << " error vero" << error_2 << std::endl;
                        h111->Delete();

			//std::cout<<hitareagood->size() <<" THIRD FIT "<<std::endl;


		      TGraphErrors *gr4 = new TGraphErrors(hitareagood->size(),tempo,nologarea,ex,ey);
		      TFitResultPtr fp = gr4->Fit("expo", "MQ");//std::cout << int(fp) <<std::endl;
		      //gr4->Fit("expo","Q");
		      TF1 *fite = gr4->GetFunction("expo");
		      slope_purity=fite->GetParameter(1);
		      intercetta_purezza=fite->GetParameter(0);
                      float mean_hit_area=0;
                      float size_hit_area=hitareagood->size();
                      int quanti_in_h111e=0;
		      TH1F *h111e = new TH1F("h111e","delta aree",100,-1000.,1000.);
		      for(int k=0;k<(int)hitareagood->size();k++)
			{
			  h111e->Fill(nologarea[k]-exp(slope_purity*tempo[k]+intercetta_purezza));
                          mean_hit_area+=nologarea[k]/size_hit_area;
                          if((nologarea[k]-exp(slope_purity*tempo[k]+intercetta_purezza))>-1000. && (nologarea[k]-exp(slope_purity*tempo[k]+intercetta_purezza))<1000)quanti_in_h111e+=1; 
			  //cout << nologarea[k]-exp(slope_purity*tempo[k]+intercetta_purezza) << endl;
			}
                      
                      float error_expo=1000.;
		      //std::cout << quanti_in_h111e << " quanti h111e " << error << std::endl;
                     if(quanti_in_h111e>1)
                      { 
			TFitResultPtr fp111e = h111e->Fit("gaus", "MQ");//std::cout << int(fp111e) <<std::endl;
		      //h111e->Fit("gaus","Q");
		      TF1 *fitge = h111e->GetFunction("gaus");
		      error_expo=fitge->GetParameter(2);
                      }
		     //std::cout << " errors " << error << " " << error_expo << std::endl;
		      h_errors->Fill(error_expo);
		      h111e->Delete();//fitge->Delete();fite->Delete();

                        for(int k=0;k<(int)hitareagood->size();k++)
                          {
                                ek[k]=error_expo;
                                ez[k]=error_expo;
                                ey[k]=error;
                          }
			////std::cout<<""<<std::endl;
			////std::cout<<"HERE line 906"<<std::endl;
			////std::cout<<""<<std::endl;
			//std::cout<<hitareagood->size() <<" FOURTH FIT "<<std::endl;

                        //TGraphErrors *gr32 = new TGraphErrors(hitareagood->size(),tempo,area,ex,ey);
			//TFitResultPtr fp2 = gr32->Fit("pol1", "MQ");std::cout << int(fp2) <<std::endl;
                        //gr32->Fit("pol1","Q");
                  
                        //TF1 *fit2 = gr32->GetFunction("pol1");
                        //float slope_purity_2=fit2->GetParameter(1);
                        //float error_slope_purity_2=fit2->GetParError(1);
                        //float intercetta_purezza_2=fit2->GetParameter(0);
                        //float chiquadro=fit2->GetChisquare()/(hitareagood->size()-2);
			//std::ofstream goodpuro("purity_results.out",std::ios::app);
                        std::ofstream goodpuro2("purity_results2.out",std::ios::app);
			
                        ////std::cout << -1/slope_purity_2 << std::endl;
                        ////std::cout << -1/(slope_purity_2+error_slope_purity_2)+1/slope_purity_2 << std::endl;
                        ////std::cout << 1/slope_purity_2-1/(slope_purity_2-error_slope_purity_2) << std::endl;
 			//std::cout<<hitareagood->size() <<" FIFTH FIT "<<std::endl;
                       TGraphAsymmErrors *gr41 = new TGraphAsymmErrors (hitareagood->size(),tempo,nologarea,ex,ex,ez,ek);
		       TFitResultPtr fp41 = gr41->Fit("expo", "MQ");//std::cout << int(fp41) <<std::endl;
		       //gr41->Fit("expo","Q");
                        TF1 *fitexo = gr41->GetFunction("expo");
                        float slope_purity_exo=fitexo->GetParameter(1);
                        float error_slope_purity_exo=fitexo->GetParError(1);
                        //fRunSubPurity2->Fill(evt.run(),evt.subRun(),-slope_purity_exo*1000.);
                        //fRunSubPurity->Fill(evt.run(),evt.subRun(),-slope_purity_2*1000.);
                        ////std::cout << -1/slope_purity_exo << std::endl;
                        //std::cout << slope_purity_exo << " " << error_slope_purity_exo<< std::endl;
                        ////std::cout << 1/slope_purity_exo-1/(slope_purity_exo-error_slope_purity_exo) << std::endl;
                        ////std::cout << fitexo->GetChisquare()/(hitareagood->size()-2) << std::endl;
			
			
                        if(/*fabs(slope_purity_2)<0.01 || */fabs(slope_purity_exo)<0.01)
			  {
			    //if(fabs(slope_purity_exo)<0.01)purityvalues->Fill(-slope_purity_exo*1000.);
                            //if(fabs(slope_purity_2)<0.01)goodpuro << evt.run() << " " << evt.subRun() << "  " << evt.event() << "  " << tpc_number << "  " << slope_purity_2 << "  " << error_slope_purity_2 << " " << chiquadro << " " << clusters_dw[icl] << " " << clusters_ds[icl] << std::endl;
			    
			    if(fabs(slope_purity_exo)<0.01)goodpuro2<< evt.run() << " " << evt.subRun() << " " << evt.event() << " " << tpc_number+fcryofcl*10 << " " << slope_purity_exo << " " << error_slope_purity_exo << " " << fitexo->GetChisquare()/(hitareagood->size()-2) << " " << clusters_dw[icl] << " " << clusters_ds[icl] << " " << clusters_mw[icl] << " " << clusters_ms[icl] << " " << delta_wire_selected << " " << delta_sample_selected << " " << sample_minimo << " " << sample_massimo << " " << wire_del_minimo << " " << wire_del_massimo << " " << wire_minimo << " " << wire_massimo << " " << whc->size() << " " << hitarea->size() << " " << quante_hit_nel_range_tempo << " " << error_expo << " " << clusters_mintime[icl] << " " << mean_hit_area << std::endl;

                  
                  if(fabs(slope_purity_exo)<0.01)
                  {
                      run_tree=evt.run();
                      //subrun_tree=evt.subRun();
                      event_tree=evt.event();
                      tpc_tree=tpc_number+fcryofcl*10;
                      slope_tree=-slope_purity_exo*1000;
                      errslope_tree=error_slope_purity_exo*1000;
                      chi_tree=fitexo->GetChisquare()/(hitareagood->size()-2);
                      dtime_tree=delta_sample_selected;
                      dwire_tree=wire_massimo-wire_minimo+1;
                      awire_tree=(wire_massimo+wire_minimo)*0.5;
                      earea_tree=error_expo;
                      marea_tree=mean_hit_area;
                      qhits_tree=hitarea->size();
		      time_tree=evt.time().timeHigh()-1600000000;
                      fpurTree->Fill();
                  }
                  if(fabs(slope_purity_exo)<0.01)purityvalues2->Fill(-slope_purity_exo*1000.);

			    if(fabs(slope_purity_exo)<0.01 && tpc_number==0)puritytpc0->Fill(-slope_purity_exo*1000.);
			    if(fabs(slope_purity_exo)<0.01 && tpc_number==1)puritytpc1->Fill(-slope_purity_exo*1000.);
			    if(fabs(slope_purity_exo)<0.01 && tpc_number==2)puritytpc2->Fill(-slope_purity_exo*1000.);
			    if(fabs(slope_purity_exo)<0.01 && tpc_number==3)puritytpc3->Fill(-slope_purity_exo*1000.);
			    
			  }
		        if(fabs(slope_purity_exo)<0.01)
                        {	

			  
			anab::TPCPurityInfo purity_info;
			purity_info.Run = evt.run();
			purity_info.Subrun = evt.subRun();
			purity_info.Event = evt.event();
			purity_info.TPC = tpc_number;
			purity_info.Wires = delta_wire_selected;
                        purity_info.Ticks = delta_sample_selected;
			
			//if(purity_info.TPC<2) purity_info.Cryostat=0;
			//else purity_info.Cryostat=1;
			
			purity_info.Cryostat=fcryofcl;

			// near/far from cathode tracks                                                                                        
/*
			  if(delta_wire_selected< 0){
			    purity_info.AttenuationNEAR = slope_purity_exo*-1000.;
			  }
			  else purity_info.AttenuationNEAR = 0;
			  if(delta_wire_selected>0){
                            purity_info.AttenuationFAR = slope_purity_exo*-1000.;
                          }
			  else purity_info.AttenuationFAR = 0;

*/

       			purity_info.Attenuation = slope_purity_exo*-1.;
			purity_info.FracError = error_slope_purity_exo / slope_purity_exo;
                        //purity_info.Attenuation_2 = slope_purity_2*-1.;
                        //purity_info.FracError_2 = error_slope_purity_2 / slope_purity_2;

			//purity_info.Wires = clusters_dw[icl];
			//                        purity_info.Ticks = clusters_ds[icl];

			if(fFillAnaTuple)
			  purityTuple->Fill(purity_info.Run,purity_info.Event,purity_info.TPC,purity_info.Wires,purity_info.Ticks,purity_info.Attenuation);

			////std::cout << "Calling again after filling attenuation … " << std::endl;
			purity_info.Print();
			outputPtrVector->push_back(purity_info);
                        }
			////std::cout << ts << " is time event " << std::endl;
                        //goodpur << -1/slope_purity_exo << std::endl;
                        //goodpur << -1/(slope_purity_exo+error_slope_purity_exo)+1/slope_purity_exo << std::endl;
                        //goodpur << 1/slope_purity_exo-1/(slope_purity_exo-error_slope_purity_exo) << std::endl;
                        //goodpur << timeevent << " is time event " << std::endl;
			  }
		      }
		    ////std::cout << "Delete hit stuff." << std::endl;

		    delete hittime;
		    delete hitarea;
		    delete hittimegood;
		    delete hitareagood;
		    delete hitwiregood;
		    delete hitwire;
		  }
		
		////std::cout << "Delete cluster stuff." << std::endl;
		
		delete shc;
		delete ahc;
		delete fahc;
		delete whc;
	       
	      }//fine if ananlisi
	  }
	}

	////std::cout << "Delete big stuff." << std::endl;

	delete www0;
	delete sss0;
	delete hhh0;
	delete ehh0;
	delete ccc0;
	delete www2;
	delete sss2;
	delete hhh2;
	delete ehh2;
	delete ccc2;
	delete aaa0;
	delete aaa2;

      }

    //std::cout << "Checking everything in the output..." << std::endl;
    //std::cout << "There are " << outputPtrVector->size() << " objects in the output vector." << std::endl;
    /* //don't need this printed here.    
    for (size_t i_info = 0; i_info<outputPtrVector->size(); ++i_info){
      auto info = outputPtrVector->at(i_info);
      info.Print();
    }
    */
    
   //put info onto the event
    evt.put(std::move(outputPtrVector));
  } // produces
  
} //end namespace

namespace icarus{

  DEFINE_ART_MODULE(ICARUSPurityDQM)
  
}

