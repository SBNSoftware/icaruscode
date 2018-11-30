////////////////////////////////////////////////////////////////////////
// Class:       CRTSimpleTrackProducer
// Module Type: producer
// File:        CRTSimpleTrackProducer_module.cc
// Description: Module for constructiong over-simplified CRT tracks.
// Copied from CRTTrackProducer by David Lorca Galindo 
//  Edited by Michelle Stancari April 3, 2018
//  Added some SBND specific stuff - Tom Brooks
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"

#include "sbndcode/CRT/CRTProducts/CRTHit.hh"
#include "sbndcode/CRT/CRTProducts/CRTTrack.hh"
#include "sbndcode/CRT/CRTProducts/CRTTzero.hh"

#include "TTree.h"
#include "TVector3.h"

#include <iostream>
#include <stdio.h>
#include <sstream>
#include <vector>
#include <map>
#include <utility>
#include <cmath> 
#include <memory>

namespace sbnd {


class CRTTrackProducer : public art::EDProducer {

public:

  explicit CRTTrackProducer(fhicl::ParameterSet const & p);

  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.
  // Plugins should not be copied or assigned.
  CRTTrackProducer(CRTTrackProducer const &) = delete;
  CRTTrackProducer(CRTTrackProducer &&) = delete;
  CRTTrackProducer & operator = (CRTTrackProducer const &) = delete;
  CRTTrackProducer & operator = (CRTTrackProducer &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

  // Function to make createing CRTHits easier
  crt::CRTHit FillCrtHit(std::vector<uint8_t> tfeb_id, std::map<uint8_t, 
                         std::vector<std::pair<int,float>>> tpesmap, float peshit, 
                         double time, double x, double ex, double y, double ey, 
                         double z, double ez, std::string tagger); 
  
  // Function to make creating CRTTracks easier
  crt::CRTTrack FillCrtTrack(crt::CRTHit hit1, crt::CRTHit hit2, bool complete);

  // Function to average hits within a certain distance of each other
  std::vector<crt::CRTHit> AverageHits(std::vector<art::Ptr<crt::CRTHit>> hits);

  // Take a list of hits and find average parameters
  crt::CRTHit DoAverage(std::vector<art::Ptr<crt::CRTHit>> hits);

  // Create CRTTracks from list of hits
  std::vector<crt::CRTTrack> CreateTracks(std::vector<crt::CRTHit> hits);

  // Calculate the tagger crossing point of CRTTrack candidate
  TVector3 CrossPoint(crt::CRTHit hit, TVector3 start, TVector3 diff);

private:

  // Declare member data here.
  //  art::ServiceHandle<art::TFileService> tfs;
  std::string  fDataLabelHits;
  std::string  fDataLabelTZeros;
  int          fTrackMethodType;
  int          fStoreTrack;
  double       fTimeLimit;
  double       fAverageHitDistance;
  double       fDistanceLimit;
  bool         fUseTopPlane;

}; // class CRTTrackProducer


// Function to calculate average and rms from a vector of values
void vmanip(std::vector<float> v, float* ave, float* rms);


// Average crt hit structure
struct CRTavehit{
  uint32_t ts0_ns;
  uint16_t ts0_ns_err;
  int32_t ts1_ns; 
  uint16_t ts1_ns_err;                                                        
  
  float x_pos;
  float x_err;
  float y_pos;
  float y_err;
  float z_pos;
  float z_err;
  float pe;
  int plane;
  std::string tagger;
} tempah;


// Function to make filling average hit easier
CRTavehit fillme(uint32_t i, uint16_t j, int32_t k, uint16_t l, float a,
                 float b, float c, float d, float e, float f, float g, 
                 int p, std::string t);


// Function to copy average hits
CRTavehit copyme(crt::CRTHit myhit);


// Function to make creating CRTTracks easier
crt::CRTTrack shcut(CRTavehit ppA,CRTavehit ppb, uint32_t time0s,uint16_t terr);


// Constructor
CRTTrackProducer::CRTTrackProducer(fhicl::ParameterSet const & p)
{  

  // Initialize member data here.
  fDataLabelHits      = p.get<std::string>("DataLabelHits");      // CRTHit producer module name
  fDataLabelTZeros    = p.get<std::string>("DataLabelTZeros");    // CRTTzero producer module name
  fStoreTrack         = p.get<int>        ("StoreTrack");         // method 1 = all, method 2 = ave, method 3 = pure, method 4 = top plane
  fTrackMethodType    = p.get<int>        ("TrackMethodType");    // Print stuff
  fTimeLimit          = p.get<double>     ("TimeLimit");          // Maximum time difference to combine hits into a tzero [us]
  fAverageHitDistance = p.get<double>     ("AverageHitDistance"); // Maximum distance to avarage hits on the same tagger over [cm]
  fDistanceLimit      = p.get<double>     ("DistanceLimit");      // Maximum distance from CRTTrack to absorb a new CRTHit [cm]
  fUseTopPlane        = p.get<bool>       ("UseTopPlane");        // Use hits from the top plane (SBND specific)
  
  // Call appropriate produces<>() functions here.
  if(fStoreTrack == 1) 
    produces< std::vector<crt::CRTTrack>   >();

} // CRTTrackProducer()


void CRTTrackProducer::produce(art::Event & evt)
{

  // For validation
  int nTrack = 0;
  int nCompTrack = 0;
  int nIncTrack = 0;

  // CRTTrack collection on this event                                                                         
  std::unique_ptr<std::vector<crt::CRTTrack> > CRTTrackCol(new std::vector<crt::CRTTrack>);

  // Implementation of required member function here.
  art::Handle< std::vector<crt::CRTHit> > rawHandle;
  evt.getByLabel(fDataLabelHits, rawHandle); //what is the product instance name? no BernZMQ

  // Check to make sure the data we asked for is valid                                                                         
  if(!rawHandle.isValid()){
    mf::LogWarning("CRTTrackProducer")
      <<"No CRTHits from producer module "<<fDataLabelHits;
    return;
  }

  // Track method 4 = SBND method with top plane (doesn't use CRTTzero)
  if(fTrackMethodType == 4){

    //Get the CRT hits from the event
    std::vector<art::Ptr<crt::CRTHit> > hitlist;
    if (evt.getByLabel(fDataLabelHits, rawHandle))
      art::fill_ptr_vector(hitlist, rawHandle);

    std::vector<std::vector<art::Ptr<crt::CRTHit>>> CRTTzeroVect;
    std::vector<int> npvec;
    int iflag[2000] = {};

    // Loop over crt hits
    for(size_t i = 0; i<hitlist.size(); i++){
      if(iflag[i] == 0){
        std::vector<art::Ptr<crt::CRTHit>> CRTTzero;
        std::map<std::string, int> nPlanes;
        double time_ns_A = hitlist[i]->ts1_ns;
        iflag[i]=1;
        CRTTzero.push_back(hitlist[i]);
        nPlanes[hitlist[i]->tagger]++;

        // Sort into a Tzero collection
        // Loop over all the other CRT hits
        for(size_t j = i+1; j<hitlist.size(); j++){
          if(iflag[j] == 0){
            // If ts1_ns - ts1_ns < diff then put them in a vector
            double time_ns_B = hitlist[j]->ts1_ns;
            double diff = std::abs(time_ns_B - time_ns_A) * 1e-3; // [us]
            if(diff < fTimeLimit){
              iflag[j] = 1;
              CRTTzero.push_back(hitlist[j]);
              nPlanes[hitlist[j]->tagger]++;
            }
          }
        }

        int np = 0;
        for(auto &nPlane : nPlanes){
          if(nPlane.second>0 && nPlane.first != "volTaggerTopHigh_0") np++;
        }
        CRTTzeroVect.push_back(CRTTzero);
        npvec.push_back(np);
      }
    }

    // Loop over tzeros
    for(size_t i = 0; i<CRTTzeroVect.size(); i++){

      //loop over hits for this tzero, sort by tagger
      std::map<std::string, std::vector<art::Ptr<crt::CRTHit>>> hits;
      for (size_t ah = 0; ah< CRTTzeroVect[i].size(); ++ah){        
        std::string ip = CRTTzeroVect[i][ah]->tagger;       
        hits[ip].push_back(CRTTzeroVect[i][ah]);
      } // loop over hits
      
      //loop over planes and calculate average hits
      std::vector<crt::CRTHit> allHits;
      for (auto &keyVal : hits){
        std::string ip = keyVal.first;
        std::vector<crt::CRTHit> ahits = AverageHits(hits[ip]);
        if(fUseTopPlane && ip == "volTaggerTopHigh_0"){ 
          allHits.insert(allHits.end(), ahits.begin(), ahits.end());
        }
        else if(ip != "volTaggerTopHigh_0"){ 
          allHits.insert(allHits.end(), ahits.begin(), ahits.end());
        }
      }

      //Create tracks with hits at the same tzero
      std::vector<crt::CRTTrack> trackCandidates = CreateTracks(allHits);
      nTrack += trackCandidates.size();
      for(size_t j = 0; j < trackCandidates.size(); j++){
        CRTTrackCol->emplace_back(trackCandidates[j]);
        if(trackCandidates[j].complete) nCompTrack++;
        else nIncTrack++;
      }
    }

  //Older track reconstruction methods from MicroBooNE
  }else{
    //Get list of tzeros             
    //  std::vector<crt::CRTHit> const& CRTHitCollection(*rawHandle);
    art::Handle< std::vector<crt::CRTTzero> > rawHandletzero;
    evt.getByLabel(fDataLabelTZeros, rawHandletzero); //what is the product instance name? no BernZMQ
    
    //check to make sure the data we asked for is valid                                                                            
    if(!rawHandletzero.isValid()){
      mf::LogWarning("CRTTrackProducer")
        <<"No CRTTzeros from producer module "<<fDataLabelTZeros;
      return;
    }
    
    std::vector<art::Ptr<crt::CRTTzero> > tzerolist;
    if (evt.getByLabel(fDataLabelTZeros,rawHandletzero))
      art::fill_ptr_vector(tzerolist, rawHandletzero);
 
    art::FindManyP<crt::CRTHit> fmht(rawHandletzero, evt, fDataLabelTZeros);
 
    //loop over tzeros
    for(size_t tzIter = 0; tzIter < tzerolist.size(); ++tzIter){   
      
      //count planes with hits for this tzero
      int np =0 ; //int ipflag[7] = {}; // CHANGED FROM 4 TO 7
      int tothits =0;
      for (int ip=0;ip<7;++ip) { // CHANGED FROM 4 TO 7
        if (tzerolist[tzIter]->nhits[ip]>0)  { np++; /*ipflag[ip]=1;*/ tothits+=tzerolist[tzIter]->nhits[ip];}  
      }
 
      if (np<2) continue;
      std::vector<art::Ptr<crt::CRTHit> > hitlist=fmht.at(tzIter);
      //for(size_t hit_i = 0; hit_i < hitlist.size(); hit
      if (fTrackMethodType==1) {
        double time_s_A = hitlist[0]->ts0_s;
 
        // find pairs of hits in different planes
        for (size_t ah = 0; ah< hitlist.size()-1; ++ah){        
          crt::CRTHit temphit=*hitlist[ah];
          CRTavehit Ahit = copyme(temphit);
          int planeA = hitlist[ah]->plane;
 
          for (size_t bh = ah+1; bh< hitlist.size(); ++bh){        
            int planeB = hitlist[bh]->plane;
 
            if (planeB!=planeA && !((planeA==3&&planeB==4)||(planeA==4&&planeB==3))) {  // make a track               
              temphit=*hitlist[bh];
              CRTavehit Bhit = copyme(temphit);
              crt::CRTTrack CRTcanTrack=shcut(Ahit,Bhit,time_s_A,0);
              CRTTrackCol->emplace_back(CRTcanTrack);
            }
 
          }
 
        }
 
      }
 
      else if ((fTrackMethodType==2) || (fTrackMethodType==3 && np==2 && tothits==2)) {        
 
        //loop over hits and get average x,y,z,pe for each plane CHANGED FROM 4 TO 7
        std::vector<float> thittime0[7];
        std::vector<float> thittime1[7];
        std::vector<float> tx[7];
        std::vector<float> ty[7];
        std::vector<float> tz[7];
        std::vector<float> pe[7];
        
        double time_s_A = hitlist[0]->ts0_s;
        //      double time_s_err = hitlist[0]->ts0_s_err;
        double time_s_err = 0.;
        double time1_ns_A = hitlist[0]->ts1_ns;
        double time0_ns_A = hitlist[0]->ts0_ns;
       
        //loop over hits for this tzero, sort by plane
        for (size_t ah = 0; ah< hitlist.size(); ++ah){        
          int ip = hitlist[ah]->plane;       
          thittime0[ip].push_back(hitlist[ah]->ts0_ns-time0_ns_A);
          thittime1[ip].push_back(hitlist[ah]->ts1_ns-time1_ns_A);
          tx[ip].push_back(hitlist[ah]->x_pos);
          ty[ip].push_back(hitlist[ah]->y_pos);
          tz[ip].push_back(hitlist[ah]->z_pos);
          pe[ip].push_back(hitlist[ah]->peshit);        
        } // loop over hits
 
        CRTavehit aveHits[7];
        //loop over planes and calculate average hits
        for (int ip = 0; ip < 7; ip++){
          if (tx[ip].size()>0){
            uint32_t at0; int32_t at1; uint16_t rt0,rt1;
            float totpe=0.0;
            float avet1=0.0; float rmst1 =0.0; 
            float avet0=0.0; float rmst0 =0.0; 
            float avex=0.0; float rmsx =0.0; 
            float avey=0.0; float rmsy =0.0; 
            float avez=0.0; float rmsz =0.0;
            vmanip(thittime0[ip],&avet0,&rmst0);
            vmanip(thittime1[ip],&avet1,&rmst1);
            at0 = (uint32_t)(avet0+time0_ns_A); rt0 = (uint16_t)rmst0;   
            at1 = (int32_t)(avet1+time1_ns_A); rt1 = (uint16_t)rmst1;
            vmanip(tx[ip],&avex,&rmsx);
            vmanip(ty[ip],&avey,&rmsy);
            vmanip(tz[ip],&avez,&rmsz);
            totpe=std::accumulate(pe[ip].begin(), pe[ip].end(), 0.0);
            CRTavehit aveHit = fillme(at0,rt0,at1,rt1,avex,rmsx,avey,rmsy,avez,rmsz,totpe,ip,"");
            aveHits[ip] = aveHit;
          }
          else {
            CRTavehit aveHit = fillme(0,0,0,0,-99999,-99999,-99999,-99999,-99999,-99999,-99999,ip,"");
            aveHits[ip] = aveHit;
          }
        }  
 
        // find pairs of hits in different planes
        for (size_t ah = 0; ah< 6; ++ah){        
          CRTavehit Ahit = aveHits[ah];
          if( Ahit.x_pos==-99999 ) continue;
 
          for (size_t bh = ah+1; bh< 7; ++bh){        
            if ( aveHits[bh].x_pos==-99999 ) continue;
 
            if (ah!=bh && !(ah==3&&bh==4)) {  // make a track               
              CRTavehit Bhit = aveHits[bh];
              crt::CRTTrack CRTcanTrack=shcut(Ahit,Bhit,time_s_A,time_s_err);
              CRTTrackCol->emplace_back(CRTcanTrack);
              nTrack++;
            }
 
          }
 
        }
 
      }
      
    }// loop over tzeros
  }
  
  //store track collection into event
  if(fStoreTrack == 1)
    evt.put(std::move(CRTTrackCol));

  mf::LogInfo("CRTTrackProducer")
    <<"Number of tracks            = "<<"\n"
    <<"Number of complete tracks   = "<<"\n"
    <<"Number of incomplete tracks = "<<nIncTrack;
  
} // CRTTrackProducer::produce()


void CRTTrackProducer::beginJob()
{

} // CRTTrackProducer::beginJob()


void CRTTrackProducer::endJob()
{

} // CRTTrackProducer::endJob()


// Function to calculate the mean and rms from a vector of values
void vmanip(std::vector<float> v, float* ave, float* rms)
{

  *ave=0.0; *rms =0.0;
  if (v.size()>0) {

    //  find the mean and *rms of all the vector elements
    double sum = std::accumulate(v.begin(), v.end(), 0.0);
    double mean = sum / v.size();
    *ave=mean;
    
    if (v.size()>1) {
    double sq_sum = std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
    double stdev = std::sqrt(sq_sum / v.size() - mean * mean);
    *rms=stdev;
    }

  }

} // vmanip()


// Function to make creating average crt hits easier
CRTavehit fillme(uint32_t ts0_ns, uint16_t ts0_ns_err, int32_t ts1_ns, uint16_t ts1_ns_err, 
                 float x_pos, float x_err, float y_pos, float y_err, float z_pos, float z_err, 
                 float pe, int plane, std::string tagger)
{

  CRTavehit h;

  h.ts0_ns     = ts0_ns;
  h.ts0_ns_err = ts0_ns_err;
  h.ts1_ns     = ts1_ns; 
  h.ts1_ns_err = ts1_ns_err;                                                        
  
  h.x_pos      = x_pos;
  h.x_err      = x_err;
  h.y_pos      = y_pos;
  h.y_err      = y_err;
  h.z_pos      = z_pos;
  h.z_err      = z_err;
  h.pe         = pe;
  h.plane      = plane;
  h.tagger     = tagger;

  return(h);

} // fillme()


// Function to copy average CRT hits
CRTavehit copyme(crt::CRTHit myhit)
{

  CRTavehit h;

  h.ts0_ns     = myhit.ts0_ns;
  h.ts0_ns_err = 0;
  h.ts1_ns     = myhit.ts1_ns;; 
  h.ts1_ns_err = 0;       
  h.x_pos      = myhit.x_pos;
  h.x_err      = myhit.x_err;
  h.y_pos      = myhit.y_pos;
  h.y_err      = myhit.y_err;
  h.z_pos      = myhit.z_pos;
  h.z_err      = myhit.z_err;
  h.pe         = myhit.peshit;
  h.plane      = myhit.plane;
  h.tagger     = myhit.tagger;

  return(h);

} // copyme()


// Function to make CRTTrack
crt::CRTTrack shcut(CRTavehit ppA, CRTavehit ppB, uint32_t time0s, uint16_t terr)
{

  crt::CRTTrack newtr;

  newtr.ts0_s         = time0s;
  newtr.ts0_s_err     = terr;
  newtr.ts0_ns_h1     = ppA.ts0_ns;
  newtr.ts0_ns_err_h1 = ppA.ts0_ns_err;
  newtr.ts0_ns_h2     = ppB.ts0_ns;
  newtr.ts0_ns_err_h2 = ppB.ts0_ns_err;
  newtr.ts0_ns        = (uint32_t)(0.5*(ppA.ts0_ns+ppB.ts0_ns));
  newtr.ts0_ns_err    = (uint16_t)(0.5*sqrt(ppA.ts0_ns_err*ppA.ts0_ns_err+ppB.ts0_ns_err*ppB.ts0_ns_err));
  newtr.ts1_ns        = (int32_t)(0.5*(ppA.ts1_ns+ppB.ts1_ns));
  newtr.ts1_ns_err    = (uint16_t)(0.5*sqrt(ppA.ts0_ns_err*ppA.ts0_ns_err+ppB.ts0_ns_err*ppB.ts0_ns_err));
  newtr.peshit        = ppA.pe+ppB.pe;
  newtr.x1_pos        = ppA.x_pos;
  newtr.x1_err        = ppA.x_err;
  newtr.y1_pos        = ppA.y_pos;
  newtr.y1_err        = ppA.y_err;
  newtr.z1_pos        = ppA.z_pos;
  newtr.z1_err        = ppA.z_err;
  newtr.x2_pos        = ppB.x_pos;
  newtr.x2_err        = ppB.x_err;
  newtr.y2_pos        = ppB.y_pos;
  newtr.y2_err        = ppB.y_err;
  newtr.z2_pos        = ppB.z_pos;
  newtr.z2_err        = ppB.z_err;
  float deltax        = ppA.x_pos-ppB.x_pos;
  float deltay        = ppA.y_pos-ppB.y_pos;
  float deltaz        = ppA.z_pos-ppB.z_pos;
  newtr.length        = sqrt(deltax*deltax+deltay*deltay+deltaz*deltaz);
  newtr.thetaxy       = atan2(deltax,deltay);
  newtr.phizy         = atan2(deltaz,deltay);
  newtr.plane1        = ppA.plane;
  newtr.plane2        = ppB.plane;

  return(newtr);

} // shcut()


// Function to make creating CRTHits easier
crt::CRTHit CRTTrackProducer::FillCrtHit(std::vector<uint8_t> tfeb_id, std::map<uint8_t, 
                                         std::vector<std::pair<int,float>>> tpesmap, float peshit, double time, 
                                         double x, double ex, double y, double ey, double z, double ez, std::string tagger)
{

    crt::CRTHit crtHit;

    crtHit.feb_id      = tfeb_id;
    crtHit.pesmap      = tpesmap;
    crtHit.peshit      = peshit;
    crtHit.ts0_s_corr  = 0;
    crtHit.ts0_ns      = time * 1e3;
    crtHit.ts0_ns_corr = 0;
    crtHit.ts1_ns      = time * 1e3;
    crtHit.ts0_s       = time * 1e-6;
    crtHit.x_pos       = x;
    crtHit.x_err       = ex;
    crtHit.y_pos       = y; 
    crtHit.y_err       = ey;
    crtHit.z_pos       = z;
    crtHit.z_err       = ez;
    crtHit.tagger      = tagger;

    return crtHit;

} // CRTTrackProducer::FillCrtHit()


// Function to make creating CRTTracks easier
crt::CRTTrack CRTTrackProducer::FillCrtTrack(crt::CRTHit hit1, crt::CRTHit hit2, bool complete)
{

  crt::CRTTrack newtr;

  newtr.ts0_s         = (hit1.ts0_s + hit2.ts0_s)/2.;
  newtr.ts0_s_err     = (uint32_t)((hit1.ts0_s - hit2.ts0_s)/2.);
  newtr.ts0_ns_h1     = hit1.ts0_ns;
  newtr.ts0_ns_err_h1 = hit1.ts0_ns_corr;
  newtr.ts0_ns_h2     = hit2.ts0_ns;
  newtr.ts0_ns_err_h2 = hit2.ts0_ns_corr;
  newtr.ts0_ns        = (uint32_t)((hit1.ts0_ns + hit2.ts0_ns)/2.);
  newtr.ts0_ns_err    = (uint16_t)(sqrt(hit1.ts0_ns_corr*hit1.ts0_ns_corr + hit2.ts0_ns_corr*hit2.ts0_ns_corr)/2.);
  newtr.ts1_ns        = (int32_t)(((double)(int)hit1.ts1_ns + (double)(int)hit2.ts1_ns)/2.);
  newtr.ts1_ns_err    = (uint16_t)(sqrt(hit1.ts0_ns_corr*hit1.ts0_ns_corr + hit2.ts0_ns_corr*hit2.ts0_ns_corr)/2.);
  newtr.peshit        = hit1.peshit+hit2.peshit;
  newtr.x1_pos        = hit1.x_pos;
  newtr.x1_err        = hit1.x_err;
  newtr.y1_pos        = hit1.y_pos;
  newtr.y1_err        = hit1.y_err;
  newtr.z1_pos        = hit1.z_pos;
  newtr.z1_err        = hit1.z_err;
  newtr.x2_pos        = hit2.x_pos;
  newtr.x2_err        = hit2.x_err;
  newtr.y2_pos        = hit2.y_pos;
  newtr.y2_err        = hit2.y_err;
  newtr.z2_pos        = hit2.z_pos;
  newtr.z2_err        = hit2.z_err;
  float deltax        = hit1.x_pos - hit2.x_pos;
  float deltay        = hit1.y_pos - hit2.y_pos;
  float deltaz        = hit1.z_pos - hit2.z_pos;
  newtr.length        = sqrt(deltax*deltax + deltay*deltay+deltaz*deltaz);
  newtr.thetaxy       = atan2(deltax,deltay);
  newtr.phizy         = atan2(deltaz,deltay);
  newtr.plane1        = hit1.plane;
  newtr.plane2        = hit2.plane;
  newtr.complete      = complete;

  return(newtr);

} // CRTTrackProducer::FillCrtTrack()


// Function to average hits within a certain distance of each other
std::vector<crt::CRTHit> CRTTrackProducer::AverageHits(std::vector<art::Ptr<crt::CRTHit>> hits)
{

  std::vector<crt::CRTHit> returnHits;
  std::vector<art::Ptr<crt::CRTHit>> aveHits;
  std::vector<art::Ptr<crt::CRTHit>> spareHits;

  if (hits.size()>0){
    // loop over size of tx
    bool first = true;
    TVector3 middle(0., 0., 0.);
    for (size_t i = 0; i < hits.size(); i++){
      // Get the position of the hit
      TVector3 pos(hits[i]->x_pos, hits[i]->y_pos, hits[i]->z_pos);
      // If first then set average = hit pos
      if(first){
        middle = pos;
        first = false;
      }
      // If distance from average < limit then add to average
      if((pos-middle).Mag() < fAverageHitDistance){
        aveHits.push_back(hits[i]);
      }
      // Else add to another vector
      else{
        spareHits.push_back(hits[i]);
      }
    }

    crt::CRTHit aveHit = DoAverage(aveHits);
    returnHits.push_back(aveHit);

    //Do this recursively
    std::vector<crt::CRTHit> moreHits = AverageHits(spareHits);
    returnHits.insert(returnHits.end(), moreHits.begin(), moreHits.end());
    return returnHits;
  }
  else {
    return returnHits;
  }

} // CRTTrackProducer::AverageHits()

  
// Take a list of hits and find average parameters
crt::CRTHit CRTTrackProducer::DoAverage(std::vector<art::Ptr<crt::CRTHit>> hits)
{

  // Initialize values
  std::string tagger = hits[0]->tagger;
  double xpos = 0.; 
  double ypos = 0.;
  double zpos = 0.;
  double xmax = -99999; double xmin = 99999;
  double ymax = -99999; double ymin = 99999;
  double zmax = -99999; double zmin = 99999;
  double ts1_ns = 0.;
  int nhits = 0;

  // Loop over hits
  for( auto& hit : hits ){
    // Get the mean x,y,z and times
    xpos += hit->x_pos;
    ypos += hit->y_pos;
    zpos += hit->z_pos;
    ts1_ns += (double)(int)hit->ts1_ns;
    // For the errors get the maximum limits
    if(hit->x_pos + hit->x_err > xmax) xmax = hit->x_pos + hit->x_err;
    if(hit->x_pos - hit->x_err < xmin) xmin = hit->x_pos - hit->x_err;
    if(hit->y_pos + hit->y_err > ymax) ymax = hit->y_pos + hit->y_err;
    if(hit->y_pos - hit->y_err < ymin) ymin = hit->y_pos - hit->y_err;
    if(hit->z_pos + hit->z_err > zmax) zmax = hit->z_pos + hit->z_err;
    if(hit->z_pos - hit->z_err < zmin) zmin = hit->z_pos - hit->z_err;
    // Add all the unique IDs in the vector
    nhits++;
  }

  // Create a hit
  crt::CRTHit crtHit = FillCrtHit(hits[0]->feb_id, hits[0]->pesmap, hits[0]->peshit, 
                                  (ts1_ns/nhits)*1e-3, xpos/nhits, (xmax-xmin)/2,
                                  ypos/nhits, (ymax-ymin)/2., zpos/nhits, (zmax-zmin)/2., tagger);

  return crtHit;

} // CRTTrackProducer::DoAverage()


// Function to create tracks from tzero hit collections
std::vector<crt::CRTTrack> CRTTrackProducer::CreateTracks(std::vector<crt::CRTHit> hits)
{

  std::vector<crt::CRTTrack> returnTracks;

  //Store list of hit pairs with distance between them
  std::vector<std::pair<std::pair<size_t, size_t>, double>> hitPairDist;
  std::vector<std::pair<size_t, size_t>> usedPairs;

  //Calculate the distance between all hits on different planes
  for(size_t i = 0; i < hits.size(); i++){
    crt::CRTHit hit1 = hits[i];

    for(size_t j = 0; j < hits.size(); j++){
      crt::CRTHit hit2 = hits[j];
      std::pair<size_t, size_t> hitPair = std::make_pair(i, j);
      std::pair<size_t, size_t> rhitPair = std::make_pair(j, i);

      //Only compare hits on different taggers and don't reuse hits
      if(hit1.tagger!=hit2.tagger && std::find(usedPairs.begin(), usedPairs.end(), rhitPair)==usedPairs.end()){
        //Calculate the distance between hits and store
        TVector3 pos1(hit1.x_pos, hit1.y_pos, hit1.z_pos);
        TVector3 pos2(hit2.x_pos, hit2.y_pos, hit2.z_pos);
        double dist = (pos1 - pos2).Mag();
        usedPairs.push_back(hitPair);
        hitPairDist.push_back(std::make_pair(hitPair, dist));
      }

    }

  }

  //Sort map by distance
  std::sort(hitPairDist.begin(), hitPairDist.end(), [](auto& left, auto& right){
            return left.second > right.second;});

  //Store potential hit collections + distance along 1D hit
  std::vector<std::pair<std::vector<size_t>, double>> tracks;
  for(size_t i = 0; i < hitPairDist.size(); i++){
    size_t hit_i = hitPairDist[i].first.first;
    size_t hit_j = hitPairDist[i].first.second;
    //Make sure bottom plane hit is always hit_i
    if(hits[hit_j].tagger=="volTaggerBot_0") std::swap(hit_i, hit_j);
    crt::CRTHit ihit = hits[hit_i];
    crt::CRTHit jhit = hits[hit_j];

    //If the bottom plane hit is a 1D hit
    if(ihit.x_err>100. || ihit.z_err>100.){
      double facMax = 1;
      std::vector<size_t> nhitsMax;
      double minDist = 99999;

      //Loop over the length of the 1D hit
      for(int i = 0; i<21; i++){
        double fac = (i)/10.;
        std::vector<size_t> nhits;
        double totalDist = 0.;
        TVector3 start(ihit.x_pos-(1.-fac)*ihit.x_err, ihit.y_pos, ihit.z_pos-(1.-fac)*ihit.z_err);
        TVector3 end(jhit.x_pos, jhit.y_pos, jhit.z_pos);
        TVector3 diff = start - end;

        //Loop over the rest of the hits
        for(size_t k = 0; k < hits.size(); k++){
          if(k == hit_i || k == hit_j || hits[k].tagger == ihit.tagger || hits[k].tagger == jhit.tagger) continue;
          //Calculate the distance between the track crossing point and the true hit
          TVector3 mid(hits[k].x_pos, hits[k].y_pos, hits[k].z_pos);
          TVector3 cross = CrossPoint(hits[k], start, diff);
          double dist = (cross-mid).Mag();

          //If the distance is less than some limit add the hit to the track and record the distance
          if(dist < fDistanceLimit){
            nhits.push_back(k);
            totalDist += dist;
          }

        }

        //If the distance down the 1D hit means more hits are included and they are closer to the track record it
        if(nhits.size()>=nhitsMax.size() && totalDist/nhits.size() < minDist){
          nhitsMax = nhits;
          facMax = fac;
          minDist = totalDist/nhits.size();
        }

        nhits.clear();
      }

      //Record the track candidate
      std::vector<size_t> trackCand;
      trackCand.push_back(hit_i);
      trackCand.push_back(hit_j);
      trackCand.insert(trackCand.end(), nhitsMax.begin(), nhitsMax.end());
      tracks.push_back(std::make_pair(trackCand, facMax));
    }

    //If there is no 1D hit
    else{
      TVector3 start(ihit.x_pos, ihit.y_pos, ihit.z_pos);
      TVector3 end(jhit.x_pos, jhit.y_pos, jhit.z_pos);
      TVector3 diff = start - end;
      std::vector<size_t> trackCand;
      trackCand.push_back(hit_i);
      trackCand.push_back(hit_j);

      //Loop over all the other hits
      for(size_t k = 0; k < hits.size(); k++){
        if(k == hit_i || k == hit_j || hits[k].tagger == ihit.tagger || hits[k].tagger == jhit.tagger) continue;
        //Calculate distance to other hits not on the planes of the track hits
        TVector3 mid(hits[k].x_pos, hits[k].y_pos, hits[k].z_pos);
        TVector3 cross = CrossPoint(hits[k], start, diff);
        double dist = (cross-mid).Mag();

        //Record any within a certain distance
        if(dist < fDistanceLimit){
          trackCand.push_back(k);
        }

      }

      tracks.push_back(std::make_pair(trackCand, 1));
    }

  }

  //Sort track candidates by number of hits
  std::sort(tracks.begin(), tracks.end(), [](auto& left, auto& right){
            return left.first.size() > right.first.size();});

  //Record used hits
  std::vector<size_t> usedHits;
  //Loop over candidates
  for(auto& track : tracks){
    size_t hit_i = track.first[0];
    size_t hit_j = track.first[1];

    // Make sure the first hit is the top high tagger if there are only two hits
    if(hits[hit_j].tagger=="volTaggerTopHigh_0") std::swap(hit_i, hit_j);
    crt::CRTHit ihit = hits[track.first[0]];
    crt::CRTHit jhit = hits[track.first[1]];

    //Check no hits in track have been used
    bool used = false;

    //Loop over hits in track candidate
    for(size_t i = 0; i < track.first.size(); i++){
      //Check if any of the hits have been used
      if(std::find(usedHits.begin(), usedHits.end(), track.first[i]) != usedHits.end()) used=true;
    }

    //If any of the hits have already been used skip this track
    if(used) continue;
    ihit.x_pos -= (1.-track.second)*ihit.x_err;
    ihit.z_pos -= (1.-track.second)*ihit.z_err;
    //Create track
    crt::CRTTrack crtTrack = FillCrtTrack(ihit, jhit, true);

    //If only the top two planes are hit create an incomplete/stopping track
    if(track.first.size()==2 && ihit.tagger == "volTaggerTopHigh_0" && jhit.tagger == "volTaggerTopLow_0"){ 
      crtTrack.complete = false;
    }

    returnTracks.push_back(crtTrack);

    //Record which hits were used only if the track has more than two hits
    //If there are multiple 2 hit tracks there is no way to distinguish between them
    //TODO: Add charge matching for ambiguous cases
    for(size_t i = 0; i < track.first.size(); i++){
      if(track.first.size()>2) usedHits.push_back(track.first[i]);
    }

  }

  return returnTracks;

} // CRTTrackProducer::CreateTracks()


// Function to calculate the crossing point of a track and tagger
TVector3 CRTTrackProducer::CrossPoint(crt::CRTHit hit, TVector3 start, TVector3 diff)
{
  TVector3 cross;
  // Use the error to get the fixed coordinate of a tagger
  // FIXME: can this be done better?
  if(hit.x_err > 0.39 && hit.x_err < 0.41){
    double xc = hit.x_pos;
    TVector3 crossp(xc, 
                    ((xc - start.X()) / (diff.X()) * diff.Y()) + start.Y(), 
                    ((xc - start.X()) / (diff.X()) * diff.Z()) + start.Z());
    cross = crossp;
  }

  else if(hit.y_err > 0.39 && hit.y_err < 0.41){
    double yc = hit.y_pos;
    TVector3 crossp(((yc - start.Y()) / (diff.Y()) * diff.X()) + start.X(), 
                    yc, 
                    ((yc - start.Y()) / (diff.Y()) * diff.Z()) + start.Z());
    cross = crossp;
  }

  else if(hit.z_err > 0.39 && hit.z_err < 0.41){
    double zc = hit.z_pos;
    TVector3 crossp(((zc - start.Z()) / (diff.Z()) * diff.X()) + start.X(), 
                    ((zc - start.Z()) / (diff.Z()) * diff.Y()) + start.Y(), 
                    zc);
    cross = crossp;
  }

  return cross;

} // CRTTrackProducer::CrossPoint()

DEFINE_ART_MODULE(CRTTrackProducer)

}// namespace sbnd
