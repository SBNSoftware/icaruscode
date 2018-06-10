////////////////////////////////////////////////////////////////////////
/// \file   CandHitICARUS.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include "larreco/HitFinder/HitFinderTools/ICandidateHitFinder.h"

#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"

#include <cmath>
#include <fstream>

namespace reco_tool
{

class CandHitICARUS : ICandidateHitFinder
{
public:
    explicit CandHitICARUS(const fhicl::ParameterSet& pset);
    
    ~CandHitICARUS();
    
    void configure(const fhicl::ParameterSet& pset) override;
    
    void findHitCandidates(const std::vector<float>&,
                           size_t,
                           double,
                           HitCandidateVec&) const override;
    
    void findHitCandidates(std::vector<float>::const_iterator,
                           std::vector<float>::const_iterator,
                           size_t,
                           double,
                           HitCandidateVec&) const override;
    
    void MergeHitCandidates(const Waveform&,
                            const HitCandidateVec&,
                            MergeHitCandidateVec&) const override;
    
private:
      void expandHit(HitCandidate& h, std::vector<float> holder, HitCandidateVec how);
    void prova() {return ;}
    
    int          fInd1Width;            //INITIAL WIDTH FOR INDUCTION HITFINDING.
    int              fInd2Width;            //INITIAL WIDTH FOR INDUCTION HITFINDING.
    int              fColWidth;            //INITIAL WIDTH FOR COLLECTION HITFINDING.
    unsigned int              fInd1Window;            //INITIAL WINDOW FOR INDUCTION HITFINDING.
    unsigned int              fInd2Window;            //INITIAL WINDOW FOR INDUCTION HITFINDING.
    unsigned int              fColWindow;            //INITIAL WINDOW FOR COLLECTION HITFINDING.
    int              fInd1Threshold;            //THRESHOLD FOR INDUCTION HITFINDING.
    int              fInd2Threshold;            //THRESHOLD FOR INDUCTION HITFINDING.
    int              fColThreshold;            //THRESHOLD FOR COLLECTION HITFINDING.
    int              fInd1Above;            //MINIMAL NUMBER OF TICKS ABOVE THRESHOLD FOR INDUCTION HITFINDING.
    int              fInd2Above;
    int              fColAbove;
    int              fInd1Fall;
    int              fInd2Fall;
    int              fColFall;
    
    const geo::GeometryCore*  fGeometry = lar::providerFrom<geo::Geometry>();
};
    
//----------------------------------------------------------------------
// Constructor.
CandHitICARUS::CandHitICARUS(const fhicl::ParameterSet& pset)
{
    configure(pset);
}
    
CandHitICARUS::~CandHitICARUS()
{
}
    
void CandHitICARUS::configure(const fhicl::ParameterSet& pset)
{
    // Start by recovering the parameters

    fInd1Width = pset.get< int  >("Ind1Width");
    fInd2Width = pset.get< int  >("Ind2Width");
    fColWidth = pset.get< int  >("ColWidth");
    fInd1Window = pset.get< unsigned int  >("Ind1Window");
    fInd2Window = pset.get< unsigned int  >("Ind2Window");
    fColWindow = pset.get< unsigned int  >("ColWindow");
    fInd1Threshold = pset.get< int  >("Ind1Threshold");
    fInd2Threshold = pset.get< int  >("Ind2Threshold");
    fColThreshold = pset.get< int  >("ColThreshold");
    fInd1Above = pset.get< int  >("Ind1Above");
    fInd2Above = pset.get< int  >("Ind2Above");
    fColAbove = pset.get< int  >("ColAbove");
    fInd1Fall = pset.get< int  >("Ind1Fall");
    fInd2Fall = pset.get< int  >("Ind2Fall");
    fColFall = pset.get< int  >("ColFall");
    
    return;
}
    
void CandHitICARUS::findHitCandidates(const std::vector<float>& waveform,
                                        size_t plane, double wire,
                                        HitCandidateVec& hits) const
{
    if(plane!=2) return;
   // if(wire>2639) return;
//for( int j=0;j<4096;j++)
  //  std::cout << " j " << " waveform " << waveform[j] << std::endl;
int iflag;
    int localbellow,rising;
    int begin;
    // int nSamp=digitVec->Samples();
    int lastcomputedmean,count;
    int peakheight;
    int localminidx,localmin;
    // int jj, jcount;
    // float area;
    HitCandidate h;
    //std::vector<HitCandidate> hits;
    unsigned int i;
    // Hit finding parameters
    const int rise=5;
    
    double threshold=0;
    unsigned int window=0;
    int abovecut(0),fall(0);
    unsigned int width(0);
    
    
    // initialize parameters
    iflag=0;                // equal to one if we are within a hit candidate
    peakheight=-9999;       // last found hit maximum
    begin=-1;               // last found hit initial sample
    localbellow=0;          // number of times we are bellow peakheight
    lastcomputedmean=0;     // the value we compare with to know if we have a hit
    count=0;
    for(unsigned int j=0;j<window;j++)
        count+=waveform[j];
    lastcomputedmean=(count>=0)? 0 : count/window;
    //      negpeakwidth=0;
    //      if(lastcomputedmean<0) negpeakwidth++;
    localmin=9999;
    localminidx=-1;
    //     if ( typ == BasicData::ISBASICPLANE_PMT && iw==1 ){
    //cout << "window, mean " <<  window <<" "<<lastcomputedmean<<endl;
    //cout << "wire, nsamp " <<  iw <<" "<<nSamp<<endl;
    //}
    // loop on the selected samples
    
            //std::cout << "candhitfinderICARUS " << std::endl;
    
    if(plane == 0) threshold=fInd1Threshold;
    else if(plane == 1) threshold=fInd2Threshold;
    else if(plane == 2) threshold=fColThreshold;
    if(plane == 0) window=fInd1Window;
    else if(plane == 1) window=fInd2Window;
    else if(plane == 2) window=fColWindow;
    if(plane == 0) abovecut=fInd1Above;
    else if(plane == 1) abovecut=fInd2Above;
    else if(plane == 2) abovecut=fColAbove;
    if(plane == 0) fall=fInd1Fall;
    else if(plane == 1) fall=fInd2Fall;
    else if(plane == 2) fall=fColFall;
    if(plane == 0) width=fInd1Width;
    else if(plane == 1) width=fInd2Width;
    else if(plane == 2) width=fColWidth;
    
    for( i=0;i<4096;i++)
    { //2
        //std::cout << " i " << i << " waveform " << waveform[i] << std::endl;
        // skip sharp peaks
        //if(abs(waveform[i]-lastcomputedmean)>100)
        //  continue;
        
        if(!iflag)
            lastcomputedmean=0;
        
        if(waveform[i]-lastcomputedmean>threshold) // we're within a hit OR hit group
        { //3
            //      if(iwire==4526&&plane==2&&cryostat==0&&tpc==0)
             //   std::cout << " over threshold bin " << i << std::endl;
            iflag=1;
            
            // we're in the beginning of the hit
            if(begin<0) {
                
                begin=i;     // hit starting point
            }
            
            // keep peak info
            if(waveform[i]-lastcomputedmean>peakheight)
            {
                peakheight=waveform[i]-lastcomputedmean;
     //           h.hitCenter=waveform[i];
                h.hitHeight=peakheight;
//                h.iWire=iw;
                h.hitCenter=i;
                std::cout << " hitcenter " << i << std::endl;
                localbellow=0;
            }
            
            // resolve close hits
            //          if(*(pf+i)-(peakheight+lastcomputedmean)<-threshold) // we're in the slope down
            if(waveform[i]-(peakheight+lastcomputedmean)<-1) // we're in the slope down
            {
                localbellow++;
                std::cout << " i " << i << " localbellow " << localbellow << " abovecut " << abovecut << std::endl;
            }
            if(localbellow>abovecut)
            { //4
                // keep local minimum as border between consecutive hits
                if(waveform[i]<localmin) {localmin=waveform[i];localminidx=i;}
                
                // count the number of rising samples within the following n=<rise> samples
                rising=0;
                for(int l=0;l<rise;l++)
                {
                    if(i+l<4096) {
                        if(waveform[i+l+1]-waveform[i+l]>0) rising++;
                        else if(waveform[i+l+1]-waveform[i+l]<0) rising--;
                    }
                }
                std::cout << " rising " << rising << " abovecut " << abovecut << std::endl;
                // if after a slope down there's a slope up save the previous hit
                if(rising>abovecut)
                { //5
                    h.startTick=begin;
                    //              h.finDrift=i+iniSamp;
                    h.stopTick=localminidx;

                    std::cout << "  saving previous hit " << begin << " fin " << localminidx << std::endl;
                    
                    if((h.stopTick-h.hitCenter)>=fall && h.stopTick-h.startTick>width)
                    { //6
                        h.minTick=h.startTick;
                        h.maxTick=h.stopTick;
                        // std::cout << " adding hit case 1" << std::endl;
                        h.hitSigma      = 0.5*(h.stopTick-h.startTick);

                        hits.push_back(h);
                        peakheight=-9999;
                        h.hitHeight=0;
                        //begin=i+iniSamp+1;
                        begin=localminidx+1;
                        localbellow=0;
                        localminidx=-1;
                        localmin=9999;
                    }
                }
            }
        }
        else // outside the hit
        { //3
            //if (  typ == BasicData::ISBASICPLANE_PMT  && iw==1){
            //    cout << "out hit " << iflag << "iadc " << h.iAdcMax <<endl;
            //    cout << " in, fin " << begin << " " << i+iniSamp << "peak "<<h.iDrift <<endl;}
            if (iflag==1 && h.hitHeight) //just getting out of the latest hit
            { //4
                h.startTick=begin;
                h.stopTick=i;

                if((h.stopTick-h.hitCenter)>=fall && (h.stopTick-h.startTick)>width)
                { //5
                    h.minTick=h.startTick;
                    h.maxTick=h.stopTick;
                    //if(iwire==4526&&plane==1&&cryostat==0&&tpc==0)
                    //std::cout << " adding hit case 2, tick " << i << std::endl;
                 //   if(cryostat==0&&tpc==0&&plane==2&&h.iWire==2656)
                   //     std::cout << "  before expand ini " << h.iniDrift << " fin " << h.finDrift << std::endl;
                //expandHit(h,waveform,hits);
                   
                   // if(cryostat==0&&tpc==0&&plane==2&&h.iWire==2656)
                     //   std::cout << "  after expand ini " << h.iniDrift << " fin " << h.finDrift << std::endl;
                    h.hitSigma      = 0.5*(h.stopTick-h.startTick);

                    hits.push_back(h);
                    //InsertHit(&h);
                    //  if ( typ != BasicData::ISBASICPLANE_PMT)
                    
                }
                
                peakheight=-9999;
                begin=-1;
                iflag=0;
                localbellow=0;
            }
            
        }
        
        //keep the last mean value to which we have to come back after the hit
        if(i>=window && window) {
            if(waveform[i]-count/window>-10)
            {
                count+=waveform[i];
                count-=waveform[i-window];
                //        if((float) count/window<0) negpeakwidth++;
                //        else negpeakwidth=0;
            }
        }
    } //end loop on samples
    
    //if we were within a hit while reaching last sample, keep it
    if(iflag==1 && h.hitHeight) //just getting out of the latest hit
    { //3
        h.startTick=begin;
        h.stopTick=i-1;
        
        if((h.stopTick-h.hitCenter)>=fall && (h.stopTick-h.startTick)>width)
        {
          //  if(cryostat==0&&tpc==0&&plane==2&&h.iWire==2656)
//                std::cout << "  before expand ini " << h.iniDrift << " fin " << h.finDrift << std::endl;
       // expandHit(h,waveform,hits);
          
  //          if(cryostat==0&&tpc==0&&plane==2&&h.iWire==2656)
    //            std::cout << "  after expand ini " << h.iniDrift << " fin " << h.finDrift << std::endl;
            h.minTick=h.startTick;
            h.maxTick=h.stopTick;
            //std::cout << " adding hit case 3 " << i << std::endl;

            h.hitSigma      = 0.5*(h.stopTick-h.startTick);
            hits.push_back(h);
            // InsertHit(&h);
        }
    }
    if(hits.size())
     std::cout << " wire " << wire << " found hits " << hits.size() << std::endl;
    
    
    
    return;
}
  
void CandHitICARUS::findHitCandidates(std::vector<float>::const_iterator startItr,
                                        std::vector<float>::const_iterator stopItr,
                                        size_t                             roiStartTick,
                                        double                             roiThreshold,
                                        HitCandidateVec&                   hitCandidateVec) const
{
    // Need a minimum number of ticks to do any work here
    if (std::distance(startItr,stopItr) > 4)
    {
        // Find the highest peak in the range given
        auto maxItr = std::max_element(startItr, stopItr);
        
        float maxValue = *maxItr;
        int   maxTime  = std::distance(startItr,maxItr);
        
        if (maxValue > roiThreshold)
        {
            // backwards to find first bin for this candidate hit
            auto firstItr = std::distance(startItr,maxItr) > 2 ? maxItr - 1 : startItr;
            
            while(firstItr != startItr)
            {
                // Check for pathology where waveform goes too negative
                if (*firstItr < -roiThreshold) break;
                
                // Check both sides of firstItr and look for min/inflection point
                if (*firstItr < *(firstItr+1) && *firstItr <= *(firstItr-1)) break;
                
                firstItr--;
            }
            
            int firstTime = std::distance(startItr,firstItr);
            
            // Recursive call to find all candidate hits earlier than this peak
            findHitCandidates(startItr, firstItr + 1, roiStartTick, roiThreshold, hitCandidateVec);
            
            // forwards to find last bin for this candidate hit
            auto lastItr = std::distance(maxItr,stopItr) > 2 ? maxItr + 1 : stopItr - 1;
            
            while(lastItr != stopItr - 1)
            {
                // Check for pathology where waveform goes too negative
                if (*lastItr < -roiThreshold) break;
                
                // Check both sides of firstItr and look for min/inflection point
                if (*lastItr <= *(lastItr+1) && *lastItr < *(lastItr-1)) break;
                
                lastItr++;
            }
            
            int lastTime = std::distance(startItr,lastItr);
            
            // Now save this candidate's start and max time info
            HitCandidate hitCandidate;
            hitCandidate.startTick     = roiStartTick + firstTime;
            hitCandidate.stopTick      = roiStartTick + lastTime;
            hitCandidate.maxTick       = roiStartTick + firstTime;
            hitCandidate.minTick       = roiStartTick + lastTime;
            hitCandidate.maxDerivative = *(startItr + firstTime);
            hitCandidate.minDerivative = *(startItr + lastTime);
            hitCandidate.hitCenter     = roiStartTick + maxTime;
            hitCandidate.hitSigma      = std::max(2.,float(lastTime - firstTime) / 6.);
            hitCandidate.hitHeight     = maxValue;
            
            hitCandidateVec.push_back(hitCandidate);
            
            // Recursive call to find all candidate hits later than this peak
            findHitCandidates(lastItr + 1, stopItr, roiStartTick + std::distance(startItr,lastItr + 1), roiThreshold, hitCandidateVec);
        }
    }
    
    return;
}
    void CandHitICARUS::MergeHitCandidates(const Waveform&        signalVec,
                                               const HitCandidateVec& hitCandidateVec,
                                               MergeHitCandidateVec&  mergedHitsVec) const
    {
        if(hitCandidateVec.size())
        mergedHitsVec.push_back(hitCandidateVec);
      //  std::cout << " mergedhitsvec size" << mergedHitsVec.size() << std::endl;
    }
    
    void CandHitICARUS::expandHit(HitCandidate& h, std::vector<float> holder, HitCandidateVec how)
    {
        // Given a hit or hit candidate <hit> expand its limits to the closest minima
        int nsamp=50;
        int cut=1;
        int upordown;
        int found=0;
        
        unsigned int first=h.startTick;
        unsigned int last =h.stopTick;
        
        std::vector<HitCandidate>::iterator hiter;
        std::vector<HitCandidate> hlist;
        
        // fill list of existing hits on this wire
        for(unsigned int j=0;j<how.size();j++)
        {
            HitCandidate h2=how[j];
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
               HitCandidate h2=*hiter;
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
                HitCandidate h2=*hiter;
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

DEFINE_ART_CLASS_TOOL(CandHitICARUS)
}
