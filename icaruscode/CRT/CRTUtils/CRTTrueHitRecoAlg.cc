#ifndef ICARUS_CRTTRUEHITRECOALG_CC
#define ICARUS_CRTTRUEHITRECOALG_CC

#include "icaruscode/CRT/CRTUtils/CRTTrueHitRecoAlg.h"

using namespace icarus::crt;

//----------------------------------------------------------------------
CRTTrueHitRecoAlg::CRTTrueHitRecoAlg(const Config& config){
  this->reconfigure(config);

  fGeometryService = lar::providerFrom<geo::Geometry>();
  fFebMap = CRTCommonUtils::GetFebMap();
}

//---------------------------------------------------------------------
CRTTrueHitRecoAlg::CRTTrueHitRecoAlg(){
  fGeometryService = lar::providerFrom<geo::Geometry>();
  fFebMap = CRTCommonUtils::GetFebMap();
}


CRTTrueHitRecoAlg::~CRTTrueHitRecoAlg(){}

//---------------------------------------------------------------------
void CRTTrueHitRecoAlg::reconfigure(const Config& config){
  fUseReadoutWindow = config.UseReadoutWindow();
  return;
}


//------------------------------------------------------------------------------
vector<pair<CRTHit,vector<sim::AuxDetIDE>>> CRTTrueHitRecoAlg::CreateCRTHits(
       vector<art::Ptr<sim::AuxDetSimChannel>> adscList) 
{
    vector<pair<CRTHit,vector<sim::AuxDetIDE>>> hitCol;
    map<int,map<int,tagger>> trackTaggers; //trackID -> (moduleID -> tagger )
    
    for(auto const& adsc : adscList) {

        auto const& adGeo = fGeometryService->AuxDet(adsc->AuxDetID());
        const int layerID = CRTCommonUtils::GetLayerID(fGeometryService,adsc);
        const char type = CRTCommonUtils::GetAuxDetType(adGeo);
        const int region = CRTCommonUtils::GetAuxDetRegionNum(CRTCommonUtils::GetAuxDetRegion(adGeo));
        const int adID = adsc->AuxDetID();
        const int adsID = adsc->AuxDetSensitiveID();

        for(auto const& ide : adsc->AuxDetIDEs()) {

            tagger& tag = (trackTaggers[ide.trackID])[adID];
            tag.type = type;
            tag.region = region;
            tag.layerID.insert(layerID);
            tag.stripLayer[adsID] = layerID;
            tag.stripIDE[adsID] = ide;
        }
    }

    map <int,tagger> mTaggers;

    //loop over trackIDs
    for (auto const& trk : trackTaggers) {

        // loop over taggers: modID->strip hit info
        for (auto const& tag : trk.second) {
    
            vector<sim::AuxDetIDE> vide; //IDEs in hit
            //XYZTVector 
            TLorentzVector rHit(0.,0.,0.,0.); //average hit position
            uint8_t feb_id = (CRTCommonUtils::ADToMac(fFebMap, tag.first)).first;
            map<uint8_t,vector< pair<int,float> > > pesmap;
            float peshit = 0.;
            double xerr=0., yerr=0., zerr = 0.;

	    // if c ord typrrmodule
            if (tag.second.type=='c' || tag.second.type=='d') {

                // if "X-Y" coincidence
                if (tag.second.layerID.size()>1) {
    
                    // loop over module strips map: stripID->pos 4-vec
                    for (auto const& ide : tag.second.stripIDE) {
                        rHit += CRTCommonUtils::AvgIDEPoint(ide.second);
                        vide.push_back(ide.second);
                        peshit+=ide.second.energyDeposited*1000;
                        pesmap[feb_id].push_back(std::make_pair(ide.first,ide.second.energyDeposited*1000));
                    }
    
                    rHit*=1.0/tag.second.stripIDE.size();                    
 
                    //hit position RMS
                    for ( auto const& ide : tag.second.stripIDE) {
                        //XYZTVector 
                        TLorentzVector point = CRTCommonUtils::AvgIDEPoint(ide.second);
                        xerr += pow(point.X()-rHit.X(),2);
                        yerr += pow(point.Y()-rHit.Y(),2);
                        zerr += pow(point.Z()-rHit.Z(),2);
                    }
                    xerr = sqrt(xerr/(tag.second.stripIDE.size()-1));
                    yerr = sqrt(yerr/(tag.second.stripIDE.size()-1));
                    zerr = sqrt(zerr/(tag.second.stripIDE.size()-1));
    
                } //if coincidence
                else continue;
   
                hitCol.push_back(std::make_pair( 
                            FillCrtHit({feb_id}, pesmap, peshit, rHit.T(), rHit.T(), 0, rHit.X(), xerr, 
                              rHit.Y(), yerr, rHit.Z(), zerr, CRTCommonUtils::GetRegionNameFromNum(tag.second.region)),
                            vide) ); 

            }//if c or d type
    
            if ( tag.second.type=='m' ) {
                mTaggers[tag.first] = tag.second;
            }
    
        } //loop over taggers

        set <int> mPairs;
        int nmisspair = 0;
        for (auto const& tag : mTaggers) {
   
            vector<sim::AuxDetIDE> vide; //IDEs in hit
            //XYZTVector 
            TLorentzVector rHit(0.,0.,0.,0.); //average hit position
            uint8_t mac1 = (CRTCommonUtils::ADToMac(fFebMap, tag.first)).first;
            vector<uint8_t> feb_id = {mac1};
            map<uint8_t,vector< pair<int,float> > > pesmap;
            float peshit = 0.; 
            bool pairFound = false;
            double xerr=0., yerr=0., zerr = 0.;   
 
            if (mPairs.find(tag.first) != mPairs.end()) continue; //don't double count
            for (auto const& tag2 : mTaggers) {
                if ( tag.first == tag2.first ) continue; //not the same module
                if ( mPairs.find(tag2.first) != mPairs.end()) continue; //not aleady counted
                if ( tag.second.region != tag2.second.region ) continue; //modules in same region
                if ( tag.second.layerID == tag2.second.layerID ) continue; //modules in opposite layers
    
                //mark modules as counted
                mPairs.insert(tag.first);
                mPairs.insert(tag2.first);

                uint8_t mac2 = (CRTCommonUtils::ADToMac(fFebMap, tag2.first)).first;
                feb_id.push_back(mac2);    

                // loop over module strips map: stripID->pos 4-vec
                if(pesmap.find(mac1)==pesmap.end()) 
                    for (auto const& ide : tag.second.stripIDE) {
                        rHit += CRTCommonUtils::AvgIDEPoint(ide.second);
                        vide.push_back(ide.second);
                        peshit+=ide.second.energyDeposited*1000;
                        pesmap[mac1].push_back(std::make_pair(ide.first,ide.second.energyDeposited*1000));
                    }
    
                for (auto const& ide : tag2.second.stripIDE) {
                    rHit += CRTCommonUtils::AvgIDEPoint(ide.second);
                    vide.push_back(ide.second);
                    peshit+=ide.second.energyDeposited*1000;
                    pesmap[mac2].push_back(std::make_pair(ide.first,ide.second.energyDeposited*1000)); 
                }//Hfor xyzt in second tagger (module)

                rHit*=1.0/(tag.second.stripIDE.size()+tag2.second.stripIDE.size());
    
                //hit position RMS
                for ( auto const& ide : tag.second.stripIDE) {
                    //XYZTVector 
                    TLorentzVector point = CRTCommonUtils::AvgIDEPoint(ide.second);
                    xerr += pow(point.X()-rHit.X(),2);
                    yerr += pow(point.Y()-rHit.Y(),2);
                    zerr += pow(point.Z()-rHit.Z(),2);
                }

                for ( auto const& ide : tag2.second.stripIDE) {
                    //XYZTVector 
                    TLorentzVector point = CRTCommonUtils::AvgIDEPoint(ide.second);
                    xerr += pow(point.X()-rHit.X(),2);
                    yerr += pow(point.Y()-rHit.Y(),2);
                    zerr += pow(point.Z()-rHit.Z(),2);
                }

                xerr = sqrt(xerr/(tag.second.stripIDE.size()+tag2.second.stripIDE.size()-1));
                yerr = sqrt(yerr/(tag.second.stripIDE.size()+tag2.second.stripIDE.size()-1));
                zerr = sqrt(zerr/(tag.second.stripIDE.size()+tag2.second.stripIDE.size()-1));   
 
                pairFound = true;
            }//inner loop over taggers
    
            if (pairFound) {
                hitCol.push_back(std::make_pair(
                            FillCrtHit({feb_id}, pesmap, peshit, rHit.T(), rHit.T(), 0, rHit.X(), xerr,
                              rHit.Y(), yerr, rHit.Z(), zerr, CRTCommonUtils::GetRegionNameFromNum(tag.second.region)),
                            vide) );

            } 
            else nmisspair++;
    
        } // outer loop over taggers
    } //loop over trackTaggers
    
    return hitCol;
}

//--------------------------------------------------------------------------------------------
// Function to make filling a CRTHit a bit faster
CRTHit CRTTrueHitRecoAlg::FillCrtHit(vector<uint8_t> tfeb_id, map<uint8_t,vector<pair<int,float>>> tpesmap, 
                            float peshit, double time0, double time1, int plane,
                            double x, double ex, double y, double ey, double z, double ez, string tagger){
    CRTHit crtHit;
    crtHit.feb_id      = tfeb_id;
    crtHit.pesmap      = tpesmap;
    crtHit.peshit      = peshit;
    crtHit.ts0_s_corr  = time0*1e-9;
    crtHit.ts0_ns      = time0;
    crtHit.ts0_ns_corr = time0;
    crtHit.ts1_ns      = time1;
    crtHit.ts0_s       = time0 * 1e-9;
    crtHit.plane       = plane;
    crtHit.x_pos       = x;
    crtHit.x_err       = ex;
    crtHit.y_pos       = y;
    crtHit.y_err       = ey;
    crtHit.z_pos       = z;
    crtHit.z_err       = ez;
    crtHit.tagger      = tagger;

    return crtHit;

} //FillCrtHit()

#endif
