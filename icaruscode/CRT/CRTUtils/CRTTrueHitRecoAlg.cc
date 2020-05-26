#ifndef ICARUS_CRTTRUEHITRECOALG_CC
#define ICARUS_CRTTRUEHITRECOALG_CC

#include "icaruscode/CRT/CRTUtils/CRTTrueHitRecoAlg.h"

using namespace icarus::crt;

//----------------------------------------------------------------------
CRTTrueHitRecoAlg::CRTTrueHitRecoAlg(const Config& config){
  this->reconfigure(config);

  fGeometryService = lar::providerFrom<geo::Geometry>();
  fCrtutils = new CRTCommonUtils();
}

//---------------------------------------------------------------------
CRTTrueHitRecoAlg::CRTTrueHitRecoAlg(){
  fGeometryService = lar::providerFrom<geo::Geometry>();
  fCrtutils = new CRTCommonUtils();
}


CRTTrueHitRecoAlg::~CRTTrueHitRecoAlg(){}

//---------------------------------------------------------------------
void CRTTrueHitRecoAlg::reconfigure(const Config& config){
  fUseReadoutWindow = config.UseReadoutWindow();
  fEDepMin = config.EDepMin();
  fRollupUnusedIds = config.RollupUnusedIds();
  fGlobalT0Offset = config.GlobalT0Offset();
  return;
}


//------------------------------------------------------------------------------
vector<pair<CRTHit,vector<sim::AuxDetIDE>>> CRTTrueHitRecoAlg::CreateCRTHits(
       vector<art::Ptr<sim::AuxDetSimChannel>> adscList) 
{
    vector<pair<CRTHit,vector<sim::AuxDetIDE>>> hitCol;
    map<int,map<int,tagger>> trackTaggers; //trackID -> (moduleID -> tagger )
   
    //fill taggers
    //loop over AuxDetSimChannels 
    for(auto const& adsc : adscList) {

        const int    adID    = adsc->AuxDetID();
        const int    layerID = fCrtutils->GetLayerID(adsc);
        const char   type    = fCrtutils->GetAuxDetType(adID);
        const string region  = fCrtutils->GetAuxDetRegion(adID);
        const int    adsID   = adsc->AuxDetSensitiveID();

        //loop over AuxDetIDEs
        for(auto const& ide : adsc->AuxDetIDEs()) {

            //FIX ME: for now, ignoring negative IDs, impliment fRollupUnusedIds
            if(/*ide.trackID>-1 &&*/ ide.energyDeposited*1000 < fEDepMin)
                continue;
            if(!fRollupUnusedIds && ide.energyDeposited*1000 < fEDepMin)
                continue;
            tagger& tag = (trackTaggers[ide.trackID])[adID];
            tag.type = type;
            tag.region = region;
            tag.layerID.insert(layerID);
            tag.stripLayer[adsID] = layerID;
            tag.stripIDE[adsID] = ide;
        }//AuxDetIDEs
    }//AuxDetSimChannels

    map <int,tagger> mTaggers; //MINOS module ID -> tagger
    int nmisscd=0, nmisspair=0;

    //apply logic to form hits
    //loop over trackIDs
    for (auto const& trk : trackTaggers) {

        // loop over taggers: modID->strip hit info
        for (auto const& tag : trk.second) {
    
            vector<sim::AuxDetIDE> vide; //IDEs in hit
            //XYZTVector 
            TLorentzVector rHit(0.,0.,0.,0.); //average hit position
            vector<uint8_t> feb_id = {(fCrtutils->ADToMac( tag.first)).first};
            map<uint8_t,vector< pair<int,float> > > pesmap; 
            float peshit = 0.;
            double xerr=0., yerr=0., zerr = 0.;

	    // if c or d type module
            if (tag.second.type=='c' || tag.second.type=='d') {

                // if "X-Y" coincidence
                if (tag.second.layerID.size()>1) {
    
                    // loop over module strips map: stripID->pos 4-vec
                    for (auto const& ide : tag.second.stripIDE) {
                        rHit += fCrtutils->AvgIDEPoint(ide.second);
                        vide.push_back(ide.second);
                        peshit+=ide.second.energyDeposited*1000;
                        pesmap[feb_id[0]].push_back(std::make_pair(ide.first,ide.second.energyDeposited*1000));
                    }
    
                    rHit*=1.0/tag.second.stripIDE.size();                    
                    rHit.SetT(rHit.T()+fGlobalT0Offset); 

                    //hit position RMS
                    for ( auto const& ide : tag.second.stripIDE) {
                        //XYZTVector 
                        TLorentzVector point = fCrtutils->AvgIDEPoint(ide.second);
                        xerr += pow(point.X()-rHit.X(),2);
                        yerr += pow(point.Y()-rHit.Y(),2);
                        zerr += pow(point.Z()-rHit.Z(),2);
                    }
                    xerr = sqrt(xerr/(tag.second.stripIDE.size()-1));
                    yerr = sqrt(yerr/(tag.second.stripIDE.size()-1));
                    zerr = sqrt(zerr/(tag.second.stripIDE.size()-1));
    
                } //if coincidence
                else{
                    nmisscd++;
                    continue;
                }
   
                hitCol.push_back(std::make_pair( 
                            FillCrtHit(feb_id, pesmap, peshit, rHit.T(), rHit.T(), 0, rHit.X(), xerr, 
                              rHit.Y(), yerr, rHit.Z(), zerr, tag.second.region),
                            vide) ); 

            }//if c or d type
    
            if ( tag.second.type=='m' ) {
                mTaggers[tag.first] = tag.second;
            }
    
        } //loop over taggers

        set <int> mPairs; //keep track of used auxdetIDs
        for (auto const& tag : mTaggers) {

           if (mPairs.find(tag.first) != mPairs.end()) continue; //don't double count   

            vector<sim::AuxDetIDE> vide; //IDEs in hit
            //XYZTVector 
            TLorentzVector rHit(0.,0.,0.,0.); //average hit position
            auto macpair = fCrtutils->ADToMac(tag.first);
            uint8_t mac1 = macpair.first;
            uint8_t mac11 = macpair.second;
            vector<uint8_t> feb_id = {mac1};
            if(mac1!=mac11) feb_id.push_back(mac11);
            map<uint8_t,vector< pair<int,float> > > pesmap;
            float peshit = 0.; 
            bool pairFound = false;
            double xerr=0., yerr=0., zerr = 0.;   
            set<int> layers;
            layers.insert(*(tag.second.layerID.begin()));
            mPairs.insert(tag.first);

            for (auto const& ide : tag.second.stripIDE) {
                rHit += fCrtutils->AvgIDEPoint(ide.second);
                vide.push_back(ide.second);
                peshit+=ide.second.energyDeposited*1000;
                pesmap[mac1].push_back(std::make_pair(ide.first,ide.second.energyDeposited*1000));
            }

            //inner loop over mTaggers (try to find a coincidence match) 
            for (auto const& tag2 : mTaggers) {
                if ( tag.first == tag2.first ) continue; //not the same module
                if ( mPairs.find(tag2.first) != mPairs.end()) continue; //not aleady counted
                if ( tag.second.region != tag2.second.region ) continue; //modules in same region

                layers.insert(*(tag2.second.layerID.begin()));
    
                //mark module as counted
                mPairs.insert(tag2.first);

                //mac5's for 2nd module in pair
                macpair = fCrtutils->ADToMac(tag2.first);
                uint8_t mac2 = macpair.first;
                uint8_t mac22 = macpair.second;
                feb_id.push_back(mac2);
                if(mac2!=mac22) feb_id.push_back(mac22);    

                // loop over module strips map: stripID->pos 4-vec
                for (auto const& ide : tag2.second.stripIDE) {
                    rHit += fCrtutils->AvgIDEPoint(ide.second);
                    vide.push_back(ide.second);
                    peshit+=ide.second.energyDeposited*1000;
                    pesmap[mac2].push_back(std::make_pair(ide.first,ide.second.energyDeposited*1000));
                }//for xyzt in second tagger (module)
 
                if(layers.size()==2)
                    pairFound = true;
            }//inner loop over taggers
    
            if (pairFound) {

                rHit*=1.0/vide.size();
                rHit.SetT(rHit.T()+fGlobalT0Offset);

                //hit position RMS
                for ( auto const& ide : vide) {
                    //XYZTVector 
                    TLorentzVector point = fCrtutils->AvgIDEPoint(ide);
                    xerr += pow(point.X()-rHit.X(),2);
                    yerr += pow(point.Y()-rHit.Y(),2);
                    zerr += pow(point.Z()-rHit.Z(),2);
                }

                xerr = sqrt(xerr/(vide.size()-1));
                yerr = sqrt(yerr/(vide.size()-1));
                zerr = sqrt(zerr/(vide.size()-1));

                hitCol.push_back(std::make_pair(
                            FillCrtHit(feb_id, pesmap, peshit, rHit.T(), rHit.T(), 0, rHit.X(), xerr,
                              rHit.Y(), yerr, rHit.Z(), zerr, tag.second.region),
                            vide) );

            } 
            else nmisspair++;
    
        } // outer loop over minos taggers
    } //loop over trackTaggers
 
    std::cout << "CRTTrueHitRecoAlg: nmisscd=" << nmisscd << ", nmissm=" << nmisspair << std::endl;
   
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
