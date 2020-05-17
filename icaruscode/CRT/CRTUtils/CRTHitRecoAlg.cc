#include "CRTHitRecoAlg.h"

using namespace icarus::crt;

//----------------------------------------------------------------------
CRTHitRecoAlg::CRTHitRecoAlg(const Config& config){
    this->reconfigure(config);
  
    fGeometryService = lar::providerFrom<geo::Geometry>();
    fDetectorClocks = lar::providerFrom<detinfo::DetectorClocksService>();
    fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
    fTrigClock = fDetectorClocks->TriggerClock();
    fFebMap = CRTCommonUtils::GetFebMap();
}

//---------------------------------------------------------------------
CRTHitRecoAlg::CRTHitRecoAlg(){
    fGeometryService = lar::providerFrom<geo::Geometry>();
    fDetectorClocks = lar::providerFrom<detinfo::DetectorClocksService>();
    fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
    fTrigClock = fDetectorClocks->TriggerClock();
    fFebMap = crt::CRTCommonUtils::GetFebMap();
}


CRTHitRecoAlg::~CRTHitRecoAlg(){}

//---------------------------------------------------------------------
void CRTHitRecoAlg::reconfigure(const Config& config){
    fVerbose = config.Verbose();
    fUseReadoutWindow = config.UseReadoutWindow(); 
    fQPed = config.QPed();
    fQSlope = config.QSlope();
    fPropDelay = config.PropDelay();
    fPEThresh = config.PEThresh();
    fCoinWindow = config.CoinWindow();
    return;
}
//---------------------------------------------------------------------------------------
vector<pair<CRTHit, vector<int>>> CRTHitRecoAlg::CreateCRTHits(vector<art::Ptr<CRTData>> crtList) {

    vector<pair<CRTHit, vector<int>>> returnHits;
    vector<int> dataIds;
  
    uint16_t nMissC = 0, nMissD = 0, nMissM = 0, nHitC = 0, nHitD = 0, nHitM = 0;
    if (fVerbose) mf::LogInfo("CRT") << "Found " << crtList.size() << " FEB events" << '\n';
  
    map<int,int> regCounts;
    std::set<int> regs;
    int febdat_last = -1;
  
    //loop over time-ordered CRTData
    for (size_t febdat_i=0; febdat_i<crtList.size(); febdat_i++) {
  
        uint8_t mac = crtList[febdat_i]->fMac5;
        char type = CRTCommonUtils::MacToType(mac);
        int region = CRTCommonUtils::MacToRegion(mac);
        CRTHit hit;
 
        std::cout << "found data with mac5 = " << (int)mac << ", " << string(1,type) 
                  << " type, region " << region << std::endl;
 
        dataIds.clear();
  
        if ((regs.insert(region)).second) regCounts[region] = 1;
        else regCounts[region]++;
  
        //CERN modules (intramodule coincidence)
        if ( type == 'c' ) {
            hit = MakeTopHit(crtList[febdat_i]);
            if(IsEmptyHit(hit))
                nMissC++;
            else {
                std::cout << "CERN hit produced" << std::endl;
                dataIds.push_back(febdat_i);
                returnHits.push_back(std::make_pair(hit,dataIds));
                nHitC++;
            }
        }
  
        //DC modules (intramodule coincidence)
        if ( type == 'd' ) {
            hit = MakeBottomHit(crtList[febdat_i]);
            if(IsEmptyHit(hit))
                nMissD++;
            else {
                std::cout << "DC hit produced" << std::endl;
                dataIds.push_back(febdat_i);
                returnHits.push_back(std::make_pair(hit,dataIds));
                nHitD++;
            }
        }
 
        //MINOS modules (intermodule coincidence) 
        if ( type == 'm' && (int)febdat_i>febdat_last) {
  
            vector<art::Ptr<CRTData>> coinData;
  
            for (size_t febdat_j=febdat_i; febdat_j<crtList.size(); febdat_j++) {
                std::cout << "i: " << febdat_j << ", region: " << 
                          CRTCommonUtils::MacToRegion(crtList[febdat_j]->fMac5) 
                          << ", ts0: " << crtList[febdat_j]->fTs0 << std::endl;
                if(CRTCommonUtils::MacToRegion(crtList[febdat_j]->fMac5)!=region) //not same region
                    continue;
                if(crtList[febdat_j]->fTs0 > crtList[febdat_i]->fTs0 + fCoinWindow) { //out of coinWindow

                    std::cout << "attempting to produce MINOS hit from " << coinData.size() 
                              << " data products..." << std::endl;
                    hit = MakeSideHit(coinData);
                    if(IsEmptyHit(hit))
                        nMissM++;
                    else {
                        std::cout << "MINOS hit produced" << std::endl;
                        returnHits.push_back(std::make_pair(hit,dataIds));
                        nHitM++;
                    }
                       
                    break;
                }
  
                febdat_last = (int)febdat_j;
                coinData.push_back(crtList[febdat_j]);
                dataIds.push_back(febdat_j);
  
            }// for febdat         
        } //if m type
  
    }//loop over CRTData products
  
    if(fVerbose) {
          mf::LogInfo("CRT") << returnHits.size() << " CRT hits produced!" << '\n'
              << "  nHitC: " << nHitC << " , nHitD: " << nHitD << " , nHitM: " << nHitM << '\n'
              << "    " << nMissC << " CRT hits missed!" << '\n';
          std::map<int,int>::iterator cts = regCounts.begin();
          mf::LogInfo("CRT") << " CRT Hits by region" << '\n';
          while (cts != regCounts.end()) {
              mf::LogInfo("CRT") << "reg: " << (*cts).first << " , hits: " << (*cts).second << '\n';
              cts++;
          }
    }//if Verbose
  
    return returnHits;

}
//--------------------------------------------------------------------------------------------
// Function to make filling a CRTHit a bit faster
CRTHit CRTHitRecoAlg::FillCRTHit(vector<uint8_t> tfeb_id, map<uint8_t,vector<pair<int,float>>> tpesmap,
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

} // CRTHitRecoAlg::FillCRTHit()

//------------------------------------------------------------------------------------------
CRTHit CRTHitRecoAlg::MakeTopHit(art::Ptr<CRTData> data){

    uint8_t mac = data->fMac5;
    map< uint8_t, vector< pair<int,float> > > pesmap;
    string region = CRTCommonUtils::MacToRegionName(mac);
    int adid  = CRTCommonUtils::MacToAuxDetID(mac,0); //module ID
    auto const& adGeo = fGeometryService->AuxDet(adid); //module

    double hitpoint[3], hitpointerr[3], hitlocal[3];
    TVector3 hitpos (0.,0.,0.);
    float petot = 0.;
    int adsid_max = -1, nabove=0;
    uint16_t adc_max = 0;
    TVector3 postrig;

    for(int chan=0; chan<32; chan++) {

        if(data->fAdc[chan]<=fPEThresh) continue;
        nabove++;
        int adsid = CRTCommonUtils::ChannelToAuxDetSensitiveID(mac,chan);
        petot += data->fAdc[chan];
        pesmap[mac].push_back(std::make_pair(chan,data->fAdc[chan]));

        TVector3 postmp = CRTCommonUtils::ChanToLocalCoords(fGeometryService,mac,chan);
        //strip along z-direction
        if(adsid < 8){
            hitpos.SetX(data->fAdc[chan]*postmp.X()+hitpos.X());
        }
        //strip along x-direction
        else {
            hitpos.SetZ(data->fAdc[chan]*postmp.Z()+hitpos.Z());
        }
        //identify trigger channel
        if(data->fAdc[chan]>adc_max) {
            adsid_max = chan;
            adc_max = data->fAdc[chan];
            postrig = postmp;
        }
    }

    //no channels above threshold? return empty hit
    if(nabove==0)
        return FillCRTHit({},{},0,0,0,0,0,0,0,0,0,0,"");

    hitpos*=1.0/petot; //hit position weighted by deposited charge
    hitlocal[0] = hitpos.X();
    hitlocal[1] = 0.;
    hitlocal[2] = hitpos.Z();
    
    auto const& adsGeo = adGeo.SensitiveVolume(adsid_max); //trigger strip
    double thit = data->fTs0;
    if(adsid_max<8)
        thit -= (adsGeo.HalfLength() - hitpos.Z() )*fPropDelay;
    else
        thit -= (adsGeo.HalfLength() - hitpos.X() )*fPropDelay;

    adGeo.LocalToWorld(hitlocal,hitpoint); //tranform from module to world coords

    hitpointerr[0] = adsGeo.HalfWidth1()*2/sqrt(12);
    hitpointerr[1] = adGeo.HalfHeight();
    hitpointerr[2] = adsGeo.HalfWidth1()*2/sqrt(12);

    CRTHit hit = FillCRTHit({mac},pesmap,petot,thit,thit,0,hitpoint[0],hitpointerr[0],
                            hitpoint[1],hitpointerr[1],hitpoint[2],hitpointerr[2],region);

    return hit;

} // CRTHitRecoAlg::MakeTopHit

//------------------------------------------------------------------------------------------
CRTHit CRTHitRecoAlg::MakeBottomHit(art::Ptr<CRTData> data){

    uint8_t mac = data->fMac5;
    map< uint8_t, vector< pair<int,float> > > pesmap;
    string region = CRTCommonUtils::MacToRegionName(mac);
    int adid  = CRTCommonUtils::MacToAuxDetID(mac,0); //module ID
    auto const& adGeo = fGeometryService->AuxDet(adid); //module

    double hitpoint[3], hitpointerr[3], hitlocal[3];
    TVector3 hitpos (0.,0.,0.);
    float petot = 0.;
    int adsid_max = -1, nabove=0;
    uint16_t adc_max = 0;
    TVector3 postrig;
    double xmin=0.,xmax=0.;

    for(int chan=0; chan<64; chan++) {

        if(data->fAdc[chan]<=fPEThresh) continue;
        nabove++;
        int adsid = CRTCommonUtils::ChannelToAuxDetSensitiveID(mac,chan);
        petot += data->fAdc[chan];
        pesmap[mac].push_back(std::make_pair(chan,data->fAdc[chan]));

        TVector3 postmp = CRTCommonUtils::ChanToLocalCoords(fGeometryService,mac,chan);
        //all strips along z-direction
        hitpos.SetX(data->fAdc[chan]*postmp.X()+hitpos.X());
        if(postmp.X()<xmin)
            xmin = postmp.X();
        if(postmp.X()>xmax)
            xmax = postmp.X();

        //identify trigger channel
        if(data->fAdc[chan]>adc_max) {
            adsid_max = adsid;
            adc_max = data->fAdc[chan];
            postrig = postmp;
        }
    }

    //no channels above threshold? return empty hit
    if(nabove==0)
        return FillCRTHit({},{},0,0,0,0,0,0,0,0,0,0,"");

    hitpos*=1.0/petot; //hit position weighted by deposited charge
    hitlocal[0] = hitpos.X();
    hitlocal[1] = 0.;
    hitlocal[2] = 0;

    auto const& adsGeo = adGeo.SensitiveVolume(adsid_max); //trigger strip
    double thit = data->fTs0 - adsGeo.HalfLength()*fPropDelay;

    adGeo.LocalToWorld(hitlocal,hitpoint); //tranform from module to world coords

    hitpointerr[0] = (xmax-xmin+2*adsGeo.HalfWidth1()*2)/sqrt(12);
    hitpointerr[1] = adGeo.HalfHeight();
    hitpointerr[2] = adsGeo.Length()/sqrt(12);

    CRTHit hit = FillCRTHit({mac},pesmap,petot,thit,thit,0,hitpoint[0],hitpointerr[0],
                            hitpoint[1],hitpointerr[1],hitpoint[2],hitpointerr[2],region);

    return hit;

} // CRTHitRecoAlg::MakeBottomHit

//-----------------------------------------------------------------------------------
CRTHit CRTHitRecoAlg::MakeSideHit(vector<art::Ptr<CRTData>> coinData) {

    vector<uint8_t> macs;
    map< uint8_t, vector< pair<int,float> > > pesmap;
    string region = CRTCommonUtils::MacToRegionName(coinData[0]->fMac5);
    int adid  = CRTCommonUtils::MacToAuxDetID(coinData[0]->fMac5,0); //module ID
    auto const& adGeo = fGeometryService->AuxDet(adid); //module

    double hitpoint[3], hitpointerr[3];
    TVector3 hitpos (0.,0.,0.);
    float petot = 0.;
    int adsid_max = -1, nabove=0;
    uint16_t adc_max = 0;
    TVector3 postrig;
    vector<double> ttrigs;
    double zmin=DBL_MAX, zmax = -DBL_MAX;
    double ymin=DBL_MAX, ymax = -DBL_MAX;
    double xmin=DBL_MAX, xmax = -DBL_MAX;
    std::set<int> layID;

    std::cout << "makeing MINOS hit...looping over coinData..." << std::endl;
    //loop over FEBs
    for(auto const& data : coinData) {

        macs.push_back(data->fMac5);

        //loop over channels
        for(int chan=0; chan<32; chan++) {

            if(data->fAdc[chan]<=fPEThresh) continue;
            nabove++;

            int adsid = CRTCommonUtils::ChannelToAuxDetSensitiveID(macs.back(),chan);
            petot += data->fAdc[chan];
            pesmap[macs.back()].push_back(std::make_pair(chan,data->fAdc[chan]));

            //inner or outer layer
            int layer = CRTCommonUtils::GetMINOSLayerID(fGeometryService,adGeo);
            layID.insert(layer);    
            TVector3 postmp = CRTCommonUtils::ChanToWorldCoords(fGeometryService,macs.back(),chan);

            //East/West Walls (all strips along z-direction) or
            // North/South inner walls (all strips along x-direction)
            if(!(region=="South" && layer==1)) {
                hitpos.SetY(data->fAdc[chan]*postmp.Y()+hitpos.Y());
                if(postmp.Y()<ymin)
                    ymin = postmp.Y();
                if(postmp.Y()>ymax)
                    ymax = postmp.Y();
                if(region!="South") {
                    hitpos.SetX(data->fAdc[chan]*postmp.X()+hitpos.X());
                    if(postmp.X()<xmin)
                        xmin = postmp.X();
                    if(postmp.X()>xmax)
                        xmax = postmp.X();
                }
            } 
            else { //else vertical strips in South wall
                hitpos.SetX(data->fAdc[chan]*postmp.X()+hitpos.X());
                if(postmp.X()<xmin)
                    xmin = postmp.X();
                if(postmp.X()>xmax)
                    xmax = postmp.X();
            }

            hitpos.SetZ(data->fAdc[chan]*postmp.Z()+hitpos.Z());
            if(postmp.X()<xmin)
                zmin = postmp.X();
            if(postmp.X()>xmax)
                zmax = postmp.X();

            //identify trigger channel
            if(data->fAdc[chan]>adc_max) {
                adsid_max = adsid;
                adc_max = data->fAdc[chan];
                postrig = postmp;
            }


        }//loop over channels

        auto const& adsGeo = adGeo.SensitiveVolume(adsid_max); //trigger strip
        ttrigs.push_back(data->fTs0 - adsGeo.HalfLength()*fPropDelay);

    }//loop over FEBs

    std::cout << "done...used " << nabove << " charge amplitudes" << std::endl;

    //no channels above threshold or no intermodule coincidences? return empty hit
    if(nabove==0 || layID.size()!=2)
        return FillCRTHit({},{},0,0,0,0,0,0,0,0,0,0,"");

    auto const& adsGeo = adGeo.SensitiveVolume(adsid_max);
    hitpos*=1.0/petot; //hit position weighted by deposited charge
    hitpoint[0] = hitpos.X();
    hitpoint[1] = hitpos.Y();
    hitpoint[2] = hitpos.Z();

    //error estimates (likely need to be revisted)
    if(region!="North" && region!="South"){
        hitpointerr[0] = (xmax-xmin)/sqrt(12);
        hitpointerr[1] = (ymax-ymin)/sqrt(12);
        hitpointerr[2] = adsGeo.Length()/sqrt(12);
    }

    if(region=="North"){
        hitpointerr[0] = (xmax-xmin)/sqrt(12);
        hitpointerr[1] = (ymax-ymin)/sqrt(12);
        hitpointerr[2] = (zmax-zmin)/sqrt(12);
    }

    if(region=="South"){
        hitpointerr[0] = adsGeo.HalfWidth1()*2/sqrt(12);
        hitpointerr[1] = adsGeo.HalfWidth1()*2/sqrt(12);
        hitpointerr[2] = (zmax-zmin)/sqrt(12);
    }

    //time stamp
    double thit = 0.;
    for(double const t : ttrigs)
        thit += t;
    thit*=1.0/coinData.size();

    std::cout << "generating CRTHit..." << std::endl;
    //generate hit
    CRTHit hit = FillCRTHit(macs,pesmap,petot,thit,thit,0,hitpoint[0],hitpointerr[0],
                            hitpoint[1],hitpointerr[1],hitpoint[2],hitpointerr[2],region);

    return hit;
}

//-----------------------------------------------------------------------------
bool CRTHitRecoAlg::IsEmptyHit(CRTHit hit) {

    if ( hit.feb_id.empty() && hit.pesmap.empty() && hit.peshit == 0
      && hit.ts0_ns == 0 && hit.ts1_ns == 0 && hit.plane == 0
      && hit.x_pos == 0 && hit.x_err == 0 && hit.y_pos == 0
      && hit.y_err == 0 && hit.z_pos == 0 && hit.z_err == 0 && hit.tagger == "") 
        return true;

    return false;
}
