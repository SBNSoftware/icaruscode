#include "CRTHitRecoAlg.h"

using namespace icarus::crt;

//----------------------------------------------------------------------
CRTHitRecoAlg::CRTHitRecoAlg(const Config& config){
    this->reconfigure(config);
  
    fGeometryService  = lar::providerFrom<geo::Geometry>();
    fDetectorClocks   = lar::providerFrom<detinfo::DetectorClocksService>();
    fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
    fTrigClock = fDetectorClocks->TriggerClock();
    fCrtutils = new CRTCommonUtils();
}

//---------------------------------------------------------------------
CRTHitRecoAlg::CRTHitRecoAlg(){
    fGeometryService = lar::providerFrom<geo::Geometry>();
    fDetectorClocks = lar::providerFrom<detinfo::DetectorClocksService>();
    fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
    fTrigClock = fDetectorClocks->TriggerClock();
    fCrtutils = new CRTCommonUtils();
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
vector<pair<sbn::crt::CRTHit, vector<int>>> CRTHitRecoAlg::CreateCRTHits(vector<art::Ptr<CRTData>> crtList) {

    vector<pair<CRTHit, vector<int>>> returnHits;
    vector<int> dataIds;
  
    uint16_t nMissC = 0, nMissD = 0, nMissM = 0, nHitC = 0, nHitD = 0, nHitM = 0;
    if (fVerbose) mf::LogInfo("CRT") << "Found " << crtList.size() << " FEB events" << '\n';
  
    map<string,int> regCounts;
    std::set<string> regs;
    map<string,vector<size_t>> sideRegionToIndices;

    //loop over time-ordered CRTData
    for (size_t febdat_i=0; febdat_i<crtList.size(); febdat_i++) {
  
        uint8_t mac = crtList[febdat_i]->fMac5;
        int adid  = fCrtutils->MacToAuxDetID(mac,0); //module ID
        string region = fCrtutils->GetAuxDetRegion(adid);
        char type = fCrtutils->GetAuxDetType(adid);
        CRTHit hit;

        //if(fVerbose) 
        //    std::cout << "found data with mac5 = " << (int)mac << ", " << string(1,type) 
        //          << " type, region " << region << std::endl;
 
        dataIds.clear();
  
        //CERN modules (intramodule coincidence)
        if ( type == 'c' ) {
            hit = MakeTopHit(crtList[febdat_i]);
            if(IsEmptyHit(hit))
                nMissC++;
            else {
                //if(fVerbose)
                //    std::cout << "CERN hit produced" << std::endl;
                dataIds.push_back(febdat_i);
                returnHits.push_back(std::make_pair(hit,dataIds));
                if ((regs.insert(region)).second) regCounts[region] = 1;
                else regCounts[region]++;

                nHitC++;
            }
        }
  
        //DC modules (intramodule coincidence)
        if ( type == 'd' ) {
            hit = MakeBottomHit(crtList[febdat_i]);
            if(IsEmptyHit(hit))
                nMissD++;
            else {
                //if(fVerbose)
                //    std::cout << "DC hit produced" << std::endl;
                dataIds.push_back(febdat_i);
                returnHits.push_back(std::make_pair(hit,dataIds));
                if ((regs.insert(region)).second) regCounts[region] = 1;
                else regCounts[region]++;

                nHitD++;
            }
        }
 
        if ( type == 'm' )
            sideRegionToIndices[region].push_back(febdat_i);

    }//loop over CRTData products

    vector<size_t> unusedDataIndex;
    for(auto const& regIndices : sideRegionToIndices) {

        if(fVerbose) 
            std::cout << "searching for side CRT hits in region, " << regIndices.first << std::endl;
        vector<size_t> indices = regIndices.second;

        for(size_t index_i=0; index_i < indices.size(); index_i++) {

            dataIds.clear();
            dataIds.push_back(indices[index_i]);
            vector<art::Ptr<CRTData>> coinData = {crtList[indices[index_i]]};

 
            //inner loop over data after data_i in time
            for (size_t index_j=index_i+1; index_j<indices.size(); index_j++) {

                if(crtList[indices[index_j]]->fTs0 < crtList[indices[index_i]]->fTs0)
                    mf::LogError("CRTHitRecoAlg::CreateCRTHits") <<
                        "bad time ordering!" << '\n';

                if(crtList[indices[index_j]]->fTs0 <= crtList[indices[index_i]]->fTs0 + fCoinWindow) {
                    coinData.push_back(crtList[indices[index_j]]);
                    dataIds.push_back(indices[index_j]);
                }

                //out of coinWindow
                if(crtList[indices[index_j]]->fTs0 > crtList[indices[index_i]]->fTs0 + fCoinWindow
                   || index_j==indices.size()-1) 
                { 

                    if(fVerbose)
                        std::cout << "attempting to produce MINOS hit from " << coinData.size() 
                              << " data products..." << std::endl;

                    CRTHit hit = MakeSideHit(coinData);

                    if(IsEmptyHit(hit)){
                        unusedDataIndex.push_back(indices[index_i]);
                        nMissM++;
                    }
                    else {
                        if(fVerbose)
                            std::cout << "MINOS hit produced" << std::endl;

                        returnHits.push_back(std::make_pair(hit,dataIds));

                        if ((regs.insert(regIndices.first)).second) 
                            regCounts[regIndices.first] = 1;
                        else 
                            regCounts[regIndices.first]++;

                        nHitM++;
                    }
                    index_i = index_j-1;
                    if(index_j==indices.size()-1)
                        index_i++;
                    
                    break;
                }//if jth data out of coinc window

            }//inner loop over data
        }// outer loop over data
    }//loop over side CRTData products
  
    if(fVerbose) {
          mf::LogInfo("CRT") << returnHits.size() << " CRT hits produced!" << '\n'
              << "  nHitC: " << nHitC  << " , nHitD: " << nHitD  << " , nHitM: " << nHitM  << '\n'
              << "  nMisC: " << nMissC << " , nMisD: " << nMissD << " , nMisM: " << nMissM << '\n';
          auto cts = regCounts.begin();
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
sbn::crt::CRTHit CRTHitRecoAlg::FillCRTHit(vector<uint8_t> tfeb_id, map<uint8_t,vector<pair<int,float>>> tpesmap,
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
sbn::crt::CRTHit CRTHitRecoAlg::MakeTopHit(art::Ptr<CRTData> data){

    uint8_t mac = data->fMac5;
    if(fCrtutils->MacToType(mac)!='c')
        mf::LogError("CRTHitRecoAlg::MakeTopHit") 
            << "CRTUtils returned wrong type!" << '\n';

    map< uint8_t, vector< pair<int,float> > > pesmap;
    int adid  = fCrtutils->MacToAuxDetID(mac,0); //module ID
    auto const& adGeo = fGeometryService->AuxDet(adid); //module
    string region = fCrtutils->GetAuxDetRegion(adid);

    double hitpoint[3], hitpointerr[3], hitlocal[3];
    TVector3 hitpos (0.,0.,0.);
    float petot = 0., pemax=0., pemaxx=0., pemaxz=0.;
    int adsid_max = -1, nabove=0;
    TVector3 postrig;
    bool findx = false, findz = false;
    int maxx=0, maxz=0;

    for(int chan=0; chan<32; chan++) {

        float pe = (data->fAdc[chan]-fQPed)/fQSlope;
        if(pe<=fPEThresh) continue;
        nabove++;
        int adsid = fCrtutils->ChannelToAuxDetSensitiveID(mac,chan);
        petot += pe;
        pesmap[mac].push_back(std::make_pair(chan,pe));

        //TVector3 postmp = fCrtutils->ChanToLocalCoords(mac,chan);
        //strip along z-direction
        if(adsid < 8 && adsid > -1){
            //hitpos.SetX(pe*postmp.X()+hitpos.X());
            //hitpos.SetX(postmp.X()+hitpos.X());
            if(pe>pemaxx){
                pemaxx = pe;
                maxx = adsid;
            }
            findx = true;   
        }
        //strip along x-direction
        else if(adsid > -1 && adsid < 16 ){
            //hitpos.SetZ(pe*postmp.Z()+hitpos.Z());
            //hitpos.SetZ(postmp.Z()+hitpos.Z());
            if(pe > pemaxz) {
                pemaxz = pe;
                maxz = adsid;
            }
            findz = true;
        }
        else {
            mf::LogError("CRTHitRecoAlg::MakeTopHit")
                << "auxDetSensitive ID out of range!" << '\n';
        }
        //identify trigger channel
        if(pe>pemax) {
            TVector3 postmp = fCrtutils->ChanToLocalCoords(mac,chan);
            adsid_max = chan;
            pemax = pe;
            postrig = postmp;
        }
    }

    TVector3 postmp = fCrtutils->ChanToLocalCoords(mac,maxx*2);
    hitpos.SetX(postmp.X());
    postmp = fCrtutils->ChanToLocalCoords(mac,maxz*2); 
    hitpos.SetZ(postmp.Z());

    if(!findx)
        mf::LogWarning("CRTHitRecoAlg::MakeTopHit") << " no interlayer coincidence found! Missing X coord." << '\n';
    if(!findz)
        mf::LogWarning("CRTHitRecoAlg::MakeTopHit") << " no interlayer coincidence found! Missing Z coord." << '\n';

    //no channels above threshold? return empty hit
    if(nabove==0||!findx||!findz)
        return FillCRTHit({},{},0,0,0,0,0,0,0,0,0,0,"");

    //hitpos*=1.0/petot; //hit position weighted by deposited charge
    //hitpos*=1.0/nabove;
    hitlocal[0] = hitpos.X();
    hitlocal[1] = 0.;
    hitlocal[2] = hitpos.Z();
    
    auto const& adsGeo = adGeo.SensitiveVolume(adsid_max); //trigger strip
    double thit = data->fTs0;
    if(adsid_max<8)
        thit -= hitpos.Z()*fPropDelay;
    else
        thit -= hitpos.X()*fPropDelay;

    adGeo.LocalToWorld(hitlocal,hitpoint); //tranform from module to world coords

    hitpointerr[0] = adsGeo.HalfWidth1()*2/sqrt(12);
    hitpointerr[1] = adGeo.HalfHeight();
    hitpointerr[2] = adsGeo.HalfWidth1()*2/sqrt(12);

    CRTHit hit = FillCRTHit({mac},pesmap,petot,thit,thit,0,hitpoint[0],hitpointerr[0],
                            hitpoint[1],hitpointerr[1],hitpoint[2],hitpointerr[2],region);

    return hit;

} // CRTHitRecoAlg::MakeTopHit

//------------------------------------------------------------------------------------------
sbn::crt::CRTHit CRTHitRecoAlg::MakeBottomHit(art::Ptr<CRTData> data){

    uint8_t mac = data->fMac5;
    map< uint8_t, vector< pair<int,float> > > pesmap;
    int adid  = fCrtutils->MacToAuxDetID(mac,0); //module ID
    auto const& adGeo = fGeometryService->AuxDet(adid); //module
    string region = fCrtutils->GetAuxDetRegion(adid);

    double hitpoint[3], hitpointerr[3], hitlocal[3];
    TVector3 hitpos (0.,0.,0.);
    float petot = 0., pemax=0.;
    int adsid_max = -1, nabove=0;
    TVector3 postrig;
    double xmin=0.,xmax=0.;

    for(int chan=0; chan<64; chan++) {

        float pe = (data->fAdc[chan]-fQPed)/fQSlope;
        if(pe<=fPEThresh) continue;
        nabove++;
        int adsid = fCrtutils->ChannelToAuxDetSensitiveID(mac,chan);
        petot += pe;
        pesmap[mac].push_back(std::make_pair(chan,pe));

        TVector3 postmp = fCrtutils->ChanToLocalCoords(mac,chan);
        //all strips along z-direction
        hitpos.SetX(pe*postmp.X()+hitpos.X());
        if(postmp.X()<xmin)
            xmin = postmp.X();
        if(postmp.X()>xmax)
            xmax = postmp.X();

        //identify trigger channel
        if(pe>pemax) {
            adsid_max = adsid;
            pemax = pe;
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
sbn::crt::CRTHit CRTHitRecoAlg::MakeSideHit(vector<art::Ptr<CRTData>> coinData) {

    vector<uint8_t> macs;
    map< uint8_t, vector< pair<int,float> > > pesmap;
    int adid  = fCrtutils->MacToAuxDetID(coinData[0]->fMac5,0); //module ID
    auto const& adGeo = fGeometryService->AuxDet(adid); //module
    string region = fCrtutils->GetAuxDetRegion(adid);

    double hitpoint[3], hitpointerr[3];
    TVector3 hitpos (0.,0.,0.);
    float petot = 0., pemax = 0., pex=0., pey=0.;
    int adsid_max = -1, nabove=0, nx=0, ny=0;
    TVector3 postrig;
    //map<int,vector<double>> ttrigs;
    vector<int> ttrigs;
    double zmin=DBL_MAX, zmax = -DBL_MAX;
    double ymin=DBL_MAX, ymax = -DBL_MAX;
    double xmin=DBL_MAX, xmax = -DBL_MAX;
    std::set<int> layID;

    //loop over FEBs
    for(auto const& data : coinData) {

        macs.push_back(data->fMac5);
        adid  = fCrtutils->MacToAuxDetID(macs.back(),0);

        //loop over channels
        for(int chan=0; chan<32; chan++) {

            float pe = (data->fAdc[chan]-fQPed)/fQSlope;
            if(pe<=fPEThresh) continue;
            nabove++;

            int adsid = fCrtutils->ChannelToAuxDetSensitiveID(macs.back(),chan);
            petot += pe;
            pesmap[macs.back()].push_back(std::make_pair(chan,pe));

            //inner or outer layer
            int layer = fCrtutils->GetMINOSLayerID(adid);
            layID.insert(layer);    
            TVector3 postmp = fCrtutils->ChanToWorldCoords(macs.back(),chan);

            //East/West Walls (all strips along z-direction) or
            // North/South inner walls (all strips along x-direction)
            if(!(region=="South" && layer==1)) {
                hitpos.SetY(pe*postmp.Y()+hitpos.Y());
                ny++;
                pey+=pe;
                if(postmp.Y()<ymin)
                    ymin = postmp.Y();
                if(postmp.Y()>ymax)
                    ymax = postmp.Y();
                if(region!="South") { //region is E/W/N
                    hitpos.SetX(pe*postmp.X()+hitpos.X());
                    nx++;
                    pex+=pe;
                    if(postmp.X()<xmin)
                        xmin = postmp.X();
                    if(postmp.X()>xmax)
                        xmax = postmp.X();
                }
            } 
            else { //else vertical strips in South wall
                hitpos.SetX(pe*postmp.X()+hitpos.X());
                nx++;
                pex+=pe;
                if(postmp.X()<xmin)
                    xmin = postmp.X();
                if(postmp.X()>xmax)
                    xmax = postmp.X();
            }

            //nz = ny
            hitpos.SetZ(pe*postmp.Z()+hitpos.Z());
            if(postmp.X()<xmin)
                zmin = postmp.X();
            if(postmp.X()>xmax)
                zmax = postmp.X();

            //identify trigger channel
            if(pe>pemax) {
                adsid_max = adsid;
                pemax = pe;
                postrig = postmp;
            }


        }//loop over channels

        //correct trigger time for propegation delay
        auto const& adsGeo = adGeo.SensitiveVolume(adsid_max); //trigger stripi
        
       // ttrigs[layer].push_back(data->fTs0);// - adsGeo.HalfLength()*fPropDelay);
        ttrigs.push_back(data->fTs0 - adsGeo.HalfLength()*fPropDelay);

    }//loop over FEBs

    //no channels above threshold or no intermodule coincidences? return empty hit
    if(nabove==0 || layID.size()!=2) {
        if(nabove==0) std::cout << "no channels above threshold!" << std::endl;
        if(layID.size()<2) std::cout << "no coincidence found" << std::endl;
        return FillCRTHit({},{},0,0,0,0,0,0,0,0,0,0,"");
    }

    //finish averaging and fill hit point array
    if(region=="South") {
        hitpos.SetX(hitpos.X()*1.0/pex);
        hitpos.SetY(hitpos.Y()*1.0/pey);
        hitpos.SetZ(hitpos.Z()*1.0/petot);
    }
    else
        hitpos*=1.0/petot; //hit position weighted by deposited charge

    hitpoint[0] = hitpos.X();
    hitpoint[1] = hitpos.Y();
    hitpoint[2] = hitpos.Z();

    //time stamp averaged over all FEBs
    double thit = 0.;//, thit_0 = 0., thit_1 = 0.;
    //for(double const t : ttrigs[0]) 
    //    thit_0 += t;
    for(double const t : ttrigs)
        thit += t;
    //thit_0
    thit*=1.0/ttrigs.size();

    //error estimates (likely need to be revisted)
    auto const& adsGeo = adGeo.SensitiveVolume(adsid_max);
    if(region!="North" && region!="South"){
        hitpointerr[0] = (xmax-xmin)/sqrt(12);
        hitpointerr[1] = (ymax-ymin)/sqrt(12);
        hitpointerr[2] = adsGeo.Length()/sqrt(12);
        //thit=(thit_0
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
