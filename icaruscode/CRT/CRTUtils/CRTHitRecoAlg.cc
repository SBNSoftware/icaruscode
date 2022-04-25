#include "CRTHitRecoAlg.h"
#include <algorithm>
using namespace icarus::crt;

//----------------------------------------------------------------------
CRTHitRecoAlg::CRTHitRecoAlg(const Config& config){
  this->reconfigure(config);
  fChannelMap = art::ServiceHandle<icarusDB::IICARUSChannelMap const>{}.get();  
  fGeometryService  = lar::providerFrom<geo::Geometry>();
  fCrtutils = new CRTCommonUtils();
}

//---------------------------------------------------------------------
CRTHitRecoAlg::CRTHitRecoAlg(){
  fChannelMap = art::ServiceHandle<icarusDB::IICARUSChannelMap const>{}.get();
  fGeometryService = lar::providerFrom<geo::Geometry>();
  fCrtutils = new CRTCommonUtils();
}

//---------------------------------------------------------------------
void CRTHitRecoAlg::reconfigure(const Config& config){
    fVerbose = config.Verbose();
    fUseReadoutWindow = config.UseReadoutWindow(); 
    fQPed = config.QPed();
    fQSlope = config.QSlope();
    fPropDelay = config.PropDelay(); 
    fPEThresh = config.PEThresh();
    fCoinWindow = config.CoinWindow();
    fCrtWindow = config.CrtWindow();
    foutCSVFile = config.outCSVFile();
    fCSVFile = config.CSVFile();
    fData = config.Data();
    if (foutCSVFile) 
      filecsv.open(fCSVFile.c_str());
    return;
}

//---------------------------------------------------------------------------------------
vector<art::Ptr<CRTData>> CRTHitRecoAlg::PreselectCRTData(vector<art::Ptr<CRTData>>& crtList, uint64_t trigger_timestamp){
  if (fVerbose)  mf::LogInfo("CRTHitRecoAlg: ") << "In total " << crtList.size() << " CRTData found in an event" << '\n';
  vector<art::Ptr<CRTData>> crtdata;

  for (size_t febdat_i=0; febdat_i<crtList.size(); febdat_i++) {
    
    uint8_t mac = crtList[febdat_i]->fMac5;
    int adid    = fCrtutils->MacToAuxDetID(mac,0);
    char type   = fCrtutils->GetAuxDetType(adid);

    /// Looking for data within +/- 3ms within trigger time stamp
    /// Here t0 - trigger time -ve, only adding 1s makes the value +ve or -ve
    if (fData && (std::fabs(int64_t(crtList[febdat_i]->fTs0 - trigger_timestamp) + 1'000'000'000) > fCrtWindow)) continue;

    //==== Check if any channel has pe>fPEThresh
    bool passPECut = false;

    if ( type == 'm'){
      for(int chan=0; chan<32; chan++) {
        std::pair<double,double> const chg_cal = fChannelMap->getSideCRTCalibrationMap((int)crtList[febdat_i]->fMac5,chan);
        float pe = (crtList[febdat_i]->fAdc[chan]-chg_cal.second)/chg_cal.first;
        if(pe>fPEThresh){
          passPECut = true;
          continue;
        }
        if(fVerbose)
      mf::LogInfo("CRTHitRecoAlg: ") << "\nfebP (mac5, channel, gain, pedestal, adc, pe) = (" << (int)crtList[febdat_i]->fMac5 << ", " << chan << ", " << chg_cal.first << ", " << chg_cal.second << "," << crtList[febdat_i]->fAdc[chan] << "," << pe << ")\n";
      }
    }else if ( type == 'c' ) {
      for(int chan=0; chan<32; chan++) {
        float pe = (crtList[febdat_i]->fAdc[chan]-fQPed)/fQSlope;
        if(pe>fPEThresh){
          passPECut = true;
          continue;
        }
      }
    }else if ( type == 'd'){
      for(int chan=0; chan<64; chan++) {
        float pe = (crtList[febdat_i]->fAdc[chan]-fQPed)/fQSlope;
        if(pe>fPEThresh){
          passPECut = true;
          continue;
        }
      }
    }
  
    if (passPECut) crtdata.push_back(crtList[febdat_i]);
    
  }
  mf::LogInfo("CRTHitRecoAlg:") << "Found " << crtdata.size() << " after preselection "<< '\n';
  return crtdata;  
}

//---------------------------------------------------------------------------------------
vector<pair<sbn::crt::CRTHit, vector<int>>> CRTHitRecoAlg::CreateCRTHits(vector<art::Ptr<CRTData>> &crtList) {
  
  vector<pair<CRTHit, vector<int>>> returnHits;
  vector<int> dataIds;
  
    uint16_t nMissC = 0, nMissD = 0, nMissM = 0, nHitC = 0, nHitD = 0, nHitM = 0;
    if (fVerbose) mf::LogInfo("CRTHitRecoAlg: ") << "Found " << crtList.size() << " FEB events" << '\n';

    map<string,int> regCounts;
    std::set<string> regs;
    map<string,vector<size_t>> sideRegionToIndices;
    
    // sort by the time 
    std::sort(crtList.begin(), crtList.end(), compareBytime);        
    
    //loop over time-ordered CRTData
    for (size_t febdat_i=0; febdat_i<crtList.size(); febdat_i++) {
  
      uint8_t mac = crtList[febdat_i]->fMac5;
      int adid  = fCrtutils->MacToAuxDetID(mac,0); //module ID
      
      string region = fCrtutils->GetAuxDetRegion(adid);
      char type = fCrtutils->GetAuxDetType(adid);
      CRTHit hit;

      dataIds.clear();

      //CERN modules (intramodule coincidence)
      if ( type == 'c' ) {
        hit = MakeTopHit(crtList[febdat_i]);
        if(IsEmptyHit(hit))
            nMissC++;
        else {
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
          dataIds.push_back(febdat_i);
          returnHits.push_back(std::make_pair(hit,dataIds));
          if ((regs.insert(region)).second) regCounts[region] = 1;
          else regCounts[region]++;
        
          nHitD++;
        }
      }

      if ( type == 'm' ){
        sideRegionToIndices[region].push_back(febdat_i);
      }

    }//loop over CRTData products

    //==== Looping over Side CRT regions
    for(auto const& regIndices : sideRegionToIndices) {

      if(fVerbose) 
        mf::LogInfo("CRTHitRecoAlg: ") << "searching for side CRT hits in region, " << regIndices.first << '\n';
      
      const vector<size_t>& indices = regIndices.second;
      
      if(fVerbose)
        mf::LogInfo("CRTHitRecoAlg: ") << "number of hits associated to this region : " << indices.size() << '\n';

      for(size_t index_i=0; index_i < indices.size(); index_i++) {

        //uint8_t mac = crtList[indices[index_i]]->fMac5;
        //int adid  = fCrtutils->MacToAuxDetID(mac,0); //module ID

        dataIds.clear();
        dataIds.push_back(indices[index_i]);
        vector<art::Ptr<CRTData>> coinData = {crtList[indices[index_i]]};

        if(fVerbose)
          mf::LogInfo("CRTHitRecoAlg: ") << "size ..  " << coinData.size()
            << " data products enetring to time ordring" << '\n';
    
        //inner loop over data after data_i in time
        for (size_t index_j=index_i+1; index_j<indices.size(); index_j++) {
      
          if(fVerbose)
            mf::LogInfo("CRTHitRecoAlg: ") << "i \t"<<index_i << ", j \t"<<index_j << "\t"<<crtList[indices[index_j]]->fTs0 << "\t"<<crtList[indices[index_i]]->fTs0 
            << "\t"<<crtList[indices[index_i]]->fTs0+ fCoinWindow <<'\n';
        
          if(crtList[indices[index_j]]->fTs0 < crtList[indices[index_i]]->fTs0)
            mf::LogError("CRTHitRecoAlg::CreateCRTHits") <<
            "bad time ordering!" << '\n';
        
          if(fVerbose)
            mf::LogInfo("CRTHitRecoAlg: ") << "size ..  " << coinData.size()
                          << " data products before coincidence" << '\n';

          if( (crtList[indices[index_j]]->fTs0 >= crtList[indices[index_i]]->fTs0 && 
              (crtList[indices[index_j]]->fTs0 - crtList[indices[index_i]]->fTs0) < fCoinWindow) ||
              (crtList[indices[index_j]]->fTs0 < crtList[indices[index_i]]->fTs0 && 
              (crtList[indices[index_i]]->fTs0 - crtList[indices[index_j]]->fTs0) < fCoinWindow)) {
            if(fVerbose)
            mf::LogInfo("CRTHitRecoAlg: ") <<  " in coincidence: i \t " << index_i << " ,j: \t" << index_j <<",i mac: \t" 
              << (int)crtList[indices[index_i]]->fMac5 << ", j mac: \t" <<(int)crtList[indices[index_j]]->fMac5<< '\n';
          
            coinData.push_back(crtList[indices[index_j]]);
            dataIds.push_back(indices[index_j]);
          }
            
          //out of coinWindow
          if( (crtList[indices[index_j]]->fTs0 - crtList[indices[index_i]]->fTs0) > fCoinWindow || index_j==indices.size()-1){
            if(fVerbose)
              mf::LogInfo("CRTHitRecoAlg: ") <<  "out of coincidence  " << index_j << "\t" << indices.size() <<"\t" <<indices.size()-1
                << " data products..." << '\n';
            if(fVerbose)
              mf::LogInfo("CRTHitRecoAlg: ") << "attempting to produce MINOS hit from " << coinData.size() 
                << " data products..." << '\n';

            CRTHit hit = MakeSideHit(crtList, dataIds);

            if(IsEmptyHit(hit)){
              nMissM++;
            }
            else {
              if(fVerbose)
                mf::LogInfo("CRTHitRecoAlg: ") << "MINOS hit produced" << '\n';
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
                            float peshit, uint64_t time0, uint64_t time1, int plane, 
                            double x, double ex, double y, double ey, double z, double ez, string tagger){
  CRTHit crtHit;
  crtHit.feb_id      = tfeb_id;
  crtHit.pesmap      = tpesmap;
  crtHit.peshit      = peshit;
  crtHit.ts0_s_corr  = time0 / 1'000'000'000; 
  crtHit.ts0_ns      = time0;
  crtHit.ts0_ns_corr = time0; 
  crtHit.ts1_ns      = time1;
  crtHit.ts0_s       = time0 / 1'000'000'000;
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
    int plane = fCrtutils->AuxDetRegionNameToNum(region);
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
    uint64_t thit = data->fTs0;

    if(adsid_max<8)
        thit -= hitpos.Z()*fPropDelay;
    else
        thit -= hitpos.X()*fPropDelay;

    adGeo.LocalToWorld(hitlocal,hitpoint); //tranform from module to world coords

    hitpointerr[0] = adsGeo.HalfWidth1()*2/sqrt(12);
    hitpointerr[1] = adGeo.HalfHeight();
    hitpointerr[2] = adsGeo.HalfWidth1()*2/sqrt(12);

    CRTHit hit = FillCRTHit({mac},pesmap,petot,thit,thit,plane,hitpoint[0],hitpointerr[0],
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
    int plane =fCrtutils->AuxDetRegionNameToNum(region);
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
    uint64_t thit = data->fTs0 - adsGeo.HalfLength()*fPropDelay;
    
    adGeo.LocalToWorld(hitlocal,hitpoint); //tranform from module to world coords

    hitpointerr[0] = (xmax-xmin+2*adsGeo.HalfWidth1()*2)/sqrt(12);
    hitpointerr[1] = adGeo.HalfHeight();
    hitpointerr[2] = adsGeo.Length()/sqrt(12);

    CRTHit hit = FillCRTHit({mac},pesmap,petot,thit,thit,plane,hitpoint[0],hitpointerr[0],
                            hitpoint[1],hitpointerr[1],hitpoint[2],hitpointerr[2],region);

    return hit;

} // CRTHitRecoAlg::MakeBottomHit

//-----------------------------------------------------------------------------------

//==== idxList : vector of indeces that 1) share same (Side)CRT region 2) time coincidence (150ns)
sbn::crt::CRTHit CRTHitRecoAlg::MakeSideHit(vector<art::Ptr<CRTData>>& crtList, vector<int>& idxList){

  std::map<int, vector<int>> map_ADID_to_Indicies;

  for(auto const& idx : idxList) {
    int this_mac5 = crtList[idx]->fMac5;
    int this_adid  = fCrtutils->MacToAuxDetID(this_mac5,0); // but using channel==0
    map_ADID_to_Indicies[this_adid].push_back(idx);
  }

  bool HitFromLayer0(false);
  float pemaxFromLayer0(-DBL_MAX);
  CRTHit hitFromLayer0;
  bool HitFromLayer1(false);
  float pemaxFromLayer1(-DBL_MAX);
  CRTHit hitFromLayer1;
  for(std::map<int, vector<int>>::const_iterator it=map_ADID_to_Indicies.begin(); it!=map_ADID_to_Indicies.end(); ++it){
    vector<int> idxList = it->second;
    CRTHit this_hit = MakeSideHitPerModule(crtList, idxList);
    int this_layer = fCrtutils->GetMINOSLayerID(it->first);

    if(this_layer==0){
      HitFromLayer0 = true;
      if( this_hit.peshit > pemaxFromLayer0 ){
        pemaxFromLayer0 = this_hit.peshit;
        hitFromLayer0 = this_hit;
      }
    }

    if(this_layer==1){
      HitFromLayer1 = true;
      if( this_hit.peshit > pemaxFromLayer1 ){
        pemaxFromLayer1 = this_hit.peshit;
        hitFromLayer1 = this_hit;
      }
    }

  }

  if(!HitFromLayer0 || !HitFromLayer1){
    return FillCRTHit({},{},0,0,0,0,0,0,0,0,0,0,"");
  }

  //==== merge per-layer information here
  sbn::crt::CRTHit out_hit = MergeSideHits( {hitFromLayer0,hitFromLayer1} );

  //==== per-layer info
  out_hit.nLayer = 2;
  out_hit.layerHits = new sbn::crt::CRTLayerHit[2];
  out_hit.layerHits[0] = {
    hitFromLayer0.peshit,
    hitFromLayer0.ts0_ns,
    hitFromLayer0.ts1_ns,
    hitFromLayer0.x_pos,
    hitFromLayer0.x_err,
    hitFromLayer0.y_pos,
    hitFromLayer0.y_err,
    hitFromLayer0.z_pos,
    hitFromLayer0.z_err
  };
  out_hit.layerHits[1] = {
    hitFromLayer1.peshit,
    hitFromLayer1.ts0_ns,
    hitFromLayer1.ts1_ns,
    hitFromLayer1.x_pos,
    hitFromLayer1.x_err,
    hitFromLayer1.y_pos,
    hitFromLayer1.y_err,
    hitFromLayer1.z_pos,
    hitFromLayer1.z_err
  };

  //==== direction : outside-in vector
  //==== East wall : layer0 = outer
  //==== Other : layer1 = outer
  TVector3 vec_dir(
    hitFromLayer0.x_pos - hitFromLayer1.x_pos,
    hitFromLayer0.y_pos - hitFromLayer1.y_pos,
    hitFromLayer0.z_pos - hitFromLayer1.z_pos
  ); // 1->0 vector
  vec_dir = vec_dir.Unit();
  //==== If East, flip the direction to make it outside-in
  if(hitFromLayer0.plane>=43 && hitFromLayer0.plane<=45) vec_dir *= -1.;
  out_hit.x_dir = vec_dir.X();
  out_hit.y_dir = vec_dir.Y();
  out_hit.z_dir = vec_dir.Z();

  return out_hit;

}

sbn::crt::CRTHit CRTHitRecoAlg::MakeSideHitPerModule(vector<art::Ptr<CRTData>>& crtList, vector<int>& idxList){

  const int nFEBs = idxList.size();
  //==== mac5
  vector<uint8_t> mac5s;
  //==== pe
  float pe_sum(0.); // sum of pe above fPEThresh
  map<uint8_t, vector<pair<int,float>>> pe_map; // key = mac5, value = vector of pair, (channel, pe)
  vector<int> vec_maxPEChanNumber; // channel number that gives maximum pe for each FEB
  //==== timing variables
  uint64_t t0_avg(0), t1_avg(0);
  //==== plane
  string regionName = fCrtutils->GetAuxDetRegion(fCrtutils->MacToAuxDetID(crtList[idxList[0]]->fMac5,0)); // = region
  int regionNum = fCrtutils->AuxDetRegionNameToNum(regionName); // = plane
  //==== position
  int nChAboveThsh(0); // number of channels used for the position average
  TVector3 hitpos(0.,0.,0.);
  //==== position err
  double x_err(-DBL_MAX), y_err(-DBL_MAX), z_err(-DBL_MAX);
  //==== Detector lengths
  //==== https://github.com/LArSoft/larcorealg/blob/9758633b30fcaa5586c4a28befb02c8d2d3494e6/larcorealg/Geometry/AuxDetGeo.h#L102-L105
  double geoHalfWidth1(-DBL_MAX);
  double geoHalfHeight(-DBL_MAX);
  double geoLength(-DBL_MAX);

  //==== loop over FEBs
  for(auto const& idx : idxList){

    auto const& data = crtList[idx];

    //==== mac5
    uint8_t this_mac5 = data->fMac5;
    mac5s.push_back(this_mac5);
    //==== adid (module id)
    int this_adid  = fCrtutils->MacToAuxDetID(this_mac5,0); //module ID
    auto const& adGeo = fGeometryService->AuxDet(this_adid); //module

    //==== loop over channels
    float tmp_maxPE(-DBL_MAX);
    int tmp_maxPEChan(-1);
    for(int chan=0; chan<32; chan++){
      //==== Calib. db.
      std::pair<double,double> const chg_cal = fChannelMap->getSideCRTCalibrationMap((int)this_mac5,chan);
      float pe = (data->fAdc[chan]-chg_cal.second)/chg_cal.first;
      //==== pe cut
      if(pe<=fPEThresh) continue;

      nChAboveThsh++;

      //==== pe sum
      pe_sum += pe;
      pe_map[this_mac5].push_back(std::make_pair(chan,pe));
      //==== position
      TVector3 this_pos = fCrtutils->ChanToWorldCoords(this_mac5,chan);
      hitpos += pe*this_pos; // pe-weighted

      if(tmp_maxPE < pe){
        tmp_maxPE = pe;
        tmp_maxPEChan = chan;

        int this_adsid = fCrtutils->ChannelToAuxDetSensitiveID(this_mac5,chan);
        auto const& adsGeo = adGeo.SensitiveVolume(this_adsid);
        geoHalfWidth1 = adsGeo.HalfWidth1();
        geoHalfHeight = adsGeo.HalfHeight();
        geoLength = adsGeo.Length();
      }

    }

    vec_maxPEChanNumber.push_back(tmp_maxPEChan);

    //==== timing variables
    t0_avg += data->fTs0;
    t1_avg += data->fTs1;

  }

  //==== When no channel with pe>thrh, return an empty hit
  if(nChAboveThsh==0 || pe_sum==0.){
    return FillCRTHit({},{},0,0,0,0,0,0,0,0,0,0,"");
  }

  hitpos *= 1./pe_sum; // pe-weighted
  t0_avg /= nFEBs;
  t1_avg /= nFEBs;

  //==== z-position calculation for East/West module
  if(40<=regionNum && regionNum<=45){

    if(pe_map.size()==2){ // i.e., checking number of unique mac5s; i.e., the two readouts

      int prev_mac5_int(-999);
      double relZPosFromCeneter(-999.);
      for(auto const& idx : idxList){
        auto ii = &idx - idxList.data();
        int this_mac5_int = (int)crtList[idx]->fMac5;
        if(abs(prev_mac5_int-this_mac5_int)==1){
          uint64_t t0_ro0 = crtList[idx-1]->fTs0;
          uint64_t t0_ro1 = crtList[idx]->fTs0;
          float zaxixpos = 0.5*(int64_t(t0_ro1 - t0_ro0)/fPropDelay);
          relZPosFromCeneter = zaxixpos;

          int this_adid  = fCrtutils->MacToAuxDetID(crtList[idx]->fMac5,0);
          auto const& adGeo = fGeometryService->AuxDet(this_adid);
          int adsid = fCrtutils->ChannelToAuxDetSensitiveID(crtList[idx-1]->fMac5,vec_maxPEChanNumber[ii-1]);
          auto const& adsGeo = adGeo.SensitiveVolume(adsid);

          geo::Point_t tmp_pos = adsGeo.GetCenter() + geo::Zaxis() * zaxixpos;

          hitpos.SetZ(tmp_pos.Z());

          break;
        }
        else{
          prev_mac5_int = this_mac5_int;
        }
      }

      x_err = 2.*geoHalfWidth1; // Depth of the strip
      y_err = 2.*2.*geoHalfHeight/std::sqrt(12); // TODO two strips are coupled right..?
      z_err = relZPosFromCeneter/std::sqrt(12); //TODO this one is just random.. should be fixed

    }
    else{

      //==== TODO set to center

      for(auto const& idx : idxList){
        auto ii = &idx - idxList.data();
        int this_adid  = fCrtutils->MacToAuxDetID(crtList[idx]->fMac5,0);
        auto const& adGeo = fGeometryService->AuxDet(this_adid); 
        int adsid = fCrtutils->ChannelToAuxDetSensitiveID(crtList[idx]->fMac5,vec_maxPEChanNumber[ii]);
        auto const& adsGeo = adGeo.SensitiveVolume(adsid); 
        geo::Point_t tmp_pos = adsGeo.GetCenter();

        hitpos.SetZ(tmp_pos.Z());

        break;
      }

      x_err = 2.*geoHalfWidth1; // Depth of the strip
      y_err = 2.*2.*geoHalfHeight/std::sqrt(12); // TODO two strips are coupled right..?
      z_err = -geoLength/std::sqrt(12); // intended to be negative.. look error treament at MergeSideHits()

    }
  }
  //==== South wall
  else if(regionNum==46){
    x_err = 2.*2.*geoHalfHeight/std::sqrt(12); // TODO two strips are coupled right..?
    y_err = 2.*2.*geoHalfHeight/std::sqrt(12); // TODO two strips are coupled right..?
    z_err = 2.*geoHalfWidth1; // Depth of the strip
  }
  //==== North
  else if(regionNum==47){
    x_err = geoLength;
    y_err = 2.*2.*geoHalfHeight/std::sqrt(12); // TODO two strips are coupled right..?
    z_err = 2.*geoHalfWidth1; // Depth of the strip
  }
  else{
    mf::LogInfo("CRTHitRecoAlg:") << "This is side CRT, but regionNum = " << regionNum << "\n";
    abort();
  }

  return FillCRTHit(
    mac5s, // vector<uint8_t> tfeb_id
    pe_map, // map<uint8_t, vector<pair<int,float>>> tpesmap
    pe_sum, // float peshi
    t0_avg, // uint64_t time0
    t1_avg, // uint64_t time1
    regionNum, // int plane
    hitpos.X(), // double x
    x_err, // double ex
    hitpos.Y(), // double y
    y_err, // double ey
    hitpos.Z(), // double z
    z_err, // double ez
    regionName // string tagger
  );

}

sbn::crt::CRTHit CRTHitRecoAlg::MergeSideHits(const vector<sbn::crt::CRTHit>& crtHits){

  //==== For now, crtHits has two elements
  //==== 0th : Layer 0
  //==== 1th : Layer 1

  vector<uint8_t> mac5s;
  map<uint8_t, vector<pair<int,float>>> pe_map;
  float pe_sum(0.);
  uint64_t t0(0), t1(0);
  int regionNum;
  double pos_x(0.), pos_y(0.), pos_z(0.);
  double pos_x_err(0.), pos_y_err(0.), pos_z_err(0.);
  std::string regionName;

  for(auto const& crtHit : crtHits){
    mac5s.insert(mac5s.end(), crtHit.feb_id.begin(), crtHit.feb_id.end());
    pe_map.insert(crtHit.pesmap.begin(), crtHit.pesmap.end());
    pe_sum += crtHit.peshit;
    t0 += crtHit.ts0_ns;
    t1 += crtHit.ts1_ns;
    regionNum = crtHit.plane;
    regionName = crtHit.tagger;
  }
  t0 /= crtHits.size();
  t1 /= crtHits.size();

  //==== Position determination
  //==== - East or West: both layers are horizontal
  if(regionNum>=40 &&regionNum<=45){

    pos_x = (crtHits[0].x_pos + crtHits[1].x_pos)/2.;
    pos_y = (crtHits[0].y_pos + crtHits[1].y_pos)/2.;
    //==== TODO how should we add the errs from two layers?
    pos_x_err = std::sqrt( crtHits[0].x_err*crtHits[0].x_err + crtHits[1].x_err*crtHits[1].x_err );
    pos_y_err = std::sqrt( crtHits[0].y_err*crtHits[0].y_err + crtHits[1].y_err*crtHits[1].y_err );

    //==== From MakeSideHitPerModule(), z_err is set to negative when the two-readout method failed
    if(crtHits[0].z_err>0 && crtHits[1].z_err>0){
      pos_z = (crtHits[0].z_pos + crtHits[1].z_pos)/2.;
      //==== TODO how should we add the errs from two layers?
      pos_z_err = std::sqrt( crtHits[0].z_err*crtHits[0].z_err + crtHits[1].z_err*crtHits[1].z_err );
    }
    else if(crtHits[0].z_err<=0 || crtHits[1].z_err<=0){
      pos_z = crtHits[0].z_err<=0 ? crtHits[1].z_pos : crtHits[0].z_pos;
      pos_z_err = crtHits[0].z_err<=0 ? crtHits[1].z_err : crtHits[0].z_err;
    }
    else{
      pos_z = (crtHits[0].z_pos + crtHits[1].z_pos)/2.;
      //==== TODO how should we add the errs from two layers?
      pos_z_err = std::sqrt( crtHits[0].z_err*crtHits[0].z_err + crtHits[1].z_err*crtHits[1].z_err );
    }
  }
  //==== South
  //==== layer 0 : Inner, horizontal -> use y position
  //==== layer 1 : Outer, vertical -> use x position
  else if(regionNum==46){
    pos_x = crtHits[1].x_pos;
    pos_y = crtHits[0].y_pos;
    pos_z = (crtHits[0].z_pos + crtHits[1].z_pos)/2.;

    pos_x_err = crtHits[1].x_err;
    pos_y_err = crtHits[0].y_err;
    pos_z_err = std::sqrt( crtHits[0].z_err*crtHits[0].z_err + crtHits[1].z_err*crtHits[1].z_err );
  }
  //==== North
  //==== both horizontal
  else if(regionNum==47){
    pos_x = (crtHits[0].x_pos + crtHits[1].x_pos)/2.;
    pos_y = (crtHits[0].y_pos + crtHits[1].y_pos)/2.;
    pos_z = (crtHits[0].z_pos + crtHits[1].z_pos)/2.;

    pos_x_err = std::sqrt( crtHits[0].x_err*crtHits[0].x_err + crtHits[1].x_err*crtHits[1].x_err );
    pos_y_err = std::sqrt( crtHits[0].y_err*crtHits[0].y_err + crtHits[1].y_err*crtHits[1].y_err );
    pos_z_err = std::sqrt( crtHits[0].z_err*crtHits[0].z_err + crtHits[1].z_err*crtHits[1].z_err );
  }
  else{
    mf::LogInfo("CRTHitRecoAlg:") << "This is side CRT, but regionNum = " << regionNum << "\n";
    abort();
  }


  return FillCRTHit(
    mac5s, // vector<uint8_t> tfeb_id
    pe_map, // map<uint8_t, vector<pair<int,float>>> tpesmap
    pe_sum, // float peshi
    t0, // uint64_t time0
    t1, // uint64_t time1
    regionNum, // int plane
    pos_x, // double x
    pos_x_err, // double ex
    pos_y, // double y
    pos_y_err, // double ey
    pos_z, // double z
    pos_z_err, // double ez
    regionName // string tagger
  );

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
