#include "CRTHitRecoAlg.h"

using namespace icarus::crt;

//string filename = "filename.txt";

//----------------------------------------------------------------------
CRTHitRecoAlg::CRTHitRecoAlg(const Config& config){
    this->reconfigure(config);
  
    fGeometryService  = lar::providerFrom<geo::Geometry>();
    fCrtutils = new CRTCommonUtils();
}

//---------------------------------------------------------------------
CRTHitRecoAlg::CRTHitRecoAlg(){
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
    foutCSVFile = config.outCSVFile();
    fCSVFile = config.CSVFile();
    if (foutCSVFile) 
      filecsv.open(fCSVFile.c_str());
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
        //      std::cout << "In CRTHitRecoAlg::CreateCRTHits functions mac is " << (int)mac 
        //        << "  with module number " << adid <<std::endl; 
        string region = fCrtutils->GetAuxDetRegion(adid);
        char type = fCrtutils->GetAuxDetType(adid);
        CRTHit hit;

        //if(fVerbose) 
        //std::cout << "found data with mac5 = " << (int)mac << ", " << string(1,type) 
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

      //      if(!(regIndices.first=="East North")) continue;
        if(fVerbose) 
            std::cout << "searching for side CRT hits in region, " << regIndices.first << std::endl;
    
        vector<size_t> indices = regIndices.second;
        
        for(size_t index_i=0; index_i < indices.size(); index_i++) {
          
          dataIds.clear();
          dataIds.push_back(indices[index_i]);
          vector<art::Ptr<CRTData>> coinData = {crtList[indices[index_i]]};
          
          if(fVerbose)
            std::cout << "size ..  " << coinData.size()
                        << " data products enetring to time ordring" << std::endl;
          
          //inner loop over data after data_i in time
          for (size_t index_j=index_i+1; index_j<indices.size(); index_j++) {
            
            if(fVerbose)
              std::cout << "\t"<<index_i << "\t"<<index_j << "\t"<<crtList[indices[index_j]]->fTs0 << "\t"<<crtList[indices[index_i]]->fTs0 
                        << "\t"<<crtList[indices[index_i]]->fTs0+ fCoinWindow <<std::endl;
            
            if(crtList[indices[index_j]]->fTs0 < crtList[indices[index_i]]->fTs0)
              mf::LogError("CRTHitRecoAlg::CreateCRTHits") <<
                "bad time ordering!" << '\n';
            //          if(fVerbose) std::cout<< "should not enter this line" << std::endl;
            //}
            
            if(fVerbose)
              std::cout << "size ..  " << coinData.size()
                        << " data products before coincidence" << std::endl;
            // in coincidence
            //      if(crtList[indices[index_j]]->fTs0 <= crtList[indices[index_i]]->fTs0 + fCoinWindow) {
	    //            if(std::llabs(crtList[indices[index_j]]->fTs0 - crtList[indices[index_i]]->fTs0) < fCoinWindow) {
	    if( (crtList[indices[index_j]]->fTs0>=crtList[indices[index_i]]->fTs0 && 
		 (crtList[indices[index_j]]->fTs0 - crtList[indices[index_i]]->fTs0) < fCoinWindow) ||
		(crtList[indices[index_j]]->fTs0<crtList[indices[index_i]]->fTs0 && 
		 (crtList[indices[index_i]]->fTs0 - crtList[indices[index_j]]->fTs0) < fCoinWindow)) {
              if(fVerbose)
                std::cout <<  " in coincidence: i \t " << index_i << " ,j: \t" << index_j <<",i mac: \t" 
                          << (int)crtList[indices[index_i]]->fMac5 << ", j mac: \t" <<(int)crtList[indices[index_j]]->fMac5<< std::endl;
              //                  if(fVerbose) std::cout<< "should not enter this line....." << std::endl;
              coinData.push_back(crtList[indices[index_j]]);
              dataIds.push_back(indices[index_j]);
            }
            
            //out of coinWindow
            if(crtList[indices[index_j]]->fTs0 > crtList[indices[index_i]]->fTs0 + fCoinWindow
                 || index_j==indices.size()-1) 
                { 
                  
                  if(fVerbose)
                    std::cout <<  "out of coincidence  " << index_j << "\t" << indices.size() <<"\t" <<indices.size()-1
                              << " data products..." << std::endl;
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
		    //      std::cout << "++++++++++++++++++++++++++++ line 169" <<std::endl;
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
	    //  std::cout << "++++++++++++++++++++++++++++ line 183 after break" <<std::endl;
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
    //    std::cout << "++++++++++++++++++++++++++++ line 201 after return hit" <<std::endl;
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
    double thit = data->fTs0;
    std::cout << "double thit: " << thit << "\n";
    uint64_t thit_64 = data->fTs0;
    std::cout << "uint64_t thit: " << thit_64 << "\n";

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
    double thit_dub = data->fTs0 - adsGeo.HalfLength()*fPropDelay;
    std::cout << "thit_dub: " << thit_dub << "\n";
    uint64_t thit = data->fTs0 - adsGeo.HalfLength()*fPropDelay;
    std::cout << "thit: " << thit << "\n";
    
    adGeo.LocalToWorld(hitlocal,hitpoint); //tranform from module to world coords

    hitpointerr[0] = (xmax-xmin+2*adsGeo.HalfWidth1()*2)/sqrt(12);
    hitpointerr[1] = adGeo.HalfHeight();
    hitpointerr[2] = adsGeo.Length()/sqrt(12);

    CRTHit hit = FillCRTHit({mac},pesmap,petot,thit,thit,plane,hitpoint[0],hitpointerr[0],
                            hitpoint[1],hitpointerr[1],hitpoint[2],hitpointerr[2],region);

    return hit;

} // CRTHitRecoAlg::MakeBottomHit

//-----------------------------------------------------------------------------------
sbn::crt::CRTHit CRTHitRecoAlg::MakeSideHit(vector<art::Ptr<CRTData>> coinData) {

  //  std::ofstream filecsv;
  //filecsv.open("filename.txt");

    vector<uint8_t> macs;
    map< uint8_t, vector< pair<int,float> > > pesmap;
    map< uint8_t, vector< pair<int,TVector3> > > chantopos;

    struct infoA {
      uint8_t mac5s;
      int channel;
      uint64_t t0;
      TVector3 pos;
      int strip;
    };

    struct infoB {
      uint8_t mac5s;
      int channel;
      uint64_t t0;
      TVector3 pos;
      int strip;
    };

    vector<infoA> informationA;
    vector<infoA> informationB;

    int adid  = fCrtutils->MacToAuxDetID(coinData[0]->fMac5,0); //module ID
    //    std::cout << "In CRTHitRecoAlg::MakeSideHit functions mac is " << (int)coinData[0]->fMac5 
    //        << "  with module number " << adid <<std::endl; 
    auto const& adGeo = fGeometryService->AuxDet(adid); //module
    string region = fCrtutils->GetAuxDetRegion(adid);
    int plane =fCrtutils->AuxDetRegionNameToNum(region);

    //int nfeb = -1;
    double hitpoint[3], hitpointerr[3];
    TVector3 hitpos (0.,0.,0.);

    float petot = 0., pemax = 0., pex=0., pey=0.;
    int adsid_max = -1, nabove=0, nx=0, ny=0, nz = 0, ntrig = 0;
    TVector3 postrig;
    //map<int,vector<double>> ttrigs;
    vector<uint64_t> ttrigs;
    vector<TVector3> tpos;
    double zmin=DBL_MAX, zmax = -DBL_MAX;
    double ymin=DBL_MAX, ymax = -DBL_MAX;
    double xmin=DBL_MAX, xmax = -DBL_MAX;
    //std::set<int> layID;
    std::vector<int> layID;
    //vector<pair<uint8_t,uint8_t> > mac5;
    //std::vector<int> nfebs;
    //    map<int, int> febtoadid;
    std::vector<int> febA;
    std::vector<int> febB;
    //    std::vector<int> modid;
    
    //loop over FEBs
    for(auto const& data : coinData) {
      //      macs.push_back(data->fMac5);
      //adid  = fCrtutils->MacToAuxDetID(macs.back(),0);
      //      std ::cout << (int)data->fMac5 << "\t"<< (int)fCrtutils->MacToAuxDetID(52,0) <<"\t"<< adid << std::endl; 
      // std ::cout << (int)data->fMac5 << "\t"<< fCrtutils->MacToAuxDetID((int)data->fMac5,0) <<"\t"<< adid << std::endl; 
      if (adid == (int)fCrtutils->MacToAuxDetID((int)data->fMac5,0)){
	//std ::cout << "entering if condition\t" << std::endl;
	febA.push_back(data->fMac5);
      }else {
	//std ::cout << "not going to else........\t" << std::endl;
	febB.push_back(data->fMac5);
      }
      //std ::cout << "end of for loop\t" << std::endl;
      // modid.push_back(fCrtutils->MacToAuxDetID(feb.back(),0));

    }

    if(fVerbose)
      std ::cout << "line 451: size of febA: \t" << (int)febA.size() 
		 << " size of febB: " << (int)febB.size() << std::endl;
    
    
    //   if ((int)febA.size()==0 or (int)febB.size()==0) return FillCRTHit({},{},0,0,0,0,0,0,0,0,0,0,"");
    /*
    if ((int)febA.size()!=0 or (int)febB.size()!= 0){// return FillCRTHit({},{},0,0,0,0,0,0,0,0,0,0,"");
      
      for (int i = 0; i < (int)febA.size(); i++){
	
	if(fVerbose)
	  std::cout << "febA:  " << (int)febA[i]
		    << " ,febB : " << (int)febB[i]<<std::endl;
      }
    }
    */

    //    std ::cout << "line 458 \t" << std::endl;
    //if (feb.size() != 4 || feb.size() != 2) return 0;

    //if (modid[0] == modid[1]) 


    //loop over FEBs
    for(auto const& data : coinData) {

      // if(!(region=="North")) continue;
        macs.push_back(data->fMac5);
        adid  = fCrtutils->MacToAuxDetID(macs.back(),0);

	//        nfebs.push_back(fCrtutils->NFeb(adid));
	// mac5.push_back(fCrtutils->ADToMac(adid));

        //      pair<uint8_t,uint8_t> feb =  
        int layer = fCrtutils->GetMINOSLayerID(adid);
        layID.push_back(layer);   

	auto idx = &data - coinData.data();
	if(fVerbose)	std :: cout << "index: " << idx << std::endl;

	if ((int)febA.size() == 0 or (int)febB.size() == 0) continue;

	for (int i = 0; i < (int)febA.size(); i++){
	  if  (macs.back() == (int)febA[i]) {
	    //loop over channels
	    for(int chan=0; chan<32; chan++) {

	      float pe = (data->fAdc[chan]-fQPed)/fQSlope;
	      if(pe<=fPEThresh) continue;
	      nabove++;

	      int adsidA = fCrtutils->ChannelToAuxDetSensitiveID(macs.back(),chan);
	      TVector3 postmp = fCrtutils->ChanToWorldCoords(macs.back(),chan);

	      informationA.push_back ({macs.back(), chan, data->fTs0, postmp, adsidA});
	    }

	    if(fVerbose)
	      std::cout << "looking for mac " << (int)macs.back()
			<< " and matching with febA : " << (int)febA[i]<<std::endl; 
	  }else if ( macs.back() == (int)febB[i]){

	    //loop over channels
            for(int chan=0; chan<32; chan++) {

              float pe = (data->fAdc[chan]-fQPed)/fQSlope;
              if(pe<=fPEThresh) continue;
              nabove++;

	      int adsidB = fCrtutils->ChannelToAuxDetSensitiveID(macs.back(),chan);
              TVector3 postmp = fCrtutils->ChanToWorldCoords(macs.back(),chan);

              informationB.push_back ({macs.back(), chan, data->fTs0, postmp, adsidB});
            }
	    if(fVerbose)
	      std::cout << "else if looking for mac "<< (int)macs.back()
			<< " and matching with febB : " << (int)febB[i]<<std::endl;
	  }
	}

	if(fVerbose)
	  std::cout << "In CRTHitRecoAlg::MakeSideHit 1st feb is " << (int)fCrtutils->ADToMac(adid).first 
		    << " ,2nd feb :"<< (int)fCrtutils->ADToMac(adid).second  
		    << ", time "<< data->fTs0 
                    << "  with module number " << adid << ", no. of FEB " <<std::endl;

	if(fVerbose)
	  std::cout << "In CRTHitRecoAlg::MakeSideHit functions mac is " << (int)macs.back()
		    << "  with module number " << adid << ", no. of FEB " <<std::endl;

        //loop over channels
        for(int chan=0; chan<32; chan++) {
	  
	  float pe = (data->fAdc[chan]-fQPed)/fQSlope;
	  if(pe<=fPEThresh) continue;
	  nabove++;
	  
	  int adsid = fCrtutils->ChannelToAuxDetSensitiveID(macs.back(),chan);
	  petot += pe;
	  pesmap[macs.back()].push_back(std::make_pair(chan,pe));
	  
	  //inner or outer layer
	  //	  int layer = fCrtutils->GetMINOSLayerID(adid);
	  //layID.insert(layer);    

	  TVector3 postmp = fCrtutils->ChanToWorldCoords(macs.back(),chan);
	  
	  // for each feb store position to channel
	  chantopos[macs.back()].push_back(std::make_pair(chan,postmp));
	  //	  if(fVerbose)
	  //std::cout << " chan:\t "<< chan 
	  //	      << " ,pe: \t"<< pe << ", adc:\t" << data->fAdc[chan] << ", time: \t"<< data->fTs0 << std::endl;
	    
	  if(fVerbose)
	    std::cout <<  "feb: " << (int)macs.back() << " ,chan : \t" << chan
		      << " ,pe: \t"<< pe << ", adc:\t" << data->fAdc[chan] << ", time: \t"<< data->fTs0
                      << " ,x: \t"<< postmp.X() <<" ,y: \t" << postmp.Y()  <<" ,z: \t" << postmp.Z()<< std::endl;

	  //East/West Walls (all strips along z-direction) or
	  // North/South inner walls (all strips along x-direction)
	  // All the horizontal layers measure Y first,
	  if(!(region=="South" && layer==1)) {
	    hitpos.SetY(pe*postmp.Y()+hitpos.Y());
	    // hitpos.SetY(1.0*postmp.Y()+hitpos.Y());
	    if (fVerbose){
	      std::cout << "!(region==South && layer==1) : \t" << " feb: " << (int)macs.back() << " ,chan : \t" << chan
			<< " ,pe: \t"<< pe << ", adc:\t" << data->fAdc[chan] << ", time: \t"<< data->fTs0
			<< " ,x: \t"<< postmp.X() <<" ,y: \t" << postmp.Y()  <<" ,z: \t" << postmp.Z()<< std::endl;
	    }
	    ny++;
	    pey+=pe;
	    if(postmp.Y()<ymin)
	      ymin = postmp.Y();
	    if(postmp.Y()>ymax)
	      ymax = postmp.Y();
	    if(region!="South") { //region is E/W/N
	      hitpos.SetX(pe*postmp.X()+hitpos.X());
	      //hitpos.SetX(1.0*postmp.X()+hitpos.X());
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
	    //  hitpos.SetX(1.0*postmp.X()+hitpos.X());
	    if (fVerbose){
	    std::cout << "vertical strips in South wall : \t" << " feb: " << (int)macs.back() << " ,chan : \t" << chan
                      << " ,pe: \t"<< pe << ", adc:\t" << data->fAdc[chan] << ", time: \t"<< data->fTs0
                      << " ,x: \t"<< postmp.X() <<" ,y: \t" << postmp.Y()  <<" ,z: \t" << postmp.Z()<< std::endl;
	    }
	    nx++;
	    pex+=pe;
	    if(postmp.X()<xmin)
	      xmin = postmp.X();
	    if(postmp.X()>xmax)
	      xmax = postmp.X();
	  }
	  
	  //nz = ny
	  hitpos.SetZ(pe*postmp.Z()+hitpos.Z());
	  //  hitpos.SetZ(1.0*postmp.Z()+hitpos.Z());
	  nz++;
	  if (fVerbose){
	    if(region =="South")
	      std::cout << " South wall z: \t" << " feb: " << (int)macs.back() << " ,chan : \t" << chan
			<< " ,pe: \t"<< pe << ", adc:\t" << data->fAdc[chan] << ", time: \t"<< data->fTs0
			<< " ,x: \t"<< postmp.X() <<" ,y: \t" << postmp.Y()  <<" ,z: \t" << postmp.Z()<< " petotal : "<< petot<< std::endl;
	  }
	  if(postmp.Z()<zmin)
	    zmin = postmp.Z();
	  if(postmp.Z()>zmax)
	    zmax = postmp.Z();
	  
	  //identify trigger channel
	  if(pe>pemax) {
	    adsid_max = adsid;
	    pemax = pe;
	    postrig = postmp;
	  }

	  //	  if(fVerbose)
	  //std::cout <<  "feb: " << (int)macs.back() << " ,chan : \t" << chan 
	  //	      << " ,x: \t"<< hitpos[0] <<" ,y: \t" << hitpos[1]  <<" ,z: \t" << hitpos[2]<< std::endl;	  
        }//loop over channels
	
	
        //correct trigger time for propegation delay
        auto const& adsGeo = adGeo.SensitiveVolume(adsid_max); //trigger stripi

	// Timing determination
	// T1 = t_hit + x/v_propagation;
	// T1 = t_hit + (L-x)/v_propagation;
	//double deltat = (adsGeo.Length()-2*hitpos[2])/fPropDelay;
         
	// ttrigs[layer].push_back(data->fTs0);// - adsGeo.HalfLength()*fPropDelay);
        ttrigs.push_back(data->fTs0 - uint64_t(adsGeo.HalfLength()*fPropDelay));
	tpos.push_back(postrig);
	ntrig++;
	    
	if(fVerbose)
	  std::cout <<  "raw T0: " << data->fTs0 << " ,half lenth : \t" << adsGeo.HalfLength() << " ,delay: \t" <<fPropDelay 
		    << " ,corrected time: " << data->fTs0 - uint64_t(adsGeo.HalfLength()*fPropDelay) << std::endl;
    }//loop over FEBs
    

    
    for(auto const&  chan2pos: chantopos) {
      int moduleid =-5;
      for (auto const&  chan2posvec: chan2pos.second){
	moduleid  = fCrtutils->MacToAuxDetID(chan2pos.first, chan2posvec.first); //module ID
      //      if (moduleid )
      if(fVerbose)
	std::cout<< "moduleid: " << moduleid << " ,macs "<< (int)chan2pos.first  //<< std::endl;
		 << " ,chal "<<chan2posvec.first 
		 << " ,  position "<<chan2posvec.second[0] << std::endl;
      }
    }
    
    //if(fVerbose)
    //std::cout <<  "size of nfebs : \t"  <<nfebs.size() <<  " size of mac5 : \t"  << mac5.size()<< std::endl;

    if(fVerbose)
      std::cout <<  "size of layer ID : \t"  <<layID.size() << std::endl;


    
    //no channels above threshold or no intermodule coincidences? return empty hit
    //    if(nabove==0 || layID.size()!=2) {
    if(nabove==0 || layID.size() < 2) {
      if(nabove==0 && fVerbose) std::cout << "no channels above threshold!" << std::endl;
      if(layID.size()<2 && fVerbose) std::cout << "no coincidence found" << std::endl;
      return FillCRTHit({},{},0,0,0,0,0,0,0,0,0,0,"");
    }
    

    //-----------------------------------------------------------
    // Matching between both end readouts for full length module
    //-----------------------------------------------------------
    //1. The axis you are measuring on, with the correct sign (d)
    //2. The centre of the detector (c)
    //3. The absolute position is c + d · (x - L/2), 
    //   that is a displacement from the centre c by a distance x - L/2 (relative to the centre position L/2).
    //   we are determining along the z which is here d = z (aligned with the long axis of the ICARUS detector).
    //   The centre of the module you are reading is adsGeo.GetCenter()
    //   so adsGeo.GetCenter() + geo::Zaxis() * (x - adsGeo.HalfLength()) may work.
    //   We can define a different x, from the middle of the detector:
    
    // T1 = t_hit + (L/2 + x)/v_propagation;
    // T2 = t_hit + (L/2 - x)/v_propagation;
    // double deltat(T1-T2) = 2*x/v_propagation;
    //---------------------------------------------------------------------
    //std::ofstream o(filename.c_str());
    
    geo::Point_t posA; geo::Point_t posB;
    bool layer1 = false; 
    bool layer2 = false;

    geo::Point_t crossfebpos;

    uint64_t t0_1=-5; uint64_t t0_2=-5;
    uint64_t t1_1=-5; uint64_t t1_2=-5;
    uint64_t t2_1=-5; uint64_t t2_2=-5;
    int mac5_1=-5; int mac5_2 = -5;
    //float halflength = 400.;
    geo::Point_t center;

    for(auto const&  infn: informationA) {    
      auto i = &infn - informationA.data();
      auto const& adsGeo = adGeo.SensitiveVolume(infn.strip); //trigger stripi
      //mac5_1 = (int)infn.mac5s;
      //t0_1 = (long int)infn.t0;
      center = adsGeo.GetCenter();

      if(fVerbose)
	std::cout<< "A type ----------> time: " << (long long int)infn.t0 << " ,macs "<< (int)infn.mac5s  //<< std::endl;
                 << " ,chal "<< infn.channel
                 << " ,  position "<<infn.pos[2] << std::endl;
    
      if ((int)infn.mac5s != (int)informationA[i+1].mac5s and i < (int)informationA.size()-1){
	layer1 = true;

	if ((int)infn.mac5s % 2 == 0) t1_1 = infn.t0;
	else t1_1 = informationA[i+1].t0;
	if ((int)informationA[i+1].mac5s % 2 != 0) t1_2 = informationA[i+1].t0;
	else t1_2 = infn.t0;
	if(fVerbose)
	  std::cout<< "t1: " << t1_1 << ", t2:"<< t1_2 << ", deltat : "<< int64_t(t1_1 - t1_2) << std::endl;

	//if (foutCSVFile)
	//if (!filecsv.is_open()) throw cet::exception("CRTHitRecoAlg") << "Failed to create CSV file!\n";
	//std::cout<<"line 775: "<< plane << "\t"<<  int64_t(t1_1 - t1_2) << "\n";
	//if (foutCSVFile)
	//if (!filecsv) throw cet::exception("CRTHitRecoAlg") << "CSV file is in a bad state before writing!\n"; 
	if (foutCSVFile) filecsv << plane << "\t"<<  int64_t(t1_1 - t1_2) << "\n";
	//if (foutCSVFile)
	//if (!filecsv) throw cet::exception("CRTHitRecoAlg") << "CSV file is in a bad state after writing!\n"; 
	//std::cout<<"line 777: "<< plane << "\t"<<  int64_t(t1_1 - t1_2) << "\n";
 	float zaxixpos = 0.5*(int64_t(t1_1 - t1_2)/fPropDelay);
	//float zaxixpos = 0.5*(std::llabs((long int)infn.t0 - (long int)informationA[i+1].t0)*fPropDelay);
	//float zaxixpos = 0.5*((long int)infn.t0 - (long int)informationA[i+1].t0)/fPropDelay;
	/*
	if(fVerbose)
	  std::cout<< i << " ,1st mac5: "<< (int)infn.mac5s << " 1st time: " << (long int)infn.t0 
		   << " ,2nd mac5: "<<(int)informationA[i+1].mac5s << ", 2nd time " << (long int)informationA[i+1].t0 << " , deltaT: "
		   << std::llabs((long int)infn.t0 - (long int)informationA[i+1].t0) << " , length: " << adsGeo.Length() 
		   << " ,propagation delay: "<< fPropDelay << ", pos z: "
		   << 0.5*(adsGeo.Length()-(std::llabs((long int)infn.t0 - (long int)informationA[i+1].t0)*fPropDelay)) << std::endl;
	*/


	/*	if (zaxixpos > 0){
	  posA = adsGeo.GetCenter() + geo::Zaxis() * (zaxixpos - adsGeo.HalfLength());
	  std::cout<< "posA (>0): "<< posA<< std::endl;	
	}else if (zaxixpos < 0){
	  posA = adsGeo.GetCenter() + geo::Zaxis() * (zaxixpos + adsGeo.HalfLength());
	  std::cout<< "posA (<0): "<< posA<< std::endl;
	}else {
	*/
	posA = adsGeo.GetCenter() + geo::Zaxis() * zaxixpos;
	if (fVerbose) std::cout<< "posA (==0): "<< posA<< std::endl;
	  //}

	if(fVerbose)
	  std::cout<< i << " ,1st mac5: "<< (int)infn.mac5s << " 1st time: " << (long long int)infn.t0 
		   << " ,2nd mac5: "<<(int)informationA[i+1].mac5s << ", 2nd time " << (long long int)informationA[i+1].t0 << " , deltaT: "
		   << int64_t(t1_1 - t1_2) << " , length: " << adsGeo.Length()
	    //	   << (long int)infn.t0 - (long int)informationA[i+1].t0 << " , length: " << adsGeo.Length() 
		   << " ,propagation delay: "<< fPropDelay << " , pos z: " << zaxixpos
	    //<< 0.5*((long long int)infn.t0 - (long long int)informationA[i+1].t0)/fPropDelay 
		   << " , center: " << adsGeo.GetCenter() << " , zaxis: "<< geo::Zaxis() <<  " , half length:  " << adsGeo.HalfLength()
		   << " , actual pos w.rt. z: " << posA << std::endl;
	    //<< adsGeo.GetCenter() + geo::Zaxis() * (zaxixpos - adsGeo.HalfLength()) << std::endl;
      }
    }

    for(auto const&  infn: informationB) {

     
      if(fVerbose)
	std::cout<< " B type ----------> time: " << infn.t0 << " ,macs "<< (int)infn.mac5s  //<< std::endl;
                 << " ,chal "<< infn.channel
                 << " ,  position "<<infn.pos[2] << std::endl;

      auto i = &infn - informationB.data();
      auto const& adsGeo = adGeo.SensitiveVolume(infn.strip); //trigger stripi

      // Timing determination
      // T1 = t_hit + x/v_propagation;
      // T1 = t_hit + (L-x)/v_propagation;
      //double deltat = (adsGeo.Length()-2*hitpos[2])/fPropDelay;
      if ((int)infn.mac5s != (int)informationB[i+1].mac5s and i < (int)informationB.size()-1){
	//float zaxixpos = 0.5*(std::llabs((long int)infn.t0 - (long int)informationB[i+1].t0)*fPropDelay);
	//      float zaxixpos = 0.5*((long int)infn.t0 - (long int)informationB[i+1].t0)/fPropDelay;



      layer2 = true;
      mac5_2 = (int)infn.mac5s;
      t0_2 = (long int)infn.t0;
      std::cout << "t0_2 (line 849) : " << t0_2 << "\n";
      std::cout << "t0 (64?): " << uint64_t(infn.t0) << "\n";

      if ((int)infn.mac5s % 2 == 0) t2_1 = infn.t0;
      else t2_1 = informationB[i+1].t0;
      if ((int)informationB[i+1].mac5s % 2 != 0) t2_2 = informationB[i+1].t0;
      else t2_2 = infn.t0;

      if(fVerbose)
	std::cout<< "t1: " << t2_1 <<", t2:"<< t2_2 << ", deltat : "<< int64_t(t2_1 - t2_2) << std::endl;
      if (foutCSVFile) filecsv << plane << "\t"<<  int64_t(t2_1 - t2_2) << "\n";
      float zaxixpos = 0.5*(int64_t(t2_1 - t2_2)/fPropDelay);

      //if (foutCSVFile) filecsv.close();

      /*

      if(fVerbose)
	std::cout<< i << " ,1st mac5: "<< (int)infn.mac5s << " 1st time: " << (long int)infn.t0
                 << " ,2nd mac5: "<<(int)informationB[i+1].mac5s << ", 2nd time " << (long int)informationB[i+1].t0 << " , deltaT: "
                 << std::llabs((long int)infn.t0 - (long int)informationB[i+1].t0) << " , length: " << adsGeo.Length()
                 << " ,propagation delay: "<< fPropDelay << ", pos z: "
                 << 0.5*(adsGeo.Length()-(std::llabs((long int)infn.t0 - (long int)informationB[i+1].t0)*fPropDelay)) << std::endl;
      */

      /*
      if (zaxixpos > 0){
	posB = adsGeo.GetCenter() + geo::Zaxis() * (zaxixpos - adsGeo.HalfLength());
	std::cout<< "posB (>0): "<< posB<< std::endl;
      }else if (zaxixpos < 0){
	posB = adsGeo.GetCenter() + geo::Zaxis() * (zaxixpos + adsGeo.HalfLength());
	std::cout<< "posB (< 0): "<< posB<< std::endl;
      }else {
      */
	posB = adsGeo.GetCenter() + geo::Zaxis() * zaxixpos;
	if (fVerbose) std::cout<< "posB (== 0): "<< posB<< std::endl;
	//}
      //posB = adsGeo.GetCenter() + geo::Zaxis() * (zaxixpos - adsGeo.HalfLength());

      if(fVerbose)
	std::cout<< i << " ,1st mac5: "<< (int)infn.mac5s << " 1st time: " << (long long int)infn.t0
		 << " ,2nd mac5: "<<(int)informationB[i+1].mac5s << ", 2nd time " << (long long int)informationB[i+1].t0 << " , deltaT: "
		 << int64_t(t2_1 - t2_2) << " , length: " << adsGeo.Length()
	  //	 << (long int)infn.t0 - (long int)informationB[i+1].t0 << " , length: " << adsGeo.Length()
		 << " ,propagation delay: "<< fPropDelay << " , pos z: " << zaxixpos
	  //	 << 0.5*((long long int)infn.t0 - (long long int)informationB[i+1].t0)/fPropDelay
		 << " , center: " << adsGeo.GetCenter() << " , zaxis: "<< geo::Zaxis() <<  " , half length:  " << adsGeo.HalfLength()
		 << " , actual pos w.rt. z: " << posB << std::endl;
	  //<< adsGeo.GetCenter() + geo::Zaxis() * (zaxixpos - adsGeo.HalfLength()) << std::endl;

      }
      
    }


    for (auto const&  infna: informationA) {
      for (auto const&  infnb: informationB) {
	if(fVerbose)
	  std::cout<< "macs a "<< (int)infna.mac5s  //<< std::endl;
		   << " ,chal a "<< infna.channel
		   << " ,macs b "<< (int)infnb.mac5s
		   << " ,chal b "<< infnb.channel
		   << " ,  position a ("<<infna.pos[0] << ", " << infna.pos[1] << " , " << infna.pos[2]<< ")"
		   << " ,  position b ("<<infnb.pos[0] << ", " << infnb.pos[1] << " , " << infnb.pos[2]<< ")" 
		   << " , slope: "<< (infnb.pos[1] - infna.pos[1])/(infnb.pos[0] - infna.pos[0])<< std::endl;
      }
    }


    int crossfeb = std::abs(mac5_1 - mac5_2);


    // side crt and match the both layers
    if (layer1 && layer2 && region!="South" && region!="North" ){//&& nx==4){
      float avg = 0.5*(posA.Z() + posB.Z());
      hitpos.SetZ(avg);
      //hitpos.SetX(hitpos.X()*1.0/nx);
      //      hitpos.SetY(hitpos.Y()*1.0/nx);
      hitpos.SetX(hitpos.X()*1.0/petot);
      hitpos.SetY(hitpos.Y()*1.0/petot);
      if(fVerbose)
	std::cout << "z position in layer 1: "<< posA.Z() << " and in layer 2 "<< posB.Z() 
		  << " average is "<< (posA.Z()+ posB.Z())/2. << " ,hitpos z " << hitpos[2] << std::endl;

      
    }else if ((int)informationA.size()==1 and (int)informationB.size()==1
	      and (crossfeb == 7 or crossfeb == 5) and 
	      region!="South" && region!="North"){
      //      int z_pos =  0.5*(std::llabs(t0_1 - t0_2)*fPropDelay);
      int z_pos =  int64_t(t0_1 - t0_2)/(uint64_t(2*fPropDelay));
      /*
      if (z_pos > 0){
	crossfebpos =  center + geo::Zaxis()*(z_pos - halflength);
	
      }else if (z_pos < 0){
	crossfebpos =  center + geo::Zaxis()*(z_pos + halflength);

      }else {
      */
      crossfebpos =  center + geo::Zaxis()*z_pos;
	//}

      // crossfebpos =  center + geo::Zaxis()*(z_pos - halflength);
      hitpos.SetZ(crossfebpos.Z());
      // hitpos.SetX(hitpos.X()*1.0/nx);
      //      hitpos.SetY(hitpos.Y()*1.0/nx);
      hitpos.SetX(hitpos.X()*1.0/petot);
      hitpos.SetY(hitpos.Y()*1.0/petot);
      if(fVerbose)
	std::cout << "hello hi namaskar,  hitpos z " << hitpos[2] << std::endl;
      // side crt and only single layer match
    }else if (layer1 && region!="South" && region!="North"){// && nx==1){
      hitpos.SetZ(posA.Z());
      //      hitpos.SetX(hitpos.X()*1.0/nx);
      //hitpos.SetY(hitpos.Y()*1.0/nx);
      hitpos.SetX(hitpos.X()*1.0/petot);
      hitpos.SetY(hitpos.Y()*1.0/petot);
      if(fVerbose)
	std::cout << " same layer coincidence:  z position in layer 1: "<< posA.Z() << " ,hitpos z " << hitpos[2] << std::endl;
      
      // side crt and only single layer match
    }else if (layer2 && region!="South" && region!="North" ){//&& nx==1){
      hitpos.SetZ(posB.Z());
      //hitpos.SetX(hitpos.X()*1.0/nx);
      //hitpos.SetY(hitpos.Y()*1.0/nx);
      hitpos.SetX(hitpos.X()*1.0/petot);
      hitpos.SetY(hitpos.Y()*1.0/petot);
      if(fVerbose)
	std::cout << " same layer coincidence: z position in layer 2 "<< posB.Z() << " ,hitpos z " << hitpos[2] << std::endl;
      
    }else if (region!="South" && region!="North" ){//&& nx==2){
      //  hitpos*=1.0/nx;
      hitpos.SetX(hitpos.X()*1.0/petot);
      hitpos.SetY(hitpos.Y()*1.0/petot);
      hitpos.SetZ(hitpos.Z()*1.0/petot);
      if (fVerbose) std::cout << " In side CRTs [E/W] x: \t"<< hitpos[0] <<" ,y: \t" << hitpos[1]  <<" ,z: \t" << hitpos[2]<< std::endl;
    }/*else {
      hitpos.SetX(-99999);
      hitpos.SetY(-99999);
      hitpos.SetZ(-99999);
    }*/
    

    //finish averaging and fill hit point array
    if(region=="South") {
      
      hitpos.SetX(hitpos.X()*1.0/pex);
      hitpos.SetY(hitpos.Y()*1.0/pey);
      hitpos.SetZ(hitpos.Z()*1.0/petot);
      /*
	hitpos.SetX(hitpos.X()/nx);
	hitpos.SetY(hitpos.Y()/ny);
	hitpos.SetZ(hitpos.Z()/nz);
        */
	// }else
      //hitpos*=1.0/petot; //hit position weighted by deposited charge
   
    }else if (region=="North"){
      hitpos*=1.0/petot;
      //hitpos*=1.0/nz;

      //}else if (region!="South" && region!="North"){
      // hitpos.SetX(hitpos.X()*1.0/petot);
      //hitpos.SetY(hitpos.Y()*1.0/petot);
    }
   
    //else
    //hitpos*=1.0/petot; //hit position weighted by deposited charge
    //hitpos*=1.0; //hit position weighted by deposited charge

    hitpoint[0] = hitpos.X();
    hitpoint[1] = hitpos.Y();
    hitpoint[2] = hitpos.Z();
    //if (fVerbose){
    if (region=="South" && hitpoint[0] >= 366. && hitpoint[1] > 200. && fVerbose)
      std::cout << "I am looking for south wall :   macs " << (int)macs.back() << " x: \t"<< hitpoint[0] 
		<<" ,y: \t" << hitpoint[1]  <<" ,z: \t" << hitpoint[2] << std::endl;
    //}
    if (fVerbose){
      if (region=="North") std::cout << "north wall x: \t"<< hitpoint[0] <<" ,y: \t" << hitpoint[1]  <<" ,z: \t" << hitpoint[2]<< std::endl;
    } 
    if (fVerbose) std::cout << " nx: \t"<< nx <<" ,ny: \t" << ny  <<" ,nz: \t" << nz<< std::endl;
    if (fVerbose) std::cout << " x: \t"<< hitpoint[0] <<" ,y: \t" << hitpoint[1]  <<" ,z: \t" << hitpoint[2]<< std::endl;
    
    //time stamp averaged over all FEBs
    uint64_t thit = 0.;//, thit_0 = 0., thit_1 = 0.;
    //uint64_t thit_64 = 0.;
    //for(double const t : ttrigs[0]) 
    //    thit_0 += t;
    for(uint64_t const t : ttrigs){
      std::cout << "t in ttrigs " << (long long int)t << "\n";
      thit += t;
      //thit +=   uint64_t(t);
      std::cout << "in loop      thit: " << thit << "\n";
    }
      //thit_0
    thit*=1.0/uint64_t(ttrigs.size());
    std::cout << "thit: " << (long long int)thit << "\n";

    /*for(uint64_t const t_64 : ttrigs){
      std::cout << "t_64 in ttrigs " << t_64 << "\n";
      thit_64 +=   uint64_t(t_64);
      //thit_0                                                                                                                            
    }
      thit_64*=1.0/uint64_t(ttrigs.size());
      std::cout << "thit_64: " << thit_64 << "--------\n";*/
    if (fVerbose) std::cout << " <time>: \t"<< thit <<" size ttrig: \t" <<   ttrigs.size()<< std::endl;    

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
    CRTHit hit = FillCRTHit(macs,pesmap,petot,thit,thit,plane,hitpoint[0],hitpointerr[0],
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
