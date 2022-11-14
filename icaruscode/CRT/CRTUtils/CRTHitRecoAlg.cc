#include "CRTHitRecoAlg.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include <algorithm>
using namespace icarus::crt;
/*
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
*/
CRTHitRecoAlg::CRTHitRecoAlg(const fhicl::ParameterSet& pset){
  this->reconfigure(pset);
  fChannelMap = art::ServiceHandle<icarusDB::IICARUSChannelMap const>{}.get();
  fGeometryService  = lar::providerFrom<geo::Geometry>();
  fCrtutils = new CRTCommonUtils();
  return;
}

CRTHitRecoAlg::CRTHitRecoAlg() = default;
/*
//---------------------------------------------------------------------
void CRTHitRecoAlg::reconfigure(const Config& config){
    fVerbose = config.Verbose();
    fUseReadoutWindow = config.UseReadoutWindow(); 
    fQPed = config.QPed();
    fQSlope = config.QSlope();
    fPropDelay = config.PropDelay(); 
    fPEThresh = config.PEThresh();
    ftopGain = config.topGain();
    ftopPed = config.topPed();
    fSiPMtoFEBdelay = config.SiPMtoFEBdelay();
    fCoinWindow = config.CoinWindow();
    fCrtWindow = config.CrtWindow();
    foutCSVFile = config.outCSVFile();
    fCSVFile = config.CSVFile();
    fData = config.Data();
    if (foutCSVFile) 
      filecsv.open(fCSVFile.c_str());
    return;
}
*/
void CRTHitRecoAlg::reconfigure(const fhicl::ParameterSet& pset){
  fVerbose          = pset.get<bool>("Verbose", false);
  fUseReadoutWindow = pset.get<bool>("UseReadoutWindow", false);
  fQPed             = pset.get<double>("QPed", 0.);
  fQSlope           = pset.get<double>("QSlope", 0.);
  fPropDelay        = pset.get<double>("PropDelay",0.);
  fPEThresh         = pset.get<double>("PEThresh",0.);
  ftopGain          = pset.get<double>("topGain",0.);
  ftopPed           = pset.get<double>("topPed",0.);
  fSiPMtoFEBdelay   = pset.get<uint64_t>("SiPMtoFEBdelay", 0.);
  fCoinWindow       = pset.get<uint64_t>("CoinWindow", 0.);
  fCrtWindow        = pset.get<uint64_t>("CrtWindow", 0.);
  foutCSVFile       = pset.get<bool>("outCSVFile", false);
  fCSVFile          = pset.get<std::string>("CSVFile", "");
  fData             = pset.get<bool>("Data", false);
  if (!fCSVFile.empty())  filecsv.open(fCSVFile);
  {
    std::vector<std::vector<int32_t> > T1delays =  pset.get<std::vector<std::vector<int32_t> > >("FEB_T1delay_side");
    std::vector<std::vector<int32_t> > T0delays =  pset.get<std::vector<std::vector<int32_t> > >("FEB_T0delay_side");
    for(auto & feb : T1delays) {
      int32_t & mac = feb[0];
      int32_t & d   = feb[1];
      FEB_T1delay_side[mac] = d;
    }
    for(auto & feb : T0delays) {
      int32_t & mac = feb[0];
      int32_t & d   = feb[1];
      FEB_T0delay_side[mac] = d;
    }
  }

  return;
}

//---------------------------------------------------------------------------------------
vector<art::Ptr<CRTData>> CRTHitRecoAlg::PreselectCRTData(const vector<art::Ptr<CRTData>> &crtList, uint64_t trigger_timestamp){
  if (fVerbose)  mf::LogInfo("CRTHitRecoAlg: ") << "In total " << crtList.size() << " CRTData found in an event" << '\n';
  vector<art::Ptr<CRTData>> crtdata;
  bool presel = false;

  for (size_t febdat_i=0; febdat_i<crtList.size(); febdat_i++) {
    
    uint8_t mac = crtList[febdat_i]->fMac5;
    int adid    = fCrtutils->MacToAuxDetID(mac,0);
    char type   = fCrtutils->GetAuxDetType(adid);

    /// Looking for data within +/- 3ms within trigger time stamp
    /// Here t0 - trigger time -ve
    if (fData && (std::fabs(int64_t(crtList[febdat_i]->fTs0 - trigger_timestamp)) > fCrtWindow)) continue;
    
    if ( type == 'm'){
      for(int chan=0; chan<32; chan++) {
	std::pair<double,double> const chg_cal = fChannelMap->getSideCRTCalibrationMap((int)crtList[febdat_i]->fMac5,chan);
	float pe = (crtList[febdat_i]->fAdc[chan]-chg_cal.second)/chg_cal.first;
	if(pe<=fPEThresh && !crtList[febdat_i]->IsReference_TS1()) continue; // filter out low PE and non-T1 ref values
	presel = true;
	if(fVerbose)
	  mf::LogInfo("CRTHitRecoAlg: ") << "\nfebP (mac5, channel, gain, pedestal, adc, pe) = (" << (int)crtList[febdat_i]->fMac5 << ", " << chan << ", " 
					 << chg_cal.first << ", " << chg_cal.second << "," << crtList[febdat_i]->fAdc[chan] << "," << pe << ")\n";
      }
    }else if ( type == 'c' ) {
      for(int chan=0; chan<32; chan++) {
	//float pe = (crtList[febdat_i]->fAdc[chan]-fQPed)/fQSlope;
	//if(pe<=fPEThresh) continue;
	presel = true;
      }
    }else if ( type == 'd'){
      for(int chan=0; chan<64; chan++) {
	float pe = (crtList[febdat_i]->fAdc[chan]-fQPed)/fQSlope;
	if(pe<=fPEThresh) continue;
	presel = true;
      }
    }
  
    if (presel) crtdata.push_back(crtList[febdat_i]);
    presel = false;
    
  }
  mf::LogInfo("CRTHitRecoAlg:") << "Found " << crtdata.size() << " after preselection "<< '\n';
  //std::cout<<trigger_timestamp<<std::endl;
  return crtdata;  
}

//---------------------------------------------------------------------------------------
vector<pair<sbn::crt::CRTHit, vector<int>>> CRTHitRecoAlg::CreateCRTHits(vector<art::Ptr<CRTData>> crtList, uint64_t trigger_timestamp) {
  
  vector<pair<CRTHit, vector<int>>> returnHits;
  vector<int> dataIds;
  
    uint16_t nMissC = 0, nMissD = 0, nMissM = 0, nHitC = 0, nHitD = 0, nHitM = 0;
    if (fVerbose) mf::LogInfo("CRTHitRecoAlg: ") << "Found " << crtList.size() << " FEB events" << '\n';

    map<string,int> regCounts;
    std::set<string> regs;
    map<string,vector<size_t>> sideRegionToIndices;
    
    // sort by the time 
    std::sort(crtList.begin(), crtList.end(), compareBytime);        

    //Load Delays map for Top CRT
    CRT_delay_map const FEB_delay_map = LoadFEBMap();
    //vectors to store CRT Resets. Note: West and North Side CRTs are on the West timing rack, East and South Side CRTs are on East timing rack
    std::vector<std::pair<int,ULong64_t>> CRTReset_top;
    std::vector<std::pair<int,ULong64_t>> CRTReset_side_west;
    std::vector<std::pair<int,ULong64_t>> CRTReset_side_east;
    ULong64_t TriggerArray_top[305]={0};
    ULong64_t TriggerArray_side_west[305]={0};
    ULong64_t TriggerArray_side_east[305]={0};
    
   for (size_t crtdat_i=0; crtdat_i<crtList.size(); crtdat_i++) {
	uint8_t mac = crtList[crtdat_i]->fMac5;
	int adid  = fCrtutils->MacToAuxDetID(mac,0);
	char type = fCrtutils->GetAuxDetType(adid);
	string region = fCrtutils->GetAuxDetRegion(adid);
        int plane =fCrtutils->AuxDetRegionNameToNum(region);
	//For the time being, Only Top CRT delays are loaded, nothing to do for Side CRT yet
	if (type == 'c' && crtList[crtdat_i]->IsReference_TS1()) {
	    ULong64_t Ts0T1ResetEvent = crtList[crtdat_i]->fTs0 + FEB_delay_map.at((int)mac+73).T0_delay - FEB_delay_map.at((int)mac+73).T1_delay;
	    TriggerArray_top[(int) mac]=Ts0T1ResetEvent;
	    CRTReset_top.emplace_back((int) mac,Ts0T1ResetEvent);
	}
	if (type == 'm' && crtList[crtdat_i]->IsReference_TS1()) {
	  try
	    {
	      ULong64_t Ts0T1ResetEvent = crtList[crtdat_i]->fTs0 + FEB_T0delay_side.at(mac) - FEB_T1delay_side.at(mac);
	      if(plane == 40 || plane == 41 || plane == 42 || plane == 47){//CRT T1 reset on West/North 
		TriggerArray_side_west[(int) mac]=Ts0T1ResetEvent;
                CRTReset_side_west.emplace_back((int) mac,Ts0T1ResetEvent);
	      }
	      else{
		TriggerArray_side_east[(int) mac]=Ts0T1ResetEvent;
                CRTReset_side_east.emplace_back((int) mac,Ts0T1ResetEvent);
	      }
	    }
	  catch(std::out_of_range& e)
	    {
	      std::cout << "not found in the FEB_delay array!!! Please update FEB_delay FHiCL file \n";
            }
	}
	
   }
   //std::cout << "------\nCreateCRTHits::Trigger timestamp = "<<trigger_timestamp<<"\n------\n";
   //Global Trigger Calculation for Top CRT
   const int trigger_offset= 60; //Average distance between Global Trigger and Trigger_timestamp (ns)
   ULong64_t GlobalTrigger_top= trigger_timestamp;
   if (!CRTReset_top.empty()) GlobalTrigger_top = GetMode(CRTReset_top);
   //Add average difference between trigger_timestamp and Global trigger
   else GlobalTrigger_top=GlobalTrigger_top-trigger_offset;// In this event, the T1 Reset was probably "vetoed" by the T0 Reset
   for (int i=0; i<305; i++){
     if (TriggerArray_top[i]==0) TriggerArray_top[i]=GlobalTrigger_top;
   }
   //std::cout<<"------\nCreateCRTHits::Mode Top CRT Global Trigger ="<<GlobalTrigger_top<<"\n------\n";

   //Global Trigger Calculation for West/North Side CRTs
   ULong64_t GlobalTrigger_side_west= trigger_timestamp;
   if (!CRTReset_side_west.empty()) GlobalTrigger_side_west = GetMode(CRTReset_side_west);
   //Add average difference between trigger_timestamp and Global trigger    
   else GlobalTrigger_side_west=GlobalTrigger_side_west-trigger_offset;// In this event, the T1 Reset was probably "vetoed" by the T0 Reset
   for (int i=0; i<305; i++){
     if (TriggerArray_side_west[i]==0) TriggerArray_side_west[i]=GlobalTrigger_side_west;
   }
   //std::cout<<"------\nCreateCRTHits::Mode Side West CRT Global Trigger ="<<GlobalTrigger_side_west<<"\n------\n";

   //Global Trigger Calculation for East/South Side CRTs
   ULong64_t GlobalTrigger_side_east= trigger_timestamp;
   if (!CRTReset_side_east.empty()) GlobalTrigger_side_east = GetMode(CRTReset_side_east);
   //Add average difference between trigger_timestamp and Global trigger
   else GlobalTrigger_side_east=GlobalTrigger_side_east-trigger_offset;// In this event, the T1 Reset was probably "vetoed" by the T0 Reset
   for (int i=0; i<305; i++){
     if (TriggerArray_side_east[i]==0) TriggerArray_side_east[i]=GlobalTrigger_side_east;
   }
   //std::cout<<"------\nCreateCRTHits::Mode Side East CRT Global Trigger ="<<GlobalTrigger_side_east<<"\n------\n";


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
            hit = MakeTopHit(crtList[febdat_i], TriggerArray_top);
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
 
        if ( type == 'm' )
            sideRegionToIndices[region].push_back(febdat_i);

    }//loop over CRTData products

    vector<size_t> unusedDataIndex;
    for(auto const& regIndices : sideRegionToIndices) {
      
      if(fVerbose) 
	mf::LogInfo("CRTHitRecoAlg: ") << "searching for side CRT hits in region, " << regIndices.first << '\n';
      
      vector<size_t> indices = regIndices.second;
      
      if(fVerbose)
	mf::LogInfo("CRTHitRecoAlg: ") << "number of hits associated to this region : " << indices.size() << '\n';
      
      for(size_t index_i=0; index_i < indices.size(); index_i++) {
          
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
	  // in coincidence
	  //      if(crtList[indices[index_j]]->fTs0 <= crtList[indices[index_i]]->fTs0 + fCoinWindow) {
	  //            if(std::llabs(crtList[indices[index_j]]->fTs0 - crtList[indices[index_i]]->fTs0) < fCoinWindow) {
	  if( (crtList[indices[index_j]]->fTs0>=crtList[indices[index_i]]->fTs0 && 
	       (crtList[indices[index_j]]->fTs0 - crtList[indices[index_i]]->fTs0) < fCoinWindow) ||
	      (crtList[indices[index_j]]->fTs0<crtList[indices[index_i]]->fTs0 && 
	       (crtList[indices[index_i]]->fTs0 - crtList[indices[index_j]]->fTs0) < fCoinWindow)) {
	    if(fVerbose)
	      mf::LogInfo("CRTHitRecoAlg: ") <<  " in coincidence: i \t " << index_i << " ,j: \t" << index_j <<",i mac: \t" 
			<< (int)crtList[indices[index_i]]->fMac5 << ", j mac: \t" <<(int)crtList[indices[index_j]]->fMac5<< '\n';
	    
	    coinData.push_back(crtList[indices[index_j]]);
	    dataIds.push_back(indices[index_j]);
	  }
          
	  //out of coinWindow
	  
	  if  ((crtList[indices[index_j]]->fTs0 - crtList[indices[index_i]]->fTs0) > fCoinWindow || index_j==indices.size()-1){
	    if(fVerbose)
	      mf::LogInfo("CRTHitRecoAlg: ") <<  "out of coincidence  " << index_j << "\t" << indices.size() <<"\t" <<indices.size()-1
			<< " data products..." << '\n';
	    if(fVerbose)
	      mf::LogInfo("CRTHitRecoAlg: ") << "attempting to produce MINOS hit from " << coinData.size() 
			<< " data products..." << '\n';
	    uint8_t imac = (int)crtList[indices[index_i]]->fMac5;
            int adid  = fCrtutils->MacToAuxDetID(imac,0);
            string region = fCrtutils->GetAuxDetRegion(adid);
            int plane =fCrtutils->AuxDetRegionNameToNum(region);
	    //CRTHit hit = MakeSideHit(coinData, TriggerArray);
            CRTHit hit;
	    if(plane == 40 || plane == 41 || plane == 42 || plane == 47){//CRT T1 reset on West/North
	      hit = MakeSideHit(coinData, TriggerArray_side_west);
	    }
	    else{
	      hit = MakeSideHit(coinData, TriggerArray_side_east);
	    }
	      
	    if(IsEmptyHit(hit)){
	      unusedDataIndex.push_back(indices[index_i]);
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
                            float peshit, uint64_t time0, Long64_t time1, int plane, 
                            double x, double ex, double y, double ey, double z, double ez, string tagger){
    CRTHit crtHit;
    crtHit.feb_id      = tfeb_id;
    crtHit.pesmap      = tpesmap;
    crtHit.peshit      = peshit;
    crtHit.ts0_s_corr  = time0 / 1'000'000'000; 
    crtHit.ts0_ns      = time0 % 1'000'000'000;
    crtHit.ts0_ns_corr = time0; 
    crtHit.ts1_ns      = time1 /*% 1'000'000'000*/; //TODO: Update the CRTHit data product /sbnobj/common/CRT . Discussion with SBND people needed
    crtHit.ts0_s       = time0 / 1'000'000'000;//'
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
sbn::crt::CRTHit CRTHitRecoAlg::MakeTopHit(art::Ptr<CRTData> data, ULong64_t GlobalTrigger[305]){

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
    double sum=0;
    for(int chan=0; chan<32; chan++) {
	sum=sum+data->fAdc[chan];
        float pe = (data->fAdc[chan]-ftopPed)/ftopGain;
//      if(pe<=fPEThresh) continue;
        nabove++;
        int adsid = fCrtutils->ChannelToAuxDetSensitiveID(mac,chan);
        petot += pe;
        pesmap[mac].push_back(std::make_pair(chan,pe));

        //TVector3 postmp = fCrtutils->ChanToLocalCoords(mac,chan);
        //strip along z-direction
        if(adsid >= 0 && adsid < 8){
            //hitpos.SetX(pe*postmp.X()+hitpos.X());
            //hitpos.SetX(postmp.X()+hitpos.X());
            if(pe>pemaxx){
                pemaxx = pe;
                maxx = adsid;
            }
            findx = true;   
        }
        //strip along x-direction
        else if(adsid >= 8 && adsid < 16 ){
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
    Long64_t thit1 = data->fTs1;
    thit -= fSiPMtoFEBdelay;
    thit1 -= fSiPMtoFEBdelay;

    //92.0 is the middle of one of the Top CRT modules (each of them is 184 cm)    
    //TO DO: Move hardcoded numbers to parameter fcl files.
    //Another possibility is using object values from GDML, but I (Francesco Poppi) found weird numbers some months ago and needed to double check.
    double const corrPos = std::max(-hitpos.X(), hitpos.Z());
    uint64_t const corr = (uint64_t)round(abs((92.0+corrPos)*fPropDelay));
    thit -= corr;
    thit1 -= corr;

    adGeo.LocalToWorld(hitlocal,hitpoint); //tranform from module to world coords

    hitpointerr[0] = adsGeo.HalfWidth1()*2/sqrt(12);
    hitpointerr[1] = adGeo.HalfHeight();
    hitpointerr[2] = adsGeo.HalfWidth1()*2/sqrt(12);
    //thit1 = (Long64_t)(thit-GlobalTrigger[(int)mac+73]);
    thit1 = (Long64_t)(thit-GlobalTrigger[(int)mac]);

    //Remove T1 Reset event not correctly flagged, remove T1 reset events, remove T0 reset events
    if((sum<10000 && thit1<2'001'000 && thit1>2'000'000)||data->IsReference_TS1() || data->IsReference_TS0()) return FillCRTHit({},{},0,0,0,0,0,0,0,0,0,0,"");

    CRTHit hit = FillCRTHit({mac},pesmap,petot,thit,thit1,plane,hitpoint[0],hitpointerr[0],
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
sbn::crt::CRTHit CRTHitRecoAlg::MakeSideHit(vector<art::Ptr<CRTData>> coinData, ULong64_t GlobalTrigger[305]) {

    vector<uint8_t> macs;
    map< uint8_t, vector< pair<int,float> > > pesmap;


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
    auto const& adGeo = fGeometryService->AuxDet(adid); //module
    string region = fCrtutils->GetAuxDetRegion(adid);
    int plane =fCrtutils->AuxDetRegionNameToNum(region);


    double hitpoint[3], hitpointerr[3];
    TVector3 hitpos (0.,0.,0.);


    float petot = 0., pemax = 0., pex=0., pey=0.;
    int adsid_max = -1, nabove=0, nx=0, ny=0, nz = 0, ntrig = 0;

    TVector3 postrig;

    vector<uint64_t> ttrigs;
    vector<uint64_t> t1trigs;
    vector<TVector3> tpos;
    double zmin=DBL_MAX, zmax = -DBL_MAX;
    double ymin=DBL_MAX, ymax = -DBL_MAX;
    double xmin=DBL_MAX, xmax = -DBL_MAX;
    std::vector<int> layID;
    std::vector<int> febA;
    std::vector<int> febB;


    uint64_t southt0_v = -999, southt0_h =-999;
    
    //loop over FEBs
    for(auto const& data : coinData) {

      if (adid == (int)fCrtutils->MacToAuxDetID((int)data->fMac5,0)){
	febA.push_back(data->fMac5);
      }else {
	febB.push_back(data->fMac5);
      }
    }

    if(fVerbose)
      std ::cout << "line 451: size of febA: \t" << (int)febA.size() 
		 << " size of febB: " << (int)febB.size() << '\n';
    
    
    //loop over FEBs
    for(auto const& data : coinData) {

      //if(!(region=="South")) continue;
        macs.push_back(data->fMac5);
        adid  = fCrtutils->MacToAuxDetID(macs.back(),0);


        int layer = fCrtutils->GetMINOSLayerID(adid);
        layID.push_back(layer);   

	auto idx = &data - coinData.data();
	if(fVerbose)    std :: cout << "index: " << idx <<" , mac5: " << (int)macs.back() << " ,adid: " << adid<< '\n';


	if ((int)febA.size() == 0 or (int)febB.size() == 0) continue;

	// deciding to use largest size of FEBs in the loop
	int size = 0;
	if (febA.size() < febB.size())  size = febB.size();
	else size = febA.size();


	for (int i = 0; i < size; i++){
	  if  (macs.back() == (int)febA[i]) {
	    //loop over channels

	    for(int chan=0; chan<32; chan++) {
	      std::pair<double,double> const chg_cal = fChannelMap->getSideCRTCalibrationMap((int)data->fMac5,chan);
	      float pe = (data->fAdc[chan]-chg_cal.second)/chg_cal.first;

	      //float pe = (data->fAdc[chan]-fQPed)/fQSlope;
	      if(pe<=fPEThresh) continue;
	      nabove++;
	      if(fVerbose)
		mf::LogInfo("CRTHitRecoAlg: ") << "\nfebA (mac5, channel, gain, pedestal, adc, pe) = (" << (int)data->fMac5 << ", " << chan << ", "
			  << chg_cal.first << ", " << chg_cal.second << "," << data->fAdc[chan] << "," << pe << ")\n";
	      int adsidA = fCrtutils->ChannelToAuxDetSensitiveID(macs.back(),chan);
	      TVector3 postmp = fCrtutils->ChanToWorldCoords(macs.back(),chan);

	      informationA.push_back ({macs.back(), chan, data->fTs0, postmp, adsidA});
	    }

	    if(fVerbose)
	      mf::LogInfo("CRTHitRecoAlg: ") << "looking for mac " << (int)macs.back()
			<< " and matching with febA : " << (int)febA[i]<<'\n'; 
	  }else if ( macs.back() == (int)febB[i]){

	    //loop over channels
            for(int chan=0; chan<32; chan++) {
	      std::pair<double,double> const chg_cal = fChannelMap->getSideCRTCalibrationMap((int)data->fMac5,chan);
              float pe = (data->fAdc[chan]-chg_cal.second)/chg_cal.first;
              //float pe = (data->fAdc[chan]-fQPed)/fQSlope;
              if(pe<=fPEThresh) continue;
              nabove++;
	      if(fVerbose)
		mf::LogInfo("CRTHitRecoAlg: ") << "\nfebB (mac5, channel, gain, pedestal, adc, pe) = (" << (int)data->fMac5 << ", " << chan << ", "
			  << chg_cal.first << ", " << chg_cal.second << "," << data->fAdc[chan] << "," << pe << ")\n";
	      int adsidB = fCrtutils->ChannelToAuxDetSensitiveID(macs.back(),chan);
              TVector3 postmp = fCrtutils->ChanToWorldCoords(macs.back(),chan);

              informationB.push_back ({macs.back(), chan, data->fTs0, postmp, adsidB});
            }
	    if(fVerbose)
	      mf::LogInfo("CRTHitRecoAlg: ") << "else if looking for mac "<< (int)macs.back()
			<< " and matching with febB : " << (int)febB[i]<<'\n';
	  }
	}

	if(fVerbose)
	  mf::LogInfo("CRTHitRecoAlg: ") << "In CRTHitRecoAlg::MakeSideHit 1st feb is " << (int)fCrtutils->ADToMac(adid).first 
		    << " ,2nd feb :"<< (int)fCrtutils->ADToMac(adid).second  
		    << ", time "<< data->fTs0 
                    << "  with module number " << adid << ", no. of FEB " <<'\n';

	if(fVerbose)
	  mf::LogInfo("CRTHitRecoAlg: ") << "In CRTHitRecoAlg::MakeSideHit functions mac is " << (int)macs.back()
		    << "  with module number " << adid << ", no. of FEB " <<'\n';

        //loop over channels
        for(int chan=0; chan<32; chan++) {
	  // std :: cout << "chan: ---------------- " << chan << " , "<< fCrtutils->ChannelToAuxDetSensitiveID(macs.back(),chan) << '\n';	  
	  std::pair<double,double> const chg_cal = fChannelMap->getSideCRTCalibrationMap((int)data->fMac5,chan);
	  float pe = (data->fAdc[chan]-chg_cal.second)/chg_cal.first;
	  //float pe = (data->fAdc[chan]-fQPed)/fQSlope;
	  if(pe<=fPEThresh) continue;
	  nabove++;
	  if(fVerbose)
	    mf::LogInfo("CRTHitRecoAlg: ") << "\nfebC (mac5, channel, gain, pedestal, adc, pe) = (" << (int)data->fMac5 << ", " << chan << ", "
		      << chg_cal.first << ", " << chg_cal.second << "," << data->fAdc[chan] << "," << pe << ")\n";
	  int adsid = fCrtutils->ChannelToAuxDetSensitiveID(macs.back(),chan);
	  petot += pe;
	  pesmap[macs.back()].push_back(std::make_pair(chan,pe));

	  TVector3 postmp = fCrtutils->ChanToWorldCoords(macs.back(),chan);
	    
	  if(fVerbose)
	    mf::LogInfo("CRTHitRecoAlg: ") <<  "feb: " << (int)macs.back() << " ,chan : \t" << chan
		      << " ,pe: \t"<< pe << ", adc:\t" << data->fAdc[chan] << ", time: \t"<< data->fTs0
                      << " ,x: \t"<< postmp.X() <<" ,y: \t" << postmp.Y()  <<" ,z: \t" << postmp.Z() << " petotal : "<< petot << '\n';

	  //East/West Walls (all strips along z-direction) or
	  // North/South inner walls (all strips along x-direction)
	  // All the horizontal layers measure Y first,
	  if(!(region=="South" && layer==1)) {
	    //hitpos.SetY(pe*postmp.Y()+hitpos.Y());
	    //southvertypos.SetX(pe*postmp.X()+southvertypos.X());
	    hitpos.SetY(1.0*postmp.Y()+hitpos.Y());
	    if (fVerbose){
	      mf::LogInfo("CRTHitRecoAlg: ") << "!(region==South && layer==1) : \t" << " feb: " << (int)macs.back() << " ,chan : \t" << chan
			<< " ,pe: \t"<< pe << ", adc:\t" << data->fAdc[chan] << ", time: \t"<< data->fTs0
			<< " ,x: \t"<< postmp.X() <<" ,y: \t" << postmp.Y()  <<" ,z: \t" << postmp.Z() << " petotal : "<< petot<< '\n';
	    }
	    ny++;
	    pey+=pe;
	    if(postmp.Y()<ymin)
	      ymin = postmp.Y();
	    if(postmp.Y()>ymax)
	      ymax = postmp.Y();
	    if(region!="South") { //region is E/W/N
	      //    hitpos.SetX(pe*postmp.X()+hitpos.X());
	      hitpos.SetX(1.0*postmp.X()+hitpos.X());
	      nx++;
	      pex+=pe;
	      if(postmp.X()<xmin)
		xmin = postmp.X();
	      if(postmp.X()>xmax)
		xmax = postmp.X();
	    }
	  } 
	  else { //else vertical strips in South wall
	    // hitpos.SetX(pe*postmp.X()+hitpos.X());
	    // hitpos.SetY(pe*postmp.Y()+hitpos.Y());
	    // southvertypos.SetY(pe*postmp.Y()+southvertypos.Y());
	    hitpos.SetX(1.0*postmp.X()+hitpos.X());
	    if (fVerbose){
	    mf::LogInfo("CRTHitRecoAlg: ") << "vertical strips in South wall : \t" << " feb: " << (int)macs.back() << " ,chan : \t" << chan
                      << " ,pe: \t"<< pe << ", adc:\t" << data->fAdc[chan] << ", time: \t"<< data->fTs0
                      << " ,x: \t"<< postmp.X() <<" ,y: \t" << postmp.Y()  <<" ,z: \t" << postmp.Z() << " petotal : "<< petot<< '\n';
	    }
	    nx++;
	    pex+=pe;
	    if(postmp.X()<xmin)
	      xmin = postmp.X();
	    if(postmp.X()>xmax)
	      xmax = postmp.X();
	  }
	  
	  //nz = ny
	  //hitpos.SetZ(pe*postmp.Z()+hitpos.Z());
	  hitpos.SetZ(1.0*postmp.Z()+hitpos.Z());
	  nz++;
	  if (fVerbose){
	    if(region =="South")
	      mf::LogInfo("CRTHitRecoAlg: ") << " South wall z: \t" << " feb: " << (int)macs.back() << " ,chan : \t" << chan
			<< " ,pe: \t"<< pe << ", adc:\t" << data->fAdc[chan] << ", time: \t"<< data->fTs0
			<< " ,x: \t"<< postmp.X() <<" ,y: \t" << postmp.Y()  <<" ,z: \t" << postmp.Z()<< " petotal : "<< petot<< '\n';
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

        }//loop over channels
	


        //correct trigger time for propegation delay
	auto const& adsGeo = adGeo.SensitiveVolume(adsid_max); //trigger stripi
	//auto const& adsGeo = adGeo.SensitiveVolume(last_chan); //trigger strip

	// Timing determination
	// ttrigs[layer].push_back(data->fTs0);// - adsGeo.HalfLength()*fPropDelay);


	//ttrigs.push_back(data->fTs0 - uint64_t(adsGeo.HalfLength()*fPropDelay));
        //t1trigs.push_back(data->fTs1 - uint64_t(adsGeo.HalfLength()*fPropDelay));
	ttrigs.push_back(data->fTs0);
	t1trigs.push_back(data->fTs1); 
	tpos.push_back(postrig);
	ntrig++;
	    
	if(fVerbose)
	  mf::LogInfo("CRTHitRecoAlg: ") <<  "raw T0: " << data->fTs0 << " ,half lenth : \t" << adsGeo.HalfLength() << " ,delay: \t" <<fPropDelay 
		    << " ,corrected time: " << data->fTs0 - uint64_t(adsGeo.HalfLength()*fPropDelay) << '\n';


	if(region=="South" && layer==1) {
	  southt0_h = data->fTs0;
	  if(fVerbose)
	    mf::LogInfo("CRTHitRecoAlg: ") << "southt0_h : " << layer << "\t" << southt0_h << '\n';
	}else if (region=="South" && layer!=1){
	  southt0_v = data->fTs0;
	  if(fVerbose)
	    mf::LogInfo("CRTHitRecoAlg: ") << "southt0_v : " << layer << "\t"<< southt0_v << '\n';
	}else {
	  southt0_h = -999;
	  southt0_v = -999;
	}
	

    }//loop over FEBs
    


    std::sort(layID.begin(), layID.end());
    layID.resize(std::distance(layID.begin(), std::unique(layID.begin(), layID.end())));

    if(fVerbose)
      mf::LogInfo("CRTHitRecoAlg: ") <<  "size of layer ID : \t"  <<layID.size() << '\n';


    
    //no channels above threshold or no intermodule coincidences? return empty hit
    //    if(nabove==0 || layID.size()!=2) {
    if(nabove==0 || layID.size() < 2) {
      if(nabove==0 && fVerbose) mf::LogInfo("CRTHitRecoAlg: ") << "no channels above threshold!" << '\n';
      if(layID.size()<2 && fVerbose) mf::LogInfo("CRTHitRecoAlg: ") << "no coincidence found" << '\n';
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
    geo::Point_t center;

    for(auto const&  infn: informationA) {    
      auto i = &infn - informationA.data();
      auto const& adsGeo = adGeo.SensitiveVolume(infn.strip); //trigger stripi
      center = adsGeo.GetCenter();

      if(fVerbose)
	mf::LogInfo("CRTHitRecoAlg: ")<< "A type ----------> time: " << (long long int)infn.t0 << " ,macs "<< (int)infn.mac5s  //<< '\n';
                 << " ,chal "<< infn.channel
                 << " ,  position "<<infn.pos[2] << '\n';
    
      if ((int)infn.mac5s != (int)informationA[i+1].mac5s and i < (int)informationA.size()-1){
	layer1 = true;

	if ((int)infn.mac5s % 2 == 0) t1_1 = infn.t0;
	else t1_1 = informationA[i+1].t0;
	if ((int)informationA[i+1].mac5s % 2 != 0) t1_2 = informationA[i+1].t0;
	else t1_2 = infn.t0;
	if(fVerbose)
	  mf::LogInfo("CRTHitRecoAlg: ")<< "t1: " << t1_1 << ", t2:"<< t1_2 << ", deltat : "<< int64_t(t1_1 - t1_2) << '\n';

 	float zaxixpos = 0.5*(int64_t(t1_1 - t1_2)/fPropDelay);

	posA = adsGeo.GetCenter() + geo::Zaxis() * zaxixpos;
	if (fVerbose) mf::LogInfo("CRTHitRecoAlg: ")<< "posA (==0): "<< posA<< '\n';


	if(fVerbose)
	  mf::LogInfo("CRTHitRecoAlg: ")<< i << " ,1st mac5: "<< (int)infn.mac5s << " 1st time: " << (long long int)infn.t0 
		   << " ,2nd mac5: "<<(int)informationA[i+1].mac5s << ", 2nd time " << (long long int)informationA[i+1].t0 << " , deltaT: "
		   << int64_t(t1_1 - t1_2) << " , length: " << adsGeo.Length()
		   << " ,propagation delay: "<< fPropDelay << " , pos z: " << zaxixpos
		   << " , center: " << adsGeo.GetCenter() << " , zaxis: "<< geo::Zaxis() <<  " , half length:  " << adsGeo.HalfLength()
		   << " , actual pos w.rt. z: " << posA << '\n';
      }
    }

    for(auto const&  infn: informationB) {

     
      if(fVerbose)
	mf::LogInfo("CRTHitRecoAlg: ")<< " B type ----------> time: " << infn.t0 << " ,macs "<< (int)infn.mac5s  //<< '\n';
                 << " ,chal "<< infn.channel
                 << " ,  position "<<infn.pos[2] << '\n';

      auto i = &infn - informationB.data();
      auto const& adsGeo = adGeo.SensitiveVolume(infn.strip); //trigger stripi


      if ((int)infn.mac5s != (int)informationB[i+1].mac5s and i < (int)informationB.size()-1){

      layer2 = true;
      mac5_2 = (int)infn.mac5s;
      t0_2 = (uint64_t)infn.t0;

      if ((int)infn.mac5s % 2 == 0) t2_1 = infn.t0;
      else t2_1 = informationB[i+1].t0;
      if ((int)informationB[i+1].mac5s % 2 != 0) t2_2 = informationB[i+1].t0;
      else t2_2 = infn.t0;

      if(fVerbose)
	mf::LogInfo("CRTHitRecoAlg: ")<< "t1: " << t2_1 <<", t2:"<< t2_2 << ", deltat : "<< int64_t(t2_1 - t2_2) << '\n';
      //if (foutCSVFile) filecsv << plane << "\t"<<  int64_t(t2_1 - t2_2) << "\n";
      float zaxixpos = 0.5*(int64_t(t2_1 - t2_2)/fPropDelay);

      posB = adsGeo.GetCenter() + geo::Zaxis() * zaxixpos;
      if (fVerbose) mf::LogInfo("CRTHitRecoAlg: ")<< "posB (== 0): "<< posB<< '\n';


      if(fVerbose)
	mf::LogInfo("CRTHitRecoAlg: ")<< i << " ,1st mac5: "<< (int)infn.mac5s << " 1st time: " << (long long int)infn.t0
		 << " ,2nd mac5: "<<(int)informationB[i+1].mac5s << ", 2nd time " << (long long int)informationB[i+1].t0 << " , deltaT: "
		 << int64_t(t2_1 - t2_2) << " , length: " << adsGeo.Length()
		 << " ,propagation delay: "<< fPropDelay << " , pos z: " << zaxixpos
		 << " , center: " << adsGeo.GetCenter() << " , zaxis: "<< geo::Zaxis() <<  " , half length:  " << adsGeo.HalfLength()
		 << " , actual pos w.rt. z: " << posB << '\n';

      }
      
    }

    int crossfeb = std::abs(mac5_1 - mac5_2);


    // side crt and match the both layers
    if (layer1 && layer2 && region!="South" && region!="North" ){//&& nx==4){
      float avg = 0.5*(posA.Z() + posB.Z());
      hitpos.SetZ(avg);
      hitpos.SetX(hitpos.X()*1.0/nx);
      hitpos.SetY(hitpos.Y()*1.0/nx);
      //hitpos.SetX(hitpos.X()*1.0/petot);
      //hitpos.SetY(hitpos.Y()*1.0/petot);
      if(fVerbose)
	mf::LogInfo("CRTHitRecoAlg: ") << "z position in layer 1: "<< posA.Z() << " and in layer 2 "<< posB.Z() 
		  << " average is "<< (posA.Z()+ posB.Z())/2. << " ,hitpos z " << hitpos[2] << '\n';

      
    }else if ((int)informationA.size()==1 and (int)informationB.size()==1
	      and (crossfeb == 7 or crossfeb == 5) and 
	      region!="South" && region!="North"){

      int z_pos =  int64_t(t0_1 - t0_2)/(uint64_t(2*fPropDelay));
      crossfebpos =  center + geo::Zaxis()*z_pos;

      hitpos.SetZ(crossfebpos.Z());
      hitpos.SetX(hitpos.X()*1.0/nx);
      hitpos.SetY(hitpos.Y()*1.0/nx);
      //hitpos.SetX(hitpos.X()*1.0/petot);
      //hitpos.SetY(hitpos.Y()*1.0/petot);
      if(fVerbose)
	mf::LogInfo("CRTHitRecoAlg: ") << "hello hi namaskar,  hitpos z " << hitpos[2] << '\n';
      // side crt and only single layer match
    }else if (layer1 && region!="South" && region!="North"){// && nx==1){
      hitpos.SetZ(posA.Z());
      hitpos.SetX(hitpos.X()*1.0/nx);
      hitpos.SetY(hitpos.Y()*1.0/nx);
      //hitpos.SetX(hitpos.X()*1.0/petot);
      //hitpos.SetY(hitpos.Y()*1.0/petot);
      if(fVerbose)
	mf::LogInfo("CRTHitRecoAlg: ") << " same layer coincidence:  z position in layer 1: "<< posA.Z() << " ,hitpos z " << hitpos[2] << '\n';
      
      // side crt and only single layer match
    }else if (layer2 && region!="South" && region!="North" ){//&& nx==1){
      hitpos.SetZ(posB.Z());
      hitpos.SetX(hitpos.X()*1.0/nx);
      hitpos.SetY(hitpos.Y()*1.0/nx);
      //hitpos.SetX(hitpos.X()*1.0/petot);
      //hitpos.SetY(hitpos.Y()*1.0/petot);
      if(fVerbose)
	mf::LogInfo("CRTHitRecoAlg: ") << " same layer coincidence: z position in layer 2 "<< posB.Z() << " ,hitpos z " << hitpos[2] << '\n';
      
    }else if (region!="South" && region!="North" ){//&& nx==2){
      hitpos*=1.0/nx;
      //hitpos.SetX(hitpos.X()*1.0/petot);
      //hitpos.SetY(hitpos.Y()*1.0/petot);
      //hitpos.SetZ(hitpos.Z()*1.0/petot);
      if (fVerbose) mf::LogInfo("CRTHitRecoAlg: ") << " In side CRTs [E/W] x: \t"<< hitpos[0] <<" ,y: \t" << hitpos[1]  <<" ,z: \t" << hitpos[2]<< '\n';
    }/*else {
      hitpos.SetX(-99999);
      hitpos.SetY(-99999);
      hitpos.SetZ(-99999);
    }*/
    

    //finish averaging and fill hit point array
    if(region=="South") {
      /*     
      hitpos.SetX(hitpos.X()*1.0/pex);
      hitpos.SetZ(hitpos.Z()*1.0/petot);
      */
      ///*
	hitpos.SetX(hitpos.X()/nx);
	hitpos.SetY(hitpos.Y()/ny);
	hitpos.SetZ(hitpos.Z()/nz);
	//  */
	// }else
      //hitpos*=1.0/petot; //hit position weighted by deposited charge
   
    }else if (region=="North"){
      //hitpos*=1.0/petot;
      hitpos*=1.0/nz;

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
    
    if (region=="South" && hitpoint[0] >= 366. && hitpoint[1] > 200. && fVerbose)
      mf::LogInfo("CRTHitRecoAlg: ") << "I am looking for south wall :   macs " << (int)macs.back() << " x: \t"<< hitpoint[0] 
		<<" ,y: \t" << hitpoint[1]  <<" ,z: \t" << hitpoint[2] << '\n';
    
    if (fVerbose){
      if (region=="North") mf::LogInfo("CRTHitRecoAlg: ") << "north wall x: \t"<< hitpoint[0] <<" ,y: \t" << hitpoint[1]  <<" ,z: \t" << hitpoint[2]<< '\n';
    } 
    if (fVerbose) mf::LogInfo("CRTHitRecoAlg: ") << " nx: \t"<< nx <<" ,ny: \t" << ny  <<" ,nz: \t" << nz<< '\n';
    if (fVerbose) mf::LogInfo("CRTHitRecoAlg: ") << " x: \t"<< hitpoint[0] <<" ,y: \t" << hitpoint[1]  <<" ,z: \t" << hitpoint[2]<< '\n';
    
    //time stamp averaged over all FEBs 
	// NOTE: in cases where ttrigs>12, then adding the unix timestamps @ line 1136-1137 can overload the max value of an uint64_t, giving a bad T0 timestamp. 
	// for these cases(ttrigss.size()>10), average over the ns portion of the full unix timestamps and then add thit_s to thit_ns get the averaged full unix timestamp 
    uint64_t thit = 0., t1hit = 0.;// thit_1 = 0.;
    uint64_t thit_ns =0., thit_s=0.;//'thit_s': s portion of Unix Timestamp, 'thit_ns': ns portion of Unix Timestamp.

    ttrigs.resize(std::distance(ttrigs.begin(), std::unique(ttrigs.begin(), ttrigs.end())));
    //thit = std::reduce(ttrigs.begin(), ttrigs.end()) / (uint64_t)ttrigs.size();
    //thit = std::accumulate(ttrigs.begin(), ttrigs.end()) / (uint64_t)ttrigs.size();
   
    t1trigs.resize(std::distance(t1trigs.begin(), std::unique(t1trigs.begin(), t1trigs.end())));
    //    t1hit = std::reduce(t1trigs.begin(), t1trigs.end()) / (uint64_t)t1trigs.size();
    //    t1hit = std::accumulate(t1trigs.begin(), t1trigs.end()) / (uint64_t)t1trigs.size();
    if(ttrigs.size()>10){
      for(uint64_t const t : ttrigs){
        thit_s = t - t%1'000'000'000;
	if (region=="North" || region=="South") thit_ns += t%1'000'000'000-uint64_t(200.*fPropDelay)-fSiPMtoFEBdelay;
	else thit_ns += t%1'000'000'000 -uint64_t(400.*fPropDelay)-fSiPMtoFEBdelay; 
      } 
      if(fVerbose)
	mf::LogInfo("CRTHitRecoAlg: ")  << "Average: thit_ns / ttrigs.size = " << thit_ns << " / " << uint64_t(ttrigs.size()) << " = ";
      thit_ns = thit_ns /uint64_t(ttrigs.size()); //average ns portion of timestamp of CRT data products in the CRT Hit 
      thit = thit_s + thit_ns; //add averaged ns portion of timestamp to full s portion of unix timestamp
       if(fVerbose)
	mf::LogInfo("CRTHitRecoAlg: ") << "add on full second: thit = thit_s + thit_ns = " << thit_s << " + " << thit_ns << " = " << thit << "\n";
	    
    }
    else{
      for(uint64_t const t : ttrigs){
        if (region=="North" || region=="South") /*thit += t;*/ thit += t-uint64_t(200.*fPropDelay)-fSiPMtoFEBdelay;
        else /*thit += t;*/thit += t-uint64_t(400.*fPropDelay)-fSiPMtoFEBdelay;
        if (fVerbose) 
	  mf::LogInfo("CRTHitRecoAlg: ") << " Inside the loop: t: \t"<< t << "\t" << uint64_t(200.*fPropDelay) 
				       << " t0hit : "<< thit <<" size ttrig: \t" <<   ttrigs.size()<< '\n';
     
      }
      thit=thit/uint64_t(ttrigs.size());
    }
    for(double const t1 : t1trigs)
      t1hit +=   t1-uint64_t(400.*fPropDelay)-fSiPMtoFEBdelay;
         

    t1hit=t1hit/uint64_t(t1trigs.size());

    if (region=="South" && fVerbose) mf::LogInfo("CRTHitRecoAlg: ") << "..................... Hello ....Welcome to Beam............"<< '\n';
    if (fVerbose) mf::LogInfo("CRTHitRecoAlg: ") << " <time>: T0: \t"<< thit << " T1 : "<< t1hit <<" size ttrig: \t" <<   ttrigs.size()<< '\n';

    if(region=="South" && fVerbose)
      mf::LogInfo("CRTHitRecoAlg: ") <<  "southt0_h: " << southt0_h << " ,southt0_v : " << southt0_v 
		<< " ,deltaT: \t" << int64_t(southt0_h - southt0_v) << '\n';
    if (foutCSVFile) filecsv << southt0_h << "\t" << southt0_v << "\t"<<  int64_t(southt0_h - southt0_v) << "\n";

    //error estimates (likely need to be revisted)
    auto const& adsGeo = adGeo.SensitiveVolume(adsid_max);
    if(region!="North" && region!="South"){
      hitpointerr[0] = (xmax-xmin)/sqrt(12);
      hitpointerr[1] = (ymax-ymin)/sqrt(12);
      hitpointerr[2] = (zmax-zmin)/sqrt(12);
      //      hitpointerr[2] = adsGeo.Length()/sqrt(12);
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
   
    Long64_t thit1=(Long64_t)(thit-GlobalTrigger[(int)macs.at(0)]);
    
    //generate hit
    CRTHit hit = FillCRTHit(macs,pesmap,petot,thit,thit1,plane,hitpoint[0],hitpointerr[0],
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

//-----------------------------------------------------------------------------
ULong64_t icarus::crt::GetMode(std::vector<std::pair<int, ULong64_t>> vector) {

        sort(vector.begin(), vector.end(), icarus::crt::sortbytime);

        int modecounter = 0;
        int isnewmodecounter = 0;
        ULong64_t Mode = 0;
        ULong64_t isnewMode = 0;
        bool isFirst = true;
        for (auto i : vector) {
                if (!isFirst) {
                        if (i.second == Mode) modecounter++;
                        else if (i.second !=isnewMode) {
                                isnewMode = i.second;
                                isnewmodecounter = 1;
                        }
                        else if (i.second == isnewMode) {
                                isnewmodecounter++;
                                if (isnewmodecounter > modecounter) {
                                        Mode = isnewMode;
                                        modecounter = isnewmodecounter;
                                }
                        }
                }
                else {
                        isFirst = false;
                        Mode = i.second;
                        modecounter++;
                }
        }
        return Mode;
}
