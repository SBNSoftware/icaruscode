#include "CRTHitRecoAlg.h"
#include <algorithm>
#include "larcore/CoreUtils/ServiceUtil.h"  // lar::providerFrom()
using namespace icarus::crt;

CRTHitRecoAlg::CRTHitRecoAlg(const fhicl::ParameterSet& pset)
    : CRTHitRecoAlg() {
  this->reconfigure(pset);
}

CRTHitRecoAlg::CRTHitRecoAlg()
    : fGeometryService(lar::providerFrom<geo::Geometry>()),
      fChannelMap(
          art::ServiceHandle<icarusDB::IICARUSChannelMap const>{}.get()) {}

//---------------------------------------------------------------------
void CRTHitRecoAlg::reconfigure(const fhicl::ParameterSet& pset) {
  fVerbose = pset.get<bool>("Verbose", false);
  fUseReadoutWindow = pset.get<bool>("UseReadoutWindow", false);
  fQPed = pset.get<double>("QPed", 0.);
  fQSlope = pset.get<double>("QSlope", 0.);
  fPropDelay = pset.get<double>("PropDelay", 0.);
  fPEThresh = pset.get<double>("PEThresh", 0.);
  ftopGain = pset.get<double>("topGain", 0.);
  ftopPed = pset.get<double>("topPed", 0.);
  fSiPMtoFEBdelay = pset.get<uint64_t>("SiPMtoFEBdelay", 0.);
  fCoinWindow = pset.get<uint64_t>("CoinWindow", 0.);
  fCrtWindow = pset.get<uint64_t>("CrtWindow", 0.);
  fCSVFile = pset.get<std::string>("CSVFile", "");
  fData = pset.get<bool>("Data", false);
  fGlobalT0Offset = pset.get<double>("GlobalT0Offset",0.);
  if (!fCSVFile.empty()) filecsv.open(fCSVFile);
  //if(foutCSVFile) filecsv << "0, 1 \n";
  return;
}

//---------------------------------------------------------------------------------------
// Preselect CRTData to be constructed into a CRT Hit 
// 1. look for CRTData within +/- 3ms (fCrtWindow) of trigger timestamp
// 2. filtering 
// -- Side CRTs filter out low PEs (values < fPEThresh) and hits that are not a 
//    T0 or T1 reference hit (if its a T0 or T1 reset event, they can have low PE..)
// -- Bottom CRT Hit reco also using fPEThresh, probably need to be revisited
// -- Top CRT does not require additional filtering here
// 
vector<art::Ptr<CRTData>> CRTHitRecoAlg::PreselectCRTData(
    const vector<art::Ptr<CRTData>>& crtList, uint64_t trigger_timestamp) {
  if (fVerbose)
    mf::LogInfo("CRTHitRecoAlg: ") << "In total " << crtList.size()
                                   << " CRTData found in an event" << '\n';
  vector<art::Ptr<CRTData>> crtdata;
  bool presel = false;

  for (size_t febdat_i = 0; febdat_i < crtList.size(); febdat_i++) {
    uint8_t mac = crtList[febdat_i]->fMac5;
    int adid = fCrtutils.MacToAuxDetID(mac, 0);
    char type = fCrtutils.GetAuxDetType(adid);

    /// Looking for data within +/- 3ms within trigger time stamp
    /// Here t0 - trigger time -ve
    if (fData && (std::fabs(int64_t(crtList[febdat_i]->fTs0 -
                                    trigger_timestamp)) > fCrtWindow))
      continue;

    if (type == 'm') { // 'm' = MINOS, Side CRTs
      for (int chan = 0; chan < 32; chan++) {
        std::pair<double, double> const chg_cal =
            fChannelMap->getSideCRTCalibrationMap((int)crtList[febdat_i]->fMac5,
                                                  chan);
        float pe =
            (crtList[febdat_i]->fAdc[chan] - chg_cal.second) / chg_cal.first;
        if (pe <= fPEThresh && !crtList[febdat_i]->IsReference_TS1() &&
            !crtList[febdat_i]->IsReference_TS0())
          continue;  // filter out low PE and non-T1 and non-T0 ref values
        presel = true;
        /*if (fVerbose)
          mf::LogInfo("CRTHitRecoAlg: ")
              << "\nfebP (mac5, channel, gain, pedestal, adc, pe) = ("
              << (int)crtList[febdat_i]->fMac5 << ", " << chan << ", "
              << chg_cal.first << ", " << chg_cal.second << ","
              << crtList[febdat_i]->fAdc[chan] << "," << pe << ")\n";*/
      }
    } else if (type == 'c') { // 'c' = CERN, Top CRTs
      for (int chan = 0; chan < 32; chan++) {
        // float pe = (crtList[febdat_i]->fAdc[chan]-fQPed)/fQSlope;
        // if(pe<=fPEThresh) continue;
        presel = true;
      }
    } else if (type == 'd') { //'d' = Double Chooz, Bottom CRTs
      for (int chan = 0; chan < 64; chan++) {
        float pe = (crtList[febdat_i]->fAdc[chan] - fQPed) / fQSlope;
        if (pe <= fPEThresh) continue;
        presel = true;
      }
    }

    if (presel) crtdata.push_back(crtList[febdat_i]);
    presel = false;
  }
  mf::LogInfo("CRTHitRecoAlg:")
      << "Found " << crtdata.size() << " after preselection " << '\n';
  return crtdata;
}

//---------------------------------------------------------------------------------------
vector<pair<sbn::crt::CRTHit, vector<int>>> CRTHitRecoAlg::CreateCRTHits(
    vector<art::Ptr<CRTData>> crtList, uint64_t trigger_timestamp) {
  vector<pair<CRTHit, vector<int>>> returnHits;
  vector<int> dataIds;

  uint16_t nMissC = 0, nMissD = 0, nMissM = 0, nHitC = 0, nHitD = 0, nHitM = 0;
  if (fVerbose)
    mf::LogInfo("CRTHitRecoAlg: ")
        << "Found " << crtList.size() << " FEB events" << '\n';

  map<string, int> regCounts;
  std::set<string> regs;
  map<string, vector<size_t>> sideRegionToIndices;

  // sort by the time
  std::sort(crtList.begin(), crtList.end(), compareBytime);
  // Below is the calculation of the Global Trigger timestamp from the T1 reset
  // currently only done for top CRT T1 reset values
  // TODO: Validate side CRT global trigger timestamp reconstruction to be used 
  // for the reference for side CRT hits, currently a couple ns off. -AH 01/19/2024

  // Load Delays map for Top CRT
  CRT_delay_map const FEB_delay_map = LoadFEBMap();
  std::vector<std::pair<int, ULong64_t>> CRTReset;
  ULong64_t TriggerArray[305] = {0};
  
  for (size_t crtdat_i = 0; crtdat_i < crtList.size(); crtdat_i++) {
    uint8_t mac = crtList[crtdat_i]->fMac5;
    int adid = fCrtutils.MacToAuxDetID(mac, 0);
    char type = fCrtutils.GetAuxDetType(adid);
    string region = fCrtutils.GetAuxDetRegion(adid);
    
    // For the time being, Only Top CRT delays are loaded, nothing to do for
    // Side CRT yet
    if (type == 'c' && crtList[crtdat_i]->IsReference_TS1()) {
      ULong64_t Ts0T1ResetEvent = crtList[crtdat_i]->fTs0 +
                                  FEB_delay_map.at((int)mac + 73).T0_delay -
                                  FEB_delay_map.at((int)mac + 73).T1_delay;
      TriggerArray[(int)mac] = Ts0T1ResetEvent;
      CRTReset.emplace_back((int)mac, Ts0T1ResetEvent);  // single GT
    }
  }

  const int trigger_offset =
      60;  // Average distance between Global Trigger and Trigger_timestamp (ns)
           // TODO: Make configurable parameter
  ULong64_t GlobalTrigger = trigger_timestamp;
  if (!CRTReset.empty()) GlobalTrigger = GetMode(CRTReset);
  // Add average difference between trigger_timestamp and Global trigger
  else
    GlobalTrigger =
        GlobalTrigger - trigger_offset;  // In this event, the T1 Reset was
                                         // probably "vetoed" by the T0 Reset
  for (int i = 0; i < 305; i++) {
    if (TriggerArray[i] == 0) TriggerArray[i] = GlobalTrigger;
  }
  
  // loop over time-ordered CRTData
  for (size_t febdat_i = 0; febdat_i < crtList.size(); febdat_i++) {
    uint8_t mac = crtList[febdat_i]->fMac5;
    int adid = fCrtutils.MacToAuxDetID(mac, 0);  // module ID

    string region = fCrtutils.GetAuxDetRegion(adid);
    char type = fCrtutils.GetAuxDetType(adid);
    CRTHit hit;

    dataIds.clear();

    // CERN modules (intramodule coincidence)
    if (type == 'c') {
      hit = MakeTopHit(crtList[febdat_i], TriggerArray);  // single GT
      // hit = MakeTopHit(crtList[febdat_i], TriggerArray_top); //3 seperate GT
      if (IsEmptyHit(hit))
        nMissC++;
      else {
        dataIds.push_back(febdat_i);
        returnHits.push_back(std::make_pair(hit, dataIds));
        if ((regs.insert(region)).second)
          regCounts[region] = 1;
        else
          regCounts[region]++;

        nHitC++;
      }
    }

    // DC modules (intramodule coincidence)
    if (type == 'd') {
      hit = MakeBottomHit(crtList[febdat_i]);
      if (IsEmptyHit(hit))
        nMissD++;
      else {
        dataIds.push_back(febdat_i);
        returnHits.push_back(std::make_pair(hit, dataIds));
        if ((regs.insert(region)).second)
          regCounts[region] = 1;
        else
          regCounts[region]++;

        nHitD++;
      }
    }

    if (type == 'm') sideRegionToIndices[region].push_back(febdat_i);

  }  // End loop over time-ordered CRTData products

  vector<size_t> unusedDataIndex;
  
  for (auto const& regIndices : sideRegionToIndices) {
    if (fVerbose)
      mf::LogInfo("CRTHitRecoAlg: ")
          << "searching for side CRT hits in region, " << regIndices.first
          << '\n';

    vector<size_t> indices = regIndices.second;

    if (fVerbose)
      mf::LogInfo("CRTHitRecoAlg: ")  
	<< "\n-------------------------\nCreateCRTHits: found " 
	<< indices.size() << " side CRT hits in region " << regIndices.first 
	<< "\n----------\n";

    for (size_t index_i = 0; index_i < indices.size(); index_i++) {
      dataIds.clear();
      dataIds.push_back(indices[index_i]);

      vector<art::Ptr<CRTData>> coinData = {crtList[indices[index_i]]};

      if (fVerbose)
        mf::LogInfo("CRTHitRecoAlg: ")
            << "  " << coinData.size()
            << " data products entering to time ordering" << '\n';

      // inner loop over data after data_i in time
      for (size_t index_j = index_i + 1; index_j < indices.size(); index_j++) {
        /*if (fVerbose)
          mf::LogInfo("CRTHitRecoAlg: ")
              << "index_i = " << index_i << ", index_j = " << index_j << ", ts0_i = "
              << crtList[indices[index_i]]->fTs0 << ", ts0_j = "
	      << crtList[indices[index_j]]->fTs0 << "\nts0_i+fCoinWindow = "  
	      << crtList[indices[index_i]]->fTs0 << " + " << fCoinWindow 
	      << " = "
              << crtList[indices[index_i]]->fTs0 + fCoinWindow << '\n';*/

        if (crtList[indices[index_j]]->fTs0 < crtList[indices[index_i]]->fTs0)
          mf::LogError("CRTHitRecoAlg::CreateCRTHits")
              << "bad time ordering!" << '\n';

        /*if (fVerbose)
          mf::LogInfo("CRTHitRecoAlg: ")
              << "size ..  " << coinData.size()
              << " data products before coincidence" << '\n';*/
        
        if ((crtList[indices[index_j]]->fTs0 >=
                 crtList[indices[index_i]]->fTs0 &&
             (crtList[indices[index_j]]->fTs0 -
              crtList[indices[index_i]]->fTs0) < fCoinWindow) ||
            (crtList[indices[index_j]]->fTs0 <
                 crtList[indices[index_i]]->fTs0 &&
             (crtList[indices[index_i]]->fTs0 -
              crtList[indices[index_j]]->fTs0) < fCoinWindow)) {
          if (fVerbose)
            mf::LogInfo("CRTHitRecoAlg: ")
                << " in coincidence: i: " << index_i << " ,j: " << index_j
                << ", i mac = " << (int)crtList[indices[index_i]]->fMac5
		<< " @ ts0_i = " << crtList[indices[index_i]]->fTs0 
                << ", j mac = " << (int)crtList[indices[index_j]]->fMac5
		<< " @ ts0_j = " << crtList[indices[index_j]]->fTs0 
                << '\n';

          coinData.push_back(crtList[indices[index_j]]);
          dataIds.push_back(indices[index_j]);
        }

        // out of coinWindow

        if ((crtList[indices[index_j]]->fTs0 -
             crtList[indices[index_i]]->fTs0) > fCoinWindow ||
            index_j == indices.size() - 1) {
          /*if (fVerbose)
            mf::LogInfo("CRTHitRecoAlg: ")
                << "out of coincidence: j: " << index_j << ", indices.size = " << indices.size()
                << ", indices.size - 1 = " << indices.size() - 1 << " data products..." << '\n';*/
	  if (coinData.size() < 2) continue; //if coinData<2, no coincidence found
          if (fVerbose)
            mf::LogInfo("CRTHitRecoAlg: ")
                << "attempting to produce MINOS hit from " << coinData.size()
                << " data products..." << '\n';
          uint8_t imac = (int)crtList[indices[index_i]]->fMac5;
          int adid = fCrtutils.MacToAuxDetID(imac, 0);
          string region = fCrtutils.GetAuxDetRegion(adid);
          CRTHit hit = MakeSideHit(coinData, TriggerArray);  // using top CRT GT

          if (IsEmptyHit(hit)) {
	    if(fVerbose) mf::LogInfo("CRTHitRecoAlg") << "no hit produced!!\n-------------------------------\n";
            unusedDataIndex.push_back(indices[index_i]);
            nMissM++;
          } else {
            if (fVerbose)
              mf::LogInfo("CRTHitRecoAlg: ") << "MINOS hit produced" << "\n-------------------------------\n";

            returnHits.push_back(std::make_pair(hit, dataIds));

            if ((regs.insert(regIndices.first)).second)
              regCounts[regIndices.first] = 1;
            else
              regCounts[regIndices.first]++;

            nHitM++;
          }
          index_i = index_j - 1;
          if (index_j == indices.size() - 1) index_i++;

          break;
        }  // if jth data out of coinc window
      }    // inner loop over data
    }      // outer loop over data
  }        // loop over side CRTData products
  
  if (fVerbose) {
    mf::LogInfo("CRTHitRecoAlg") << returnHits.size() << " CRT hits produced!" << '\n'
                       << "  nHitC: " << nHitC << " , nHitD: " << nHitD
                       << " , nHitM: " << nHitM << '\n'
                       << "  nMisC: " << nMissC << " , nMisD: " << nMissD
                       << " , nMisM: " << nMissM << '\n';
    auto cts = regCounts.begin();
    mf::LogInfo("CRT") << " CRT Hits by region" << '\n';
    while (cts != regCounts.end()) {
      if(fVerbose) mf::LogInfo("CRTHitRecoAlg") << "reg: " << (*cts).first
		<< " , hits: " << (*cts).second << '\n';
      cts++;
    }
  }  // if Verbose

  return returnHits;
}
//--------------------------------------------------------------------------------------------
// Function to make filling a CRTHit a bit faster
sbn::crt::CRTHit CRTHitRecoAlg::FillCRTHit(
    vector<uint8_t> tfeb_id, map<uint8_t, vector<pair<int, float>>> tpesmap,
    float peshit, uint64_t time0, Long64_t time1, int plane, double x,
    double ex, double y, double ey, double z, double ez, string tagger) {
  CRTHit crtHit;
  crtHit.feb_id = tfeb_id;
  crtHit.pesmap = tpesmap;
  crtHit.peshit = peshit;
  crtHit.ts0_s_corr = time0 / 1'000'000'000;
  crtHit.ts0_ns = time0 % 1'000'000'000;
  crtHit.ts0_ns_corr = time0;
  crtHit.ts1_ns =
      time1 /*% 1'000'000'000*/;  // TODO: Update the CRTHit data product
                                  // /sbnobj/common/CRT . Discussion with SBND
                                  // people needed
  crtHit.ts0_s = time0 / 1'000'000'000;//'
  crtHit.plane = plane;
  crtHit.x_pos = x;
  crtHit.x_err = ex;
  crtHit.y_pos = y;
  crtHit.y_err = ey;
  crtHit.z_pos = z;
  crtHit.z_err = ez;
  crtHit.tagger = tagger;

  return crtHit;

}  // CRTHitRecoAlg::FillCRTHit()

//------------------------------------------------------------------------------------------
int64_t CRTHitRecoAlg::RegionDelay(std::string const& region) const {
  return fSiPMtoFEBdelay +
         uint64_t(((region == "North" || region == "South") ? 200. : 400) *
                  fPropDelay);
}
//------------------------------------------------------------------------------------------

sbn::crt::CRTHit CRTHitRecoAlg::MakeTopHit(
    art::Ptr<CRTData> data,
    ULong64_t GlobalTrigger[305]) {  // single GT: GlobalTrigger[305], 3
                                     // seperate GT: GlobalTrigger[232]
  uint8_t mac = data->fMac5;
  if (fCrtutils.MacToType(mac) != 'c')
    mf::LogError("CRTHitRecoAlg::MakeTopHit")
        << "CRTUtils returned wrong type!" << '\n';

  map<uint8_t, vector<pair<int, float>>> pesmap;
  int adid = fCrtutils.MacToAuxDetID(mac, 0);          // module ID
  auto const& adGeo = fGeometryService->AuxDet(adid);  // module
  string region = fCrtutils.GetAuxDetRegion(adid);
  int plane = fCrtutils.AuxDetRegionNameToNum(region);
  double hitpointerr[3];
  TVector3 hitpos(0., 0., 0.);
  float petot = 0., pemax = 0., pemaxx = 0., pemaxz = 0.;
  int adsid_max = -1, nabove = 0;
  TVector3 postrig;
  bool findx = false, findz = false;
  int maxx = 0, maxz = 0;
  double sum = 0;
  for (int chan = 0; chan < 32; chan++) {
    sum = sum + data->fAdc[chan];
    std::pair<double, double> const chg_cal =
              fChannelMap->getSideCRTCalibrationMap((int)data->fMac5, chan);
    float pe = (data->fAdc[chan] - chg_cal.second) / chg_cal.first;
    //float pe = (data->fAdc[chan] - ftopPed) / ftopGain;
    //      if(pe<=fPEThresh) continue;
    if (pe<0) pe=0;
    nabove++;
    int adsid = fCrtutils.ChannelToAuxDetSensitiveID(mac, chan);
    petot += pe;
    pesmap[mac].push_back(std::make_pair(chan, pe));

    // TVector3 postmp = fCrtutils.ChanToLocalCoords(mac,chan);
    // strip along z-direction
    if (adsid >= 0 && adsid < 8) {
      // hitpos.SetX(pe*postmp.X()+hitpos.X());
      // hitpos.SetX(postmp.X()+hitpos.X());
      if (pe > pemaxx) {
        pemaxx = pe;
        maxx = adsid;
      }
      findx = true;
    }
    // strip along x-direction
    else if (adsid >= 8 && adsid < 16) {
      // hitpos.SetZ(pe*postmp.Z()+hitpos.Z());
      // hitpos.SetZ(postmp.Z()+hitpos.Z());
      if (pe > pemaxz) {
        pemaxz = pe;
        maxz = adsid;
      }
      findz = true;
    } else {
      mf::LogError("CRTHitRecoAlg::MakeTopHit")
          << "auxDetSensitive ID out of range!" << '\n';
    }
    // identify trigger channel
    if (pe > pemax) {
      TVector3 postmp = fCrtutils.ChanToLocalCoords(mac, chan);
      adsid_max = chan;
      pemax = pe;
      postrig = postmp;
    }
  }

  TVector3 postmp = fCrtutils.ChanToLocalCoords(mac, maxx * 2);
  hitpos.SetX(postmp.X());
  postmp = fCrtutils.ChanToLocalCoords(mac, maxz * 2);
  hitpos.SetZ(postmp.Z());
  int sector = -1;
  if (findz == true && findx == true) sector = (maxz - 8) * 8 + maxx;
  if (!findx)
    mf::LogWarning("CRTHitRecoAlg::MakeTopHit")
        << " no interlayer coincidence found! Missing X coord." << '\n';
  if (!findz)
    mf::LogWarning("CRTHitRecoAlg::MakeTopHit")
        << " no interlayer coincidence found! Missing Z coord." << '\n';

  // no channels above threshold? return empty hit
  if (nabove == 0 || !findx || !findz)
    return FillCRTHit({}, {}, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, "");

  // hitpos*=1.0/petot; //hit position weighted by deposited charge
  // hitpos*=1.0/nabove;
  geo::AuxDetGeo::LocalPoint_t const hitlocal{hitpos.X(), 0., hitpos.Z()};

  auto const& adsGeo = adGeo.SensitiveVolume(adsid_max);  // trigger strip
  uint64_t thit = data->fTs0;
  Long64_t thit1 = data->fTs1;
  thit -= fSiPMtoFEBdelay;
  thit1 -= fSiPMtoFEBdelay;

  // 92.0 is the middle of one of the Top CRT modules (each of them is 184 cm)
  // TO DO: Move hardcoded numbers to parameter fcl files.
  // Another possibility is using object values from GDML, but I (Francesco
  // Poppi) found weird numbers some months ago and needed to double check.
  // double const corrPos = std::max(-hitpos.X(), hitpos.Z());
  // uint64_t const corr = (uint64_t)round(abs((92.0+corrPos)*fPropDelay));
  // //Obsolete
  double corr = 0;
  if (findz == true && findx == true) corr = TopCRT_TimingCorr[sector];
  thit -= (uint64_t)round(corr);
  thit1 -= (int64_t)round(corr);

  auto const hitpoint =
      adGeo.toWorldCoords(hitlocal);  // tranform from module to world coords

  hitpointerr[0] = adsGeo.HalfWidth1() * 2 / sqrt(12);
  hitpointerr[1] = adGeo.HalfHeight();
  hitpointerr[2] = adsGeo.HalfWidth1() * 2 / sqrt(12);
  // thit1 = (Long64_t)(thit-GlobalTrigger[(int)mac+73]);
  if (fData) thit1 = (Long64_t)(thit - GlobalTrigger[(int)mac]);
  else thit1 = thit - fGlobalT0Offset;

  // Remove T1 Reset event not correctly flagged, remove T1 reset events, remove
  // T0 reset events
  if ((sum < 10000 && thit1 < 2'001'000 && thit1 > 2'000'000) ||
      data->IsReference_TS1() || data->IsReference_TS0())
    return FillCRTHit({}, {}, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, "");

  CRTHit hit = FillCRTHit({mac}, pesmap, petot, thit, thit1, plane,
                          hitpoint.X(), hitpointerr[0], hitpoint.Y(),
                          hitpointerr[1], hitpoint.Z(), hitpointerr[2], region);

  return hit;

}  // CRTHitRecoAlg::MakeTopHit

//------------------------------------------------------------------------------------------
sbn::crt::CRTHit CRTHitRecoAlg::MakeBottomHit(art::Ptr<CRTData> data) {
  uint8_t mac = data->fMac5;
  map<uint8_t, vector<pair<int, float>>> pesmap;
  int adid = fCrtutils.MacToAuxDetID(mac, 0);          // module ID
  auto const& adGeo = fGeometryService->AuxDet(adid);  // module
  string region = fCrtutils.GetAuxDetRegion(adid);
  int plane = fCrtutils.AuxDetRegionNameToNum(region);
  double hitpointerr[3];
  TVector3 hitpos(0., 0., 0.);
  float petot = 0., pemax = 0.;
  int adsid_max = -1, nabove = 0;
  TVector3 postrig;
  double xmin = 0., xmax = 0.;

  for (int chan = 0; chan < 64; chan++) {
    float pe = (data->fAdc[chan] - fQPed) / fQSlope;
    if (pe <= fPEThresh) continue;
    nabove++;
    int adsid = fCrtutils.ChannelToAuxDetSensitiveID(mac, chan);
    petot += pe;
    pesmap[mac].push_back(std::make_pair(chan, pe));

    TVector3 postmp = fCrtutils.ChanToLocalCoords(mac, chan);
    // all strips along z-direction
    hitpos.SetX(pe * postmp.X() + hitpos.X());
    if (postmp.X() < xmin) xmin = postmp.X();
    if (postmp.X() > xmax) xmax = postmp.X();

    // identify trigger channel
    if (pe > pemax) {
      adsid_max = adsid;
      pemax = pe;
      postrig = postmp;
    }
  }

  // no channels above threshold? return empty hit
  if (nabove == 0) return FillCRTHit({}, {}, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, "");

  hitpos *= 1.0 / petot;  // hit position weighted by deposited charge
  geo::AuxDetGeo::LocalPoint_t const hitlocal{hitpos.X(), 0., 0.};

  auto const& adsGeo = adGeo.SensitiveVolume(adsid_max);  // trigger strip
  uint64_t thit = data->fTs0 - adsGeo.HalfLength() * fPropDelay;
  Long64_t thit1;
  if (fData) thit1 = data->fTs0 - adsGeo.HalfLength() * fPropDelay; // this should probably be revisited, but we dont currently reconstruct bottom CRT entries anyway
  else thit1 = data->fTs0 - fGlobalT0Offset;
  auto const hitpoint =
      adGeo.toWorldCoords(hitlocal);  // tranform from module to world coords

  hitpointerr[0] = (xmax - xmin + 2 * adsGeo.HalfWidth1() * 2) / sqrt(12);
  hitpointerr[1] = adGeo.HalfHeight();
  hitpointerr[2] = adsGeo.Length() / sqrt(12);

  CRTHit hit = FillCRTHit({mac}, pesmap, petot, thit, thit1, plane, hitpoint.X(),
                          hitpointerr[0], hitpoint.Y(), hitpointerr[1],
                          hitpoint.Z(), hitpointerr[2], region);

  return hit;

}  // CRTHitRecoAlg::MakeBottomHit

//-----------------------------------------------------------------------------------
sbn::crt::CRTHit CRTHitRecoAlg::MakeSideHit(
  vector<art::Ptr<CRTData>> coinData,
  ULong64_t GlobalTrigger[305]) {
  
  std::vector<uint8_t> layerA, layerB; //vectors to sort macs on inner or outer layer
  int adid = fCrtutils.MacToAuxDetID(coinData[0]->fMac5, 0);  // 1st module ID (1st modID in coindata collection)
  int layer = fCrtutils.GetMINOSLayerID(adid); // 1st layerID //modID_1-->adid
  string region = fCrtutils.GetAuxDetRegion(adid);            //region name
  int plane = fCrtutils.AuxDetRegionNameToNum(region);        //region code (ranges from 30-50)
  double hitpoint[3], hitpointerr[3];
  TVector3 hitpos(0., 0., 0.);

  //-----------------------------------------
  // Sort coinData into inner or outer layers
  //-----------------------------------------
  for (auto const& data : coinData) {
    int modID_x = (int)fCrtutils.MacToAuxDetID((int)data->fMac5, 0);
    int layer_x = fCrtutils.GetMINOSLayerID(modID_x); //layerID of data entry in coinData
    if (layer_x == layer){
      std::cout << "layerA: mac " << (int)data->fMac5 << ", modID " << modID_x << ", layer " << layer_x << "\n";
      layerA.push_back(data->fMac5); 
    }
    else {
      std::cout << "layerB: mac " << (int)data->fMac5 << ", modID " << modID_x << ", layer " << layer_x << "\n";
      layerB.push_back(data->fMac5); 
    }
  }

  std::unordered_map<int,std::pair<uint8_t,uint8_t>> layer1_map, layer2_map;
  groupByLayer(layerA, layer1_map, "layer1");
  groupByLayer(layerB, layer2_map, "layer2");
  std::cout << "layerA.size = " << layerA.size() << ", layerB.size " << layerB.size() <<  "\n";
  std::cout << "layer1_map.size = " << layer1_map.size() << ", layer2_map.size " << layer2_map.size() << "\n";

  vector<uint8_t> macs;
  vector<uint64_t> ttrigs, t1trigs; 
  vector<info> informationA, informationB;
  auto const& adGeo = fGeometryService->AuxDet(adid);         // module
  float petot = 0., pemax = 0.;
  map<uint8_t, vector<pair<int, float>>> pesmap;
  int nx = 0, ny = 0, nz = 0, adsid_max = -1;
  TVector3 postrig;

  // 2nd loop over coinData
  for (auto const& data : coinData) {
    if((int)layerA.size()==0 || (int)layerB.size() == 0){
      std::cout << "one layer has no hits, return empty hit!..\n";
      //continue;
      return FillCRTHit({}, {}, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, "");
    }
    // here do i want to check if mac is in the layer1_map?????
    int modID = fCrtutils.MacToAuxDetID(data->fMac5, 0);
    auto in_layer1 = layer1_map.find(modID);
    auto in_layer2 = layer2_map.find(modID);
    /*std::cout << "is modID " << modID << " in layer1_map?\n";
    if (in_layer1 != layer1_map.end()){
      std::cout << "found modID " << modID << " in layer1_map!\n";
    }
    else{
      std::cout << "no modID " << modID << " in layer1_map...\n";
    }
    
    std::cout << "is modID " << modID << " in layer2_map?\n";
    if (in_layer2 != layer2_map.end()){
      std::cout << "found modID " << modID << " in layer2_map!\n";
    }
    else{
      std::cout << "no modID " << modID << " in layer2_map...\n";
      }*/
	//std::cout << "in_layer1, in_layer2 = " << (int)in_layer1 << ", " << (int)in_layer2 << "\n";
    if ((layer1_map.size() > 0 and in_layer1 != layer1_map.end())
	or  (layer2_map.size() > 0 and in_layer2 != layer2_map.end())
	or (layer1_map.size() == 0 and layerA.size() > 0)
	or (layer2_map.size() == 0 and layerB.size() > 0)){
      //std::cout << "modID found or single ended readout!\n";
    }
    else {
      std::cout << "modID " << modID << " not in either layer1 or layer2 map! contining..\n";
      continue;
    }

    macs.push_back(data->fMac5);
    ttrigs.push_back(data->fTs0);
    t1trigs.push_back(data->fTs1);

    size_t maxSize = std::max(layerA.size(), layerB.size());
    for (size_t i = 0; i < maxSize; i++) {
      if (i < layerA.size() && macs.back() == (int)layerA[i]) {
	std::cout << "layerA loop over chan:\n";
	for (int chan = 0; chan < 32; chan++) {
	  processChannel(chan, macs.back(), data, informationA, 
			 petot, pesmap, 
			 hitpos, nx, ny, nz, 
			 //xmin, xmax, ymin, ymax, zmin, zmax,
			 adsid_max, pemax, postrig);
	}
      }
      else if(i < layerB.size() && macs.back() == (int)layerB[i]){
	std::cout << "layerB loop over chan:\n";
	for (int chan = 0; chan < 32; chan++) {
	  processChannel(chan, macs.back(), data, informationB, 
			 petot, pesmap, 
			 hitpos, nx, ny, nz, 
			 //xmin, xmax, ymin, ymax, zmin, zmax,
			 adsid_max, pemax, postrig);
	}
      }
    }
  } // end loop over coindata
  std::cout << "pesmap.size = " << pesmap.size() << "\n";
  if (informationA.size() == 0 and informationB.size() == 0) {
    std::cout << "infoA and infoB are both 0! return empty hit..\n";
    return FillCRTHit({}, {}, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, "");
  }
  std::cout << "adsid_max = " << adsid_max << ", pemax = " << pemax << ", postrig.x = " << postrig.X() << "\n";
    
  std::cout << "macs.size = " << macs.size() << /*", nabove = " << nabove <<*/ "\n";
  std::cout << "infoA.size = " << informationA.size() << ", infoB.size = " << informationB.size() << "\n";
  bool layer1 = false, layer2 = false;
  geo::Point_t posA, posB;
  geo::Point_t center;
  float zposA = 0, zposB = 0;
  float zpos_local = 0, zpos_local_pe=0;
  vector<float> zposA_vec = {}, zposB_vec = {};
  int mac5_1 = -555, mac5_2 = -555; //used to check for hits that cross module
  uint64_t t0_1 = -5, t0_2 = -5; //T0 used to check for hits that cross module
  bool is_infoA_deltaT_over50ns = false, is_infoB_deltaT_over50ns = false, is_crossfeb_deltaT_over50ns = false;
  for (int i = 0; i < (int)informationA.size(); ++i) {
    const auto& infn = informationA[i];
    const auto& adsGeo = adGeo.SensitiveVolume(infn.strip);
    center = adsGeo.GetCenter();
    //// if double ended readout
    if (layer1_map.size()>0){
      if (((int)infn.mac5s != (int)informationA[i+1].mac5s
	   and i < (int)informationA.size()-1
	   and layerA.size()>0)){
	 //and layer1_map.size()>0)
	std::cout << "double ended readout on layer1\n";
	int modID_x = fCrtutils.MacToAuxDetID(infn.mac5s, 0);
	int modID_a = fCrtutils.MacToAuxDetID(informationA[i+1].mac5s, 0);
	std::cout << "attempting to match mac " << (int)infn.mac5s << " on modID " << modID_x << " to mac " << (int)informationA[i+1].mac5s << " on modID " << modID_a << "\n";
	if(modID_x == modID_a){
	  recoZwithTiming(infn.mac5s, informationA[i+1].mac5s, infn.t0, informationA[i + 1].t0, fPropDelay, zpos_local, is_infoA_deltaT_over50ns);
	  layer1 = true;
	  std::cout << "zpos = zpos_local + center = " << zpos_local << " + " << center.Z() << " = " << zpos_local+center.Z() << "\n";
	  posA = geo::Zaxis()*(zpos_local+center.Z());
	  zposA_vec.push_back(zpos_local+center.Z());//normal z reco from timing 
	  recoZwithPE(infn.mac5s, informationA[i+1].mac5s, pesmap, zpos_local_pe);
	  std::cout << "zpos_PE = zpos_local_pe + center = " << zpos_local_pe << " + " << center.Z() << " = " << zpos_local_pe+center.Z() << "\n";
	  //posA = geo::Zaxis()*(zpos_local_pe+center.Z());
	  //zposA_vec.push_back(zpos_local_pe+center.Z());
	  //if (is_infoA_deltaT_over50ns) std::cout << "infoA delta_t > 50 ns!!!\n";
	  //could make a function to do recoZ with PE
	  //recoZwithPE(infn.mac5s, informationA[i+1].mac5s,pesmap, zpos_local_pe)
	  /*if(mac % 2){// if mac odd-->North side
	    for (auto &p : pesmap[(int)mac]){
	      pe_n = p.second;
	    }
	  }
	  else { // if mac even-->South side
	    for (auto &p : pesmap[(int)mac]){
	      pe_s = p.second;
	      }*/
	    
	    //pesmap[
	    //	   recoZwithPE(infn.mac5s, informationA[i+1].mac5s,)
	}
	else{
	  std::cout << "infoA: modIDs dont match.. look for rest in list\n";// index = "<< index <<"\n";
	  int matched_index = checkNextModID(layer1_map, i, informationA);
	  std::cout << "matched_index = "<< matched_index <<"\n";
	  if (matched_index == 0){
	    std::cout << "matched_index = 0 in infoA! continuing...\n";
	    continue;
	  }
	  if (matched_index!=-1) {
	    recoZwithTiming(infn.mac5s, informationA[matched_index].mac5s, infn.t0, informationA[matched_index].t0, fPropDelay, zpos_local, is_infoA_deltaT_over50ns);
	    layer1 = true;
	    std::cout << "zpos = zpos_local + center = " << zpos_local << " + " << center.Z() << " = " << zpos_local+center.Z() << "\n";
	    posA = geo::Zaxis()*(zpos_local+center.Z());
	    zposA_vec.push_back(zpos_local+center.Z());
	    recoZwithPE(infn.mac5s, informationA[matched_index].mac5s, pesmap, zpos_local_pe);
	    std::cout << "zpos_PE = zpos_local_pe + center = " << zpos_local_pe << " + " << center.Z() << " = " << zpos_local_pe+center.Z() << "\n";
	    //posA = geo::Zaxis()*(zpos_local_pe+center.Z());
	    //zposA_vec.push_back(zpos_local_pe+center.Z());
	  }
	}
      }
    }
    //// else if single ended readout on East or West
    else if(region != "South" && region != "North"){
      std::cout << "single ended readout on layer1! infoA.size = " << informationA.size() << "\n";
      TVector3 postmp = fCrtutils.ChanToWorldCoords(infn.mac5s, infn.channel);
      mac5_1 = infn.mac5s;
      t0_1 = infn.t0;
    }
    
  }
  // if mutliple FEB coincidences on single layer, take average of all reco'd z-positions                                                      
  if (zposA_vec.size() > 1){
    zposA = 0;
    std::cout << "ZposA = ";
    
    for(float i_zpos : zposA_vec){
      zposA += i_zpos;
      std::cout << i_zpos << " + ";
    }
    if (fVerbose)
      mf::LogInfo("CRTHitRecoAlg") << " = " << zposA << " / " << zposA_vec.size() << " = " << zposA /zposA_vec.size() << "\n";
    zposA = zposA / zposA_vec.size();
    posA = geo::Zaxis()*zposA;
  }

  // for infoB
  for (int i = 0; i < (int)informationB.size(); ++i) {
    const auto& infn = informationB[i];
    const auto& adsGeo = adGeo.SensitiveVolume(infn.strip);
    center = adsGeo.GetCenter();
    //// if double ended readout
    if (layer2_map.size()>0){
      if (((int)infn.mac5s != (int)informationB[i+1].mac5s
	   and i < (int)informationB.size()-1
	   and layerB.size()>0)){
	 //and layer1_map.size()>0)
	std::cout << "double ended readout on layer2\n";
	int modID_x = fCrtutils.MacToAuxDetID(infn.mac5s, 0);
	int modID_b = fCrtutils.MacToAuxDetID(informationB[i+1].mac5s, 0);
	std::cout << "attempting to match mac " << (int)infn.mac5s << " on modID " << modID_x << " to mac " << (int)informationB[i+1].mac5s << " on modID " << modID_b << "\n";
	if (modID_x == modID_b){
	  recoZwithTiming(infn.mac5s, informationB[i+1].mac5s, infn.t0, informationB[i + 1].t0, fPropDelay, zpos_local, is_infoB_deltaT_over50ns);
	  layer2 = true;
	  std::cout << "zpos = zpos_local + center = " << zpos_local << " + " << center.Z() << " = " << zpos_local+center.Z() << "\n";
	  posB = geo::Zaxis()*(zpos_local+center.Z());
	  zposB_vec.push_back(zpos_local+center.Z());
	  recoZwithPE(infn.mac5s, informationB[i+1].mac5s, pesmap, zpos_local_pe);
	  std::cout << "zpos_PE = zpos_local_pe + center = " << zpos_local_pe << " + " << center.Z() << " = " << zpos_local_pe+center.Z() << "\n";
	  //posB = geo::Zaxis()*(zpos_local_pe+center.Z());
	  //zposB_vec.push_back(zpos_local_pe+center.Z());
	  
	}
	else{
	  //int index = -5;
	  //check for matching modIDs in rest of list
	  std::cout << "infoB: modIDs dont match or modID not in layer_map.. look for rest in list\n";// index = "<< index <<"\n";
	  //if(layer2_map.size()>0) {
	  /*for (auto modID : layer2_map){
	    std::cout << "modID " << modID.first << "\n";
	    }*/
	  //checkNextModID(layer1_map, i, informationA, index);
	  
	  int matched_index = checkNextModID(layer2_map, i, informationB);
	  std::cout << "matched_index = "<< matched_index <<"\n";
	  if (matched_index == 0){
	    std::cout << "matched_index = 0 in infoB! continuing...\n";
	    continue;
	  }
	  if (matched_index!=-1) {
	    recoZwithTiming(infn.mac5s, informationB[matched_index].mac5s, infn.t0, informationB[matched_index].t0, fPropDelay, zpos_local, is_infoB_deltaT_over50ns);
	    layer2=true;
	    std::cout << "zpos = zpos_local + center = " << zpos_local << " + " << center.Z() << " = " << zpos_local+center.Z() << "\n";
	    posB = geo::Zaxis()*(zpos_local+center.Z());
	    zposB_vec.push_back(zpos_local+center.Z());
	    recoZwithPE(infn.mac5s, informationB[matched_index].mac5s, pesmap, zpos_local_pe);
	    std::cout << "zpos_PE = zpos_local_pe + center = " << zpos_local_pe << " + " << center.Z() << " = " << zpos_local_pe+center.Z() << "\n";
	    //posB = geo::Zaxis()*(zpos_local_pe+center.Z());
	    //zposB_vec.push_back(zpos_local_pe+center.Z());
	  }
	}
	/*std::cout << "zposB = " << zposB << "\n";
	posB = zposB * geo::Zaxis();
	zposB_vec.push_back(zposB);*/

      }
    }
    //// else if single ended readout
    //else{
    else if(region != "South" && region != "North"){
      std::cout << "single ended readout on layer2! infoB.size = " << informationB.size() << "\n";
      TVector3 postmp = fCrtutils.ChanToWorldCoords(infn.mac5s, infn.channel);
      mac5_2 = infn.mac5s;
      t0_2 = infn.t0;
    }
  }
  // if mutliple FEB coincidences on single layer, take average of all reco'd z-positions                                                      
  if (zposB_vec.size() > 1){
    zposB = 0;
    std::cout << "zposB = " << zposB << " (before loop over zposB_vect)\n";
    std::cout << "ZposB = (";
    for(float i_zpos : zposB_vec){
      zposB += i_zpos;
      std::cout << i_zpos << " + ";
    }
    if (fVerbose)
      mf::LogInfo("CRTHitRecoAlg") << ")/size = " << zposB << " / " << zposB_vec.size()
                << " = " << (zposB /zposB_vec.size()) << "\n";
    zposB = zposB / zposB_vec.size();
    posB =geo::Zaxis()*zposB;
  }
  std::cout << "hitpos.Z() / nz = " << hitpos.Z() << "/" << nz << " = " << hitpos.Z() / nz << "\n";
  std::cout << "hitpos.Z()*1.0/petot = " << hitpos.Z() << "*1.0/" << petot << " = " << hitpos.Z()*1.0/petot << "\n";
  auto const& adsGeo = adGeo.SensitiveVolume(adsid_max);
  geo::Point_t crossfebpos;
  // assign final Hit Positions based on coin and region
  int crossfeb = std::abs(mac5_1 - mac5_2);
  if((crossfeb == 7 or crossfeb == 5)
     and layerA.size()>0 and layerB.size()>0 
     and region != "South" and region != "North"){
    std::cout << "crossfeb hit!\n";
    hitpos.SetX(hitpos.X() * 1.0 / nx);//nominal
    //hitpos.SetX(hitpos.X() * 1.0 / petot);
    hitpos.SetY(hitpos.Y() * 1.0 / ny);//nominal
    //std::cout << "setY = hitpos.Y/petot = " << hitpos.Y()*1.0 << "/" << petot << " = " << hitpos.Y()*1.0/petot << "\n";
    //hitpos.SetY(hitpos.Y()*1.0/petot);
    auto const& adsGeo = adGeo.SensitiveVolume(adsid_max);
    geo::Point_t center = adsGeo.GetCenter();
    recoZwithTiming(mac5_1, mac5_2, t0_1, t0_2, fPropDelay, zpos_local, is_crossfeb_deltaT_over50ns);
    std::cout
      << "crossfeb = abs(mac5_1 - mac5_2) = " << mac5_1
      << " - " << mac5_2 << " = " << crossfeb << "\n";
    std::cout << "zpos = zpos_local + center = " << zpos_local << " + " << center.Z() << " = " << zpos_local+center.Z() << "\n";
    std::cout << " OR zpos = hitpos.Z / petot = " << hitpos.Z() << " / " << petot << " = " << hitpos.Z() / petot << "\n";
    //hitpos.SetY(hitpos.Y() / petot);
    crossfebpos = geo::Zaxis()*zpos_local + center;
    std::cout << "vs. crossfeb pos.Z() = " << crossfebpos.Z() << "\n";
    hitpos.SetZ(crossfebpos.Z());
    //hitpos.SetZ(hitpos.Z() / petot);
    //posB = geo::Zaxis()*(zpos_local+center.Z());
    //zposB_vec.push_back(zpos_local+center.Z());
    if(is_crossfeb_deltaT_over50ns) // if delta_t > 50, add a dummy value of 100 to error to seperate. REVISIT!
      hitpointerr[2] = int(2.3/fPropDelay) + 100;
    else hitpointerr[2] = int(2.3/fPropDelay);
    std::cout << "zerr = timeRes/propDelay = " << 2.3 << " / " << fPropDelay << " = " << 2.3/fPropDelay << " cm \n";
  }  
  else if (layer1 and layer2 and // if 2 layer readout on E/W modules
      region != "South" and region != "North" ){
    hitpos.SetX(hitpos.X() * 1.0 / nx); //nominal
    //hitpos.SetX(hitpos.X() * 1.0 / petot);
    hitpos.SetY(hitpos.Y() * 1.0 / ny); //nominal 
    std::cout << "setY = hitpos.Y/petot = " << hitpos.Y()*1.0 << "/" << petot << " = " << hitpos.Y()*1.0/petot << "\n";
    //hitpos.SetY(hitpos.Y() * 1.0 / petot);
    //hitpos.SetY(hitpos.Y()*1.0/petot);
    float avg = 0.5 * (posA.Z() + posB.Z());
    std::cout << "avg (posA.Z + posB.Z)/2 = (" << posA.Z() << " + " << posB.Z() << ")/2 = " << avg << "\n";
    std::cout << " OR zpos = hitpos.Z / petot = " << hitpos.Z() << " / " << petot << " = " << hitpos.Z() / petot << "\n";
    hitpos.SetZ(avg); //nominal
    /*float pe_s = -1, pe_n = -1;
    for (int i=0; i<(int)macs.size(); i++){
      int mac = (int)macs.at(i);
      if(mac % 2){// if mac odd-->North side
	for (auto &p : pesmap[(int)mac]){
	  pe_n = p.second;
	}
      }
      else { // if mac even-->South side
	for (auto &p : pesmap[(int)mac]){
	  pe_s = p.second;
	}
      }
    }*/
    //std::cout << "PE_s = " 

    //std::cout << " OR z = z_PE + center.Z = L_atten*ln(PE_s/PE_n)/2 = 1028 cm * ln (" << pe_s << " / " << 
    //hitpos.SetZ(hitpos.Z() / petot);
    //fCrtutils.MacToAuxDetID(mac, 0);
    if (is_infoA_deltaT_over50ns || is_infoB_deltaT_over50ns){
      if(avg > (center.Z() - adsGeo.HalfLength()) && avg < (center.Z() + adsGeo.HalfLength())){
	// if avg Z is still w/in bounds, add 100
	std::cout << "delta_t > 50 on one of the modules, but final Z still w/in bounds\n";
	hitpointerr[2] = int(sqrt(2)*2.3/fPropDelay) + 100; //add a dummy value of 100 to error to seperate. REVISIT!
      }
      else{
	//if oob, add 200
	hitpointerr[2] = int(sqrt(2)*2.3/fPropDelay) + 200; //add a dummy value of 200 to error to seperate. REVISIT!
      }
    }else{
      hitpointerr[2] = int(sqrt(2)*2.3/fPropDelay);
      std::cout << "zerr = sqrt(2)*timeRes/propDelay = sqrt(2)*" << 2.3 << " / " << fPropDelay << " = " << sqrt(2)*2.3/fPropDelay << " cm \n";

      // for good hits (2 layers and deltaT<50), save dist. between reco Z and end of module and PE value for that channel
      // could make a function to do this, pass the mac, L, pesmap and possibly fOutCSVFile
      // this is a temporary thing so maybe making a function isnt necessary? but might want to do it by indv. layer???
      std::cout << "global Z = " << avg << "\n";// for macs: ";
      for (int i=0; i<(int)macs.size(); i++){
	//std::cout << (int)macs.at(i) << " , ";
	int mac = (int)macs.at(i);
	//std::cout << " OR z = z_PE + center.Z = L_atten*ln(PE_s/PE_n)/2 = 1028 cm * ln (" << 
	
	if (pesmap[(int)mac].size() > 1){
	    std::cout << "multiple pesmap entries for mac " << (int)mac << "!!\n";
	  }
	if(mac % 2){// if mac odd
	  for (auto &p : pesmap[(int)mac]){
	    //zrange_max = center.Z() + adsGeo.HalfLength();
	    std::cout << "odd mac: " << (int)mac << "; L = center.Z + halfLength - recoZ = " << center.Z() << " + " << adsGeo.HalfLength() << " - " << avg << " = " << center.Z() + adsGeo.HalfLength() - avg << ",\t";
	    std::cout << "pe = " << p.second << "\n";
	    //if(foutCSVFile) filecsv << center.Z() + adsGeo.HalfLength() - avg << ", " << p.second << "\n";
	  } 
	}
	else{ //if mac even
	  for (auto &p : pesmap[(int)mac]){
	    //zrange_min = center.Z() - adsGeo.HalfLength();
	    std::cout << "even mac: " << (int)mac << "; L = recoZ - (center.Z - HalfLength) = " 
		      << avg << " - (" << center.Z() << " - " << adsGeo.HalfLength() << ") = " 
		      << avg - (center.Z() - adsGeo.HalfLength()) << ",\t";
	    std::cout << "pe = " << p.second << "\n";

	    //if(foutCSVFile) filecsv << avg - (center.Z() - adsGeo.HalfLength()) << ", " << p.second << "\n";
	  }
	}
      }
      std::cout << "\n";
    }
  }
  else if(layer1 and region != "South" and region != "North"){
    hitpos.SetX(hitpos.X() * 1.0 / nx);
    hitpos.SetY(hitpos.Y() * 1.0 / ny);
    //hitpos.SetX(hitpos.X() * 1.0 / petot);
    //hitpos.SetY(hitpos.Y() * 1.0 / petot);
    std::cout << " OR zpos = hitpos.Z / petot = " << hitpos.Z() << " / " << petot << " = " << hitpos.Z() / petot << "\n";
    hitpos.SetZ(posA.Z());
    //hitpos.SetZ(hitpos.Z() / petot);
    //hitpointerr[2] hitpos.Y() / petot= sqrt(2)*2*adsGeo.HalfLength() / sqrt(12);
    if (is_infoA_deltaT_over50ns) 
      hitpointerr[2] = 234 + 100; // dummy value for single ended readouts on one layer
    else
      hitpointerr[2] = 234; // dummy value for single ended readouts on one layer
    //std::cout << "zerr = sqrt(2)*800cm/sqrt(12) = sqrt(2)*" << 2*adsGeo.HalfLength() << "/" << sqrt(12) << " = " << sqrt(2)*2*adsGeo.HalfLength() / sqrt(12) << "\n"; 
  }
  else if(layer2 && region != "South" && region != "North"){
    hitpos.SetX(hitpos.X() * 1.0 / nx);
    hitpos.SetY(hitpos.Y() * 1.0 / ny);
    //hitpos.SetX(hitpos.X() * 1.0 / petot);
    //hitpos.SetY(hitpos.Y() * 1.0 / petot);
    std::cout << " OR zpos = hitpos.Z / petot = " << hitpos.Z() << " / " << petot << " = " << hitpos.Z() / petot << "\n";
    //hitpos.SetZ(hitpos.Z() / petot);
    hitpos.SetZ(posB.Z()); //nominal z
    //std::cout << "zerr = sqrt(2)*800cm/sqrt(12) = sqrt(2)*" << 2*adsGeo.HalfLength() << "/" << sqrt(12) << " = " << sqrt(2)*2*adsGeo.HalfLength() / sqrt(12) << "\n"; 
    //hitpointerr[2] = sqrt(2)*2*adsGeo.HalfLength() / sqrt(12);
    if (is_infoB_deltaT_over50ns) 
      hitpointerr[2] = 234 + 100; // dummy value for single ended readouts on one layer
    else
      hitpointerr[2] = 234; // dummy value for single ended readouts on one layer
  }
  else if(region != "South" && region != "North"){
    //hitpos *= 1.0 / nx; //nominal 
    //hitpos *= 1.0 / petot;
    //std::cout << " OR zpos = hitpos.Z / petot = " << hitpos.Z() << " / " << petot << " = " << hitpos.Z() / petot << "\n";
    hitpos.SetX(hitpos.X() / nx);
    //hitpos.SetX(hitpos.X() * 1.0 / petot);
    hitpos.SetY(hitpos.Y() / nx);    
    //hitpos.SetX(hitpos.X() * 1.0 / nx);
    //hitpos.SetY(hitpos.Y()*1.0/petot);
    //hitpos.SetZ(hitpos.Z() * 1.0 / nx);
    
    if((hitpos.Z() / petot) > (center.Z() - adsGeo.HalfLength()) && (hitpos.Z() / petot) < (center.Z() + adsGeo.HalfLength())){
      std::cout << "single ended readout z within module bounds: z = hitpos.Z / petot = " << hitpos.Z() << " / " << petot << " = " << hitpos.Z() / petot << "\n";
      hitpointerr[2] = int(sqrt(2)*2*adsGeo.HalfLength() / sqrt(12));
      std::cout << "zerr = sqrt(2)*800cm/sqrt(12) = sqrt(2)*" << 2*adsGeo.HalfLength() << "/" << sqrt(12) << " = " << sqrt(2)*2*adsGeo.HalfLength() / sqrt(12) << "\n"; 
    }
    else{
      std::cout << "single ended readout z out of module bounds: z = hitpos.Z / petot = " << hitpos.Z() << " / " << petot << " = " << hitpos.Z() / petot << "\n";
      hitpointerr[2] = int(sqrt(2)*2*adsGeo.HalfLength() / sqrt(12)) + 100; // if avg oob, add 100 to seperate 
      std::cout << "zerr = sqrt(2)*800cm/sqrt(12) + 100 = sqrt(2)*" << 2*adsGeo.HalfLength() << "/" << sqrt(12) << " + 100 = " << sqrt(2)*2*adsGeo.HalfLength() / sqrt(12) + 100 << "\n"; 
    }
    std::cout << "hitpos.Z/petot = " << hitpos.Z() << " / " << petot << " = " << hitpos.Z() / petot << "\n";
    hitpos.SetZ(hitpos.Z() / petot);
    std::cout << " In side CRTs [E/W] (x,y,z) = (" << hitpos[0] << ", " << hitpos[1]
	      << ", " << hitpos[2] <<  ")\n";
  }
  else if(region == "South"){ //for south, only use PE for Z pos 
    std::cout << " OR zpos = hitpos.Z / petot = " << hitpos.Z() << " / " << petot << " = " << hitpos.Z() / petot << "\n";
    hitpos.SetX(hitpos.X() / nx);
    hitpos.SetY(hitpos.Y() / ny);
    //hitpos.SetX(hitpos.X() * 1.0 / petot);
    //hitpos.SetY(hitpos.Y() / petot);
    hitpos.SetZ(hitpos.Z() / petot);
    hitpointerr[0] = int(sqrt(2)*4*adsGeo.HalfHeight()/sqrt(12));
    hitpointerr[2] = int(sqrt(2)*2*adsGeo.HalfWidth1()/sqrt(12)); //zerr should be small 
    
    std::cout << "xerr = sqrt(2)*4*HalfHeight/sqrt(12) = sqrt(2)*" << 4*adsGeo.HalfHeight() << "/" << sqrt(12) << " = " << sqrt(2)*4*adsGeo.HalfHeight()/sqrt(12) << "\n"; //strips are 4cm tall, but each SiPM reads out 2 strips so each chan. is 8cm. 2*halfHeight-->4*halfHeight 
    std::cout << "zerr = sqrt(2)*2*HalfWidth/sqrt(12) = sqrt(2)*" << 2*adsGeo.HalfWidth1() << "/" <<  sqrt(12) << " = " << sqrt(2)*2*adsGeo.HalfWidth1()  / sqrt(12) << "\n";
  }
  else if(region == "North"){
    //hitpos *= 1.0 / nz;
    //hitpos*=1.0/petot;
    std::cout << " OR zpos = hitpos.Z / petot = " << hitpos.Z() << " / " << petot << " = " << hitpos.Z() / petot << "\n";
    //hitpos.SetX(hitpos.X() / nx);
    hitpos.SetX(hitpos.X() * 1.0 / petot);
    hitpos.SetY(hitpos.Y() / petot);
    hitpos.SetZ(hitpos.Z() / petot);
    std::cout << " North region! (x,y,z) = (" << hitpos[0] << ", " << hitpos[1]
	      << ", " << hitpos[2] << '\n';
    hitpointerr[0] = int(sqrt(2)*2*adsGeo.HalfLength()  / sqrt(12)); 
    std::cout << "xerr = sqrt(2)*2*HalfLength/sqrt(12) = sqrt(2)*2*" << adsGeo.HalfLength() << "/" <<  sqrt(12) << " = " << sqrt(2)*2*adsGeo.HalfLength()  / sqrt(12) << "\n";
    hitpointerr[2] = int(sqrt(2)*2*adsGeo.HalfWidth1() / sqrt(12));
    std::cout << "zerr = sqrt(2)*2*HalfWidth/sqrt(12) = sqrt(2)*" << 2*adsGeo.HalfWidth1() << "/" <<  sqrt(12) << " = " << sqrt(2)*2*adsGeo.HalfWidth1()  / sqrt(12) << "\n";
  }
  // i think y-err can be estimated by below for all walls 
  hitpointerr[1] = int(sqrt(2)*4*adsGeo.HalfHeight()/sqrt(12));
  std::cout << "yerr = sqrt(2)*4*HalfHeight/sqrt(12) = sqrt(2)*4*" << adsGeo.HalfHeight() << "/" << sqrt(12) << " = " << sqrt(2)*4*adsGeo.HalfHeight()/sqrt(12) << "\n"; //strips are 4cm tall, but each SiPM reads out 2 strips so each chan. is 8cm. 2*halfHeight-->4*halfHeight 
  if(region!="South" && region != "North"){
    hitpointerr[0] = int(sqrt(2)*2*adsGeo.HalfWidth1()/sqrt(12));
    std::cout << "xerr = sqrt(2)*2*halfWidth/sqrt(12) = sqrt(2)*" << 2*adsGeo.HalfWidth1() << "/" << sqrt(12) << " = " << sqrt(2)*2*adsGeo.HalfWidth1()/sqrt(12) << "\n"; 
  }

  hitpoint[0] = hitpos.X();
  hitpoint[1] = hitpos.Y();
  hitpoint[2] = hitpos.Z();
  // ------ End Hit position reconstruction ------

  // ------ Timing reco ------
  // time stamp averaged over all FEBs
  uint64_t thit = 0., t1hit = 0.;
  ttrigs.resize(std::distance(ttrigs.begin(), std::unique(ttrigs.begin(), ttrigs.end())));
  t1trigs.resize(std::distance(t1trigs.begin(), std::unique(t1trigs.begin(), t1trigs.end())));
  if (!ttrigs.empty()) {
    uint64_t const offset = ttrigs[0];
    int64_t rel_thit = 0;
    for (uint64_t const t : ttrigs) {
      rel_thit += t - offset;
      rel_thit -= RegionDelay(region);
    }
    if (fVerbose)
      mf::LogVerbatim("CRTHitRecoAlg: ")
	<< "Average: offset + rel_thit / ttrigs.size = " << offset << " + "
	<< rel_thit << " / " << uint64_t(ttrigs.size()) << " = ";

    
    thit = offset + // average relative portion of timestamp of CRT data products in the CRT Hit
      (int64_t)(rel_thit / int64_t(ttrigs.size())); 
  } else thit = 0;

  for (double const t1 : t1trigs)
    t1hit += t1 - uint64_t(400. * fPropDelay) - fSiPMtoFEBdelay;

  t1hit = t1hit / uint64_t(t1trigs.size());
  Long64_t thit1;
  if (fData) thit1=(Long64_t)(thit-GlobalTrigger[(int)macs.at(0)]);
  else thit1 = thit - fGlobalT0Offset;
 
  std::cout  << "End of make side hit, z = " << hitpoint[2] << " is saved from macs: ( ";
  for (int i=0; i<(int)macs.size(); i++){
    if(fVerbose) std::cout << (int)macs.at(i) << " , ";
  }
  if(fVerbose) std::cout << ") in region " << region << "\n";
  
  std::cout << "is_infoA_deltaT_over50ns = " << is_infoA_deltaT_over50ns << ", is_infoB_deltaT_over50ns = " << is_infoB_deltaT_over50ns << ", is_crossfeb_deltaT_over50ns = " << is_crossfeb_deltaT_over50ns << "\n";
  //// if double ended readout
  //// else if single ended readout

  // else if single ended readout, assign Z
  std::cout << "(nx,ny,nz) = (" << nx<< "," << ny <<"," << nz << ")\n";
  std::cout << "Filling CRT hit with";
  std::cout << "\n\tmacs.size = " << macs.size() 
	    << "\n\tpesmap.size = " << pesmap.size() 
	    << "\n\tpetot = " << petot 
	    << "\n\tplane, region = " << plane << ", " << region << ", center = " << center 
	    << "\n\t(x,y,z) = (" <<  hitpoint[0] << ", " << hitpoint[1] << ", " << hitpoint[2] 
	    << ")\n\t(xerr, yerr, zerr) = (" 
	    << hitpointerr[0] << ", " << hitpointerr[1] << ", " << hitpointerr[2] 
	    << ")\n\tthit, T0 = " << thit << ", T1 = " << thit1 << "\n";


  // generate hit
  CRTHit hit = FillCRTHit(macs, pesmap, petot, thit, thit1, plane, hitpoint[0],
                          hitpointerr[0], hitpoint[1], hitpointerr[1],
                          hitpoint[2], hitpointerr[2], region);

  return hit;
}
//-----------------------------------------------------------------------------------

//-----------------------------------------------------------------------------
bool CRTHitRecoAlg::IsEmptyHit(CRTHit hit) {
  if (hit.feb_id.empty() && hit.pesmap.empty() && hit.peshit == 0 &&
      hit.ts0_ns == 0 && hit.ts1_ns == 0 && hit.plane == 0 && hit.x_pos == 0 &&
      hit.x_err == 0 && hit.y_pos == 0 && hit.y_err == 0 && hit.z_pos == 0 &&
      hit.z_err == 0 && hit.tagger == "")
    return true;

  return false;
}
//-----------------------------------------------------------------------------
bool containsPair(const vector<pair<int, float>>& vec, const pair<int, float>& p) {
  for (const auto& item : vec) {
    if (item == p) {
      return true;
    }
  }
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
      if (i.second == Mode)
        modecounter++;
      else if (i.second != isnewMode) {
        isnewMode = i.second;
        isnewmodecounter = 1;
      } else if (i.second == isnewMode) {
        isnewmodecounter++;
        if (isnewmodecounter > modecounter) {
          Mode = isnewMode;
          modecounter = isnewmodecounter;
        }
      }
    } else {
      isFirst = false;
      Mode = i.second;
      modecounter++;
    }
  }
  return Mode;
}
//-----------------------------------------------------------------------------
// Function to sort and match FEBs on opposite ends of same module
void CRTHitRecoAlg::groupByLayer(const std::vector<uint8_t>& layer, std::unordered_map<int, std::pair<uint8_t, uint8_t>>& layer_map, const char* layer_name) {
  std::unordered_map<int, std::vector<uint8_t>> modID_to_macs;
  //std::set<std::pair<uint8_t, uint8_t>> unique_pairs;
  // Populate the modID to mac addresses mapping
  for (uint8_t mac : layer) {
    int modID = fCrtutils.MacToAuxDetID(mac, 0);
    modID_to_macs[modID].push_back(mac);
  }

  // Process each modID to find pairs of mac addresses
  for (const auto& entry : modID_to_macs) {
    int modID = entry.first;
    const std::vector<uint8_t>& macs = entry.second;

    // Check all pairs within the same modID
    for (size_t i = 0; i < macs.size(); ++i) {
      for (size_t j = i + 1; j < macs.size(); ++j) {
	uint8_t mac_i = macs[i];
	uint8_t mac_j = macs[j];
	if (mac_i != mac_j) { // Different macs on the same module
	  layer_map[modID] = std::make_pair(mac_i, mac_j);
	  /*auto pair = std::make_pair(std::min(mac_i, mac_j), std::max(mac_i, mac_j));
	  if (unique_pairs.find(pair) == unique_pairs.end()) {
	    std::cout << "macs " << (int)mac_i << ", " << (int)mac_j << " are both on modID " << modID << ", save to " << layer_name << "_map \n";
	    layer_map[modID] = pair;
	    unique_pairs.insert(pair);
	  }
	  else{
	    std::cout << "found a duplicate pair!!\n";
	    }*/
	}
      }
    }
  }
}
//-----------------------------------------------------------------------------
void CRTHitRecoAlg::processChannel(int chan, uint8_t mac, art::Ptr<CRTData> data, std::vector<info>& information, 
				   float& petot, map<uint8_t, vector<pair<int, float>>>& pesmap, 
				   TVector3& hitpos, int& nx, int& ny, int& nz, 
				   int& adsid_max, float& pemax, TVector3& postrig) {
  std::pair<double, double> const chg_cal = fChannelMap->getSideCRTCalibrationMap(mac, chan);
  float pe = (data->fAdc[chan] - chg_cal.second) / chg_cal.first;
  if (pe <= fPEThresh) return;
  TVector3 postmp = fCrtutils.ChanToWorldCoords(mac, chan);
  int adsid = fCrtutils.ChannelToAuxDetSensitiveID(mac, chan);
  
  //check for duplicate (chan,pe) pairs before filling pesmap 
  pair<int, float> newPair = std::make_pair(chan, pe);
  if (pesmap.find(mac) == pesmap.end() || !containsPair(pesmap[mac], newPair)) {
    
    //save unique (chan,PE) pairs to pesmap 
    pesmap[mac].push_back(newPair);
    // add PE info
    petot += pe;

    std::cout<< "pe = (adc[chan]-QPed) / QSlope = (" << data->fAdc[chan] << " - " << chg_cal.second << ")/" << chg_cal.first << " = " << pe << "\n";
    std::cout
      << " (mac5, channel, gain, pedestal, adc, pe, ts0) = ("
      << (int)mac << ", " << chan << ", " << chg_cal.first << ", " << chg_cal.second << ","
      << data->fAdc[chan] << "," << pe << ", " << data->fTs0 << ")\n";
    
    /*bool isNewPair = true;
  if (pesmap.find(mac) != pesmap.end() && containsPair(pesmap[mac], newPair)) {
    isNewPair = false;
    std::cout << "Pair (chan, pe) already exists for this MAC." << std::endl;
    }*/

  //pesmap[mac].push_back(std::make_pair(chan,pe));
    information.push_back({mac, chan, data->fTs0, postmp, adsid}); //orginally push back info here, but try above instead

    // save channel position 
    //int adid = fCrtutils.MacToAuxDetID(mac, chan);
    int adid = fCrtutils.MacToAuxDetID(mac, 0); // passing "0" instead of "chan" is needed for the south wall..
    auto const& adGeo = fGeometryService->AuxDet(adid);  
    auto const& adsGeo = adGeo.SensitiveVolume(adsid); 
    string region = fCrtutils.GetAuxDetRegion(adid);
    int layer = fCrtutils.GetMINOSLayerID(adid); // layer ID
    std::cout<< " adid(mac, chan) =  adid(" << (int)mac << "," << chan << ") = " << adid << "\n";
    std::cout<< " postmp(x,y,z) = (" << postmp.X() << ", " << postmp.Y() << ", " << postmp.Z() << ")\n";
    
    // East/West Walls (all strips along z-direction) or                                                                                     
    // North/South inner walls (all strips along x-direction)                                               
    // All the horizontal layers measure Y first, 
    if (!(region == "South" && layer == 1)) {
      if (region != "South") {  // region is E/W/N
	std::cout<< " setX = postmp.X + hitpos.X = " << 1.0 * postmp.X() << " + " << hitpos.X() << " = " << 1.0 * postmp.X() + hitpos.X() << "\n";
	std::cout<< " setX = pe*postmp.X + hitpos.X = " << pe << "*" << postmp.X() << " + " <<  hitpos.X() << " = " << pe * postmp.X() << " + " << hitpos.X() << " = " << pe * postmp.X() + hitpos.X() << "\n";
	std::cout<< " setY = pe*postmp.Y + hitpos.Y = " << pe << "*" << postmp.Y() << " + " << hitpos.Y() <<  " = "<<  pe * postmp.Y() << " + " << hitpos.Y() << " = " << pe * postmp.Y() + hitpos.Y() << "\n";
	hitpos.SetX(1.0 * postmp.X() + hitpos.X()); // nominal
	hitpos.SetY(1.0 * postmp.Y() + hitpos.Y()); //nominal 
	//hitpos.SetX(pe * postmp.X() + hitpos.X());
	//hitpos.SetY(pe * postmp.Y()+hitpos.Y());
	nx++;
	ny++;
    }
      else{
      //auto const& adsGeo = adGeo.SensitiveVolume(adsid); 
	std::cout<< " horizonal strips in South wall: \n";
	std::cout<< " halfWidth = " << adsGeo.HalfWidth1() << ", halfLength = " << adsGeo.HalfLength() << ", halfHeight = " << adsGeo.HalfHeight() << "\n";
	std::cout<< " setY = postmp.Y + hitpos.Y = " << 1.0 * postmp.Y() << " + " << hitpos.Y() << " = " << 1.0 * postmp.Y() + hitpos.Y() << "\n";
	hitpos.SetY(1.0 * postmp.Y() + hitpos.Y());
	ny++;
      }
    }
    else{  // else vertical strips in South wall   
      //auto const& adsGeo = adGeo.SensitiveVolume(adsid); 
      std::cout<< " vertical strips in South wall \n";
      std::cout<< " halfWidth = " << adsGeo.HalfWidth1() << ", halfLength = " << adsGeo.HalfLength() << ", halfHeight = " << adsGeo.HalfHeight() << "\n";
      std::cout<< " setX = postmp.X + hitpos.X = " << 1.0 * postmp.X() << " + " << hitpos.X() << " = " << 1.0 * postmp.X() + hitpos.X() << "\n";
      //std::cout<< " setX = pe*postmp.X + hitpos.X = " << pe << "*" << postmp.X() << " + " <<  hitpos.X() << " = " << pe * postmp.X() << " + " << hitpos.X() << " = " << pe * postmp.X() + hitpos.X() << "\n";
      hitpos.SetX(1.0 * postmp.X() + hitpos.X());
    //hitpos.SetX(pe * postmp.X() + hitpos.X());
      nx++;
    }
    std::cout<< " setZ = postmp.Z + hitpos.Z = " << 1.0 * postmp.Z() << " + " << hitpos.Z() << " = " << 1.0 * postmp.Z() + hitpos.Z() << "\n";
    //hitpos.SetZ(1.0 * postmp.Z() + hitpos.Z());
    nz++;
    if(region == "South" || region == "North"){
      std::cout<< " setZ = pe*postmp.Z + hitpos.Z = " << pe << "*" << postmp.Z() << " + " << hitpos.Z() <<  " = "<<  pe * postmp.Z() << " + " << hitpos.Z() << " = " << pe * postmp.Z() + hitpos.Z() << "\n";
      hitpos.SetZ(pe * postmp.Z() + hitpos.Z()); // use this for North and South 
    }
    else{
      hitpos.SetZ(1.0 * postmp.Z() + hitpos.Z());
      // for east and west walls, try using PE to assign a Z based on what side that module is on 
      //recoZwithPE(mac,pe,zpos);
      if(mac % 2 == 0){ // even mac, calculate Z from south end 
	std::cout<< " even mac " << (int)mac << ", postmp.Z() = " << postmp.Z() << "\n";
	std::cout<< " FEB is located at (center - halfLength = " << postmp.Z() << " - " << adsGeo.HalfLength() << ") = " << postmp.Z() - adsGeo.HalfLength() << "\n";
	std::cout<< "L from fit: L = ln(" << pe << ") - " << 3.90439 << " / " << -9.72707e-04 << " = (" << log(pe) << " - " << 3.90439 << ")/ " <<-9.72707e-04 <<" = " << (log(pe) - 3.90439)/-9.72707e-04 << "\n";
	std::cout<< "z from fit: evenEnd + L = " << postmp.Z() - adsGeo.HalfLength() << " + " << (log(pe) - 3.90439)/-9.72707e-04 << " = " 
		 << postmp.Z() - adsGeo.HalfLength() + (log(pe) - 3.90439)/-9.72707e-04 << " cm.\n";
	float zpos_fromfit = postmp.Z() - adsGeo.HalfLength() + (log(pe) - 3.90439)/-9.72707e-04;
	hitpos.SetZ(pe * zpos_fromfit + hitpos.Z()); // make sure to use petot to setZ if using this variable
      }
      else{
	std::cout<< "odd mac " << (int)mac << ", postmp.Z() = " << postmp.Z() << "\n";
	std::cout<< "FEB is located at (center + halfLength = " << postmp.Z() << " + " << adsGeo.HalfLength() << ") = " << postmp.Z() + adsGeo.HalfLength() << "\n";
	std::cout<< "L from fit: L = ln(" << pe << ") - " << 3.90439 << " / " << -9.72707e-04 << " = (" << log(pe) << " - " << 3.90439 << ")/ " <<-9.72707e-04 <<" = " << (log(pe) - 3.90439)/-9.72707e-04 << "\n";
	std::cout<< "z from fit: oddEnd - L = " << postmp.Z() + adsGeo.HalfLength() << " - " << (log(pe) - 3.90439)/-9.72707e-04 << " = " 
		 << postmp.Z() + adsGeo.HalfLength() - (log(pe) - 3.90439)/-9.72707e-04 << " cm.\n";
	float zpos_fromfit = postmp.Z() + adsGeo.HalfLength() - (log(pe) - 3.90439)/-9.72707e-04;
	hitpos.SetZ(pe * zpos_fromfit + hitpos.Z());
	
      }
    }
    if(pe>pemax){
      adsid_max = adsid;
      pemax = pe;
      postrig = postmp; // postrig seems unneeded throughout all CRTHitReco? even makeTop/makeBottom
    }
  }
  else{
    std::cout << "Pair (chan, pe) already exists for this MAC." << std::endl; 
  }
}

//-----------------------------------------------------------------------------
void CRTHitRecoAlg::recoZwithTiming(uint8_t mac1, uint8_t mac2, uint64_t t0_1, uint64_t t0_2, float fPropDelay, float& zpos, bool& is_deltaT_over50ns){
  //-----------------------------------------------------------
  // Matching between both end readouts for full length module
  //-----------------------------------------------------------
  // 1. The axis you are measuring on, with the correct sign (d)
  // 2. The centre of the detector (c)
  // 3. The absolute position is c + d · (x - L/2),
  //   that is a displacement from the centre c by a distance x - L/2 (relative
  //   to the centre position L/2). we are determining along the z which is here
  //   d = z (aligned with the long axis of the ICARUS detector). The centre of
  //   the module you are reading is adsGeo.GetCenter() so adsGeo.GetCenter() +
  //   geo::Zaxis() * (x - adsGeo.HalfLength()) may work. We can define a
  //   different x, from the middle of the detector:

  // T1 = t_hit + (L/2 + x)/v_propagation;
  // T2 = t_hit + (L/2 - x)/v_propagation;
  // double deltat(T1-T2) = 2*x/v_propagation;
  //---------------------------------------------------------------------
  uint64_t t0_s = 0, t0_n = 0;
  if ((int)mac1 % 2 == 0) //if mac is even, south side
    t0_s = t0_1;
  else 
    t0_s = t0_2;

  if((int)mac2 % 2 != 0) //if mac is odd, north side
    t0_n = t0_2;
  else 
    t0_n = t0_1;

  zpos = 0.5 * (int64_t(t0_s - t0_n) / fPropDelay);
  if(fVerbose) std::cout
    << "---\n\tFEB : (mac_1,mac_2 = " << (int)mac1 << ","
    << (int)mac2<< "), delta_t = t1 - t2 = "
    << t0_s << " - " << t0_n << " = " << int64_t(t0_s - t0_n) 
    << "\n\tlocal Zpos = .5*(delta_t)/prop = .5*(" 
    <<  int64_t(t0_s - t0_n) << ")/" << fPropDelay << " = " << zpos << "\n";
  
  // Check if deltaT > 50
  is_deltaT_over50ns = false;
  if(std::abs(int64_t(t0_s - t0_n))>=50){
    is_deltaT_over50ns = true;
    std::cout<< "\tdelta_t >= 50! delta_t = " 
	      << int64_t(t0_s - t0_n) << "\n";
  }
}
//-----------------------------------------------------------------------------
int CRTHitRecoAlg::checkNextModID(std::unordered_map<int,std::pair<uint8_t,uint8_t>> layer_map, int i, const std::vector<info>& information){
  int index = -1;
  // Check against layer1_map if no direct match is found
  const auto& infn = information[i];
  int modID = fCrtutils.MacToAuxDetID(infn.mac5s, 0);
  /*for (const auto& key : layer_map) {
    std::cout << "does modID = " << modID << " from mac " << (int)infn.mac5s<< " match modID in layer_map "<< key.first << "?\n";
    std::cout << "macs( " << (int)key.second.first << ", " << (int)key.second.second << ") on modID " << key.first << "\n";
    if (modID==key.first and infn.mac5s != key.second.first){
      std::cout << "different macs on same modID!\n";
    }
    else{
      std::cout << "same mac or different modIDs, break..!\n";
      break;
    }
    // if (infn.mac5s == key.second.first) {
      std::cout << "mac " << (int)infn.mac5s << " matches layer1_map 1st mac " << (int)key.second.first << "\n";
      break;
      }//
      // doesnt seem like the above actually does anything?
  }*/

  // Find matching modID in the informationA vector
  for (size_t j = 0; j < information.size(); ++j) {
    int modID_j = fCrtutils.MacToAuxDetID(information[j].mac5s, 0);
    if (modID == modID_j && infn.mac5s != information[j].mac5s) {
      std::cout<< "found a match!!!! modID(infn.macs= " << (int)infn.mac5s << ") = " << modID << " = modID(info[j=" << j << "].macs=" << (int)information[j].mac5s << " )= " << modID_j << "\n";
      index = j;
      break;
    }
  }
  return index;
}
//-----------------------------------------------------------------------------
void CRTHitRecoAlg::recoZwithPE(uint8_t mac1, uint8_t mac2, map<uint8_t, vector<pair<int, float>>> pesmap, float& zpos){
  float pe_s, pe_n;
  if ((int)mac1 % 2 == 0) //if mac is even, south side
    for (auto &p : pesmap[(int)mac1]){
      pe_s = p.second;
      std::cout<< "\t mac " << (int)mac1 << " on south side w PE=" << p.second << "\n";
    }
  else
    for (auto &p : pesmap[(int)mac1]){
      pe_s = p.second;
      std::cout<< "\t mac " << (int)mac1 << " on south side w PE=" << p.second << "\n";
    }
  if((int)mac2 % 2 != 0) //if mac is odd, north side
    for (auto &p : pesmap[(int)mac2]){
      pe_n = p.second;
      std::cout<< "\t mac " << (int)mac2 << " on north side w PE=" << p.second << "\n";
    }
  else
    for (auto &p : pesmap[(int)mac2]){
      pe_n = p.second;
      std::cout<< "\t mac " << (int)mac2 << " on north side w PE=" << p.second << "\n";
    }
  std::cout<< "---\n\tFEB : (mac_1,mac_2 = " << (int)mac1 << ","
	    << (int)mac2<< "), (PE_s, PE_n = " << pe_s << ", " << pe_n << "\n";
  std::cout<< "\n\t local Zpos with PE = L_atten*ln(PE_s/PE_n)/2 = 1028 cm * ln (" << pe_s << " / " << pe_n << ")/2 = 1028cm/2 * " << log(pe_s/pe_n) << " = "<< 1028/2 * log(pe_s/pe_n) << " cm. \n";
  zpos = 1028/2 * log(pe_s/pe_n);
    //pe_s = 

}
//-----------------------------------------------------------------------------
/*void CRTHitRecoAlg::processInformation(const std::vector<info>& information, fGeometryService& adGeo, int fVerbose, bool& layer, uint64_t& t1, uint64_t& t2, TVector3& pos, float fPropDelay) {
  for (size_t i = 0; i < information.size(); ++i) {
    const auto& infn = information[i];
    const auto& adsGeo = adGeo.SensitiveVolume(infn.strip);
    TVector3 center = adsGeo.GetCenter();

    if (i < information.size() - 1 && (int)infn.mac5s != (int)information[i + 1].mac5s) {
      layer = true;
      if ((int)infn.mac5s % 2 == 0)
	t1 = infn.t0;
      else
	t1 = information[i + 1].t0;

      if ((int)information[i + 1].mac5s % 2 != 0)
	t2 = information[i + 1].t0;
      else
	t2 = infn.t0;

      float zaxixpos = 0.5 * (int64_t(t1 - t2) / fPropDelay);
      pos = adsGeo.GetCenter() + geo::Zaxis() * zaxixpos;

      std::cout
	<< "---\n\tFEB A: (mac_1,mac_2 = " << (int)infn.mac5s << ","
	<< (int)informationA[index].mac5s<< "), delta_t = t1_1 - t1_2 = "
	<< t1_1 << " - " << t1_2 << " = " << int64_t(t1_1 - t1_2) << "\n"
	<< "\tFEB A z hit pos = .5*(delta_t)/prop + center = .5*("
	<< int64_t(t1_1 - t1_2) << ")/ " << fPropDelay << " + "
	<< 1.0*adsGeo.GetCenter().Z() << " = " << posA.Z() << "\n";
    }
  }
  }*/
//-----------------------------------------------------------------------------
/*void CRTHitRecoAlg::processHitInformation(const std::vector<info>& informationA, const std::vector<info>& informationB, CRTGeo& adGeo, int fVerbose, bool& layer1, uint64_t& t1_1, uint64_t& t1_2, TVector3& posA, bool& layer2, uint64_t& t2_1, uint64_t& t2_2, TVector3& posB, float fPropDelay) {
  processInformation(informationA, adGeo, fVerbose, layer1, t1_1, t1_2, posA, fPropDelay);
  processInformation(informationB, adGeo, fVerbose, layer2, t2_1, t2_2, posB, fPropDelay);
  }*/
//-----------------------------------------------------------------------------
/*void CRTHitRecoAlg::processChannels(int mac, int adid, int i, const std::string& layer, 
				    std::unordered_map<int, std::pair<uint8_t, uint8_t>>& layer_map,
				    std::vector<uint8_t>& layer_data, std::vector<info>& information) {
  // Check if MAC is on modID saved in layer_map
  std::cout << "check: is mac " << mac << ", modID " << adid << " on same modID saved in " << layer << "? " << layer << "_map.size = "
              << layer_map.size() << "\n";
  auto processChannelData = [&](int chan) {
    std::pair<double, double> const chg_cal = fChannelMap->getSideCRTCalibrationMap(mac, chan);
    float pe = (data->fAdc[chan] - chg_cal.second) / chg_cal.first;
    if (pe > fPEThresh) {
      nabove++;
      std::cout << "\ni=" << i << ", feb" << layer << " (mac5, channel, gain, pedestal, adc, pe, t0) = ("
		<< mac << ", " << chan << ", " << chg_cal.first << ", " << chg_cal.second << ","
		<< data->fAdc[chan] << "," << pe << ", " << data->fTs0 << ") (using " << layer << ")\n";
      int adsid = fCrtutils.ChannelToAuxDetSensitiveID(mac, chan);
      TVector3 postmp = fCrtutils.ChanToWorldCoords(mac, chan);
      information.push_back({mac, chan, data->fTs0, postmp, adsid});
    }
  };
  // Process if layer_map is not empty
  if (!layer_map.empty()) {
    for (const auto& key : layer_map) {
      if (adid == key.first) {
	std::cout << "modID = " << adid << " matches modID = " << key.first << " in " << layer << "!\n";
	for (int chan = 0; chan < 32; chan++) {
	  processChannelData(chan);
	}
            } else {
                std::cout << "modID = " << adid << " doesn't match modID = " << key.first << " in map\n";
            }
        }
    } else {
        std::cout << layer << "_map.size = 0! " << layer << ".size = " << layer_data.size() << "\n";
        if (!layer_data.empty()) {
            for (int chan = 0; chan < 32; chan++) {
                processChannelData(chan);
            }
        }
    }
    }
//-----------------------------------------------------------------------------
void CRTHitRecoAlg::processLayerData(int i, std::vector<uint8_t>& layerA, std::vector<uint8_t>& layerB, int adid, 
                      std::vector<std::pair<int, std::pair<uint8_t, uint8_t>>>& layer1_map, 
                      std::vector<std::pair<int, std::pair<uint8_t, uint8_t>>>& layer2_map, 
				     std::vector<CRTHitRecoAlg::info>& informationA, std::vector<CRTHitRecoAlg::info>& informationB) {
  int mac = static_cast<int>(macs.back());
  if (mac == static_cast<int>(layerA[i])) {
    processChannels(mac, adid, i, "layer1", layer1_map, layerA, informationA);
  } else if (mac == static_cast<int>(layerB[i])) {
    processChannels(mac, adid, i, "layer2", layer2_map, layerB, informationB);
  }
  }*/
//-----------------------------------------------------------------------------
