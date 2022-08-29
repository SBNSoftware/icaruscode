#ifndef CRTHITRECOALG_H_SEEN
#define CRTHITRECOALG_H_SEEN
////////////////////////////////////////////////////////////////////
// CRTHitRecoAlg.h
//
// Functions for CRT hit reconstruction
// Chris Hilgenberg (Chris.Hilgenberg@colostate.edu), June 2020
////////////////////////////////////////////////////////////////////

// framework
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "canvas/Persistency/Common/FindManyP.h"
#include "larcorealg/CoreUtils/NumericUtils.h"

// LArSoft
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "larcore/Geometry/AuxDetGeometry.h"

// Utility libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"
#include "cetlib/pow.h" // cet::sum_of_squares()

//icaruscode includes
#include "sbnobj/Common/CRT/CRTHit.hh"
#include "sbnobj/ICARUS/CRT/CRTData.hh"
#include "icaruscode/CRT/CRTUtils/CRTCommonUtils.h"

#include "icaruscode/Decode/ChannelMapping/IChannelMapping.h"
#include "icaruscode/Decode/ChannelMapping/IICARUSChannelMap.h"
#include "icaruscode/Decode/ChannelMapping/ICARUSChannelMapProvider.h"
// c++
#include <iostream>
#include <stdio.h>
#include <sstream>
#include <vector>
#include <map>
#include <utility>
#include <cmath> 
#include <memory>
#include <cstdint>
// ROOT
#include "TVector3.h"
#include "TGeoManager.h"

using std::vector;
using std::map;
using std::pair;
using std::string;

namespace icarus {
 namespace crt {
    class CRTHitRecoAlg;
 }
}

class icarus::crt::CRTHitRecoAlg {

 public:
    
  using CRTData = icarus::crt::CRTData;
  using CRTHit = sbn::crt::CRTHit;

  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    fhicl::Atom<bool> UseReadoutWindow {
      Name("UseReadoutWindow"),
	Comment("Only reconstruct hits within readout window")
	};
    fhicl::Atom<double> QPed {
      Name("QPed"),
	Comment("Pedestal offset [ADC]")
	};
    fhicl::Atom<double> QSlope {
      Name("QSlope"),
	Comment("Pedestal slope [ADC/photon]")
	};
    fhicl::Atom<bool> Verbose {
      Name("Verbose"),
	Comment("Output verbosity")
	};
    fhicl::Atom<bool> Data {
      Name("Data"),
	Comment("choose data/mc")
	};
    fhicl::Atom<bool> outCSVFile {
      Name("outCSVFile"),
	Comment("Output a csv file")
	};
    fhicl::Atom<std::string> CSVFile {
      Name("CSVFile"),
	Comment("file name for output CSV File")
        };
    fhicl::Atom<double> PropDelay {
      Name("PropDelay"),
	Comment("group velocity in WLS fiber [ns/cm]")
	};
    fhicl::Atom<double> PEThresh {
      Name("PEThresh"),
	Comment("threshold in photoelectrons above which charge amplitudes used in hit reco")
	};
    fhicl::Atom<double> topGain {
      Name("topGain"),
        Comment("Dummy Gain value for Top CRT")
        };
    fhicl::Atom<double> topPed {
      Name("topPed"),
        Comment("Dummy Pedestal value for Top CRT")
        };
    fhicl::Atom<uint64_t> SiPMtoFEBdelay {
      Name("SiPMtoFEBdelay"),
	Comment("Delay for SiPM to FEB signal correction 11.6 [ns]")
	};
    fhicl::Atom<uint64_t> CoinWindow {
      Name("CoinWindow"),
	Comment("window for finding side CRT trigger coincidences [ns]")
	};
    fhicl::Atom<uint64_t> CrtWindow {
      Name("CrtWindow"),
	Comment("window for looking data [ns]")
	};
  };//Config

  //constructors
  CRTHitRecoAlg(const Config& config);
  CRTHitRecoAlg(const fhicl::ParameterSet& pset) :
    CRTHitRecoAlg{fhicl::Table<Config>{pset}()} {}
  CRTHitRecoAlg();

  //configure module from fcl file
  void reconfigure(const Config& config);

  //produce CRTHits with associated data indices from input vector of CRTData
  vector<pair<CRTHit, vector<int>>> CreateCRTHits(vector<art::Ptr<CRTData>> crtList);

  //preselection based on charge in a CRTData
  vector<art::Ptr<CRTData>> PreselectCRTData(vector<art::Ptr<CRTData>> crtList, uint64_t trigger_timestamp);

  // Function to make filling a CRTHit a bit faster
  CRTHit FillCRTHit(vector<uint8_t> tfeb_id, map<uint8_t, vector<pair<int,float>>> tpesmap,
                    float peshit, uint64_t time0, Long64_t time1, int plane,
                    double x, double ex, double y, double ey, double z, double ez, string tagger);


 private:

  geo::GeometryCore const* fGeometryService;

  CRTCommonUtils* fCrtutils;

  //Params from fcl file
  bool fVerbose;          ///< print info
  bool foutCSVFile;       ///<FCL input: Write a CSV File?
  std::string fCSVFile;   ///<FCL input: file name for output CSV File
  bool fUseReadoutWindow; ///< Only reconstruct hits within TPC readout window
  double fQPed;           ///< Pedestal offset of SiPMs [ADC]
  double fQSlope;         ///< Pedestal slope of SiPMs [ADC/photon]
  double fPropDelay;      ///< propegation time [ns/cm]
  double fPEThresh;       ///< threshold[PE] above which charge amplitudes used in hit reco
  double ftopGain;        ///< Dummy Top CRT Gain Value
  double ftopPed;         ///< Dummy Top CRT Pedestal Value
  uint64_t fSiPMtoFEBdelay; ///< SiPM to FEB cable induced delay: 11.6 [ns]
  uint64_t fCoinWindow;   ///< Coincidence window used for grouping side CRT triggers [ns]
  uint64_t fCrtWindow;    ///< Looking data window within trigger timestamp [ns]
  std::ofstream filecsv;
  bool fData;             ///< look for only data
  const icarusDB::IICARUSChannelMap* fChannelMap = nullptr;

  //Given top CRTData product, produce CRTHit
  CRTHit MakeTopHit(art::Ptr<CRTData> data, ULong64_t GlobalTrigger);
  //Given bottom CRTData product, produce CRTHit
  CRTHit MakeBottomHit(art::Ptr<CRTData> data);
  //Given vector of side CRTData products, produce CRTHit
  CRTHit MakeSideHit(vector<art::Ptr<CRTData>> coinData, ULong64_t GlobalTrigger);
  // Check if a hit is empty
  bool IsEmptyHit(CRTHit hit);

  static  bool compareBytime(art::Ptr<CRTData> const &a, art::Ptr<CRTData> const &b){
    return a->fTs0 < b->fTs0;
  }
}; //class CRTHitRecoAlg


inline ULong64_t GetMode(std::vector<ULong64_t> vector) {

	sort(vector.begin(), vector.end(), std::greater<int>());

	int modecounter = 0;
	int isnewmodecounter = 0;
	ULong64_t Mode = 0;
	ULong64_t isnewMode = 0;
	bool isFirst = true;
	for (auto i : vector) {
		if (!isFirst) {
			if (i == Mode) modecounter++;
			else if (i!=isnewMode) {
				isnewMode = i;
				isnewmodecounter = 1;
			}
			else if (i == isnewMode) {
				isnewmodecounter++;
				if (isnewmodecounter > modecounter) {
					Mode = isnewMode;
					modecounter = isnewmodecounter;
				}
			}
		}
		else {
			isFirst = false;
			Mode = i;
			modecounter++;
		}
	}
	return Mode;
}

typedef struct FEB_delay {
	int HW_mac;
	int SW_mac;
	int SW_modID;
	ULong64_t T0_delay;
	ULong64_t T1_delay;
} FEB_delay;

inline std::map<int, FEB_delay> LoadFEBMap() {

	std::map<int, struct FEB_delay> FEBs;

	struct FEB_delay FEB;
	FEB.HW_mac = 81;
	FEB.SW_mac = 198;
	FEB.SW_modID = 271;
	FEB.T0_delay = 283;
	FEB.T1_delay = 2000309;
	FEBs.insert({ 271, FEB });
	FEB.HW_mac = 119;
	FEB.SW_mac = 197;
	FEB.SW_modID = 270;
	FEB.T0_delay = 298;
	FEB.T1_delay = 2000324;
	FEBs.insert({ 270, FEB });
	FEB.HW_mac = 87;
	FEB.SW_mac = 196;
	FEB.SW_modID = 269;
	FEB.T0_delay = 313;
	FEB.T1_delay = 2000339;
	FEBs.insert({ 269, FEB });
	FEB.HW_mac = 92;
	FEB.SW_mac = 195;
	FEB.SW_modID = 268;
	FEB.T0_delay = 329;
	FEB.T1_delay = 2000355;
	FEBs.insert({ 268, FEB });
	FEB.HW_mac = 180;
	FEB.SW_mac = 194;
	FEB.SW_modID = 267;
	FEB.T0_delay = 344;
	FEB.T1_delay = 2000370;
	FEBs.insert({ 267, FEB });
	FEB.HW_mac = 97;
	FEB.SW_mac = 193;
	FEB.SW_modID = 266;
	FEB.T0_delay = 359;
	FEB.T1_delay = 2000385;
	FEBs.insert({ 266, FEB });
	FEB.HW_mac = 174;
	FEB.SW_mac = 192;
	FEB.SW_modID = 265;
	FEB.T0_delay = 374;
	FEB.T1_delay = 2000400;
	FEBs.insert({ 265, FEB });
	FEB.HW_mac = 238;
	FEB.SW_mac = 178;
	FEB.SW_modID = 251;
	FEB.T0_delay = 390;
	FEB.T1_delay = 2000416;
	FEBs.insert({ 251, FEB });
	FEB.HW_mac = 234;
	FEB.SW_mac = 164;
	FEB.SW_modID = 237;
	FEB.T0_delay = 405;
	FEB.T1_delay = 2000431;
	FEBs.insert({ 237, FEB });
	FEB.HW_mac = 189;
	FEB.SW_mac = 224;
	FEB.SW_modID = 297;
	FEB.T0_delay = 420;
	FEB.T1_delay = 2000446;
	FEBs.insert({ 297, FEB });
	FEB.HW_mac = 190;
	FEB.SW_mac = 223;
	FEB.SW_modID = 296;
	FEB.T0_delay = 436;
	FEB.T1_delay = 2000462;
	FEBs.insert({ 296, FEB });
	FEB.HW_mac = 80;
	FEB.SW_mac = 222;
	FEB.SW_modID = 295;
	FEB.T0_delay = 451;
	FEB.T1_delay = 2000477;
	FEBs.insert({ 295, FEB });
	FEB.HW_mac = 162;
	FEB.SW_mac = 221;
	FEB.SW_modID = 294;
	FEB.T0_delay = 466;
	FEB.T1_delay = 2000492;
	FEBs.insert({ 294, FEB });
	FEB.HW_mac = 64;
	FEB.SW_mac = 220;
	FEB.SW_modID = 293;
	FEB.T0_delay = 482;
	FEB.T1_delay = 2000508;
	FEBs.insert({ 293, FEB });
	FEB.HW_mac = 172;
	FEB.SW_mac = 182;
	FEB.SW_modID = 255;
	FEB.T0_delay = 298;
	FEB.T1_delay = 2000324;
	FEBs.insert({ 255, FEB });
	FEB.HW_mac = 114;
	FEB.SW_mac = 181;
	FEB.SW_modID = 254;
	FEB.T0_delay = 313;
	FEB.T1_delay = 2000339;
	FEBs.insert({ 254, FEB });
	FEB.HW_mac = 100;
	FEB.SW_mac = 180;
	FEB.SW_modID = 253;
	FEB.T0_delay = 328;
	FEB.T1_delay = 2000355;
	FEBs.insert({ 253, FEB });
	FEB.HW_mac = 150;
	FEB.SW_mac = 179;
	FEB.SW_modID = 252;
	FEB.T0_delay = 344;
	FEB.T1_delay = 2000370;
	FEBs.insert({ 252, FEB });
	FEB.HW_mac = 176;
	FEB.SW_mac = 165;
	FEB.SW_modID = 238;
	FEB.T0_delay = 359;
	FEB.T1_delay = 2000385;
	FEBs.insert({ 238, FEB });
	FEB.HW_mac = 67;
	FEB.SW_mac = 151;
	FEB.SW_modID = 224;
	FEB.T0_delay = 374;
	FEB.T1_delay = 2000400;
	FEBs.insert({ 224, FEB });
	FEB.HW_mac = 138;
	FEB.SW_mac = 150;
	FEB.SW_modID = 223;
	FEB.T0_delay = 390;
	FEB.T1_delay = 2000416;
	FEBs.insert({ 223, FEB });
	FEB.HW_mac = 170;
	FEB.SW_mac = 136;
	FEB.SW_modID = 209;
	FEB.T0_delay = 405;
	FEB.T1_delay = 2000431;
	FEBs.insert({ 209, FEB });
	FEB.HW_mac = 101;
	FEB.SW_mac = 122;
	FEB.SW_modID = 195;
	FEB.T0_delay = 420;
	FEB.T1_delay = 2000446;
	FEBs.insert({ 195, FEB });
	FEB.HW_mac = 142;
	FEB.SW_mac = 108;
	FEB.SW_modID = 181;
	FEB.T0_delay = 435;
	FEB.T1_delay = 2000462;
	FEBs.insert({ 181, FEB });
	FEB.HW_mac = 139;
	FEB.SW_mac = 206;
	FEB.SW_modID = 279;
	FEB.T0_delay = 451;
	FEB.T1_delay = 2000477;
	FEBs.insert({ 279, FEB });
	FEB.HW_mac = 185;
	FEB.SW_mac = 207;
	FEB.SW_modID = 280;
	FEB.T0_delay = 466;
	FEB.T1_delay = 2000492;
	FEBs.insert({ 280, FEB });
	FEB.HW_mac = 6;
	FEB.SW_mac = 109;
	FEB.SW_modID = 182;
	FEB.T0_delay = 481;
	// TEMPORARY
//	FEB.T1_delay = 2000447;
	FEB.T1_delay = 2000508;
	FEBs.insert({ 182, FEB });
	FEB.HW_mac = 177;
	FEB.SW_mac = 123;
	FEB.SW_modID = 196;
	FEB.T0_delay = 497;
	FEB.T1_delay = 2000523;
	FEBs.insert({ 196, FEB });
	FEB.HW_mac = 61;
	FEB.SW_mac = 137;
	FEB.SW_modID = 210;
	FEB.T0_delay = 512;
	FEB.T1_delay = 2000538;
	FEBs.insert({ 210, FEB });
	FEB.HW_mac = 125;
	FEB.SW_mac = 183;
	FEB.SW_modID = 256;
	FEB.T0_delay = 298;
	FEB.T1_delay = 2000325;
	FEBs.insert({ 256, FEB });
	FEB.HW_mac = 116;
	FEB.SW_mac = 169;
	FEB.SW_modID = 242;
	FEB.T0_delay = 314;
	FEB.T1_delay = 2000340;
	FEBs.insert({ 242, FEB });
	FEB.HW_mac = 104;
	FEB.SW_mac = 168;
	FEB.SW_modID = 241;
	FEB.T0_delay = 329;
	FEB.T1_delay = 2000355;
	FEBs.insert({ 241, FEB });
	FEB.HW_mac = 91;
	FEB.SW_mac = 167;
	FEB.SW_modID = 240;
	FEB.T0_delay = 344;
	FEB.T1_delay = 2000371;
	FEBs.insert({ 240, FEB });
	FEB.HW_mac = 88;
	FEB.SW_mac = 166;
	FEB.SW_modID = 239;
	FEB.T0_delay = 360;
	FEB.T1_delay = 2000386;
	FEBs.insert({ 239, FEB });
	FEB.HW_mac = 120;
	FEB.SW_mac = 152;
	FEB.SW_modID = 225;
	FEB.T0_delay = 375;
	FEB.T1_delay = 2000401;
	FEBs.insert({ 225, FEB });
	FEB.HW_mac = 132;
	FEB.SW_mac = 138;
	FEB.SW_modID = 211;
	FEB.T0_delay = 390;
	FEB.T1_delay = 2000417;
	FEBs.insert({ 211, FEB });
	FEB.HW_mac = 95;
	FEB.SW_mac = 124;
	FEB.SW_modID = 197;
	FEB.T0_delay = 405;
	FEB.T1_delay = 2000432;
	FEBs.insert({ 197, FEB });
	FEB.HW_mac = 232;
	FEB.SW_mac = 110;
	FEB.SW_modID = 183;
	FEB.T0_delay = 421;
	FEB.T1_delay = 2000447;
	FEBs.insert({ 183, FEB });
	FEB.HW_mac = 165;
	FEB.SW_mac = 208;
	FEB.SW_modID = 281;
	FEB.T0_delay = 436;
	FEB.T1_delay = 2000463;
	FEBs.insert({ 281, FEB });
	FEB.HW_mac = 148;
	FEB.SW_mac = 209;
	FEB.SW_modID = 282;
	FEB.T0_delay = 451;
	FEB.T1_delay = 2000478;
	FEBs.insert({ 282, FEB });
	FEB.HW_mac = 237;
	FEB.SW_mac = 111;
	FEB.SW_modID = 184;
	FEB.T0_delay = 467;
	FEB.T1_delay = 2000493;
	FEBs.insert({ 184, FEB });
	FEB.HW_mac = 102;
	FEB.SW_mac = 125;
	FEB.SW_modID = 198;
	FEB.T0_delay = 482;
	FEB.T1_delay = 2000508;
	FEBs.insert({ 198, FEB });
	FEB.HW_mac = 94;
	FEB.SW_mac = 139;
	FEB.SW_modID = 212;
	FEB.T0_delay = 497;
	FEB.T1_delay = 2000524;
	FEBs.insert({ 212, FEB });
	FEB.HW_mac = 130;
	FEB.SW_mac = 153;
	FEB.SW_modID = 226;
	FEB.T0_delay = 513;
	FEB.T1_delay = 2000539;
	FEBs.insert({ 226, FEB });
	FEB.HW_mac = 181;
	FEB.SW_mac = 184;
	FEB.SW_modID = 257;
	FEB.T0_delay = 284;
	FEB.T1_delay = 2000310;
	FEBs.insert({ 257, FEB });
	FEB.HW_mac = 124;
	FEB.SW_mac = 170;
	FEB.SW_modID = 243;
	FEB.T0_delay = 299;
	FEB.T1_delay = 2000325;
	FEBs.insert({ 243, FEB });
	FEB.HW_mac = 152;
	FEB.SW_mac = 156;
	FEB.SW_modID = 229;
	FEB.T0_delay = 314;
	FEB.T1_delay = 2000341;
	FEBs.insert({ 229, FEB });
	FEB.HW_mac = 98;
	FEB.SW_mac = 155;
	FEB.SW_modID = 228;
	FEB.T0_delay = 329;
	FEB.T1_delay = 2000356;
	FEBs.insert({ 228, FEB });
	FEB.HW_mac = 173;
	FEB.SW_mac = 154;
	FEB.SW_modID = 227;
	FEB.T0_delay = 345;
	FEB.T1_delay = 2000371;
	FEBs.insert({ 227, FEB });
	FEB.HW_mac = 169;
	FEB.SW_mac = 140;
	FEB.SW_modID = 213;
	FEB.T0_delay = 360;
	FEB.T1_delay = 2000387;
	FEBs.insert({ 213, FEB });
	FEB.HW_mac = 144;
	FEB.SW_mac = 126;
	FEB.SW_modID = 199;
	FEB.T0_delay = 375;
	FEB.T1_delay = 2000402;
	FEBs.insert({ 199, FEB });
	FEB.HW_mac = 239;
	FEB.SW_mac = 112;
	FEB.SW_modID = 185;
	FEB.T0_delay = 391;
	FEB.T1_delay = 2000417;
	FEBs.insert({ 185, FEB });
	FEB.HW_mac = 147;
	FEB.SW_mac = 210;
	FEB.SW_modID = 283;
	FEB.T0_delay = 406;
	FEB.T1_delay = 2000433;
	FEBs.insert({ 283, FEB });
	FEB.HW_mac = 105;
	FEB.SW_mac = 211;
	FEB.SW_modID = 284;
	FEB.T0_delay = 421;
	FEB.T1_delay = 2000448;
	FEBs.insert({ 284, FEB });
	FEB.HW_mac = 231;
	FEB.SW_mac = 113;
	FEB.SW_modID = 186;
	FEB.T0_delay = 437;
	FEB.T1_delay = 2000463;
	FEBs.insert({ 186, FEB });
	FEB.HW_mac = 117;
	FEB.SW_mac = 127;
	FEB.SW_modID = 200;
	FEB.T0_delay = 452;
	FEB.T1_delay = 2000478;
	FEBs.insert({ 200, FEB });
	FEB.HW_mac = 126;
	FEB.SW_mac = 141;
	FEB.SW_modID = 214;
	FEB.T0_delay = 467;
	FEB.T1_delay = 2000494;
	FEBs.insert({ 214, FEB });
	FEB.HW_mac = 90;
	FEB.SW_mac = 142;
	FEB.SW_modID = 215;
	FEB.T0_delay = 482;
	FEB.T1_delay = 2000509;
	FEBs.insert({ 215, FEB });
	FEB.HW_mac = 183;
	FEB.SW_mac = 128;
	FEB.SW_modID = 201;
	FEB.T0_delay = 498;
	FEB.T1_delay = 2000524;
	FEBs.insert({ 201, FEB });
	FEB.HW_mac = 241;
	FEB.SW_mac = 114;
	FEB.SW_modID = 187;
	FEB.T0_delay = 513;
	FEB.T1_delay = 2000540;
	FEBs.insert({ 187, FEB });
	FEB.HW_mac = 113;
	FEB.SW_mac = 212;
	FEB.SW_modID = 285;
	FEB.T0_delay = 528;
	FEB.T1_delay = 2000555;
	FEBs.insert({ 285, FEB });
	FEB.HW_mac = 233;
	FEB.SW_mac = 185;
	FEB.SW_modID = 258;
	FEB.T0_delay = 283;
	FEB.T1_delay = 2000310;
	FEBs.insert({ 258, FEB });
	FEB.HW_mac = 164;
	FEB.SW_mac = 171;
	FEB.SW_modID = 244;
	FEB.T0_delay = 299;
	FEB.T1_delay = 2000325;
	FEBs.insert({ 244, FEB });
	FEB.HW_mac = 161;
	FEB.SW_mac = 157;
	FEB.SW_modID = 230;
	FEB.T0_delay = 314;
	FEB.T1_delay = 2000341;
	FEBs.insert({ 230, FEB });
	FEB.HW_mac = 203;
	FEB.SW_mac = 158;
	FEB.SW_modID = 231;
	FEB.T0_delay = 329;
	FEB.T1_delay = 2000356;
	FEBs.insert({ 231, FEB });
	FEB.HW_mac = 122;
	FEB.SW_mac = 159;
	FEB.SW_modID = 232;
	FEB.T0_delay = 345;
	FEB.T1_delay = 2000371;
	FEBs.insert({ 232, FEB });
	FEB.HW_mac = 2;
	FEB.SW_mac = 145;
	FEB.SW_modID = 218;
	FEB.T0_delay = 360;
	FEB.T1_delay = 2000387;
	FEBs.insert({ 218, FEB });
	FEB.HW_mac = 112;
	FEB.SW_mac = 131;
	FEB.SW_modID = 204;
	FEB.T0_delay = 375;
	FEB.T1_delay = 2000402;
	FEBs.insert({ 204, FEB });
	FEB.HW_mac = 62;
	FEB.SW_mac = 117;
	FEB.SW_modID = 190;
	FEB.T0_delay = 391;
	FEB.T1_delay = 2000417;
	FEBs.insert({ 190, FEB });
	FEB.HW_mac = 133;
	FEB.SW_mac = 215;
	FEB.SW_modID = 288;
	FEB.T0_delay = 406;
	FEB.T1_delay = 2000432;
	FEBs.insert({ 288, FEB });
	FEB.HW_mac = 168;
	FEB.SW_mac = 214;
	FEB.SW_modID = 287;
	FEB.T0_delay = 421;
	FEB.T1_delay = 2000448;
	FEBs.insert({ 287, FEB });
	FEB.HW_mac = 182;
	FEB.SW_mac = 116;
	FEB.SW_modID = 189;
	FEB.T0_delay = 436;
	FEB.T1_delay = 2000463;
	FEBs.insert({ 189, FEB });
	FEB.HW_mac = 107;
	FEB.SW_mac = 130;
	FEB.SW_modID = 203;
	FEB.T0_delay = 452;
	FEB.T1_delay = 2000478;
	FEBs.insert({ 203, FEB });
	FEB.HW_mac = 252;
	FEB.SW_mac = 144;
	FEB.SW_modID = 217;
	FEB.T0_delay = 467;
	FEB.T1_delay = 2000494;
	FEBs.insert({ 217, FEB });
	FEB.HW_mac = 141;
	FEB.SW_mac = 143;
	FEB.SW_modID = 216;
	FEB.T0_delay = 482;
	FEB.T1_delay = 2000509;
	FEBs.insert({ 216, FEB });
	FEB.HW_mac = 160;
	FEB.SW_mac = 129;
	FEB.SW_modID = 202;
	FEB.T0_delay = 498;
	// TEMPORARY
	//FEB.T1_delay = 2000424;
	FEB.T1_delay = 2000524;
	FEBs.insert({ 202, FEB });
	FEB.HW_mac = 137;
	FEB.SW_mac = 115;
	FEB.SW_modID = 188;
	FEB.T0_delay = 513;
	FEB.T1_delay = 2000540;
	FEBs.insert({ 188, FEB });
	FEB.HW_mac = 179;
	FEB.SW_mac = 213;
	FEB.SW_modID = 286;
	FEB.T0_delay = 528;
	FEB.T1_delay = 2000555;
	FEBs.insert({ 286, FEB });
	FEB.HW_mac = 66;
	FEB.SW_mac = 186;
	FEB.SW_modID = 259;
	FEB.T0_delay = 298;
	FEB.T1_delay = 2000325;
	FEBs.insert({ 259, FEB });
	FEB.HW_mac = 247;
	FEB.SW_mac = 172;
	FEB.SW_modID = 245;
	FEB.T0_delay = 314;
	FEB.T1_delay = 2000340;
	FEBs.insert({ 245, FEB });
	FEB.HW_mac = 198;
	FEB.SW_mac = 173;
	FEB.SW_modID = 246;
	FEB.T0_delay = 329;
	FEB.T1_delay = 2000356;
	FEBs.insert({ 246, FEB });
	FEB.HW_mac = 243;
	FEB.SW_mac = 174;
	FEB.SW_modID = 247;
	FEB.T0_delay = 344;
	FEB.T1_delay = 2000371;
	FEBs.insert({ 247, FEB });
	FEB.HW_mac = 72;
	FEB.SW_mac = 175;
	FEB.SW_modID = 248;
	FEB.T0_delay = 360;
	FEB.T1_delay = 2000386;
	FEBs.insert({ 248, FEB });
	FEB.HW_mac = 250;
	FEB.SW_mac = 161;
	FEB.SW_modID = 234;
	FEB.T0_delay = 375;
	FEB.T1_delay = 2000401;
	FEBs.insert({ 234, FEB });
	FEB.HW_mac = 249;
	FEB.SW_mac = 147;
	FEB.SW_modID = 220;
	FEB.T0_delay = 390;
	FEB.T1_delay = 2000417;
	FEBs.insert({ 220, FEB });
	FEB.HW_mac = 248;
	FEB.SW_mac = 133;
	FEB.SW_modID = 206;
	FEB.T0_delay = 405;
	FEB.T1_delay = 2000432;
	FEBs.insert({ 206, FEB });
	FEB.HW_mac = 60;
	FEB.SW_mac = 119;
	FEB.SW_modID = 192;
	FEB.T0_delay = 421;
	FEB.T1_delay = 2000447;
	FEBs.insert({ 192, FEB });
	FEB.HW_mac = 145;
	FEB.SW_mac = 217;
	FEB.SW_modID = 290;
	FEB.T0_delay = 436;
	FEB.T1_delay = 2000463;
	FEBs.insert({ 290, FEB });
	FEB.HW_mac = 110;
	FEB.SW_mac = 216;
	FEB.SW_modID = 289;
	FEB.T0_delay = 451;
	FEB.T1_delay = 2000478;
	FEBs.insert({ 289, FEB });
	FEB.HW_mac = 59;
	FEB.SW_mac = 118;
	FEB.SW_modID = 191;
	FEB.T0_delay = 467;
	FEB.T1_delay = 2000493;
	FEBs.insert({ 191, FEB });
	FEB.HW_mac = 202;
	FEB.SW_mac = 132;
	FEB.SW_modID = 205;
	FEB.T0_delay = 482;
	FEB.T1_delay = 2000509;
	FEBs.insert({ 205, FEB });
	FEB.HW_mac = 135;
	FEB.SW_mac = 146;
	FEB.SW_modID = 219;
	FEB.T0_delay = 497;
	FEB.T1_delay = 2000524;
	FEBs.insert({ 219, FEB });
	FEB.HW_mac = 246;
	FEB.SW_mac = 160;
	FEB.SW_modID = 233;
	FEB.T0_delay = 513;
	FEB.T1_delay = 2000539;
	FEBs.insert({ 233, FEB });
	FEB.HW_mac = 253;
	FEB.SW_mac = 187;
	FEB.SW_modID = 260;
	FEB.T0_delay = 342;
	FEB.T1_delay = 2000369;
	FEBs.insert({ 260, FEB });
	FEB.HW_mac = 245;
	FEB.SW_mac = 188;
	FEB.SW_modID = 261;
	FEB.T0_delay = 358;
	FEB.T1_delay = 2000384;
	FEBs.insert({ 261, FEB });
	FEB.HW_mac = 65;
	FEB.SW_mac = 189;
	FEB.SW_modID = 262;
	FEB.T0_delay = 373;
	FEB.T1_delay = 2000400;
	FEBs.insert({ 262, FEB });
	FEB.HW_mac = 57;
	FEB.SW_mac = 190;
	FEB.SW_modID = 263;
	FEB.T0_delay = 388;
	FEB.T1_delay = 2000415;
	FEBs.insert({ 263, FEB });
	FEB.HW_mac = 63;
	FEB.SW_mac = 176;
	FEB.SW_modID = 249;
	FEB.T0_delay = 404;
	FEB.T1_delay = 2000430;
	FEBs.insert({ 249, FEB });
	FEB.HW_mac = 251;
	FEB.SW_mac = 177;
	FEB.SW_modID = 250;
	FEB.T0_delay = 419;
	FEB.T1_delay = 2000445;
	FEBs.insert({ 250, FEB });
	FEB.HW_mac = 70;
	FEB.SW_mac = 163;
	FEB.SW_modID = 236;
	FEB.T0_delay = 434;
	FEB.T1_delay = 2000461;
	FEBs.insert({ 236, FEB });
	FEB.HW_mac = 155;
	FEB.SW_mac = 149;
	FEB.SW_modID = 222;
	FEB.T0_delay = 449;
	FEB.T1_delay = 2000476;
	FEBs.insert({ 222, FEB });
	FEB.HW_mac = 154;
	FEB.SW_mac = 135;
	FEB.SW_modID = 208;
	FEB.T0_delay = 465;
	FEB.T1_delay = 2000491;
	FEBs.insert({ 208, FEB });
	FEB.HW_mac = 85;
	FEB.SW_mac = 121;
	FEB.SW_modID = 194;
	FEB.T0_delay = 480;
	FEB.T1_delay = 2000507;
	FEBs.insert({ 194, FEB });
	FEB.HW_mac = 134;
	FEB.SW_mac = 219;
	FEB.SW_modID = 292;
	FEB.T0_delay = 495;
	FEB.T1_delay = 2000522;
	FEBs.insert({ 292, FEB });
	FEB.HW_mac = 129;
	FEB.SW_mac = 218;
	FEB.SW_modID = 291;
	FEB.T0_delay = 511;
	FEB.T1_delay = 2000537;
	FEBs.insert({ 291, FEB });
	FEB.HW_mac = 115;
	FEB.SW_mac = 120;
	FEB.SW_modID = 193;
	FEB.T0_delay = 526;
	FEB.T1_delay = 2000553;
	FEBs.insert({ 193, FEB });
	FEB.HW_mac = 204;
	FEB.SW_mac = 134;
	FEB.SW_modID = 207;
	FEB.T0_delay = 541;
	FEB.T1_delay = 2000568;
	FEBs.insert({ 207, FEB });
	FEB.HW_mac = 244;
	FEB.SW_mac = 148;
	FEB.SW_modID = 221;
	FEB.T0_delay = 557;
	FEB.T1_delay = 2000583;
	FEBs.insert({ 221, FEB });
	FEB.HW_mac = 82;
	FEB.SW_mac = 162;
	FEB.SW_modID = 235;
	FEB.T0_delay = 572;
	FEB.T1_delay = 2000598;
	FEBs.insert({ 235, FEB });
	FEB.HW_mac = 186;
	FEB.SW_mac = 199;
	FEB.SW_modID = 272;
	FEB.T0_delay = 284;
	FEB.T1_delay = 2000310;
	FEBs.insert({ 272, FEB });
	FEB.HW_mac = 83;
	FEB.SW_mac = 200;
	FEB.SW_modID = 273;
	FEB.T0_delay = 299;
	FEB.T1_delay = 2000326;
	FEBs.insert({ 273, FEB });
	FEB.HW_mac = 254;
	FEB.SW_mac = 201;
	FEB.SW_modID = 274;
	FEB.T0_delay = 314;
	FEB.T1_delay = 2000341;
	FEBs.insert({ 274, FEB });
	FEB.HW_mac = 166;
	FEB.SW_mac = 202;
	FEB.SW_modID = 275;
	FEB.T0_delay = 330;
	FEB.T1_delay = 2000356;
	FEBs.insert({ 275, FEB });
	FEB.HW_mac = 178;
	FEB.SW_mac = 203;
	FEB.SW_modID = 276;
	FEB.T0_delay = 345;
	FEB.T1_delay = 2000371;
	FEBs.insert({ 276, FEB });
	FEB.HW_mac = 136;
	FEB.SW_mac = 204;
	FEB.SW_modID = 277;
	FEB.T0_delay = 360;
	FEB.T1_delay = 2000387;
	FEBs.insert({ 277, FEB });
	FEB.HW_mac = 184;
	FEB.SW_mac = 205;
	FEB.SW_modID = 278;
	FEB.T0_delay = 375;
	FEB.T1_delay = 2000402;
	FEBs.insert({ 278, FEB });
	FEB.HW_mac = 187;
	FEB.SW_mac = 191;
	FEB.SW_modID = 264;
	FEB.T0_delay = 391;
	FEB.T1_delay = 2000417;
	FEBs.insert({ 264, FEB });
	FEB.HW_mac = 240;
	FEB.SW_mac = 231;
	FEB.SW_modID = 304;
	FEB.T0_delay = 406;
	FEB.T1_delay = 2000433;
	FEBs.insert({ 304, FEB });
	FEB.HW_mac = 242;
	FEB.SW_mac = 230;
	FEB.SW_modID = 303;
	FEB.T0_delay = 421;
	FEB.T1_delay = 2000448;
	FEBs.insert({ 303, FEB });
	FEB.HW_mac = 188;
	FEB.SW_mac = 229;
	FEB.SW_modID = 302;
	FEB.T0_delay = 437;
	FEB.T1_delay = 2000463;
	FEBs.insert({ 302, FEB });
	FEB.HW_mac = 58;
	FEB.SW_mac = 228;
	FEB.SW_modID = 301;
	FEB.T0_delay = 452;
	FEB.T1_delay = 2000479;
	FEBs.insert({ 301, FEB });
	FEB.HW_mac = 143;
	FEB.SW_mac = 227;
	FEB.SW_modID = 300;
	FEB.T0_delay = 467;
	FEB.T1_delay = 2000494;
	FEBs.insert({ 300, FEB });
	FEB.HW_mac = 235;
	FEB.SW_mac = 226;
	FEB.SW_modID = 299;
	FEB.T0_delay = 483;
	FEB.T1_delay = 2000509;
	FEBs.insert({ 299, FEB });

	return FEBs;
}

#endif
