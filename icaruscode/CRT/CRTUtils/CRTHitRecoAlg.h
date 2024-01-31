#ifndef CRTHITRECOALG_H_SEEN
#define CRTHITRECOALG_H_SEEN
////////////////////////////////////////////////////////////////////
// CRTHitRecoAlg.h
//
// Functions for CRT hit reconstruction
// Chris Hilgenberg (Chris.Hilgenberg@colostate.edu), June 2020
////////////////////////////////////////////////////////////////////

// framework
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "fhiclcpp/ParameterSet.h"
#include "larcorealg/CoreUtils/NumericUtils.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft
#include "larcore/Geometry/AuxDetGeometry.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"

// Utility libraries
#include "cetlib/pow.h"  // cet::sum_of_squares()
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// icaruscode includes
#include "icaruscode/CRT/CRTUtils/CRTCommonUtils.h"
#include "sbnobj/Common/CRT/CRTHit.hh"
#include "sbnobj/ICARUS/CRT/CRTData.hh"

#include "icaruscode/Decode/ChannelMapping/ICARUSChannelMapProvider.h"
#include "icaruscode/Decode/ChannelMapping/IChannelMapping.h"
#include "icaruscode/Decode/ChannelMapping/IICARUSChannelMap.h"
// c++
#include <stdio.h>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <utility>
#include <vector>
// ROOT
#include "TGeoManager.h"
#include "TVector3.h"

using std::map;
using std::pair;
using std::string;
using std::vector;

namespace icarus::crt {
class CRTHitRecoAlg;
inline bool sortbytime(const pair<int, ULong64_t>& a,
                       const pair<int, ULong64_t>& b) {
  return (a.second < b.second);
}
}  // namespace icarus::crt

class icarus::crt::CRTHitRecoAlg {
 public:
  using CRTData = icarus::crt::CRTData;
  using CRTHit = sbn::crt::CRTHit;
  
  explicit CRTHitRecoAlg(const fhicl::ParameterSet& pset);
  CRTHitRecoAlg();
  void reconfigure(const fhicl::ParameterSet& pset);

  // produce CRTHits with associated data indices from input vector of CRTData
  vector<pair<CRTHit, vector<int>>> CreateCRTHits(
      vector<art::Ptr<CRTData>> crtList, uint64_t trigger_timestamp);

  // preselection based on charge in a CRTData
  vector<art::Ptr<CRTData>> PreselectCRTData(
      const vector<art::Ptr<CRTData>>& crtList, uint64_t trigger_timestamp);

  // Function to make filling a CRTHit a bit faster
  CRTHit FillCRTHit(vector<uint8_t> tfeb_id,
                    map<uint8_t, vector<pair<int, float>>> tpesmap,
                    float peshit, uint64_t time0, Long64_t time1, int plane,
                    double x, double ex, double y, double ey, double z,
                    double ez, string tagger);

 private:
  geo::GeometryCore const* fGeometryService;

  CRTCommonUtils fCrtutils;

  // Params from fcl file
  bool fVerbose;           ///< print info
  bool foutCSVFile;        ///< FCL input: Write a CSV File?
  std::string fCSVFile;    ///< FCL input: file name for output CSV File
  bool fUseReadoutWindow;  ///< Only reconstruct hits within TPC readout window
  double fQPed;            ///< Pedestal offset of SiPMs [ADC]
  double fQSlope;          ///< Pedestal slope of SiPMs [ADC/photon]
  double fPropDelay;       ///< propegation time [ns/cm]
  double fPEThresh;  ///< threshold[PE] above which charge amplitudes used in
                     ///< hit reco
  double ftopGain;   ///< Dummy Top CRT Gain Value
  double ftopPed;    ///< Dummy Top CRT Pedestal Value
  uint64_t fSiPMtoFEBdelay;  ///< SiPM to FEB cable induced delay: 11.6 [ns]
  uint64_t fCoinWindow;      ///< Coincidence window used for grouping side CRT
                             ///< triggers [ns]
  uint64_t fCrtWindow;  ///< Looking data window within trigger timestamp [ns]
  std::ofstream filecsv;
  bool fData;  ///< look for only data
  double fGlobalT0Offset;  ///< Offset to be applied to MC CRT hit times.
  // See DocDB 34763 for CRT Time distributions before and after this 
  // fGlobalT0Offset parameter was moved to the  CRTHitRecoAlg.

  const icarusDB::IICARUSChannelMap* fChannelMap = nullptr;

  // Given top CRTData product, produce CRTHit
  CRTHit MakeTopHit(art::Ptr<CRTData> data, ULong64_t GlobalTrigger[]);
  // Given bottom CRTData product, produce CRTHit
  CRTHit MakeBottomHit(art::Ptr<CRTData> data);
  // Given vector of side CRTData products, produce CRTHit
  CRTHit MakeSideHit(vector<art::Ptr<CRTData>> coinData,
                     ULong64_t GlobalTrigger[]);
  // Check if a hit is empty
  bool IsEmptyHit(CRTHit hit);
  // function to appply appropriate prop delay for Side full vs cut modules
  // (North and South walls are cut modules)
  int64_t RegionDelay(std::string const& region) const;

  std::map<uint8_t, int32_t> FEB_T1delay_side;  //<mac5, delay in ns>
  std::map<uint8_t, int32_t> FEB_T0delay_side;  //<mac5, delay in ns>
  // Following: Data driven timing correction of the light propagation
  // to the SiPM readout sector-by-sector for the Top CRT.
  // These numbers refer to DocDB 29405 slide number 9
  double const TopCRT_TimingCorr[64] = {  // [ns]
      8.2, 7.0,  5.8,  4.5, 3.2,  1.9,  1.1,  0.7,  8.3,  7.3,  5.8,  4.7,  3.4,
      2.3, 1.5,  1.2,  9.0, 8.0,  6.6,  5.6,  4.6,  3.5,  3.4,  3.1,  10.0, 8.9,
      7.4, 6.6,  5.9,  5.0, 5.4,  4.7,  10.7, 9.9,  8.4,  7.7,  7.2,  6.6,  6.7,
      6.4, 11.7, 10.6, 9.2, 8.5,  8.4,  7.8,  8.1,  7.7,  12.5, 11.6, 10.2, 9.5,
      9.2, 8.9,  9.0,  8.8, 13.0, 12.3, 10.9, 10.5, 10.4, 9.8,  10.0, 9.7};
  static bool compareBytime(art::Ptr<CRTData> const& a,
                            art::Ptr<CRTData> const& b) {
    return a->fTs0 < b->fTs0;
  }
};  // class CRTHitRecoAlg

namespace icarus::crt {
ULong64_t GetMode(std::vector<std::pair<int, ULong64_t>> vector);

struct FEB_delay {
  int HW_mac = -1;
  int SW_mac = -1;
  int SW_modID = -1;
  ULong64_t T0_delay = 0;  //[ns]
  ULong64_t T1_delay = 0;  //[ns]
};

typedef int feb_index;
typedef std::map<feb_index, FEB_delay> CRT_delay_map;

CRT_delay_map LoadFEBMap();
}  // namespace icarus::crt

inline icarus::crt::CRT_delay_map icarus::crt::LoadFEBMap() {
  CRT_delay_map FEBs;

  FEBs = {{271, {81, 198, 271, 283ull, 2000309ull}},
          {270, {119, 197, 270, 298ull, 2000324ull}},
          {269, {87, 196, 269, 313ull, 2000339ull}},
          {268, {92, 195, 268, 329ull, 2000355ull}},
          {267, {180, 194, 267, 344ull, 2000370ull}},
          {266, {97, 193, 266, 359ull, 2000385ull}},
          {265, {174, 192, 265, 374ull, 2000400ull}},
          {251, {238, 178, 251, 390ull, 2000416ull}},
          {237, {234, 164, 237, 405ull, 2000431ull}},
          {297, {189, 224, 297, 420ull, 2000446ull}},
          {296, {190, 223, 296, 436ull, 2000462ull}},
          {295, {80, 222, 295, 451ull, 2000477ull}},
          {294, {162, 221, 294, 466ull, 2000492ull}},
          {293, {64, 220, 293, 482ull, 2000508ull}},
          {255, {172, 182, 255, 298ull, 2000324ull}},
          {254, {114, 181, 254, 313ull, 2000339ull}},
          {253, {100, 180, 253, 328ull, 2000355ull}},
          {252, {150, 179, 252, 344ull, 2000370ull}},
          {238, {176, 165, 238, 359ull, 2000385ull}},
          {224, {67, 151, 224, 374ull, 2000400ull}},
          {223, {138, 150, 223, 390ull, 2000416ull}},
          {209, {170, 136, 209, 405ull, 2000431ull}},
          {195, {101, 122, 195, 420ull, 2000446ull}},
          {181, {142, 108, 181, 435ull, 2000462ull}},
          {279, {139, 206, 279, 451ull, 2000477ull}},
          {280, {185, 207, 280, 466ull, 2000492ull}},
          {182, {6, 109, 182, 481ull, 2000508ull}},
          {196, {177, 123, 196, 497ull, 2000523ull}},
          {210, {61, 137, 210, 512ull, 2000538ull}},
          {256, {123, 183, 256, 298ull, 2000325ull}},
          {242, {116, 169, 242, 314ull, 2000340ull}},
          {241, {104, 168, 241, 329ull, 2000355ull}},
          {240, {91, 167, 240, 344ull, 2000371ull}},
          {239, {88, 166, 239, 360ull, 2000386ull}},
          {225, {120, 152, 225, 375ull, 2000401ull}},
          {211, {132, 138, 211, 390ull, 2000417ull}},
          {197, {95, 124, 197, 405ull, 2000432ull}},
          {183, {232, 110, 183, 421ull, 2000447ull}},
          {281, {165, 208, 281, 436ull, 2000463ull}},
          {282, {148, 209, 282, 451ull, 2000478ull}},
          {184, {237, 111, 184, 467ull, 2000493ull}},
          {198, {102, 125, 198, 482ull, 2000508ull}},
          {212, {94, 139, 212, 497ull, 2000524ull}},
          {226, {130, 153, 226, 513ull, 2000539ull}},
          {257, {181, 184, 257, 284ull, 2000310ull}},
          {243, {124, 170, 243, 299ull, 2000325ull}},
          {229, {152, 156, 229, 314ull, 2000341ull}},
          {228, {98, 155, 228, 329ull, 2000356ull}},
          {227, {173, 154, 227, 345ull, 2000371ull}},
          {213, {169, 140, 213, 360ull, 2000387ull}},
          {199, {144, 126, 199, 375ull, 2000402ull}},
          {185, {239, 112, 185, 391ull, 2000417ull}},
          {283, {147, 210, 283, 406ull, 2000433ull}},
          {284, {105, 211, 284, 421ull, 2000448ull}},
          {186, {231, 114, 186, 437ull, 2000463ull}},
          {200, {117, 127, 200, 452ull, 2000478ull}},
          {214, {126, 141, 214, 467ull, 2000494ull}},
          {215, {90, 142, 215, 482ull, 2000509ull}},
          {201, {183, 128, 201, 498ull, 2000524ull}},
          {187, {241, 114, 187, 513ull, 2000540ull}},
          {285, {113, 212, 285, 528ull, 2000555ull}},
          {258, {233, 185, 258, 283ull, 2000310ull}},
          {244, {164, 171, 244, 299ull, 2000325ull}},
          {230, {161, 157, 230, 314ull, 2000341ull}},
          {231, {203, 158, 231, 329ull, 2000356ull}},
          {232, {122, 159, 232, 345ull, 2000371ull}},
          {218, {2, 145, 218, 360ull, 2000387ull}},
          {204, {112, 131, 204, 375ull, 2000402ull}},
          {190, {62, 117, 190, 391ull, 2000417ull}},
          {288, {133, 215, 288, 406ull, 2000432ull}},
          {287, {168, 214, 287, 421ull, 2000448ull}},
          {189, {182, 116, 189, 436ull, 2000463ull}},
          {203, {107, 130, 203, 452ull, 2000478ull}},
          {217, {252, 144, 217, 467ull, 2000494ull}},
          {216, {141, 143, 216, 482ull, 2000509ull}},
          {202, {160, 129, 202, 498ull, 2000524ull}},
          {188, {137, 115, 188, 513ull, 2000540ull}},
          {286, {179, 213, 286, 528ull, 2000555ull}},
          {259, {66, 186, 259, 298ull, 2000325ull}},
          {245, {247, 172, 245, 314ull, 2000340ull}},
          {246, {198, 173, 246, 329ull, 2000356ull}},
          {247, {243, 174, 247, 344ull, 2000371ull}},
          {248, {72, 175, 248, 360ull, 2000386ull}},
          {234, {250, 161, 234, 375ull, 2000401ull}},
          {220, {249, 147, 220, 390ull, 2000417ull}},
          {206, {248, 133, 206, 405ull, 2000432ull}},
          {192, {60, 119, 192, 421ull, 2000447ull}},
          {290, {145, 217, 290, 436ull, 2000463ull}},
          {289, {110, 216, 289, 451ull, 2000478ull}},
          {191, {59, 118, 191, 467ull, 2000493ull}},
          {205, {202, 132, 205, 482ull, 2000509ull}},
          {219, {135, 146, 219, 497ull, 2000524ull}},
          {233, {246, 160, 233, 513ull, 2000539ull}},
          {260, {253, 187, 260, 342ull, 2000369ull}},
          {261, {245, 188, 261, 358ull, 2000384ull}},
          {262, {65, 189, 262, 373ull, 2000400ull}},
          {263, {57, 190, 263, 388ull, 2000415ull}},
          {249, {63, 176, 249, 404ull, 2000430ull}},
          {250, {251, 177, 250, 419ull, 2000445ull}},
          {236, {70, 163, 236, 434ull, 2000461ull}},
          {222, {155, 149, 222, 449ull, 2000476ull}},
          {208, {154, 135, 208, 465ull, 2000491ull}},
          {194, {85, 121, 194, 480ull, 2000507ull}},
          {292, {134, 219, 292, 495ull, 2000522ull}},
          {291, {129, 218, 291, 511ull, 2000537ull}},
          {193, {115, 120, 193, 526ull, 2000553ull}},
          {207, {204, 134, 207, 541ull, 2000568ull}},
          {221, {244, 148, 221, 557ull, 2000583ull}},
          {235, {82, 162, 235, 572ull, 2000598ull}},
          {272, {186, 199, 272, 284ull, 2000310ull}},
          {273, {83, 200, 273, 299ull, 2000326ull}},
          {274, {254, 201, 274, 314ull, 2000341ull}},
          {275, {166, 202, 275, 330ull, 2000356ull}},
          {276, {178, 203, 276, 345ull, 2000371ull}},
          {277, {136, 204, 277, 360ull, 2000387ull}},
          {278, {184, 205, 278, 375ull, 2000402ull}},
          {264, {187, 191, 264, 391ull, 2000417ull}},
          {304, {240, 231, 304, 406ull, 2000433ull}},
          {303, {242, 230, 303, 421ull, 2000448ull}},
          {302, {188, 229, 302, 437ull, 2000463ull}},
          {301, {58, 228, 301, 452ull, 2000479ull}},
          {300, {143, 227, 300, 467ull, 2000494ull}},
          {299, {235, 226, 299, 483ull, 2000509ull}}};
  return FEBs;
}

#endif
