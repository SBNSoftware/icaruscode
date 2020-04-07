#ifndef CRTHITRECOALG_H_SEEN
#define CRTHITRECOALG_H_SEEN
///////////////////////////////////////////////
// CRTHitRecoAlg.h
//
// Functions for CRT hit reconstruction
// T Brooks (tbrooks@fnal.gov), November 2018
///////////////////////////////////////////////

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
#include "icaruscode/CRT/CRTProducts/CRTHit.hh"
#include "icaruscode/CRT/CRTProducts/CRTData.hh"
#include "icaruscode/CRT/CRTUtils/CRTCommonUtils.h"

// c++
#include <iostream>
#include <stdio.h>
#include <sstream>
#include <vector>
#include <map>
#include <utility>
#include <cmath>
#include <memory>

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
    fhicl::Atom<double> PropDelay {
      Name("PropDelay"),
      Comment("group velocity in WLS fiber [ns/cm]")
    };
    fhicl::Atom<double> PEThresh {
      Name("PEThresh"),
      Comment("threshold in photoelectrons above which charge amplitudes used in hit reco")
    };
    fhicl::Atom<uint32_t> CoinWindow {
      Name("CoinWindow"),
      Comment("window for finding side CRT trigger coincidences [ns]")
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
  // Function to make filling a CRTHit a bit faster
  CRTHit FillCRTHit(vector<uint8_t> tfeb_id, map<uint8_t, vector<pair<int,float>>> tpesmap,
                    float peshit, double time0, double time1, int plane,
                    double x, double ex, double y, double ey, double z, double ez, string tagger);


 private:

  geo::GeometryCore const* fGeometryService;

  //Params from fcl file
  bool fVerbose;          ///< print info
  bool fUseReadoutWindow; ///< Only reconstruct hits within TPC readout window
  double fQPed;           ///< Pedestal offset of SiPMs [ADC]
  double fQSlope;         ///< Pedestal slope of SiPMs [ADC/photon]
  double fPropDelay;      ///< propegation time [ns/cm]
  double fPEThresh;       ///< threshold[PE] above which charge amplitudes used in hit reco
  uint32_t fCoinWindow;   ///< Coincidence window used for grouping side CRT triggers [ns]

  //Given top CRTData product, produce CRTHit
  CRTHit MakeTopHit(art::Ptr<CRTData> data);
  //Given bottom CRTData product, produce CRTHit
  CRTHit MakeBottomHit(art::Ptr<CRTData> data);
  //Given vector of side CRTData products, produce CRTHit
  CRTHit MakeSideHit(vector<art::Ptr<CRTData>> coinData);
  // Check if a hit is empty
  bool IsEmptyHit(CRTHit hit);

}; //class CRTHitRecoAlg

#endif
