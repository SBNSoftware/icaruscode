
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
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "canvas/Persistency/Common/FindManyP.h"
#include "larcorealg/CoreUtils/NumericUtils.h"

// LArSoft
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "larcore/Geometry/AuxDetGeometry.h"

// Utility libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"
#include "cetlib/pow.h" // cet::sum_of_squares()
#include "icaruscode/CRT/CRTProducts/CRTHit.hh"
#include "icaruscode/CRT/CRTProducts/CRTData.hh"
#include "icaruscode/CRT/CRTProducts/CRTChannelData.h"

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

namespace icarus {

  class CRTHitRecoAlg {
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

    };

    CRTHitRecoAlg(const Config& config);

    CRTHitRecoAlg(const fhicl::ParameterSet& pset) :
      CRTHitRecoAlg(fhicl::Table<Config>(pset, {})()) {}

    CRTHitRecoAlg();

    ~CRTHitRecoAlg();

    void reconfigure(const Config& config);

    char MacToType(int mac);
    int MacToRegion(int mac);
    std::string MacToRegionName(int mac);
    int MacToAuxDetID(int mac, int chan);
    int ChannelToAuxDetSensitiveID(int mac, int chan);

    std::vector<std::pair<crt::CRTHit, std::vector<int>>> CreateCRTHits(std::vector<art::Ptr<crt::CRTData>> crtList);

    // Function to make filling a CRTHit a bit faster
    icarus::crt::CRTHit FillCrtHit(std::vector<uint8_t> tfeb_id, std::map<uint8_t, 
                           std::vector<std::pair<int,float>>> tpesmap, float peshit, double time0, double time1, int plane, 
                           double x, double ex, double y, double ey, double z, double ez, std::string tagger); 
  private:

    geo::GeometryCore const* fGeometryService;
    detinfo::DetectorClocks const* fDetectorClocks;
    detinfo::DetectorProperties const* fDetectorProperties;
    detinfo::ElecClock fTrigClock;
    //art::ServiceHandle<geo::AuxDetGeometry> fAuxDetGeoService;
    //const geo::AuxDetGeometry* fAuxDetGeo;
    //const geo::AuxDetGeometryCore* fAuxDetGeoCore;

    //Params from fcl file
    bool fVerbose;          ///< print info
    bool fUseReadoutWindow; ///< Only reconstruct hits within readout window
    double fQPed;           ///< Pedestal offset of SiPMs [ADC]
    double fQSlope;         ///< Pedestal slope of SiPMs [ADC/photon]
    double fPropDelay;      ///< propegation time [ns/cm]

  };

}

#endif
