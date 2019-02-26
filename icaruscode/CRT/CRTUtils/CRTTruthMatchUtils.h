#ifndef CRTTRUTHMATCHUTILS_H_SEEN
#define CRTTRUTHMATCHUTILS_H_SEEN

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

// LArSoft
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "icaruscode/CRT/CRTProducts/CRTData.hh"
#include "icaruscode/CRT/CRTProducts/CRTChannelData.h"
#include "icaruscode/CRT/CRTProducts/CRTHit.hh"

// c++
#include <iostream>
#include <stdio.h>
#include <sstream>
#include <vector>
#include <map>
#include <utility>
#include <cmath> 
#include <memory>

namespace icarus {
namespace CRTTruthMatchUtils{

  std::vector<int> AllTrueIds(art::Handle<std::vector<crt::CRTHit>> hitHandle, 
                              const art::Event& event, art::InputTag hitLabel, int hit_i);

  std::vector<int> AllTrueIds(art::Handle<std::vector<crt::CRTData>> dataHandle, 
                              const art::Event& event, art::InputTag dataLabel, int data_i);

  std::vector<int> AllTrueIds(art::Handle<std::vector<crt::CRTHit>> hitHandle, 
                              const art::Event& event, art::InputTag hitLabel, art::InputTag dataLabel, int hit_i);

  int TrueIdFromTotalEnergy(art::Handle<std::vector<crt::CRTData>> dataHandle, 
                            const art::Event& event, art::InputTag dataLabel, int data_i);

  int TrueIdFromTotalEnergy(art::Handle<std::vector<crt::CRTHit>> hitHandle, 
                            const art::Event& event, art::InputTag hitLabel, art::InputTag dataLabel, int hit_i);

}//namespace CRTTruthMatchUtils
}//namespace icarus

#endif
