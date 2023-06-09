#include "canvas/Persistency/Common/Wrapper.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/Assns.h"

#include "icaruscode/IcarusObj/SimEnergyDepositSummary.h"
#include "icaruscode/IcarusObj/CRTTPCMatchingInfo.h"
#include "icaruscode/IcarusObj/OpDetWaveformMeta.h"
#include "icaruscode/IcarusObj/PMTWaveformTimeCorrection.h"
#include "icaruscode/IcarusObj/ChannelROI.h"
#include "icaruscode/IcarusObj/CRTPMTMatching.h"

#include "sbnobj/ICARUS/PMT/Trigger/Data/OpticalTriggerGate.h"
#include "sbnobj/ICARUS/PMT/Trigger/Data/TriggerGateData.h"

#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/Simulation/BeamGateInfo.h"
#include "sbnobj/Common/CRT/CRTHit.hh"
#include "lardataobj/RawData/OpDetWaveform.h"
#include <vector>

namespace {
  icarus::SimEnergyDepositSummary EDepSum;
  
  icarus::crt::CRTPMTMatchingInfo meta1;
  
  art::Assns<sbn::crt::CRTHit, recob::OpFlash, icarus::crt::CRTPMTMatchingInfo> assn11;
  art::Assns<recob::OpFlash, sbn::crt::CRTHit, icarus::crt::CRTPMTMatchingInfo> assn12;
  
  art::Assns<recob::OpFlash, icarus::crt::CRTPMTMatching> assn21;
  art::Assns<icarus::crt::CRTPMTMatching, recob::OpFlash> assn22;
  
  art::Assns<sbn::crt::CRTHit, icarus::crt::CRTPMTMatching> assn31;
  art::Assns<icarus::crt::CRTPMTMatching, sbn::crt::CRTHit> assn32;
  
}
