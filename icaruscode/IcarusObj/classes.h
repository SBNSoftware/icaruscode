#include "canvas/Persistency/Common/Wrapper.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/Assns.h"
#include "icaruscode/IcarusObj/SimEnergyDepositSummary.h"
#include "icaruscode/IcarusObj/CRTTPCMatchingInfo.h"
#include "icaruscode/IcarusObj/OpDetWaveformMeta.h"
#include "sbnobj/ICARUS/PMT/Trigger/Data/OpticalTriggerGate.h"
#include "sbnobj/ICARUS/PMT/Trigger/Data/TriggerGateData.h"
<<<<<<< HEAD
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/Simulation/BeamGateInfo.h"
=======
#include "sbnobj/Common/CRT/CRTHit.hh"
>>>>>>> SBNSoftware/feature/gp_CRTmatchMetadata
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/AnalysisBase/T0.h"
#include <vector>

namespace {
  icarus::SimEnergyDepositSummary EDepSum;
  icarus::CRTTPCMatchingInfo CRTTPCMatchingInfo;
}
