#include "canvas/Persistency/Common/Wrapper.h"
#include "canvas/Persistency/Common/Assns.h"
#include "icaruscode/IcarusObj/SimEnergyDepositSummary.h"
#include "icaruscode/IcarusObj/CRTTPCMatchingInfo.h"
#include "icaruscode/IcarusObj/OpDetWaveformMeta.h"
#include "sbnobj/ICARUS/PMT/Trigger/Data/OpticalTriggerGate.h"
#include "sbnobj/ICARUS/PMT/Trigger/Data/TriggerGateData.h"
#include "sbnobj/Common/CRT/CRTHit.hh"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/AnalysisBase/T0.h"
#include <vector>

namespace {
  icarus::SimEnergyDepositSummary EDepSum;
  icarus::CRTTPCMatchingInfo CRTTPCMatchingInfo;
}
