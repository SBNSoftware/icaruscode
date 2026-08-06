#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Persistency/Common/Wrapper.h"

#include "icaruscode/TPC/Tracking/MCS/MCSFitResultGS.h"

namespace recob {
  class MCSFitResultGS;
}

template class art::Wrapper<std::vector<recob::MCSFitResultGS>>;