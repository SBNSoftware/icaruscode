#ifndef CRT_TIMING_H
#define CRT_TIMING_H

#include <TTree.h>

#include <vector>
#include <utility>
#include <map>

//#include "./CRTRawTree.h"
//#include "./CRTPreProcessTree.h"

#include "icaruscode/CRT/CRTDecoder/CRTRawTree.h"
#include "icaruscode/CRT/CRTDecoder/CRTPreProcessTree.h"

using std::map;
using std::vector;
using std::pair;
using std::make_pair;

namespace icarus {
  namespace crt {
	class CRTTiming;
  }
}

class icarus::crt::CRTTiming {

  public:
	explicit CRTTiming(CRTPreProcessTree &raw);
	explicit CRTTiming(CRTRawTree &raw);
	void TimeOrder();
	const map<size_t,size_t>* GetRawToOrderedMap();
	const map<size_t,size_t>* GetOrderedToRawMap();
	void DumpSortedTimes(size_t nmax);
	void DumpRawTimes(size_t nmax);

  private:
	const CRTRawTree* fRaw;
	const CRTPreProcessTree* fPre;
	map<size_t,size_t> fRawToOrdered;	
	map<size_t,size_t> fOrderedToRaw;
	bool fHasSort;
	char fType;
	size_t fNEntries;
};

#endif
