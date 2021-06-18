#ifndef CRT_TIMING_CC
#define CRT_TIMING_CC

//#include "./CRTTiming.h"
#include "icaruscode/CRT/CRTDecoder/CRTTiming.h"


using namespace icarus::crt;

CRTTiming::CRTTiming(CRTPreProcessTree& pre) {
        fPre = &pre;
	fRaw = 0;
        fHasSort = false;
	fType='p';
	fNEntries = fPre->GetNEntries();
}

CRTTiming::CRTTiming(CRTRawTree& raw) {
	fRaw = &raw;
	fPre = 0;
	fHasSort = false;
	fType = 'r';
	fNEntries = fRaw->GetNEntries();
}

void CRTTiming::TimeOrder() {

	vector<pair<size_t,uint64_t>> pairs;

	for(size_t ientry=0; ientry<fNEntries; ientry++) {
		uint64_t t0=0;
		if(fType=='p') t0 = fPre->GetAbsTime(ientry);
		if(fType=='r') t0 = fRaw->GetAbsTime(ientry);
		pairs.push_back(make_pair(ientry,t0));
	}

	std::sort(pairs.begin(),pairs.end(),[](pair<size_t,uint64_t> a, pair<size_t,uint64_t> b) {
        	return a.second < b.second;   
        });

	for(size_t ientry=0; ientry<fNEntries; ientry++) {
		fRawToOrdered[pairs[ientry].first] = ientry;
		fOrderedToRaw[ientry] = pairs[ientry].first;
	}

	fHasSort = true;

	return;
}

const map<size_t,size_t>* CRTTiming::GetRawToOrderedMap() {

	if(!fHasSort) TimeOrder();

	return &fRawToOrdered;
}

const map<size_t,size_t>* CRTTiming::GetOrderedToRawMap() {

        if(!fHasSort) TimeOrder();

        return &fOrderedToRaw;
}

void CRTTiming::DumpSortedTimes(size_t nmax=0) {
	for(size_t i=0; i<fNEntries; i++) {
		if(nmax!=0 && i==nmax) return;

		if(fType=='p') std::cout << fPre->GetAbsTime(fOrderedToRaw[i]) << std::endl;
		if(fType=='r') std::cout << fRaw->GetAbsTime(fOrderedToRaw[i]) << std::endl;

	}
}

void CRTTiming::DumpRawTimes(size_t nmax=0) {
        for(size_t i=0; i<fNEntries; i++) {
                if(nmax!=0 && i==nmax) return;

		if(fType=='p') std::cout << fPre->GetAbsTime(i) << std::endl;
                if(fType=='r') std::cout << fRaw->GetAbsTime(i) << std::endl;

        }
}
#endif
