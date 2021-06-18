#ifndef CRT_RAW_TREE_H
#define CRT_RAW_TREE_H

#include <TTree.h>
#include <TBranch.h>

#include <iostream>
#include <cfloat>

namespace icarus{
 namespace crt {
	class CRTRawTree;
 }
}


class icarus::crt::CRTRawTree{

 public:
	explicit CRTRawTree(TTree* tr);
	size_t   GetNEntries() const;
	uint64_t GetAbsTime(size_t ientry) const;
	uint8_t  GetMac(size_t ientry) const;
	uint16_t GetADC(size_t ientry, uint8_t chan) const;
	float    GetInstRate(size_t ientry_prev, size_t ientry_next) const;
	float    GetPollRate(size_t ientry) const;
	uint64_t GetInterval(size_t ientry1, size_t ientry2) const;

 private:
	TTree *fTree;
	uint8_t   fMac5; //last 8 bits of FEB mac5 address
	uint16_t  fFlags;
	uint16_t  fLostcpu;
	uint16_t  fLostfpga;
	uint32_t  fTs0;
	uint32_t  fTs1;
	uint16_t  fAdc[32];
	uint32_t  fCoinc;		
	ULong64_t  fRun_start_time;
	ULong64_t  fThis_poll_start;
	ULong64_t  fThis_poll_end;
	ULong64_t  fLast_poll_start;
	ULong64_t  fLast_poll_end;
	int32_t   fSystem_clock_deviation;
	uint32_t  fFeb_per_poll;
	uint32_t  fFeb_event_number;
	uint32_t  fSequence_id;
	ULong64_t  fFragment_timestamp;

	TBranch* b_Mac5; //last 8 bits of FEB mac5 address
	TBranch* b_Flags;
	TBranch* b_Lostcpu;
	TBranch* b_Lostfpga;
	TBranch* b_Ts0;
	TBranch* b_Ts1;
	TBranch* b_Adc;
	TBranch* b_Coinc;
	TBranch* b_Run_start_time;
	TBranch* b_This_poll_start;
	TBranch* b_This_poll_end;
	TBranch* b_Last_poll_start;
	TBranch* b_Last_poll_end;
	TBranch* b_System_clock_deviation;
	TBranch* b_Feb_per_poll;
	TBranch* b_Feb_event_number;
	TBranch* b_Sequence_id;
	TBranch* b_Fragment_timestamp;

	void Init();

};

#endif
