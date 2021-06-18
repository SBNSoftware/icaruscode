#ifndef CRT_RAW_TREE_CC
#define CRT_RAW_TREE_CC

#include "icaruscode/CRT/CRTDecoder/CRTRawTree.h"

using namespace icarus::crt;

CRTRawTree::CRTRawTree(TTree* tr){

	fTree = tr;
	//fTree->SetMaxVirtualSize(16e9);
	fTree->LoadBaskets(16e9);
	Init();
}

void CRTRawTree::Init(){

  fTree->SetBranchAddress("mac5",                   &fMac5,                   &b_Mac5);
  fTree->SetBranchAddress("flags",                  &fFlags,                  &b_Flags);
  fTree->SetBranchAddress("lostcpu",                &fLostcpu,                &b_Lostcpu);
  fTree->SetBranchAddress("lostfpga",               &fLostfpga,               &b_Lostfpga);
  fTree->SetBranchAddress("ts0",                    &fTs0,                    &b_Ts0);
  fTree->SetBranchAddress("ts1",                    &fTs1,                    &b_Ts1);
  fTree->SetBranchAddress("adc",                    &fAdc,                    &b_Adc);
  fTree->SetBranchAddress("coinc",                  &fCoinc,                  &b_Coinc);
  fTree->SetBranchAddress("run_start_time",         &fRun_start_time,         &b_Run_start_time);
  fTree->SetBranchAddress("this_poll_start",        &fThis_poll_start,        &b_This_poll_start);
  fTree->SetBranchAddress("this_poll_end",          &fThis_poll_end,          &b_This_poll_end);
  fTree->SetBranchAddress("last_poll_start",        &fLast_poll_start,        &b_Last_poll_start);
  fTree->SetBranchAddress("last_poll_end",          &fLast_poll_end,          &b_Last_poll_end);
  fTree->SetBranchAddress("system_clock_deviation", &fSystem_clock_deviation, &b_System_clock_deviation);
  fTree->SetBranchAddress("feb_per_poll",           &fFeb_per_poll,           &b_Feb_per_poll);
  fTree->SetBranchAddress("feb_event_number",       &fFeb_event_number,       &b_Feb_event_number);
  fTree->SetBranchAddress("sequence_id",            &fSequence_id,            &b_Sequence_id);
  fTree->SetBranchAddress("fragment_timestamp",     &fFragment_timestamp,     &b_Fragment_timestamp);


}

size_t CRTRawTree::GetNEntries() const {
	return (size_t)fTree->GetEntriesFast();
}

uint64_t CRTRawTree::GetAbsTime(size_t ientry) const {
	if(ientry >= GetNEntries()) {
		std::cout << "ERROR in CRTRawTree::GetAbsTime: entry out of range" << std::endl;
		return UINT64_MAX;
	}

	fTree->GetEntry(ientry);
	return (uint64_t)fFragment_timestamp;
}

uint8_t CRTRawTree::GetMac(size_t ientry) const {
        if(ientry >= GetNEntries()) {
                std::cout << "ERROR in CRTRawTree::GetMac: entry out of range" << std::endl;
                return UINT8_MAX;
        }

        fTree->GetEntry(ientry);
        return fMac5;	
}

uint16_t CRTRawTree::GetADC(size_t ientry, uint8_t chan) const {
        if(ientry >= GetNEntries()) {
                std::cout << "ERROR in CRTRawTree::GetADC: entry out of range" << std::endl;
                return UINT16_MAX;
        }

        fTree->GetEntry(ientry);
        return fAdc[chan];
}

float CRTRawTree::GetInstRate(size_t ientry_prev, size_t ientry_next) const {

	
        return 1.0/GetInterval(ientry_prev, ientry_next);
}

float CRTRawTree::GetPollRate(size_t ientry) const {
        if(ientry >= GetNEntries()) {
                std::cout << "ERROR in CRTRawTree::GetPollRate: entry out of range" << std::endl;
                return FLT_MAX;
        }

        fTree->GetEntry(ientry);
        return fFeb_per_poll/(fThis_poll_start-fLast_poll_start);
}

uint64_t CRTRawTree::GetInterval(size_t ientry1, size_t ientry2) const {
        if(ientry1 >= GetNEntries() || ientry2 >= GetNEntries()) {
                std::cout << "ERROR in CRTRawTree::GetInterval: entry out of range" << std::endl;
                return UINT64_MAX;
        }

	uint64_t tprev, tnext;
        fTree->GetEntry(ientry1);
        tprev = fFragment_timestamp;

        fTree->GetEntry(ientry2);
	tnext = fFragment_timestamp;

	if(tnext < tprev){
		std::cout << "ERROR in CRTRawTree::GetInterval: reverse ordered!" << std::endl;
		return UINT64_MAX;
	}

	return tnext - tprev;
}

#endif
