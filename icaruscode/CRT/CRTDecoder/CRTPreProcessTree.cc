#ifndef CRT_PREPROCESS_TREE_CC
#define CRT_PREPROCESS_TREE_CC

//#include "./CRTPreProcessTree.h"
#include "icaruscode/CRT/CRTDecoder/CRTPreProcessTree.h"

using namespace icarus::crt;

CRTPreProcessTree::CRTPreProcessTree(TTree* tr) {
	fTree = tr;
        //fTree->SetMaxVirtualSize(16e9);
        fTree->LoadBaskets(16e9);
	Init();
}

void CRTPreProcessTree::Init() {

  fTree->SetBranchAddress("mac5",           &fMac5,        &b_Mac5);
  fTree->SetBranchAddress("pe",             fPE,           &b_PE);
  fTree->SetBranchAddress("active",         fActive,       &b_Active);
  fTree->SetBranchAddress("maxChan",        &fMaxChan,     &b_MaxChan);
  fTree->SetBranchAddress("maxPE",          &fMaxPE,       &b_MaxPE);
  fTree->SetBranchAddress("totPE",          &fTotPE,       &b_TotPE);
  fTree->SetBranchAddress("nAbove",         &fNChanAbove,  &b_NChanAbove);
  fTree->SetBranchAddress("above",          fAbove,        &b_Above);
  fTree->SetBranchAddress("isNoise",        &fIsNoise,     &b_IsNoise);
  fTree->SetBranchAddress("region",         &fRegion,      &b_Region);
  fTree->SetBranchAddress("layer",          &fLayer,       &b_Layer);
  fTree->SetBranchAddress("t0",             &fT0,          &b_T0);
  fTree->SetBranchAddress("pollRate",       &fPollRate,    &b_PollRate);
  fTree->SetBranchAddress("instRate",       &fInstRate,    &b_InstRate);


}

size_t  CRTPreProcessTree::GetNEntries() const {
        return (size_t)fTree->GetEntriesFast();
}

uint64_t CRTPreProcessTree::GetAbsTime(size_t ientry) const {
        if(ientry >= GetNEntries()) {
                std::cout << "ERROR in CRTRawTree::GetAbsTime: entry out of range" << std::endl;
                return UINT64_MAX;
        }

        fTree->GetEntry(ientry);
        return (uint64_t)fT0;
}

uint64_t CRTPreProcessTree::GetAbsTime() const {
        return (uint64_t)fT0;
}

void    CRTPreProcessTree::Load(size_t ientry) const {
	fTree->GetEntry(ientry);	
}
uint8_t CRTPreProcessTree::Mac5() const {
	return fMac5;
}
bool    CRTPreProcessTree::IsNoise() const {
	return fIsNoise;
}
uint8_t CRTPreProcessTree::MaxChan() const {
	return fMaxChan;
}
float   CRTPreProcessTree::MaxPE() const {
	return fMaxPE;
}
float   CRTPreProcessTree::TotPE() const {
	return fTotPE;
}
float   CRTPreProcessTree::PE(uint8_t chan) const {
	return fPE[chan];
}
uint8_t CRTPreProcessTree::NChanAbove() const {
	return fNChanAbove;
}
bool    CRTPreProcessTree::Above(uint8_t chan) const {
	return fAbove[chan];
}
bool    CRTPreProcessTree::Active(uint8_t chan) const {
	return fActive[chan];
}
int     CRTPreProcessTree::Region() const {
	return fRegion;
}
int     CRTPreProcessTree::Layer() const {
	return fLayer;
}
float   CRTPreProcessTree::PollRate() const {
	return fPollRate;
}
float   CRTPreProcessTree::InstRate() const {
	return fInstRate;
}
uint64_t CRTPreProcessTree::T0() const {
	return (uint64_t)fT0;
}

#endif 
