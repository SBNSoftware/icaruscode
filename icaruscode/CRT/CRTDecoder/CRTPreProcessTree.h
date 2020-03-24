#ifndef CRT_PREPROCESS_TREE_H
#define CRT_PREPROCESS_TREE_H

#include <TTree.h>
#include <TBranch.h>

#include <iostream>

namespace icarus {
 namespace crt {
	class CRTPreProcessTree;
 }
}

class icarus::crt::CRTPreProcessTree {

 public:
	CRTPreProcessTree(TTree* tr);
	void      Init();
	size_t    GetNEntries() const;
	uint64_t  GetAbsTime(size_t ientry) const;
	void      Load(size_t ientry) const;
	uint8_t   Mac5() const;
	bool      IsNoise() const;
	uint8_t   MaxChan() const;
	float     MaxPE() const;
	float     TotPE() const;
	float     PE(uint8_t chan) const;
	uint8_t   NChanAbove() const;
	bool      Above(uint8_t chan) const;
	bool      Active(uint8_t chan) const;
	int       Region() const;
	int       Layer() const;
	float     PollRate() const;
	float     InstRate() const;
	uint64_t  T0() const;	

 private:
	TTree*    fTree;

	uint8_t   fMac5;
	bool      fIsNoise;
	uint8_t   fMaxChan;
	float     fMaxPE;
	float     fTotPE;
	float     fPE[32];
	uint8_t   fNChanAbove;
	bool      fAbove[32];
	bool	  fActive[32];
	ULong64_t fT0;
	int       fRegion;
	int       fLayer;
	float     fPollRate;
	float     fInstRate;

	TBranch*  b_Mac5;
	TBranch*  b_IsNoise;
	TBranch*  b_MaxChan;
	TBranch*  b_MaxPE;
	TBranch*  b_TotPE;
	TBranch*  b_PE;
	TBranch*  b_NChanAbove;
	TBranch*  b_Above;
	TBranch*  b_Active;
	TBranch*  b_T0;
	TBranch*  b_Region;
	TBranch*  b_Layer;
	TBranch*  b_PollRate;
	TBranch*  b_InstRate;


};

#endif
