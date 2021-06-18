#ifndef CRTCALTREE_H
#define CRTCALTREE_H

#include <TDirectory.h>
#include <TBranch.h>
#include <TTree.h>
#include <TFile.h>

#include <map>
#include <string>
#include <iostream>
#include <vector>

using std::map;
using std::string;
using std::vector;

namespace icarus{
  namespace crt {
	class CrtCalTree;
  }
}

class icarus::crt::CrtCalTree {


  public:

	CrtCalTree(TTree* tree);
	CrtCalTree(string fname);
	~CrtCalTree();

	float GetGain(uint8_t mac, uint8_t channel);
	float GetPed(uint8_t mac, uint8_t channel);
	bool  GetActive(uint8_t mac, uint8_t channel);
	vector<uint8_t> GetMacs();
	float GetGainXsqr(uint8_t mac, uint8_t channel);
	short GetGainNdf(uint8_t mac, uint8_t channel);
	float GetGainErr(uint8_t mac, uint8_t channel);
	float GetPedErr(uint8_t mac, uint8_t channel);
	float GetPedXsqr(uint8_t mac, uint8_t channel);
	short GetPedNdf(uint8_t mac, uint8_t channel);
	float GetPedSigma(uint8_t mac, uint8_t channel);
	float GetPedSigmaErr(uint8_t mac, uint8_t channel);

	void  Dump();
  private:
	void  Init();
	
	TTree* fTree;
	std::map<std::pair<uint8_t,uint8_t>,float> fMacChanToGain;
	std::map<std::pair<uint8_t,uint8_t>,float> fMacChanToPed;
	std::map<std::pair<uint8_t,uint8_t>,bool> fMacChanToActive;
	std::vector<uint8_t> fMacs;
	std::map<uint8_t,int> fMacToEntry;

        uint8_t mac5;
        float   gain[32];
        float   ped[32];
        bool    active[32];
        float   gainXsqr[32];
        UShort_t   gainNdf[32];
        float   gainErr[32];
        float   pedXsqr[32];
        UShort_t   pedNdf[32];
        float   pedErr[32];
        float   pedSigma[32];
        float   pedSigmaErr[32];

        TBranch* b_mac5;
        TBranch* b_gain;
        TBranch* b_ped;
        TBranch* b_active;
        TBranch* b_gainXsqr;
        TBranch* b_gainNdf;
        TBranch* b_gainErr;
        TBranch* b_pedXsqr;
        TBranch* b_pedNdf;
        TBranch* b_pedErr;
        TBranch* b_pedSigma;
        TBranch* b_pedSigmaErr;


};

#endif 
