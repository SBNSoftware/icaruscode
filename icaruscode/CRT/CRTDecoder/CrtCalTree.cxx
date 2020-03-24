#ifndef CRTCALTREE_CXX
#define CRTCALTREE_CXX

#include "icaruscode/CRT/CRTDecoder/CrtCalTree.h"

using namespace icarus::crt;

CrtCalTree::CrtCalTree(TTree* tree){
	fTree = tree;
	Init();
}

CrtCalTree::CrtCalTree(string fname){
        TFile fin(fname.c_str());//,"READ");
	fTree = (TTree*)fin.FindObjectAny("calTree");
	Init();
}

CrtCalTree::~CrtCalTree(){}

void CrtCalTree::Init(){

        uint8_t mac5=0;
        float   gain[32];
        float   ped[32];
        bool    active[32];

        TBranch* b_mac5;
        TBranch* b_gain;
        TBranch* b_ped;
        TBranch* b_active;

        fTree->SetBranchAddress("mac5",   &mac5,   &b_mac5);
        fTree->SetBranchAddress("gain",   &gain,   &b_gain);
        fTree->SetBranchAddress("ped",    &ped,    &b_ped);
        fTree->SetBranchAddress("active", &active, &b_active);

        const size_t nmac = fTree->GetEntriesFast();
	std::cout << "retreiving calibration data for " << nmac << " FEBs..." << std::endl;

        for(size_t imac=0; imac<nmac; imac++){
                fTree->GetEntry(imac);

		fMacs.push_back(mac5);
                for(size_t chan=0; chan < 32; chan++){
                        fMacChanToGain[std::make_pair(mac5,chan)] = gain[chan];
                        fMacChanToPed[std::make_pair(mac5,chan)] = ped[chan];
                        fMacChanToActive[std::make_pair(mac5,chan)] = active[chan];
                }
        }	
}

float CrtCalTree::GetGain(uint8_t mac, uint8_t channel){
	return fMacChanToGain[std::make_pair(mac,channel)];
}

float CrtCalTree::GetPed(uint8_t mac, uint8_t channel){
        return fMacChanToPed[std::make_pair(mac,channel)];
}

bool CrtCalTree::GetActive(uint8_t mac, uint8_t channel){
        return fMacChanToActive[std::make_pair(mac,channel)];
}

std::vector<uint8_t> CrtCalTree::GetMacs() {
	return fMacs;
}

void CrtCalTree::Dump() {
	std::cout << "preparing to dump " << fMacChanToGain.size() << " channels." << std::endl;

	std::cout << "mac5, channel -> gain, ped (active)" << std::endl;
	for(auto const& mc : fMacChanToGain) {
		std::cout << "   " << (int)mc.first.first << ", " 
                          << (int)mc.first.second << "-> " << mc.second 
			  << ", " << fMacChanToPed[mc.first] << " ("
			  << fMacChanToActive[mc.first] << ")" << std::endl;
	}
}

#endif
