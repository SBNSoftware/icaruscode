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

        fTree->SetBranchAddress("mac5",        &mac5,       &b_mac5);
        fTree->SetBranchAddress("gain",        gain,        &b_gain);
        fTree->SetBranchAddress("ped",         ped,         &b_ped);
        fTree->SetBranchAddress("active",      active,      &b_active);
	fTree->SetBranchAddress("gainXsqr",    gainXsqr,    &b_gainXsqr);
	fTree->SetBranchAddress("gainNdf",     gainNdf,     &b_gainNdf);
	fTree->SetBranchAddress("gainErr",     gainErr,     &b_gainErr);
	fTree->SetBranchAddress("pedXsqr",     pedXsqr,     &b_pedXsqr);
	fTree->SetBranchAddress("pedNdf",      pedNdf,      &b_pedNdf);
	fTree->SetBranchAddress("pedErr",      pedErr,      &b_pedErr);
	fTree->SetBranchAddress("pedSigma",    pedSigma,    &b_pedSigma);
	fTree->SetBranchAddress("pedSigmaErr", pedSigmaErr, &b_pedSigmaErr);

        const size_t nmac = fTree->GetEntriesFast();
	std::cout << "retreiving calibration data for " << nmac << " FEBs..." << std::endl;

        for(size_t imac=0; imac<nmac; imac++){
                fTree->GetEntry(imac);

		fMacs.push_back(mac5);
		fMacToEntry[mac5] = imac;
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

float CrtCalTree::GetGainXsqr(uint8_t mac, uint8_t channel){
	float val = -1.;
	fTree->GetEntry(fMacToEntry[mac]);
	val = gainXsqr[channel];
	return val;
}
short CrtCalTree::GetGainNdf(uint8_t mac, uint8_t channel){
        float val = -1.;
        fTree->GetEntry(fMacToEntry[mac]);
        val = gainNdf[channel];
        return val;
}
float CrtCalTree::GetGainErr(uint8_t mac, uint8_t channel){
        float val = -1.;
        fTree->GetEntry(fMacToEntry[mac]);
        val = gainErr[channel];
        return val;
}
float CrtCalTree::GetPedErr(uint8_t mac, uint8_t channel){
        float val = -1.;
        fTree->GetEntry(fMacToEntry[mac]);
        val = pedErr[channel];
        return val;
}
float CrtCalTree::GetPedXsqr(uint8_t mac, uint8_t channel){
        float val = -1.;
        fTree->GetEntry(fMacToEntry[mac]);
        val = pedXsqr[channel];
        return val;
}
short CrtCalTree::GetPedNdf(uint8_t mac, uint8_t channel){
        float val = -1.;
        fTree->GetEntry(fMacToEntry[mac]);
        val = pedNdf[channel];
        return val;
}
float CrtCalTree::GetPedSigma(uint8_t mac, uint8_t channel){
        float val = -1.;
        fTree->GetEntry(fMacToEntry[mac]);
        val = pedSigma[channel];
        return val;
}
float CrtCalTree::GetPedSigmaErr(uint8_t mac, uint8_t channel){
        float val = -1.;
        fTree->GetEntry(fMacToEntry[mac]);
        val = pedSigmaErr[channel];
        return val;
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
