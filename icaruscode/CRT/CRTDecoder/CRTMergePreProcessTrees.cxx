#ifndef CRTMergePreProcessTrees_CXX
#define CRTMergePreProcessTrees_CXX

#include <TTree.h>
#include <TFile.h>

#include <map>

#include "./CRTPreProcessTree.h"
//#include "./CRTRawTree.h"
#include "./CRTTiming.h"

//#include "icaruscode/CRT/CRTDecoder/CRTPreProcessTree.h"
//#include "icaruscode/CRT/CRTDecoder/CRTRawTree.h"
//#include "icaruscode/CRT/CRTDecoder/CRTTiming.h"


using namespace std;

//takes 2 args
//to run: ./CRTMergePreProcessTrees <input file name>.root <output file name>.root
int main(int argc, char *argv[]) {

	cout << "Opening file for reading with name, '" << argv[1] << "'" << endl;	

	cout << "Opening new file for writing with name, '" << argv[2] << "'" << endl;

        uint8_t   fMac5;
        bool      fIsNoise;
        uint8_t   fMaxChan;
        float     fMaxPE;
        float     fTotPE;
        float     fPE[32];
        uint8_t   fNChanAbove;
        bool      fAbove[32];
        bool      fActive[32];
        uint64_t  fT0;
        int       fRegion;
        int       fLayer;
        float     fPollRate;
        float     fInstRate;

	TTree* fAnaTree = new TTree("anaTree","calibrated charge, trigger 'flag', and time ordered entries");

	fAnaTree->Branch("mac5",           &fMac5,                "mac5/b");
	fAnaTree->Branch("pe",             fPE,                   "pe[32]/F");
	fAnaTree->Branch("active",         fActive,               "active[32]/O");
	fAnaTree->Branch("maxChan",        &fMaxChan,             "maxChan/b");
	fAnaTree->Branch("maxPE",          &fMaxPE,               "maxPE/F");
	fAnaTree->Branch("totPE",          &fTotPE,               "totPE/F");
	fAnaTree->Branch("nAbove",         &fNChanAbove,          "nAbove/b");
	fAnaTree->Branch("above",          fAbove,                "above[32]/O");
	fAnaTree->Branch("isNoise",        &fIsNoise,             "isNoise/O");
	fAnaTree->Branch("region",         &fRegion,              "region/I");
	fAnaTree->Branch("layer",          &fLayer,               "layer/I");
	fAnaTree->Branch("t0",             &fT0,                  "t0/l");
	fAnaTree->Branch("pollRate",       &fPollRate,            "pollRate/F");
	fAnaTree->Branch("instRate",       &fInstRate,            "instRate/F");

	TFile fin(argv[1],"READ");
	TTree* intree = (TTree*)fin.FindObjectAny("anaTree");
	icarus::crt::CRTPreProcessTree cpt(intree);
	icarus::crt::CRTTiming ct(cpt);

	const map<size_t,size_t>* sortedToRaw = ct.GetOrderedToRawMap();
	std::cout << "done sorting " << sortedToRaw->size() << " entries" << std::endl;


	TFile fout(argv[2],"RECREATE");

	std::cout << "filling new, sorted tree..." << std::endl;
	const size_t nentries = cpt.GetNEntries();
	for(size_t i=0; i<nentries; i++) {

		if(i%100000==0)
			std::cout << 100.0*i/nentries << " % complete" << std::endl;

		cpt.Load(sortedToRaw->at(i));
		
        	fMac5 = cpt.Mac5();
        	fIsNoise = cpt.IsNoise();
        	fMaxChan = cpt.MaxChan();
        	fMaxPE = cpt.MaxPE();
        	fTotPE = cpt.TotPE();
        	fNChanAbove = cpt.NChanAbove();
        	fT0 = cpt.T0();
        	fRegion = cpt.Region();
        	fLayer = cpt.Layer();
        	fPollRate = cpt.PollRate();
        	fInstRate = cpt.InstRate();

		for(size_t ch=0; ch<32; ch++){
			fPE[ch] = cpt.PE(ch);
			fAbove[ch] = cpt.Above(ch);
			fActive[ch] = cpt.Active(ch);
		}

		fAnaTree->Fill();
	}
	fin.Close();

	fAnaTree->Write();
	fout.Close();

	return 0;
}

#endif
