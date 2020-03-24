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

	void  Dump();
  private:
	void  Init();

	TTree* fTree;
	std::map<std::pair<uint8_t,uint8_t>,float> fMacChanToGain;
	std::map<std::pair<uint8_t,uint8_t>,float> fMacChanToPed;
	std::map<std::pair<uint8_t,uint8_t>,bool> fMacChanToActive;
	std::vector<uint8_t> fMacs;
};

#endif 
