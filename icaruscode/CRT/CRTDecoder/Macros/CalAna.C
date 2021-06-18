#define CalAna_cxx

#include <TTree.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <map>

#include "/icarus/app/users/chilgenb/dev_areas/icaruscode_v08_45_00_prof/srcs/icaruscode/icaruscode/CRT/CRTDecoder/CrtCalTree.cxx"
#include "TimeUtils.C"

using namespace std;

const bool  WRITE_PED_FILE = false;
const float GAIN_XSQR_MAX = 2.;
const float PED_XSQR_MAX = 1.e5;

int macColors[20] = {1,  2,  3,  4,  5,  6,  7,  8,  9,  25,
		     28, 32, 37, 46, 52, 41, 70, 87, 91, 96};
void CalAna(){

	//get all files from list
	ifstream fin;
	fin.open("./caltrees.list");

	vector<string> files;
	string line;
	while(getline(fin,line)){
		cout << line << endl;
		files.push_back(line);
	}

	fin.close();

	const size_t nfile = files.size();
	vector<uint8_t> febs; //all mac5 addresses found
	map<int,map<uint8_t,vector<float>>> peds;  //file index -> (mac5 -> pedestali[32])
	map<int,map<uint8_t,vector<float>>> gains; //file index -> (mac5 -> gain[32])
	map<int,map<uint8_t,vector<int>>>   chans; //file index -> (mac5 -> channel[32])
	map<int,map<uint8_t,vector<bool>>>  gainFails; //file index -> (mac5 -> gain failures [32])
	map<int,map<uint8_t,vector<bool>>>  pedFails; //file index -> (mac5 -> ped failures [32])
	map<int,Date>                       dates; //file index -> Date

	TH1F* hngainfail_perun = new TH1F("hngainfail_perun","number of gain fit failures per run",560,0,560);
        TH1F* hnpedfail_perun  = new TH1F("hnpedfail_perun","number of pedestal fit failures per run",560,0,560);

	//loop over files
	for(size_t ifile=0; ifile<nfile; ifile++){
		cout << "reading from " << files[ifile] << endl;
		TFile f(files[ifile].c_str(),"READ");
		TTree* tree = (TTree*)f.FindObjectAny("calTree");
		if(tree==nullptr) 
			tree=(TTree*)f.FindObjectAny("calAnalyzerTree");
		const size_t nfeb = tree->GetEntriesFast(); //1 entry per FEB
		cout << "found " << nfeb << "FEBs" << endl;
		CrtCalTree cct(tree);

		vector<uint8_t> febstmp = cct.GetMacs();
		if(ifile==0) febs=febstmp;
		else if (febstmp!=febs) cout << "different FEBs between runs!" << endl;

		dates[ifile] = PathToDate(files[ifile]);

		size_t nfail_gain = 0;
		size_t nfail_ped  = 0;

		//loop over mac5 addresses, fill maps
		for(size_t ifeb=0; ifeb<nfeb; ifeb++){
			int mac5 = febs[ifeb];
			for(size_t ch=0; ch<32; ch++){

                                if(cct.GetPedXsqr(mac5,ch)>PED_XSQR_MAX){
                                        pedFails[ifile][mac5].push_back(false);
                                        nfail_ped++;
                                }
                                else{
                                        pedFails[ifile][mac5].push_back(true);
                                }

                                peds[ifile][mac5].push_back(cct.GetPed(mac5,ch));

				if(!cct.GetActive(mac5,ch)) continue;
				chans[ifile][mac5].push_back(ch);
				gains[ifile][mac5].push_back(cct.GetGain(mac5,ch));

				if(cct.GetGainNdf(mac5,ch)==0 || cct.GetGainXsqr(mac5,ch)>GAIN_XSQR_MAX) {
					gainFails[ifile][mac5].push_back(false);
					nfail_gain++;
				}
				else {
					gainFails[ifile][mac5].push_back(true);
				}
			}
			cout << "found " << chans[ifile][mac5].size() << " active channels" << endl;
		}

		hngainfail_perun->Fill(nfail_gain);
		hnpedfail_perun->Fill(nfail_ped);

		f.Close();
	}//end loop over files

	const size_t nfeb = febs.size();
	float day[nfile];
	float gainArr[nfeb][32][nfile];
	float pedArr[nfeb][32][nfile];

	float dayChanGain[nfeb][32][nfile];
	float dayChanPed[nfeb][32][nfile];
	float gainArrCut[nfeb][32][nfile];
	float pedArrCut[nfeb][32][nfile];
	size_t npassGain[nfeb][32];
	size_t npassPed[nfeb][32];

	float pedAvg[nfeb][32];
	float pedAvgMac[nfeb];

	TH1F* hpedavg_chan = new TH1F("hpedavg_chan","average channel pedestal values",45,50,275);
	TH1F* hpedavg_feb = new TH1F("hpedavg_feb","average FEB pedestal values",45,50,275);
	TH1F* hpassgain = new TH1F("hpassgain","fraction of runs with passing gains per channel",100,0,1);

	for(size_t ifeb=0; ifeb<nfeb; ifeb++){
	
		pedAvgMac[ifeb] = 0.;

		for(size_t ichan=0; ichan<32; ichan++) {

			pedAvg[ifeb][ichan] = 0.;

			for(size_t ifile=0; ifile<nfile; ifile++){

                                pedArr[ifeb][ichan][ifile] = peds[ifile][febs[ifeb]][ichan];
                                //pedAvg[ifeb][ichan] += pedArr[ifeb][ichan][ifile];
				if(pedFails[ifile][febs[ifeb]][ichan]) {
					dayChanPed[ifeb][ichan][npassPed[ifeb][ichan]] = DaysPassed(dates[0],dates[ifile]);
					pedArrCut[ifeb][ichan][npassPed[ifeb][ichan]]  = peds[ifile][febs[ifeb]][ichan];
					pedAvg[ifeb][ichan] += peds[ifile][febs[ifeb]][ichan];
					npassPed[ifeb][ichan]++;
				}

				if(ichan<chans[0][febs[ifeb]].size()) {
					gainArr[ifeb][ichan][ifile] = gains[ifile][febs[ifeb]][ichan];
	
					if(gainFails[ifile][febs[ifeb]][ichan] && pedFails[ifile][febs[ifeb]][chans[ifile][febs[ifeb]][ichan]]) {
						dayChanGain[ifeb][ichan][npassGain[ifeb][ichan]]    = DaysPassed(dates[0],dates[ifile]);
						gainArrCut[ifeb][ichan][npassGain[ifeb][ichan]] = gains[ifile][febs[ifeb]][ichan];
						npassGain[ifeb][ichan]++;
					}
				}
			}

			//pedAvg[ifeb][ichan] *= 1.0/nfile;
			pedAvg[ifeb][ichan] *= 1.0/npassPed[ifeb][ichan];
			pedAvgMac[ifeb] += pedAvg[ifeb][ichan];

			hpedavg_chan->Fill(pedAvg[ifeb][ichan]);

			if(ichan<chans[0][febs[ifeb]].size())
				hpassgain->Fill(1.0*npassGain[ifeb][ichan]/nfile);
		}


		pedAvgMac[ifeb] *= 1.0/32;
		hpedavg_feb->Fill(pedAvgMac[ifeb]);

	}

        ofstream fout;
	if(WRITE_PED_FILE)
	        fout.open("sbn_fd_crt_pedestals_2020mar30.csv");

	//calculate deviations from the mean
	float peddev[nfeb][32][nfile];
        for(size_t ifeb=0; ifeb<nfeb; ifeb++){

		if(WRITE_PED_FILE)
			fout << (int)febs[ifeb] << "," << pedAvgMac[ifeb] << ",";

                for(size_t ichan=0; ichan<32; ichan++) {

			if(WRITE_PED_FILE){
				fout << pedAvg[ifeb][ichan];
				if(ichan<31) fout << ",";
				else fout << '\n';
			}

			for(size_t ifile=0; ifile<nfile; ifile++){
				peddev[ifeb][ichan][ifile] = pedArr[ifeb][ichan][ifile] - pedAvg[ifeb][ichan];
			}
		}
	}

	if(WRITE_PED_FILE)
		fout.close();

	//setup for plotting
        vector<TH1F*> hGain, hPed;
	for(size_t ifile=0; ifile<nfile; ifile++){
		day[ifile] = DaysPassed(dates[0],dates[ifile]);
		string obj="h"+to_string(ifile);
		hGain.push_back(new TH1F(obj.c_str(),"gains",30,30,90));

	}

	vector<vector<TGraph*>> gGain;
	vector<vector<TGraph*>> gPed;
	vector<vector<TGraph*>> gGainCut;
	vector<vector<TGraph*>> gPedCut;
	vector<vector<TGraph*>> gPedDev;
        for(size_t ifeb=0; ifeb<nfeb; ifeb++){
		gGain.push_back({});
		gPed.push_back({});
                gGainCut.push_back({});
                gPedCut.push_back({});
		gPedDev.push_back({});
		cout << "filling for mac5 " << (int)febs[ifeb] << endl;
		cout << "adding new TGraphs for " << chans[0][febs[ifeb]].size() << " channels" << endl;

                for(size_t ichan=0; ichan<32; ichan++) {

			gPed[ifeb].push_back(new TGraph(nfile,day,pedArr[ifeb][ichan]));
			gPed[ifeb].back()->SetMarkerStyle(8);
			gPed[ifeb].back()->SetLineWidth(1);

			gPedDev[ifeb].push_back(new TGraph(nfile,day,peddev[ifeb][ichan]));

			gPedCut[ifeb].push_back(new TGraph(npassPed[ifeb][ichan],dayChanPed[ifeb][ichan],pedArrCut[ifeb][ichan]));
                        gPedCut[ifeb].back()->SetMarkerStyle(8);
                        gPedCut[ifeb].back()->SetLineWidth(1);

			if(ichan<chans[0][febs[ifeb]].size()) {
				for(size_t ifile=0; ifile<nfile; ifile++){
					hGain[ifile]->Fill(gainArr[ifeb][ichan][ifile]);
				}

                		gGain[ifeb].push_back(new TGraph(nfile,day,gainArr[ifeb][ichan]));
				gGain[ifeb].back()->SetMarkerStyle(8);
				gGain[ifeb].back()->SetLineWidth(1);

                                gGainCut[ifeb].push_back(new TGraph(npassGain[ifeb][ichan],dayChanGain[ifeb][ichan],gainArrCut[ifeb][ichan]));
                                gGainCut[ifeb].back()->SetMarkerStyle(8);
                                gGainCut[ifeb].back()->SetLineWidth(1);

			}

                }
        }

	//gGain[0][0]->GetYaxis()->SetRangeUser(30.0,70.0);
	//gGain[0][0]->GetXaxis()->SetTitle("days since first measurement");
	//gGain[0][0]->GetYaxis()->SetTitle("gain [ADC/PE]");
	//gGain[0][0]->SetTitle("gains over time, all FEBs/channels");

	/*TCanvas* cg = new TCanvas("cg","all gain measurements");
	cg->Divide(5,4);
	//gGain[0][0]->Draw("APL");
	for(size_t ifeb=0; ifeb<nfeb; ifeb++){
		cg->cd(ifeb+1);
		string title = "gains over time: FEB "+to_string(febs[ifeb]);
		gGain[ifeb][0]->SetTitle(title.c_str());
		gGain[ifeb][0]->GetYaxis()->SetRangeUser(30.0,70.0);
		gGain[ifeb][0]->GetXaxis()->SetTitle("days since first measurement");
		gGain[ifeb][0]->GetYaxis()->SetTitle("gain [ADC/PE]");
		gGain[ifeb][0]->Draw("APL");
		for(size_t ichan=1; ichan<chans[0][febs[ifeb]].size(); ichan++) {
			gGain[ifeb][ichan]->Draw("samepl");
		}
	}*/

	//pedestal vs time
	/*gPed[0][0]->GetYaxis()->SetRangeUser(50,275);
        gPed[0][0]->GetXaxis()->SetTitle("days since first measurement");
        gPed[0][0]->GetYaxis()->SetTitle("pedestal [ADC]");
        gPed[0][0]->SetTitle("pedestals over time, all FEBs/channels");

        TCanvas* cp = new TCanvas();
        gPed[0][0]->Draw("APL");
        for(size_t ifeb=0; ifeb<nfeb; ifeb++){
                for(size_t ichan=0; ichan<32; ichan++) {
                        gPed[ifeb][ichan]->Draw("samepl");
                }
        }*/

	//gGainCut[0][0]->GetYaxis()->SetRangeUser(30.0,70.0);
        //gGainCut[0][0]->GetXaxis()->SetTitle("days since first measurement");
        //gGainCut[0][0]->GetYaxis()->SetTitle("gain [ADC/PE]");
        //gGainCut[0][0]->SetTitle("gains over time (quality cuts), all FEBs/channels");

        TCanvas* cgc = new TCanvas("cgc","all gain measurements passing quality cuts");
	cgc->Divide(5,4);
        //gGainCut[0][0]->Draw("APL");
        for(size_t ifeb=0; ifeb<nfeb; ifeb++){
                cgc->cd(ifeb+1);
                string title = "gains over time: FEB "+to_string(febs[ifeb]);
                gGainCut[ifeb][0]->SetTitle(title.c_str());
                gGainCut[ifeb][0]->GetYaxis()->SetRangeUser(45.0,70.0);
                gGainCut[ifeb][0]->GetXaxis()->SetTitle("days since first measurement");
                gGainCut[ifeb][0]->GetYaxis()->SetTitle("gain [ADC/PE]");
                gGainCut[ifeb][0]->Draw("APL");

                for(size_t ichan=0; ichan<chans[0][febs[ifeb]].size(); ichan++) {
                        gGainCut[ifeb][ichan]->Draw("samepl");
                }
        }

        //pedestal vs time
        /*gPedCut[0][0]->GetYaxis()->SetRangeUser(50,275);
        gPedCut[0][0]->GetXaxis()->SetTitle("days since first measurement");
        gPedCut[0][0]->GetYaxis()->SetTitle("pedestal [ADC]");
        gPedCut[0][0]->SetTitle("pedestals over time (quality cuts), all FEBs/channels");

        TCanvas* cpc = new TCanvas();
        gPedCut[0][0]->Draw("APL");
        for(size_t ifeb=0; ifeb<nfeb; ifeb++){
                for(size_t ichan=0; ichan<32; ichan++) {
                        gPedCut[ifeb][ichan]->Draw("samepl");
                }
        }

	new TCanvas();
	hGain[0]->Draw();

	new TCanvas();
	hpedavg_chan->Draw();

        new TCanvas();
        hpedavg_feb->Draw();

	new TCanvas();
	hngainfail_perun->Draw();

	new TCanvas();
	hnpedfail_perun->Draw();

	new TCanvas();
        gPedDev[0][0]->Draw("APL");
        for(size_t ifeb=0; ifeb<nfeb; ifeb++){
                for(size_t ichan=0; ichan<chans[0][febs[ifeb]].size(); ichan++) {
                        gPedDev[ifeb][ichan]->Draw("samepl");
                }
        }*/


	new TCanvas();
	hpassgain->Draw();

	return;
}
