#define noise_ana_cxx

#include <TGraph.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TFile.h>

#include "TimeUtils.C"

const float MAXRATE = 12.0; //khz

size_t macNI[4] = {1,3,6,7};
size_t macNO[4] = {4,5,8,9};
size_t macWI[6] = {17,19,21,22,23,24};
size_t macWO[6] = {10,11,13,14,15,16};

vector<string> locN = {"W Mezz","W Pit","E Mezz", "E Pit"};
vector<string> locW = {"N Mezz","N MezzPit","N Pit","S Pit","S MezzPit","S Mezz"};

const char mode ='a'; //'n', 'w' or 'a'

void noise_ana(){

	vector<string> fwo, fwi, fno, fni;

        //get all files from list
        ifstream flist;
        flist.open("./noisefiles.list");

        string line;
        while(getline(flist,line)){
		if(line.find("north_inner")!=UINT64_MAX){
	                fni.push_back(line);
		}
                if(line.find("north_outer")!=UINT64_MAX){
                        fno.push_back(line);
                }
                if(line.find("west_inner")!=UINT64_MAX){
                        fwi.push_back(line);
                }
                if(line.find("west_outer")!=UINT64_MAX){
                        fwo.push_back(line);
                }

        }

        flist.close();

	cout << "--- found files ---" << endl;
	cout << "  north, inner: " << fni.size() << endl;
        cout << "  north, outer: " << fno.size() << endl;
        cout << "  west, inner:  " << fwi.size() << endl;
        cout << "  west, outer:  " << fwo.size() << endl;


	const size_t n_wo = fwo.size();
        const size_t n_wi = fwi.size();
	const size_t n_no = fno.size();
	const size_t n_ni = fni.size();
        float rateWO[6][n_wo];
        float rateWI[6][n_wi];
	float rateNO[4][n_no];
	float rateNI[4][n_ni];
        float timeWO[n_wo];
        float timeWI[n_wi];
	float timeNO[n_no];
	float timeNI[n_ni];

	Date date0 = PathToDate(fno[0]);

	if(mode=='n'||mode=='a'){
	  for(size_t i=0; i<n_no; i++){
		cout << "opening file " << fno[i] << endl;
		TFile fin(fno[i].c_str());
		timeNO[i] = (float)DaysPassed(date0,PathToDate(fno[i]));

		for(size_t imac=0; imac<4; imac++){

			string hname = "h_"+to_string(macNO[imac]);
			TH1F* h = (TH1F*)fin.FindObjectAny(hname.c_str());

			rateNO[imac][i] = h->GetMean();
		}

		fin.Close();

	  }

          for(size_t i=0; i<n_ni; i++){
                cout << "opening file " << fni[i] << endl;
                TFile fin(fni[i].c_str());
                timeNI[i] = (float)DaysPassed(date0,PathToDate(fni[i]));

                for(size_t imac=0; imac<4; imac++){

                        string hname = "h_"+to_string(macNI[imac]);
                        TH1F* h = (TH1F*)fin.FindObjectAny(hname.c_str());

                        rateNI[imac][i] = h->GetMean();
                }

                fin.Close();

          }
	}

        if(mode=='w'||mode=='a'){
          for(size_t i=0; i<n_wo; i++){
                cout << "opening file " << fwo[i] << endl;
                TFile fin(fwo[i].c_str());
                timeWO[i] = (float)DaysPassed(date0,PathToDate(fwo[i]));

                for(size_t imac=0; imac<6; imac++){

                        string hname = "h_"+to_string(macWO[imac]);
                        TH1F* h = (TH1F*)fin.FindObjectAny(hname.c_str());

                        rateWO[imac][i] = h->GetMean();
                }

                fin.Close();

          }

          for(size_t i=0; i<n_wi; i++){
                cout << "opening file " << fwi[i] << endl;
                TFile fin(fwi[i].c_str());
                timeWI[i] = (float)DaysPassed(date0,PathToDate(fwi[i]));

                for(size_t imac=0; imac<6; imac++){

                        string hname = "h_"+to_string(macWI[imac]);
                        TH1F* h = (TH1F*)fin.FindObjectAny(hname.c_str());

                        rateWI[imac][i] = h->GetMean();
                }

                fin.Close();

          }
        }

	// --- plotting ---

	TLegend* leg_ni = new TLegend(0.3,0.65,0.6,0.9);
        TLegend* leg_no = new TLegend(0.3,0.65,0.6,0.9);
        TLegend* leg_wi = new TLegend(0.3,0.55,0.6,0.9);
        TLegend* leg_wo = new TLegend(0.3,0.55,0.6,0.9);

	leg_ni->SetHeader("mac5 (location)", "C");
        leg_no->SetHeader("mac5 (location)", "C");
        leg_wi->SetHeader("mac5 (location)", "C");
        leg_wo->SetHeader("mac5 (location)", "C");
	
	vector<TGraph*> gO;

	if(mode=='a' || mode=='n') {
	  TCanvas* cno = new TCanvas();
    	  for(size_t imac=0; imac<4; imac++){
		gO.push_back(new TGraph(n_no,timeNO,rateNO[imac]));
		gO.back()->SetMarkerColor(imac+1);
		gO.back()->SetMarkerStyle(8);
		string text = to_string(macNO[imac])+" ("+locN[imac]+")";
		leg_no->AddEntry(gO.back(),text.c_str(),"p");
		if(imac==0){
			gO.back()->GetXaxis()->SetTitle("days since start");
			gO.back()->GetYaxis()->SetTitle("rate [kHz]");
			gO.back()->GetXaxis()->SetTitleSize(0.05);
			gO.back()->GetXaxis()->SetLabelSize(0.05);
			gO.back()->GetYaxis()->SetTitleSize(0.05);
			gO.back()->GetYaxis()->SetLabelSize(0.05);
			gO.back()->GetXaxis()->SetTitleOffset(0.9);
			gO.back()->GetYaxis()->SetTitleOffset(0.9);
                        gO.back()->GetYaxis()->SetRangeUser(0,MAXRATE);
			gO.back()->SetTitle("CRT Noise Rate Over Time: North Wall, Out");
			gO.back()->Draw("APL");
		}
		else{
			gO.back()->Draw("samepl");
		}
	  }
	  leg_no->Draw();
	  cno->SaveAs("noise_north_outer.png");
	}	

        vector<TGraph*> gI;

	if(mode=='a' || mode=='n') {
	  TCanvas* cni = new TCanvas();
          for(size_t imac=0; imac<4; imac++){
                gI.push_back(new TGraph(n_ni,timeNI,rateNI[imac]));
                gI.back()->SetMarkerColor(imac+1);
                gI.back()->SetMarkerStyle(8);
                string text = to_string(macNI[imac])+" ("+locN[imac]+")";
                leg_ni->AddEntry(gI.back(),text.c_str(),"p");

                if(imac==0){
                        gI.back()->GetXaxis()->SetTitle("days since start");
                        gI.back()->GetYaxis()->SetTitle("rate [kHz]");
                        gI.back()->GetXaxis()->SetTitleSize(0.05);
                        gI.back()->GetXaxis()->SetLabelSize(0.05);
                        gI.back()->GetYaxis()->SetTitleSize(0.05);
                        gI.back()->GetYaxis()->SetLabelSize(0.05);
                        gI.back()->GetXaxis()->SetTitleOffset(0.9);
                        gI.back()->GetYaxis()->SetTitleOffset(0.9);
                        gI.back()->GetYaxis()->SetRangeUser(0,MAXRATE);
                        gI.back()->SetTitle("CRT Noise Rate Over Time: North Wall, In");
                        gI.back()->Draw("APL");
                }
                else{
                        gI.back()->Draw("samepl");
                }
          }
	  leg_ni->Draw();
	  cni->SaveAs("noise_north_inner.png");

	  /*new TCanvas();
	  gO[0]->GetYaxis()->SetRangeUser(-0.1,30);
	  //gO[0]->GetXaxis()->SetRangeUser(-0.5,7.5);
	  gO[0]->Draw("apl");
	  gI[0]->Draw("samepl");
	  for(size_t imac=1; imac<4; imac++){
		gO[imac]->Draw("samepl");
		gI[imac]->Draw("samepl");
	  }*/
	}

        vector<TGraph*> gWO;
	vector<TGraph*> gWI;

        if(mode=='a' || mode=='w') {
	  TCanvas* cwo = new TCanvas();
          for(size_t imac=0; imac<6; imac++){
                gWO.push_back(new TGraph(n_wo,timeWO,rateWO[imac]));
                gWO.back()->SetMarkerColor(imac+1);
                gWO.back()->SetMarkerStyle(8);
                string text = to_string(macWO[imac])+" ("+locW[imac]+")";
                leg_wo->AddEntry(gWO.back(),text.c_str(),"p");
                if(imac==0){
                        gWO.back()->GetXaxis()->SetTitle("days since start");
                        gWO.back()->GetYaxis()->SetTitle("rate [kHz]");
                        gWO.back()->GetXaxis()->SetTitleSize(0.05);
                        gWO.back()->GetXaxis()->SetLabelSize(0.05);
                        gWO.back()->GetYaxis()->SetTitleSize(0.05);
                        gWO.back()->GetYaxis()->SetLabelSize(0.05);
                        gWO.back()->GetXaxis()->SetTitleOffset(0.9);
                        gWO.back()->GetYaxis()->SetTitleOffset(0.9);
			gWO.back()->GetYaxis()->SetRangeUser(0,MAXRATE);
                        gWO.back()->SetTitle("CRT Noise Rate Over Time: West Wall, Out");
                        gWO.back()->Draw("APL");
                }
                else{
                        gWO.back()->Draw("samepl");
                }
          }
	  leg_wo->Draw();
	  cwo->SaveAs("noise_west_outer.png");

	  TCanvas* cwi = new TCanvas();
	  for(size_t imac=0; imac<6; imac++){
                gWI.push_back(new TGraph(n_wi,timeWI,rateWI[imac]));
                gWI.back()->SetMarkerColor(imac+1);
                gWI.back()->SetMarkerStyle(8);
                string text = to_string(macWI[imac])+" ("+locW[imac]+")";
                leg_wi->AddEntry(gWI.back(),text.c_str(),"p");
                if(imac==0){
                        gWI.back()->GetXaxis()->SetTitle("days since start");
                        gWI.back()->GetYaxis()->SetTitle("rate [kHz]");
                        gWI.back()->GetXaxis()->SetTitleSize(0.05);
                        gWI.back()->GetXaxis()->SetLabelSize(0.05);
                        gWI.back()->GetYaxis()->SetTitleSize(0.05);
                        gWI.back()->GetYaxis()->SetLabelSize(0.05);
                        gWI.back()->GetXaxis()->SetTitleOffset(0.9);
                        gWI.back()->GetYaxis()->SetTitleOffset(0.9);
                        gWI.back()->GetYaxis()->SetRangeUser(0,MAXRATE);
                        gWI.back()->SetTitle("CRT Noise Rate Over Time: West Wall, In");
                        gWI.back()->Draw("APL");
                }
                else{
                        gWI.back()->Draw("samepl");
                }
	  }
	  leg_wi->Draw();
	  cwi->SaveAs("noise_west_inner.png");

          /*new TCanvas();
          gWO[0]->GetYaxis()->SetRangeUser(-0.1,30);
          gWO[0]->Draw("apl");
          gWI[0]->Draw("samepl");
          for(size_t imac=1; imac<4; imac++){
                gWO[imac]->Draw("samepl");
                gWI[imac]->Draw("samepl");
          }*/	

	}
}
