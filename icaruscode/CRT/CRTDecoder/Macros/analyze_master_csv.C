#define analyze_master_csv_C

#include <TGraph.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TFile.h>

void analyze_master_csv(){

	vector<vector<float>> master_array;
	string filename = "/icarus/app/users/tboone/decoder/cal/adc_histos/csvs/crt_cal_db_reco_filtered.csv";
	ifstream f_1;
	f_1.open(filename.c_str());

	//cut off the top line which has column header information
	string dummyline_1; getline(f_1,dummyline_1);

	string line_1, field_1; 

	vector<float> myv_1;

	while(getline(f_1,line_1)){

		myv_1.clear();
		stringstream ss_1(line_1);

		while(getline(ss_1,field_1,',')){
		
			myv_1.push_back(stof(field_1));

		}//end getting individual columns

		master_array.push_back(myv_1);

	}//end getline(f, line)

	f_1.close();

	TH1F *h_gain = new TH1F("h_gain","Side CRT Gains",201,39.5,80.5);
	h_gain->GetXaxis()->SetTitle("gain (ADC/PE)");
	h_gain->GetYaxis()->SetTitle("counts");

	TH1F *h_ped  = new TH1F("h_ped", "Side CRT Pedestals",201,-0.5,400.5);
	h_ped->GetXaxis()->SetTitle("pedestal (ADC)");
	h_ped->GetYaxis()->SetTitle("counts");

	TH1F *h_gain_flags = new TH1F("h_gain_flags","Number of Side CRT gain fit failures out of five measurements",8,-1.5,6.5);
	h_gain_flags->GetXaxis()->SetTitle("# gain fit failures");
	h_gain_flags->GetYaxis()->SetTitle("counts");

//	cout <<"master_array.size()= " << master_array.size();
//	cout << "gain,ped,gain_flags\n";
	for(int i=0; i<master_array.size(); i++){

		h_gain->Fill(master_array[i][3]);
		h_ped->Fill(master_array[i][4]);
		h_gain_flags->Fill(master_array[i][5]);

//		cout << master_array[i][3] << "," << master_array[i][4] << "," << master_array[i][5] << endl;

	}//end loop over array

	TCanvas *c1 = new TCanvas();
	h_gain->Draw();
	c1->SaveAs("gains.png");

	TCanvas *c2 = new TCanvas();
	h_ped->Draw();
	c2->SaveAs("peds.png");

	TCanvas *c3 = new TCanvas();
	h_gain_flags->Draw();
	c3->SaveAs("gainflags.png");

}//end analyze_master_csv function definition

