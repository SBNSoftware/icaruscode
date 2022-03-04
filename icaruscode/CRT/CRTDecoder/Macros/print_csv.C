#define print_csv_C

#include <TGraph.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TFile.h>

void print_csv(){

        TFile* fin = new TFile("/icarus/app/users/tboone/decoder/cal/adc_histos/ftest.root");
        TTree *tr = (TTree*)fin->FindObjectAny("calAnalyzerTree");

	ofstream myfile;
	myfile.open("run_7262_cal_data.csv");
	myfile << "channel,mac5,localchannel,gain,pedestal,stats\n";


	int numentries = tr->GetEntries();

//	int master_chan = 0;
	for(int i=0; i<numentries; i++){

		tr->GetEntry(i);

		for(int j=0; j<32; j++){

			int master_chan;
			int temp_mac5 = tr->GetBranch("mac5")->GetLeaf("mac5")->GetValue();
			double temp_gain = tr->GetBranch("gain")->GetLeaf("gain")->GetValue(j);
			double temp_ped  = tr->GetBranch("ped")->GetLeaf("ped")->GetValue(j);
			double temp_stats = tr->GetBranch("stats")->GetLeaf("stats")->GetValue(j);

			master_chan = temp_mac5 + i + (temp_mac5-1)*31;

			myfile << master_chan << "," << temp_mac5 << "," << j << "," << temp_gain << "," << temp_ped << "," << temp_stats << "\n";
//			master_chan++;

		}


	}//end loop over entries in calAnalyzerTree

	myfile.close();

}

