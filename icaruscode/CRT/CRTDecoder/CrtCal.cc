#ifndef CRT_CAL_CC
#define CRT_CAL_CC

#include "icaruscode/CRT/CRTDecoder/CrtCal.h"

using namespace icarus::crt;

CrtCal::CrtCal(const vector<TH1F*>* histos) : fHistos(histos) {

	fMac5 = 0;
        this->IndexToMacChan();

	fHasActive = false;
	fHasThresh = false;
	fHasPedCal = false;
	fHasGainCal = false;	

        fActive       = new bool[32];

        fGain         = new float[32];
        fGainErr      = new float[32];
        fGainXsqr     = new float[32];
	fGainNdf      = new short[32];
        fGainPed      = new float[32];
        fGainPedErr   = new float[32];
        fNpeak        = new short[32];

        fPed          = new float[32];
        fPedErr       = new float[32];
        fPedXsqr      = new float[32];
        fPedNdf       = new short[32];
        fPedSigma     = new float[32];
        fPedSigmaErr  = new float[32];
        fPedNorm      = new float[32];
        fPedNormErr   = new float[32];
        fThreshADC    = new int[32];
        fThreshPE     = new float[32];
        fNabove       = new int[32];

        fPeakNorm     = new float*[32];
        fPeakNormErr  = new float*[32];
        fPeakSigma    = new float*[32];
        fPeakSigmaErr = new float*[32];
        fPeakMean     = new float*[32]; 
        fPeakMeanErr  = new float*[32];
        fPeakXsqr     = new float*[32];
	fPeakNdf      = new short*[32];

	for(size_t i=0; i<32; i++){
        	fPeakNorm[i] = new float[5];
        	fPeakNormErr[i] = new float[5];
        	fPeakSigma[i] = new float[5];
        	fPeakSigmaErr[i] = new float[5];
        	fPeakMean[i] = new float[5];
        	fPeakMeanErr[i] = new float[5];
        	fPeakXsqr[i] = new float[5];
		fPeakNdf[i]  = new short[5];
	}


	for(size_t ch=0; ch<32; ch++){
	        fActive[ch] = false;

 		fGain[ch]         =  FLT_MAX;
 		fGainErr[ch]      =  FLT_MAX;
 		fGainXsqr[ch]     =  FLT_MAX;
 		fGainNdf[ch]      =  SHRT_MAX;
 		fGainPed[ch]      =  FLT_MAX;
 		fGainPedErr[ch]   =  FLT_MAX;
 		fNpeak[ch]        =  SHRT_MAX;

 		fPed[ch]          =  FLT_MAX;
 		fPedErr[ch]       =  FLT_MAX;
 		fPedXsqr[ch]      =  FLT_MAX;
 		fPedNdf[ch]       =  SHRT_MAX;
 		fPedSigma[ch]     =  FLT_MAX;
 		fPedSigmaErr[ch]  =  FLT_MAX;
 		fPedNorm[ch]      =  FLT_MAX;
 		fPedNormErr[ch]   =  FLT_MAX;
 		fThreshADC[ch]    =  INT_MAX;
 		fThreshPE[ch]     =  FLT_MAX;
 		fNabove[ch]       =  INT_MAX;

		for(size_t p=0; p<5; p++){
			fPeakNorm[ch][p]     =  FLT_MAX;
			fPeakNormErr[ch][p]  =  FLT_MAX;
			fPeakSigma[ch][p]    =  FLT_MAX;
			fPeakSigmaErr[ch][p] =  FLT_MAX;
			fPeakMean[ch][p]     =  FLT_MAX;
			fPeakMeanErr[ch][p]  =  FLT_MAX;
			fPeakXsqr[ch][p]     =  FLT_MAX;
			fPeakNdf[ch][p]      =  SHRT_MAX;
		}

	}
}

CrtCal::~CrtCal(){

        delete[] fActive; 

        delete[] fGain;    
        delete[] fGainErr;
        delete[] fGainXsqr;   
        delete[] fGainNdf;
        delete[] fGainPed;   
        delete[] fGainPedErr; 
        delete[] fNpeak;

        delete[] fPed;
        delete[] fPedErr;
        delete[] fPedXsqr;    
        delete[] fPedNdf;  
        delete[] fPedSigma;
        delete[] fPedSigmaErr;
        delete[] fPedNorm;
        delete[] fPedNormErr;
        delete[] fThreshADC;
        delete[] fThreshPE;
        delete[] fNabove;

        for(size_t i=0; i<32; i++){
                delete[] fPeakNorm[i];
                delete[] fPeakNormErr[i];
                delete[] fPeakSigma[i];
                delete[] fPeakSigmaErr[i];
                delete[] fPeakMean[i];
                delete[] fPeakMeanErr[i];
                delete[] fPeakXsqr[i];
                delete[] fPeakNdf[i];
        }

        delete[] fPeakNorm;
        delete[] fPeakNormErr;
        delete[] fPeakSigma;
        delete[] fPeakSigmaErr;
        delete[] fPeakMean;
        delete[] fPeakMeanErr;
        delete[] fPeakXsqr;
        delete[] fPeakNdf;
}

void CrtCal::IndexToMacChan(){

        if(fChanMap.size()!=0){
                return;
        }

        for(size_t i=0; i<fHistos->size(); i++){
                //hname = hadc_mac_chan
                const char* hname = fHistos->at(i)->GetName();
                string name = "";
                name+=hname;
		if(i==0) {
			const uint8_t nDigitMac = name.find("_",name.find("_")+1)-name.find("_")-1;
			fMac5 = atoi(name.substr(name.find("_")+1,nDigitMac).c_str());
		}
                const uint8_t nDigitChan = name.size()-name.find("_",name.find("_")+1)-1;
                uint8_t chan = atoi(name.substr(name.size()-nDigitChan,nDigitChan).c_str());
		if(chan>31) std::cout << "ERROR in CrtCal::IndexToMacCHan: channel out of range (" << (short)chan << ")" << std::endl;
                fChanMap[i] = chan;
        }

	std::cout << "found " << fChanMap.size() << " histograms/channels" << std::endl;

        return;
}

void CrtCal::Cal(){

	size_t nactive = 0;
	std::cout << "ActiveChannel scan..." << std::endl;
        for(size_t i=0; i<fHistos->size(); i++){
                TH1F* h1 = (TH1F*)fHistos->at(i)->Clone();
		TH1F* h2 = (TH1F*)fHistos->at(i)->Clone();
		TH1F* h3 = (TH1F*)fHistos->at(i)->Clone();
		const size_t chan = fChanMap[i];
		fActive[chan] = IsActive(h1);
		delete h1;
		if(fActive[chan]) {
			nactive++;
			fThreshADC[chan] = FindThreshADC(h2);
			delete h2;
			std::cout << "ch. " << chan << " ADC thresh: " << fThreshADC[chan] << std::endl;
			fNabove[chan] = FindNabove(h3,fThreshADC[chan]);
			delete h3;
		} 
	}
	std::cout << "found " << nactive << " active channels" << std::endl;

	std::cout << "PedCal..." << std::endl;
	PedCal();
	std::cout << "GainCal..." << std::endl;
	GainCal();

	for(size_t ch=0; ch<32; ch++){
		if(!fActive[ch]) continue;
		fThreshPE[ch] = 1.0*(fThreshADC[ch]-fPed[ch])/fGain[ch];
	}	
}

bool CrtCal::IsActive(TH1F* h){

	const size_t cutoff = 600;

	if(h->Integral(h->FindBin(cutoff),h->GetNbinsX())==0) {
		return false;
	}

	fHasActive = true;

	return true;
}

int CrtCal::FindThreshADC(TH1F* h){
	const size_t low = 300;
	const size_t high = 1000;
	int thresh=0;
	
	h->GetXaxis()->SetRangeUser(low,high);
	thresh = (int)h->GetBinLowEdge(h->GetMinimumBin());
	fHasThresh = true;

	return thresh;
}

int CrtCal::FindNabove(TH1F* h, int thresh){
	return h->Integral(h->FindBin(thresh),h->GetNbinsX()+1);
}

void CrtCal::PedCal(){

	if(fHasPedCal) return;

        float* statarr = new float[8];
        for(size_t i=0; i<8; i++){
        	statarr[i] = FLT_MAX;
        }

	for(size_t i=0; i<fHistos->size(); i++){
		const size_t chan = fChanMap[i];
		TH1F* h = (TH1F*)fHistos->at(i)->Clone();
		PedFit(h,chan,statarr,true);
		ParsePedStats(statarr, fPed[chan], fPedErr[chan], fPedNorm[chan], fPedNormErr[chan],
                               fPedSigma[chan], fPedSigmaErr[chan], fPedXsqr[chan], fPedNdf[chan] );
		delete h;
	}

	delete statarr;
	fHasPedCal = true;

	return;
}

void CrtCal::GainCal(){

	if(!fHasPedCal) {
		std::cout << "need ped cal before gain cal!" << std::endl;
		return;
	}
	if(fHasGainCal) return;

	float** statarr = new float*[15];
	for(size_t i=0; i<15; i++){
                statarr[i] = new float[5];
                for(size_t j=0; j<5; j++){
                        statarr[i][j] = FLT_MAX;
                }
        }

	for(size_t i=0; i<fHistos->size(); i++){
		const size_t chan = fChanMap[i];
		if(!fActive[chan]) continue;
		//std::cout << "channel " << chan << std::endl;
		TH1F* h = (TH1F*)fHistos->at(i)->Clone();
		//std::cout << "call to GainFit()" << std::endl;
        	GainFit(h,chan,statarr,true);
		//std::cout << "ParseGainStats()..." << std::endl;
		ParseGainStats(statarr, fGainXsqr[chan], fGainNdf[chan], fGain[chan], fGainErr[chan],
                            fGainPed[chan], fGainPedErr[chan], fNpeak[chan], fPeakXsqr[chan],
                            fPeakMean[chan], fPeakMeanErr[chan], fPeakNorm[chan], fPeakNormErr[chan],
                            fPeakSigma[chan], fPeakSigmaErr[chan], fPeakNdf[chan]);
		delete h;
	}

	delete statarr;
	fHasGainCal = true;
	return;
}

void CrtCal::PedFit(TH1F* h, size_t chan, float* statsarr, bool save=false){

	const int low = 1, high = 300;

	h->GetXaxis()->SetRangeUser(low,high);

        const int max_bin = h->GetMaximumBin();
        const int max_val = h->GetMaximum();
        const int max_adc = h->GetBinLowEdge(max_bin);

        //Define fitting function:
        TF1 *gausfit = new TF1("gausfit","[0]*exp(-0.5*((x-[1])/[2])^2)",max_adc-12,max_adc+12);

        //Set parameter names:
        gausfit->SetParName(0,"Constant");
        gausfit->SetParName(1,"Peak value");
        gausfit->SetParName(2,"sigma");

        //Initial guesses for the parameters:
        gausfit->SetParameter(0,max_val);
        gausfit->SetParLimits(0,0.5*max_val,1000*max_val);
        gausfit->SetParameter(1,max_adc);
        gausfit->SetParLimits(1,max_adc-20,max_adc+20);
        gausfit->SetParameter(2,50.0);//h->GetStdDev());
        gausfit->SetParLimits(2,1,50);

        TCanvas *c = new TCanvas();
        c->SetGrid();

        gStyle->Reset("Modern");
        if(h->GetMean()>(high-low)/2)
                gStyle->SetStatX(0.4);
        else
                gStyle->SetStatX(0.9);

        gStyle->SetStatY(0.89);
        gStyle->SetStatH(0.45);
        gStyle->SetStatW(0.28);
        h->Draw("e0samehist");
        //h->Draw("E0same");//hist");
        h->Fit(gausfit,"RMQ");

        float csq = gausfit->GetChisquare(); //chi-squared
        float ndf = gausfit->GetNDF(); //NDF
        float rcsq = csq/ndf; //reduced chi-squared

     if(rcsq>10.0)
        {
                //cout << "X^2/NDF > 10...recursive refit..." << endl;
                gausfit->SetRange(gausfit->GetParameter(1)-10,gausfit->GetParameter(1)+10);
                gausfit->SetParameter(0,gausfit->GetParameter(0));
                gausfit->SetParLimits(0,gausfit->GetParameter(0)-50,gausfit->GetParameter(0)+50);
                gausfit->SetParameter(1,gausfit->GetParameter(1));
                gausfit->SetParLimits(1,gausfit->GetParameter(1)-10,gausfit->GetParameter(1)+10);
                gausfit->SetParameter(2,gausfit->GetParameter(2));
                gausfit->SetParLimits(2,gausfit->GetParameter(2)-10,gausfit->GetParameter(2)+10);
                h->Fit(gausfit,"RMQ");
        }

        gausfit->Draw("same");

	//float* statsarr = new float[8];
        statsarr[0]=gausfit->GetParameter(0); //const
        statsarr[1]=gausfit->GetParError(0);  //const err
        statsarr[2]=gausfit->GetParameter(1); //mean
        statsarr[3]=gausfit->GetParError(1);  //mean err
        statsarr[4]=gausfit->GetParameter(2); //sigma
        statsarr[5]=gausfit->GetParError(2);  //sigma err
	statsarr[6]=gausfit->GetChisquare();
	statsarr[7]=gausfit->GetNDF(); 

	//std::cout << "   ped value = " << statarr[2] << std::endl;

        //Adding chi^2/ndf and fit parameters to stat box:
        gStyle->SetOptStat(1111);
        gStyle->SetOptFit(111);
        c->Update();

        csq = gausfit->GetChisquare();
        ndf = gausfit->GetNDF();
        rcsq = csq/ndf;
        if (rcsq>200.0) std::cout << "warning: possibly bad ped fit mac5 " << fMac5 << ", ch. "
                        << chan << " X^2/NDF=" << rcsq << " ADC" << std::endl;
        //save plot to file if desired
        if ( save )
        {
                int ctr = 1;
                TString suff=".png";
                TString fname = "./";//calPlotDir+"/";
                fname+=fMac5;
		fname+="_";
		fname+="chan"; //FormatMacString(mac)+ "_ch";
                if(chan<10) fname+="0";
                fname+=to_string(chan)+"_pedfit_"+to_string(ctr)+suff;

                while (!gSystem->AccessPathName(fname))
                {
                        fname.Remove(fname.Length()-suff.Length(),fname.Length());
                        fname.Remove(fname.Length()-1,fname.Length());
                        ctr++;
                        fname+=to_string(ctr)+suff;
                }

                //write image to file
                TImage *img = TImage::Create();
                img->FromPad(c);
                img->WriteImage(fname);

                delete img;
                delete gausfit;
                delete c;
        }//end save opt

        return;// statarr;

}

void CrtCal::GainFit(TH1F* h, size_t chan, float** statsarr, bool save){

	const float gainSeed = 55.0;
	const Double_t tSpectrumSigma = 1.75;
        const Double_t tSpectrumThreshold = 0.18;
	const size_t nPeakMax = 5;
	
	h->Rebin(4);
	h->GetXaxis()->SetRangeUser(350,700);
	h->SetStats(kFALSE);

	TCanvas* cspec = new TCanvas();
	h->Draw("e0hist");

        TSpectrum *s = new TSpectrum();
        //args=source histo, sigma of searched peaks, options, threshold
        int nPeak = s->Search(h,tSpectrumSigma,"",tSpectrumThreshold);
	//std::cout << "found " << nPeak << " peaks to fit" << std::endl;

        //int ctr = 2;

	//std::cout << "declaring in function arrays" << std::endl;
        float* peakXsqr     = new float[5];
        float* peakNdf      = new float[5];
        float* peakSigma    = new float[5];
        float* peakSigmaErr = new float[5];
	float* peakNorm     = new float[5];
	float* peakNormErr  = new float[5];
	float* peakMean     = new float[5];
	float* peakMeanErr  = new float[5];
        /*float ** statarr    = new float*[15];
	for(size_t i=0; i<15; i++){
		statarr[i] = new float[5];*/
		for(size_t j=0; j<5; j++){
			//statarr[i][j] = FLT_MAX;
			peakXsqr[j] = FLT_MAX;
			peakNdf[j] = FLT_MAX;
			peakSigma[j] = FLT_MAX;
			peakSigmaErr[j] = FLT_MAX;
			peakNorm[j] = FLT_MAX;
			peakNormErr[j] = FLT_MAX;
			peakMean[j] = FLT_MAX;
			peakMeanErr[j] = FLT_MAX;
		}
	//}

	//std::cout << "done." << std::endl;

        float *peds    = GetPed();
        float *pwidths = GetPedSigma();
	//std::cout << "retreived ped info" << std::endl;

        Double_t *peaks = s->GetPositionX(); //candidate photopeaks ADC position

	//std::cout << "sorting peaks" << std::endl;
        //Sort peaks
        double temp = 0.;
        int nchanges = 0;
        do
        {
		nchanges=0;
                for(int p=0;p<nPeak-1;p++)
                        if(peaks[p]>peaks[p+1])
                        {
                                temp=peaks[p];
                                peaks[p]=peaks[p+1];
                                peaks[p+1]=temp;
                                nchanges++;
                        }
        } while( nchanges != 0 );

	//std::cout << "done." << std::endl;

        //ascending list of peak number(x) vs. ADC value(y) from TSpectrum
        float x[10];
        float y[10];//, ey[10];

        //same list with "bad" peaks excluded (what eventually goes into the gain fit)
        float gx[11];
        float gy[11], gey[11];
        int peak_offset = round((peaks[0]-peds[chan])/gainSeed); //estimate peak number

       for (int j=0; j<nPeak&&j<10; j++)
        {
		if(peaks[j]<0) std::cout << "bad peak value: " << peaks[j] << std::endl;
                x[j] = j+peak_offset;
                y[j] = peaks[j];
        }

        int gg = 1; //index of passing peak array
        int nplow = 0; //no. peaks close to low hist edge
        int nphigh = 0; //no. peaks close to high hist edge

        //first peak in passing array is pedestal
        gx[0] = 0;
        gy[0] = peds[chan];
        gey[0] = pwidths[chan];

        //first and refit chi-squareds
        double chisqr = 0.;//, chisqr0;
       // int rnum = 0; //refit number

        //fit peaks about TSpectrum values
        for (int g=0 ; g<nPeak&&g<(int)nPeakMax; g++)
        {
                //initial gaus fit to peak from TSpectrum
                TF1 *gfit = new TF1("gfit","gaus",y[g]-20,y[g]+20);
                gfit->SetParameter(0,h->GetBinContent(h->FindBin(y[g])));//200);
                gfit->SetParameter(1,y[g]);
                gfit->SetParameter(2,12);
                gfit->SetParLimits(0,0,20000);
                gfit->SetParLimits(1,y[g]-15,y[g]+15);
                gfit->SetParLimits(2,8,40);
                //std::cout << "fitting gaussian to peak detected at " << y[g] << " ADC..." << std::endl;

                h->Fit(gfit,"MQLR"); //fit peak
                gfit->Draw("same");
		//std::cout << "done" << std::endl;
                //chisqr0 = gfit->GetChisquare()/gfit->GetNDF(); //get reduced chi-square from fit

                //check for peaks near the edges of the histogam
                if( y[g]<h->GetBinLowEdge(1)+15) nplow++;
                if( y[g]>h->GetBinLowEdge(h->GetNbinsX())-15) nphigh++;

                //ignore edge peaks and peaks with fit mean > 15 from peak
                if( y[g]>h->GetBinLowEdge(1)+15 && y[g]<h->GetBinLowEdge(h->GetNbinsX())-15&&abs(y[g] - gfit->GetParameter(1)) < 15)
                {
                        //skip false peaks (may need improvement - currently relies on init values)
                        if(g!=0 && g!=nPeak-1 //check it's not first or last peaks
                         && (y[g+1]-y[g]<gainSeed*0.7 || y[g]-y[g-1]<gainSeed*0.7))
                        { //check peak within 30% of expected gain w.r.t adj.
                                //std::cout << "determined peak " << g << " (" << y[g] << " ADC) is false peak. removing..." << std::endl;
                                for (int j=g; j<nPeak; j++) x[j] = x[j]-1;
                        }//overwrite current value, shifting higher values down by 1 index
                        //if not missed peak, could there have been a skipped peak?
                        else
                        {
                                if(g!=0&&g!=nPeak-1&&y[g+1]-y[g]<gainSeed*1.2&&y[g]-y[g-1]>gainSeed*1.5) //missed peak adjust
                                {
                                        for (int j=g; j<nPeak; j++) x[j] = x[j]+1;
                                        //std::cout << "missed peak detected. shifting peaks " << g << " and higher up by 1..." << std::endl;
                                }
                                //if last peak likely occuring after skipped peak
                                if(g==nPeak-1 && (y[g]-y[g-1])/gainSeed >1.5) {
                                        x[g]+=(int)((y[g]-y[g-1])/gainSeed);
                                        //std::cout << "missed peak(s) detected before last peak...shifting last peak" << std::endl;
                                }

                                peakXsqr[gg] = gfit->GetChisquare();
                                peakNdf[gg] = gfit->GetNDF();
				peakMean[gg] = gfit->GetParameter(1);
				peakMeanErr[gg] = gfit->GetParError(1);
				peakNorm[gg] = gfit->GetParameter(0);
				peakNormErr[gg] = gfit->GetParError(0);
                                peakSigma[gg] =  gfit->GetParameter(2);
                                peakSigmaErr[gg] = gfit->GetParError(2);

                                gx[gg] = x[g];
                                gy[gg] = gfit->GetParameter(1);
                                gey[gg] = sqrt(gx[gg]+gey[0]*gey[0]);//gfit->GetParError(1);
                                gg++;
                        }//end else not missed but was there skip
                }//end if not edge peak, not far from TSpectrum value
        }//end for peaks

	//std::cout << "done looping over peaks" << std::endl;

        if (nplow>0) for (int i=1; i<gg; i++) gx[i] = gx[i]-(nplow-1);

        //graph of adc(y) vs. photo-peak number (x)
        TGraphErrors* gr_mean = new TGraphErrors(gg,gx,gy,0,gey);

        //linear fit function
        TF1 *fit = new TF1("fit","[0] + [1]*x",gx[0]-0.25,gx[gg-1]+0.25);

        //name, initialize gain fit parameters
        fit->SetParName(1,"Gain");
        fit->SetParName(0, "Pedestal");
        fit->SetParameter(1,gainSeed);
        fit->SetParameter(0,peds[chan]);
        fit->SetParLimits(1,gainSeed-20,gainSeed+20);
        fit->SetParLimits(0,peds[chan]*0.8,peds[chan]*1.2);

	//std::cout << "gain fit..." << std::endl;
        //perform gain fit
        gr_mean->Fit(fit, "QR");

	//std::cout << "check chi-square..." << std::endl;
        //check if gain fit is bad according to chi-square
        if (fit->GetChisquare()/fit->GetNDF()>5.0)
        {
                //std::cout << "gain fit X^2 too large...shifting all peaks by 1" << std::endl;
                chisqr=fit->GetChisquare();
                for(int i=1; i<gg; i++) gx[i]+=1;
                gr_mean = new TGraphErrors(gg,gx,gy,0,gey);
                fit->SetRange(gx[0]-0.25,gx[gg-1]+0.25);
                gr_mean->Fit(fit, "QR");
                if (fit->GetChisquare()<chisqr) chisqr=fit->GetChisquare();
                else
                {
                        for(int i=1; i<gg; i++) gx[i]-=2;
                        gr_mean = new TGraphErrors(gg,gx,gy,0,gey);
                        fit->SetRange(gx[0]-0.25,gx[gg-1]+0.25);
                        gr_mean->Fit(fit, "QR");
                        if (fit->GetChisquare()<chisqr) chisqr=fit->GetChisquare();
                        else
                        {
                                for(int i=1; i<gg; i++) gx[i]+=1;
                                gr_mean = new TGraphErrors(gg,gx,gy,0,gey);
                                fit->SetRange(gx[0]-0.25,gx[gg-1]+0.25);
                                gr_mean->Fit(fit, "QR");
                        }
                }
        }

        string grtitle = "Mac5 "+std::to_string(fMac5)+", Ch. "+std::to_string(chan)+" Gain Fit";

	//std::cout << "style and range opts" << std::endl;

        gr_mean->SetTitle(grtitle.c_str());
        gr_mean->GetXaxis()->SetTitle("Peak #");
        gr_mean->GetYaxis()->SetTitle("ADC value");
        gr_mean->SetMarkerColor(4);
        gr_mean->SetMarkerStyle(20);
        gr_mean->SetFillColor(0);
	//std::cout << "set range..." << std::endl;
        gr_mean->GetXaxis()->SetRangeUser(gx[0]-0.5,gx[gg-1]+0.5);

	//std::cout << "done" << std::endl;
        //if (fit->GetChisquare()/fit->GetNDF()>30.0) cout << "Poor X^2! mac5 " << mac << "ch. " << chan << endl;

        gStyle->SetOptStat(0100);
        gStyle->SetOptFit(1111);

	//std::cout << "drawing..." << std::endl;

        TCanvas *c2 = new TCanvas();
        c2->cd();
        c2->SetGrid();

        gStyle->SetStatX(0.5);
        gStyle->SetStatY(0.9);
        gStyle->SetStatH(0.15);
        gStyle->SetStatW(0.2);

	//std::cout << "draw gr_mean..." << std::endl;
        gr_mean->Draw("ALP");
	//std::cout << "draw fit" << std::endl;
        fit->Draw("L,SAME");

	//std::cout << "filling statarr" << std::endl;

        statsarr[0][0] = fit->GetParameter(1); //gain
        statsarr[1][0] = fit->GetParError(1);  //gain error
        statsarr[2][0] = fit->GetParameter(0); //pedestal mean
        statsarr[3][0] = fit->GetParError(0);  //pedestal mean error
        statsarr[4][0] = fit->GetChisquare();  //X^2
        statsarr[5][0] = fit->GetNDF();        //NDF
        statsarr[6][0] = nPeak;
	for(size_t i=0; i<5; i++) {
	        statsarr[7][i]    = peakXsqr[i];
	        statsarr[8][i]    = peakMean[i];
	        statsarr[9][i]    = peakMeanErr[i];
	       	statsarr[10][i]   = peakNorm[i];
	        statsarr[11][i]   = peakNormErr[i];
	        statsarr[12][i]   = peakSigma[i];
	        statsarr[13][i]   = peakSigmaErr[i];
		statsarr[14][i]   = peakNdf[i];
	}
	
	delete peakXsqr;
	delete peakMean;
	delete peakMeanErr;
	delete peakNorm;
	delete peakNormErr;
	delete peakSigma;
	delete peakSigmaErr;
	delete peakNdf;

        if ( save )
        {
		//std::cout << "saving plots..."<< std::endl;
                int ctr=1, ctr2=1;
                TString suff = ".png";
                TString fname="./";
		TString fname2 = fname;
                fname+=fMac5;
		fname2+=fMac5;
		fname+="_ch";
		fname2+="_ch";
                if(chan<10) fname+="0";
		if(chan<10) fname2+="0";
		fname2+=chan;
		fname2+="_peak_fit_spec";
                fname+=to_string(chan)+"_fit-gain";
                fname+="_1"+suff;
		fname2+="_1"+suff;
                while (!gSystem->AccessPathName(fname))
                {
                        fname.Remove(fname.Length()-suff.Length(),fname.Length());
                        fname.Remove(fname.Length()-1,fname.Length());
                        ctr++;
                        fname+=to_string(ctr)+suff;
                }

                while (!gSystem->AccessPathName(fname2))
                {
                        fname2.Remove(fname2.Length()-suff.Length(),fname2.Length());
                        fname2.Remove(fname2.Length()-1,fname2.Length());
                        ctr2++;
                        fname2+=to_string(ctr2)+suff;
                }

                //write image to file
                TImage *img = TImage::Create();
                img->FromPad(c2);
                img->WriteImage(fname);

                //TImage *img2 = TImage::Create();
                img->FromPad(cspec);
                img->WriteImage(fname2);

                //deallocate memory
                delete img;
                delete c2;
		delete cspec;
		//delete img2;
                delete gr_mean;
                delete fit;
                delete s;
        }

	//std::cout << "done" << std::endl;

        return;// statarr;

}

//methods for parsing stats arrays
void CrtCal::ParsePedStats(const float* statarr, float& ped, float& pedErr, 
                           float& pedNorm, float& pedNormErr, float& pedSigma, float& pedSigmaErr, 
                           float& pedXsqr, short& pedNdf ){

        pedNorm     = statarr[0]; //const
        pedNormErr  = statarr[1]; //const err
        ped         = statarr[2]; //mean
        pedErr      = statarr[3]; //mean err
        pedSigma    = statarr[4]; //sigma
        pedSigmaErr = statarr[5]; //sigma err
	pedXsqr     = statarr[6];
        pedNdf      = (short)statarr[7];

	return;
}

void CrtCal::ParseGainStats(float** statarr, float& gainXsqr, short& gainNdf, float& gain, float& gainErr, 
                            float& gainPed, float& gainPedErr, short& nPeak, float* peakXsqr,
                            float* peakMean, float* peakMeanErr, float* peakNorm, float* peakNormErr,
                            float* peakSigma, float* peakSigmaErr, short* peakNdf) {

        gain         = statarr[0][0];
        gainErr      = statarr[1][0];
        gainPed      = statarr[2][0];
        gainPedErr   = statarr[3][0];
        gainXsqr     = statarr[4][0];
        gainNdf      = (short)statarr[5][0]; 
	nPeak        = statarr[6][0]; 
	for(size_t i=0; i<5; i++){
		//std::cout << "i="<< i << std::endl;
                peakXsqr[i]     = statarr[7][i];
                peakMean[i]     = statarr[8][i];
                peakMeanErr[i]  = statarr[9][i];
                peakNorm[i]     = statarr[10][i];
                peakNormErr[i]  = statarr[11][i];
                peakSigma[i]    = statarr[12][i];
                peakSigmaErr[i] = statarr[13][i];
		peakNdf[i]      = (short)statarr[14][i];
	}

	return;
}

//Lots of getters
//
bool*  CrtCal::GetActive() const {
	if(!fHasActive)
		return nullptr;

        return fActive;
}

float* CrtCal::GetPed() const{
	if(!fHasPedCal) 
		return nullptr;

	return fPed;
}

float*  CrtCal::GetPedErr() const{
        if(!fHasPedCal)
                return nullptr;

        return fPedErr;

}
float*  CrtCal::GetPedXsqr() const{
        if(!fHasPedCal)
                return nullptr;

        return fPedXsqr;
}
short*  CrtCal::GetPedNdf() const {
        if(!fHasPedCal)
                return nullptr;

        return fPedNdf;
}
float*  CrtCal::GetPedNorm() const{
        if(!fHasPedCal)
                return nullptr;

        return fPedNorm;
}
float*  CrtCal::GetPedNormErr() const{
        if(!fHasPedCal)
                return nullptr;

        return fPedNormErr;
}
float*  CrtCal::GetPedSigma() const{
        if(!fHasPedCal)
                return nullptr;

        return fPedSigma;
}
float*  CrtCal::GetPedSigmaErr() const{
        if(!fHasPedCal)
                return nullptr;

        return fPedSigmaErr;
}
float*  CrtCal::GetGain() const {
        if(!fHasGainCal)
                return nullptr;

        return fGain;
}
float*  CrtCal::GetGainErr() const {
        if(!fHasGainCal)
                return nullptr;

        return fGainErr;
}
float*  CrtCal::GetGainXsqr() const{
        if(!fHasGainCal)
                return nullptr;

        return fGainXsqr;
}
short* CrtCal::GetGainNdf() const{
        if(!fHasGainCal)
                return nullptr;

        return fGainNdf;
}

float*  CrtCal::GetGainPed() const{
        if(!fHasGainCal)
                return nullptr;

        return fGainPed;
}
float*  CrtCal::GetGainPedErr () const{
        if(!fHasGainCal)
                return nullptr;

        return fGainPedErr;
}
short*  CrtCal::GetNpeak() const{
        if(!fHasGainCal)
                return nullptr;

        return fNpeak;
}
int*  CrtCal::GetThreshADC() const{
	if(!fHasThresh)
		return nullptr;

        return fThreshADC;
}
float*  CrtCal::GetThreshPE() const{
        if(!fHasGainCal)
                return nullptr;

        return fThreshPE;
}
int*  CrtCal::GetNabove() const{
        if(!fHasThresh)
                return nullptr;

        return fNabove;
}
float** CrtCal::GetPeakNorm() const{
        if(!fHasGainCal)
                return nullptr;

        return fPeakNorm;
}
float** CrtCal::GetPeakNormErr() const{
        if(!fHasGainCal)
                return nullptr;

        return fPeakNormErr;
}
float** CrtCal::GetPeakSigma() const{
        if(!fHasGainCal)
                return nullptr;

        return fPeakSigma;
}
float** CrtCal::GetPeakSigmaErr() const{
        if(!fHasGainCal)
                return nullptr;

        return fPeakSigmaErr;
}
float** CrtCal::GetPeakMean() const{
        if(!fHasGainCal)
                return nullptr;

        return fPeakMean;
}
float** CrtCal::GetPeakMeanErr() const{
        if(!fHasGainCal)
                return nullptr;

        return fPeakMeanErr;
}
float** CrtCal::GetPeakXsqr() const{
        if(!fHasGainCal)
                return nullptr;

        return fPeakXsqr;
}
short** CrtCal::GetPeakNdf() const{
        if(!fHasGainCal)
                return nullptr;

        return fPeakNdf;
}

#endif
