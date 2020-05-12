#ifndef CRT_CAL_H
#define CRT_CAL_H

//c++ includes
#include <vector>
#include <map>
#include <string>
#include <climits>

//ROOT includes
#include<TH1F.h>
#include <TTree.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TF1.h>
#include <TImage.h>
#include <TSystem.h>
#include <TSpectrum.h>
#include <TGraphErrors.h>

namespace icarus {
 namespace crt {
  class CrtCal;
 }
}

using std::to_string;
using std::string;
using std::vector;
using std::map;

class icarus::crt::CrtCal {

  public:
	//CrtCal(string runName);
	CrtCal(const vector<TH1F*>* histos);
        ~CrtCal();

	void    Cal();
	void    PedCal();
	void    GainCal();
        void    GainFit(TH1F* h, size_t chan, float** statsarr, bool save);
        void    PedFit(TH1F* h, size_t chan, float* statsarr, bool save);
        bool    IsActive(TH1F*);

        void    ParsePedStats(const float* statarr, float& pedXsqr, float& ped, float& pedErr,
                           float& pedNorm, float& pedNormErr, float& pedSigma, float& pedSigmaErr,
			   short& pedNdf);

	void    ParseGainStats(float** statarr, float& gainXsqr, short& gainNdf, float& gain, float& gainErr,
                            float& gainPed, float& gainPedErr, short& nPeak, float* peakXsqr,
                            float* peakMean, float* peakMeanErr, float* peakNorm, float* peakNormErr,
                            float* peakSigma, float* peakSigmaErr,short* peakNdf);

	int  FindThreshADC(TH1F* h);
	float   FindThreshPE(TH1F* h);
	int  FindNabove(TH1F *h, int thresh);

	bool*   GetActive() const;
	float*  GetPed() const;
	float*  GetPedErr() const;
	float*  GetPedXsqr() const;
	short*  GetPedNdf() const;
	float*  GetPedNorm() const;
	float*  GetPedNormErr() const;
	float*  GetPedSigma() const;
	float*  GetPedSigmaErr() const;
	float*  GetGain() const;
	float*  GetGainErr() const;
	float*  GetGainXsqr() const;
	short*  GetGainNdf() const;
	float*  GetGainPed() const;
	float*  GetGainPedErr () const;
	short*  GetNpeak() const;
	int*    GetThreshADC() const;
	float*  GetThreshPE() const;
	int*    GetNabove() const;
	float** GetPeakNorm() const;
        float** GetPeakNormErr() const;
	float** GetPeakSigma() const;
        float** GetPeakSigmaErr() const;
        float** GetPeakMean() const;
        float** GetPeakMeanErr() const;
        float** GetPeakXsqr() const;
	short** GetPeakNdf() const;

  private:
        void IndexToMacChan();

        const vector<TH1F*>* fHistos;
        //string fRunName;
        //string fOutDir;
        map<uint8_t,uint8_t> fChanMap;

	bool fHasPedCal;
	bool fHasGainCal;
	bool fHasActive;
	bool fHasThresh;

	uint8_t  fMac5;
	bool*    fActive;
	float*   fGain;
	float*   fGainErr;
	float*   fGainXsqr;
	short*   fGainNdf;
	float*   fGainPed;
	float*   fGainPedErr;
	short*   fNpeak;
	float**  fPeakNorm;
	float**  fPeakNormErr;
	float**  fPeakSigma;
	float**  fPeakSigmaErr;
	float**  fPeakMean;
	float**  fPeakMeanErr;
	float**  fPeakXsqr;
	short**  fPeakNdf;
	float*   fPed;
	float*   fPedErr;
	float*   fPedXsqr;
        short*   fPedNdf;
	float*   fPedSigma;
	float*   fPedSigmaErr;
	float*   fPedNorm;
	float*   fPedNormErr;
	int*  fThreshADC;
	float*   fThreshPE;
	int*  fNabove;


};

#endif

