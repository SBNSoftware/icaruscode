#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLine.h"

#include "CaloTools/SPEChargeResponse.h"
#include "CaloTools/SPEChargeIdealResponse.h"
#include "CaloTools/SPEChargeLogisticResponse.h"
#include "CaloTools/FitBackgroundPhotons.h"

#include <stdio.h>
#include <fstream>
#include <iostream>

// ---------------------------------------------------------------------
/*  This macro computes the fits to the PMT pulse charge distribution
 *  to extract the gain for each channel. It takes as input an output
 *  file produced by PMTBackgroundphotonsCalibration_module.cc
 *
 *  The fitting itself is managed by CaloTools/FitBackgroundPhotons.cc
 *  Two SPE charge response functions are currently available (and more
 *  can be defined by deriving from CaloTools/SPEChargeResponse.h).
 *  All fitting options can be set via CLI arguments.
 *
 *  List of available arguments:
 *  -i inputFilename 
 *  -d destinationFolder (default: calibrationdb)
 *  -w writeOutut (default: 1)
 *  -v verbosity (default: 0)
 *  -s start channel (default: 0)
 *  -e end channel (default: 0)
 *  -r useIdealResponse (default: 0);     
 *  -p useNoisePedestal, only for ideal response (default: 0)
 *  -c cutOnAmplitude, always true for logistic response (default: 1)
 *  -n number of poisson-gaussians in sum (default: 2)
 *  -l fitRangeLow (default: 0.0)
 *  -h fitRangeHigh (default: 3.0)
 *  -al ampCutLowLimit (default: 1.22)
 *  -ah ampCutHighLimit (default: 4.0) 
 *  -qmin gain low fit limit (default: 0.15)
 *  -qmax gain high fit limit (default: 1.0)
 *
 *  The default settings use the amplitude distribution to find a
 *  selection thresholds between 1.22 mV and 4. mV to remove the noise
 *  pedestal. The ratio between the charge distribution before and after
 *  the cut is used to determine the parameters of the logistic function.
 *  The fit function is then an ideal response (sum of two poisson-gaussians)
 *  multiplied by the pre-determined logistic function.
 *
 *  The output is saved as .csv file in the destination folder.
 *  If debug>0, it also produces a .pdf file to visually inspect all fits.
*/
// ---------------------------------------------------------------------

void computeLogisticFunction( TH1D* hist, TH1D* histcut, double &x0, double &k ){

  TH1D *hratio = (TH1D*)hist->Clone("hratio");
  hratio->Reset();

  for(int ibin=1; ibin<=hratio->GetNbinsX(); ibin++)
  {  
    int total = hist->GetBinContent(ibin);
    int cut = histcut->GetBinContent(ibin);
    if( total < 1 ) continue;
    hratio->SetBinContent(ibin, float(cut)/total );
  }
 
  double lw = 0.;
  double hg = x0+(x0-lw);

  TF1 *fitf = new TF1( "logistic", "1./(1.+TMath::Exp(-[1]*(x-[0])))", lw, hg);
  fitf->SetParameters(x0,k);

  int status = hratio->Fit("logistic", "RQ0", "", lw, hg);
  x0 = fitf->GetParameter(0);
  k = fitf->GetParameter(1);

  if( status < 0 ) std::cout << "Logistic function fit FAILED for " << hist->GetTitle() << std::endl;

}

void saveToPdf(TH1D* hist, std::string name, double threshold){

  char cname[200];
  sprintf( cname, "c_%s", hist->GetName() );
  TCanvas *c = new TCanvas(cname,"fit",800.,600.); 	
        
  hist->Draw();
  // add line that shows threshold
  TLine *l = new TLine(threshold,0.,threshold,hist->GetMaximum()); 
  l->SetLineColor(kRed);
  l->SetLineStyle(kDashed);
  l->Draw("SAME");

  // add canvas to pdf
  c->Print(name.c_str(),"pdf");
}

// ---------------------------------------------------------------------

int main( int argc, char **argv ){

  // Read the CLI arguments
  std::string inputFilename;
  std::string  destinationFolder="calibrationdb";
  
  int debug =0;
  int startch = 0; 
  int endch = 359; 
  int writeOutput = 1;

  int useIdealResponse = 0;
  int usePedestal=0;
  int cutOnAmplitude=1;
  float fitRangeLow=0., fitRangeHigh=3.0;
  int nSum = 2;  
  float ampCutLow=1.22, ampCutHigh=4.;
  float qmin = 0.15, qmax = 1.0;

  for ( int i=1; i<argc; i=i+2 )
  {
    // scripts parameters
    if      ( std::string(argv[i]) == "-i" ) inputFilename = argv[i+1];
    else if ( std::string(argv[i]) == "-d" ) destinationFolder = argv[i+1];
    else if ( std::string(argv[i]) == "-v" ) debug = atoi(argv[i+1]);
    else if ( std::string(argv[i]) == "-w" ) writeOutput = atoi(argv[i+1]);
    // channel selection
    else if ( std::string(argv[i]) == "-s" ) startch = atoi(argv[i+1]);
    else if ( std::string(argv[i]) == "-e" ) endch = atoi(argv[i+1]);
    // fitting options
    else if ( std::string(argv[i]) == "-p" ) usePedestal = atoi(argv[i+1]);
    else if ( std::string(argv[i]) == "-r" ) useIdealResponse = atoi(argv[i+1]);     
    else if ( std::string(argv[i]) == "-c" ) cutOnAmplitude = atoi(argv[i+1]);
    else if ( std::string(argv[i]) == "-n" ) nSum = atoi(argv[i+1]);
    else if ( std::string(argv[i]) == "-l" ) fitRangeLow = atof(argv[i+1]);
    else if ( std::string(argv[i]) == "-h" ) fitRangeHigh = atof(argv[i+1]);
    else if ( std::string(argv[i]) == "-al") ampCutLow = atof(argv[i+1]);
    else if ( std::string(argv[i]) == "-ah") ampCutHigh = atof(argv[i+1]);
    else if ( std::string(argv[i]) == "-qmin") qmin = atof(argv[i+1]);
    else if ( std::string(argv[i]) == "-qmax") qmax = atof(argv[i+1]);
    else {
      std::cout << "Unknown option \n" << argv[i+1] << std::endl;
      return 1;
    }
  }
  
  // The function to fit is initialized here according to the inputs
  // SPEChargeIdealResponde (w/ or w/o pedestal) or SPEChargeLogisticResponse
  // in both cases nparameters=6, but the last two are used differently
  
  int nparameters = 6;
  pmtcalo::SPEChargeResponse *speResponseFunct = nullptr;

  bool hasExponential=false;  
  if( usePedestal == 1 ) hasExponential=true;
  
  if( useIdealResponse == 1 ){
    speResponseFunct = new pmtcalo::SPEChargeIdealResponse(1, nSum, hasExponential);
  } else {
    speResponseFunct = new pmtcalo::SPEChargeLogisticResponse(1, nSum, hasExponential);
    cutOnAmplitude = 1; //force it
  }
 
  // preparing fitting manager
  pmtcalo::FitBackgroundPhotons fitPMTResponse(nparameters, speResponseFunct);
  fitPMTResponse.setFitRange( fitRangeLow, fitRangeHigh );

  if ( debug>0 ){
  
    if ( useIdealResponse == 1 ){
      std::cout << "Use the SPEChargeIdealResponse function (nSum=" << nSum << ") between [" << fitRangeLow << " : " << fitRangeHigh << "]";
      if(hasExponential) std::cout<< " with exponential background model" << std::endl;
      else std::cout << " without exponential background model" << std::endl;  
    
    } else 
      std::cout << "Use the SPEChargeLogisticResponse function (nSum=" << nSum << ") between [" << fitRangeLow << " : " << fitRangeHigh << "]" << std::endl;
  
    if( cutOnAmplitude == 1 ) 
      std::cout << "Selection using amplitude cut between [" << ampCutLow << " : " << ampCutHigh << "]" << std::endl;

   }

  // Here we open the input file and get the first timestamp to give a time reference
  // We also prepare to extract the ophits from the tree (if needed later)
  
  TFile *tfile = new TFile(inputFilename.c_str(), "UPDATE");

  if( !tfile->IsOpen() ){
    std::cout << inputFilename << " not found!" << std::endl;
    return 1;
  }

  int run = 0;
  int timestamp = 0;
  int channel = -1;
  double amplitude = 0;
  double integral = 0;

  TTree *ttree = (TTree*)tfile->Get("bkgcalibration/ophits");
  if( !ttree ){
    std::cout << "TTree 'ophits' not found!" << std::endl;
    return 1;
  }

  ttree->SetBranchAddress("run", &run);
  ttree->SetBranchAddress("timestamp", &timestamp);
  ttree->SetBranchAddress("channel", &channel);
  ttree->SetBranchAddress("amplitude", &amplitude);
  ttree->SetBranchAddress("integral", &integral);
  ttree->GetEntry(0);
      
  if( debug>0 ) std::cout << "Run " << run << " timestamp " << timestamp << std::endl;

  // if cutOnAmplitude is set to true, we need to build a new
  // histogram using only the ophits that survive the amplitude cut

  std::map<unsigned int, TH1D*> hintegral;
  std::map<unsigned int, TH1D*> hintegralcut;
  std::map<unsigned int, double> thresholds;

  for( int pmt=startch; pmt<=endch; pmt++ )
  {

    char histname[100]; 
    char ampname[100]; 
    sprintf(histname, "bkgcalibration/hintegral%d", pmt);
    sprintf(ampname, "bkgcalibration/hamplitude%d", pmt);
    TH1D *hamplitude;

    // Just ignore the keys not found
    try{
      hintegral[pmt] = (TH1D*)tfile->Get( histname );
      hamplitude = (TH1D*)tfile->Get( ampname );
    }
    catch (...) {
      std::cout << histname << " or " << ampname << " is a key not found! Ignore" << std::endl;
      continue;
    }

    // skip if empty
    if( !hintegral[pmt] || !hamplitude ){ continue; }
    int nentries = hintegral[pmt]->GetEntries();
    if( nentries < 100 ){ continue; }
      
    // Find the threshold for the amplituce cut: look for the minimum of the
    // amplitude histogram within tunable boundaries
   
    hamplitude->GetXaxis()->SetRangeUser(ampCutLow,ampCutHigh);
    int minbin = hamplitude->GetMinimumBin();
    thresholds[pmt] = hamplitude->GetXaxis()->GetBinCenter(minbin);
    if( debug>1 ) std::cout << "Amplitude threshold for " << pmt << " is " << thresholds[pmt] << std::endl;

    if( cutOnAmplitude ){
    
	// same binning and limits, but reset
        hintegralcut[pmt] = (TH1D*)hintegral[pmt]->Clone("hcut");
        hintegralcut[pmt]->Reset();
	
        //change name to later save it in the file
    	sprintf(histname, "hintegralcut%d", pmt);
    	hintegralcut[pmt]->SetName(histname); 
    
        if( debug>0 ){
          std::string fname = "run" + std::to_string(run) + "_thresholds.pdf";
          std::string sfname = fname + "(";
          std::string efname = fname + ")";
          if( pmt == startch ) saveToPdf(hamplitude,sfname,thresholds[pmt]);
          else if( pmt == endch ) saveToPdf(hamplitude,efname,thresholds[pmt]);
          else saveToPdf(hamplitude,fname,thresholds[pmt]);
        }
    }    
  }

  // Build new integral histogram with cut (if needed) 
  // we already reset them, just loop the tree
  if( cutOnAmplitude )
  {    
    for( int jj=0; jj<ttree->GetEntries(); jj++){
    
      ttree->GetEntry(jj);
      if( channel < startch || channel > endch ) continue; 
      if( amplitude > thresholds[channel] ) hintegralcut[channel]->Fill(integral);
    
    }
  }
  
  // Prepare the output file by setting the header
  char csvFilename[100];
  sprintf( csvFilename, "%s/backgroundphotons_run%d_%d.csv", destinationFolder.c_str(), run, timestamp );

  std::ofstream myfile; 
  std::string line = "pmt,nentries,mu,emu,q,eq,sigma,esigma,amplitude,eamplitude,chi2,ndf,fitstatus\n";
  
  if ( writeOutput ){
    myfile.open( csvFilename );
    myfile << line;
  }
  if( debug>0 ) std::cout << line << std::endl;
   
  // Now we fit the histograms
  for( int pmt=startch; pmt<=endch; pmt++ )
  {

    // Set the starting values for the parameters
    // the last two will need to be FIXED later if using the logistic function
    std::vector<double> params = {  0.1,(qmax+qmin)/2.,0.3,hintegral[pmt]->Integral()*0.2,hintegral[pmt]->Integral()*0.8,-3. };
    fitPMTResponse.setFitParameters( params, qmin, qmax );

    // if we created a new histogram for the fit, save it for book-keeping
    // find the logistic function parameters (if requested) and update the fitting function
    // then do the fit

    int nentries = -1;

    if( cutOnAmplitude ){

      tfile->cd("bkgcalibration");
      if( writeOutput ) hintegralcut[pmt]->Write();
      tfile->cd();

      double x0 = 2., k = 1.;
      
      if( useIdealResponse == 0 ) {  // fit a logitstic function on the ratio before/after cut
        computeLogisticFunction( hintegral[pmt], hintegralcut[pmt], x0, k );
        fitPMTResponse.fixFitParameter(4, x0);
        fitPMTResponse.fixFitParameter(5, k);
      }
    
      // Do the fit 
      fitPMTResponse.fitHistogram( hintegralcut[pmt], "RQ0" );
      nentries = hintegralcut[pmt]->GetEntries();

    } else {
     
      // if not threshold on amplitude is needed, jut fit the original histograms
      // note that it is not possible to use the logistic function for the fit in this case

      fitPMTResponse.fitHistogram( hintegral[pmt], "RQ0" );
      nentries = hintegral[pmt]->GetEntries();

    }
   
    // Write the results to file
    std::string line = std::to_string(pmt) + "," + std::to_string(nentries) + "," ;

    for( int i=0; i<4; i++ )
      line += std::to_string(fitPMTResponse.getParameter(i)) + "," + std::to_string(fitPMTResponse.getParError(i)) + "," ;

    line += std::to_string(fitPMTResponse.getChi2()) + "," + std::to_string(fitPMTResponse.getNDF()) + "," ;
    line += std::to_string(fitPMTResponse.getFitStatus()) + "\n" ;

    if( writeOutput ) myfile << line;

    if( debug>0 ){
      std::cout << line << std::endl;
      std::string fname = "run" + std::to_string(run) + "_plots.pdf";
      std::string sfname = fname + "(";
      std::string efname = fname + ")";
      if( pmt == startch ) fitPMTResponse.SaveToPdf(hintegral[pmt],hintegralcut[pmt],sfname,cutOnAmplitude);
      else if( pmt == endch ) fitPMTResponse.SaveToPdf(hintegral[pmt],hintegralcut[pmt],efname,cutOnAmplitude);
      else fitPMTResponse.SaveToPdf(hintegral[pmt],hintegralcut[pmt],fname,cutOnAmplitude);
    }

  }

  myfile.close();
  tfile->Close();
  return 0;

 } // End main
