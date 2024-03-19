#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"

#include "CaloTools/IdealPmtResponse.h"
#include "CaloTools/FitBackgroundPhotons.h"

#include <stdio.h>
#include <fstream>
#include <iostream>

// ---------------------------------------------------------------------

int main( int argc, char **argv ){

  // Read the CLI arguments
  std::string inputFilename;
  std::string  destinationFolder="calibrationdb";

  int debug =0;
  int startch = 0; 
  int endch = 359; 

  float fitRangeLow=0.3, fitRangeHigh=2.0;
  int modelPedestal=0;
  int nSum = 4;  

  float ampCutLow=1.22, ampCutHigh=4.;
  int cutOnAmplitude=0;

  for ( int i=1; i<argc; i=i+2 )
  {
    if      ( std::string(argv[i]) == "-i" ) inputFilename = argv[i+1];
    else if ( std::string(argv[i]) == "-d" ) destinationFolder = argv[i+1];
    else if ( std::string(argv[i]) == "-l" ) fitRangeLow = atof(argv[i+1]);
    else if ( std::string(argv[i]) == "-h" ) fitRangeHigh = atof(argv[i+1]);
    else if ( std::string(argv[i]) == "-al") ampCutLow = atof(argv[i+1]);
    else if ( std::string(argv[i]) == "-ah") ampCutHigh = atof(argv[i+1]);
    else if ( std::string(argv[i]) == "-p" ) modelPedestal = atoi(argv[i+1]);
    else if ( std::string(argv[i]) == "-n" ) nSum = atoi(argv[i+1]);
    else if ( std::string(argv[i]) == "-c" ) cutOnAmplitude = atoi(argv[i+1]);
    else if ( std::string(argv[i]) == "-v" ) debug = atoi(argv[i+1]);
    else if ( std::string(argv[i]) == "-s" ) startch = atoi(argv[i+1]);
    else if ( std::string(argv[i]) == "-e" ) endch = atoi(argv[i+1]);
    else {
      std::cout << "Unknown option \n" << argv[i+1] << std::endl;
      return 1;
    }
  }

  // The function to fit is initialized here
  bool hasExponential=false;  
  if( modelPedestal == 1 ) 
    hasExponential=true;

  IdealPmtResponse idealPmtResponseFunct(1, nSum, hasExponential);
  FitBackgroundPhotons fitPMTResponse(6, idealPmtResponseFunct);
  fitPMTResponse.setFitRange( fitRangeLow, fitRangeHigh );

  if ( debug>0 && hasExponential )
    std::cout << "Use the IdealPMTResponse function between [" << fitRangeLow << " : " << fitRangeHigh << "] with exponential background model" << std::endl; 
  else if ( debug>0 && !hasExponential )
    std::cout << "Use the IdealPMTResponse function between [" << fitRangeLow << " : " << fitRangeHigh << "] without exponential background model" << std::endl; 


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
  ttree->SetBranchAddress("run", &run);
  ttree->SetBranchAddress("timestamp", &timestamp);
  ttree->SetBranchAddress("channel", &channel);
  ttree->SetBranchAddress("amplitude", &amplitude);
  ttree->SetBranchAddress("integral", &integral);
  ttree->GetEntry(0);

  std::cout << "Perform fit for run: " << run << ", Timestamp: " << timestamp << std::endl;

  // Prepare the output file by setting the header
  char csvFilename[100];
  sprintf( csvFilename, "%s/backgroundphotons_run%d_%d.csv", destinationFolder.c_str(), run, timestamp );

  std::ofstream myfile; myfile.open( csvFilename );
  std::string line = "pmt,nentries,mu,emu,q,eq,sigma,esigma,amplitude,eamplitude,chi2,ndf,fitstatus\n";

  myfile << line;

  if( debug>0 )
    std::cout << line << std::endl;

  // if cutOnAmplitude is set to true, we need to build a new
  // histogram using only the ophits that survive the amplitude cut

  std::map<unsigned int, TH1D*> hintegral;
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

    if( cutOnAmplitude ){
    
	// same binning and limits, but reset
    	hintegral[pmt]->Reset();
    
	//change name to later save it in the file
    	sprintf(histname, "hintegralcut%d", pmt);
    	hintegral[pmt]->SetName(histname); 
    }  
  }

  // Build new integral histogram if needed 
  // we already reset them, just loop the tree
  if( cutOnAmplitude )
  {    
    for( int jj=0; jj<ttree->GetEntries(); jj++){
    
      ttree->GetEntry(jj);
      if( channel < startch || channel > endch ) continue; 
      if( amplitude > thresholds[channel] ) hintegral[channel]->Fill(integral);
    
    }
  }
   
  // Now we fit the histograms
  for( int pmt=startch; pmt<=endch; pmt++ )
  {

    // if we created a new one, save it for bookkeeping
    if( cutOnAmplitude ){
      tfile->cd("bkgcalibration");
      hintegral[pmt]->Write();
      tfile->cd();
    }
    
    if( debug>1 )
      std::cout << "Amplitude threshold for " << pmt << " is " << thresholds[pmt] << std::endl;
 
    // Do the fit 
    fitPMTResponse.fitHistogram( hintegral[pmt], "RQ0" );

    // Write the results to file
    int nentries = hintegral[pmt]->GetEntries();
    std::string line = std::to_string(pmt) + "," + std::to_string(nentries) + "," ;

    for( int i=0; i<4; i++ )
      line += std::to_string(fitPMTResponse.getParameter(i)) + "," + std::to_string(fitPMTResponse.getParError(i)) + "," ;

    line += std::to_string(fitPMTResponse.getChi2()) + "," + std::to_string(fitPMTResponse.getNDF()) + "," ;
    line += std::to_string(fitPMTResponse.getFitStatus()) + "\n" ;

    myfile << line;

    if( debug>0 ){
      std::cout << line << std::endl;
          
      if( pmt == startch ) fitPMTResponse.SaveToPdf(hintegral[pmt],"plots.pdf(");
      else if( pmt == endch ) fitPMTResponse.SaveToPdf(hintegral[pmt],"plots.pdf)");
      else fitPMTResponse.SaveToPdf(hintegral[pmt],"plots.pdf");
    }

  }

  myfile.close();
  tfile->Close();
  return 0;

 } // End main
