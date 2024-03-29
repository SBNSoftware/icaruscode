#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"

#include "CaloTools/IdealPmtResponse.h"
#include "CaloTools/FitBackgroundPhotons.h"

#include <stdio.h>
#include <fstream>
#include <iostream>


int main( int argc, char **argv ){


  // Read the CLI arguments

  std::string inputFilename;
  std::string  destinationFolder="calibrationdb";

  int debug =0;
  int startch = 0; 
  int endch = 359; 

  float rangeLow=0.3, rangeHigh=2.0;
  int modelPedestal=0;

  for ( int i=1; i<argc; i=i+2 )
  {
    if      ( std::string(argv[i]) == "-i" ) inputFilename = argv[i+1];
    else if ( std::string(argv[i]) == "-d" ) destinationFolder = argv[i+1];
    else if ( std::string(argv[i]) == "-l" ) rangeLow = atof(argv[i+1]);
    else if ( std::string(argv[i]) == "-h" ) rangeHigh = atof(argv[i+1]);
    else if ( std::string(argv[i]) == "-p" ) modelPedestal = atoi(argv[i+1]);
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

  IdealPmtResponse idealPmtResponseFunct(1, 4, hasExponential);
  FitBackgroundPhotons fitPMTResponse(6, idealPmtResponseFunct);
  fitPMTResponse.setFitRange( rangeLow, rangeHigh );

  if ( debug>0 && hasExponential )
    std::cout << "Use the IdealPMTResponse function between [" << rangeLow << " : " << rangeHigh << "] with exponential background model" << std::endl; 
  else if ( debug>0 && !hasExponential )
    std::cout << "Use the IdealPMTResponse function between [" << rangeLow << " : " << rangeHigh << "] without exponential background model" << std::endl; 


  // Here we open the input file and get the first timestamp to give a time reference
  TFile *tfile = new TFile(inputFilename.c_str(), "READ");

  if( !tfile->IsOpen() ){
    std::cout << inputFilename << " not found!" << std::endl;
    return 1;
  }

  int run = 0;
  int timestamp = 0;

  TTree *ttree = (TTree*)tfile->Get("bkgcalibration/event");
  ttree->SetBranchAddress("run", &run);
  ttree->SetBranchAddress("timestamp", &timestamp);
  ttree->GetEntry(0);


  std::cout << "Perform fit for run: " << run << ", Timestamp: " << timestamp << std::endl;

  // Here there is the output file
  char csvFilename[100];
  sprintf( csvFilename, "%s/backgroundphotons_run%d_%d.csv", destinationFolder.c_str(), run, timestamp );

  std::ofstream myfile; myfile.open( csvFilename );
  std::string line = "pmt,nentries,mu,emu,q,eq,sigma,esigma,amplitude,eamplitude,chi2,ndf,fitstatus\n";

  myfile << line;

  if( debug>0 )
    std::cout << line << std::endl;


  // Now we fit the histograms
  for( int pmt=startch; pmt<=endch; pmt++ )
  {

    char histname[100]; 
    sprintf(histname, "bkgcalibration/hintegral%d", pmt);

    TH1D *hintegral;

    // Just ignore the keys not found
    try{
      hintegral = (TH1D*)tfile->Get( histname );
    }
    catch (...) {
      std::cout << histname << "is a key not found! Ignore" << std::endl;
      continue;
    }

    if( !hintegral ){ continue; }
    
    int nentries = hintegral->GetEntries();
    if( nentries < 100 ){ continue; }

    // Do the fit 
    fitPMTResponse.fitHistogram( hintegral, "RQ0" );

    // Write the results to file
    std::string line = std::to_string(pmt) + "," + std::to_string(nentries) + "," ;

    for( int i=0; i<4; i++ )
      line += std::to_string(fitPMTResponse.getParameter(i)) + "," + std::to_string(fitPMTResponse.getParError(i)) + "," ;

    line += std::to_string(fitPMTResponse.getChi2()) + "," + std::to_string(fitPMTResponse.getNDF()) + "," ;
    line += std::to_string(fitPMTResponse.getFitStatus()) + "\n" ;

    myfile << line;

    if( debug>0 )
      std::cout << line << std::endl;

  }


  myfile.close();
  tfile->Close();

	return 0;

 } // End main
