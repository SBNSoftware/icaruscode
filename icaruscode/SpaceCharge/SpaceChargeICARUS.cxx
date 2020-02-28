///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// SpaceChargeICARUS.cxx; brief implementation of class for storing/accessing space charge distortions for ICARUS
// Based largely on SpaceChargeSBND and SpaceChargeProtoDUNE
// rlazur@colostate.edu
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// C++ language includes
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <stdio.h>
// LArSoft includes
#include "icaruscode/SpaceCharge/SpaceChargeICARUS.h"

// Framework includes
#include "cetlib_except/exception.h"

spacecharge::SpaceChargeICARUS::SpaceChargeICARUS(fhicl::ParameterSet const& pset)
{
  Configure(pset);
}

bool spacecharge::SpaceChargeICARUS::Configure(fhicl::ParameterSet const& pset)
{
  fEnableSimSpatialSCE = pset.get<bool>("EnableSimSpatialSCE");
  fEnableSimEfieldSCE = pset.get<bool>("EnableSimEfieldSCE");
  //fEnableCorrSCE = pset.get<bool>("EnableCorrSCE");
  fEnableCalSpatialSCE = pset.get<bool>("EnableCalSpatialSCE");
  fEnableCalEfieldSCE = pset.get<bool>("EnableCalEfieldSCE");

  std::cout << "Configuring SpaceCharge..." << std::endl;

  if((fEnableSimSpatialSCE == true) || (fEnableSimEfieldSCE == true)){
    fRepresentationType = pset.get<std::string>("RepresentationType");
    fInputFilename = pset.get<std::string>("InputFilename");

    std::string fname;
    cet::search_path sp("FW_SEARCH_PATH");
    sp.find_file(fInputFilename, fname);

    std::unique_ptr<TFile> infile(new TFile(fname.c_str(), "READ"));
    if(!infile->IsOpen()){
      throw cet::exception("SpaceChargeICARUS") << "Could not find the space charge input file '" << fInputFilename << "'!\n";
    }

    if(fRepresentationType == "Voxelized_TH3"){
      std::cout << "begin loading voxelized TH3s..." << std::endl;

      //Load in histograms
      TH3F* hTrueFwdX = (TH3F*) infile->Get("TrueFwd_Displacement_X");
      TH3F* hTrueFwdY = (TH3F*) infile->Get("TrueFwd_Displacement_Y");
      TH3F* hTrueFwdZ = (TH3F*) infile->Get("TrueFwd_Displacement_Z");
      TH3F* hTrueBkwdX = (TH3F*) infile->Get("TrueBkwd_Displacement_X");
      TH3F* hTrueBkwdY = (TH3F*) infile->Get("TrueBkwd_Displacement_Y");
      TH3F* hTrueBkwdZ = (TH3F*) infile->Get("TrueBkwd_Displacement_Z");
      TH3F* hTrueEFieldX = (TH3F*) infile->Get("True_ElecField_X");
      TH3F* hTrueEFieldY = (TH3F*) infile->Get("True_ElecField_Y");
      TH3F* hTrueEFieldZ = (TH3F*) infile->Get("True_ElecField_Z");

      //https://root.cern.ch/doc/master/classTH1.html#a0367fe04ae8709fd4b82795d0a5462c3
      //Set hist directories so they can be referenced elsewhere
      //This needs to be done because they were read in from ext file
      //Note this is not a property of the TH3F, so does't survive copying
      hTrueFwdX->SetDirectory(0);
      hTrueFwdY->SetDirectory(0);
      hTrueFwdZ->SetDirectory(0);
      hTrueBkwdX->SetDirectory(0);
      hTrueBkwdY->SetDirectory(0);
      hTrueBkwdZ->SetDirectory(0);
      hTrueEFieldX->SetDirectory(0);
      hTrueEFieldY->SetDirectory(0);
      hTrueEFieldZ->SetDirectory(0);

      //SCEhistograms can be accessed globally in this script
      SCEhistograms = {hTrueFwdX, hTrueFwdY, hTrueFwdZ,
		       hTrueBkwdX, hTrueBkwdY, hTrueBkwdZ,
		       hTrueEFieldX, hTrueEFieldY, hTrueEFieldZ};


      std::cout << "...finished loading TH3s" << std::endl;
    }
    infile->Close();
  }
  if(fEnableCorrSCE == true){
    //keeping here (for now) for historic reasons
  }
  return true;
}

bool spacecharge::SpaceChargeICARUS::Update(uint64_t ts)
{
  if (ts == 0){return false;}
  return true;
}

// Whether or not to turn simulation of SCE on for spatial distortions
bool spacecharge::SpaceChargeICARUS::EnableSimSpatialSCE() const
{
  return fEnableSimSpatialSCE;
}

// Whether or not to turn simulation of SCE on for E-field distortions
bool spacecharge::SpaceChargeICARUS::EnableSimEfieldSCE() const
{
  return fEnableSimEfieldSCE;
}

// Return boolean indicating whether 
bool spacecharge::SpaceChargeICARUS::EnableCalSpatialSCE() const
{
  return fEnableCalSpatialSCE;
}

// Return boolean indicating whether or not to apply SCE corrections
bool spacecharge::SpaceChargeICARUS::EnableCalEfieldSCE() const
{
  return fEnableCalEfieldSCE;
}

/////////////////////////////////////////////////////////////////////////////
// BELOW ARE THE WORKHORSE FUNCTIONS
////////////////////////////////////////////////////////////////////////////

// Primary working method of service that provides position offsets
geo::Vector_t spacecharge::SpaceChargeICARUS::GetPosOffsets(geo::Point_t const& point) const
{
  std::vector<double> thePosOffsets;
  double xx=point.X(), yy=point.Y(), zz=point.Z();
  double corr=1.;

  if(fRepresentationType == "Voxelized_TH3"){
    //handle OOAV by projecting edge cases
    //also only have map for positive cryostat (assume symmetry)
    //need to invert coordinates for cryo0

    //in larsim, this is how the offsets are used in DriftElectronstoPlane_module
    // DriftDistance += -1.0 * thePosOffsets[0]
    // transversePos1/2 = xyz[1/2] + thePosOffsets[1/2]
    if(xx>0){
      corr=1.0;
    }else{
      corr=-1.0;
    }
    fixCoords(&xx, &yy, &zz); //bring into AV and x = abs(x)
    double offset_x=0., offset_y=0., offset_z=0.;
    offset_x = corr*SCEhistograms.at(0)->Interpolate(xx,yy,zz);
    offset_y = SCEhistograms.at(1)->Interpolate(xx,yy,zz);
    offset_z = SCEhistograms.at(2)->Interpolate(xx,yy,zz);
    thePosOffsets = {offset_x, offset_y, offset_z};
  }else{
    thePosOffsets.resize(3, 0.0);
  }

  return { thePosOffsets[0], thePosOffsets[1], thePosOffsets[2] };
}

// Returns the SCE correction at a specific point in the AV
  geo::Vector_t spacecharge::SpaceChargeICARUS::GetCalPosOffsets(geo::Point_t const& point, int const& TPCid) const
{
  std::vector<double> theCalPosOffsets;
  //make copies of const vars to modify
  double xx=point.X(), yy=point.Y(), zz=point.Z();
  int tpcid = TPCid;

  if(fRepresentationType == "Voxelized_TH3"){
    //handle OOAV by projecting edge cases
    //also only have map for positive cryostat (assume symmetry)
    //need to invert coordinates for cryo0
    double corr=1.;
    if(xx<0){
      corr=-1.0;
    }
    fixCoords(&xx, &yy, &zz, &tpcid); //bring into AV and x = abs(x)
    //handle the depositions that was reconstructed in the wrong TPC   
    //hard code in the cathode faces (got from dump_icarus_geometry.fcl)
    if ( tpcid == 0 && xx > 220.14 ) { xx = 220.14; }
    if ( tpcid == 1 && xx < 220.29 ) { xx = 220.29; }
    double offset_x=0., offset_y=0., offset_z=0.;
    offset_x = corr*SCEhistograms.at(3)->Interpolate(xx,yy,zz);
    offset_y = SCEhistograms.at(4)->Interpolate(xx,yy,zz);
    offset_z = SCEhistograms.at(5)->Interpolate(xx,yy,zz);
    theCalPosOffsets = {offset_x, offset_y, offset_z};
  }else{
    theCalPosOffsets.resize(3, 0.0);
  }
  
  return { theCalPosOffsets[0], theCalPosOffsets[1], theCalPosOffsets[2] };
}

geo::Vector_t spacecharge::SpaceChargeICARUS::GetCalPosOffsets(geo::Point_t const& point, geo::TPCID const& TPCid ) const
{
  return GetCalPosOffsets(point, TPCid.TPC);
}

// Primary working method of service that provides E field offsets
geo::Vector_t spacecharge::SpaceChargeICARUS::GetEfieldOffsets(geo::Point_t const& point) const
{
  //chiefly utilized by larsim, ISCalculationSeparate
  //the magnitude of the Efield is most important
  std::vector<double> theEfieldOffsets;
  double xx=point.X(), yy=point.Y(), zz=point.Z();
  double offset_x=0., offset_y=0., offset_z=0.;

  if(fRepresentationType == "Voxelized_TH3"){
    //handle OOAV by projecting edge cases
    //also only have map for positive cryostat (assume symmetry)
    fixCoords(&xx, &yy, &zz);
    offset_x = SCEhistograms.at(6)->Interpolate(xx, yy, zz);
    offset_y = SCEhistograms.at(7)->Interpolate(xx, yy, zz);
    offset_z = SCEhistograms.at(8)->Interpolate(xx, yy, zz);

    theEfieldOffsets = {offset_x, offset_y, offset_z};
  }
  return { theEfieldOffsets[0], theEfieldOffsets[1], theEfieldOffsets[2] };
}

void spacecharge::SpaceChargeICARUS::fixCoords(double* xx, double* yy, double* zz) const{
  //handle the edge cases by projecting SCE corrections onto boundaries
  *xx = abs(*xx);
  if(*xx<71.94){*xx=71.94;}
  if(*xx>368.49){*xx=368.489;}
  if(*yy<-181.86){*yy=-181.86;}
  if(*yy>134.96){*yy=134.96;}
  if(*zz<-894.951){*zz=-894.951;}
  if(*zz>894.951){*zz=894.9509;}
}

void spacecharge::SpaceChargeICARUS::fixCoords(double* xx, double* yy, double* zz, int* tpcid) const{
  //handle the edge cases by projecting SCE corrections onto boundaries
  //using tpcid to disambiguate hits that cross cathode due to SCE
  //need to do some fancy flipping of the TPC ids as well because the use of abs()
  if (xx < 0) { *tpcid = abs(*tpcid - 1); }
  *xx = abs(*xx);
  if(*xx<71.94){*xx=71.94;}
  if(*xx>368.49){*xx=368.489;}
  if(*yy<-181.86){*yy=-181.86;}
  if(*yy>134.96){*yy=134.96;}
  if(*zz<-894.951){*zz=-894.951;}
  if(*zz>894.951){*zz=894.9509;}
}
