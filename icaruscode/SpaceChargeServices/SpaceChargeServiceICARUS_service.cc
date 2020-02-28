/////////////////////////////////////////////////////////////////////////////
// SpaceChargeServiceICARUS_service.cc; brief implementation of class for storing/accessing space charge distortions for ICARUS
// Based on SpaceChargeServiceSBND_service.cc
// rlazur@fnal.gov
/////////////////////////////////////////////////////////////////////////////

// C++ language includes                                                                  
#include <iostream>

// LArSoft includes                                                                       
#include "icaruscode/SpaceChargeServices/SpaceChargeServiceICARUS.h"

// ROOT includes                                                                          
#include "TMath.h"

// Framework includes                                                                     
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

//-----------------------------------------------                                          
spacecharge::SpaceChargeServiceICARUS::SpaceChargeServiceICARUS(fhicl::ParameterSet const& pset, art::ActivityRegistry &reg)
{
  fProp.reset(new spacecharge::SpaceChargeICARUS(pset));
  
  reg.sPreBeginRun.watch(this, &SpaceChargeServiceICARUS::preBeginRun);
}

//----------------------------------------------                                          
void spacecharge::SpaceChargeServiceICARUS::preBeginRun(const art::Run& run)
{
  fProp->Update(run.id().run());
}

//------------------------------------------------                                        
void spacecharge::SpaceChargeServiceICARUS::reconfigure(fhicl::ParameterSet const& pset)
{
  fProp->Configure(pset);
  return;
}

//------------------------------------------------                                         
DEFINE_ART_SERVICE_INTERFACE_IMPL(spacecharge::SpaceChargeServiceICARUS, spacecharge::SpaceChargeService)
