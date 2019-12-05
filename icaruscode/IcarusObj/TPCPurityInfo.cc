#include "icaruscode/IcarusObj/TPCPurityInfo.hh"
#include <iostream>

anab::TPCPurityInfo::TPCPurityInfo()
{
  Run=0;
  Subrun=0;
  Event=0;
  Cryostat=999999;
  TPC=999999;
  Attenuation = -9999999;
}

void anab::TPCPurityInfo::Print()
{
  std::cout << "TPCPurityInfo:" << std::endl;
  std::cout << "\tRun,Subrun,Event: " << Run << "," << Subrun  << "," << Event;
  std::cout << "\n\tTPC: " << TPC;
  std::cout << "\n\tAttenuation = " << Attenuation;
  std::cout << std::endl;
}
