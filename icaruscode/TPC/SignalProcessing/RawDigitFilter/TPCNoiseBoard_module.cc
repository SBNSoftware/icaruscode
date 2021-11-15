////////////////////////////////////////////////////////////////////////
// Class:       TPCNoiseBoard
// Plugin Type: producer (art v3_05_01)
// File:        TPCNoiseBoard_module.cc
//
// Generated at Thu Jan 21 16:22:15 2021 by Justin Mueller using cetskelgen
// from cetlib version v3_10_00.
////////////////////////////////////////////////////////////////////////

// LArSoft Includes
#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalService.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalProvider.h"
#include "cetlib/cpu_timer.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib_except/exception.h"

// icaruscode Includes
#include "icaruscode/TPC/SignalProcessing/RawDigitFilter/Algorithms/RawDigitCharacterizationAlg.h"
#include "icaruscode/Decode/ChannelMapping/IICARUSChannelMap.h"
// icarus_signal_processing Includes
#include "icarus_signal_processing/Filters/ICARUSFFT.h"

// ROOT Includes
#include "TTree.h"
#include "TFile.h"
// C++ Includes
#include <map>
#include <vector>
#include <tuple>
#include <algorithm>
#include <iostream>
#include <string>
#include <cmath>
#include <iostream>
#include <fstream>
#include <memory>

#include <numeric> // std::accumulate

namespace TPCNoiseBoard {
  class TPCNoiseBoard;
}


class TPCNoiseBoard::TPCNoiseBoard : public art::EDAnalyzer {
public:
  explicit TPCNoiseBoard(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  //TPCNoiseBoard(TPCNoiseBoard const&) = delete;
  //TPCNoiseBoard(TPCNoiseBoard&&) = delete;
  //TPCNoiseBoard& operator=(TPCNoiseBoard const&) = delete;
  //TPCNoiseBoard& operator=(TPCNoiseBoard&&) = delete;

  void analyze(const art::Event& e);
  void reconfigure(fhicl::ParameterSet const& pset);
  void beginJob();
//  void beginRun();
  void endJob();
 // void endRun();

private:
  // Various FHiCL parameters.
  std::string fRawDigitModuleLabel;
  std::string fRawDigitProcess;
  std::string fRawInstance;
  std::string fIntInstance;
std::string fCohInstance;
std::string fHistoFileName;


  // FFT calculation.
  using FFTPointer = std::unique_ptr<icarus_signal_processing::ICARUSFFT<double>>;
  FFTPointer fFFT;
  int NumberTimeSamples;

  // FFT variables.
  std::vector< std::vector<float> > fRawPower;
  std::vector< std::vector<float> > fIntPower;
std::vector< std::vector<float> > fCohPower;

int ct;

std::vector<TH1D*> fRawPowerHistos;
std::vector<TH1D*> fCohPowerHistos;
std::vector<TH1D*> fIntPowerHistos;
std::vector<TH1D*> fRawRMSHistos;
std::vector<TH1D*> fCohRMSHistos;
std::vector<TH1D*> fIntRMSHistos;
std::vector<TH1D*> fMediaHistos;


  // The variables that will go into the n-tuple.
  int fEvent;
  int fRun;
  int fSubRun;
  std::vector<float> fPed;
  std::vector<float> fRawMean;
  std::vector<double> fRawRMS;
  std::vector<double> fRawRMSTrim;
  std::vector<float> fIntMean;
  std::vector<double> fIntRMS;
  std::vector<double> fIntRMSTrim;
  std::vector<float> fCohMean;
  std::vector<double> fCohRMS;
  std::vector<double> fCohRMSTrim;
  std::vector<unsigned short int> fChannel;

  // The variables that will go into the power n-tuple.
  //std::vector< std::vector<float> > fPower;
  unsigned int NEvents;
 
  // The output trees.
  TTree* fNoiseTree;
 /* TTree* fRawPowerTree;
  TTree* fIntrinsicPowerTree;
 TTree* fCoherentPowerTree;*/

    const icarusDB::IICARUSChannelMap*      fChannelMap;
};


TPCNoiseBoard::TPCNoiseBoard::TPCNoiseBoard(fhicl::ParameterSet const& p)
  : EDAnalyzer(p),
    NEvents(0),
 fChannelMap(art::ServiceHandle<icarusDB::IICARUSChannelMap const>{}.get())
{
  // Retrieve proper number of time samples.
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob();
  NumberTimeSamples = detProp.NumberTimeSamples();

  std::cout << "Number of time samples: " << NumberTimeSamples << std::endl;
ct=0;

  fRawPower.resize(1000);
  for(auto &it : fRawPower)
    {
      it.resize(NumberTimeSamples);
    }
std::cout << " after resizing " << std::endl;
  fIntPower.resize(1000);
  for(auto &it : fIntPower)
    {
      it.resize(NumberTimeSamples);
    }
std::cout << " after resizing " << std::endl;
 fCohPower.resize(1000);
  for(auto &it : fCohPower)
    {
      it.resize(NumberTimeSamples);
    }
std::cout << " after intrinsic resizing " << std::endl;

  fFFT = std::make_unique<icarus_signal_processing::ICARUSFFT<double>>(NumberTimeSamples);
  this->reconfigure(p);

 double freqBin=0.6103515625;

int maxBoard=1000;

for(int jb=0;jb<maxBoard;jb++) {
std::string rawPowerHistoName="rawFFT"+std::to_string(jb);
std::string rawRMSHistoName="rawRMS"+std::to_string(jb);
std::string cohRMSHistoName="cohRMS"+std::to_string(jb);
std::string intRMSHistoName="intRMS"+std::to_string(jb);
TH1D* fRawPowerHisto=new TH1D(rawPowerHistoName.c_str(),rawPowerHistoName.c_str(),2048,0.,2048*freqBin);
fRawPowerHistos.push_back(fRawPowerHisto);
std::string cohPowerHistoName="cohFFT"+std::to_string(jb);
TH1D* fCohPowerHisto=new TH1D(cohPowerHistoName.c_str(),cohPowerHistoName.c_str(),2048,0.,2048*freqBin);
fCohPowerHistos.push_back(fCohPowerHisto);
std::string intPowerHistoName="intFFT"+std::to_string(jb);
TH1D* fIntPowerHisto=new TH1D(intPowerHistoName.c_str(),intPowerHistoName.c_str(),2048,0.,2048*freqBin);
fIntPowerHistos.push_back(fIntPowerHisto);
TH1D* fRawRMSHisto=new TH1D(rawRMSHistoName.c_str(),rawRMSHistoName.c_str(),100,0.,30.);
fRawRMSHistos.push_back(fRawRMSHisto);
fRawRMS.push_back(0);
TH1D* fCohRMSHisto=new TH1D(cohRMSHistoName.c_str(),cohRMSHistoName.c_str(),100,0.,30.);
fCohRMSHistos.push_back(fCohRMSHisto);
fCohRMS.push_back(0);
TH1D* fIntRMSHisto=new TH1D(intRMSHistoName.c_str(),intRMSHistoName.c_str(),100,0.,30.);
fIntRMSHistos.push_back(fIntRMSHisto);
fIntRMS.push_back(0);
}


  fRawRMSTrim.clear();

  fIntMean.clear();

  fIntRMSTrim.clear();
    fChannelMap = art::ServiceHandle<icarusDB::IICARUSChannelMap const>{}.get();

std::cout << " end constructor " << std::endl;
}

void TPCNoiseBoard::TPCNoiseBoard::analyze(const art::Event& e)
{
std::cout << " begin analyze " << std::endl;
   art::ServiceHandle<geo::Geometry> geom;

  // Clear vectors before filling for this event.
  fChannel.clear();
std::cout << " after clearing " << std::endl;
  fPed.clear();
  fRawMean.clear();

std::cout << " after clearing " << std::endl;
  
  // Fill event level variables.
  fEvent = e.id().event();
std::cout << " event " << fEvent << std::endl;
  fRun = e.id().run();
std::cout << " run " << fRun << std::endl;
  fSubRun = e.id().subRun();
  std::cout << " subrun " << fSubRun << std::endl;
  std::cout << "Processing event " << NEvents+1 << " for TPC Noise Analysis " << " Run " << fRun << ", " << "Event " << fEvent << "." << std::endl;

  
  ///////////////////////////
  // "Raw" RawDigits.
  ///////////////////////////

  art::Handle< std::vector<raw::RawDigit> > RawDigitHandle;
  e.getByLabel(fRawDigitModuleLabel, fRawInstance, fRawDigitProcess, RawDigitHandle);

  std::cout << " raw instance " << fRawInstance << std::endl;
 /* std::vector<int,std::vector<short>* > RawPowers;
int maxBoard=1000;
for(int jm=0;jm<maxBoard;jm++) {
 std::vector<short>* v;
v->clear();
RawPowers.push_back(v);
}

*/
  for(const auto& RawDigit : *RawDigitHandle)
    {
std::vector<geo::WireID> widVec = geom->ChannelToWire(RawDigit.Channel());

 //       size_t                   plane  = widVec[0].Plane;
 //size_t                   wire  = widVec[0].Wire;
        size_t                   cryo  = widVec[0].Cryostat;
 size_t                   tpc  = widVec[0].TPC;
std::cout << "channel " << RawDigit.Channel() << " cryo " << cryo << " tpc " << tpc << std::endl;
    const icarusDB::TPCReadoutBoardToChannelMap& readoutBoardToChannelMap = fChannelMap->getReadoutBoardToChannelMap();
int boardFound=-1;
 for(const auto& boardPair : readoutBoardToChannelMap)
    {
 // A little song and dance to make sure this board is in our TPC group
        // What we need is a proper mapping but because we can have channels which are not valid for the geometry
        // service we have a little thing we have to go through here...
        std::vector<geo::WireID> wireIDVec;

        for(const auto& channelPair : boardPair.second.second)
        {
if(channelPair.first==RawDigit.Channel())
boardFound=boardPair.first;
        }
}
std::cout << " board found " << boardFound << std::endl;

//std::vector<short>* fRawPower=RawPowers[boardFound];
 std::vector<short> RawADC;
      // Grab raw waveform, ensuring that the size is set appropriately.
      unsigned int DataSize = RawDigit.Samples();

      RawADC.resize(DataSize);
      raw::Uncompress(RawDigit.ADCs(), RawADC, RawDigit.Compression());

      // We need a sorted waveform (by absolute value) for the truncated RMS and median calculation.
      std::vector<short> SortedADC(RawADC);
      std::sort(SortedADC.begin(),SortedADC.end(),[](const auto& left, const auto& right){return std::fabs(left) < std::fabs(right);});
      float median(SortedADC.at(SortedADC.size()/2));

      // Calculate mean values.
      float mean(float(std::accumulate(SortedADC.begin(),SortedADC.end(),0))/float(SortedADC.size()));





      // Remove pedestal of waveform.
      std::vector<float> ADCLessPed;
      ADCLessPed.resize(SortedADC.size());
      std::transform(SortedADC.begin(),SortedADC.end(),ADCLessPed.begin(),std::bind(std::minus<float>(),std::placeholders::_1,median));

      // Calculate full RMS.
      float rms(std::sqrt(std::inner_product(ADCLessPed.begin(), ADCLessPed.end(), ADCLessPed.begin(), 0.) / float(ADCLessPed.size())));

      // Calculate the truncated RMS.
      unsigned int MinBins((1.0 - 0.01)*ADCLessPed.size());
      //unsigned int BinsToKeep;
      //for(BinsToKeep = 0; BinsToKeep < ADCLessPed.size(); ++BinsToKeep)
      //{
      //  if(std::fabs(ADCLessPed.at(BinsToKeep)) >= 3*rms) break;
      //}
      //float truncrms(std::sqrt(std::inner_product(ADCLessPed.begin(), ADCLessPed.begin() + BinsToKeep, ADCLessPed.begin(), 0.) / float(BinsToKeep)));
      float truncrms(std::sqrt(std::inner_product(ADCLessPed.begin(), ADCLessPed.begin() + MinBins, ADCLessPed.begin(), 0.) / float(MinBins)));

      // Calculate the power.
      std::vector<double> power(DataSize);
      std::vector<double> RawLessPed;
      RawLessPed.resize(RawADC.size());
      std::transform(RawADC.begin(),RawADC.end(),RawLessPed.begin(),std::bind(std::minus<double>(),std::placeholders::_1,median));
    fFFT->getFFTPower(RawLessPed, power);


    
std::transform(fRawPower.at(boardFound).begin(), fRawPower.at(boardFound).end(), power.begin(), fRawPower.at(boardFound).begin(), std::plus<float>());  ct++; 


      fPed.push_back(median);
  fRawMean.push_back(mean);
std::cout << " before filling frawrms " << std::endl;
fRawRMS.at(boardFound)=rms;  
      fRawRMSTrim.push_back(truncrms);
    
      fChannel.push_back(RawDigit.Channel());

    }
float integ=0;
std::vector<float> rp=fRawPower.at(0);

for (unsigned int j=0; j< rp.size(); j++) integ+=rp.at(j); 

std::cout << " integral " << integ << std::endl;


  ///////////////////////////
  // "Coh" RawDigits.
  ///////////////////////////

  art::Handle< std::vector<raw::RawDigit> > CohDigitHandle;
  e.getByLabel(fRawDigitModuleLabel, fCohInstance, fRawDigitProcess, CohDigitHandle);

  std::cout << " coh instance " << fCohInstance << std::endl;
 
  for(const auto& RawDigit : *CohDigitHandle)
    {
std::vector<geo::WireID> widVec = geom->ChannelToWire(RawDigit.Channel());

 //       size_t                   plane  = widVec[0].Plane;
 //size_t                   wire  = widVec[0].Wire;
std::cout << " before channel mapping " << std::endl;
    const icarusDB::TPCReadoutBoardToChannelMap& readoutBoardToChannelMap = fChannelMap->getReadoutBoardToChannelMap();
int boardFound=-1;
 for(const auto& boardPair : readoutBoardToChannelMap)
    {
 // A little song and dance to make sure this board is in our TPC group
        // What we need is a proper mapping but because we can have channels which are not valid for the geometry
        // service we have a little thing we have to go through here...
        std::vector<geo::WireID> wireIDVec;

        for(const auto& channelPair : boardPair.second.second)
        {
if(channelPair.first==RawDigit.Channel())
boardFound=boardPair.first;
        }
}
std::cout << " board found " << boardFound << std::endl;

//std::vector<short>* fRawPower=RawPowers[boardFound];
 std::vector<short> RawADC;
      // Grab raw waveform, ensuring that the size is set appropriately.
      unsigned int DataSize = RawDigit.Samples();

      RawADC.resize(DataSize);
      raw::Uncompress(RawDigit.ADCs(), RawADC, RawDigit.Compression());

      // We need a sorted waveform (by absolute value) for the truncated RMS and median calculation.
      std::vector<short> SortedADC(RawADC);
      std::sort(SortedADC.begin(),SortedADC.end(),[](const auto& left, const auto& right){return std::fabs(left) < std::fabs(right);});
      float median(SortedADC.at(SortedADC.size()/2));

      // Calculate mean values.
      float mean(float(std::accumulate(SortedADC.begin(),SortedADC.end(),0))/float(SortedADC.size()));





      // Remove pedestal of waveform.
      std::vector<float> ADCLessPed;
      ADCLessPed.resize(SortedADC.size());
      std::transform(SortedADC.begin(),SortedADC.end(),ADCLessPed.begin(),std::bind(std::minus<float>(),std::placeholders::_1,median));

      // Calculate full RMS.
      float rms(std::sqrt(std::inner_product(ADCLessPed.begin(), ADCLessPed.end(), ADCLessPed.begin(), 0.) / float(ADCLessPed.size())));

      // Calculate the truncated RMS.
      unsigned int MinBins((1.0 - 0.01)*ADCLessPed.size());
      //unsigned int BinsToKeep;
      //for(BinsToKeep = 0; BinsToKeep < ADCLessPed.size(); ++BinsToKeep)
      //{
      //  if(std::fabs(ADCLessPed.at(BinsToKeep)) >= 3*rms) break;
      //}
      //float truncrms(std::sqrt(std::inner_product(ADCLessPed.begin(), ADCLessPed.begin() + BinsToKeep, ADCLessPed.begin(), 0.) / float(BinsToKeep)));
      float truncrms(std::sqrt(std::inner_product(ADCLessPed.begin(), ADCLessPed.begin() + MinBins, ADCLessPed.begin(), 0.) / float(MinBins)));

      // Calculate the power.
      std::vector<double> power(DataSize);
      std::vector<double> RawLessPed;
      RawLessPed.resize(RawADC.size());
      std::transform(RawADC.begin(),RawADC.end(),RawLessPed.begin(),std::bind(std::minus<double>(),std::placeholders::_1,median));
    fFFT->getFFTPower(RawLessPed, power);


    
std::transform(fCohPower.at(boardFound).begin(), fCohPower.at(boardFound).end(), power.begin(), fCohPower.at(boardFound).begin(), std::plus<float>());  ct++; 


      fPed.push_back(median);
  fCohMean.push_back(mean);
fCohRMS.at(boardFound)=rms;  
      fCohRMSTrim.push_back(truncrms);
    
    }

  ///////////////////////////
  // "Int" RawDigits.
  ///////////////////////////

  art::Handle< std::vector<raw::RawDigit> > IntDigitHandle;
  e.getByLabel(fRawDigitModuleLabel, fIntInstance, fRawDigitProcess, IntDigitHandle);

  std::cout << " int instance " << fIntInstance << std::endl;
 /* std::vector<int,std::vector<short>* > RawPowers;
int maxBoard=1000;
for(int jm=0;jm<maxBoard;jm++) {
 std::vector<short>* v;
v->clear();
RawPowers.push_back(v);
}

*/
  for(const auto& RawDigit : *IntDigitHandle)
    {
std::vector<geo::WireID> widVec = geom->ChannelToWire(RawDigit.Channel());

 //       size_t                   plane  = widVec[0].Plane;
 //size_t                   wire  = widVec[0].Wire;
std::cout << " before channel mapping " << std::endl;
    const icarusDB::TPCReadoutBoardToChannelMap& readoutBoardToChannelMap = fChannelMap->getReadoutBoardToChannelMap();
int boardFound=-1;
 for(const auto& boardPair : readoutBoardToChannelMap)
    {
 // A little song and dance to make sure this board is in our TPC group
        // What we need is a proper mapping but because we can have channels which are not valid for the geometry
        // service we have a little thing we have to go through here...
        std::vector<geo::WireID> wireIDVec;

        for(const auto& channelPair : boardPair.second.second)
        {
if(channelPair.first==RawDigit.Channel())
boardFound=boardPair.first;
        }
}
std::cout << " board found " << boardFound << std::endl;

//std::vector<short>* fRawPower=RawPowers[boardFound];
 std::vector<short> RawADC;
      // Grab raw waveform, ensuring that the size is set appropriately.
      unsigned int DataSize = RawDigit.Samples();

      RawADC.resize(DataSize);
      raw::Uncompress(RawDigit.ADCs(), RawADC, RawDigit.Compression());

      // We need a sorted waveform (by absolute value) for the truncated RMS and median calculation.
      std::vector<short> SortedADC(RawADC);
      std::sort(SortedADC.begin(),SortedADC.end(),[](const auto& left, const auto& right){return std::fabs(left) < std::fabs(right);});
      float median(SortedADC.at(SortedADC.size()/2));

      // Calculate mean values.
      float mean(float(std::accumulate(SortedADC.begin(),SortedADC.end(),0))/float(SortedADC.size()));





      // Remove pedestal of waveform.
      std::vector<float> ADCLessPed;
      ADCLessPed.resize(SortedADC.size());
      std::transform(SortedADC.begin(),SortedADC.end(),ADCLessPed.begin(),std::bind(std::minus<float>(),std::placeholders::_1,median));

      // Calculate full RMS.
      float rms(std::sqrt(std::inner_product(ADCLessPed.begin(), ADCLessPed.end(), ADCLessPed.begin(), 0.) / float(ADCLessPed.size())));

      // Calculate the truncated RMS.
      unsigned int MinBins((1.0 - 0.01)*ADCLessPed.size());
      //unsigned int BinsToKeep;
      //for(BinsToKeep = 0; BinsToKeep < ADCLessPed.size(); ++BinsToKeep)
      //{
      //  if(std::fabs(ADCLessPed.at(BinsToKeep)) >= 3*rms) break;
      //}
      //float truncrms(std::sqrt(std::inner_product(ADCLessPed.begin(), ADCLessPed.begin() + BinsToKeep, ADCLessPed.begin(), 0.) / float(BinsToKeep)));
      float truncrms(std::sqrt(std::inner_product(ADCLessPed.begin(), ADCLessPed.begin() + MinBins, ADCLessPed.begin(), 0.) / float(MinBins)));

      // Calculate the power.
      std::vector<double> power(DataSize);
      std::vector<double> RawLessPed;
      RawLessPed.resize(RawADC.size());
      std::transform(RawADC.begin(),RawADC.end(),RawLessPed.begin(),std::bind(std::minus<double>(),std::placeholders::_1,median));
    fFFT->getFFTPower(RawLessPed, power);


    
std::transform(fIntPower.at(boardFound).begin(), fIntPower.at(boardFound).end(), power.begin(), fIntPower.at(boardFound).begin(), std::plus<float>());  ct++; 


      fPed.push_back(median);
  fIntMean.push_back(mean);
fIntRMS.at(boardFound)=rms;  
      fIntRMSTrim.push_back(truncrms);
    
 

    }


    
  
   
  
  fNoiseTree->Fill();


  ++NEvents;

  return;

}

void TPCNoiseBoard::TPCNoiseBoard::reconfigure(fhicl::ParameterSet const& p)
{
  fRawDigitModuleLabel = p.get< std::string >("RawDigitModuleLabel", std::string("daqTPC"));
  //std::cerr << "fRawDigitModuleLabel: " << fRawDigitModuleLabel << std::endl;
  fRawDigitProcess = p.get< std::string >("RawDigitProcess", std::string("decode"));
  //std::cerr << "fRawDigitProcess: " << fRawDigitProcess << std::endl;
  fRawInstance = p.get< std::string >("RawInstance", std::string("RAW"));
  //std::cerr << "fRawInstance: " << fRawInstance << std::endl;
  fIntInstance = p.get< std::string >("IntrinsicInstance", std::string("."));
  //std::cerr << "fIntrinsicInstance: " << fIntrinsicInstance << std::endl;
fCohInstance = p.get< std::string >("CoherentInstance", std::string("Cor"));
 // std::cout << "fCoherenInstance: " << fCoherentInstance << std::endl;
fHistoFileName = p.get< std::string >("HistoFileName");

  return;

}

void TPCNoiseBoard::TPCNoiseBoard::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;

  fNoiseTree = tfs->makeAndRegister<TTree>("TPCNoiseBoard_t", "TPC Noise");
 /* fRawPowerTree = tfs->makeAndRegister<TTree>("TPCRawPower_t", "TPC Raw Power");
  fIntrinsicPowerTree = tfs->makeAndRegister<TTree>("TPCIntrinsicPower_t", "TPC Intrinsic Power");
fCoherentPowerTree = tfs->makeAndRegister<TTree>("TPCCoherentPower_t", "TPC Coherent Power");*/
  
  // Initialize branches of the noise TTree.
  fNoiseTree->Branch("Event", &fEvent, "Event/I");
  fNoiseTree->Branch("Run", &fRun, "Run/I");
  fNoiseTree->Branch("SubRun", &fSubRun, "SubRun/I");
  fNoiseTree->Branch("Channel", &fChannel);
  fNoiseTree->Branch("Pedestal", &fPed);
  fNoiseTree->Branch("RawMean", &fRawMean);
  fNoiseTree->Branch("RawRMS", &fRawRMS);

  fNoiseTree->Branch("RawTrimmedRMS", &fRawRMSTrim);
  fNoiseTree->Branch("IntrinsicMean", &fIntMean);
  fNoiseTree->Branch("IntrinsicRMS", &fIntRMS);
  fNoiseTree->Branch("IntrinsicTrimmedRMS", &fIntRMSTrim);

  return;
}

void TPCNoiseBoard::TPCNoiseBoard::endJob()
{
double freqBin=0.6103515625;




  // Average the power vectors over boards
  std::cout << "Rawpower size  " << fRawPower.size() << std::endl;
  std::cout << "fChannel size  " << fChannel.size() << std::endl;
  std::vector<float> TMPVect;
 // fRawPowerTree->Branch("Power", &TMPVect);
int count=0;
  for(long unsigned int jrp=0;jrp<fRawPower.size();jrp++)
    {
for(unsigned int j=0;j<(fRawPower.at(jrp)).size();j++) {(fRawPower.at(jrp)).at(j)/=64;}
std::cout << " in loop before channel " << std::endl;
std::cout << " count " << count << " channel " << fChannel[count] << std::endl;
count++;

      std::transform(fRawPower[jrp].begin(), fRawPower[jrp].end(), fRawPower[jrp].begin(), std::bind(std::divides<float>(), std::placeholders::_1, NEvents));
      TMPVect = fRawPower[jrp];
   //   fRawPowerTree->Fill();
     std::cout << " jrp " << jrp << " frawrms " << fRawRMS.at(jrp) << std::endl;
     std::cout << " frawrmshistos size " << fRawRMSHistos.size() << std::endl;
std::string rawPowerHisto="rawFFT"+std::to_string(jrp);

for(unsigned int jv=0;jv<TMPVect.size();jv++) {
 //std::cout << " filling raw power " << jv*freqBin << std::endl;
//if(jv==100)
  fRawPowerHistos[jrp]->Fill(jv*freqBin,TMPVect[jv]);
    }

fRawRMSHistos[jrp]->Fill(fRawRMS.at(jrp));
}

std::cout << " raw power count " << count << std::endl;
//fRawPowerHisto->Scale(1./count);
std::cout << " after filling raw power histo " << std::endl;
  


 TFile *f = new TFile(fHistoFileName.c_str(),"UPDATE");
    
std::cout << " coherent rms histo entries "<< std::endl;

 for(long unsigned int jrp=0;jrp<fRawPower.size();jrp++)
    {
std::cout << " jrp " << jrp << " raw mean " << fRawRMSHistos[jrp]->GetMean() << std::endl;
if(fRawPowerHistos[jrp]->Integral())
  fRawPowerHistos[jrp]->Write();
if(fRawRMSHistos[jrp]->GetMean())
  fRawRMSHistos[jrp]->Write();
}
  // Average the power vectors over boards
  std::cout << "Cohpower size  " << fCohPower.size() << std::endl;
  std::cout << "fChannel size  " << fChannel.size() << std::endl;

 // fCohPowerTree->Branch("Power", &TMPVect);
 count=0;
  for(long unsigned int jrp=0;jrp<fCohPower.size();jrp++)
    {
for(unsigned int j=0;j<(fCohPower.at(jrp)).size();j++) {(fCohPower.at(jrp)).at(j)/=64;}
std::cout << " in loop before channel " << std::endl;
std::cout << " count " << count << " channel " << fChannel[count] << std::endl;
count++;

      std::transform(fCohPower[jrp].begin(), fCohPower[jrp].end(), fCohPower[jrp].begin(), std::bind(std::divides<float>(), std::placeholders::_1, NEvents));
      TMPVect = fCohPower[jrp];
   //   fRawPowerTree->Fill();
     std::cout << " tmpvect size " << TMPVect.size() << std::endl;
std::string cohPowerHisto="cohFFT"+std::to_string(jrp);


for(unsigned int jv=0;jv<TMPVect.size();jv++) {
 //std::cout << " filling raw power " << jv*freqBin << std::endl;
//if(jv==100)
  fCohPowerHistos[jrp]->Fill(jv*freqBin,TMPVect[jv]);
fCohRMSHistos[jrp]->Fill(fCohRMS.at(jrp));
    }
}

std::cout << " coh power count " << count << std::endl;
//fRawPowerHisto->Scale(1./count);
std::cout << " after filling coh power histo " << std::endl;
  


 
std::cout << " coherent rms histo entries "<< std::endl;

 for(long unsigned int jrp=0;jrp<fCohPower.size();jrp++)
    {
std::cout << " jrp " << jrp << " coh mean " << fCohRMSHistos[jrp]->GetMean() << std::endl;
if(fCohPowerHistos[jrp]->Integral())
  fCohPowerHistos[jrp]->Write();
if(fCohRMSHistos[jrp]->GetMean())
  fCohRMSHistos[jrp]->Write();
}

  // Average the power vectors over boards
  std::cout << "Rawpower size  " << fIntPower.size() << std::endl;
  std::cout << "fChannel size  " << fChannel.size() << std::endl;

  //fIntPowerTree->Branch("Power", &TMPVect);
 count=0;
  for(long unsigned int jrp=0;jrp<fIntPower.size();jrp++)
    {
for(unsigned int j=0;j<(fIntPower.at(jrp)).size();j++) {(fIntPower.at(jrp)).at(j)/=64;}
std::cout << " in loop before channel " << std::endl;
std::cout << " count " << count << " channel " << fChannel[count] << std::endl;
count++;

      std::transform(fIntPower[jrp].begin(), fIntPower[jrp].end(), fIntPower[jrp].begin(), std::bind(std::divides<float>(), std::placeholders::_1, NEvents));
      TMPVect = fIntPower[jrp];
   //   fRawPowerTree->Fill();
     std::cout << " tmpvect size " << TMPVect.size() << std::endl;
std::string intPowerHisto="intFFT"+std::to_string(jrp);


for(unsigned int jv=0;jv<TMPVect.size();jv++) {
 //std::cout << " filling raw power " << jv*freqBin << std::endl;
//if(jv==100)
  fIntPowerHistos[jrp]->Fill(jv*freqBin,TMPVect[jv]);
fIntRMSHistos[jrp]->Fill(fIntRMS.at(jrp));
    }
}

std::cout << " int power count " << count << std::endl;
//fRawPowerHisto->Scale(1./count);
std::cout << " after filling int power histo " << std::endl;
 

 for(long unsigned int jrp=0;jrp<fIntPower.size();jrp++)
    {
std::cout << " jrp " << jrp << " int mean " << fIntRMSHistos[jrp]->GetMean() << std::endl;
if(fIntPowerHistos[jrp]->Integral())
  fIntPowerHistos[jrp]->Write();
if(fIntRMSHistos[jrp]->GetMean())
  fIntRMSHistos[jrp]->Write();
}


  f->Close();
        f->Delete();
std::cout << " after filling c histos " << std::endl;
  std::cout << "Ending job..." << std::endl;
//exit(11);
  return;
}

DEFINE_ART_MODULE(TPCNoiseBoard::TPCNoiseBoard)
