////////////////////////////////////////////////////////////////////////
// Class:       TPCNoise
// Plugin Type: producer (art v3_05_01)
// File:        TPCNoise_module.cc
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

namespace tpcnoise {
  class TPCNoise;
}


class tpcnoise::TPCNoise : public art::EDAnalyzer {
public:
  explicit TPCNoise(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  //TPCNoise(TPCNoise const&) = delete;
  //TPCNoise(TPCNoise&&) = delete;
  //TPCNoise& operator=(TPCNoise const&) = delete;
  //TPCNoise& operator=(TPCNoise&&) = delete;

  void analyze(const art::Event& e);
  void reconfigure(fhicl::ParameterSet const& pset);
  void beginJob();
  void beginRun();
  void endJob();
  void endRun();

private:
  // Various FHiCL parameters.
  std::string fRawDigitModuleLabel;
  std::string fRawDigitProcess;
  std::string fRawInstance;
  std::string fIntrinsicInstance;
std::string fCoherentInstance;
std::string fHistoFileName;


  // FFT calculation.
  using FFTPointer = std::unique_ptr<icarus_signal_processing::ICARUSFFT<double>>;
  FFTPointer fFFT;
  int NumberTimeSamples;

  // FFT variables.
  std::vector< std::vector<float> > fRawPowerC;
  std::vector< std::vector<float> > fIntrinsicPowerC;
std::vector< std::vector<float> > fCoherentPowerC;
 // FFT variables.
  std::vector< std::vector<float> > fRawPowerI1;
  std::vector< std::vector<float> > fIntrinsicPowerI1;
std::vector< std::vector<float> > fCoherentPowerI1;
 // FFT variables.
  std::vector< std::vector<float> > fRawPowerI2;
  std::vector< std::vector<float> > fIntrinsicPowerI2;
std::vector< std::vector<float> > fCoherentPowerI2;

// average FFT histos
TH1D* fRawPowerHistoI1;
 TH1D* fIntrinsicPowerHistoI1;
TH1D* fCoherentPowerHistoI1;
TH1D* fRawRMSHistoI1;
TH1D* fIntrinsicRMSHistoI1;
TH1D* fCoherentRMSHistoI1;
TH1D* fMediaHistoI1;
// average FFT histos
TH1D* fRawPowerHistoI2;
 TH1D* fIntrinsicPowerHistoI2;
TH1D* fCoherentPowerHistoI2;
TH1D* fRawRMSHistoI2;
TH1D* fIntrinsicRMSHistoI2;
TH1D* fCoherentRMSHistoI2;
TH1D* fMediaHistoI2;
// average FFT histos
TH1D* fRawPowerHistoC;
 TH1D* fIntrinsicPowerHistoC;
TH1D* fCoherentPowerHistoC;
TH1D* fRawRMSHistoC;
TH1D* fIntrinsicRMSHistoC;
TH1D* fCoherentRMSHistoC;
TH1D* fMediaHistoC;
  // The variables that will go into the n-tuple.
  int fEvent;
  int fRun;
  int fSubRun;
  std::vector<float> fPed;
  std::vector<float> fRawMeanC;
 std::vector<float> fRawMeanI1;
 std::vector<float> fRawMeanI2;
  std::vector<double> fRawRMSC;
 std::vector<double> fRawRMSI1;
 std::vector<double> fRawRMSI2;
  std::vector<double> fRawRMSTrim;
  std::vector<float> fIntrinsicMean;
  std::vector<double> fIntrinsicRMSC;
 std::vector<double> fIntrinsicRMSI1;
 std::vector<double> fIntrinsicRMSI2;
  std::vector<double> fIntrinsicRMSTrim;
  std::vector<float> fCoherentMean;
  std::vector<double> fCoherentRMSC;
 std::vector<double> fCoherentRMSI1;
 std::vector<double> fCoherentRMSI2;
  std::vector<double> fCoherentRMSTrim;
  std::vector<unsigned short int> fChannel;

  // The variables that will go into the power n-tuple.
  //std::vector< std::vector<float> > fPower;
  unsigned int NEvents;
 
  // The output trees.
  TTree* fNoiseTree;
  TTree* fRawPowerTree;
  TTree* fIntrinsicPowerTree;
 TTree* fCoherentPowerTree;
};


tpcnoise::TPCNoise::TPCNoise(fhicl::ParameterSet const& p)
  : EDAnalyzer(p),
    NEvents(0)
{
  // Retrieve proper number of time samples.
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob();
  NumberTimeSamples = detProp.NumberTimeSamples();

  std::cout << "Number of time samples: " << NumberTimeSamples << std::endl;

  fRawPowerC.resize(1);
  for(auto &it : fRawPowerC)
    {
      it.resize(NumberTimeSamples);
    }
std::cout << " after resizing " << std::endl;
  fIntrinsicPowerC.resize(1);
  for(auto &it : fIntrinsicPowerC)
    {
      it.resize(NumberTimeSamples);
    }
std::cout << " after resizing " << std::endl;
 fCoherentPowerC.resize(1);
  for(auto &it : fCoherentPowerC)
    {
      it.resize(NumberTimeSamples);
    }

std::cout << " after resizing " << std::endl;
 fRawPowerI1.resize(1);
  for(auto &it : fRawPowerI1)
    {
      it.resize(NumberTimeSamples);
    }
std::cout << " after resizing " << std::endl;
  fIntrinsicPowerI1.resize(1);
  for(auto &it : fIntrinsicPowerI1)
    {
      it.resize(NumberTimeSamples);
    }
std::cout << " after resizing " << std::endl;
 fCoherentPowerI1.resize(1);
  for(auto &it : fCoherentPowerI1)
    {
      it.resize(NumberTimeSamples);
    } 
fCoherentPowerI2.resize(1);
  for(auto &it : fCoherentPowerI2)
    {
      it.resize(NumberTimeSamples);
    }
std::cout << " after resizing " << std::endl;
 fRawPowerI2.resize(1);
  for(auto &it : fRawPowerI2)
    {
      it.resize(NumberTimeSamples);
    }
std::cout << " after resizing " << std::endl;
  fIntrinsicPowerI2.resize(1);
  for(auto &it : fIntrinsicPowerI2)
    {
      it.resize(NumberTimeSamples);
    }
std::cout << " after intrinsic i2 resizing " << std::endl;

  fFFT = std::make_unique<icarus_signal_processing::ICARUSFFT<double>>(NumberTimeSamples);
  this->reconfigure(p);
double freqBin=0.6103515625;
fRawPowerHistoC=new TH1D("rawpowerC","rawpowerC",2048,0.,2048*freqBin);
fIntrinsicPowerHistoC=new TH1D("intpowerC","intpowerC",2048,0.,2048*freqBin);
fCoherentPowerHistoC=new TH1D("cohpowerC","cohpowerC",2048,0.,2048*freqBin);
fRawRMSHistoC=new TH1D("rawRMSC","rawRMSC",100,0.,30.);
fMediaHistoC=new TH1D("mediaC","mediaC",100,-10.,10.);
fIntrinsicRMSHistoC=new TH1D("intRMSC","intRMSC",100,0.,30.);
fCoherentRMSHistoC=new TH1D("cohRMSC","cohRMSC",100,0.,30.);
fRawPowerHistoI1=new TH1D("rawpowerI1","rawpowerI1",2048,0.,2048*freqBin);
fIntrinsicPowerHistoI1=new TH1D("intpowerI1","intpowerI1",2048,0.,2048*freqBin);
fCoherentPowerHistoI1=new TH1D("cohpowerI1","cohpowerI1",2048,0.,2048*freqBin);
fRawRMSHistoI1=new TH1D("rawRMSI1","rawRMSI1",100,0.,30.);
fMediaHistoI1=new TH1D("mediaI1","mediaI1",100,-10.,10.);
fIntrinsicRMSHistoI1=new TH1D("intRMSI1","intRMSI1",100,0.,30.);
fCoherentRMSHistoI1=new TH1D("cohRMSI1","cohRMSI1",100,0.,30.);
fRawPowerHistoI2=new TH1D("rawpowerI2","rawpowerI2",2048,0.,2048*freqBin);
fIntrinsicPowerHistoI2=new TH1D("intpowerI2","intpowerI2",2048,0.,2048*freqBin);
fCoherentPowerHistoI2=new TH1D("cohpowerI2","cohpowerI2",2048,0.,2048*freqBin);
fRawRMSHistoI2=new TH1D("rawRMSI2","rawRMSI2",100,0.,30.);
fMediaHistoI2=new TH1D("mediaI2","mediaI2",100,-10.,10.);
fIntrinsicRMSHistoI2=new TH1D("intRMSI2","intRMSI2",100,0.,30.);
fCoherentRMSHistoI2=new TH1D("cohRMSI2","cohRMSI2",100,0.,30.);

std::cout << " end constructor " << std::endl;
}

void tpcnoise::TPCNoise::analyze(const art::Event& e)
{
std::cout << " begin analyze " << std::endl;
   art::ServiceHandle<geo::Geometry> geom;

  // Clear vectors before filling for this event.
  fChannel.clear();
std::cout << " after clearing " << std::endl;
  fPed.clear();
  fRawMeanC.clear();
fRawMeanI1.clear();
fRawMeanI2.clear();
  fRawRMSC.clear();
  fCoherentRMSC.clear();
  fRawRMSTrim.clear();
fRawRMSI1.clear();
  fCoherentRMSI1.clear();
fRawRMSI2.clear();
  fCoherentRMSI2.clear();
  fIntrinsicMean.clear();
  fIntrinsicRMSC.clear();
fIntrinsicRMSI1.clear();
fIntrinsicRMSI2.clear();
  fIntrinsicRMSTrim.clear();
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
  int ctI1=0; int ctI2=0; int ctC=0;
  for(const auto& RawDigit : *RawDigitHandle)
    {
      // Grab raw waveform, ensuring that the size is set appropriately.
      unsigned int DataSize = RawDigit.Samples();
      std::vector<short> RawADC;
      RawADC.resize(DataSize);
      raw::Uncompress(RawDigit.ADCs(), RawADC, RawDigit.Compression());

      // We need a sorted waveform (by absolute value) for the truncated RMS and median calculation.
      std::vector<short> SortedADC(RawADC);
      std::sort(SortedADC.begin(),SortedADC.end(),[](const auto& left, const auto& right){return std::fabs(left) < std::fabs(right);});
      float median(SortedADC.at(SortedADC.size()/2));

      // Calculate mean values.
      float mean(float(std::accumulate(SortedADC.begin(),SortedADC.end(),0))/float(SortedADC.size()));
std::vector<geo::WireID> widVec = geom->ChannelToWire(RawDigit.Channel());
        size_t                   plane  = widVec[0].Plane;
 size_t                   wire  = widVec[0].Wire;


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

      if(plane==0) std::cout << " Induction 1 channel " <<std::endl;
if(plane==0&&wire!=32&&wire<3000)  {std::cout << " i1 power wire " << wire << " channel " << RawDigit.Channel() << " ind1 raw RMS " << rms << " ctI1 " << ctI1 <<std::endl; 
std::transform(fRawPowerI1.at(0).begin(), fRawPowerI1.at(0).end(), power.begin(), fRawPowerI1.at(0).begin(), std::plus<float>());  ctI1++; }
if(plane==1)  { std::transform(fRawPowerI2.at(0).begin(), fRawPowerI2.at(0).end(), power.begin(), fRawPowerI2.at(0).begin(), std::plus<float>()); ctI2++; }
if(plane==2)  { std::transform(fRawPowerC.at(0).begin(), fRawPowerC.at(0).end(), power.begin(), fRawPowerC.at(0).begin(), std::plus<float>()); ctC++; }  


      fPed.push_back(median);
   if(plane==2)   fRawMeanC.push_back(mean);
if(plane==1)   fRawMeanI2.push_back(mean);
if(plane==0)   fRawMeanI1.push_back(mean);
  if(plane==2)   { fRawRMSC.push_back(rms); std::cout << " coll raw RMS " << rms << " rawrmsc size " << fRawRMSC.size() <<  std::endl; }
if(plane==0&&wire!=32&&wire<3000) { std::cout << " filling I1 rms " << std::endl; fRawRMSI1.push_back(rms);   }
 if(plane==1) { fRawRMSI2.push_back(rms);std::cout << " ind2 raw RMS " << rms << std::endl; }
      fRawRMSTrim.push_back(truncrms);
    
      fChannel.push_back(RawDigit.Channel());

    }
float intI1=0, intI2=0, intC=0;
std::vector<float> rpi1=fRawPowerI1.at(0);
std::vector<float> rpi2=fRawPowerI2.at(0);
std::vector<float> rpc=fRawPowerC.at(0);
for (unsigned int j=0; j< rpi1.size(); j++) intI1+=rpi1.at(j); 
for (unsigned int j=0; j< rpi2.size(); j++) intI2+=rpi2.at(j); 
for (unsigned int j=0; j< rpc.size(); j++) intC+=rpc.at(j); 

std::cout << " integral i1 " << intI1 << " integral i2 " << intI2 << " integral c " << intC << std::endl;
  ///////////////////////////
  // "Intrinsic" RawDigits.
  ///////////////////////////
   art::Handle< std::vector<raw::RawDigit> > IntrinsicHandle;
  e.getByLabel(fRawDigitModuleLabel, fIntrinsicInstance, fRawDigitProcess, IntrinsicHandle);

  //std::cerr << RawDigitHandle << std::endl;
 //ctI1=0; ctI2=0; ctC=0;
std::cout << " intrinsic instance " << fIntrinsicInstance << std::endl;
  for(const auto& RawDigit : *IntrinsicHandle)
    {
      // Grab raw waveform, ensuring that the size is set appropriately.
      unsigned int DataSize = RawDigit.Samples();
      std::vector<short> RawADC;
      RawADC.resize(DataSize);
      raw::Uncompress(RawDigit.ADCs(), RawADC, RawDigit.Compression());

      // We need a sorted waveform (by absolute value) for the truncated RMS and median calculation.
      std::vector<short> SortedADC(RawADC);
      std::sort(SortedADC.begin(),SortedADC.end(),[](const auto& left, const auto& right){return std::fabs(left) < std::fabs(right);});
      float median(SortedADC.at(SortedADC.size()/2));

      // Calculate mean values.
      float mean(float(std::accumulate(SortedADC.begin(),SortedADC.end(),0))/float(SortedADC.size()));

      // Remove pedestal of waveform.
      std::vector<short> ADCLessPed;
      ADCLessPed.resize(SortedADC.size());
      std::transform(SortedADC.begin(),SortedADC.end(),ADCLessPed.begin(),std::bind(std::minus<float>(),std::placeholders::_1,median));

      // Calculate full RMS.
      float rms(std::sqrt(std::inner_product(ADCLessPed.begin(), ADCLessPed.end(), ADCLessPed.begin(), 0.) / float(ADCLessPed.size())));

      // Calculate the truncated RMS.
      unsigned int MinBins((1.0 - 0.2)*ADCLessPed.size());
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
std::vector<geo::WireID> widVec = geom->ChannelToWire(RawDigit.Channel());
        size_t                   plane  = widVec[0].Plane;

if(plane==0)  { std::transform(fIntrinsicPowerI1.at(0).begin(), fIntrinsicPowerI1.at(0).end(), power.begin(), fIntrinsicPowerI1.at(0).begin(), std::plus<float>());  }
if(plane==1)  { std::transform(fIntrinsicPowerI2.at(0).begin(), fIntrinsicPowerI2.at(0).end(), power.begin(), fIntrinsicPowerI2.at(0).begin(), std::plus<float>()); }
if(plane==2)  { std::transform(fIntrinsicPowerC.at(0).begin(), fIntrinsicPowerC.at(0).end(), power.begin(), fIntrinsicPowerC.at(0).begin(), std::plus<float>());  }   
     // std::transform(fIntrinsicPower.at(RawDigit.Channel()).begin(), fIntrinsicPower.at(RawDigit.Channel()).end(), power.begin(), fIntrinsicPower.at(RawDigit.Channel()).begin(), std::plus<float>());

      fIntrinsicMean.push_back(mean);
  if(plane==2)  {  fIntrinsicRMSC.push_back(rms);std::cout << " coll intrinsic RMS " << rms << std::endl; }
 if(plane==0)fIntrinsicRMSI1.push_back(rms);
 if(plane==1) {fIntrinsicRMSI2.push_back(rms);std::cout << " ind2 intrinsic RMS " << rms << std::endl; }

      fIntrinsicRMSTrim.push_back(truncrms);
    }

  ///////////////////////////
  // "Coherent" RawDigits.
  ///////////////////////////
  
  //RawDigitHandle.clear();
  //e.getByLabel(art::InputTag(fRawDigitModuleLabel), RawDigitHandle);

   art::Handle< std::vector<raw::RawDigit> > CoherentHandle;
  e.getByLabel(fRawDigitModuleLabel, fCoherentInstance, fRawDigitProcess, CoherentHandle);
std::cout << " coherent instance " << fCoherentInstance << std::endl;
  //std::cerr << RawDigitHandle << std::endl;
//ctI1=0; ctI2=0; ctC=0;
  for(const auto& RawDigit : *CoherentHandle)
    {

      // Grab raw waveform, ensuring that the size is set appropriately.
      unsigned int DataSize = RawDigit.Samples();
      std::vector<short> RawADC;
      RawADC.resize(DataSize);
      raw::Uncompress(RawDigit.ADCs(), RawADC, RawDigit.Compression());

      // We need a sorted waveform (by absolute value) for the truncated RMS and median calculation.
      std::vector<short> SortedADC(RawADC);
      std::sort(SortedADC.begin(),SortedADC.end(),[](const auto& left, const auto& right){return std::fabs(left) < std::fabs(right);});
      float median(SortedADC.at(SortedADC.size()/2));
//for(unsigned int jadc=0;jadc<SortedADC.size();jadc++)
 //if(SortedADC.at(jadc))
  //std::cout << " sorted adc " << SortedADC.at(jadc) << std::endl;

      // Calculate mean values.
      float mean(float(std::accumulate(SortedADC.begin(),SortedADC.end(),0))/float(SortedADC.size()));

      // Remove pedestal of waveform.
      std::vector<short> ADCLessPed;
      ADCLessPed.resize(SortedADC.size());
      std::transform(SortedADC.begin(),SortedADC.end(),ADCLessPed.begin(),std::bind(std::minus<float>(),std::placeholders::_1,median));

      // Calculate full RMS.
      float rms(std::sqrt(std::inner_product(ADCLessPed.begin(), ADCLessPed.end(), ADCLessPed.begin(), 0.) / float(ADCLessPed.size())));

      // Calculate the truncated RMS.
      unsigned int MinBins((1.0 - 0.2)*ADCLessPed.size());
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
std::vector<geo::WireID> widVec = geom->ChannelToWire(RawDigit.Channel());
        size_t                   plane  = widVec[0].Plane;
     
if(plane==0)  { std::transform(fCoherentPowerI1.at(0).begin(), fCoherentPowerI1.at(0).end(), power.begin(), fCoherentPowerI1.at(0).begin(), std::plus<float>());  }
if(plane==1)  { std::transform(fCoherentPowerI2.at(0).begin(), fCoherentPowerI2.at(0).end(), power.begin(), fCoherentPowerI2.at(0).begin(), std::plus<float>()); }
if(plane==2)  { std::transform(fCoherentPowerC.at(0).begin(), fCoherentPowerC.at(0).end(), power.begin(), fCoherentPowerC.at(0).begin(), std::plus<float>());  }

//if(rms) std::cout << " coherent mean " << mean << " rms " << rms << std::endl;
      fCoherentMean.push_back(mean);
if(rms) {
  if(plane==2)  {  fCoherentRMSC.push_back(rms);std::cout << " coll coherent RMS " << rms << std::endl; }
 if(plane==0)fCoherentRMSI1.push_back(rms);
 if(plane==1) { fCoherentRMSI2.push_back(rms);    std::cout << " ind2 coherent RMS " << rms << std::endl; }
      fCoherentRMSTrim.push_back(truncrms);
}
    }
//std::cout << " cohrms size " << fCoherentRMSC.size() << std::endl;
//for(int j=0;j<10;j++) std::cout << " j " << j << " cohrms C " << fCoherentRMSC.at(j) << " cohrms I2 " << fCoherentRMSI2.at(j) << std::endl;
  fNoiseTree->Fill();

std::cout << " ctI1 " << ctI1 << " ctI2 " << ctI2 << " ctC " << ctC << std::endl;
for(unsigned int j=0;j<(fRawPowerI1.at(0)).size();j++) {(fRawPowerI1.at(0)).at(j)/=ctI1;}
for(unsigned int j=0;j<(fCoherentPowerI1.at(0)).size();j++) {(fCoherentPowerI1.at(0)).at(j)/=ctI1;}
for(unsigned int j=0;j<(fIntrinsicPowerI1.at(0)).size();j++) {(fIntrinsicPowerI1.at(0)).at(j)/=ctI1;}
for(unsigned int j=0;j<(fRawPowerI2.at(0)).size();j++) {(fRawPowerI2.at(0)).at(j)/=ctI2;}
for(unsigned int j=0;j<(fCoherentPowerI2.at(0)).size();j++) {(fCoherentPowerI2.at(0)).at(j)/=ctI2;}
for(unsigned int j=0;j<(fIntrinsicPowerI2.at(0)).size();j++) {(fIntrinsicPowerI2.at(0)).at(j)/=ctI2;}
for(unsigned int j=0;j<(fRawPowerC.at(0)).size();j++) {(fRawPowerC.at(0)).at(j)/=ctC;}
for(unsigned int j=0;j<(fCoherentPowerC.at(0)).size();j++) {(fCoherentPowerC.at(0)).at(j)/=ctC;}
for(unsigned int j=0;j<(fIntrinsicPowerC.at(0)).size();j++) {(fIntrinsicPowerC.at(0)).at(j)/=ctC;}

  ++NEvents;

  return;

}

void tpcnoise::TPCNoise::reconfigure(fhicl::ParameterSet const& p)
{
  fRawDigitModuleLabel = p.get< std::string >("RawDigitModuleLabel", std::string("daqTPC"));
  //std::cerr << "fRawDigitModuleLabel: " << fRawDigitModuleLabel << std::endl;
  fRawDigitProcess = p.get< std::string >("RawDigitProcess", std::string("decode"));
  //std::cerr << "fRawDigitProcess: " << fRawDigitProcess << std::endl;
  fRawInstance = p.get< std::string >("RawInstance", std::string("RAW"));
  //std::cerr << "fRawInstance: " << fRawInstance << std::endl;
  fIntrinsicInstance = p.get< std::string >("IntrinsicInstance", std::string("."));
  //std::cerr << "fIntrinsicInstance: " << fIntrinsicInstance << std::endl;
fCoherentInstance = p.get< std::string >("CoherentInstance", std::string("Cor"));
 // std::cout << "fCoherentInstance: " << fCoherentInstance << std::endl;
fHistoFileName = p.get< std::string >("HistoFileName");

  return;

}

void tpcnoise::TPCNoise::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;

  fNoiseTree = tfs->makeAndRegister<TTree>("TPCNoise_t", "TPC Noise");
  fRawPowerTree = tfs->makeAndRegister<TTree>("TPCRawPower_t", "TPC Raw Power");
  fIntrinsicPowerTree = tfs->makeAndRegister<TTree>("TPCIntrinsicPower_t", "TPC Intrinsic Power");
fCoherentPowerTree = tfs->makeAndRegister<TTree>("TPCCoherentPower_t", "TPC Coherent Power");
  
  // Initialize branches of the noise TTree.
  fNoiseTree->Branch("Event", &fEvent, "Event/I");
  fNoiseTree->Branch("Run", &fRun, "Run/I");
  fNoiseTree->Branch("SubRun", &fSubRun, "SubRun/I");
  fNoiseTree->Branch("Channel", &fChannel);
  fNoiseTree->Branch("Pedestal", &fPed);
  fNoiseTree->Branch("RawMeanC", &fRawMeanC);
  fNoiseTree->Branch("RawRMSC", &fRawRMSC);

  fNoiseTree->Branch("RawTrimmedRMS", &fRawRMSTrim);
  fNoiseTree->Branch("IntrinsicMean", &fIntrinsicMean);
  fNoiseTree->Branch("IntrinsicRMSC", &fIntrinsicRMSC);
  fNoiseTree->Branch("IntrinsicTrimmedRMS", &fIntrinsicRMSTrim);

  return;
}

void tpcnoise::TPCNoise::endJob()
{
double freqBin=0.6103515625;
  // Average the power vectors.
  std::cout << "Averaging power vectors..." << std::endl;
  std::vector<float> TMPVect;
  fRawPowerTree->Branch("Power", &TMPVect);
int count=0;
  for(auto &it : fRawPowerC)
    {
count++;
      std::transform(it.begin(), it.end(), it.begin(), std::bind(std::divides<float>(), std::placeholders::_1, NEvents));
      TMPVect = it;
   //   fRawPowerTree->Fill();
     std::cout << " tmpvect size " << TMPVect.size() << std::endl;
for(unsigned int jv=0;jv<TMPVect.size();jv++) {
 std::cout << " filling raw power " << jv*freqBin << std::endl;
//if(jv==100)
  fRawPowerHistoC->Fill(jv*freqBin,TMPVect[jv]);

    }
}
std::cout << " raw power entries " << fRawPowerHistoC->GetEntries() << std::endl;
std::cout << " raw power count " << count << std::endl;
fRawPowerHistoC->Scale(1./count);
std::cout << " after filling raw power histo " << std::endl;
  fIntrinsicPowerTree->Branch("Power", &TMPVect);
count=0;
  for(auto &it : fIntrinsicPowerC)
    {
count++;
      std::transform(it.begin(), it.end(), it.begin(), std::bind(std::divides<float>(), std::placeholders::_1, NEvents));
      TMPVect = it;
  std::cout << " tmpvect size " << TMPVect.size() << std::endl;
   //   fIntrinsicPowerTree->Fill();
  std::cout << " after fill tmpvect size " << TMPVect.size() << std::endl;
for(unsigned int jv=0;jv<TMPVect.size();jv++) {
//std::cout << " filling jv " << jv << std::endl; 

 fIntrinsicPowerHistoC->Fill(jv*freqBin,TMPVect[jv]);
}
fIntrinsicPowerHistoC->Scale(1./count);
std::cout << " after filling intrinsic power histo " << std::endl;
    }

  fCoherentPowerTree->Branch("Power", &TMPVect);
count=0;
  for(auto &it : fCoherentPowerC)
    {
count++;
      std::transform(it.begin(), it.end(), it.begin(), std::bind(std::divides<float>(), std::placeholders::_1, NEvents));
      TMPVect = it;

    //  fCoherentPowerTree->Fill();

for(unsigned int jv=0;jv<TMPVect.size();jv++)
  fCoherentPowerHistoC->Fill(jv*freqBin,TMPVect[jv]);
    }
fCoherentPowerHistoC->Scale(1./count);
std::cout << " frawrmsc size " << fRawRMSC.size() << std::endl;
 for(unsigned int jj=0;jj<fRawRMSC.size();jj++)
  fRawRMSHistoC->Fill(fRawRMSC.at(jj));
std::cout << " after filling raw rms histo " << std::endl;
for(unsigned int j=0;j<fRawMeanC.size(); j++) {
std::cout << " filling media " << std::endl;
  fMediaHistoC->Fill(fRawMeanC.at(j));
}
std::cout << " after filling media histo " << std::endl;

for(unsigned int jj=0;jj<fIntrinsicRMSC.size();jj++)
  fIntrinsicRMSHistoC->Fill(fIntrinsicRMSC.at(jj));

std::cout << " fillhisto cohrms size " << fCoherentRMSC.size() << std::endl;

for(unsigned int j=0;j<fCoherentRMSC.size();j++) {
std::cout << " filling coherent RMS " << fCoherentRMSC.at(j) << std::endl;
  fCoherentRMSHistoC->Fill(fCoherentRMSC.at(j));
}

  // Average the power vectors.
  std::cout << "Averaging power vectors..." << std::endl;
  //std::vector<float> TMPVect;
  fRawPowerTree->Branch("Power", &TMPVect);
count=0;
  for(auto &it : fRawPowerI1)
    {
     std::cout << " i1 counter " << count++ << std::endl;
      std::transform(it.begin(), it.end(), it.begin(), std::bind(std::divides<float>(), std::placeholders::_1, NEvents));
      TMPVect = it;
     // fRawPowerTree->Fill();
    // std::cout << " tmpvect size " << TMPVect.size() << std::endl;
for(unsigned int jv=0;jv<TMPVect.size();jv++)
//if(jv==100)
  fRawPowerHistoI1->Fill(jv*freqBin,TMPVect[jv]);
    }
fRawPowerHistoI1->Scale(1./count);
std::cout << " after filling raw power histo " << std::endl;
  fIntrinsicPowerTree->Branch("Power", &TMPVect);
count=0;
  for(auto &it : fIntrinsicPowerI1)
    {
 std::cout << " i1 counter " << count++ << std::endl;
      std::transform(it.begin(), it.end(), it.begin(), std::bind(std::divides<float>(), std::placeholders::_1, NEvents));
      TMPVect = it;
  //astd::cout << " tmpvect size " << TMPVect.size() << std::endl;
    //  fIntrinsicPowerTree->Fill();
 // std::cout << " after fill tmpvect size " << TMPVect.size() << std::endl;
for(unsigned int jv=0;jv<TMPVect.size();jv++) {
//std::cout << " filling jv " << jv << std::endl; 
 fIntrinsicPowerHistoI1->Fill(jv*freqBin,TMPVect[jv]);
}
fIntrinsicPowerHistoI1->Scale(1./count);
std::cout << " after filling intrinsic power histo " << std::endl;
    }

  fCoherentPowerTree->Branch("Power", &TMPVect);
count=0;
  for(auto &it : fCoherentPowerI1)
    {
 std::cout << " i1 counter " << count++ << std::endl;
      std::transform(it.begin(), it.end(), it.begin(), std::bind(std::divides<float>(), std::placeholders::_1, NEvents));
      TMPVect = it;

     // fCoherentPowerTree->Fill();

for(unsigned int jv=0;jv<TMPVect.size();jv++)
  fCoherentPowerHistoI1->Fill(jv*freqBin,TMPVect[jv]);
    }
fCoherentPowerHistoI1->Scale(1./count);
std::cout << " rawrmsi1 size " << fRawRMSI1.size() << std::endl;
 for(unsigned int jj=0;jj<fRawRMSI1.size();jj++) {
std::cout << " filling rawrmsi1 value " << fRawRMSI1.at(jj) << std::endl;
  fRawRMSHistoI1->Fill(fRawRMSI1.at(jj));
}
std::cout << " after filling raw rms histo entries " << fRawRMSHistoI1->GetEntries() <<  std::endl;
for(unsigned int j=0;j<fRawMeanI1.size(); j++) {
//std::cout << " filling media " << std::endl;
  fMediaHistoI1->Fill(fRawMeanI1.at(j));
}
std::cout << " after filling media histo " << std::endl;

for(unsigned int jj=0;jj<fIntrinsicRMSI1.size();jj++)
  fIntrinsicRMSHistoI1->Fill(fIntrinsicRMSI1.at(jj));

std::cout << " fillhisto cohrms size " << fCoherentRMSI1.size() << std::endl;

for(unsigned int j=0;j<fCoherentRMSI1.size();j++) {
std::cout << " filling coherent RMS " << fCoherentRMSI1.at(j) << std::endl;
  fCoherentRMSHistoI1->Fill(fCoherentRMSI1.at(j));
}

  // Average the power vectors.
  std::cout << "Averaging power vectors..." << std::endl;
  //std::vector<float> TMPVect;
  fRawPowerTree->Branch("Power", &TMPVect);
count=0;
  for(auto &it : fRawPowerI2)
    {
count++;
      std::transform(it.begin(), it.end(), it.begin(), std::bind(std::divides<float>(), std::placeholders::_1, NEvents));
      TMPVect = it;
     // fRawPowerTree->Fill();
     std::cout << " tmpvect size " << TMPVect.size() << std::endl;
for(unsigned int jv=0;jv<TMPVect.size();jv++)
//if(jv==100)
  fRawPowerHistoI2->Fill(jv*freqBin,TMPVect[jv]);
    }
std::cout << " after filling raw power histo " << std::endl;
fRawPowerHistoI2->Scale(1./count);
  fIntrinsicPowerTree->Branch("Power", &TMPVect);
count=0;
  for(auto &it : fIntrinsicPowerI2)
    {
count++;
      std::transform(it.begin(), it.end(), it.begin(), std::bind(std::divides<float>(), std::placeholders::_1, NEvents));
      TMPVect = it;
  std::cout << " tmpvect size " << TMPVect.size() << std::endl;
   //   fIntrinsicPowerTree->Fill();
  std::cout << " after fill tmpvect size " << TMPVect.size() << std::endl;
for(unsigned int jv=0;jv<TMPVect.size();jv++) {
//std::cout << " filling jv " << jv << std::endl; 
 fIntrinsicPowerHistoI2->Fill(jv*freqBin,TMPVect[jv]);
}
fIntrinsicPowerHistoI2->Scale(1./count);
std::cout << " after filling intrinsic power histo " << std::endl;
    }

  fCoherentPowerTree->Branch("Power", &TMPVect);
count=0;
  for(auto &it : fCoherentPowerI2)
    {
count++;
      std::transform(it.begin(), it.end(), it.begin(), std::bind(std::divides<float>(), std::placeholders::_1, NEvents));
      TMPVect = it;

    //  fCoherentPowerTree->Fill();

for(unsigned int jv=0;jv<TMPVect.size();jv++)
  fCoherentPowerHistoI2->Fill(jv*freqBin,TMPVect[jv]);
    }
fCoherentPowerHistoI2->Scale(1./count);
std::cout << " after filling coherent power histo " << std::endl;
 for(unsigned int jj=0;jj<fRawRMSI2.size();jj++)
  fRawRMSHistoI2->Fill(fRawRMSI2.at(jj));
std::cout << " after filling raw rms histo " << std::endl;
for(unsigned int j=0;j<fRawMeanI2.size(); j++) {
//std::cout << " filling media " << std::endl;
  fMediaHistoI2->Fill(fRawMeanI2.at(j));
}
std::cout << " after filling media histo " << std::endl;

for(unsigned int jj=0;jj<fIntrinsicRMSI2.size();jj++)
  fIntrinsicRMSHistoI2->Fill(fIntrinsicRMSI2.at(jj));

std::cout << " fillhisto cohrms size " << fCoherentRMSI2.size() << std::endl;

for(unsigned int j=0;j<fCoherentRMSI2.size();j++) {
std::cout << " filling coherent RMS " << fCoherentRMSI2.at(j) << std::endl;
  fCoherentRMSHistoI2->Fill(fCoherentRMSI2.at(j));
}

 TFile *f = new TFile(fHistoFileName.c_str(),"RECREATE");
    
std::cout << " coherent rms histo entries "<< std::endl;

  fRawPowerHistoI2->Write();
fIntrinsicPowerHistoI2->Write();
fCoherentPowerHistoI2->Write();
 fRawRMSHistoI2->Write();
fMediaHistoI2->Write();
fIntrinsicRMSHistoI2->Write();
fCoherentRMSHistoI2->Write();

std::cout << " after filling i2 histos " << std::endl;

  fRawPowerHistoI1->Write();
fIntrinsicPowerHistoI1->Write();
fCoherentPowerHistoI1->Write();
 fRawRMSHistoI1->Write();
fMediaHistoI1->Write();
fIntrinsicRMSHistoI1->Write();
fCoherentRMSHistoI1->Write();

std::cout << " after filling i1 histos " << std::endl;

  fRawPowerHistoC->Write();
fIntrinsicPowerHistoC->Write();
fCoherentPowerHistoC->Write();
 fRawRMSHistoC->Write();
fMediaHistoC->Write();
fIntrinsicRMSHistoC->Write();
fCoherentRMSHistoC->Write();
  f->Close();
        f->Delete();
std::cout << " after filling c histos " << std::endl;
  std::cout << "Ending job..." << std::endl;
//exit(11);
  return;
}

DEFINE_ART_MODULE(tpcnoise::TPCNoise)
