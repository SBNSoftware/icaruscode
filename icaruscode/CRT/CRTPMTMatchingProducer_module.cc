// icaruscode includes
#include "sbnobj/ICARUS/CRT/CRTData.hh"
#include "sbnobj/Common/CRT/CRTHit.hh"
#include "icaruscode/CRT/CRTUtils/CRTHitRecoAlg.h"
#include "sbnobj/Common/Trigger/ExtraTriggerInfo.h"
#include "icaruscode/IcarusObj/CRTPMTMatching.h"

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Event.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileService.h" 
#include "art_root_io/TFileDirectory.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Utilities/Exception.h"

// C++ includes
#include <memory>
#include <iostream>
#include <map>
#include <iterator>
#include <algorithm>
#include <vector>

// LArSoft
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/AuxDetGeometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RawData/ExternalTrigger.h"
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

// ROOT
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TVector3.h"
#include "TGeoManager.h"

using std::vector;

namespace icarus::crt {

bool flashInTime(double flashTime, int gateType, double gateDiff,
                 double gateWidth) {
  // As a reminder, I will leave here the commented part of the vetoOffset, in
  // case something changes in the future
  /*double vetoOffset = 3.5;*/
  double activeGate = gateWidth /*- vetoOffset*/;

  double relFlashTime = flashTime + gateDiff / 1000. /*- vetoOffset*/;
  mf::LogDebug("CRTPMTMatchingProducer FlashInTime")
      << "Gate Diff " << gateDiff / 1000 << " Ftime+gateDiff "
      << flashTime + gateDiff / 1000. << " " << activeGate;

  return ((relFlashTime > 0) && (relFlashTime < activeGate));
};

enum hasCRTHit {
  /**
   * Type of the matching.
   *
   * Secret decoder ring:
   *  * `0`: No CRT match
   *  * `1`: matched with Top CRT hit before optical flash
   *  * `2`: matched with Side CRT hit before optical flash
   *  * `3`: entering from Top and exiting from side
   *  * `4`: matched with a Top CRT after the optical flash
   *  * `5`: matched with a Side CRT after the optical flash
   *  * `6`: matched with multiple Top CRT hits before the optical flash
   *  * `7`: matched with multiple Top CRT hits before the optical flash and
   * more then 1 side CRT hits after the optical flash
   *  * `9`: all the other cases
   *
   */
  noMatch = 0,
  enTop = 1,
  enSide = 2,
  enTop_exSide = 3,
  exTop = 4,
  exSide = 5,
  enTop_mult = 6,
  enTop_exSide_mult = 7,
  others = 9
};

struct CRTPMT {
  double 	tof;    	///< Time difference between CRT Hit and optical flash [ns]
  double 	distance;    	///< Distance between CRT Hit and optical flash centroid [cm]
  art::Ptr<sbn::crt::CRTHit> CRTHit;
};

struct CRTMatches {
  // vector of pairs where first is the matched Time of Flight and second is the
  // matched CRTHit
  std::vector<CRTPMT> 	entering;
  std::vector<CRTPMT> 	exiting;
  hasCRTHit matchType;
};

struct MatchedCRT {
  geo::Point_t CRTHitPos;
  double CRTPMTTimeDiff_ns;
  double CRTTime_us;
  int CRTSys;
  int CRTRegion;
};

struct FlashType {
  geo::Point_t FlashPos;
  double FlashTime_us;
  double FlashGateTime_ns;
  bool inBeam;
  bool inGate;
  hasCRTHit Classification;
  std::vector<MatchedCRT> CRTmatches;
};

CRTMatches CRTHitmatched(
    double flashTime, geo::Point_t const& flashpos,
    std::vector<art::Ptr<sbn::crt::CRTHit>>& crtHits, double interval) {

  	std::vector<icarus::crt::CRTPMT> enteringCRTHits;
  	std::vector<icarus::crt::CRTPMT> exitingCRTHits;
  	hasCRTHit MatchType;
  	int topen = 0, topex = 0, sideen = 0, sideex = 0;
  	for (auto const& crtHit : crtHits) {
    		double tof = crtHit->ts1_ns - flashTime * 1e3;
    		double distance =
		  (flashpos - geo::Point_t{crtHit->x_pos, crtHit->y_pos, crtHit->z_pos})
            	    .R();
    		if (abs(tof) >= interval) continue;
    		if (tof < 0) {
      			if (crtHit->plane > 36)
        			sideen++;
      			else
        			topen++;
      			CRTPMT m_match = {tof, distance, crtHit};
      			enteringCRTHits.push_back(m_match);
    		} else if (tof >= 0) {
      			if (crtHit->plane > 36)
        			sideex++;
      			else
        			topex++;
      			CRTPMT m_match = {tof, distance, crtHit};
      			exitingCRTHits.push_back(m_match);
    		}
  	}
  	if (topen == 0 && sideen == 0 && topex == 0 && sideex == 0)
    		MatchType = noMatch;
  	else if (topen == 1 && sideen == 0 && topex == 0 && sideex == 0)
    		MatchType = enTop;
  	else if (topen == 0 && sideen == 1 && topex == 0 && sideex == 0)
    		MatchType = enSide;
  	else if (topen == 1 && sideen == 0 && topex == 0 && sideex == 1)
    		MatchType = enTop_exSide;
  	else if (topen == 0 && sideen == 0 && topex == 1 && sideex == 0)
    		MatchType = exTop;
  	else if (topen == 0 && sideen == 0 && topex == 0 && sideex == 1)
    		MatchType = exSide;
  	else if (topen >= 1 && sideen >= 1 && topex == 0 && sideex == 0)
    		MatchType = enTop_mult;
  	else if (topen >= 1 && sideen >= 1 && topex == 0 && sideex >= 1)
    		MatchType = enTop_exSide_mult;
  	else
    		MatchType = others;

  	return {std::move(enteringCRTHits), std::move(exitingCRTHits), MatchType};
};

  class CRTPMTMatchingProducer : public art::EDProducer {
  public:
 
    using CRTHit = sbn::crt::CRTHit;
    
    explicit CRTPMTMatchingProducer(fhicl::ParameterSet const & p);

    // The destructor generated by the compiler is fine for classes
    // without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    CRTPMTMatchingProducer(CRTPMTMatchingProducer const &) = delete;
    CRTPMTMathcingProducer(CRTPMTMatchingProducer &&) = delete;
    CRTPMTMatchingProducer & operator = (CRTPMTMatchingProducer const &) = delete; 
    CRTPMTMatchingProducer & operator = (CRTPMTMatchingProducer &&) = delete;

    // Required functions.
    void produce(art::Event & e) override;

    // Selected optional functions.
    void beginJob() override;

    void endJob() override;

    //void reconfigure(fhicl::ParameterSet const & p);

  private:

	// Params from fcl file.......

  	// art::InputTag fOpHitModuleLabel;
  	art::InputTag fOpFlashModuleLabel0;		///< OpticalFlashes Cryo 0.
  	art::InputTag fOpFlashModuleLabel1;		///< OpticalFlashes Cryo1.
	art::InputTag fCrtModuleLabel;		      	///< name of crt producer.
	art::InputTag fTriggerLabel;		        ///< name of trigger producer.
	art::InputTag fTriggerConfigurationLabel;	///< name of the trigger configuration.
    
	double fTimeOfFlightInterval;	 		///< CRTPMT time difference interval to find the match.
	int fPMTADCThresh;				///< ADC amplitude for a PMT to be considered above threshold.
	int fnOpHitToTrigger;				///< Number of OpHit above threshold to mimic the triggering PMT.
	double fBNBBeamGateMin;
	double fBNBBeamGateMax;
	double fBNBinBeamMin;
	double fBNBinBeamMax;
	double fNuMIBeamGateMin;
	double fNuMIBeamGateMax;
	double fNuMIinBeamMin;
	double fNuMIinBeamMax;

	// Trigger data product variables
	unsigned int m_gate_type;
	std::string m_gate_name;
	uint64_t m_trigger_timestamp;
	uint64_t m_gate_start_timestamp;
	uint64_t m_trigger_gate_diff;
	uint64_t m_gate_width;


	geo::GeometryCore const* fGeometryService;	///< pointer to Geometry provider

  }; // class CRTPMTMatchingProducer

  CRTPMTMatchingProducer::CRTPMTMatchingProducer(fhicl::ParameterSet const & p)
  : EDProducer{p},
      fOpFlashModuleLabel0(
          p.get<art::InputTag>("OpFlashModuleLabel0")),
      fOpFlashModuleLabel1(
          p.get<art::InputTag>("OpFlashModuleLabel1")),
      fCrtHitModuleLabel(p.get<art::InputTag>("CrtHitModuleLabel", "crthit")),
      fTriggerLabel(p.get<art::InputTag>("TriggerLabel", "daqTrigger")),
      fTriggerConfigurationLabel(
          p.get<art::InputTag>("TriggerConfiguration", "triggerconfig")),
      fnOpHitToTrigger(p.get<int>("nOpHitToTrigger")),
      fTimeOfFlightInterval(p.get<double>("TimeOfFlightInterval")),
      fPMTADCThresh(p.get<int>("PMTADCThresh")),
      fBNBBeamGateMin(p.get<double>("BNBBeamGateMin")),
      fBNBBeamGateMax(p.get<double>("BNBBeamGateMax")),
      fBNBinBeamMin(p.get<double>("BNBinBeamMin")),
      fBNBinBeamMax(p.get<double>("BNBinBeamMax")),
      fNuMIBeamGateMin(p.get<double>("NuMIBeamGateMin")),
      fNuMIBeamGateMax(p.get<double>("NuMIBeamGateMax")),
      fNuMIinBeamMin(p.get<double>("NuMIinBeamMin")),
      fNuMIinBeamMax(p.get<double>("NuMIinBeamMax")),
      fGeometryService(lar::providerFrom<geo::Geometry>()) {
	  fFlashLabels = { fOpFlashModuleLabel0, fOpFlashModuleLabel1 };
	}

  // Initialize member data here, if know don't want to reconfigure on the fly
  {
 
   // Call appropriate produces<>() functions here.
    //produces< vector<CRTHit> >();
    //produces< art::Assns<CRTHit, CRTData> >();
    
    //reconfigure(p);

  } // CRTPMTMatchingProducer()

  /*void CRTPMTMatchingProducer::reconfigure(fhicl::ParameterSet const & p)
  {
    fCrtModuleLabel = (p.get<art::InputTag> ("CrtModuleLabel")); 
    fTriggerLabel   = (p.get<art::InputTag> ("TriggerLabel")); 
  } // CRTPMTMatchingProducer::reconfigure() */

  void CRTPMTMatchingProducer::beginJob()
  {

  } // CRTPMTMatchingProducer::beginJob()
 
  void CRTPMTMatchingProducer::produce(art::Event & event)
  {
	mf::LogDebug("CRTPMTMatchingProducer: ") << "beginning analyis";

 	// add trigger info
  	auto const& triggerInfo = e.getProduct<sbn::ExtraTriggerInfo>(fTriggerLabel);
  	sbn::triggerSource bit = triggerInfo.sourceType;
  	m_gate_type = (unsigned int)bit;
  	m_gate_name = bitName(bit);
  	m_trigger_timestamp = triggerInfo.triggerTimestamp;
  	m_gate_start_timestamp = triggerInfo.beamGateTimestamp;
  	m_trigger_gate_diff = triggerInfo.triggerTimestamp - triggerInfo.beamGateTimestamp;
  	m_gate_width = fTriggerConfiguration.getGateWidth(m_gate_type);
	
	// add CRTHits
  	art::Handle<std::vector<CRTHit>> crtHitListHandle;
  	std::vector<art::Ptr<CRTHit>> crtHitList;
  	if (e.getByLabel(fCrtHitModuleLabel, crtHitListHandle))
    	  art::fill_ptr_vector(crtHitList, crtHitListHandle);

	// add optical flashes
 	for (art::InputTag const& flashLabel : fFlashLabels) {
    		auto const flashHandle =
        	  e.getHandle<std::vector<recob::OpFlash>>(flashLabel);
    		art::FindMany<recob::OpHit> findManyHits(flashHandle, e, flashLabel);

		std::vector<FlashType> thisEventFlashes;

		for (auto const& [iflash, flash] : util::enumerate(*flashHandle)) {
			double tflash = flash.Time();
      			vector<recob::OpHit const*> const& hits = findManyHits.at(iflash);
      			int nPMTsTriggering = 0;
      			double firstTime = 999999;
      			geo::vect::MiddlePointAccumulator flashCentroid;
      			double ampsum = 0, t_m = 0;
      			std::vector<double> fHitX, fHitY, fHitZ, fHitT, fHitA;
      			for (auto const& hit : hits) {
        			if (hit->Amplitude() > fPMTADCThresh) nPMTsTriggering++;
        			if (firstTime > hit->PeakTime()) firstTime = hit->PeakTime();
        			geo::Point_t const pos =
            		  	  fGeometryService->OpDetGeoFromOpChannel(hit->OpChannel())
                	  	  .GetCenter();
        			double amp = hit->Amplitude();
        			ampsum += amp;
        			fHitX.push_back(pos.X());
        			fHitY.push_back(pos.Y());
        			fHitZ.push_back(pos.Z());
        			fHitT.push_back(hit->PeakTime());
        			fHitA.push_back(amp);
        			flashCentroid.add(pos, amp);
        			t_m = t_m + hit->PeakTime();
      			}
      			geo::Point_t flash_pos = flashCentroid.middlePoint();
			t_m = t_m / nPMTsTriggering;
      			if (nPMTsTriggering < fnOpHitToTrigger) continue;

			bool inTime = flashInTime(firstTime, m_gate_type, m_trigger_gate_diff,
                                m_gate_width);
			double fThisRelGateTime = m_trigger_gate_diff + tflash * 1e3;
      			bool fThisInTime_gate = false;
     			bool fThisInTime_beam = false;
      			if (m_gate_type == 1 || m_gate_type == 3) {  // BNB OffBeamBNB
        			if (fThisRelGateTime > fBNBBeamGateMin &&
            			  fThisRelGateTime < fBNBBeamGateMax)
          			    fThisInTime_gate = true;
        			if (fThisRelGateTime > fBNBinBeamMin &&
            			  fThisRelGateTime < fBNBinBeamMax)
          			    fThisInTime_beam = true;
      			}
      			if (m_gate_type == 2 || m_gate_type == 4) {  // NuMI OffBeamNuMI
        			if (fThisRelGateTime > fNuMIBeamGateMin &&
            			  fThisRelGateTime < fNuMIBeamGateMax)
          			    fThisInTime_gate = true;
        			if (fThisRelGateTime > fNuMIinBeamMin &&
           			  fThisRelGateTime < fNuMIinBeamMax)
          			    fThisInTime_beam = true;
      			}
      			inTime = fThisInTime_gate;
		        icarus::crt::CRTMatches CRTmatches = CRTHitmatched(
          		  firstTime, flash_pos, crtHitList, fTimeOfFlightInterval);
			int TopEn = 0, TopEx = 0, SideEn = 0, SideEx = 0;
      			auto nCRTHits = CRTmatches.entering.size() + CRTmatches.exiting.size();
			eventType = CRTmatches.matchType;
			if (nCRTHits > 0) {
        			for (auto const& entering : CRTmatches.entering) {
          				vector<double> CRTpos {entering.CRTHit->x_pos,
                                 		entering.CRTHit->y_pos,
                                 		entering.CRTHit->z_pos};
          				geo::Point_t thisCRTpos {entering.CRTHit->x_pos,
                                   		entering.CRTHit->y_pos,
                                   		entering.CRTHit->z_pos};
          				double CRTtime = entering.CRTHit->ts1_ns / 1e3;
          				int CRTRegion = entering.CRTHit->plane;
          				int CRTSys = 0;
          				if (CRTRegion >= 36) CRTSys = 1;  	// Very lazy way to determine if the Hit is a Top or a Side.
                    					    			// Will update it when bottom CRT will be availble.
          				double CRTTof_opflash = entering.CRTHit->ts1_ns - tflash * 1e3;
          				std::vector<int> HitFebs;
          				for (auto crts : entering.CRTHit->feb_id) {
						HitFebs.emplace_back((int)crts);
					}
          				if (CRTSys == 0) TopEn++;
          				if (CRTSys == 1) SideEn++;		
					MatchedCRT thisCRTMatch = { /* .CRTHitPos = */ thisCRTpos, // C++20: restore initializers
                                     	  /* .CRTPMTTimeDiff_ns = */ CRTTof_opflash,
                                     	  /* .CRTTime_us = */ CRTtime,
                                     	  /* .CRTSys = */ CRTSys,
                                     	  /* .CRTRegion = */ CRTRegion};	
					thisFlashCRTmatches.push_back(thisCRTMatch);
        			}
				for (auto const& exiting : CRTmatches.exiting) {
          				vector<double> CRTpos {exiting.CRTHit->x_pos,
                                 		exiting.CRTHit->y_pos,
                                 		exiting.CRTHit->z_pos};
          				geo::Point_t thisCRTpos {exiting.CRTHit->x_pos,
                                   		enxiting.CRTHit->y_pos,
                                   		exiting.CRTHit->z_pos};
          				double CRTtime = exiting.CRTHit->ts1_ns / 1e3;
          				int CRTRegion = exiting.CRTHit->plane;
          				int CRTSys = 0;
          				if (CRTRegion >= 36) CRTSys = 1; 
          				double CRTTof_opflash = exiting.CRTHit->ts1_ns - tflash * 1e3;
          				std::vector<int> HitFebs;
          				for (auto crts : exiting.CRTHit->feb_id) {
						HitFebs.emplace_back((int)crts);
					}
          				if (CRTSys == 0) TopEx++;
          				if (CRTSys == 1) SideEx++;		
					MatchedCRT thisCRTMatch = { /* .CRTHitPos = */ thisCRTpos, // C++20: restore initializers
                                     	  /* .CRTPMTTimeDiff_ns = */ CRTTof_opflash,
                                     	  /* .CRTTime_us = */ CRTtime,
                                     	  /* .CRTSys = */ CRTSys,
                                     	  /* .CRTRegion = */ CRTRegion};	
					thisFlashCRTmatches.push_back(thisCRTMatch);
        			}
			} //Fine CRT
		} // Fine di questo flash
      }
    //std::unique_ptr< vector<CRTHit> > CRTHitcol( new vector<CRTHit>);
    //std::unique_ptr< art::Assns<CRTHit, CRTData> > Hitassn( new art::Assns<CRTHit, CRTData>);
    //art::PtrMaker<sbn::crt::CRTHit> makeHitPtr(event);

    if (event.getByLabel(fCrtModuleLabel, crtListHandle))
      art::fill_ptr_vector(crtList, crtListHandle);

    //add trigger info
    if( !fTriggerLabel.empty() ) {

      art::Handle<sbn::ExtraTriggerInfo> trigger_handle;
      event.getByLabel( fTriggerLabel, trigger_handle );
      if( trigger_handle.isValid() )
      	m_trigger_timestamp = trigger_handle->triggerTimestamp; 
      else
      	mf::LogError("CRTPMTMatchingProducer") << "No raw::Trigger associated to label: " << fTriggerLabel.label() << "\n" ;
    } else{ 
      std::cout  << "Trigger Data product " << fTriggerLabel.label() << " not found!\n" ;
    }

    mf::LogInfo("CRTPMTMatchingProducer")
      <<"Number of SiPM hits = "<<crtList.size();
    //m_trigger_timestamp = event.getProduct<sbn::ExtraTriggerInfo>(fTriggerLabel).triggerTimestamp;

    mf::LogInfo("CRTPMTMatchingProducer")
      <<"Number of SiPM hits = "<<crtList.size();

    vector<art::Ptr<CRTData>> crtData = hitAlg.PreselectCRTData(crtList, m_trigger_timestamp);

    vector<std::pair<CRTHit, vector<int>>> crtHitPairs = hitAlg.CreateCRTHits(crtData, m_trigger_timestamp);
    //vector<std::pair<CRTHit, vector<int>>> crtHitPairs = hitAlg.CreateCRTHits(crtList);

    mf::LogInfo("CRTPMTMatchingProducer")
      << "Number of CRTHit,data indices pairs = " << crtHitPairs.size();

    for(auto const& crtHitPair : crtHitPairs){

      CRTHitcol->push_back(crtHitPair.first);
      art::Ptr<CRTHit> hitPtr = makeHitPtr(CRTHitcol->size()-1);
      nHits++;

      for(auto const& data_i : crtHitPair.second){

        Hitassn->addSingle(hitPtr, crtList[data_i]);
      }
    }
      
    event.put(std::move(CRTHitcol));
    event.put(std::move(Hitassn));

    mf::LogInfo("CRTSimHitProducer")
      <<"Number of CRT hits produced = "<<nHits;

  } // CRTPMTMatchingProducer::produce()

  void CRTPMTMatchingProducer::endJob()
  {

  } // CRTPMTMatchingProducer::endJob()

  DEFINE_ART_MODULE(CRTPMTMatchingProducer)

}
} //end namespace
