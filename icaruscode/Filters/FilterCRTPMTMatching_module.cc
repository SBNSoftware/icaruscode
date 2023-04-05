////////////////////////////////////////////////////////////////////////
// Class:       FilterCRTPMTMatching
// Plugin Type: analyzer (Unknown Unknown)
// File:        FilterCRTPMTMatching_module.cc
//
// Generated at Thu Feb 17 13:31:00 2022 by Biswaranjan Behera using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindMany.h"


// LArSoft includes
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

// Data product includes

#include "icaruscode/CRT/CRTUtils/CRTCommonUtils.h"
#include "icaruscode/Decode/DataProducts/TriggerConfiguration.h"
#include "larcorealg/CoreUtils/enumerate.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "sbnobj/Common/CRT/CRTHit.hh"
#include "sbnobj/Common/Trigger/ExtraTriggerInfo.h"

// C++ includes
#include <vector>
#include <map>
#include <numeric>

// ROOT includes
#include "TTree.h"
#include "TVector3.h"

using std::map;
using std::vector;

namespace icarus
{
    namespace crt
    {
        class FilterCRTPMTMatching;
        struct CRTPMT
        {
            double tof;
	    double distance;
            art::Ptr<sbn::crt::CRTHit> CRTHit;
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
            *  * `7`: matched with multiple Top CRT hits before the optical flash and more then 1 side CRT hits after the optical flash
            *  * `8`: all the other cases
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
        struct MatchedCRT
        {
	    // vector of pairs where first is the matched Time of Flight and second is the
  	    // matched CRTHit
            std::vector<CRTPMT> entering;
            std::vector<CRTPMT> exiting;
	    hasCRTHit matchType;
        };
        struct TrajPoint
        {
            double X;
            double Y;
            double Z;
            double T;
        };
    }
}

using namespace icarus::crt;

bool flashInTime(double const &flashTime, int gateType, double gateDiff, double gateWidth)
{
    //   {1, 2.2},  // BNB
    //   {2, 10.1}, // NuMI
    //   {3, 2.2},  // BNB Offbeam
    //   {4, 10.1}, // NuMI Offbeam
   

    //As a reminder, I will leave here the commented part of the vetoOffset, in case something changes in the future
    /*double vetoOffset = 3.5;*/
    double activeGate = gateWidth /*- vetoOffset*/;

    double relFlashTime = flashTime + gateDiff / 1000. /*- vetoOffset*/;
    mf::LogInfo("FilterCRTPMTMatching::flashInTime") << "Gate Diff " << gateDiff / 1000 << " Ftime+gateDiff " << flashTime + gateDiff / 1000. << " " << activeGate;

    return ((relFlashTime > 0) && (relFlashTime < activeGate));
}

icarus::crt::MatchedCRT CRTHitmatched(const double& flashTime, geo::Point_t const& flashpos,
    std::vector<art::Ptr<sbn::crt::CRTHit>>& crtHits, const double& interval) {
    std::vector<icarus::crt::CRTPMT> enteringCRTHits;
    std::vector<icarus::crt::CRTPMT> exitingCRTHits;
    hasCRTHit MatchType;
    int topen = 0, topex = 0, sideen = 0, sideex = 0;
    for (auto const crtHit : crtHits)
    {
        double tof = crtHit->ts1_ns / 1e3 - flashTime;
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
}

class icarus::crt::FilterCRTPMTMatching : public art::EDFilter
{
public:
    using CRTHit = sbn::crt::CRTHit;

    explicit FilterCRTPMTMatching(fhicl::ParameterSet const& p);

    // Required functions.
    void getTriggerConf(art::Run &);
    bool filter(art::Event &) override;

private:
    // Declare member data here.

    static bool HitCompare(const art::Ptr<CRTHit> &h1, const art::Ptr<CRTHit> &h2);
    void ClearVecs();

    art::InputTag fOpFlashModuleLabel0;
    art::InputTag fOpFlashModuleLabel1;
    art::InputTag fOpFlashModuleLabel2;
    art::InputTag fOpFlashModuleLabel3;
    art::InputTag fCrtHitModuleLabel;
    art::InputTag fTriggerLabel;
    art::InputTag fTriggerConfigurationLabel;

    icarus::TriggerConfiguration fTriggerConfiguration;
    
    //double                     fFlashTimeCut;

    int                        fEvent;  ///< number of the event being processed.
    int                        fRun;    ///< number of the run being processed.
    int                        fSubRun; ///< number of the sub-run being processed.

    // add trigger data product vars
    unsigned int               m_gate_type;
    std::string                m_gate_name;
    uint64_t                   m_trigger_timestamp;
    uint64_t                   m_gate_start_timestamp;
    uint64_t                   m_trigger_gate_diff;
    uint64_t		       m_gate_width;

    std::string		       fFilterLevel;		// Filter level, default is loose 
    int                        fnOpHitToTrigger;        // Number of OpHit above threshold to mimic the trigger
    double 		       fTimeOfFlightInterval;   // Time of Flight interval to find the match
    bool                       fOutputTree;             // Output tree or not
    double fCoinWindow;
    double fOpDelay;
    double fCrtDelay;
    int fPMTADCThresh;
    int fFlashPeThresh;
    int fHitPeThresh;
    double fFlashVelocity;
    double fFlashZOffset;
    double fHitVelocityMax;
    double fHitVelocityMin;
    std::vector<int> fTopBefore;
    std::vector<int> fTopAfter;
    std::vector<int> fSideBefore;
    std::vector<int> fSideAfter;
    double fBNBBeamGateMin;
    double fBNBBeamGateMax;
    double fBNBinBeamMin;
    double fBNBinBeamMax;
    double fNuMIBeamGateMin;
    double fNuMIBeamGateMax;
    double fNuMIinBeamMin;
    double fNuMIinBeamMax;

    int  		       fFiltered;      ///< Was the event filtered?

    std::vector<art::InputTag> fFlashLabels;

    std::vector<int>           fHitsInTime;   ///< Classification of the Hits in Time

    TTree*                     fMatchTree;

    // matchTree vars
    vector<double>             fOpFlashX;     ///< X position of the optical Flash Barycenter.
    vector<double>             fOpFlashY;     ///< Y position of the optical Flash Barycenter.
    vector<double>             fOpFlashZ;     ///< Z position of the optical Flash Barycenter.
    vector<double>                     fOpFlashTime_us;         ///< Time of the optical Flash w.r.t. Global Trigger in us.
    vector<double>                     fOpFlashTimehit_us;				///< Time of the first optical hit w.r.t. Global Trigger in us.
    vector<double>                     fOpFlashTimeAbs;				///< Time of the optical Flash in absolute timestestamp.
    vector<bool>                       fInTime;                 ///< Was the OpFlash in beam spill time?
    vector<bool> 		       fInTime_gate;		///< Was the OpFlash inside the beam gate window?
    vector<bool>		       fInTime_beam;		///< Was the OpFlash inside the true beam gate window?
    vector<double>                     fOpFlashPE;              ///< Total PEs of the optical Flash.
    vector<vector<double>>             fOpHitX;                 ///< X position of the OP hit.
    vector<vector<double>>             fOpHitY;                 ///< Y position of the OP hit.
    vector<vector<double>>             fOpHitZ;                 ///< Z position of the OP hit.
    vector<vector<double>>             fOpHitT;                 ///< T of the OP hit in us.
    vector<vector<double>>             fOpHitA;                 ///< Amplitude (ADC) of the OP hit.
    vector<int>                        fNCRTmatch;		///< Number of CRT Hits matched.
    vector<vector<vector<double>>>     fMatchedCRTpos;          ///< Position of the matched CRTs.
    vector<vector<double>>             fMatchedCRTtime_us;      ///< Time of the matched CRT Hits w.r.t. Global Trigger in us.
    vector<vector<ULong64_t>>	       fMatchedCRTtime_abs;     ///< Time of the matched CRT Hits in absolute timestamp.
    vector<vector<int>>                fMatchedCRTregion;       ///< Region of the matched CRT Hits.
    vector<vector<vector<int>>>        fMatchedCRTmodID;        ///< ModID of the CRT module.
    vector<vector<int>>                fMatchedCRTsys;          ///< Subsystem of the matched CRT Hit.
    vector<vector<double>>             fMatchedCRTamplitude;    ///< Amplitude of the matched CRT Hit.
    vector<vector<double>>	       fDistance;  		///< Distance between matched CRT and light barycenter.
    vector<vector<double>>	       fTofOpHit;  		///< Time difference [ns] between matched CRT and first Optical Hit.
    vector<vector<double>>             fTofOpFlash;             ///< Time difference [ns] between matched CRT and Optical Flash.
    vector<vector<double>>             fCRTGateDiff;            ///< Difference between CRTHit and BeamGate opening.
    vector<int>                        fFlashType;              ///< Classification of the flashes.
    vector<double>		       fRelGateTime;		///< Time difference between optical Flash and Beam Gate opening.
    geo::GeometryCore const*   fGeometryService; ///< pointer to Geometry provider.

};

icarus::crt::FilterCRTPMTMatching::FilterCRTPMTMatching(fhicl::ParameterSet const &p) : EDFilter{p} 
    ,fOpFlashModuleLabel0(p.get<art::InputTag>("OpFlashModuleLabel0", "opflashTPC0"))
    ,fOpFlashModuleLabel1(p.get<art::InputTag>("OpFlashModuleLabel1", "opflashTPC1"))
    //,fOpFlashModuleLabelVec(p.get<std::vector<art::InputTag>>("OpFlashModuleLabelVec", {"opflashE","opflashW"}))
    ,fCrtHitModuleLabel(p.get<art::InputTag>("CrtHitModuleLabel", "crthit"))
    ,fTriggerLabel(p.get<art::InputTag>("TriggerLabel", "daqTrigger")) 
    ,fTriggerConfigurationLabel(p.get<art::InputTag>("TriggerConfiguration", "triggerconfig"))
    ,fFilterLevel(p.get<std::string>("FilterLevel","FilterOption"))
    ,fnOpHitToTrigger(p.get<int>("nOpHitToTrigger", 0))
    ,fTimeOfFlightInterval(p.get<double>("TimeOfFlightInterval", 0.1))
    ,fOpDelay(p.get<double>("OpDelay", 55.1))
    ,fCrtDelay(p.get<double>("CrtDelay", 1.6e6))
    ,fPMTADCThresh(p.get<int>("PMTADCThresh", 400))
    ,fFlashPeThresh(p.get<int>("FlashPeThresh", 9000))
    ,fHitPeThresh(p.get<int>("HitPeThresh", 700))
    ,fBNBBeamGateMin(p.get<double>("BNBBeamGateMin", -550))
    ,fBNBBeamGateMax(p.get<double>("BNBBeamGateMax", 2300))
    ,fBNBinBeamMin(p.get<double>("BNBinBeamMin", -300))
    ,fBNBinBeamMax(p.get<double>("BNBinBeamMax", 1300))
    ,fNuMIBeamGateMin(p.get<double>("NuMIBeamGateMin", -550))
    ,fNuMIBeamGateMax(p.get<double>("NuMIBeamGateMax", 10000))
    ,fNuMIinBeamMin(p.get<double>("NuMIinBeamMin", -300))
    ,fNuMIinBeamMax(p.get<double>("NuMIinBeamMax", 9100))
    ,fGeometryService(lar::providerFrom<geo::Geometry>())
{
    fFlashLabels.push_back(fOpFlashModuleLabel0);
    fFlashLabels.push_back(fOpFlashModuleLabel1);
    art::ServiceHandle<art::TFileService> tfs;
    if (fOutputTree==true) {

    fMatchTree = tfs->make<TTree>("matchTree", "CRTHit - OpHit/Flash matching analysis");
    fMatchTree->Branch("event", &fEvent, "event/I");
    fMatchTree->Branch("run", &fRun, "run/I");
    fMatchTree->Branch("subrun", &fSubRun, "subrun/I");
    fMatchTree->Branch("Filter", &fFiltered);
    fMatchTree->Branch("InTimeFlashes_classification", &fHitsInTime);
    fMatchTree->Branch("opFlash_X", &fOpFlashX);
    fMatchTree->Branch("opFlash_Y", &fOpFlashY);
    fMatchTree->Branch("opFlash_Z", &fOpFlashZ);
    fMatchTree->Branch("opFlash_time_us", &fOpFlashTime_us);
    fMatchTree->Branch("opFlashAbsTime", &fOpFlashTimeAbs);
    fMatchTree->Branch("opFlash_time_firsthit_us", &fOpFlashTimehit_us);
    fMatchTree->Branch("opHit_x", &fOpHitX);
    fMatchTree->Branch("opHit_y", &fOpHitY);
    fMatchTree->Branch("opHit_z", &fOpHitZ);
    fMatchTree->Branch("opHit_t", &fOpHitT);
    fMatchTree->Branch("opHit_amplitude", &fOpHitA);
    fMatchTree->Branch("inTime", &fInTime);
    fMatchTree->Branch("inTime_beam", &fInTime_beam);
    fMatchTree->Branch("inTime_gate", &fInTime_gate);
    fMatchTree->Branch("FlashType", &fFlashType);
    fMatchTree->Branch("CRT_matches", &fNCRTmatch);
    fMatchTree->Branch("CRT_pos", &fMatchedCRTpos);
    fMatchTree->Branch("CRT_Time_us", &fMatchedCRTtime_us);
    fMatchTree->Branch("CRT_absTime", &fMatchedCRTtime_abs);
    fMatchTree->Branch("CRT_region", &fMatchedCRTregion);
    fMatchTree->Branch("CRT_FEB", &fMatchedCRTmodID);
    fMatchTree->Branch("CRT_subsys", &fMatchedCRTsys);
    fMatchTree->Branch("CRT_pes", &fMatchedCRTamplitude);
    fMatchTree->Branch("distance", &fDistance);
    fMatchTree->Branch("TOF_opflash", &fTofOpFlash);
    fMatchTree->Branch("TOF_ophit", &fTofOpHit);
    fMatchTree->Branch("TopHitsBefore", &fTopBefore);
    fMatchTree->Branch("TopHitsAfter", &fTopAfter);
    fMatchTree->Branch("SideHitsBefore", &fSideBefore);
    fMatchTree->Branch("SideHitsAfter", &fSideAfter);
    fMatchTree->Branch("CRT_gate_diff", &fCRTGateDiff);
    fMatchTree->Branch("RelativeGateTime", &fRelGateTime);
    fMatchTree->Branch("gate_type", &m_gate_type, "gate_type/b");
    fMatchTree->Branch("gate_name", &m_gate_name);
    fMatchTree->Branch("trigger_timestamp", &m_trigger_timestamp, "trigger_timestamp/l");
    fMatchTree->Branch("gate_start_timestamp", &m_gate_start_timestamp, "gate_start_timestamp/l");
    fMatchTree->Branch("trigger_gate_diff", &m_trigger_gate_diff, "trigger_gate_diff/l");
    }
}

void icarus::crt::FilterCRTPMTMatching::getTriggerConf(art::Run &r) {

  fTriggerConfiguration =
    r.getProduct<icarus::TriggerConfiguration>(fTriggerConfigurationLabel);

}

bool icarus::crt::FilterCRTPMTMatching::filter(art::Event &e)
{
    mf::LogDebug("FilterCRTPMTMatching: ") << "beginning analyis";

    // Start by fetching some basic event information for our n-tuple.
    fEvent = e.id().event();
    fRun = e.run();
    fSubRun = e.subRun();

    ClearVecs();
    // add trigger info
    if (!fTriggerLabel.empty()) {
        art::Handle<sbn::ExtraTriggerInfo> trigger_handle;
        e.getByLabel(fTriggerLabel, trigger_handle);
        if (trigger_handle.isValid()) {
            sbn::triggerSource bit = trigger_handle->sourceType;
            m_gate_type = (unsigned int)bit;
            m_gate_name = bitName(bit);
            m_trigger_timestamp = trigger_handle->triggerTimestamp;
            m_gate_start_timestamp = trigger_handle->beamGateTimestamp;
            m_trigger_gate_diff = trigger_handle->triggerTimestamp - trigger_handle->beamGateTimestamp;
            m_gate_width = fTriggerConfiguration.getGateWidth(m_gate_type);
        } else {
            mf::LogError("CRTPMTMatching:")
                << "No sbn::ExtraTriggerInfo associated to label: "
                << fTriggerLabel.encode() << "\n";
     	}
    }
    std::array<art::Handle<std::vector<recob::OpFlash>>, 2U> flashHandles;
    for (int i = 0; i < 2; i++) {
        flashHandles[i] = e.getHandle<std::vector<recob::OpFlash>>(fFlashLabels[i]);
    }
    // CRTHits
    art::Handle<std::vector<CRTHit>> crtHitListHandle;
    std::vector<art::Ptr<CRTHit>> crtHitList;
    if (e.getByLabel(fCrtHitModuleLabel, crtHitListHandle))
        art::fill_ptr_vector(crtHitList, crtHitListHandle);
    //hasCRTHit Type = others;
    //bool Cosmic = false;
    if (fFilterLevel == "loose") fFilterLevel = "loose";
    else if (fFilterLevel == "medium") fFilterLevel = "medium";
    else if (fFilterLevel == "thight") fFilterLevel = "tight";
    else fFilterLevel = "loose";
    mf::LogInfo("FilterCRTPMTMatching::FilteringLevel ") << fFilterLevel;
    for (art::InputTag const& flashLabel : fFlashLabels) {
        auto const flashHandle = e.getHandle<std::vector<recob::OpFlash>>(flashLabel);
        art::FindMany<recob::OpHit> findManyHits(flashHandle, e, flashLabel);

        for (auto const& [iflash, flash] : util::enumerate(*flashHandle)) {
            hasCRTHit eventType = others;
            double tflash = flash.Time();
            double tAbsflash = flash.AbsTime();
            vector<recob::OpHit const*> const& hits = findManyHits.at(iflash);
            int nPMTsTriggering = 0;
            double firstTime = 999999;
            geo::vect::MiddlePointAccumulator flashCentroid;
            double ampsum = 0, t_m = 0;
	    std::vector<double> fHitX, fHitY, fHitZ, fHitT, fHitA;
            for (auto const& hit : hits) {
                if (hit->Amplitude() > fPMTADCThresh) nPMTsTriggering++;
                if (firstTime > hit->PeakTime()) firstTime = hit->PeakTime();
                geo::Point_t const pos = fGeometryService->OpDetGeoFromOpChannel(hit->OpChannel()).GetCenter();
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
            if (nPMTsTriggering < fnOpHitToTrigger) {
                continue;
            }
	    fOpHitX.push_back(fHitX);
	    fOpHitY.push_back(fHitY);
	    fOpHitZ.push_back(fHitZ);
	    fOpHitT.push_back(fHitT);
	    fOpHitA.push_back(fHitA);
            bool inTime = flashInTime(firstTime, m_gate_type, m_trigger_gate_diff, m_gate_width);
            double fThisRelGateTime = m_trigger_gate_diff + (tAbsflash - 1500) * 1e3;
            bool fThisInTime_gate=false;
            bool fThisInTime_beam=false;
            if (m_gate_type == 1 || m_gate_type == 3) {  // BNB OffBeamBNB
                if (fThisRelGateTime > fBNBBeamGateMin && fThisRelGateTime < fBNBBeamGateMax)
                    fThisInTime_gate=true;
                if (fThisRelGateTime > fBNBinBeamMin && fThisRelGateTime < fBNBinBeamMax)
                    fThisInTime_beam=true;
            }
            if (m_gate_type == 2 || m_gate_type == 4) {  // NuMI OffBeamNuMI
                if (fThisRelGateTime > fNuMIBeamGateMin && fThisRelGateTime < fNuMIBeamGateMax)
                    fThisInTime_gate=true;
                if (fThisRelGateTime > fNuMIinBeamMin && fThisRelGateTime < fNuMIinBeamMax)
                    fThisInTime_beam=true;
            }
	    fRelGateTime.push_back(fThisRelGateTime);
            inTime = fThisInTime_gate;
            fInTime_gate.push_back(fThisInTime_gate);
	    fInTime_beam.push_back(fThisInTime_beam);
            icarus::crt::MatchedCRT CRTmatches =
              CRTHitmatched(firstTime, flash_pos, crtHitList, fTimeOfFlightInterval);
            int TopEn = 0, TopEx = 0, SideEn = 0, SideEx = 0;
            auto nCRTHits = CRTmatches.entering.size() + CRTmatches.exiting.size();
	    double peflash = flash.TotalPE();
	    std::vector<std::vector<int>> fThisMatchedCRTmodID;
	    std::vector<int> fThisMatchedCRTregion, fThisMatchedCRTsys;
	    std::vector<double> fThisMatchedCRTtime_us, fThisTofOpFlash, fThisMatchedCRTamplitude, fThisDistance, fThisTofOpHit, fThisCRTGateDiff;
	    std::vector<vector<double>> fThisMatchedCRTpos;
	    std::vector<ULong64_t> fThisMatchedCRTtime_abs;
            if (nCRTHits >= 0) {
                eventType = CRTmatches.matchType;
                for (auto& entering : CRTmatches.entering) {
                    vector<double> CRTpos = {entering.CRTHit->x_pos, entering.CRTHit->y_pos, entering.CRTHit->z_pos};
                    double CRTtime = entering.CRTHit->ts1_ns / 1e3;
                    ULong64_t CRTAbsTime = entering.CRTHit->ts0_s;
                    int CRTRegion = entering.CRTHit->plane;
                    int CRTSys = 0;
                    if (CRTRegion >= 36) CRTSys = 1;  // Very lazy way to determine if the Hit is a Top or a Side.
                        // Will update it when bottom CRT will be availble.
                    double CRTPe = entering.CRTHit->peshit;
                    double CRTDistance = entering.distance;
                    double CRTTof_ophit = entering.tof;
                    double CRTTof_opflash = CRTtime - tflash;
                    double CRTGateDiff = CRTAbsTime - m_gate_start_timestamp;
                    std::vector<int> HitFebs;
                    for (auto& crts : entering.CRTHit->feb_id) {
                        HitFebs.emplace_back((int)crts);
                    }
                    if (CRTSys == 0) TopEn++;
                    if (CRTSys == 1) SideEn++;
                    fThisMatchedCRTmodID.emplace_back(HitFebs);
                    fThisMatchedCRTpos.emplace_back(CRTpos);
                    fThisMatchedCRTtime_us.emplace_back(CRTtime);
                    fThisTofOpFlash.emplace_back(CRTTof_opflash);
                    fThisMatchedCRTtime_abs.emplace_back(CRTAbsTime);
                    fThisMatchedCRTregion.emplace_back(CRTRegion);
                    fThisMatchedCRTsys.emplace_back(CRTSys);
                    fThisMatchedCRTamplitude.emplace_back(CRTPe);
                    fThisDistance.emplace_back(CRTDistance);
                    fThisTofOpHit.emplace_back(CRTTof_ophit);
                    fThisCRTGateDiff.emplace_back(CRTGateDiff);
                    HitFebs.clear();
       	        }
                for (auto& exiting : CRTmatches.exiting) {
                    vector<double> CRTpos = {exiting.CRTHit->x_pos, exiting.CRTHit->y_pos, exiting.CRTHit->z_pos};
                    double CRTtime = exiting.CRTHit->ts1_ns / 1e3;
                    ULong64_t CRTAbsTime = exiting.CRTHit->ts0_s;
                    int CRTRegion = exiting.CRTHit->plane;
                    int CRTSys = 0;
                    if (CRTRegion >= 36) CRTSys = 1;
                    double CRTPe = exiting.CRTHit->peshit;
                    double CRTDistance = exiting.distance;
                    double CRTTof_ophit = exiting.tof;
                    double CRTTof_opflash = CRTtime - tflash;
                    double CRTGateDiff = CRTAbsTime - m_gate_start_timestamp;
                    std::vector<int> HitFebs;
                    for (auto& crts : exiting.CRTHit->feb_id) {
                        HitFebs.emplace_back((int)crts);
                    }
                    if (CRTSys == 0) TopEx++;
                    if (CRTSys == 1) SideEx++;
                    fThisMatchedCRTmodID.emplace_back(HitFebs);
                    fThisMatchedCRTpos.push_back(std::move(CRTpos));
                    fThisMatchedCRTtime_us.emplace_back(CRTtime);
                    fThisMatchedCRTtime_abs.emplace_back(CRTAbsTime);
                    fThisMatchedCRTregion.emplace_back(CRTRegion);
                    fThisMatchedCRTsys.emplace_back(CRTSys);
                    fThisMatchedCRTamplitude.emplace_back(CRTPe);
                    fThisDistance.emplace_back(CRTDistance);
                    fThisTofOpFlash.emplace_back(CRTTof_opflash);
                    fThisTofOpHit.emplace_back(CRTTof_ophit);
                    fThisCRTGateDiff.emplace_back(CRTGateDiff);
        	}
		
                //if ((eventType == enTop || eventType == enTop_exSide ||  eventType == enTop_mult || eventType == enTop_exSide_mult) && inTime == true) Cosmic = true;
                if (inTime == true) fHitsInTime.push_back(eventType);
            }
	    fMatchedCRTmodID.push_back(fThisMatchedCRTmodID);
            fMatchedCRTpos.push_back(fThisMatchedCRTpos);
            fMatchedCRTtime_us.emplace_back(fThisMatchedCRTtime_us);
            fMatchedCRTtime_abs.emplace_back(fThisMatchedCRTtime_abs);
            fMatchedCRTregion.emplace_back(fThisMatchedCRTregion);
            fMatchedCRTsys.emplace_back(fThisMatchedCRTsys);
            fMatchedCRTamplitude.emplace_back(fThisMatchedCRTamplitude);
            fDistance.emplace_back( fThisDistance);
            fTofOpFlash.emplace_back(fThisTofOpFlash);
            fTofOpHit.emplace_back(fThisTofOpHit);
            fCRTGateDiff.emplace_back( fThisCRTGateDiff);
            fTopBefore.push_back(TopEn);
            fTopAfter.push_back(TopEx);
            fSideBefore.push_back(SideEn);
            fSideAfter.push_back(SideEx);
            fOpFlashPE.push_back(peflash);
            fOpFlashX.push_back(flash_pos.X());
            fOpFlashY.push_back(flash_pos.Y());
            fOpFlashZ.push_back(flash_pos.Z());
            fOpFlashTime_us.push_back(tflash);
            fOpFlashTimeAbs.push_back(tAbsflash);	
            fOpFlashTimehit_us.push_back(firstTime);
            fInTime.push_back(inTime);
            fFlashType.push_back(eventType);
            fNCRTmatch.push_back(nCRTHits);
        }
    }
    bool HitsInGate=false;
    if (fFilterLevel == "loose") HitsInGate = false;  // By default, with loose Filtering, everything is "intersting", nothing tagged as clear cosmic
    else if (fFilterLevel == "medium") {
	HitsInGate=true;
	for (const auto& h :  fHitsInTime){
	    bool ThisFlashType=false;
       	    if (h == enTop || h == enTop_exSide || h == enTop_mult || h == enTop_exSide_mult) ThisFlashType = true;
	    // With Medium filter, everything (inTime) which is associated with Top CRT Hit before the Flash is filtered as clear cosmic
	    else ThisFlashType = false;
	    HitsInGate = HitsInGate&&ThisFlashType;
	}
    }
    else if (fFilterLevel == "thight"){
	//if (fEventType != 0 || fEventType != 4 || fEventType != 5) Cosmic = true;
	// With Thight filter, everything (inTime) associated with a CRT Hit before the Flash is filtered as clear cosmic
	//else Cosmic = false;
	HitsInGate=true;
	for (const auto& h :  fHitsInTime){
	    bool ThisFlashType=false;
       	    if (h != noMatch || h != exTop || h != exSide ) ThisFlashType = true;
	    // With Medium filter, everything (inTime) which is associated with Top CRT Hit before the Flash is filtered as clear cosmic
	    else ThisFlashType = false;
	    HitsInGate = HitsInGate&&ThisFlashType;
	}
    }
    fFiltered=HitsInGate;
    //feventType=fHitsInTime[0];
    fMatchTree->Fill();
    //fFiltered.clear();
    fHitsInTime.clear();
    fTopBefore.clear();
    fTopAfter.clear();
    fSideBefore.clear();
    fSideAfter.clear();
    fOpFlashPE.clear();
    fOpFlashTime_us.clear();
    fOpFlashTimeAbs.clear();
    fOpFlashTimehit_us.clear();
    fInTime.clear();
    fFlashType.clear();
    fNCRTmatch.clear();
    fOpHitX.clear();
    fOpHitY.clear();
    fOpHitZ.clear();
    fOpHitT.clear();
    fOpHitA.clear();
    fOpFlashX.clear();
    fOpFlashY.clear();
    fOpFlashZ.clear();
    fMatchedCRTmodID.clear();
    fMatchedCRTpos.clear();
    fMatchedCRTtime_us.clear();
    fMatchedCRTtime_abs.clear();
    fMatchedCRTregion.clear();
    fMatchedCRTsys.clear();
    fMatchedCRTamplitude.clear();
    fDistance.clear();
    fTofOpHit.clear();
    fTofOpFlash.clear();
    fCRTGateDiff.clear();

    return HitsInGate;
}

void icarus::crt::FilterCRTPMTMatching::ClearVecs()
{
    // matchTree
    //fFiltered.clear();
    fHitsInTime.clear();
    fTopBefore.clear();
    fTopAfter.clear();
    fSideBefore.clear();
    fSideAfter.clear();
    fOpFlashPE.clear();
    fOpFlashTime_us.clear();
    fOpFlashTimeAbs.clear();
    fOpFlashTimehit_us.clear();
    fInTime.clear();
    fFlashType.clear();
    fNCRTmatch.clear();
    fOpHitX.clear();
    fOpHitY.clear();
    fOpHitZ.clear();
    fOpHitT.clear();
    fOpHitA.clear();
    fOpFlashX.clear();
    fOpFlashY.clear();
    fOpFlashZ.clear();
    fMatchedCRTmodID.clear();
    fMatchedCRTpos.clear();
    fMatchedCRTtime_us.clear();
    fMatchedCRTtime_abs.clear();
    fMatchedCRTregion.clear();
    fMatchedCRTsys.clear();
    fMatchedCRTamplitude.clear();
    fDistance.clear();
    fTofOpHit.clear();
    fCRTGateDiff.clear();
}

DEFINE_ART_MODULE(icarus::crt::FilterCRTPMTMatching)
