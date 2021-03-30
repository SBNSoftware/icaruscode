#include "icaruscode/Analysis/tools/ITrackHistogramTool.h"

#include "fhiclcpp/ParameterSet.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileDirectory.h"
#include "canvas/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "larcore/Geometry/Geometry.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RawData/RawDigit.h"

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TProfile2D.h"

#include <cmath>
#include <algorithm>

namespace BasicTrackAnalysis
{
    ////////////////////////////////////////////////////////////////////////
    //
    // Class:       BasicTrackAnalysis
    // Module Type: producer
    // File:        BasicTrackAnalysis.h
    //
    //              The intent of this module is to provide methods for
    //              "analyzing" hits on waveforms
    //
    // Configuration parameters:
    //
    // TruncMeanFraction     - the fraction of waveform bins to discard when
    //
    // Created by Tracy Usher (usher@slac.stanford.edu) on February 19, 2016
    //
    ////////////////////////////////////////////////////////////////////////
    
// The following typedefs will, obviously, be useful
using HitPtrVec       = std::vector<art::Ptr<recob::Hit>>;
using ViewHitMap      = std::map<size_t,HitPtrVec>;
using TrackViewHitMap = std::map<int,ViewHitMap>;

class BasicTrackAnalysis : virtual public ITrackHistogramTool
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    explicit BasicTrackAnalysis(fhicl::ParameterSet const & pset);
    
    /**
     *  @brief  Destructor
     */
    ~BasicTrackAnalysis();
    
    // provide for initialization
    void configure(fhicl::ParameterSet const & pset) override;

    /**
     *  @brief Interface for initializing the histograms to be filled
     *
     *  @param TFileService   handle to the TFile service
     *  @param string         subdirectory to store the hists in
     */
    void initializeHists(art::ServiceHandle<art::TFileService>&, const std::string&) override;
    
    /**
     *  @brief Interface for method to executve at the end of run processing
     *
     *  @param int            number of events to use for normalization
     */
    void endJob(int numEvents) override;
    
    /**
     *  @brief Interface for filling histograms
     */
    void fillHistograms(const art::Event&)  const override;
    void fillHistograms (const art::Ptr<recob::Track>&) const;
private:
    
    // Fcl parameters.
    std::vector<art::InputTag> fTrackProducerLabelVec; ///< tag for finding the tracks
    std::string                fLocalDirName;          ///< Directory name for histograms
    
    // Pointers to the histograms we'll create.
    TH1D*     fHitsByWire[3];
TH1D*     fHeight[3];
TH1D*     fSumadc[3];

TH1D*     fHeight3mm[3];
TH1D*     fPitch[3];
    TH2D*     fPulseHVsDx[3];
TH2D*     fSumadcVsDx[3];
    TProfile* fPulseHVsHitNo[3];
    
    // Useful services, keep copies for now (we can update during begin run periods)
    const geo::GeometryCore*           fGeometry;             ///< pointer to Geometry service
};
    
//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// pset - Fcl parameters.
///
BasicTrackAnalysis::BasicTrackAnalysis(fhicl::ParameterSet const & pset)
{
    fGeometry           = lar::providerFrom<geo::Geometry>();
    
    configure(pset);
    
    // Report.
    mf::LogInfo("BasicTrackAnalysis") << "BasicTrackAnalysis configured\n";
}

//----------------------------------------------------------------------------
/// Destructor.
BasicTrackAnalysis::~BasicTrackAnalysis()
{}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void BasicTrackAnalysis::configure(fhicl::ParameterSet const & pset)
{
    fTrackProducerLabelVec = pset.get<std::vector<art::InputTag>>("TrackProducerLabel",   std::vector<art::InputTag>());
//fTrackProducerLabelVec = pset.get<std::vector<art::InputTag>>("TrackProducerLabel",   std::vector<art::InputTag>()={"pandoraTrackICARUS"});
    fLocalDirName          = pset.get<std::string               >(  "LocalDirName",       std::string("wow"));
}

//----------------------------------------------------------------------------
/// Begin job method.
void BasicTrackAnalysis::initializeHists(art::ServiceHandle<art::TFileService>& tfs, const std::string& dirName)
{
    // Make a directory for these histograms
    art::TFileDirectory dir = tfs->mkdir(dirName.c_str());

    fHitsByWire[0]    = dir.make<TH1D>("HitsByWire0", ";Wire #", fGeometry->Nwires(0), 0., fGeometry->Nwires(0));
    fHitsByWire[1]    = dir.make<TH1D>("HitsByWire1", ";Wire #", fGeometry->Nwires(1), 0., fGeometry->Nwires(1));
    fHitsByWire[2]    = dir.make<TH1D>("HitsByWire2", ";Wire #", fGeometry->Nwires(2), 0., fGeometry->Nwires(2));
    
    fPulseHVsDx[0] = dir.make<TH2D>("PHVsDx0",    ";Dx;PH", 50, 0.3,0.8,150,0.,150.);
    fPulseHVsDx[1] = dir.make<TH2D>("PHVsDx1",    ";Dx;PH", 50, 0.3,0.8,150,0.,150.);
    fPulseHVsDx[2] = dir.make<TH2D>("PHVsDx2",    ";Dx;PH", 50, 0.3,0.8, 150,0.,150.);
fSumadcVsDx[0] = dir.make<TH2D>("SumadcVsDx0",    ";Dx;Area", 50, 0.3,0.8,150,0.,1500.);
    fSumadcVsDx[1] = dir.make<TH2D>("SumadcVsDx1",    ";Dx;Area", 50, 0.3,0.8,150,0.,1500.);
    fSumadcVsDx[2] = dir.make<TH2D>("SumadcVsDx2",    ";Dx;Area", 50, 0.3,0.8, 150,0.,1500.);
    
    fPulseHVsHitNo[0] = dir.make<TProfile>("PHVsNo0",   ";Hit #;PH", 50, 0., 5., 0., 100.);
    fPulseHVsHitNo[1] = dir.make<TProfile>("PHVsNo1",   ";Hit #;PH", 50, 0., 5., 0., 100.);
    fPulseHVsHitNo[2] = dir.make<TProfile>("PHVsNo2",   ";Hit #;PH", 50, 0., 5., 0., 100.);

    fHeight[0]    = dir.make<TH1D>("Height0", ";PH", 150,0.,150.);
    fHeight[1]    = dir.make<TH1D>("Height1", ";PH", 150,0.,150.);
    fHeight[2]    = dir.make<TH1D>("Height2", ";PH", 150,0.,150.);
 fHeight3mm[0]    = dir.make<TH1D>("Height3mm0", ";PH", 150,0.,150.);
    fHeight3mm[1]    = dir.make<TH1D>("Height3mm1", ";PH", 150,0.,150.);
    fHeight3mm[2]    = dir.make<TH1D>("Height3mm2", ";PH", 150,0.,150.);
fPitch[0]    = dir.make<TH1D>("pitch0", ";Dx", 50,0.,5.);
    fPitch[1]    = dir.make<TH1D>("pitch1", ";Dx", 50,0.,5.);
    fPitch[2]    = dir.make<TH1D>("pitch2", ";Dx", 50,0.,5.);
fSumadc[0]    = dir.make<TH1D>("Sumadc0", ";Area", 150,0.,1500.);
    fSumadc[1]    = dir.make<TH1D>("Sumadc1", ";Area", 150,0.,1500.);
    fSumadc[2]    = dir.make<TH1D>("Sumadc2", ";Area", 150,0.,1500.);
    return;
}
    
void BasicTrackAnalysis::fillHistograms(const art::Event& event) const
{
std::cout << " test filling " << std::endl;
art::ServiceHandle<geo::Geometry const> geom;
//fHitsByWire[0]->Fill(0);
    // The game plan for this module is to look at recob::Tracks and objects associated to tracks
    // To do this we need a valid track collection for those we are hoping to look at
    for(const auto& trackLabel : fTrackProducerLabelVec)
    {
std::cout << " filling track histos " << trackLabel <<std::endl;
        art::Handle<std::vector<recob::Track> > trackHandle;
        event.getByLabel(trackLabel, trackHandle);
        
        if (trackHandle.isValid())
        {
            // Recover the collection of associations between tracks and hits
            art::FindManyP<recob::Hit> trackHitAssns(trackHandle, event, trackLabel);

            for(size_t trackIdx = 0; trackIdx < trackHandle->size(); trackIdx++)
            {
           std::vector<art::Ptr<recob::Hit>> trkHits = trackHitAssns.at(trackIdx);

                art::Ptr<recob::Track> track(trackHandle,trackIdx);
                std::cout << " track " << trackIdx << " hits " <<trkHits.size() << std::endl;
geo::Vector_t track_dir = track->DirectionAtPoint(0);
std::cout << " track dir " << track_dir.X() << " " << track_dir.Y() << " " << track_dir.Z()<< std::endl;
float norm=sqrt(track_dir.X()*track_dir.X()+track_dir.Y()*track_dir.Y()+track_dir.Z()*track_dir.Z());
std::cout << " norm " << norm << std::endl;
                for(long unsigned int jh=0;jh<trkHits.size();jh++) {
auto hitPtr=trkHits[jh];

 const geo::WireID& wireID   = hitPtr->WireID();
        /*float              chi2DOF  = std::min(hitPtr->GoodnessOfFit(),float(249.8));
        int                numDOF   = hitPtr->DegreesOfFreedom();*/
        int                hitMult  = hitPtr->Multiplicity();
        /*float              peakTime = std::min(float(3199.),hitPtr->PeakTime());
        float              charge   = hitPtr->Integral();*/
        float              sumADC   = hitPtr->SummedADC();
        float              hitPH    = hitPtr->PeakAmplitude();
        //float              hitSigma = hitPtr->RMS();
        
        size_t             plane    = wireID.Plane;
       // size_t             wire     = wireID.Wire;
//std::cout << " amplitude " << hitPH << std::endl;
      
 double angleToVert = geom->WireAngleToVertical(hitPtr->View(), hitPtr->WireID().TPC, hitPtr->WireID().Cryostat) - 0.5*::util::pi<>();
double cosgamma = std::abs(std::sin(angleToVert)*track_dir.Y() + std::cos(angleToVert)*track_dir.Z());
std::cout << " plane " << plane <<  " view " << hitPtr->View() << " atv " << angleToVert << std::endl;
 
  double pitch;
  if (cosgamma) {
    pitch = geom->WirePitch(hitPtr->View())/cosgamma;
  }
  else {
    pitch = 0.;
  }

std::cout <<  " plane " << plane << " cosgamma "<< cosgamma << " pitch " << pitch << std::endl;
if(hitMult==1) {
fHeight3mm[plane]->Fill(hitPH*(0.3/pitch));
fHeight[plane]->Fill(hitPH);
fPitch[plane]->Fill(pitch);
fSumadc[plane]->Fill(sumADC);
fPulseHVsDx[plane]->Fill(pitch,hitPH);
fSumadcVsDx[plane]->Fill(pitch,sumADC);
}
}
            }
        }
    }

    return;
}

void BasicTrackAnalysis::fillHistograms( const art::Ptr<recob::Track>& trk) const
{
 

    return;
}

// Useful for normalizing histograms
void BasicTrackAnalysis::endJob(int numEvents)
{
    // Normalize wire profiles to be hits/event
    double normFactor(1./numEvents);
    
    for(size_t idx = 0; idx < 3; idx++) fHitsByWire[idx]->Scale(normFactor);
    
    return;
}
    
DEFINE_ART_CLASS_TOOL(BasicTrackAnalysis)
}
