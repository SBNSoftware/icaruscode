//std includes
#include <vector>
#include <map>
#include <fstream>
#include <string>
#include <cmath>

//ROOT includes
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TSpline.h"

//Framework includes
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "art/Framework/Core/ResultsProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Results.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Utilities/InputTag.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"

namespace WireMod {
  class WireModMakeHists : public art::ResultsProducer {
    public:
      explicit WireModMakeHists(fhicl::ParameterSet const& pset);
      ~WireModMakeHists() override = default;
      void event(art::Event const& evt) override;
      void reconfigure(fhicl::ParameterSet const& p);
      void writeResults(art::Results& r) override;
      void clear() override;

      void endJob(); // added so histograms can be made and saved after all the events have been processed

    private:                         
      const geo::GeometryCore* fGeometry = lar::providerFrom<geo::Geometry>();
      art::InputTag fLabel;      // how the hits/wires are labeled in the input file
      art::InputTag fTrackLabel; // how the tracks are labeled in the input file
      bool fGetHits;             // are we getting hits? if false the label is for the wires
      bool fGetTracks;           // are we getting the tracks?
      bool fXCorrection;         // are we applying the X correction?
      TFile* fXCorrFile;         // the file which contains the X corrections
 };// end WireModMakeHists class
        //
       //-------------------------------------------------------------
      // Define the constructor
     // The way ART works is it will call the construct on the fhicl parameter set
    // Basically we want it to pull out the variables we defined there
   // We're doing this in the reconfigure function (because that makes it clearer what we are doing imo)
  // so just call that function
  WireModMakeHists::WireModMakeHists(fhicl::ParameterSet const& pset)
  {
    this->reconfigure(pset);
  }
    //
   //-------------------------------------------------------------
  // this function pulls out the stuff we defined in the fhicl file
  void WireModMakeHists::reconfigure(fhicl::ParameterSet const& pset)
  {
    // the first arguement is where in the fhicl to look, the second is the default value if that isn't found
    fLabel = pset.get<art::InputTag>("Label", "decon1droi");
    fTrackLabel = pset.get<art::InputTag>("TrackLabel");
    fGetHits = pset.get<bool>("GetHits", false);
    fGetTracks = pset.get<bool>("GetTracks", false);
    fXCorrection = pset.get<bool>("XCorrection", false);
    if (fXCorrection)
    {
      fXCorrFile = new TFile("XCorrFile.root", "READ");
    }
    //fPFParticleLabel = pset.get<art::InputTag>("PFParticleLabel");
  }
  //
 // Used to see if path goes through anode or cathode
 bool hasValuesAroundPoints(const std::vector<double>& vec, const std::vector<double>& points) {
  for (double point : points) {
    bool hasAbove = false;
    bool hasBelow = false;
    for (double value : vec) {
      if (value > point) {
        hasAbove = true;
      } 
      if (value < point) {
        hasBelow = true;
      }
      if (hasAbove && hasBelow) {
        return true;// Found a point with values both above and below it in the vector
      }            //
    }             //
  }              //
  return false; // No point was found with values both above and below it in the vector
 }

  //Here is a 5D map to store all of the data to be made in Histograms----------------------------------------------
 /* The Following are what the indicies correspond to:
 a = 0,1,2 are the planes
 b = 0,1 are E and W
 c = 0,1,2,3 are E, E, W, W respectively
 d = 0,1,2,3,4,5,6, 7, 8 correspond to values of x,y,z,theta,phi,maximum,fwhm, gaussian height, and gaussian chi^2 respectively
 e = 0,1: 1 indicates if the track crosses any anode or cathode
 */
 class FiveDMap {
    class FourDMap {
        class ThreeDMap {
            class TwoDMap {
                class OneDMap {
                public:
                    std::vector<double>& operator[](int index) {
                        return data[index];
                    }
                private:
                    std::map<int, std::vector<double>> data;
                };

            public:
                OneDMap& operator[](int index) {
                    return data[index];
                }
            private:
                std::map<int, OneDMap> data;
            };

        public:
            TwoDMap& operator[](int index) {
                return data[index];
            }
        private:
            std::map<int, TwoDMap> data;
        };

    public:
        ThreeDMap& operator[](int index) {
            return data[index];
        }
    private:
        std::map<int, ThreeDMap> data;
    };

  public:
    FourDMap& operator[](int index) {
        return data[index];
    }
  private:
    std::map<int, FourDMap> data;
 };
 //Declaring the class the holds all the the information-------------------
  //
 //Here are some vectors that will be used to generate the names of the histograms when later looping thorugh data
 const std::vector<std::string> planes = {"0","1","2"};
 const std::vector<std::string> tpc0s = {"E","W"};
 const std::vector<std::string> tpc1s = {"E","E","W","W"};
 const std::vector<std::string> variables = {"X","Y","Z","Theta", "Phi", "Amplitude", "FWHM","GaussWidth", "GaussHeight", "GaussFit"};
 const std::vector<std::string> hitchoices = {"ALL_HITS","SELECTED_HITS"};
  //
 // Bin information--------------------------------------------------
 const int nbins_amp = 100;
 const std::vector<double> mins = {-360.0, -182.0, -896, 0, 0, 0, 0, 0, 0, 0};
 const std::vector<double> maxs = {360.0, 182.0, 896, acos(0), acos(0), 200, 200, 20, 200, 200};


 //function for making and saving histograms
 void CreateAndSaveHist(const char* name, const char* title, int nbinsX, double minX, double maxX, std::vector<double>& x_values, std::vector<double>& y_values) {
    if (x_values.size() != y_values.size() || x_values.empty()) {
        // Handle error: sizes not equal or vectors are empty
        return;
    }
    //TFile *file = TFile::CurrentFile();
    //if (!file || !file->IsOpen() || file->IsZombie()) {
    //    // Handle error: no file is open or file is not in a good state
    //    return;
    //}
    //TH2D *hist = nullptr;
    //try {
    //    hist = new TH2D(name, title, nbinsX, minX, maxX, nbins_amp, 0, 200);
    //} catch (const std::bad_alloc&) {
    //    // Handle error: memory allocation failed
    //    return;
    //}
    //for (size_t i = 0; i < x_values.size(); i++) {
    //    hist->Fill(x_values[i], y_values[i]);
    //}
    //hist->Write();
    ///*TProfile *profX = hist->ProfileX();
    //if (profX) {
    //    profX->Write();
    //}*/
    //delete hist;
    
    // let's try this better...
    art::ServiceHandle<art::TFileService>  tfs;
    TH2D* hist = tfs->make<TH2D>(name, title, nbinsX, minX, maxX, nbins_amp, 0, 200);
    for (size_t i = 0; i < x_values.size(); i++) {
      hist->Fill(x_values[i], y_values[i]);
    }
    mf::LogVerbatim("WireModMakeHists")
      << "Made 2D histogram " << hist->GetName();
 } 

 //function for making and saving histograms that are 3D
 void CreateAndSaveHist3D(const char* name, const char* title, int nbinsX, double minX, double maxX, std::vector<double>& x_values, int nbinsY, double minY, double maxY, std::vector<double>& y_values, int nbinsZ, double minZ, double maxZ, std::vector<double>& z_values) {
    if (x_values.size() != y_values.size() || x_values.size() !=z_values.size() || x_values.empty() || z_values.empty()) {
        // Handle error: sizes not equal or vectors are empty
        return;
    }
    //TFile *file = TFile::CurrentFile();
    //if (!file || !file->IsOpen() || file->IsZombie()) {
    //    // Handle error: no file is open or file is not in a good state
    //    return;
    //}
    //TH3D *hist = nullptr;
    //try {
    //    hist = new TH3D(name, title, nbinsX, minX, maxX, nbinsY, minY, maxY, nbinsZ, minZ, maxZ);
    //} catch (const std::bad_alloc&) {
    //    // Handle error: memory allocation failed
    //    return;
    //}
    //for (size_t i = 0; i < x_values.size(); i++) {
    //    hist->Fill(x_values[i], y_values[i], z_values[i]);
    //}
    //hist->Write();
    ///*TProfile *profX = hist->ProfileX();
    //if (profX) {
    //    profX->Write();
    //}*/
    //delete hist;
    
    // let's try this better...
    art::ServiceHandle<art::TFileService>  tfs;
    TH3D* hist = tfs->make<TH3D>(name, title, nbinsX, minX, maxX, nbinsY, minY, maxY, nbinsZ, minZ, maxZ);
    for (size_t i = 0; i < x_values.size(); i++) {
      hist->Fill(x_values[i], y_values[i], z_values[i]);
    }
    mf::LogVerbatim("WireModMakeHists")
      << "Made 3D histogram " << hist->GetName();
 } 

 FiveDMap data;


      //--------------------------------------------------------------
     //-----------------------------------------------------------------------------------------------------------------
    // this function is run on every event in the art file
   // the event stores the information we want to analyze
  //----------------------------------------------------------------------------------------------------------------------------------------------------------------
 void WireModMakeHists::event(art::Event const& evt)
 {  
    // this is what will make our histograms for us
    art::ServiceHandle<art::TFileService>  tfs;
     //
    // get a unique string for this event
    std::string evtStr = std::to_string(evt.id().run()) + "_"
                       + std::to_string(evt.id().subRun()) + "_"
                       + std::to_string(evt.id().event()) + "_";


    if (fGetTracks)
    {
      // we need both the hits and the tracks separeately
      art::Handle<std::vector<recob::Hit>> hitHandle;
      evt.getByLabel(fLabel, hitHandle);
      if (not hitHandle.isValid())
      {
        //MF_LOG_VERBATIM("WireModWireModMakeHists");
          //<< "Hit handle is not valid!" << '\n'
          //<< "Tried " << fLabel << '\n'
          //<< "abort";
        return;
      }
      std::vector<art::Ptr<recob::Hit>> hitPtrVec;
      art::fill_ptr_vector(hitPtrVec, hitHandle);
       //
      // also get tracks and the hits for each track
      art::FindOneP<recob::Wire> hitToWireAssns(hitHandle, evt, fLabel);
      art::Handle<std::vector<recob::Track>> trackHandle;
      evt.getByLabel(fTrackLabel, trackHandle);
      if (not trackHandle.isValid())
      {
        //MF_LOG_VERBATIM("WireModWireModMakeHists");
          //<< "Track handle is not valid!" << '\n'
          //<< "Tried " << fTrackLabel << '\n'
          //<< "abort";
        return;
      }
      std::vector<art::Ptr<recob::Track>> trackPtrVec;
      art::fill_ptr_vector(trackPtrVec, trackHandle);
       //
      // the recob::TrackHitMeta will let us find where in the track each hit is
      art::FindManyP<recob::Hit, recob::TrackHitMeta> trackToHits(trackHandle, evt, fTrackLabel);

      //Get PFParticles
      //auto pfpListHandle = evt.getValidHandle<std::vector<recob::PFParticle> >(fPFParticleLabel[it]);
      //if (!pfpListHandle.isValid()) continue;

      for (auto const& trackPtr : trackPtrVec)
      {
        //creating a data set for this track to check for selected hits
        FiveDMap loopdata;
        // get the tracks hits and metadata for the hits
        std::vector<art::Ptr<recob::Hit>> const& trackHits = trackToHits.at(trackPtr.key());
        std::vector<const recob::TrackHitMeta*> const& trackHitMetas = trackToHits.data(trackPtr.key());
         // loop over the track hits
        // will need to also find the same hit in hitPtrVec to get the associated wire
        for (size_t hitIdx = 0; hitIdx < trackHits.size(); ++hitIdx)
        {
          // get the track hit and metadata
          art::Ptr<recob::Hit> trackHit = trackHits[hitIdx];
          const recob::TrackHitMeta& hitMeta = *trackHitMetas[hitIdx];
          if (hitMeta.Index() == std::numeric_limits<unsigned int>::max() || not trackPtr->HasValidPoint(hitMeta.Index()))
          {
            //MF_LOG_VERBATIM("WireModWireModMakeHists");
              //<< "Bad Hit, get another one";
              continue;
          }
          recob::Track::Point_t  const& hitLoc = trackPtr->LocationAtPoint(hitMeta.Index());
          recob::Track::Vector_t const& hitDir = trackPtr->DirectionAtPoint(hitMeta.Index());

          // set up a wirePtr, loop over hitPtrVec to find the right wire
          art::Ptr<recob::Wire> wirePtr;
          for (auto const& hitPtr : hitPtrVec)
          {  //
            // check channel etc to see if the hits are a match
            if (trackHit->Channel()   == hitPtr->Channel()   &&
                trackHit->StartTick() == hitPtr->StartTick() &&
                trackHit->EndTick()   == hitPtr->EndTick()   && 
                trackHit->PeakTime()  == hitPtr->PeakTime()  )
            {
              //MF_LOG_VERBATIM("WireModWireModMakeHists");
                //<< "Wire Found!";
              wirePtr = hitToWireAssns.at(hitPtr.key());
               // 
              //getting the wire plane and TPC info
              int planeID = hitPtr->WireID().Plane;
              int tpcID = (hitPtr->WireID().TPC) / 2;
              int cryoID = hitPtr->WireID().Cryostat;

              //Making a histogram of the hit information --------------------------------------------------------------------
              size_t fstTick = hitPtr->StartTick();
              size_t endTick = hitPtr->EndTick();
              size_t hitWidth = endTick - fstTick;
                // the start and end ticks were aquired assuming the hits are gaussian
               // we want to plot a bit of buffer around the region of interest to get the full shape
              // default to a hit width on either side, but if the wire doesn't have enough ticks then use what we can
              size_t nTicks = wirePtr->NSignal();
              size_t fstBuffer = (fstTick > hitWidth)          ? hitWidth : fstTick;
              size_t endBuffer = (nTicks > endTick + hitWidth) ? hitWidth : nTicks - endTick;
                 // put the waveform in the histogram
                // tfs will make a whatever is in the <>, (in this case a TH1F)
               // the agruements past to it should be the same as for the constructor for the <object>
              auto waveformHist = std::make_unique<TH1F>(("adc_"+evtStr+fLabel.label()+"_"+std::to_string(wirePtr.key())).c_str(), //> name of the object
                                             ";Sample;Arbitrary Units",                                                //> axes labels
                                             hitWidth + fstBuffer + endBuffer,                                         //> number of bins
                                             fstTick - fstBuffer,                                                      //> lowest edge
                                             endTick + endBuffer);                                                     //> upper edge
              // ROOT counts from 1, everyone else counts from 0
             //for (size_t bin = 1; bin < endTick - fstTick + 1; ++bin)
             for (size_t bin = 1; bin < hitWidth + fstBuffer + endBuffer + 1; ++bin)
             {
                float wvfmVal = (wirePtr->Signal())[fstTick - fstBuffer + bin - 1]; 
                waveformHist->SetBinContent(bin, wvfmVal);
             }
              //
             // Calculate amplitude and FWHM and adding to vector of data---------------
             double amplitude = waveformHist-> GetMaximum();
             double half_max = amplitude / 2.0;
             int bin1 = waveformHist->FindFirstBinAbove(half_max);
             int bin2 = waveformHist->FindLastBinAbove(half_max);
             double fwhm = waveformHist->GetBinCenter(bin2) - waveformHist->GetBinCenter(bin1);
             
             data[planeID][cryoID][tpcID][5][0].push_back(amplitude);
             loopdata[planeID][cryoID][tpcID][5][0].push_back(amplitude);
             data[planeID][cryoID][tpcID][6][0].push_back(fwhm);
             loopdata[planeID][cryoID][tpcID][6][0].push_back(fwhm);

             //getting the gaussian information
             double gaussheight = hitPtr->PeakAmplitude();
             double gausswidth = hitPtr->RMS();
             double gaussfit = hitPtr->GoodnessOfFit();
             if (fXCorrection)
             {
               TSpline3* heightSpline = fXCorrFile->Get<TSpline3>(("xHeight_" + std::to_string(planeID) + "_ratio_prof").c_str());
               TSpline3* widthSpline  = fXCorrFile->Get<TSpline3>(("xWidth_"  + std::to_string(planeID) + "_ratio_prof").c_str());
               gaussheight *= heightSpline->Eval(hitLoc.X());
               gausswidth  *= widthSpline ->Eval(hitLoc.X());
             }

             data[planeID][cryoID][tpcID][7][0].push_back(gausswidth);
             data[planeID][cryoID][tpcID][8][0].push_back(gaussheight);
             data[planeID][cryoID][tpcID][9][0].push_back(gaussfit);
             loopdata[planeID][cryoID][tpcID][7][0].push_back(gausswidth);
             loopdata[planeID][cryoID][tpcID][8][0].push_back(gaussheight);
             loopdata[planeID][cryoID][tpcID][9][0].push_back(gaussfit);

             //Adding values to vectors for plotting
             data[planeID][cryoID][tpcID][0][0].push_back(hitLoc.X());
             data[planeID][cryoID][tpcID][1][0].push_back(hitLoc.Y());
             data[planeID][cryoID][tpcID][2][0].push_back(hitLoc.Z());
             loopdata[planeID][cryoID][tpcID][0][0].push_back(hitLoc.X());
             loopdata[planeID][cryoID][tpcID][1][0].push_back(hitLoc.Y());
             loopdata[planeID][cryoID][tpcID][2][0].push_back(hitLoc.Z());

             double sinZ = fGeometry->WirePtr(hitPtr->WireID())->SinThetaZ();
             double cosZ = fGeometry->WirePtr(hitPtr->WireID())->CosThetaZ();
             double localXDir =      hitDir.X();
             double localYDir = cosZ*hitDir.Y() + sinZ*hitDir.Z();
             double localZDir = cosZ*hitDir.Z() - sinZ*hitDir.Y();
             double theta = abs(atan(localXDir / localZDir));
             double phi   = abs(atan(localYDir / localZDir));

             //if (planeID == 0){
             // theta = atan(abs(hitLoc.X()) / abs(0.867* hitLoc.Y() - 0.5*hitLoc.Z()));
             // phi = atan( abs((0.5*hitLoc.Y() + 0.867* hitLoc.Z())) / abs(0.867* hitLoc.Y() - 0.5*hitLoc.Z()));
             //}
             //else if (planeID == 1){
             // theta = atan(hitLoc.X() / (0.867* hitLoc.Y() + 0.5*hitLoc.Z()));
             // phi = atan( (0.5*hitLoc.Y() - 0.867* hitLoc.Z()) / (0.867* hitLoc.Y() + 0.5*hitLoc.Z()));
             //}
             //if (planeID == 2){
             // theta = atan(hitLoc.X() / (hitLoc.Z()));
             // phi = atan( (hitLoc.Y()) / (hitLoc.Z()));
             //}

             data[planeID][cryoID][tpcID][3][0].push_back(theta);
             loopdata[planeID][cryoID][tpcID][3][0].push_back(theta);
             data[planeID][cryoID][tpcID][4][0].push_back(phi);
             loopdata[planeID][cryoID][tpcID][4][0].push_back(phi);

             //delete waveformHist;
        
            }         
          }
          if (wirePtr.isNull())
          {
            //MF_LOG_VERBATIM("WireModWireModMakeHists");
              //<< "Couldn't find wire" << '\n'
             // << "Continue...";
            continue;
           }
            // get the location from the track using the hitMeta
           // get X using hitLoc.X(), similarly for Y and Z
          // more info from what you can get at https://sbnsoftware.github.io/doxygen/d9/dca/classrecob_1_1Track.html
          //MF_LOG_VERBATIM("WireModWireModMakeHists");
            //<< "Hit Pos is (" << hitLoc.X() << ", " << hitLoc.Y() << ", " << hitLoc.Z() << ")"; 
        }// 
       //finding the tracks that have selected hits
      std::vector<double> points = {-358.49,-210.15,-61.94,61.94,210.15,358.49};
      for (int a : {0, 1, 2}) {
        for (int b : {0, 1}) {
          for (int c : {0, 1, 2, 3}){ 
            if (hasValuesAroundPoints(loopdata[a][b][c][0][0],points)){
              data[a][b][c][0][1].insert(data[a][b][c][0][1].end(), loopdata[a][b][c][0][0].begin(), loopdata[a][b][c][0][0].end());
              data[a][b][c][1][1].insert(data[a][b][c][1][1].end(), loopdata[a][b][c][1][0].begin(), loopdata[a][b][c][1][0].end());
              data[a][b][c][2][1].insert(data[a][b][c][2][1].end(), loopdata[a][b][c][2][0].begin(), loopdata[a][b][c][2][0].end());
              data[a][b][c][3][1].insert(data[a][b][c][3][1].end(), loopdata[a][b][c][3][0].begin(), loopdata[a][b][c][3][0].end());
              data[a][b][c][4][1].insert(data[a][b][c][4][1].end(), loopdata[a][b][c][4][0].begin(), loopdata[a][b][c][4][0].end());
              data[a][b][c][5][1].insert(data[a][b][c][5][1].end(), loopdata[a][b][c][5][0].begin(), loopdata[a][b][c][5][0].end());
              data[a][b][c][6][1].insert(data[a][b][c][6][1].end(), loopdata[a][b][c][6][0].begin(), loopdata[a][b][c][6][0].end());
              data[a][b][c][7][1].insert(data[a][b][c][7][1].end(), loopdata[a][b][c][7][0].begin(), loopdata[a][b][c][7][0].end());
              data[a][b][c][8][1].insert(data[a][b][c][8][1].end(), loopdata[a][b][c][8][0].begin(), loopdata[a][b][c][8][0].end());
              data[a][b][c][9][1].insert(data[a][b][c][9][1].end(), loopdata[a][b][c][9][0].begin(), loopdata[a][b][c][9][0].end());
            }
          }
        }
      }
     } 
    }
  }

  //function for generating name of file to be saved to
  std::string generateUniqueFilename(const std::string& baseName, const std::string& extension) {
    std::string filename = baseName + extension;
    std::ifstream file(filename);
    int counter = 1;

    // Check if the file exists
    while (file.good()) {
        file.close(); // Close the previously opened file

        // Generate a new filename with an incremented counter
        filename = baseName + "_" + std::to_string(counter++) + extension;

        // Attempt to open the new filename
        file.open(filename);
    }

    // Filename is unique at this point
    return filename;
  }

  void WireModMakeHists::endJob(){
    std::string labelStr = fLabel.label();
    std::string instanceStr = fLabel.instance();
    auto instance = labelStr + instanceStr;
    //std::cout << instance << std::endl;
    // all of the instances produce the same histograms, so only make histograms on one of them. I picked this one  
    if (instance == "gaushitTPCWW"){
    //if (true){   
      //making a Tfile for the histograms
      auto name = generateUniqueFilename("Data_Histograms",".root");
      TFile *file = new TFile(name.c_str(),"RECREATE");
      //
      //looping over the data mapping to make histograms of all the data (if confused see top for what variables are)
      for (int a : {0, 1, 2}) {
        for (int b : {0, 1}) {
          for (int c : {0, 1}){ 
            for (int d : {0, 1, 2, 3, 4}) {
              for (int e : {0,1}){
              mf::LogVerbatim("WireModMakeHists")
                << "a: " << a << ", b: " << b << ", c: " << c << ", d: " << d << ", e: " << e;
              //amp_histograms
              CreateAndSaveHist(("Plane_" + planes[a] + "_" + tpc0s[b] + "_" + tpc0s[c] + "_" + variables[d] + "_" + variables[5] + "_" + hitchoices[e]).c_str(),
                                  (variables[d] + "_" + variables[5] +  ";" + variables[d] + ";" + variables[5]).c_str(),
                                  nbins_amp,
                                  mins[d],
                                  maxs[d],
                                  data[a][b][c][d][e],
                                  data[a][b][c][5][e]
                                );
              //FWHM Histograms
              CreateAndSaveHist(("Plane_" + planes[a] + "_" + tpc0s[b] + "_" +  tpc0s[c] + "_" + variables[d] + "_" + variables[6] + "_" + hitchoices[e]).c_str(),
                                  (variables[d] + "_" + variables[6] +  ";" + variables[d] + ";" + variables[6]).c_str(),
                                  nbins_amp,
                                  mins[d],
                                  maxs[d],
                                  data[a][b][c][d][e],
                                  data[a][b][c][6][e]
                                );
              //Gaussian Hitograms               
              CreateAndSaveHist(("Plane_" + planes[a] + "_" + tpc0s[b] + "_" +  tpc0s[c] + "_" + variables[d] + "_" + variables[7] + "_" + hitchoices[e]).c_str(),
                                  (variables[d] + "_" + variables[7] +  ";" + variables[d] + ";" + variables[7]).c_str(),
                                  nbins_amp,
                                  mins[d],
                                  maxs[d],
                                  data[a][b][c][d][e],
                                  data[a][b][c][7][e]
                                );
              CreateAndSaveHist(("Plane_" + planes[a] + "_" + tpc0s[b] + "_" +  tpc0s[c] + "_" + variables[d] + "_" + variables[8] + "_" + hitchoices[e]).c_str(),
                                  (variables[d] + "_" + variables[8] +  ";" + variables[d] + ";" + variables[8]).c_str(),
                                  nbins_amp,
                                  mins[d],
                                  maxs[d],
                                  data[a][b][c][d][e],
                                  data[a][b][c][8][e]
                                );
              CreateAndSaveHist(("Plane_" + planes[a] + "_" + tpc0s[b] + "_" +  tpc0s[c] + "_" + variables[d] + "_" + variables[9] + "_" + hitchoices[e]).c_str(),
                                  (variables[d] + "_" + variables[9] +  ";" + variables[d] + ";" + variables[9]).c_str(),
                                  nbins_amp,
                                  mins[d],
                                  maxs[d],
                                  data[a][b][c][d][e],
                                  data[a][b][c][9][e]
                                );                                                  
              }
            }
          }
        }
      }
      //now for making the Y Z (variables[1] + variables[2] 3D histograms:
      for (int a : {0, 1, 2}) {
        for (int b : {0, 1}) {
          for (int c : {0, 1}){ 
            for (int e : {0,1}){
              //amp_histograms
              CreateAndSaveHist3D(("Plane_" + planes[a] + "_" + tpc0s[b] + "_" + tpc0s[c] + "_" + variables[1] + "_" + variables[2] + "_" + variables[5] + "_" + hitchoices[e]).c_str(),
                                  (variables[1] + "_" + variables[2] + "_" + variables[5] +  ";" + variables[1] + "_" + variables[2] + ";" + variables[5]).c_str(),
                                  nbins_amp,
                                  mins[1],
                                  maxs[1],
                                  data[a][b][c][1][e],
                                  nbins_amp,
                                  mins[2],
                                  maxs[2],
                                  data[a][b][c][2][e],
                                  nbins_amp,
                                  0.,
                                  200.,
                                  data[a][b][c][5][e]
                                );
              //FWHM Histograms
              CreateAndSaveHist3D(("Plane_" + planes[a] + "_" + tpc0s[b] + "_" +  tpc0s[c] + "_" + variables[1] + "_" + variables[2] + "_" + variables[6] + "_" + hitchoices[e]).c_str(),
                                  (variables[1] + "_" + variables[2] + "_" + variables[6] +  ";" + variables[1] + "_" + variables[2] + ";" + variables[6]).c_str(),
                                  nbins_amp,
                                  mins[1],
                                  maxs[1],
                                  data[a][b][c][1][e],
                                  nbins_amp,
                                  mins[2],
                                  maxs[2],
                                  data[a][b][c][2][e],
                                  nbins_amp,
                                  0.,
                                  200.,
                                  data[a][b][c][6][e]
                                );

              //gauss amp histograms
              CreateAndSaveHist3D(("Plane_" + planes[a] + "_" + tpc0s[b] + "_" + tpc0s[c] + "_" + variables[1] + "_" + variables[2] + "_" + variables[7] + "_" + hitchoices[e]).c_str(),
                                  (variables[1] + "_" + variables[2] + "_" + variables[7] +  ";" + variables[1] + "_" + variables[2] + ";" + variables[7]).c_str(),
                                  nbins_amp,
                                  mins[1],
                                  maxs[1],
                                  data[a][b][c][1][e],
                                  nbins_amp,
                                  mins[2],
                                  maxs[2],
                                  data[a][b][c][2][e],
                                  nbins_amp,
                                  0.,
                                  200.,
                                  data[a][b][c][7][e]
                                );
              //gauss width Histograms
              CreateAndSaveHist3D(("Plane_" + planes[a] + "_" + tpc0s[b] + "_" +  tpc0s[c] + "_" + variables[1] + "_" + variables[2] + "_" + variables[8] + "_" + hitchoices[e]).c_str(),
                                  (variables[1] + "_" + variables[2] + "_" + variables[8] +  ";" + variables[1] + "_" + variables[2] + ";" + variables[8]).c_str(),
                                  nbins_amp,
                                  mins[1],
                                  maxs[1],
                                  data[a][b][c][1][e],
                                  nbins_amp,
                                  mins[2],
                                  maxs[2],
                                  data[a][b][c][2][e],
                                  nbins_amp,
                                  0.,
                                  200.,
                                  data[a][b][c][8][e]
                                );
             }
           }
        }
      }

        
      file->Close();
      delete file;
    }
 }
  //-------------------------------------------------------------
  // writeResults is currently unused, but we still need to define it
  // just have it do nothing for now
  void WireModMakeHists::writeResults(art::Results& r)
  {
  } 
  //-------------------------------------------------------------
  // clear is currently unused, but we still need to define it
  // just have it do nothing for now
  void WireModMakeHists::clear()
  {
  }
  DEFINE_ART_RESULTS_PLUGIN(WireModMakeHists)
}// end WireMod namespace
