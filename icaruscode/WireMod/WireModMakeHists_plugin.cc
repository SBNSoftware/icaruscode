// std includes
#include <vector>
#include <map>
#include <fstream>
#include <string>
#include <cmath>

// ROOT includes
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TSpline.h"

// Framework includes
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

  // 5D map to store all of the data to be made in Histograms
  /* The Following are what the indices correspond to:
  a = 0, 1, 2 are the planes
  b = 0, 1 are E and W
  c = 0, 1, 2, 3 are E, E, W, W respectively
  d = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 correspond to values of x, y, z, theta, phi, inetgral, fwhm, gaussian width, gaussian integral, and gaussian chi^2 respectively
  e = 0, 1: 1 indicates if the track crosses any anode or cathode
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

  /**
   * @brief Class for making histograms from wire modification data.
   */
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
      art::InputTag fLabel;      ///< Label for hits/wires in the input file
      art::InputTag fTrackLabel; ///< Label for tracks in the input file
      bool fGetHits;             ///< Flag for getting hits
      bool fGetTracks;           ///< Flag for getting tracks
      bool fXCorrection;         ///< Flag for applying X correction
      TFile* fXCorrFile;         ///< File containing the X corrections
      FiveDMap data; ///< 5D map to store data for histograms
  };


  /**
   * @brief Constructor for WireModMakeHists class.
   * @param pset The parameter set from which to extract the configuration variables.
   */
  WireModMakeHists::WireModMakeHists(fhicl::ParameterSet const& pset) {
      this->reconfigure(pset);
  }


  // This function pulls out the stuff we defined in the fhicl file
  /**
   * @brief Reconfigures the WireModMakeHists object.
   * @param pset The parameter set where to look for the configuration parameters.
   */
  void WireModMakeHists::reconfigure(fhicl::ParameterSet const& pset) {
      fLabel = pset.get<art::InputTag>("Label", "decon1droi");
      fTrackLabel = pset.get<art::InputTag>("TrackLabel");
      fGetHits = pset.get<bool>("GetHits", false);
      fGetTracks = pset.get<bool>("GetTracks", false);
      fXCorrection = pset.get<bool>("XCorrection", false);
      if (fXCorrection) {
          fXCorrFile = new TFile("XCorrFile.root", "READ");
      }
  }


  /**
   * @brief Checks if the path goes through anode or cathode by seeing if there are values around specified points.
   * @param vec The vector of values to be checked against the points.
   * @param points The vector of points to check if there are values around them in the 'vec' vector.
   * @return true if any point in 'points' has values in 'vec' both above and below it, false otherwise.
   */
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
          return true;
        }
      }
    }
    return false; 
  }

  
  // Used to generate unique histogram names 
  const std::vector<std::string> planes = {"0", "1", "2"};
  const std::vector<std::string> tpc0s = {"E", "W"};
  const std::vector<std::string> tpc1s = {"E", "E", "W", "W"};
  const std::vector<std::string> variables = {"X", "Y", "Z", "Theta", "Phi", "Integral", "FWHM", "GaussWidth", "GaussIntegral", "GaussFit"};
  const std::vector<std::string> hitchoices = {"ALL_HITS", "SELECTED_HITS"};

  
  // Bin information
  const int nbins_amp = 100;
  const std::vector<double> mins = {-360.0, -182.0, -896, 0, 0, 0, 0, 0, 0, 0};
  const std::vector<double> maxs = {360.0, 182.0, 896, acos(0), acos(0), 200, 200, 20, 200, 200};


  /**
   * @brief Creates and saves a 2D histogram.
   * @param name The name of the histogram.
   * @param title The title of the histogram.
   * @param nbinsX The number of bins for the X axis.
   * @param minX The minimum value for the X axis.
   * @param maxX The maximum value for the X axis.
   * @param x_values The vector of X values to fill the histogram.
   * @param y_values The vector of Y values to fill the histogram.
   */
  void CreateAndSaveHist(const char* name, const char* title, int nbinsX, double minX, double maxX, std::vector<double>& x_values, std::vector<double>& y_values) {
    if (x_values.size() != y_values.size() || x_values.empty()) {
      return;
    }
    art::ServiceHandle<art::TFileService> tfs;
    TH2D* hist = tfs->make<TH2D>(name, title, nbinsX, minX, maxX, nbins_amp, 0, 200);
    for (size_t i = 0; i < x_values.size(); i++) {
      hist->Fill(x_values[i], y_values[i]);
    }
    mf::LogVerbatim("WireModMakeHists") << "Made 2D histogram " << hist->GetName();
  }


  // Function for making and saving 3D histograms
  /**
   * @brief Creates and saves a 3D histogram.
   * @param name The name of the histogram.
   * @param title The title of the histogram.
   * @param nbinsX The number of bins for the X axis.
   * @param minX The minimum value for the X axis.
   * @param maxX The maximum value for the X axis.
   * @param x_values The vector of X values to fill the histogram.
   * @param nbinsY The number of bins for the Y axis.
   * @param minY The minimum value for the Y axis.
   * @param maxY The maximum value for the Y axis.
   * @param y_values The vector of Y values to fill the histogram.
   * @param nbinsZ The number of bins for the Z axis.
   * @param minZ The minimum value for the Z axis.
   * @param maxZ The maximum value for the Z axis.
   * @param z_values The vector of Z values to fill the histogram.
   */
  void CreateAndSaveHist3D(const char* name, const char* title, int nbinsX, double minX, double maxX, std::vector<double>& x_values, int nbinsY, double minY, double maxY, std::vector<double>& y_values, int nbinsZ, double minZ, double maxZ, std::vector<double>& z_values) {
    if (x_values.size() != y_values.size() || x_values.size() != z_values.size() || x_values.empty() || z_values.empty()) {
        return;
    }
    art::ServiceHandle<art::TFileService> tfs;
    TH3D* hist = tfs->make<TH3D>(name, title, nbinsX, minX, maxX, nbinsY, minY, maxY, nbinsZ, minZ, maxZ);
    for (size_t i = 0; i < x_values.size(); i++) {
      hist->Fill(x_values[i], y_values[i], z_values[i]);
    }
    mf::LogVerbatim("WireModMakeHists") << "Made 3D histogram " << hist->GetName();
  }


  /**
   * @brief Generates a unique filename by appending a counter if the file already exists.
   * @param baseName The base name of the file.
   * @param extension The extension of the file.
   * @return A unique filename.
   */
  std::string generateUniqueFilename(const std::string& baseName, const std::string& extension) {
    std::string filename = baseName + extension;
    std::ifstream file(filename);
    int counter = 1;
    while (file.good()) {
      file.close();
      // Generate a new filename with an incremented counter
      filename = baseName + "_" + std::to_string(counter++) + extension;
      file.open(filename);
    }
    return filename;
  }


  /**
   * @brief Creates and saves histograms after all events have been processed.
   */
  void WireModMakeHists::endJob(){
  std::string labelStr = fLabel.label();
  std::string instanceStr = fLabel.instance();
  auto instance = labelStr + instanceStr;
  //std::cout << instance << std::endl;
  // all of the instances produce the same histograms, so only make histograms on one of them. I picked this one  
  if (instance == "gaushitTPCWW"){   
    //making a Tfile for the histograms
    auto name = generateUniqueFilename("Data_Histograms",".root");
    TFile *file = new TFile(name.c_str(),"RECREATE");
    
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


  /**
   * @brief Writes results to the ART framework. TODO
   * @param r The results object to be written.
   */
  void WireModMakeHists::writeResults(art::Results& r) {}


  /**
   * @brief Clears the data. TODO
   */
  void WireModMakeHists::clear() {}


  /**
   * @brief Processes each event in the ART file.
   * @param evt The event to be processed.
   */
  void WireModMakeHists::event(art::Event const& evt)
  {  
    // this is what will make our histograms for us
    art::ServiceHandle<art::TFileService>  tfs;

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
        
      // the recob::TrackHitMeta will let us find where in the track each hit is
      art::FindManyP<recob::Hit, recob::TrackHitMeta> trackToHits(trackHandle, evt, fTrackLabel);

      //Get PFParticles
      //auto pfpListHandle = evt.getValidHandle<std::vector<recob::PFParticle> >(fPFParticleLabel[it]);
      //if (!pfpListHandle.isValid()) continue;

      for (auto const& trackPtr : trackPtrVec){
        //creating a data set for this track to check for selected hits
        FiveDMap loopdata;
        // get the tracks hits and metadata for the hits
        std::vector<art::Ptr<recob::Hit>> const& trackHits = trackToHits.at(trackPtr.key());
        std::vector<const recob::TrackHitMeta*> const& trackHitMetas = trackToHits.data(trackPtr.key());
        // loop over the track hits
        // will need to also find the same hit in hitPtrVec to get the associated wire
        for (size_t hitIdx = 0; hitIdx < trackHits.size(); ++hitIdx){
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
          {  
            // check channel etc to see if the hits are a match
            if (trackHit->Channel()   == hitPtr->Channel()   &&
                trackHit->StartTick() == hitPtr->StartTick() &&
                trackHit->EndTick()   == hitPtr->EndTick()   && 
                trackHit->PeakTime()  == hitPtr->PeakTime()  )
            {
              //MF_LOG_VERBATIM("WireModWireModMakeHists");
              //<< "Wire Found!";
              wirePtr = hitToWireAssns.at(hitPtr.key());
                
              //getting the wire plane and TPC info
              int planeID = hitPtr->WireID().Plane;
              int tpcID = (hitPtr->WireID().TPC) / 2;
              int cryoID = hitPtr->WireID().Cryostat;

              //Making a histogram of the hit information 
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
              double integral = 0;
              for (size_t bin = 1; bin < hitWidth + fstBuffer + endBuffer + 1; ++bin){
                float wvfmVal = (wirePtr->Signal())[fstTick - fstBuffer + bin - 1];
                integral += wvfmVal;
                waveformHist->SetBinContent(bin, wvfmVal);
              }
              
              // Calculate integral and FWHM and adding to vector of data
              double amplitude = waveformHist-> GetMaximum();
              double half_max = amplitude / 2.0;
              int bin1 = waveformHist->FindFirstBinAbove(half_max);
              int bin2 = waveformHist->FindLastBinAbove(half_max);
              double fwhm = waveformHist->GetBinCenter(bin2) - waveformHist->GetBinCenter(bin1);
              
              data[planeID][cryoID][tpcID][5][0].push_back(integral);
              loopdata[planeID][cryoID][tpcID][5][0].push_back(integral);
              data[planeID][cryoID][tpcID][6][0].push_back(fwhm);
              loopdata[planeID][cryoID][tpcID][6][0].push_back(fwhm);

              //getting the gaussian information
              double gaussintegral = hitPtr->Integral();
              double gausswidth = hitPtr->RMS();
              double gaussfit = hitPtr->GoodnessOfFit();
              if (fXCorrection)
              {
                TSpline3* integralSpline = fXCorrFile->Get<TSpline3>(("xIntegral_" + std::to_string(planeID) + "_ratio_prof").c_str());
                TSpline3* widthSpline  = fXCorrFile->Get<TSpline3>(("xWidth_"  + std::to_string(planeID) + "_ratio_prof").c_str());
                gaussintegral *= integralSpline->Eval(hitLoc.X());
                gausswidth  *= widthSpline ->Eval(hitLoc.X());
              }

              data[planeID][cryoID][tpcID][7][0].push_back(gausswidth);
              data[planeID][cryoID][tpcID][8][0].push_back(gaussintegral);
              data[planeID][cryoID][tpcID][9][0].push_back(gaussfit);
              loopdata[planeID][cryoID][tpcID][7][0].push_back(gausswidth);
              loopdata[planeID][cryoID][tpcID][8][0].push_back(gaussintegral);
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
        }
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
  DEFINE_ART_RESULTS_PLUGIN(WireModMakeHists)
} // end WireMod namespace
