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

namespace TrackSplit {
  class TrackSplitAngle : public art::ResultsProducer {
  public:
    explicit TrackSplitAngle(fhicl::ParameterSet const& pset);
    ~TrackSplitAngle() override = default;
    void event(art::Event const& evt) override;
    void reconfigure(fhicl::ParameterSet const& p);
    void writeResults(art::Results& r) override;
    void clear() override;

  private:
    const double fPi = std::acos(-1);
    const geo::GeometryCore* fGeometry = lar::providerFrom<geo::Geometry>(); // get the geometry
    art::InputTag fTag;
    TH1D* fHist;
    TH2D* fHist2d;
  };//end TrackSplitAngle class

  //-------------------------------------------------------------
  TrackSplitAngle::TrackSplitAngle(fhicl::ParameterSet const& pset)
  {
    this->reconfigure(pset);
  }

  //-------------------------------------------------------------
  void TrackSplitAngle::reconfigure(fhicl::ParameterSet const& pset)
  {
    fTag = pset.get<art::InputTag>("Label", "pandoraTrackGausCryoE");

    art::ServiceHandle<art::TFileService> tfs;
    fHist   = tfs->make<TH1D>(("AnglesBtwnTracks_"     +fTag.label()).c_str(), ";Opening Angle Between Tracks (rad);Pairs of Tracks" , 360, 0, fPi);
    fHist2d = tfs->make<TH2D>(("AnglesBtwnTracksVsGap_"+fTag.label()).c_str(), ";Opening Angle Between Tracks (rad);Size of Gap (cm)", 360, 0, fPi, 2000, 0, 2000);
  }

  //-------------------------------------------------------------
  void TrackSplitAngle::event(art::Event const& evt)
  {
    // get the tracks
    art::Handle<std::vector<recob::Track>> trackHandle;
    evt.getByLabel(fTag, trackHandle);
    art::FindManyP<recob::Hit, recob::TrackHitMeta> hitsFromTracks(trackHandle, evt, fTag);

    // loop over the tracks
    for (size_t idx_1 = 0; idx_1 < trackHandle->size(); ++idx_1)
    {
      // get the first track
      art::Ptr<recob::Track> track1(trackHandle, idx_1);

      // get the track info
      std::vector<art::Ptr<recob::Hit>> trackHits1 = hitsFromTracks.at(track1.key());
      std::vector<const recob::TrackHitMeta*> trackMetas1 = hitsFromTracks.data(track1.key());
      recob::Track::Vector_t startDir1 = track1->StartDirection();
      recob::Track::Vector_t   endDir1 = track1->EndDirection();
      recob::Track::Point_t  startPos1 = track1->Start();
      recob::Track::Point_t    endPos1 = track1->End();
      double                   length1 = track1->Length();

      // DEMO - check the channels the end of track1 is on (check only the first track to not print too much)
      // iterate over the plane geometries to find the channel in that plane which sees the point
      if (idx_1 == 0)
      {
        mf::LogVerbatim("Demo")
          << "****************************************************";
        for (auto const& planeGeom : fGeometry->Iterate<geo::PlaneGeo>())
        {
          // try to find the nearest channel, catch exceptions where it can't find one
          try
          {
            raw::ChannelID_t endChannel = fGeometry->NearestChannel(endPos1, planeGeom.ID());
            mf::LogVerbatim("Demo")
              << planeGeom.ID().toString() << " has channel " << endChannel << " nearest to the end of track 1";
          }
          catch(const std::exception& e)
          {
            mf::LogVerbatim("Demo")
              << planeGeom.ID().toString() << " does not have a channel which is near the end of track1" << '\n'
              << "  exception:" << '\n' << e.what();
          }
        }
        mf::LogVerbatim("Demo")
          << "****************************************************";
      }
      // END DEMO
       

      // if the track is short, skip it
      if (length1 < 20)
        continue;
      
      // get start and stop ticks for track
      // loop over hits and store min/max peak time
      float peakTimeMin1 = std::numeric_limits<float>::max();
      float peakTimeMax1 = std::numeric_limits<float>::min();
      for (size_t hitIdx = 0; hitIdx < trackHits1.size(); ++hitIdx)
      {
        art::Ptr<recob::Hit> hit = trackHits1.at(hitIdx);
        const recob::TrackHitMeta* meta = trackMetas1.at(hitIdx);
        size_t hitInTrack = meta->Index();
        if (hitInTrack == std::numeric_limits<unsigned int>::max() || not track1->HasValidPoint(hitInTrack))
          continue;
        peakTimeMin1 = (peakTimeMin1 > hit->PeakTime()) ? hit->PeakTime() : peakTimeMin1;
        peakTimeMax1 = (peakTimeMax1 < hit->PeakTime()) ? hit->PeakTime() : peakTimeMax1;
      }

      // idx_1 sets the first track, pair off with all sub
      for (size_t idx_2 = idx_1 + 1; idx_2 < trackHandle->size(); ++idx_2)
      {
        // get the second track
        art::Ptr<recob::Track> track2(trackHandle, idx_2);

        // get the track info
        std::vector<art::Ptr<recob::Hit>> trackHits2 = hitsFromTracks.at(track2.key());
        std::vector<const recob::TrackHitMeta*> trackMetas2 = hitsFromTracks.data(track2.key());
        recob::Track::Vector_t startDir2 = track2->StartDirection();
        recob::Track::Vector_t   endDir2 = track2->EndDirection();
        recob::Track::Point_t  startPos2 = track2->Start();
        recob::Track::Point_t    endPos2 = track2->End();
        double                   length2 = track2->Length();

        // if the track is short, skip it
        if (length2 < 20)
          continue;

        // get start and stop ticks for track
        // loop over hits and store min/max peak time
        float peakTimeMin2 = std::numeric_limits<float>::max();
        float peakTimeMax2 = std::numeric_limits<float>::min();
        for (size_t hitIdx = 0; hitIdx < trackHits2.size(); ++hitIdx)
        {
          art::Ptr<recob::Hit> hit = trackHits2.at(hitIdx);
          const recob::TrackHitMeta* meta = trackMetas2.at(hitIdx);
          size_t hitInTrack = meta->Index();
          if (hitInTrack == std::numeric_limits<unsigned int>::max() || not track2->HasValidPoint(hitInTrack))
            continue;
          peakTimeMin2 = (peakTimeMin2 > hit->PeakTime()) ? hit->PeakTime() : peakTimeMin2;
          peakTimeMax2 = (peakTimeMax2 < hit->PeakTime()) ? hit->PeakTime() : peakTimeMax2;
        }

        // check if the start/end dirs are in the same direction
        // eg, the dot products are positive
        bool point1to2 = (endDir1.Dot(startDir2) > 0) && (endDir1.Dot(startPos2 - endPos1) > 0);
        bool point2to1 = (endDir2.Dot(startDir1) > 0) && (endDir2.Dot(startPos1 - endPos2) > 0);

        // check if track 1 comes before track 2 or vice versa
        bool track1before2 = (peakTimeMax1 < peakTimeMin2);
        bool track2before1 = (peakTimeMax2 < peakTimeMin1);

        if (point1to2 && track1before2)
        {
          // these are useful in both cases
          // don't recall if the directions are unit vectors, but best to play it safe
          recob::Track::Vector_t unit1 = endDir1.Unit();
          recob::Track::Vector_t unit2 = startDir2.Unit();
          recob::Track::Vector_t unitC = unit1.Cross(unit2);
          recob::Track::Vector_t trackDis = (startPos2 - endPos1);
          // if the gap is greater than 100 cm, skip it
          if (std::sqrt(trackDis.Mag2()) > 100)
            continue;
          // if the the directions are parallel, take the angle between the common direction and the line connecting the start and end
          // subtracks this from Pi to approximate a difflection angle. Not perfect, but this is an edge case and the choice goes to Pi
          // as the tracks line up with each other
          double angle = std::numeric_limits<double>::max();
          if (unitC.Mag2() == 0)
          {
            angle = fPi - std::acos( endDir1.Dot(trackDis) / std::sqrt(endDir1.Mag2() * trackDis.Mag2()));
            fHist->Fill(angle);
            fHist2d->Fill(angle, std::sqrt(trackDis.Mag2()));
          } else {
            // not parallel, so we follow Ronald Goldman on p.304 of _Graphics Gems_ to find the points of closest approach on the two lines
            double t1 = ( trackDis.X() * (unit2.Y()*unitC.Z() - unit2.Z()*unitC.Y())
                        - trackDis.Y() * (unit2.X()*unitC.Z() - unit2.Z()*unitC.X())
                        + trackDis.Z() * (unit2.X()*unitC.Y() - unit2.Y()*unitC.X()) ) / unitC.Mag2();
            double t2 = ( trackDis.X() * (unit1.Y()*unitC.Z() - unit1.Z()*unitC.Y())
                        - trackDis.Y() * (unit1.X()*unitC.Z() - unit1.Z()*unitC.X())
                        + trackDis.Z() * (unit1.X()*unitC.Y() - unit1.Y()*unitC.X()) ) / unitC.Mag2();
            // the points of closest approach are endPos1 + t1*unit1 and startPos2 + t2*unit2
            // if the lines intersect these points are identical. In any case take the midpoint
            recob::Track::Point_t closest1 = (  endPos1 + unit1*t1);
            recob::Track::Point_t closest2 = (startPos2 + unit2*t2);
            recob::Track::Vector_t closestSep = (closest2 - closest1);
            recob::Track::Point_t midpoint = closest1 + 0.5*closestSep; // in a just world I could just average recob::Track::Point_t...
            recob::Track::Vector_t startToMid = (midpoint - startPos2);
            recob::Track::Vector_t   endToMid = (midpoint -   endPos1);
            angle = std::acos( startToMid.Dot(endToMid) / std::sqrt(startToMid.Mag2() * endToMid.Mag2()));
            fHist->Fill(angle);
            fHist2d->Fill(angle, std::sqrt(trackDis.Mag2()));
          }
          if (angle > (8.0 * fPi / 9.0))
            mf::LogVerbatim("TrackSplitAngle")
              << "Run " << evt.run() << ", Subrun " << evt.subRun() << ", event " << evt.event() << '\n'
              << "  Point 1 to 2: Track split with angle of " << angle << " rad, and gap " << std::sqrt(trackDis.Mag2()) << " cm" << '\n'
              << "  Track 1: Length " << length1 << " cm, Theta " << track1->Theta() << ", Phi "<< track1->Phi() << '\n'
              << "  Track 2: Length " << length2 << " cm, Theta " << track2->Theta() << ", Phi "<< track2->Phi();
        } else if (point2to1 && track2before1) {
          // these are useful in both cases
          // don't recall if the directions are unit vectors, but best to play it safe
          recob::Track::Vector_t unit1 = startDir1.Unit();
          recob::Track::Vector_t unit2 =   endDir2.Unit();
          recob::Track::Vector_t unitC = unit1.Cross(unit2);
          recob::Track::Vector_t trackDis = (endPos2 - startPos1);
          // if the gap is greater than 100 cm, skip it
          if (std::sqrt(trackDis.Mag2()) > 100)
            continue;
          // if the the directions are parallel, take the angle between the common direction and the line connecting the start and end
          // subtracks this from Pi to approximate a difflection angle. Not perfect, but this is an edge case and the choice goes to Pi
          // as the tracks line up with each other
          double angle = std::numeric_limits<double>::max();
          if (unitC.Mag2() == 0)
          {
            angle = fPi - std::acos( endDir2.Dot(trackDis) / std::sqrt(endDir2.Mag2() * trackDis.Mag2()));
            fHist->Fill(angle);
            fHist2d->Fill(angle, std::sqrt(trackDis.Mag2()));
          } else {
            // not parallel, so we follow Ronald Goldman on p.304 of _Graphics Gems_ to find the points of closest approach on the two lines
            double t1 = ( trackDis.X() * (unit2.Y()*unitC.Z() - unit2.Z()*unitC.Y())
                        - trackDis.Y() * (unit2.X()*unitC.Z() - unit2.Z()*unitC.X())
                        + trackDis.Z() * (unit2.X()*unitC.Y() - unit2.Y()*unitC.X()) ) / unitC.Mag2();
            double t2 = ( trackDis.X() * (unit1.Y()*unitC.Z() - unit1.Z()*unitC.Y())
                        - trackDis.Y() * (unit1.X()*unitC.Z() - unit1.Z()*unitC.X())
                        + trackDis.Z() * (unit1.X()*unitC.Y() - unit1.Y()*unitC.X()) ) / unitC.Mag2();
            // the points of closest approach are endPos1 + t1*unit1 and startPos2 + t2*unit2
            // if the lines intersect these points are identical. In any case take the midpoint
            recob::Track::Point_t closest1 = (startPos1 + unit1*t1);
            recob::Track::Point_t closest2 = (  endPos2 + unit2*t2);
            recob::Track::Vector_t closestSep = (closest2 - closest1);
            recob::Track::Point_t midpoint = closest1 + 0.5*closestSep; // in a just world I could just average recob::Track::Point_t...
            recob::Track::Vector_t startToMid = (midpoint - startPos1);
            recob::Track::Vector_t   endToMid = (midpoint -   endPos2);
            angle = std::acos( startToMid.Dot(endToMid) / std::sqrt(startToMid.Mag2() * endToMid.Mag2()));
            fHist->Fill(angle);
            fHist2d->Fill(angle, std::sqrt(trackDis.Mag2()));
          }
          if (angle > (8.0 * fPi / 9.0))
            mf::LogVerbatim("TrackSplitAngle")
              << "Run " << evt.run() << ", Subrun " << evt.subRun() << ", event " << evt.event() << '\n'
              << "  Point 2 to 1: Track split with angle of " << angle << " rad, and gap " << std::sqrt(trackDis.Mag2()) << " cm" << '\n'
              << "  Track 1: Length " << length1 << " cm, Theta " << track1->Theta() << ", Phi "<< track1->Phi() << '\n'
              << "  Track 2: Length " << length2 << " cm, Theta " << track2->Theta() << ", Phi "<< track2->Phi();
        }
      }
    }
  }

  //-------------------------------------------------------------
  void TrackSplitAngle::writeResults(art::Results& r)
  {
  }

  //-------------------------------------------------------------
  void TrackSplitAngle::clear()
  {
  }

  DEFINE_ART_RESULTS_PLUGIN(TrackSplitAngle)
}//end TrackSplit namespace
