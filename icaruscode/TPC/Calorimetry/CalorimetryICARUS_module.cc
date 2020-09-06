//////////////////////////////////////////////////
//
// Calorimetry class
//
// maddalena.antonello@lngs.infn.it
// ornella.palamara@lngs.infn.it
// ART port echurch@fnal.gov
//  This algorithm is designed to perform the calorimetric reconstruction
//  of the 3D reconstructed tracks
////////////////////////////////////////////////////////////////////////

#include <string>
#include <optional>
#include <cmath>
#include <limits> // std::numeric_limits<>

#include "icaruscode/TPC/Calorimetry/Algorithms/CalorimetryIcarusAlg.h"
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
#include "lardata/ArtDataHelper/TrackUtils.h" // lar::util::TrackPitchInView()
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "larcorealg/CoreUtils/NumericUtils.h" // util::absDiff()

#include "larevt/SpaceCharge/SpaceCharge.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

// ROOT includes
#include <TMath.h>
#include <TGraph.h>
#include <TF1.h>
#include <TVector3.h>

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/pow.h" // cet::sum_of_squares()

///calorimetry
namespace calo
{

/**
   * @brief Estimates the energy deposited by reconstructed tracks.
   *
   * Output
   * =======
   *
   * * `std::vector<anab::Calorimetry>`: collection of calorimetry information,
   *      per reconstructed track and per wire plane
   * * `art::Assns<recob::Track, anab::Calorimetry>` association of each track
   *      with its calorimetry information
   *
   *
   * Configuration
   * ==============
   *
   * @note This documentation is grossly incomplete.
   *
   * * **NotOnTrackZcut** (real, optional): if specified, hits associated to all
   *     trajectory points whose _z_ coordinate is below `NotOnTrackZcut` value
   *     (including electric field distortion correction if enabled)
   *     are excluded from the calorimetry. The value is specified as absolute
   *     _z_ coordinate in world reference frame, in centimeters.
   *     The legacy value of this cut was hard coded to `-100.0` cm.
   *
   *
   */
class CalorimetryICARUS : public art::EDProducer
{

public:
  explicit CalorimetryICARUS(fhicl::ParameterSet const &pset);

private:
  void produce(art::Event &evt) override;
  void ReadCaloTree();

  bool BeginsOnBoundary(art::Ptr<recob::Track> lar_track);
  bool EndsOnBoundary(art::Ptr<recob::Track> lar_track);

  void GetPitch(detinfo::DetectorPropertiesData const& detProp,
                art::Ptr<recob::Hit> hit, std::vector<double> trkx, std::vector<double> trky, std::vector<double> trkz, std::vector<double> trkw, std::vector<double> trkx0, double *xyz3d, double &pitch, double TickT0);

  std::string fTrackModuleLabel;
  std::string fSpacePointModuleLabel;
  std::string fT0ModuleLabel;
  bool fUseArea;
  bool fUseIntegral;
  bool fSCE;
  bool fFlipTrack_dQdx;                  //flip track direction if significant rise of dQ/dx at the track start
  std::optional<double> fNotOnTrackZcut; ///< Exclude trajectory points with _z_ lower than this [cm]
  CalorimetryIcarusAlg caloAlg;

  int fnsps;
  std::vector<int> fwire;
  std::vector<double> ftime;
  std::vector<double> fstime;
  std::vector<double> fetime;
  std::vector<double> fMIPs;
  std::vector<double> fdQdx;
  std::vector<double> fdEdx;
  std::vector<double> fResRng;
  std::vector<float> fpitch;
  std::vector<TVector3> fXYZ;
  std::vector<size_t> fHitIndex;

}; // class Calorimetry

} // namespace calo

//-------------------------------------------------
calo::CalorimetryICARUS::CalorimetryICARUS(fhicl::ParameterSet const &pset)
    : EDProducer{pset},
      fTrackModuleLabel(pset.get<std::string>("TrackModuleLabel")),
      fSpacePointModuleLabel(pset.get<std::string>("SpacePointModuleLabel")),
      fT0ModuleLabel(pset.get<std::string>("T0ModuleLabel")),
      fUseArea(pset.get<bool>("UseArea")),
      fUseIntegral(pset.get<bool>("UseIntegral",true)),
      fSCE(pset.get<bool>("CorrectSCE")),
      fFlipTrack_dQdx(pset.get<bool>("FlipTrack_dQdx", true)),
      caloAlg(pset.get<fhicl::ParameterSet>("CaloAlg"))
{
//std::cout <<" initializing calorimetry... " << std::endl;
  if (pset.has_key("NotOnTrackZcut"))
    fNotOnTrackZcut = pset.get<double>("NotOnTrackZcut");

  produces<std::vector<anab::Calorimetry>>();
  produces<art::Assns<recob::Track, anab::Calorimetry>>();
}

//------------------------------------------------------------------------------------//
void calo::CalorimetryICARUS::produce(art::Event &evt)
{
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clockData);
  //std::cout <<" producing calorimetry... " << std::endl;
  std::cout << " useintegral " << fUseIntegral << std::endl;
  auto const *sce = lar::providerFrom<spacecharge::SpaceChargeService>();

  art::Handle<std::vector<recob::Track>> trackListHandle;
  std::vector<art::Ptr<recob::Track>> tracklist;
  if (evt.getByLabel(fTrackModuleLabel, trackListHandle))
    art::fill_ptr_vector(tracklist, trackListHandle);

//std::cout << " track list size " << tracklist.size() << std::endl;

  // Get Geometry
  art::ServiceHandle<geo::Geometry const> geom;

  // channel quality
  lariov::ChannelStatusProvider const &channelStatus = art::ServiceHandle<lariov::ChannelStatusService const>()->GetProvider();

  size_t nplanes = geom->Nplanes();

  //create anab::Calorimetry objects and make association with recob::Track
  std::unique_ptr<std::vector<anab::Calorimetry>> calorimetrycol(new std::vector<anab::Calorimetry>);
  std::unique_ptr<art::Assns<recob::Track, anab::Calorimetry>> assn(new art::Assns<recob::Track, anab::Calorimetry>);

  //art::FindManyP<recob::SpacePoint> fmsp(trackListHandle, evt, fTrackModuleLabel);
  art::FindManyP<recob::Hit> fmht(trackListHandle, evt, fTrackModuleLabel);
  art::FindManyP<recob::Hit, recob::TrackHitMeta> fmthm(trackListHandle, evt, fTrackModuleLabel); //this has more information about hit-track association, only available in PMA for now
  art::FindManyP<anab::T0> fmt0(trackListHandle, evt, fT0ModuleLabel);

  for (size_t trkIter = 0; trkIter < tracklist.size(); ++trkIter)
  {

    decltype(auto) larEnd = tracklist[trkIter]->Trajectory().End();

    // Some variables for the hit
    float time;             //hit time at maximum
    float stime;            //hit start time
    float etime;            //hit end time
    uint32_t channel = 0;   //channel number
    unsigned int cstat = 0; //hit cryostat number
    unsigned int tpc = 0;   //hit tpc number
    unsigned int wire = 0;  //hit wire number
    unsigned int plane = 0; //hit plane number

    std::vector<art::Ptr<recob::Hit>> allHits = fmht.at(trkIter);
    double T0 = 0;
    double TickT0 = 0;
    if (fmt0.isValid())
    {
      std::vector<art::Ptr<anab::T0>> allT0 = fmt0.at(trkIter);
      if (allT0.size())
        T0 = allT0[0]->Time();
      TickT0 = T0 / sampling_rate(clockData);
    }

    std::vector<std::vector<unsigned int>> hits(nplanes);

    art::FindManyP<recob::SpacePoint> fmspts(allHits, evt, fSpacePointModuleLabel);
    for (size_t ah = 0; ah < allHits.size(); ++ah)
    {
      hits[allHits[ah]->WireID().Plane].push_back(ah);
    }
    //get hits in each plane
    for (size_t ipl = 0; ipl < nplanes; ++ipl)
    { //loop over all wire planes

      geo::PlaneID planeID; //(cstat,tpc,ipl);

      fwire.clear();
      ftime.clear();
      fstime.clear();
      fetime.clear();
      fMIPs.clear();
      fdQdx.clear();
      fdEdx.clear();
      fpitch.clear();
      fResRng.clear();
      fXYZ.clear();
      fHitIndex.clear();

      float Kin_En = 0.;
      float Trk_Length = 0.;
      std::vector<float> vdEdx;
      std::vector<float> vresRange;
      std::vector<float> vdQdx;
      std::vector<float> deadwire; //residual range for dead wires
      std::vector<TVector3> vXYZ;

      // Require at least 2 hits in this view
      if (hits[ipl].size() < 2)
      {
        if (hits[ipl].size() == 1)
        {
          mf::LogWarning("Calorimetry") << "Only one hit in plane " << ipl << " associated with track id " << trkIter;
        }
        calorimetrycol->push_back(anab::Calorimetry(util::kBogusD,
                                                    vdEdx,
                                                    vdQdx,
                                                    vresRange,
                                                    deadwire,
                                                    util::kBogusD,
                                                    fpitch,
                                                    recob::tracking::convertCollToPoint(vXYZ),
                                                    planeID));
        util::CreateAssn(*this, evt, *calorimetrycol, tracklist[trkIter], *assn);
        continue;
      }

      //range of wire signals
      unsigned int wire0 = 100000;
      unsigned int wire1 = 0;
      double PIDA = 0;
      int nPIDA = 0;

      // determine track direction. Fill residual range array
      bool GoingDS = true;
      // find the track direction by comparing US and DS charge BB
      double USChg = 0;
      double DSChg = 0;
      // temp array holding distance betweeen space points
      std::vector<double> spdelta;
      //int nht = 0; //number of hits
      fnsps = 0; //number of space points
      std::vector<double> ChargeBeg;
      std::stack<double> ChargeEnd;

      // find track pitch
      double fTrkPitch = 0;
      for (size_t itp = 0; itp < tracklist[trkIter]->NumberTrajectoryPoints(); ++itp)
      {

        const auto &pos = tracklist[trkIter]->LocationAtPoint(itp);
        const auto &dir = tracklist[trkIter]->DirectionAtPoint(itp);

        const double Position[3] = {pos.X(), pos.Y(), pos.Z()};
        geo::TPCID tpcid = geom->FindTPCAtPosition(Position);
        if (tpcid.isValid)
        {
          try
          {
            fTrkPitch = lar::util::TrackPitchInView(*tracklist[trkIter], geom->Plane(ipl).View(), itp);

            //Correct for SCE
            geo::Vector_t posOffsets = {0., 0., 0.};
            geo::Vector_t dirOffsets = {0., 0., 0.};
            if (sce->EnableCalSpatialSCE() && fSCE)
              posOffsets = sce->GetCalPosOffsets(geo::Point_t(pos), tpcid.TPC);
            if (sce->EnableCalSpatialSCE() && fSCE)
              dirOffsets = sce->GetCalPosOffsets(geo::Point_t{pos.X() + fTrkPitch * dir.X(), pos.Y() + fTrkPitch * dir.Y(), pos.Z() + fTrkPitch * dir.Z()}, tpcid.TPC);
            TVector3 dir_corr = {fTrkPitch * dir.X() - dirOffsets.X() + posOffsets.X(), fTrkPitch * dir.Y() + dirOffsets.Y() - posOffsets.Y(), fTrkPitch * dir.Z() + dirOffsets.Z() - posOffsets.Z()};

            fTrkPitch = dir_corr.Mag();
          }
          catch (cet::exception &e)
          {
            mf::LogWarning("Calorimetry") << "caught exception "
                                          << e << "\n setting pitch (C) to "
                                          << util::kBogusD;
            fTrkPitch = 0;
          }
          break;
        }
      }

      // find the separation between all space points
      double xx = 0., yy = 0., zz = 0.;

      //save track 3d points
      std::vector<double> trkx;
      std::vector<double> trky;
      std::vector<double> trkz;
      std::vector<double> trkw;
      std::vector<double> trkx0;
      for (size_t i = 0; i < hits[ipl].size(); ++i)
      {
        //Get space points associated with the hit
        std::vector<art::Ptr<recob::SpacePoint>> sptv = fmspts.at(hits[ipl][i]);
        for (size_t j = 0; j < sptv.size(); ++j)
        {

          double t = allHits[hits[ipl][i]]->PeakTime() - TickT0; // Want T0 here? Otherwise ticks to x is wrong?
          double x = detProp.ConvertTicksToX(t, allHits[hits[ipl][i]]->WireID().Plane, allHits[hits[ipl][i]]->WireID().TPC, allHits[hits[ipl][i]]->WireID().Cryostat);
          double w = allHits[hits[ipl][i]]->WireID().Wire;
          if (TickT0)
          {
            trkx.push_back(sptv[j]->XYZ()[0] - detProp.ConvertTicksToX(TickT0, allHits[hits[ipl][i]]->WireID().Plane, allHits[hits[ipl][i]]->WireID().TPC, allHits[hits[ipl][i]]->WireID().Cryostat));
          }
          else
          {
            trkx.push_back(sptv[j]->XYZ()[0]);
          }
          trky.push_back(sptv[j]->XYZ()[1]);
          trkz.push_back(sptv[j]->XYZ()[2]);
          trkw.push_back(w);
          trkx0.push_back(x);
        }
      }
      for (size_t ihit = 0; ihit < hits[ipl].size(); ++ihit)
      { //loop over all hits on each wire plane

        //std::cout<<ihit<<std::endl;

        if (!planeID.isValid)
        {
          plane = allHits[hits[ipl][ihit]]->WireID().Plane;
          tpc = allHits[hits[ipl][ihit]]->WireID().TPC;
          cstat = allHits[hits[ipl][ihit]]->WireID().Cryostat;
          planeID.Cryostat = cstat;
          planeID.TPC = tpc;
          planeID.Plane = plane;
          planeID.isValid = true;
        }

        wire = allHits[hits[ipl][ihit]]->WireID().Wire;
        time = allHits[hits[ipl][ihit]]->PeakTime(); // What about here? T0
        stime = allHits[hits[ipl][ihit]]->PeakTimeMinusRMS();
        etime = allHits[hits[ipl][ihit]]->PeakTimePlusRMS();
        const size_t &hitIndex = allHits[hits[ipl][ihit]].key();

        double charge = allHits[hits[ipl][ihit]]->PeakAmplitude();
        if (fUseArea)
        {
       //   std::cout << " integral " << allHits[hits[ipl][ihit]]->Integral() << std::endl;
//std::cout << " sumadc " << allHits[hits[ipl][ihit]]->SummedADC() << std::endl;
          if (fUseIntegral) charge = allHits[hits[ipl][ihit]]->Integral();
          else              charge = allHits[hits[ipl][ihit]]->SummedADC();
        }
//std::cout << " ipl " << ipl << " ihit " << ihit << " charge " << charge << std::endl;
        //get 3d coordinate and track pitch for the current hit
        //not all hits are associated with space points, the method uses neighboring spacepts to interpolate
        double xyz3d[3];
        double pitch;
        bool fBadhit = false;
        if (fmthm.isValid())
        {
          auto vhit = fmthm.at(trkIter);
          auto vmeta = fmthm.data(trkIter);
          for (size_t ii = 0; ii < vhit.size(); ++ii)
          {
            if (vhit[ii].key() == allHits[hits[ipl][ihit]].key())
            {
              if (vmeta[ii]->Index() == std::numeric_limits<int>::max())
              {
                fBadhit = true;
                continue;
              }
              if (vmeta[ii]->Index() >= tracklist[trkIter]->NumberTrajectoryPoints())
              {
                throw cet::exception("Calorimetry_module.cc") << "Requested track trajectory index " << vmeta[ii]->Index() << " exceeds the total number of trajectory points " << tracklist[trkIter]->NumberTrajectoryPoints() << " for track index " << trkIter << ". Something is wrong with the track reconstruction. Please contact tjyang@fnal.gov";
              }
              if (!tracklist[trkIter]->HasValidPoint(vmeta[ii]->Index()))
              {
                fBadhit = true;
                continue;
              }

              //Correct location for SCE
              geo::Point_t const loc = tracklist[trkIter]->LocationAtPoint(vmeta[ii]->Index());
              geo::Vector_t locOffsets = {
                  0.,
                  0.,
                  0.,
              };
              if (sce->EnableCalSpatialSCE() && fSCE)
                locOffsets = sce->GetCalPosOffsets(loc, vhit[ii]->WireID().TPC);
              xyz3d[0] = loc.X() - locOffsets.X();
              xyz3d[1] = loc.Y() + locOffsets.Y();
              xyz3d[2] = loc.Z() + locOffsets.Z();

              double angleToVert = geom->WireAngleToVertical(vhit[ii]->View(), vhit[ii]->WireID().TPC, vhit[ii]->WireID().Cryostat) - 0.5 * ::util::pi<>();
              const geo::Vector_t &dir = tracklist[trkIter]->DirectionAtPoint(vmeta[ii]->Index());
              double cosgamma = std::abs(std::sin(angleToVert) * dir.Y() + std::cos(angleToVert) * dir.Z());
              if (cosgamma)
              {
                pitch = geom->WirePitch(vhit[ii]->View()) / cosgamma;
              }
              else
              {
                pitch = 0;
              }

              //Correct pitch for SCE
              geo::Vector_t dirOffsets = {0., 0., 0.};
              if (sce->EnableCalSpatialSCE() && fSCE)
                dirOffsets = sce->GetCalPosOffsets(geo::Point_t{loc.X() + pitch * dir.X(), loc.Y() + pitch * dir.Y(), loc.Z() + pitch * dir.Z()}, vhit[ii]->WireID().TPC);
              const TVector3 &dir_corr = {pitch * dir.X() - dirOffsets.X() + locOffsets.X(), pitch * dir.Y() + dirOffsets.Y() - locOffsets.Y(), pitch * dir.Z() + dirOffsets.Z() - locOffsets.Z()};

              pitch = dir_corr.Mag();

              break;
            }
          }
        }
        else
          GetPitch(detProp, allHits[hits[ipl][ihit]], trkx, trky, trkz, trkw, trkx0, xyz3d, pitch, TickT0);

        if (fBadhit)
          continue;
        if (fNotOnTrackZcut && (xyz3d[2] < fNotOnTrackZcut.value()))
          continue; //hit not on track
        if (pitch <= 0)
          pitch = fTrkPitch;
        if (!pitch)
          continue;

        if (fnsps == 0)
        {
          xx = xyz3d[0];
          yy = xyz3d[1];
          zz = xyz3d[2];
          spdelta.push_back(0);
        }
        else
        {
          double dx = xyz3d[0] - xx;
          double dy = xyz3d[1] - yy;
          double dz = xyz3d[2] - zz;
          spdelta.push_back(sqrt(dx * dx + dy * dy + dz * dz));
          Trk_Length += spdelta.back();
          xx = xyz3d[0];
          yy = xyz3d[1];
          zz = xyz3d[2];
        }

        ChargeBeg.push_back(charge);
        ChargeEnd.push(charge);

        double MIPs = charge;
        double dQdx = MIPs / pitch;
        double dEdx = 0;
        if (fUseArea) {
                   if(fUseIntegral) dEdx = caloAlg.dEdx_AREA(clockData, detProp, *allHits[hits[ipl][ihit]], pitch, T0);
                   else dEdx = caloAlg.dEdx_SUMADC(clockData, detProp, allHits[hits[ipl][ihit]], pitch, T0);
//std::cout << " ipl " << ipl << " ihit " << ihit << " charge " << charge << std::endl;
}
        else
          dEdx = caloAlg.dEdx_AMP(clockData, detProp, *allHits[hits[ipl][ihit]], pitch, T0);

        Kin_En = Kin_En + dEdx * pitch;

        if (allHits[hits[ipl][ihit]]->WireID().Wire < wire0)
          wire0 = allHits[hits[ipl][ihit]]->WireID().Wire;
        if (allHits[hits[ipl][ihit]]->WireID().Wire > wire1)
          wire1 = allHits[hits[ipl][ihit]]->WireID().Wire;

        fMIPs.push_back(MIPs);
        fdEdx.push_back(dEdx);
        fdQdx.push_back(dQdx);
        fwire.push_back(wire);
        ftime.push_back(time);
        fstime.push_back(stime);
        fetime.push_back(etime);
        fpitch.push_back(pitch);
        TVector3 v(xyz3d[0], xyz3d[1], xyz3d[2]);
        //std::cout << "Adding these positions to v and then fXYZ " << xyz3d[0] << " " << xyz3d[1] << " " << xyz3d[2] << "\n" <<std::endl;
        fXYZ.push_back(v);
        fHitIndex.push_back(hitIndex);
        ++fnsps;
      }
      if (fnsps < 2)
      {
        vdEdx.clear();
        vdQdx.clear();
        vresRange.clear();
        deadwire.clear();
        fpitch.clear();
        //std::cout << "Adding the aforementioned positions..." << std::endl;
        calorimetrycol->push_back(anab::Calorimetry(util::kBogusD,
                                                    vdEdx,
                                                    vdQdx,
                                                    vresRange,
                                                    deadwire,
                                                    util::kBogusD,
                                                    fpitch,
                                                    recob::tracking::convertCollToPoint(vXYZ),
                                                    planeID));
        util::CreateAssn(*this, evt, *calorimetrycol, tracklist[trkIter], *assn);
        continue;
      }
      for (int isp = 0; isp < fnsps; ++isp)
      {
        if (isp > 3)
          break;
        USChg += ChargeBeg[isp];
      }
      int countsp = 0;
      while (!ChargeEnd.empty())
      {
        if (countsp > 3)
          break;
        DSChg += ChargeEnd.top();
        ChargeEnd.pop();
        ++countsp;
      }
      if (fFlipTrack_dQdx)
      {
        // Going DS if charge is higher at the end
        GoingDS = (DSChg > USChg);
      }
      else
      {
        // Use the track direction to determine the residual range
        if (!fXYZ.empty())
        {
          TVector3 track_start(tracklist[trkIter]->Trajectory().Vertex().X(),
                               tracklist[trkIter]->Trajectory().Vertex().Y(),
                               tracklist[trkIter]->Trajectory().Vertex().Z());
          TVector3 track_end(tracklist[trkIter]->Trajectory().End().X(),
                             tracklist[trkIter]->Trajectory().End().Y(),
                             tracklist[trkIter]->Trajectory().End().Z());

          if ((fXYZ[0] - track_start).Mag() + (fXYZ.back() - track_end).Mag() <
              (fXYZ[0] - track_end).Mag() + (fXYZ.back() - track_start).Mag())
          {
            GoingDS = true;
          }
          else
          {
            GoingDS = false;
          }
        }
      }

      // determine the starting residual range and fill the array
      fResRng.resize(fnsps);
      if (fResRng.size() < 2 || spdelta.size() < 2)
      {
        mf::LogWarning("Calorimetry") << "fResrng.size() = " << fResRng.size() << " spdelta.size() = " << spdelta.size();
      }
      if (GoingDS)
      {
        fResRng[fnsps - 1] = spdelta[fnsps - 1] / 2;
        for (int isp = fnsps - 2; isp > -1; isp--)
        {
          fResRng[isp] = fResRng[isp + 1] + spdelta[isp + 1];
        }
      }
      else
      {
        fResRng[0] = spdelta[1] / 2;
        for (int isp = 1; isp < fnsps; isp++)
        {
          fResRng[isp] = fResRng[isp - 1] + spdelta[isp];
        }
      }

      MF_LOG_DEBUG("CaloPrtHit") << " pt wire  time  ResRng    MIPs   pitch   dE/dx    Ai X Y Z\n";

      double Ai = -1;
      for (int i = 0; i < fnsps; ++i)
      { //loop over all 3D points
        vresRange.push_back(fResRng[i]);
        vdEdx.push_back(fdEdx[i]);
        vdQdx.push_back(fdQdx[i]);
        vXYZ.push_back(fXYZ[i]);
        if (i != 0 && i != fnsps - 1)
        { //ignore the first and last point
          // Calculate PIDA
          Ai = fdEdx[i] * pow(fResRng[i], 0.42);
          nPIDA++;
          PIDA += Ai;
        }

        MF_LOG_DEBUG("CaloPrtHit") << std::setw(4) << trkIter
                                   //std::cout<<std::setw(4)<< trkIter
                                   << std::setw(4) << ipl
                                   << std::setw(4) << i
                                   << std::setw(4) << fwire[i]
                                   << std::setw(6) << (int)ftime[i]
                                   << std::setiosflags(std::ios::fixed | std::ios::showpoint)
                                   << std::setprecision(2)
                                   << std::setw(8) << fResRng[i]
                                   << std::setprecision(1)
                                   << std::setw(8) << fMIPs[i]
                                   << std::setprecision(2)
                                   << std::setw(8) << fpitch[i]
                                   << std::setw(8) << fdEdx[i]
                                   << std::setw(8) << Ai
                                   << std::setw(8) << fXYZ[i].x()
                                   << std::setw(8) << fXYZ[i].y()
                                   << std::setw(8) << fXYZ[i].z()
                                   << "\n";
      } //end looping over 3D points
      if (nPIDA > 0)
      {
        PIDA = PIDA / (double)nPIDA;
      }
      else
      {
        PIDA = -1;
      }
      MF_LOG_DEBUG("CaloPrtTrk") << "Plane # " << ipl
                                 << "TrkPitch= "
                                 << std::setprecision(2) << fTrkPitch
                                 << " nhits= " << fnsps
                                 << "\n"
                                 << std::setiosflags(std::ios::fixed | std::ios::showpoint)
                                 << "Trk Length= " << std::setprecision(1)
                                 << Trk_Length << " cm,"
                                 << " KE calo= " << std::setprecision(1)
                                 << Kin_En << " MeV,"
                                 << " PIDA= " << PIDA
                                 << "\n";

      // look for dead wires
      for (unsigned int iw = wire0; iw < wire1 + 1; ++iw)
      {
        plane = allHits[hits[ipl][0]]->WireID().Plane;
        tpc = allHits[hits[ipl][0]]->WireID().TPC;
        cstat = allHits[hits[ipl][0]]->WireID().Cryostat;
        channel = geom->PlaneWireToChannel(plane, iw, tpc, cstat);
        if (channelStatus.IsBad(channel))
        {
          MF_LOG_DEBUG("Calorimetry") << "Found dead wire at Plane = " << plane
                                      << " Wire =" << iw;
          unsigned int closestwire = 0;
          unsigned int endwire = 0;
          unsigned int dwire = 100000;
          double mindis = 100000;
          double goodresrange = 0;
          //hitCtr = 0;
          for (size_t ihit = 0; ihit < hits[ipl].size(); ++ihit)
          {
            //	for(art::PtrVector<recob::Hit>::const_iterator hitIter = hitsV.begin();
            //	    hitIter != hitsV.end();
            //	    ++hitCtr, hitIter++){
            channel = allHits[hits[ipl][ihit]]->Channel();
            if (channelStatus.IsBad(channel))
              continue;
            // grab the space points associated with this hit
            std::vector<art::Ptr<recob::SpacePoint>> sppv = fmspts.at(hits[ipl][ihit]);
            if (sppv.size() < 1)
              continue;
            // only use the first space point in the collection, really each hit should
            // only map to 1 space point
            const recob::Track::Point_t xyz{sppv[0]->XYZ()[0],
                                            sppv[0]->XYZ()[1],
                                            sppv[0]->XYZ()[2]};
            double dis1 = (larEnd - xyz).Mag2();
            if (dis1)
              dis1 = std::sqrt(dis1);
            if (dis1 < mindis)
            {
              endwire = allHits[hits[ipl][ihit]]->WireID().Wire;
              mindis = dis1;
            }
            if (util::absDiff(wire, iw) < dwire)
            {
              closestwire = allHits[hits[ipl][ihit]]->WireID().Wire;
              dwire = util::absDiff(allHits[hits[ipl][ihit]]->WireID().Wire, iw);
              goodresrange = dis1;
            }
          }
          if (closestwire)
          {
            if (iw < endwire)
            {
              deadwire.push_back(goodresrange + (int(closestwire) - int(iw)) * fTrkPitch);
            }
            else
            {
              deadwire.push_back(goodresrange + (int(iw) - int(closestwire)) * fTrkPitch);
            }
          }
        }
      }
      //std::cout << "Adding at the end but still same fXYZ" << std::endl;
      calorimetrycol->push_back(anab::Calorimetry(Kin_En,
                                                  vdEdx,
                                                  vdQdx,
                                                  vresRange,
                                                  deadwire,
                                                  Trk_Length,
                                                  fpitch,
                                                  recob::tracking::convertCollToPoint(vXYZ),
                                                  fHitIndex,
                                                  planeID));
      util::CreateAssn(*this, evt, *calorimetrycol, tracklist[trkIter], *assn);

    } //end looping over planes
  }   //end looping over tracks

  evt.put(std::move(calorimetrycol));
  evt.put(std::move(assn));

  return;
}

void calo::CalorimetryICARUS::GetPitch(detinfo::DetectorPropertiesData const& detProp,
                                       art::Ptr<recob::Hit> hit, std::vector<double> trkx, std::vector<double> trky, std::vector<double> trkz, std::vector<double> trkw, std::vector<double> trkx0, double *xyz3d, double &pitch, double TickT0)
{
  //Get 3d coordinates and track pitch for each hit
  //Find 5 nearest space points and determine xyz and curvature->track pitch

  //std::cout << "Start of get pitch" << std::endl;

  // Get services
  art::ServiceHandle<geo::Geometry const> geom;
  auto const *sce = lar::providerFrom<spacecharge::SpaceChargeService>();

  //save distance to each spacepoint sorted by distance
  std::map<double, size_t> sptmap;
  //save the sign of distance
  std::map<size_t, int> sptsignmap;

  double wire_pitch = geom->WirePitch(0);

  double t0 = hit->PeakTime() - TickT0;
  double x0 = detProp.ConvertTicksToX(t0, hit->WireID().Plane, hit->WireID().TPC, hit->WireID().Cryostat);
  double w0 = hit->WireID().Wire;

  for (size_t i = 0; i < trkx.size(); ++i)
  {
    double distance = cet::sum_of_squares((trkw[i] - w0) * wire_pitch, trkx0[i] - x0);
    if (distance > 0)
      distance = sqrt(distance);
    //std::cout << "Dis " << distance << ", sqaured " << distance*distance << " = " << wire_pitch*wire_pitch <<"("<<trkw[i]<<"-"<<w0<<")^2 + ("<<trkx0[i]<<"-"<<x0<<")^2"<<std::endl;
    sptmap.insert(std::pair<double, size_t>(distance, i));
    if (w0 - trkw[i] > 0)
      sptsignmap.insert(std::pair<size_t, int>(i, 1));
    else
      sptsignmap.insert(std::pair<size_t, int>(i, -1));
  }

  //x,y,z vs distance
  std::vector<double> vx;
  std::vector<double> vy;
  std::vector<double> vz;
  std::vector<double> vs;

  double kx = 0, ky = 0, kz = 0;

  int np = 0;
  for (auto isp = sptmap.begin(); isp != sptmap.end(); isp++)
  {
    //    const double *xyz = new double[3];
    //    xyz = isp->second->XYZ();
    double xyz[3];
    xyz[0] = trkx[isp->second];
    xyz[1] = trky[isp->second];
    xyz[2] = trkz[isp->second];

    double distancesign = sptsignmap[isp->second];
    //std::cout<<np<<" "<<xyz[0]<<" "<<xyz[1]<<" "<<xyz[2]<<" "<<(*isp).first<<std::endl;
    if (np == 0 && isp->first > 30)
    { //hit not on track
      xyz3d[0] = std::numeric_limits<double>::lowest();
      xyz3d[1] = std::numeric_limits<double>::lowest();
      xyz3d[2] = std::numeric_limits<double>::lowest();
      pitch = -1;
      return;
    }
    //std::cout<<np<<" "<<xyz[0]<<" "<<xyz[1]<<" "<<xyz[2]<<" "<<(*isp).first<<" Plane " << hit->WireID().Plane << " TPC " << hit->WireID().TPC << std::endl;
    if (np < 5)
    {
      vx.push_back(xyz[0]);
      vy.push_back(xyz[1]);
      vz.push_back(xyz[2]);
      vs.push_back(isp->first * distancesign);
    }
    else
    {
      break;
    }
    np++;
    //delete [] xyz;
  }
  //std::cout<<"np="<<np<<std::endl;
  if (np >= 2)
  { //at least two points
    //std::cout << "At least 2 points.."<<std::endl;
    TGraph *xs = new TGraph(np, &vs[0], &vx[0]);
    //for (int i = 0; i<np; i++) std::cout<<i<<" "<<vs[i]<<" "<<vx[i]<<" "<<vy[i]<<" "<<vz[i]<<std::endl;
    try
    {
      if (np > 2)
      {
        xs->Fit("pol2", "Q");
      }
      else
      {
        xs->Fit("pol1", "Q");
      }
      TF1 *pol = 0;
      if (np > 2)
        pol = (TF1 *)xs->GetFunction("pol2");
      else
        pol = (TF1 *)xs->GetFunction("pol1");
      xyz3d[0] = pol->Eval(0);
      kx = pol->GetParameter(1);
      //std::cout<<"X fit "<<xyz3d[0]<<" "<<kx<<std::endl;
    }
    catch (...)
    {
      mf::LogWarning("Calorimetry::GetPitch") << "Fitter failed";
      xyz3d[0] = vx[0];
    }
    delete xs;
    TGraph *ys = new TGraph(np, &vs[0], &vy[0]);
    try
    {
      if (np > 2)
      {
        ys->Fit("pol2", "Q");
      }
      else
      {
        ys->Fit("pol1", "Q");
      }
      TF1 *pol = 0;
      if (np > 2)
        pol = (TF1 *)ys->GetFunction("pol2");
      else
        pol = (TF1 *)ys->GetFunction("pol1");
      xyz3d[1] = pol->Eval(0);
      ky = pol->GetParameter(1);
      //std::cout<<"Y fit "<<xyz3d[1]<<" "<<ky<<std::endl;
    }
    catch (...)
    {
      mf::LogWarning("Calorimetry::GetPitch") << "Fitter failed";
      xyz3d[1] = vy[0];
    }
    delete ys;
    TGraph *zs = new TGraph(np, &vs[0], &vz[0]);
    try
    {
      if (np > 2)
      {
        zs->Fit("pol2", "Q");
      }
      else
      {
        zs->Fit("pol1", "Q");
      }
      TF1 *pol = 0;
      if (np > 2)
        pol = (TF1 *)zs->GetFunction("pol2");
      else
        pol = (TF1 *)zs->GetFunction("pol1");
      xyz3d[2] = pol->Eval(0);
      kz = pol->GetParameter(1);
      //std::cout<<"Z fit "<<xyz3d[2]<<" "<<kz<<std::endl;
    }
    catch (...)
    {
      mf::LogWarning("Calorimetry::GetPitch") << "Fitter failed";
      xyz3d[2] = vz[0];
    }
    delete zs;
  }
  else if (np)
  {
    xyz3d[0] = vx[0];
    xyz3d[1] = vy[0];
    xyz3d[2] = vz[0];
  }
  else
  {
    xyz3d[0] = std::numeric_limits<double>::lowest();
    xyz3d[1] = std::numeric_limits<double>::lowest();
    xyz3d[2] = std::numeric_limits<double>::lowest();
    pitch = -1;
    return;
  }
  pitch = -1;
  if (kx * kx + ky * ky + kz * kz)
  {
    double tot = sqrt(kx * kx + ky * ky + kz * kz);
    kx /= tot;
    ky /= tot;
    kz /= tot;
    //get pitch
    double wirePitch = geom->WirePitch(hit->WireID().Plane, hit->WireID().TPC, hit->WireID().Cryostat);
    double angleToVert = geom->Plane(hit->WireID().Plane, hit->WireID().TPC, hit->WireID().Cryostat).Wire(0).ThetaZ(false) - 0.5 * TMath::Pi();
    double cosgamma = TMath::Abs(TMath::Sin(angleToVert) * ky + TMath::Cos(angleToVert) * kz);
    if (cosgamma > 0)
      pitch = wirePitch / cosgamma;

    //Correct for SCE
    geo::Vector_t posOffsets = {0., 0., 0.};
    geo::Vector_t dirOffsets = {0., 0., 0.};
    if (sce->EnableCalSpatialSCE() && fSCE)
      posOffsets = sce->GetCalPosOffsets(geo::Point_t{xyz3d[0], xyz3d[1], xyz3d[2]}, hit->WireID().TPC);
    if (sce->EnableCalSpatialSCE() && fSCE)
      dirOffsets = sce->GetCalPosOffsets(geo::Point_t{xyz3d[0] + pitch * kx, xyz3d[1] + pitch * ky, xyz3d[2] + pitch * kz}, hit->WireID().TPC);

    xyz3d[0] = xyz3d[0] - posOffsets.X();
    xyz3d[1] = xyz3d[1] + posOffsets.Y();
    xyz3d[2] = xyz3d[2] + posOffsets.Z();

    //Correct pitch for SCE
    TVector3 dir = {pitch * kx - dirOffsets.X() + posOffsets.X(), pitch * ky + dirOffsets.Y() - posOffsets.Y(), pitch * kz + dirOffsets.Z() - posOffsets.Z()};
    pitch = dir.Mag();
  }
  //std::cout << "At end of get pitch " << xyz3d[0] << " " << xyz3d[1] << " " << xyz3d[2] << " " << x0 << " " << std::endl;
}

namespace calo
{

DEFINE_ART_MODULE(CalorimetryICARUS)

} // namespace calo
