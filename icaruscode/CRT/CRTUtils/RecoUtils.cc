#include <fstream>
#include <cetlib/search_path.h>
#include "RecoUtils.h"


namespace icarus::crt::dataTools{

TopCRTCentersMap LoadTopCRTCenters()
{
  TopCRTCentersMap TopCRTCenters;
  //  The follwing numbers have been extracted from the Top CRT modules geometry.
  //  The numbers are respectively moduleID and X,Y,Z coordinates of the
  //  module center (X markes the spot in the sketch below).
  //  The CRT modules are perfect square with 8 bars per side.
  //  The fixed coordinate (e.g. Y for Top CRT Horizontal) is returned from geometry.
  //  The transverce coordinate is returned from the average position of P1 and P2
  //  the average position of P1 and P3.
  //  The transformed CRT Hits are in cm.
  //
  //          ---------------------------------------------------------
  //          |      |      |      |      |      |      |      |      |       
  //          |      |      |      |      |      |      |      |      |       
  //          |      |      |      |      |      |      |      |      |       
  //          ---------------------------------------------------------
  //          |      |      |      |      |      |      |      |      |       
  //          |      |      |      |      |      |      |      |      |       
  //          |      |      |      |      |      |      |      |      |       
  //          ---------------------------------------------------------
  //          |      |      |      |      |      |      |      |      |       
  //          |      |      |      |      |      |      |      |      |       
  //          |      |      |      |      |      |      |      |      |       
  //          ---------------------------------------------------------
  //          |      |      |      |      |      |      |      |      |       
  //          |      |      |      |  P1  |  P2  |      |      |      |       
  //          |      |      |      |      |      |      |      |      |       
  //          --------------------------- X ---------------------------
  //          |      |      |      |      |      |      |      |      |       
  //          |      |      |      |  P3  |      |      |      |      |       
  //          |      |      |      |      |      |      |      |      |       
  //          ---------------------------------------------------------
  //          |      |      |      |      |      |      |      |      |       
  //          |      |      |      |      |      |      |      |      |       
  //          |      |      |      |      |      |      |      |      |       
  //          ---------------------------------------------------------
  //          |      |      |      |      |      |      |      |      |       
  //          |      |      |      |      |      |      |      |      |       
  //          |      |      |      |      |      |      |      |      |       
  //          ---------------------------------------------------------
  //          |      |      |      |      |      |      |      |      |       
  //          |      |      |      |      |      |      |      |      |       
  //          |      |      |      |      |      |      |      |      |       
  //          ---------------------------------------------------------
  //
  TopCRTCenters = {{108, { -460.975, 617.388,-1050.61}},
                  {109, { -460.975,617.388,-866.215}},
                  {110, { -460.975,617.388,-681.825}},
                  {111, { -460.975,617.388,-497.435}},
                  {112, { -460.975,617.388,-313.045}},
                  {113, { -460.975,617.388,-128.655}},
                  {114, { -460.975,617.388,55.735}},
                  {115, { -460.975,617.388,240.125}},
                  {116, { -460.975,617.388,424.515}},
                  {117, { -460.975,617.388,608.905}},
                  {118, { -460.975,617.388,793.295}},
                  {119, { -460.975,617.388,977.685}},
                  {120, { -460.975,617.388,1162.07}},
                  {121, { -460.975,617.388,1346.46}},
                  {122, { -276.585,617.388,-1050.61}},
                  {123, { -276.585,617.388,-866.215}},
                  {124, { -276.585,617.388,-681.825}},
                  {125, { -276.585,617.388,-497.435}},
                  {126, { -276.585,617.388,-313.045}},
                  {127, { -276.585,617.388,-128.655}},
                  {128, { -276.585,617.388, 55.735}},
                  {129, { -276.585,617.388, 240.125}},
                  {130, { -276.585,617.388, 424.515}},
                  {131, { -276.585,617.388, 608.905}},
                  {132, { -276.585,617.388, 793.295}},
                  {133, { -276.585,617.388, 977.685}},
                  {134, { -276.585,617.388, 1162.07}},
                  {135, { -276.585,617.388, 1346.46}},
                  {136, { -92.195 ,617.388,-1050.61}},
                  {137, { -92.195 ,617.388,-866.215}},
                  {138, { -92.195 ,617.388,-681.825}},
                  {139, { -92.195 ,617.388,-497.435}},
                  {140, { -92.195 ,617.388,-313.045}},
                  {141, { -92.195 ,617.388,-128.655}},
                  {142, { -92.195 ,617.388, 55.735}},
                  {143, { -92.195 ,617.388, 240.125}},
                  {144, { -92.195 ,617.388, 424.515}},
                  {145, { -92.195 ,617.388, 608.905}},
                  {146, { -92.195 ,617.388, 793.295}},
                  {147, { -92.195 ,617.388, 977.685}},
                  {148, { -92.195 ,617.388, 1162.07}},
                  {149, { -92.195 ,617.388, 1346.46}},
                  {150, { 92.195 , 617.388, -1050.61}},
                  {151, { 92.195 , 617.388, -866.215}},
                  {152, { 92.195 , 617.388, -681.825}},
                  {153, { 92.195 , 617.388, -497.435}},
                  {154, { 92.195 , 617.388, -313.045}},
                  {155, { 92.195 , 617.388, -128.655}},
                  {156, { 92.195 ,617.388, 55.735}},
                  {157, { 92.195 ,617.388, 240.125}},
                  {158, { 92.195 ,617.388, 424.515}},
                  {159, { 92.195 ,617.388, 608.905}},
                  {160, { 92.195 ,617.388, 793.295}},
                  {161, { 92.195 ,617.388, 977.685}},
                  {162, { 92.195 ,617.388, 1162.07}},
                  {163, { 92.195 ,617.388, 1346.46}},
                  {164, { 276.585,617.388, -1050.61}},
                  {165, { 276.585,617.388, -866.215}},
                  {166, { 276.585,617.388, -681.825}},
                  {167, { 276.585,617.388, -497.435}},
                  {168, { 276.585,617.388, -313.045}},
                  {169, { 276.585,617.388, -128.655}},
                  {170, { 276.585,617.388, 55.735}},
                  {171, { 276.585,617.388, 240.125}},
                  {172, { 276.585,617.388, 424.515}},
                  {173, { 276.585,617.388, 608.905}},
                  {174, { 276.585,617.388, 793.295}},
                  {175, { 276.585,617.388, 977.685}},
                  {176, { 276.585,617.388, 1162.07}},
                  {177, { 276.585,617.388, 1346.46}},
                  {178, { 460.975,617.388, -1050.61}},
                  {179, { 460.975,617.388, -866.215}},
                  {180, { 460.975,617.388, -681.825}},
                  {181, { 460.975,617.388, -497.435}},
                  {182, { 460.975,617.388, -313.045}},
                  {183, { 460.975,617.388, -128.655}},
                  {184, { 460.975,617.388, 55.735}},
                  {185, { 460.975,617.388, 240.125}},
                  {186, { 460.975,617.388, 424.515}},
                  {187, { 460.975,617.388, 608.905}},
                  {188, { 460.975,617.388, 793.295}},
                  {189, { 460.975,617.388, 977.685}},
                  {190, { 460.975,617.388, 1162.07}},
                  {191, { 460.975,617.388, 1346.46}},
                  {192, { 555.265, 496.038, -1050.61}},
                  {193, { 555.265, 496.038, -866.215}},
                  {194, { 555.265, 496.038, -681.825}},
                  {195, { 555.265, 496.038, -497.435}},
                  {196, { 555.265, 496.038, -313.045}},
                  {197, { 555.265, 496.038, -128.655}},
                  {198, { 555.265, 496.038, 55.735}},
                  {199, { 555.265, 496.038, 240.125}},
                  {200, { 555.265, 496.038, 424.515}},
                  {201, { 555.265, 496.038, 608.905}},
                  {202, { 555.265, 496.038, 793.295}},
                  {203, { 555.265, 496.038, 977.685}},
                  {204, { 555.265, 496.038, 1162.07}},
                  {205, { 555.265, 496.038, 1346.46}},
                  {206, { -555.265, 496.038, -1050.61}},
                  {207, { -555.265, 496.038, -866.215}},
                  {208, { -555.265, 496.038, -681.825}},
                  {209, { -555.265, 496.038, -497.435}},
                  {210, { -555.265, 496.038, -313.045}},
                  {211, { -555.265, 496.038, -128.655}},
                  {212, { -555.265, 496.038, 55.735}},
                  {213, { -555.265, 496.038, 240.125}},
                  {214, { -555.265, 496.038, 424.515}},
                  {215, { -555.265, 496.038, 608.905}},
                  {216, { -555.265, 496.038, 793.295}},
                  {217, { -555.265, 496.038, 977.685}},
                  {218, { -555.265, 496.038, 1162.07}},
                  {219, { -555.265, 496.038, 1346.46}},
                  {220, { -460.975, 496.038, -1143.4}},
                  {221, { -276.585, 496.038, -1143.4}},
                  {222, { -92.195, 496.038, -1143.4}},
                  {223, { 92.195, 496.038, -1143.4}},
                  {224, { 276.585, 496.038, -1143.4}},
                  // This module does not exist in reality, but exists in simulation
                  {225, { 0, 0, 0}},
                  {226, { -460.975, 525.038, 1533.608}},
                  {227, { -276.585, 525.038, 1533.608}},
                  {228, { -92.195, 525.038, 1533.608}},
                  {229, { 92.195, 525.038, 1533.608}},
                  {230, { 276.585, 525.038, 1533.608}},
                  {231, { 460.975, 525.038, 1533.608}}};
  return TopCRTCenters;
}

TransformedCRTHit AffineTransformation(double DX, double DZ,AffineTrans affine)
{
  double CRTX=affine.B1+DZ*affine.A12+DX*affine.A11;
  double CRTZ=affine.B2+DZ*affine.A22+DX*affine.A21;
  return std::make_pair(CRTX, CRTZ);
}

TopCRTTransformations LoadTopCRTTransformations()
{
  std::string fullFileName;
  cet::search_path searchPath("FW_SEARCH_PATH");
  searchPath.find_file("TopCrtCorrectionByTPC.txt", fullFileName);
  /*if (!searchPath.find_file("TopCrtCorrectionByTPC.txt", fullFileName)) {
      mf::LogError("CRTMatchingUtils_LoadTopCrtTransformation")
      << "Top CRT Correction transformation file not found in FW_SEARCH_PATH.";
  }*/
  std::ifstream corrFile(fullFileName, std::ios::in);
  std::map<int, AffineTrans> Type0Corr;
  std::map<int, AffineTrans> Type1Corr;
  std::map<int, AffineTrans> Type2Corr;
  std::map<int, AffineTrans> Type3Corr;
  std::map<int, AffineTrans> Type4Corr;
  std::map<int, AffineTrans> Type5Corr;
  if (!corrFile.is_open()) {
      mf::LogError("CRTRecoUtils_LoadTopCRTTransformation")
      << "Failed to open Top CRT Correction transformation file: " << fullFileName;
  }
  std::string line;
  while (std::getline(corrFile, line)) {
      std::istringstream iss(line);
      std::vector<double> variables;
      double var;
      while (iss >> var) {
          variables.push_back(var);
      }
      int pos=0;
      AffineTrans Type0 = {variables.at(1+pos),variables.at(2+pos),variables.at(3+pos),variables.at(4+pos),variables.at(5+pos),variables.at(6+pos), variables.at(7+pos),variables.at(8+pos)};
      pos=pos+8;
      AffineTrans Type1 = {variables.at(1+pos),variables.at(2+pos),variables.at(3+pos),variables.at(4+pos),variables.at(5+pos),variables.at(6+pos), variables.at(7+pos),variables.at(8+pos)};
      pos=pos+8;
      AffineTrans Type2 = {variables.at(1+pos),variables.at(2+pos),variables.at(3+pos),variables.at(4+pos),variables.at(5+pos),variables.at(6+pos), variables.at(7+pos),variables.at(8+pos)};
      pos=pos+8;
      AffineTrans Type3 = {variables.at(1+pos),variables.at(2+pos),variables.at(3+pos),variables.at(4+pos),variables.at(5+pos),variables.at(6+pos), variables.at(7+pos),variables.at(8+pos)};
      pos=pos+8;
      AffineTrans Type4 = {variables.at(1+pos),variables.at(2+pos),variables.at(3+pos),variables.at(4+pos),variables.at(5+pos),variables.at(6+pos), variables.at(7+pos),variables.at(8+pos)};
      pos=pos+8;
      AffineTrans Type5 = {variables.at(1+pos),variables.at(2+pos),variables.at(3+pos),variables.at(4+pos),variables.at(5+pos),variables.at(6+pos), variables.at(7+pos),variables.at(8+pos)};
      Type0Corr.insert({variables.at(0), Type0});
      Type1Corr.insert({variables.at(0), Type1});
      Type2Corr.insert({variables.at(0), Type2});
      Type3Corr.insert({variables.at(0), Type3});
      Type4Corr.insert({variables.at(0), Type4});
      Type5Corr.insert({variables.at(0), Type5});
  }
  TopCRTTransformations LoadedTransformations={Type2Corr, Type1Corr, Type0Corr, Type4Corr, Type5Corr, Type3Corr};
  return LoadedTransformations;
}

geo::Point_t ApplyTransformation(const sbn::crt::CRTHit& crthit, const TopCRTCorrectionMap& TopCRTCorrections, const TopCRTCentersMap& TopCRTCenters){
  double centerDX=crthit.x_pos - TopCRTCenters.at((int)crthit.feb_id[0]).X();
  double centerDY=crthit.y_pos - TopCRTCenters.at((int)crthit.feb_id[0]).Y();
  double centerDZ=crthit.z_pos - TopCRTCenters.at((int)crthit.feb_id[0]).Z();
  AffineTrans thisAffine=TopCRTCorrections.at((int)crthit.feb_id[0]);
  std::pair<double,double> transCrt;
  double crtX=crthit.x_pos, crtY=crthit.y_pos, crtZ=crthit.z_pos;
  switch (crthit.plane){
    case 30:
      transCrt=icarus::crt::dataTools::AffineTransformation(centerDX, centerDZ, thisAffine);
      crtX=transCrt.first;
      crtZ=transCrt.second;
      break;
    case 31: case 32:
      transCrt=icarus::crt::dataTools::AffineTransformation(centerDY, centerDZ, thisAffine);
      crtY=transCrt.first;
      crtZ=transCrt.second; 
      break;
    case 33: case 34:
      transCrt=icarus::crt::dataTools::AffineTransformation(centerDX, centerDY, thisAffine);
      crtX=transCrt.first;
      crtY=transCrt.second;
      break;
    default:
      throw std::logic_error("TopCRTAffineTransformationError: CRT plane/region not valid!");
      break;
  }
  return {crtX, crtY, crtZ};
}

}

int RecoUtils::TrueParticleID(detinfo::DetectorClocksData const& clockData,
                              const art::Ptr<recob::Hit> hit, bool rollup_unsaved_ids) {
  std::map<int,double> id_to_energy_map;
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  std::vector<sim::TrackIDE> track_ides = bt_serv->HitToTrackIDEs(clockData, hit);
  for (unsigned int idIt = 0; idIt < track_ides.size(); ++idIt) {
    int id = track_ides.at(idIt).trackID;
    if (rollup_unsaved_ids) id = std::abs(id);
    double energy = track_ides.at(idIt).energy;
    id_to_energy_map[id]+=energy;
  }
  //Now loop over the map to find the maximum contributor
  double likely_particle_contrib_energy = -99999;
  int likely_track_id = 0;
  for (std::map<int,double>::iterator mapIt = id_to_energy_map.begin(); mapIt != id_to_energy_map.end(); mapIt++){
    double particle_contrib_energy = mapIt->second;
    if (particle_contrib_energy > likely_particle_contrib_energy){
      likely_particle_contrib_energy = particle_contrib_energy;
      likely_track_id = mapIt->first;
    }
  }
  return likely_track_id;
}

std::map<int, std::vector<std::pair<geo::WireID, const sim::IDE*>>> RecoUtils::PrepSimChannels(const std::vector<art::Ptr<sim::SimChannel>> &simchannels, const geo::GeometryCore &geom) {
  std::map<int, std::vector<std::pair<geo::WireID, const sim::IDE*>>> ret;

  for (const art::Ptr<sim::SimChannel>& sc : simchannels) {
    // Lookup the wire of this channel
    raw::ChannelID_t channel = sc->Channel();
    std::vector<geo::WireID> maybewire = geom.ChannelToWire(channel);
    geo::WireID thisWire; // Default constructor makes invalid wire
    if (!maybewire.empty()) thisWire = maybewire[0];
    for (const auto &item : sc->TDCIDEMap()) {
      for (const sim::IDE &ide: item.second) {
        // indexing initializes empty vector
        ret[abs(ide.trackID)].push_back({thisWire, &ide});
      }
    }
  }
  return ret;
}

std::map<int, std::vector<art::Ptr<recob::Hit>>> RecoUtils::buildTrackIDtoHitsMap(const std::vector<art::Ptr<recob::Hit>> &allHits, 
  const detinfo::DetectorClocksData &clockData, const cheat::BackTrackerService &backtracker) {
  std::map<int, std::vector<art::Ptr<recob::Hit>>> ret;
  for (const art::Ptr<recob::Hit>& h: allHits) {
    for (int ID: backtracker.HitToTrackIds(clockData, *h)) {
      ret[abs(ID)].push_back(h);
    }
  }
  return ret;
}

int RecoUtils::TrueParticleIDFromTotalTrueEnergy(detinfo::DetectorClocksData const& clockData, const std::vector<art::Ptr<recob::Hit> >& hits, bool rollup_unsaved_ids) {
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  std::map<int,double> trackIDToEDepMap;
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt) {
    art::Ptr<recob::Hit> hit = *hitIt;
    std::vector<sim::TrackIDE> trackIDs = bt_serv->HitToTrackIDEs(clockData, hit);
    for (unsigned int idIt = 0; idIt < trackIDs.size(); ++idIt) {
      int id = trackIDs[idIt].trackID;
      if (rollup_unsaved_ids) id = std::abs(id);
      trackIDToEDepMap[id] += trackIDs[idIt].energy;
    }
  }

  //Loop over the map and find the track which contributes the highest energy to the hit vector
  double maxenergy = -1;
  int objectTrack = -99999;
  for (std::map<int,double>::iterator mapIt = trackIDToEDepMap.begin(); mapIt != trackIDToEDepMap.end(); mapIt++){
    double energy = mapIt->second;
    double trackid = mapIt->first;
    if (energy > maxenergy){
      maxenergy = energy;
      objectTrack = trackid;
    }
  }

  return objectTrack;
}



int RecoUtils::TrueParticleIDFromTotalRecoCharge(detinfo::DetectorClocksData const& clockData, const std::vector<art::Ptr<recob::Hit> >& hits, bool rollup_unsaved_ids) {
  // Make a map of the tracks which are associated with this object and the charge each contributes
  std::map<int,double> trackMap;
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt) {
    art::Ptr<recob::Hit> hit = *hitIt;
    int trackID = TrueParticleID(clockData, hit, rollup_unsaved_ids);
    trackMap[trackID] += hit->Integral();
  }

  // Pick the track with the highest charge as the 'true track'
  double highestCharge = 0;
  int objectTrack = -99999;
  for (std::map<int,double>::iterator trackIt = trackMap.begin(); trackIt != trackMap.end(); ++trackIt) {
    if (trackIt->second > highestCharge) {
      highestCharge = trackIt->second;
      objectTrack  = trackIt->first;
    }
  }
  return objectTrack;
}



int RecoUtils::TrueParticleIDFromTotalRecoHits(detinfo::DetectorClocksData const& clockData,const std::vector<art::Ptr<recob::Hit> >& hits, bool rollup_unsaved_ids) {
  // Make a map of the tracks which are associated with this object and the number of hits they are the primary contributor to
  std::map<int,int> trackMap;
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt) {
    art::Ptr<recob::Hit> hit = *hitIt;
    int trackID = TrueParticleID(clockData, hit, rollup_unsaved_ids);
    trackMap[trackID]++;
  }

  // Pick the track which is the primary contributor to the most hits as the 'true track'
  int objectTrack = -99999;
  int highestCount = -1;
  int NHighestCounts = 0;
  for (std::map<int,int>::iterator trackIt = trackMap.begin(); trackIt != trackMap.end(); ++trackIt) {
    if (trackIt->second > highestCount) {
      highestCount = trackIt->second;
      objectTrack  = trackIt->first;
      NHighestCounts = 1;
    }
    else if (trackIt->second == highestCount){
      NHighestCounts++;
    }
  }
  if (NHighestCounts > 1){
    std::cout<<"RecoUtils::TrueParticleIDFromTotalRecoHits - There are " << NHighestCounts << " particles which tie for highest number of contributing hits (" << highestCount<<" hits).  Using RecoUtils::TrueParticleIDFromTotalTrueEnergy instead."<<std::endl;
    objectTrack = RecoUtils::TrueParticleIDFromTotalTrueEnergy(clockData, hits,rollup_unsaved_ids);
  }
  return objectTrack;
}



bool RecoUtils::IsInsideTPC(TVector3 position, double distance_buffer){
  bool inside = false;
  art::ServiceHandle<geo::Geometry> geom;
  geo::TPCID idtpc = geom->FindTPCAtPosition(geo::vect::toPoint(position));

  if (geom->HasTPC(idtpc))
    {
      const geo::TPCGeo& tpcgeo = geom->GetElement(idtpc);
      double minx = tpcgeo.MinX(); double maxx = tpcgeo.MaxX();
      double miny = tpcgeo.MinY(); double maxy = tpcgeo.MaxY();
      double minz = tpcgeo.MinZ(); double maxz = tpcgeo.MaxZ();

      for (auto const& tpcg : geom->Iterate<geo::TPCGeo>())
	{
	      if (tpcg.MinX() < minx) minx = tpcg.MinX();
	      if (tpcg.MaxX() > maxx) maxx = tpcg.MaxX();
	      if (tpcg.MinY() < miny) miny = tpcg.MinY();
	      if (tpcg.MaxY() > maxy) maxy = tpcg.MaxY();
	      if (tpcg.MinZ() < minz) minz = tpcg.MinZ();
	      if (tpcg.MaxZ() > maxz) maxz = tpcg.MaxZ();
	}

      //x
      double dista = fabs(minx - position.X());
      double distb = fabs(position.X() - maxx);
      if ((position.X() > minx) && (position.X() < maxx) &&
	  (dista > distance_buffer) && (distb > distance_buffer)) inside = true;
      //y
      dista = fabs(maxy - position.Y());
      distb = fabs(position.Y() - miny);
      if (inside && (position.Y() > miny) && (position.Y() < maxy) &&
	  (dista > distance_buffer) && (distb > distance_buffer)) inside = true;
      else inside = false;
      //z
      dista = fabs(maxz - position.Z());
      distb = fabs(position.Z() - minz);
      if (inside && (position.Z() > minz) && (position.Z() < maxz) &&
	  (dista > distance_buffer) && (distb > distance_buffer)) inside = true;
      else inside = false;
    }

  return inside;
}

double RecoUtils::CalculateTrackLength(const art::Ptr<recob::Track> track){
  double length = 0;
  if (track->NumberTrajectoryPoints()==1) return length; //Nothing to calculate if there is only one point

  for (size_t i_tp = 0; i_tp < track->NumberTrajectoryPoints()-1; i_tp++){ //Loop from the first to 2nd to last point
    TVector3 this_point(track->TrajectoryPoint(i_tp).position.X(),track->TrajectoryPoint(i_tp).position.Y(),track->TrajectoryPoint(i_tp).position.Z());
    if (!RecoUtils::IsInsideTPC(this_point,0)){
      std::cout<<"RecoUtils::CalculateTrackLength - Current trajectory point not in the TPC volume.  Skip over this point in the track length calculation"<<std::endl;
      continue;
    }
    TVector3 next_point(track->TrajectoryPoint(i_tp+1).position.X(),track->TrajectoryPoint(i_tp+1).position.Y(),track->TrajectoryPoint(i_tp+1).position.Z());
    if (!RecoUtils::IsInsideTPC(next_point,0)){
      std::cout<<"RecoUtils::CalculateTrackLength - Next trajectory point not in the TPC volume.  Skip over this point in the track length calculation"<<std::endl;
      continue;
    }

    length+=(next_point-this_point).Mag();
  }
  return length;
}