#include "CRTT0MatchAlg.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()

namespace icarus{


  CRTT0MatchAlg::CRTT0MatchAlg(const fhicl::ParameterSet& pset)
  {
    this->reconfigure(pset);
    return;
  }

  CRTT0MatchAlg::CRTT0MatchAlg() = default;
  

  void CRTT0MatchAlg::reconfigure(const fhicl::ParameterSet& pset){

    fMinTrackLength     = pset.get<double>("MinTrackLength", 20.0);
    fTrackDirectionFrac = pset.get<double>("TrackDirectionFrac", 0.5);
    fDistanceLimit      = pset.get<double>("DistanceLimit", 100);
    fTSMode             = pset.get<int>("TSMode", 2);
    fTimeCorrection     = pset.get<double>("TimeCorrection", 0.);
    fSCEposCorr         = pset.get<bool>("SCEposCorr", true);
    fDirMethod          = pset.get<int>("DirMethod", 1);
    fDCAuseBox          = pset.get<bool>("DCAuseBox",false);
    fDCAoverLength      = pset.get<bool>("DCAoverLength", false);
    fDoverLLimit        = pset.get<double>("DoverLLimit", 1);
    fPEcut              = pset.get<double>("PEcut", 0.0);
    fMaxUncert          = pset.get<double>("MaxUncert", 1000.);
    fTPCTrackLabel      = pset.get<std::vector<art::InputTag> >("TPCTrackLabel", {""});
    //  fDistEndpointAVedge = pset.get<double>(.DistEndpointAVedge();

    fGeometryService    = lar::providerFrom<geo::Geometry>();//GeometryService;
    fSCE                = lar::providerFrom<spacecharge::SpaceChargeService>();
    //fSCE = SCE;
    return;

  }
/* 
  matchCand makeNULLmc (){
    sbn::crt::CRTHit hit;
    matchCand null;
    null.thishit = hit;
    null.t0 = -99999;
    null.dca = -99999;
    null.extrapLen = -99999;
    null.best_DCA_pos = -1;
    null.simple_cathodecrosser = false; 
    null.driftdir = -5;
    null.t0min = -99999;
    null.t0max = -99999;
    null.crtTime = -99999;
    null.startDir.SetXYZ(-99999,-99999,-99999);
    null.endDir.SetXYZ(-99999,-99999,-99999);
    null.tpc_track_start.SetXYZ(-99999,-99999,-99999);
    null.tpc_track_end.SetXYZ(-99999,-99999,-99999);

    return null;

  }

   match_geometry makeNULLmc_geo (){
    sbn::crt::CRTHit hit;
    match_geometry null;
    null.thishit = hit;
    null.t0 = -99999;
    null.dca = -99999;
    null.extrapLen = -99999;
    null.best_DCA_pos = -1;
    null.driftdir = -2;

    null.simple_cathodecrosser = false;
    null.t0min = -99999;
    null.t0max = -99999;
    null.crtTime = -99999;

    null.hit_id = -99999;
    null.track_id = -99999;
    null.meta_track_id = -99999;

    null.startDir.SetXYZ(-99999,-99999,-99999);
    null.endDir.SetXYZ(-99999,-99999,-99999);

    null.tpc_track_start.SetXYZ(-99999,-99999,-99999);
    null.tpc_track_end.SetXYZ(-99999,-99999,-99999);

    null.crt_hit_pos.SetXYZ(-99999,-99999,-99999);

    null.simpleDCA_startDir = -99999;
    null.simpleDCA_endDir = -99999;

//    null.is_best_DCA_startDir = false;
//    null.is_best_DCA_endDir = false;

    return null;
  }
*/
  // Utility function that determines the possible t0 range of a track
  std::pair<double, double> CRTT0MatchAlg::TrackT0Range(detinfo::DetectorPropertiesData const& detProp,
							double startX, double endX, int driftDirection, 
							std::pair<double, double> xLimits) const {

    // If track is stitched return zeros
    if(driftDirection == 0) return std::make_pair(0, 0);

    //std::pair<double, double> result; // unused
    double Vd = driftDirection * detProp.DriftVelocity();

    //std::cout << " [ driftdirn, vd ] = [ " << driftDirection << " , " << Vd << " ]" << std::endl;

    // Shift the most postive end to the most positive limit
    double maxX = std::max(startX, endX);
    double maxLimit = std::max(xLimits.first, xLimits.second);
    double maxShift = maxLimit - maxX;
    // Shift the most negative end to the most negative limit
    double minX = std::min(startX, endX);
    double minLimit = std::min(xLimits.first, xLimits.second);
    double minShift = minLimit - minX;
    // Convert to time
    double t0max = maxShift/Vd;
    double t0min = minShift/Vd;

    /*    
    std::cout << "[ driftdirn, vd, startx , endx, xlimits, xlimite, maxx, minx, maxl, minl, maxs, mins, t0max, t0min ] = [ "
	      << driftDirection << " , " << Vd << " , " << startX << " ," <<  endX << " ," << xLimits.first << " ," << xLimits.second 
	      << " ," << maxX <<" ," << minX << " ," <<maxLimit << " ," << minLimit << " ," <<maxShift << " ," <<minShift
	      << " ," << t0max << " ," << t0min << " ]"<< std::endl;
    */
    //  if (t0min>2500)  std::cout << " t0 min " << t0min << " t0max " << t0max << std::endl;
    return std::make_pair(std::min(t0min, t0max), std::max(t0min, t0max));


  } // CRTT0MatchAlg::TrackT0Range()


  double CRTT0MatchAlg::DistOfClosestApproach(detinfo::DetectorPropertiesData const& detProp,
					      geo::Point_t const& track_point, TVector3 trackDir, 
					      sbn::crt::CRTHit const& crtHit, int driftDirection, double t0) const{

    //double minDist = 99999;
    TVector3 trackPos(track_point.X(),track_point.Y(),track_point.Z());

    // Convert the t0 into an x shift
    double xshift = driftDirection* t0 * detProp.DriftVelocity();
    trackPos[0] += xshift;

    if (fSCE->EnableCalSpatialSCE() && fSCEposCorr) {
      geo::Point_t temppt = {trackPos.X(),trackPos.Y(),trackPos.Z()};
      geo::TPCID tpcid = fGeometryService->PositionToTPCID(temppt);
      geo::Vector_t  fPosOffsets = fSCE->GetCalPosOffsets(temppt,tpcid.TPC);
      trackPos[0] += fPosOffsets.X();
      trackPos[1] += fPosOffsets.Y();
      trackPos[2] += fPosOffsets.Z();
    }

    TVector3 end = trackPos + trackDir;
    /*
    std::cout << "[trackPosx, y, z, trackDirx, y, z, endx, y, z ] = [ "
	      << trackPos.X() << " , " << trackPos.Y() << " , " << trackPos.Z() << " , "
	      << trackDir.X() << " , " << trackDir.Y() << " , " << trackDir.Z() << " , " 
	      << end.X()      << " , " << end.Y()      << " , " << end.Z() << " ]" << std::endl;
    */
    //-------- ADDED BY ME----------
    TVector3 pos (crtHit.x_pos, crtHit.y_pos, crtHit.z_pos);
    //double denominator = trackDir.Mag();
    //double numerator = (pos - trackPos).Cross(pos - end).Mag();
    /*
    std::cout << "[ crt x, y, z, startx, directionx, endx, num, denom, distance, altnative, altdist ] = [ "
	      << crtHit.x_pos << " , "  << crtHit.y_pos << " , " << crtHit.z_pos << " , "
	      << trackPos.X() << " , "  << trackDir.X() << " , "  << end.X() << " , "
	      << numerator << " , "  << denominator << " , "  << numerator/denominator << " , "
	      << (pos - trackPos).Cross(trackDir).Mag() << " , "  << ( pos - trackPos).Cross(trackDir).Mag()/denominator << " ] " << std::endl;
    */
    //-----------------------

    // calculate distance of closest approach (DCA)
    //  default is the distance to the point specified by the CRT hit (Simple DCA)
    //    useBox is the distance to the closest edge of the rectangle with the CRT hit at the center and the sides defined
    //   the position uncertainties on the CRT hits.
    double thisdca;

    if (fDCAuseBox) thisdca =   DistToCrtHit(crtHit, trackPos, end);
    else thisdca =  SimpleDCA(crtHit, trackPos, trackDir);
    return thisdca;

  } // CRTT0MatchAlg::DistToOfClosestApproach()


  std::pair<TVector3, TVector3> CRTT0MatchAlg::TrackDirectionAverage(recob::Track const& track, double frac) const
  {
    // Calculate direction as an average over directions
    size_t nTrackPoints = track.NumberTrajectoryPoints();
    recob::TrackTrajectory trajectory  = track.Trajectory();
    std::vector<geo::Vector_t> validDirections;
    for(size_t i = 0; i < nTrackPoints; i++){
      if(trajectory.FlagsAtPoint(i)!=recob::TrajectoryPointFlags::InvalidHitIndex) continue;
      validDirections.push_back(track.DirectionAtPoint(i));
    }

    size_t nValidPoints = validDirections.size();
    int endPoint = (int)floor(nValidPoints*frac);
    double xTotStart = 0; double yTotStart = 0; double zTotStart = 0;
    double xTotEnd = 0; double yTotEnd = 0; double zTotEnd = 0;
    for(int i = 0; i < endPoint; i++){
      geo::Vector_t dirStart = validDirections.at(i);
      geo::Vector_t dirEnd = validDirections.at(nValidPoints - (i+1));
      xTotStart += dirStart.X();
      yTotStart += dirStart.Y();
      zTotStart += dirStart.Z();
      xTotEnd += dirEnd.X();
      yTotEnd += dirEnd.Y();
      zTotEnd += dirEnd.Z();
    }
    TVector3 startDir = {-xTotStart/endPoint, -yTotStart/endPoint, -zTotStart/endPoint};
    TVector3 endDir = {xTotEnd/endPoint, yTotEnd/endPoint, zTotEnd/endPoint};

    return std::make_pair(startDir, endDir);

  } // CRTT0MatchAlg::TrackDirectionAverage()


  std::pair<TVector3, TVector3> CRTT0MatchAlg::TrackDirection(detinfo::DetectorPropertiesData const& detProp,
							      recob::Track const& track, double frac, 
							      double CRTtime, int driftDirection) const{
          
    size_t nTrackPoints = track.NPoints();
    int midPt = (int)floor(nTrackPoints*frac);
    geo::Point_t startP = track.Start();
    geo::Point_t endP = track.End();
    geo::Point_t midP = track.LocationAtPoint(midPt);

    double xshift = driftDirection * CRTtime * detProp.DriftVelocity();
    TVector3  startPoint = {startP.X()+xshift,startP.Y(),startP.Z()};
    TVector3  endPoint = {endP.X()+xshift,endP.Y(),endP.Z()};
    TVector3  midPoint = {midP.X()+xshift,midP.Y(),midP.Z()};

    //    std::cout <<"[ nTrackPoints, midPt, startP, endP, midP, xshift, CRTtime, startPoint, endPoint,  midPoint ] = [ " 
    //	      << nTrackPoints << " , "  << midPt << " , "  << startP.X() << " , "  << endP.X() << " , "  << midP.X() << " , "  
    //	      << xshift << " , "  << CRTtime<< " , "  << startPoint.X() << " , "  << endPoint.X() << " , "  << midPoint.X() << " ]"<<std::endl;

    if (fSCE->EnableCalSpatialSCE() && fSCEposCorr) {

      // Apply the shift depending on which TPC the track is in                                 
      geo::Point_t fTrackPos = startP;
      //std::cout <<" before set fTrackPos " << fTrackPos.X() <<std::endl;
      fTrackPos.SetX(startPoint.X());

      geo::TPCID tpcid = fGeometryService->PositionToTPCID(fTrackPos);                        
      geo::Vector_t fPosOffsets = fSCE->GetCalPosOffsets(geo::Point_t{fTrackPos.X(),fTrackPos.Y(),fTrackPos.Z()},tpcid.TPC);

      startPoint.SetX(fTrackPos.X() + fPosOffsets.X());                                       
      startPoint.SetY(fTrackPos.Y() + fPosOffsets.Y());                                       
      startPoint.SetZ(fTrackPos.Z() + fPosOffsets.Z());                                       
      // std::cout <<" [ after set fTrackPos, tpcid, offset x, offset y, offset z, startx, starty, startz ] = [ " << fTrackPos.X()  << " , "  << tpcid.TPC
      //	<< " , "  << fPosOffsets.X() << " , "  << fPosOffsets.Y() << " , "  << fPosOffsets.Z()
      //	<< " , "  << startPoint.X()<< " , "  << startPoint.Y() << " , "  << startPoint.Z() << " ]" <<std::endl;
      fTrackPos = endP;
      fTrackPos.SetX(endPoint.X());
      tpcid = fGeometryService->PositionToTPCID(fTrackPos);
      //      fPosOffsets = fSCE->GetCalPosOffsets(fTrackPos,tpcid.TPC);
      fPosOffsets = fSCE->GetCalPosOffsets(geo::Point_t{fTrackPos.X(),fTrackPos.Y(),fTrackPos.Z()},tpcid.TPC);
      endPoint.SetX(fTrackPos.X() + fPosOffsets.X());
      endPoint.SetY(fTrackPos.Y() + fPosOffsets.Y());
      endPoint.SetZ(fTrackPos.Z() + fPosOffsets.Z());
      //std::cout <<" [ after set end fTrackPos, tpcid, offset x, offset y, offset z, startx, starty, startz ] = [ " << fTrackPos.X() << " , "  << tpcid.TPC
      //	<< " , "  << fPosOffsets.X() << " , "  << fPosOffsets.Y() << " , "  << fPosOffsets.Z()
      //	<< " , "  << endPoint.X()<< " , "  << endPoint.Y() << " , "  << endPoint.Z() << " ]" <<std::endl;

      fTrackPos = midP;
      fTrackPos.SetX(midPoint.X());
      tpcid = fGeometryService->PositionToTPCID(fTrackPos);
      //fPosOffsets = fSCE->GetCalPosOffsets(fTrackPos,tpcid.TPC);
      fPosOffsets = fSCE->GetCalPosOffsets(geo::Point_t{fTrackPos.X(),fTrackPos.Y(),fTrackPos.Z()},tpcid.TPC);
      midPoint.SetX(fTrackPos.X() + fPosOffsets.X());
      midPoint.SetY(fTrackPos.Y() + fPosOffsets.Y());
      midPoint.SetZ(fTrackPos.Z() + fPosOffsets.Z());
      //std::cout <<" [ after set mid fTrackPos, tpcid, offset x, offset y, offset z, startx, starty, startz ] = [ " << fTrackPos.X()  << " , "  << tpcid.TPC
      //        << " , "  << fPosOffsets.X() << " , "  << fPosOffsets.Y() << " , "  << fPosOffsets.Z()
      //	<< " , "  << midPoint.X()<< " , "  << midPoint.Y() << " , "  << midPoint.Z() << " ]" <<std::endl;
    }
    
    TVector3 startDir = {midPoint.X()-startPoint.X(),midPoint.Y()-startPoint.Y(),midPoint.Z()-startPoint.Z()};
    float norm = startDir.Mag();
    if (norm>0)  startDir *=(1.0/norm);
    /*
    std::cout <<" [ startDirx, startDiry, startDirz, mag, xcap, ycap, zcap ] = [ " 
	      << midPoint.X()-startPoint.X() << " , "  << midPoint.Y()-startPoint.Y() << " , "  << midPoint.Z()-startPoint.Z()
	      << " , "  << norm << " , "  << startDir.X()<< " , "  << startDir.Y() << " , "  << startDir.Z() << " ]" <<std::endl;
    */
    TVector3 endDir = {midPoint.X()-endPoint.X(),midPoint.Y()-endPoint.Y(),midPoint.Z()-endPoint.Z()};    
    norm = endDir.Mag();
    if (norm>0)  endDir *=(1.0/norm);
    /*
    std::cout <<" [ endDirx, endDiry, endDirz, mag, xcap, ycap, zcap ] = [ "
	      << midPoint.X()-endPoint.X() << " , "  << midPoint.Y()-endPoint.Y() << " , "  << midPoint.Z()-endPoint.Z()
              << " , "  << norm<< " , "  << endDir.X()<< " , "  << endDir.Y() << " , "  << endDir.Z() << " ]" <<std::endl;
    */
    return std::make_pair(startDir, endDir);
    
  } // CRTT0MatchAlg::TrackDirection()                                                                  

  std::pair<TVector3, TVector3> CRTT0MatchAlg::TrackDirectionAverageFromPoints(recob::Track const& track, double frac) const{

    // Calculate direction as an average over directions
    size_t nTrackPoints = track.NumberTrajectoryPoints();
    recob::TrackTrajectory trajectory  = track.Trajectory();
    std::vector<TVector3> validPoints;
    for(size_t i = 0; i < nTrackPoints; i++){
      if(trajectory.FlagsAtPoint(i) != recob::TrajectoryPointFlags::InvalidHitIndex) continue;
      validPoints.push_back(track.LocationAtPoint<TVector3>(i));
    }

    size_t nValidPoints = validPoints.size();
    int endPoint = (int)floor(nValidPoints*frac);
    TVector3 startDir = validPoints.at(0) - validPoints.at(endPoint-1);
    TVector3 endDir = validPoints.at(nValidPoints - 1) - validPoints.at(nValidPoints - (endPoint));

    return std::make_pair(startDir.Unit(), endDir.Unit());

  } // CRTT0MatchAlg::TrackDirectionAverageFromPoints()


  // Keeping ClosestCRTHit function for backward compatibility only
  // *** use GetClosestCRTHit instead

  std::vector<std::pair<sbn::crt::CRTHit, double> > CRTT0MatchAlg::ClosestCRTHit(detinfo::DetectorPropertiesData const& detProp,
										 recob::Track const& tpcTrack, std::vector<sbn::crt::CRTHit> const& crtHits, 
										 const art::Event& event, uint64_t trigger_timestamp) const{
    //    matchCand newmc = makeNULLmc();
    std::vector<std::pair<sbn::crt::CRTHit, double> > crthitpair;
    
    for(const auto& trackLabel : fTPCTrackLabel){
      auto tpcTrackHandle = event.getValidHandle<std::vector<recob::Track>>(trackLabel);
      if (!tpcTrackHandle.isValid()) continue;
      
      art::FindManyP<recob::Hit> findManyHits(tpcTrackHandle, event, trackLabel);
      for (auto const& tpcTrack : (*tpcTrackHandle)){
	std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(tpcTrack.ID());
	
	crthitpair.push_back(ClosestCRTHit(detProp, tpcTrack, hits, crtHits, trigger_timestamp));
	//	return ClosestCRTHit(detProp, tpcTrack, hits, crtHits);
      }
    }

    return crthitpair;
    //for(const auto& crthit : crthitpair)
    //return std::make_pair(crthit.first, crthit.second);

    //return std::make_pair( newmc.thishit, -9999);
  }


  std::pair<sbn::crt::CRTHit, double>  CRTT0MatchAlg::ClosestCRTHit(detinfo::DetectorPropertiesData const& detProp,
								    recob::Track const& tpcTrack, std::vector<art::Ptr<recob::Hit>> const& hits, 
								    std::vector<sbn::crt::CRTHit> const& crtHits, uint64_t trigger_timestamp) const{

    auto start = tpcTrack.Vertex();
    auto end = tpcTrack.End();
    // Get the drift direction from the TPC
    int driftDirection = TPCGeoUtil::DriftDirectionFromHits(fGeometryService, hits);
    std::pair<double, double> xLimits = TPCGeoUtil::XLimitsFromHits(fGeometryService, hits);
    // Get the allowed t0 range
    std::pair<double, double> t0MinMax = TrackT0Range(detProp, start.X(), end.X(), driftDirection, xLimits);

    return ClosestCRTHit(detProp, tpcTrack, t0MinMax, crtHits, driftDirection, trigger_timestamp);
  }

  std::pair<sbn::crt::CRTHit, double> CRTT0MatchAlg::ClosestCRTHit(detinfo::DetectorPropertiesData const& detProp,
								   recob::Track const& tpcTrack, std::pair<double, double> t0MinMax, 
								   std::vector<sbn::crt::CRTHit> const& crtHits, int driftDirection, uint64_t trigger_timestamp) const{

    matchCand bestmatch = GetClosestCRTHit(detProp, tpcTrack,t0MinMax,crtHits,driftDirection, trigger_timestamp, false);
    return std::make_pair(bestmatch.thishit,bestmatch.dca);

  }


  matchCand CRTT0MatchAlg::GetClosestCRTHit(detinfo::DetectorPropertiesData const& detProp,
					    recob::Track const& tpcTrack, std::vector<art::Ptr<recob::Hit>> const& hits, 
					    std::vector<sbn::crt::CRTHit> const& crtHits, uint64_t trigger_timestamp, bool IsData) const{

    auto start = tpcTrack.Vertex();
    auto end   = tpcTrack.End();



    // Get the drift direction from the TPC
    int driftDirection = TPCGeoUtil::DriftDirectionFromHits(fGeometryService, hits);
    //std::cout << "size of hit in a track: " << hits.size() << ", driftDirection: "<< driftDirection 
    //	      << " , tpc: "<< hits[0]->WireID().TPC << std::endl; //<< " , intpc: "<< icarus::TPCGeoUtil::DetectedInTPC(hits) << std::endl;
    std::pair<double, double> xLimits = TPCGeoUtil::XLimitsFromHits(fGeometryService, hits);
    // Get the allowed t0 range
    std::pair<double, double> t0MinMax = TrackT0Range(detProp, start.X(), end.X(), driftDirection, xLimits);

    return GetClosestCRTHit(detProp, tpcTrack, t0MinMax, crtHits, driftDirection, trigger_timestamp, IsData);

  }

  std::vector<matchCand> CRTT0MatchAlg::GetClosestCRTHit(detinfo::DetectorPropertiesData const& detProp,
							 recob::Track const& tpcTrack, std::vector<sbn::crt::CRTHit> const& crtHits, 
							 const art::Event& event, uint64_t trigger_timestamp) const{
    //    matchCand nullmatch = makeNULLmc();
    std::vector<matchCand> matchcanvec;
    //std::vector<std::pair<sbn::crt::CRTHit, double> > matchedCan;
    for(const auto& trackLabel : fTPCTrackLabel){
      auto tpcTrackHandle = event.getValidHandle<std::vector<recob::Track>>(trackLabel);
      if (!tpcTrackHandle.isValid()) continue;

      art::FindManyP<recob::Hit> findManyHits(tpcTrackHandle, event, trackLabel);
      for (auto const& tpcTrack : (*tpcTrackHandle)){
	std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(tpcTrack.ID());
        matchcanvec.push_back(GetClosestCRTHit(detProp, tpcTrack, hits, crtHits, trigger_timestamp, false));
	//return ClosestCRTHit(detProp, tpcTrack, hits, crtHits);
	//matchCand closestHit = GetClosestCRTHit(detProp, tpcTrack, hits, crtHits);

      }
    }
    return matchcanvec;
    //auto tpcTrackHandle = event.getValidHandle<std::vector<recob::Track>>(fTPCTrackLabel);
    //art::FindManyP<recob::Hit> findManyHits(tpcTrackHandle, event, fTPCTrackLabel);
    //std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(tpcTrack.ID());
    //return GetClosestCRTHit(detProp, tpcTrack, hits, crtHits);
    //    for (const auto& match : matchedCan)
    //return match;
    //return nullmatch;
  }


  matchCand CRTT0MatchAlg::GetClosestCRTHit(detinfo::DetectorPropertiesData const& detProp,
					    recob::Track const& tpcTrack, std::pair<double, double> t0MinMax, 
					    std::vector<sbn::crt::CRTHit> const& crtHits, int driftDirection, uint64_t& trigger_timestamp, bool IsData) const {

    auto start = tpcTrack.Vertex();
    auto end   = tpcTrack.End();

    bool simple_cathode_crosscheck =( (std::abs(start.X()) < 210.215) != (std::abs(end.X()) < 210.215));
    // ====================== Matching Algorithm ========================== //
    //  std::vector<std::pair<sbn::crt::CRTHit, double>> t0Candidates;
    std::vector<matchCand> t0Candidates;

    //    if (crtHits.size() == 0) continue;
    // Loop over all the CRT hits
    for(auto &crtHit : crtHits){
      // Check if hit is within the allowed t0 range
      double crtTime = GetCRTTime(crtHit,trigger_timestamp,IsData);  // units are us
/*      
     if(IsData){
      if (fTSMode == 1) {
	crtTime = ((double)(int)crtHit.ts1_ns) * 1e-3; //+ fTimeCorrection;
      }//end if (fTSMode == 1)
      else {
		crtTime = double(crtHit.ts0_ns - (trigger_timestamp%1'000'000'000))/1e3;
		if(crtTime<-0.5e6) 	crtTime+=1e6;
		else if(crtTime>=0.5e6) 	crtTime-=1e6;		
      }//end else
     }//end if(IsData)
     else if(!IsData){	
	crtTime = crtHit.ts0_ns/1e3;
     }//end else if(!IsData)*/
      // If track is stitched then try all hits
      if (!((crtTime >= t0MinMax.first - 10. && crtTime <= t0MinMax.second + 10.) 
            || t0MinMax.first == t0MinMax.second)) continue;

      // cut on CRT hit PE value
      if (crtHit.peshit<fPEcut) continue;
      if (crtHit.x_err>fMaxUncert) continue;
      if (crtHit.y_err>fMaxUncert) continue;
      if (crtHit.z_err>fMaxUncert) continue;
      if (tpcTrack.Length() < fMinTrackLength) continue;

      geo::Point_t crtPoint(crtHit.x_pos, crtHit.y_pos, crtHit.z_pos);

      //Calculate Track direction
      std::pair<TVector3, TVector3> startEndDir;
      // dirmethod=2 is original algorithm, dirmethod=1 is simple algorithm for which SCE corrections are possible
      if (fDirMethod==2)  startEndDir = TrackDirectionAverage(tpcTrack, fTrackDirectionFrac);
      else startEndDir = TrackDirection(detProp, tpcTrack, fTrackDirectionFrac, crtTime, driftDirection);
      TVector3 startDir = startEndDir.first;
      TVector3 endDir = startEndDir.second;
    
      // Calculate the distance between the crossing point and the CRT hit, SCE corrections are done inside but dropped
      double startDist = DistOfClosestApproach(detProp, start, startDir, crtHit, driftDirection, crtTime);
      double endDist = DistOfClosestApproach(detProp, end, endDir, crtHit, driftDirection, crtTime);

    
      double xshift = driftDirection * crtTime * detProp.DriftVelocity();
      auto thisstart = start; 
      thisstart.SetX(start.X()+xshift);
      auto thisend = end; 
      thisend.SetX(end.X()+xshift);
      // repeat SCE correction for endpoints
      if (fSCE->EnableCalSpatialSCE() && fSCEposCorr) {
	geo::TPCID tpcid = fGeometryService->PositionToTPCID(thisstart);
	thisstart+= fSCE->GetCalPosOffsets(thisstart,tpcid.TPC);
	tpcid = fGeometryService->PositionToTPCID(thisend);
	thisend+= fSCE->GetCalPosOffsets(thisend,tpcid.TPC);

      }

      matchCand newmc;// = makeNULLmc();
      if (startDist<fDistanceLimit || endDist<fDistanceLimit) {
	double distS = (crtPoint-thisstart).R();
	double distE =  (crtPoint-thisend).R();
	if (distS <= distE && startDist<fDistanceLimit){ 
	  newmc.dca = startDist;
	  newmc.extrapLen = distS;
	  newmc.best_DCA_pos=0;
	}//end if(distS < distE)
	else if(distE<=distS && endDist<fDistanceLimit ){
	  newmc.dca = endDist;
	  newmc.extrapLen = distE;
	  newmc.best_DCA_pos=1;
	}//end else if(distE<=distS && endDist<fDistanceLimit )
	else continue;
	newmc.thishit = crtHit;
	newmc.t0= crtTime;
	newmc.simple_cathodecrosser = simple_cathode_crosscheck;
	newmc.driftdir = driftDirection;
	newmc.t0min = t0MinMax.first;
	newmc.t0max = t0MinMax.second;
	newmc.crtTime = crtTime;
	newmc.startDir = startDir;
	newmc.endDir = endDir;
	newmc.tpc_track_start.SetXYZ(thisstart.X(),thisstart.Y(),thisstart.Z());
	newmc.tpc_track_end.SetXYZ(thisend.X(),thisend.Y(),thisend.Z());
	t0Candidates.push_back(newmc);

      }//end if (startDist<fDistanceLimit || endDist<fDistanceLimit) 
    }//end loop over CRT Hits


      //std::cout << " found " << t0Candidates.size() << " candidates" << std::endl;
    matchCand bestmatch;// = makeNULLmc();
    if(t0Candidates.size() > 0){
      // Find candidate with shortest DCA or DCA/L value
      bestmatch=t0Candidates[0];
      double sin_angle = bestmatch.dca/bestmatch.extrapLen;
      if (fDCAoverLength) { // Use dca/extrapLen to judge best
	for(auto &thisCand : t0Candidates){
	  double this_sin_angle = thisCand.dca/thisCand.extrapLen;
	  if (bestmatch.dca<0 )bestmatch=thisCand;
	  else if (this_sin_angle<sin_angle && thisCand.dca>=0)bestmatch=thisCand;
	}//end for(auto &thisCand : t0Candidates)
      }//end if (fDCAoverLength)
      else { // use Dca to judge best
	for(auto &thisCand : t0Candidates){
	  if (bestmatch.dca<0 )bestmatch=thisCand;
	  else if (thisCand.dca<bestmatch.dca && thisCand.dca>=0)bestmatch=thisCand;
	}//end for(auto &thisCand : t0Candidates)
      }//end else [use DCA for best match method]
    }//end if(t0Candidates.size() > 0)

    //std::cout << "best match has dca of " << bestmatch.dca << std::endl;
    return bestmatch;

  }//end function defn


  std::vector<double> CRTT0MatchAlg::T0FromCRTHits(detinfo::DetectorPropertiesData const& detProp,
						   recob::Track const& tpcTrack, std::vector<sbn::crt::CRTHit> const& crtHits, 
						   const art::Event& event, uint64_t trigger_timestamp) const{
    std::vector<double> ftime;
    for(const auto& trackLabel : fTPCTrackLabel){
      auto tpcTrackHandle = event.getValidHandle<std::vector<recob::Track>>(trackLabel);
      if (!tpcTrackHandle.isValid()) continue;

      art::FindManyP<recob::Hit> findManyHits(tpcTrackHandle, event, trackLabel);
      for (auto const& tpcTrack : (*tpcTrackHandle)){
	std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(tpcTrack.ID());
	ftime.push_back(T0FromCRTHits(detProp, tpcTrack, hits, crtHits, trigger_timestamp));
	// return T0FromCRTHits(detProp, tpcTrack, hits, crtHits);
      }
    }
    return ftime;
    //auto tpcTrackHandle = event.getValidHandle<std::vector<recob::Track>>(fTPCTrackLabel);
    //art::FindManyP<recob::Hit> findManyHits(tpcTrackHandle, event, fTPCTrackLabel);
    //std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(tpcTrack.ID());
    //return T0FromCRTHits(detProp, tpcTrack, hits, crtHits);

    //    for(const auto& t0 : ftime) return t0;
    //return -99999;
  }

  double CRTT0MatchAlg::T0FromCRTHits(detinfo::DetectorPropertiesData const& detProp,
				      recob::Track const& tpcTrack, std::vector<art::Ptr<recob::Hit>> const& hits, 
				      std::vector<sbn::crt::CRTHit> const& crtHits, uint64_t& trigger_timestamp)  const{

    if (tpcTrack.Length() < fMinTrackLength) return -99999; 

    matchCand closestHit = GetClosestCRTHit(detProp, tpcTrack, hits, crtHits, trigger_timestamp, false);
    if(closestHit.dca <0) return -99999;

    double crtTime;
    if (fTSMode == 1) {
      crtTime = ((double)(int)closestHit.thishit.ts1_ns) * 1e-3; //+ fTimeCorrection;
    }
    else {
      crtTime = ((double)(int)closestHit.thishit.ts0_ns) * 1e-3 + fTimeCorrection;
    }
    if (closestHit.dca < fDistanceLimit && (closestHit.dca/closestHit.extrapLen) < fDoverLLimit) return crtTime;

    return -99999;

  }

  std::vector<std::pair<double, double> > CRTT0MatchAlg::T0AndDCAFromCRTHits(detinfo::DetectorPropertiesData const& detProp,
									     recob::Track const& tpcTrack, std::vector<sbn::crt::CRTHit> const& crtHits, 
									     const art::Event& event, uint64_t trigger_timestamp) const{ 
   
    std::vector<std::pair<double, double> > ft0anddca;
    for(const auto& trackLabel : fTPCTrackLabel){
      auto tpcTrackHandle = event.getValidHandle<std::vector<recob::Track>>(trackLabel);
      if (!tpcTrackHandle.isValid()) continue;

      art::FindManyP<recob::Hit> findManyHits(tpcTrackHandle, event, trackLabel);
      for (auto const& tpcTrack : (*tpcTrackHandle)){
	std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(tpcTrack.ID());
	ft0anddca.push_back(T0AndDCAFromCRTHits(detProp, tpcTrack, hits, crtHits, trigger_timestamp));
	//        return T0AndDCAFromCRTHits(detProp, tpcTrack, hits, crtHits);
      }
    }
    return ft0anddca;
    // auto tpcTrackHandle = event.getValidHandle<std::vector<recob::Track>>(fTPCTrackLabel);
    //art::FindManyP<recob::Hit> findManyHits(tpcTrackHandle, event, fTPCTrackLabel);
    //std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(tpcTrack.ID());
    //return T0AndDCAFromCRTHits(detProp, tpcTrack, hits, crtHits);
    //    for(const auto& t0anddca : ft0anddca) return std::make_pair(t0anddca.first, t0anddca.second);
    //return  std::make_pair(-9999., -9999.);
  }

  std::pair<double, double> CRTT0MatchAlg::T0AndDCAFromCRTHits(detinfo::DetectorPropertiesData const& detProp,
							       recob::Track const& tpcTrack, std::vector<art::Ptr<recob::Hit>> const& hits, 
							       std::vector<sbn::crt::CRTHit> const& crtHits, uint64_t& trigger_timestamp) const{

    if (tpcTrack.Length() < fMinTrackLength) return std::make_pair(-9999., -9999.);

    matchCand closestHit = GetClosestCRTHit(detProp, tpcTrack, hits, crtHits, trigger_timestamp, false);

    if(closestHit.dca < 0 ) return std::make_pair(-9999., -9999.);
    if (closestHit.dca < fDistanceLimit && (closestHit.dca/closestHit.extrapLen) < fDoverLLimit) return std::make_pair(closestHit.t0, closestHit.dca);

    return std::make_pair(-9999., -9999.);


  }

  // Simple distance of closest approach between infinite track and centre of hit
  double CRTT0MatchAlg::SimpleDCA(sbn::crt::CRTHit const& hit, TVector3 start, TVector3 direction) const{

    TVector3 pos (hit.x_pos, hit.y_pos, hit.z_pos);
    TVector3 end = start + direction;
    double denominator = direction.Mag();
    double numerator = (pos - start).Cross(pos - end).Mag();
    return numerator/denominator;

  }

  // Minimum distance from infinite track to CRT hit assuming that hit is a 2D square
  double CRTT0MatchAlg::DistToCrtHit(sbn::crt::CRTHit const& hit, TVector3 start, TVector3 end) const{

    // Check if track goes inside hit
    TVector3 min (hit.x_pos - hit.x_err, hit.y_pos - hit.y_err, hit.z_pos - hit.z_err);
    TVector3 max (hit.x_pos + hit.x_err, hit.y_pos + hit.y_err, hit.z_pos + hit.z_err);
    if(CubeIntersection(min, max, start, end).first.X() != -99999) return 0;

    // Calculate the closest distance to each edge of the CRT hit
    // Assume min error is the fixed position of tagger
    TVector3 vertex1 (hit.x_pos, hit.y_pos - hit.y_err, hit.z_pos - hit.z_err);
    TVector3 vertex2 (hit.x_pos, hit.y_pos + hit.y_err, hit.z_pos - hit.z_err);
    TVector3 vertex3 (hit.x_pos, hit.y_pos - hit.y_err, hit.z_pos + hit.z_err);
    TVector3 vertex4 (hit.x_pos, hit.y_pos + hit.y_err, hit.z_pos + hit.z_err);
    if(hit.y_err < hit.x_err && hit.y_err < hit.z_err){
      vertex1.SetXYZ(hit.x_pos - hit.x_err, hit.y_pos, hit.z_pos - hit.z_err);
      vertex2.SetXYZ(hit.x_pos + hit.x_err, hit.y_pos, hit.z_pos - hit.z_err);
      vertex3.SetXYZ(hit.x_pos - hit.x_err, hit.y_pos, hit.z_pos + hit.z_err);
      vertex4.SetXYZ(hit.x_pos + hit.x_err, hit.y_pos, hit.z_pos + hit.z_err);
    }
    if(hit.z_err < hit.x_err && hit.z_err < hit.y_err){
      vertex1.SetXYZ(hit.x_pos - hit.x_err, hit.y_pos - hit.y_err, hit.z_pos);
      vertex2.SetXYZ(hit.x_pos + hit.x_err, hit.y_pos - hit.y_err, hit.z_pos);
      vertex3.SetXYZ(hit.x_pos - hit.x_err, hit.y_pos + hit.y_err, hit.z_pos);
      vertex4.SetXYZ(hit.x_pos + hit.x_err, hit.y_pos + hit.y_err, hit.z_pos);
    }

    double dist1 = LineSegmentDistance(vertex1, vertex2, start, end);
    double dist2 = LineSegmentDistance(vertex1, vertex3, start, end);
    double dist3 = LineSegmentDistance(vertex4, vertex2, start, end);
    double dist4 = LineSegmentDistance(vertex4, vertex3, start, end);

    return std::min(std::min(dist1, dist2), std::min(dist3, dist4));

  }


  // Distance between infinite line (2) and segment (1)
  // http://geomalgorithms.com/a07-_distance.html
  double CRTT0MatchAlg::LineSegmentDistance(TVector3 start1, TVector3 end1, TVector3 start2, TVector3 end2) const{

    double smallNum = 0.00001;

    // 1 is segment
    TVector3 direction1 = end1 - start1;
    // 2 is infinite line
    TVector3 direction2 = end2 - start2;

    TVector3 u = direction1;
    TVector3 v = direction2;
    TVector3 w = start1 - start2;

    double a = u.Dot(u);
    double b = u.Dot(v);
    double c = v.Dot(v);
    double d = u.Dot(w);
    double e = v.Dot(w);
    double D = a * c - b * b;
    double sc, sN, sD = D; // sc = sN/sD
    double tc, tN, tD = D; // sc = sN/sD

    // Compute the line parameters of the two closest points
    if(D < smallNum){ // Lines are almost parallel
      sN = 0.0;
      sD = 1.0;
      tN = e;
      tD = c;
    }
    else{
      sN = (b * e - c * d)/D;
      tN = (a * e - b * d)/D;
      if(sN < 0.){ // sc < 0, the s = 0 edge is visible
	sN = 0.;
	tN = e;
	tD = c;
      }
      else if(sN > sD){ // sc > 1, the s = 1 edge is visible
	sN = sD;
	tN = e + b;
	tD = c;
      } 
    }

    sc = (std::abs(sN) < smallNum ? 0.0 : sN / sD);
    tc = (std::abs(tN) < smallNum ? 0.0 : tN / tD);
    // Get the difference of the two closest points
    TVector3 dP = w + (sc * u) - (tc * v);

    return dP.Mag();

  }

  // Intersection between axis-aligned cube and infinite line
  // (https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-box-intersection)
  std::pair<TVector3, TVector3> CRTT0MatchAlg::CubeIntersection(TVector3 min, TVector3 max, TVector3 start, TVector3 end) const{

    TVector3 dir = (end - start);
    TVector3 invDir (1./dir.X(), 1./dir.Y(), 1/dir.Z());

    double tmin, tmax, tymin, tymax, tzmin, tzmax;

    TVector3 enter (-99999, -99999, -99999);
    TVector3 exit (-99999, -99999, -99999);

    // Find the intersections with the X plane
    if(invDir.X() >= 0){
      tmin = (min.X() - start.X()) * invDir.X();
      tmax = (max.X() - start.X()) * invDir.X();
    }
    else{
      tmin = (max.X() - start.X()) * invDir.X();
      tmax = (min.X() - start.X()) * invDir.X();
    }

    // Find the intersections with the Y plane
    if(invDir.Y() >= 0){
      tymin = (min.Y() - start.Y()) * invDir.Y();
      tymax = (max.Y() - start.Y()) * invDir.Y();
    }
    else{
      tymin = (max.Y() - start.Y()) * invDir.Y();
      tymax = (min.Y() - start.Y()) * invDir.Y();
    }

    // Check that it actually intersects
    if((tmin > tymax) || (tymin > tmax)) return std::make_pair(enter, exit);

    // Max of the min points is the actual intersection
    if(tymin > tmin) tmin = tymin;

    // Min of the max points is the actual intersection
    if(tymax < tmax) tmax = tymax;

    // Find the intersection with the Z plane
    if(invDir.Z() >= 0){
      tzmin = (min.Z() - start.Z()) * invDir.Z();
      tzmax = (max.Z() - start.Z()) * invDir.Z();
    }
    else{
      tzmin = (max.Z() - start.Z()) * invDir.Z();
      tzmax = (min.Z() - start.Z()) * invDir.Z();
    }

    // Check for intersection
    if((tmin > tzmax) || (tzmin > tmax)) return std::make_pair(enter, exit);

    // Find final intersection points
    if(tzmin > tmin) tmin = tzmin;

    // Find final intersection points
    if(tzmax < tmax) tmax = tzmax;

    // Calculate the actual crossing points
    double xmin = start.X() + tmin * dir.X();
    double xmax = start.X() + tmax * dir.X();
    double ymin = start.Y() + tmin * dir.Y();
    double ymax = start.Y() + tmax * dir.Y();
    double zmin = start.Z() + tmin * dir.Z();
    double zmax = start.Z() + tmax * dir.Z();

    // Return pair of entry and exit points
    enter.SetXYZ(xmin, ymin, zmin);
    exit.SetXYZ(xmax, ymax, zmax);
    return std::make_pair(enter, exit);

  }
  std::vector<icarus::match_geometry> CRTT0MatchAlg::GetClosestCRTHit_geo(detinfo::DetectorPropertiesData const& detProp,
					    recob::Track const& tpcTrack, std::vector<art::Ptr<recob::Hit>> const& hits, 
					    std::vector<sbn::crt::CRTHit> const& crtHits, uint64_t trigger_timestamp, bool IsData) const{

    auto start = tpcTrack.Vertex();
    auto end   = tpcTrack.End();



    // Get the drift direction from the TPC
    int driftDirection = TPCGeoUtil::DriftDirectionFromHits(fGeometryService, hits);
    //std::cout << "size of hit in a track: " << hits.size() << ", driftDirection: "<< driftDirection 
    //	      << " , tpc: "<< hits[0]->WireID().TPC << std::endl; //<< " , intpc: "<< icarus::TPCGeoUtil::DetectedInTPC(hits) << std::endl;
    std::pair<double, double> xLimits = TPCGeoUtil::XLimitsFromHits(fGeometryService, hits);
    // Get the allowed t0 range
    std::pair<double, double> t0MinMax = TrackT0Range(detProp, start.X(), end.X(), driftDirection, xLimits);

    return GetClosestCRTHit_geo(detProp, tpcTrack, t0MinMax, crtHits, driftDirection, trigger_timestamp, IsData);

  }
/*
  std::vector<icarus::match_geometry> CRTT0MatchAlg::GetClosestCRTHit_geo(detinfo::DetectorPropertiesData const& detProp,
							 recob::Track tpcTrack, std::vector<sbn::crt::CRTHit> crtHits, 
							 const art::Event& event, uint64_t trigger_timestamp, bool IsData) {
    //    matchCand nullmatch = makeNULLmc();
    std::vector<icarus::match_geometry> matchcanvec;
    //std::vector<std::pair<sbn::crt::CRTHit, double> > matchedCan;
    for(const auto& trackLabel : fTPCTrackLabel){
      auto tpcTrackHandle = event.getValidHandle<std::vector<recob::Track>>(trackLabel);
      if (!tpcTrackHandle.isValid()) continue;

      art::FindManyP<recob::Hit> findManyHits(tpcTrackHandle, event, trackLabel);
      for (auto const& tpcTrack : (*tpcTrackHandle)){
	std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(tpcTrack.ID());
        matchcanvec.push_back(GetClosestCRTHit(detProp, tpcTrack, hits, crtHits, trigger_timestamp,IsData));
	//return ClosestCRTHit(detProp, tpcTrack, hits, crtHits);
	//matchCand closestHit = GetClosestCRTHit(detProp, tpcTrack, hits, crtHits);

      }
    }
    return matchcanvec;
    //auto tpcTrackHandle = event.getValidHandle<std::vector<recob::Track>>(fTPCTrackLabel);
    //art::FindManyP<recob::Hit> findManyHits(tpcTrackHandle, event, fTPCTrackLabel);
    //std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(tpcTrack.ID());
    //return GetClosestCRTHit(detProp, tpcTrack, hits, crtHits);
    //    for (const auto& match : matchedCan)
    //return match;
    //return nullmatch;
  }
*/

    std::vector<icarus::match_geometry> CRTT0MatchAlg::GetClosestCRTHit_geo(detinfo::DetectorPropertiesData const& detProp,
					    recob::Track const& tpcTrack, std::pair<double, double> t0MinMax, 
					    std::vector<sbn::crt::CRTHit> const& crtHits, int driftDirection, uint64_t& trigger_timestamp, bool IsData) const{

    auto start = tpcTrack.Vertex();
    auto end   = tpcTrack.End();

    bool simple_cathode_crosscheck =( (std::abs(start.X()) < 210.215) != (std::abs(end.X()) < 210.215));
    int hit_id = 0;

    // ====================== Matching Algorithm ========================== //
    //  std::vector<std::pair<sbn::crt::CRTHit, double>> t0Candidates;
    std::vector<match_geometry> t0Candidates;

    //    if (crtHits.size() == 0) continue;
    // Loop over all the CRT hits
    for(auto &crtHit : crtHits){
      // Check if hit is within the allowed t0 range
      double crtTime = GetCRTTime(crtHit,trigger_timestamp,IsData);  // units are us

     icarus::match_geometry this_candidate;// = makeNULLmc_geo();
/*     if(IsData){
      if (fTSMode == 1) {
	crtTime = ((double)(int)crtHit.ts1_ns) * 1e-3; //+ fTimeCorrection;
      }//end if (fTSMode == 1)
      else {
		crtTime = double(crtHit.ts0_ns - (trigger_timestamp%1'000'000'000))/1e3;
		if(crtTime<-0.5e6) 	crtTime+=1e6;
		else if(crtTime>=0.5e6) 	crtTime-=1e6;		
      }//end else
     }//end if(IsData)
     else if(!IsData){	
	crtTime = crtHit.ts0_ns/1e3;
     }//end else if(!IsData)*/
      // If track is stitched then try all hits
      if (!((crtTime >= t0MinMax.first - 10. && crtTime <= t0MinMax.second + 10.) 
            || t0MinMax.first == t0MinMax.second)) continue;

      // cut on CRT hit PE value
      if (crtHit.peshit<fPEcut) continue;
      if (crtHit.x_err>fMaxUncert) continue;
      if (crtHit.y_err>fMaxUncert) continue;
      if (crtHit.z_err>fMaxUncert) continue;
      if (tpcTrack.Length() < fMinTrackLength) continue;

      geo::Point_t crtPoint(crtHit.x_pos, crtHit.y_pos, crtHit.z_pos);

      //Calculate Track direction
      std::pair<TVector3, TVector3> startEndDir;
      // dirmethod=2 is original algorithm, dirmethod=1 is simple algorithm for which SCE corrections are possible
      if (fDirMethod==2)  startEndDir = TrackDirectionAverage(tpcTrack, fTrackDirectionFrac);
      else startEndDir = TrackDirection(detProp, tpcTrack, fTrackDirectionFrac, crtTime, driftDirection);
      TVector3 startDir = startEndDir.first;
      TVector3 endDir = startEndDir.second;
    
      // Calculate the distance between the crossing point and the CRT hit, SCE corrections are done inside but dropped
      double startDist = DistOfClosestApproach(detProp, start, startDir, crtHit, driftDirection, crtTime);
      double endDist = DistOfClosestApproach(detProp, end, endDir, crtHit, driftDirection, crtTime);

      double xshift = driftDirection * crtTime * detProp.DriftVelocity();
      auto thisstart = start; 
      thisstart.SetX(start.X()+xshift);
      auto thisend = end; 
      thisend.SetX(end.X()+xshift);

      // repeat SCE correction for endpoints
      if (fSCE->EnableCalSpatialSCE() && fSCEposCorr) {
//	geo::Point_t temppt = {thisstart.X(),thisstart.Y(),thisstart.Z()};
	geo::TPCID tpcid = fGeometryService->PositionToTPCID(thisstart);
//	geo::Vector_t  fPosOffsets = fSCE->GetCalPosOffsets(temppt,tpcid.TPC);
	thisstart+= fSCE->GetCalPosOffsets(thisstart,tpcid.TPC);
/*	thisstart[0] += fPosOffsets.X();
	thisstart[1] += fPosOffsets.Y();
	thisstart[2] += fPosOffsets.Z();*/
//	temppt.SetX(thisend.X());
//	temppt.SetY(thisend.Y());
//	temppt.SetZ(thisend.Z());
	tpcid = fGeometryService->PositionToTPCID(thisend);
	thisend+= fSCE->GetCalPosOffsets(thisend,tpcid.TPC);
/*	fPosOffsets = fSCE->GetCalPosOffsets(temppt,tpcid.TPC);
	thisend[0] += fPosOffsets.X();
	thisend[1] += fPosOffsets.Y();
	thisend[2] += fPosOffsets.Z();*/
      }
//      TVector3 thisstart_v(thisstart.X(),thisstart.Y(),thisstart.Z());
//      TVector3 thisend_v(thisend.X(),thisend.Y(),thisend.Z());


      if (startDist<fDistanceLimit || endDist<fDistanceLimit) {
	double distS = (crtPoint-thisstart).R();
	double distE =  (crtPoint-thisend).R();
	if (distS <= distE && startDist<fDistanceLimit){ 
	  this_candidate.dca = startDist;
	  this_candidate.extrapLen = distS;
	  this_candidate.best_DCA_pos=0;
	}//end if(distS < distE)
	else if(distE<=distS && endDist<fDistanceLimit ){
	  this_candidate.dca = endDist;
	  this_candidate.extrapLen = distE;
	  this_candidate.best_DCA_pos=1;
	}//end else if(distE<=distS && endDist<fDistanceLimit )
	else continue;
	this_candidate.thishit = crtHit;
	this_candidate.t0= crtTime;
	this_candidate.simple_cathodecrosser = simple_cathode_crosscheck;
	this_candidate.driftdir = driftDirection;
	this_candidate.t0min = t0MinMax.first;
	this_candidate.t0max = t0MinMax.second;
	this_candidate.crtTime = crtTime;
	this_candidate.startDir = startDir;
	this_candidate.endDir = endDir; 
	this_candidate.hit_id = hit_id; hit_id++;
	this_candidate.track_id = tpcTrack.ID();
	this_candidate.startDir.SetXYZ(startDir.X(),startDir.Y(),startDir.Z());
	this_candidate.endDir.SetXYZ(endDir.X(),endDir.Y(),endDir.Z());
	this_candidate.crt_hit_pos.SetXYZ(crtPoint.X(), crtPoint.Y(), crtPoint.Z());
	this_candidate.simpleDCA_startDir = startDist;
      	this_candidate.simpleDCA_endDir = endDist;
	this_candidate.tpc_track_start.SetXYZ(thisstart.X(),thisstart.Y(),thisstart.Z());
	this_candidate.tpc_track_end.SetXYZ(thisend.X(),thisend.Y(),thisend.Z());
	t0Candidates.push_back(this_candidate);

      }//end if (startDist<fDistanceLimit || endDist<fDistanceLimit) 
    }//end loop over CRT Hits

      //std::cout << " found " << t0Candidates.size() << " candidates" << std::endl;
/*    icarus::match_geometry bestmatch = makeNULLmc_geo();
    icarus::match_geometry thismatch = makeNULLmc_geo();
    if(t0Candidates.size() > 0){
      // Find candidate with shortest DCA or DCA/L value
      bestmatch=t0Candidates[0];
      double sin_angle = bestmatch.dca/bestmatch.extrapLen;
      if (fDCAoverLength) { // Use dca/extrapLen to judge best
	for(auto &thisCand : t0Candidates){
	  double this_sin_angle = thisCand.dca/thisCand.extrapLen;
	  if (bestmatch.dca<0 )bestmatch=thisCand;
	  else if (this_sin_angle<sin_angle && thisCand.dca>=0)bestmatch=thisCand;
	}//end for(auto &thisCand : t0Candidates)
      }//end if (fDCAoverLength)
      else { // use Dca to judge best
	for(auto &thisCand : t0Candidates){
	  if (bestmatch.dca<0 )bestmatch=thisCand;
	  else if (thisCand.dca<bestmatch.dca && thisCand.dca>=0)bestmatch=thisCand;
	}//end for(auto &thisCand : t0Candidates)
      }//end else [use DCA for best match method]
    }//end if(t0Candidates.size() > 0)

    //std::cout << "best match has dca of " << bestmatch.dca << std::endl;
    return bestmatch;*/
    return t0Candidates;

  }//end function defn


    double CRTT0MatchAlg::GetCRTTime(sbn::crt::CRTHit const& crthit, uint64_t trigger_timestamp, bool isdata) const {

	double crtTime =DBL_MAX;

	if(isdata){
		if (fTSMode == 1) {
			crtTime = ((double)(int)crthit.ts1_ns) * 1e-3; //+ fTimeCorrection;
	      	}//end if (fTSMode == 1)
      		else {
			crtTime = double(crthit.ts0_ns - (trigger_timestamp%1'000'000'000))/1e3;
			if(crtTime<-0.5e6) 	crtTime+=1e6;
			else if(crtTime>=0.5e6) crtTime-=1e6;		
      		}//end else
	}//end if(isdata)
     	else{	
		crtTime = (crthit.ts0() - 1600000)/1e3;
	}//end else if(!IsData)

	return crtTime;

    }//end definition of double CRTT0MatchAlg::GetCRTTime(sbn::crt:CRTHit const& crthit, uint64_t trigger_timestamp, bool isdata) const


/////////////////////////////////////////////////////////////

  //This is the function that is typically called by the user to associate CRT Hits with TPC tracks.
  //It differs from the other similar function in a few ways, this one takes in a vector of art::Ptrs to the recob::Hits 
  //inside the recob::Track as well as the track itself. 
  //It pre-processes the recob::Hits to get the TPC drift direction as well as the TPC X limits and the min/max allowed T0
  //given the X range of the recob::Hits in that TPC, which are also fed into the similar (non-vector) function to 
  //assist with the matching process. I'm not 100% why this was the choice originally, I haven't disturbed it
  //as it works for me right now. -TB 03/22/23
  matchCand_PCA CRTT0MatchAlg::GetClosestCRTHit_PCA(detinfo::DetectorPropertiesData const& detProp,
					    recob::Track const& tpcTrack, std::vector<art::Ptr<recob::Hit>> const& hits, 
					    std::vector<sbn::crt::CRTHit> const& crtHits, uint64_t trigger_timestamp, bool IsData) {

    auto start = tpcTrack.Vertex();
    auto end   = tpcTrack.End();



    // Get the drift direction from the TPC
    int driftDirection = TPCGeoUtil::DriftDirectionFromHits(fGeometryService, hits);
    //std::cout << "size of hit in a track: " << hits.size() << ", driftDirection: "<< driftDirection 
    //	      << " , tpc: "<< hits[0]->WireID().TPC << std::endl; //<< " , intpc: "<< icarus::TPCGeoUtil::DetectedInTPC(hits) << std::endl;
    std::pair<double, double> xLimits = TPCGeoUtil::XLimitsFromHits(fGeometryService, hits);
    // Get the allowed t0 range
    std::pair<double, double> t0MinMax = TrackT0Range(detProp, start.X(), end.X(), driftDirection, xLimits);

    return GetClosestCRTHit_PCA(detProp, tpcTrack, t0MinMax, crtHits, driftDirection, trigger_timestamp, IsData);

  }



/*  struct  matchCand_PCA {
    sbn::crt::CRTHit thishit;
    double t0 = DBL_MIN;
    double dca = DBL_MIN;
    double extrapLen = DBL_MIN;
    int best_DCA_pos = -1;//-1=no match; 0=startdir is best; 1=enddir is best; 2=both start+end are equally good;
    bool simple_cathodecrosser = false;
    int driftdir = -5;
    double t0min = DBL_MIN;
    double t0max = DBL_MIN;
    double crtTime = DBL_MIN;
    TVector3 startDir{DBL_MIN,DBL_MIN,DBL_MIN};
    TVector3 endDir{DBL_MIN,DBL_MIN,DBL_MIN};
    TVector3 tpc_track_start{DBL_MIN,DBL_MIN,DBL_MIN};
    TVector3 tpc_track_end{DBL_MIN,DBL_MIN,DBL_MIN};

    TVector3 PCA_start_pos{-DBL_MAX,-DBL_MAX,-DBL_MAX};
    TVector3 PCA_end_pos{-DBL_MAX,-DBL_MAX,-DBL_MAX};
    TVector3 PCA_start_dir{-DBL_MAX,-DBL_MAX,-DBL_MAX};
    TVector3 PCA_end_dir{-DBL_MAX,-DBL_MAX,-DBL_MAX};
    TVector3 PCA_start_crtplanecross{-DBL_MAX,-DBL_MAX,-DBL_MAX};
    TVector3 PCA_end_crtplanecross{-DBL_MAX,-DBL_MAX,-DBL_MAX};

    double PCA_DCA_start = -1;
    double PCA_DCA_end = -1;
    double PCA_planedist_start = -1;
    double PCA_planedist_end = -1;

  };*/

  //This is the main workhorse function, in that it takes in a TPC track and a vector of CRT Hits, then searches for the best match candidate. 
  //Today I'm going to try to implement PCA into this to differ from the results found with the old method. 
  matchCand_PCA CRTT0MatchAlg::GetClosestCRTHit_PCA(detinfo::DetectorPropertiesData const& detProp,
					    recob::Track const& tpcTrack, std::pair<double, double> t0MinMax, 
					    std::vector<sbn::crt::CRTHit> const& crtHits, int driftDirection, uint64_t& trigger_timestamp, bool IsData)  {

    bool useVerbose = false;
    auto start = tpcTrack.Vertex();
    auto end   = tpcTrack.End();

    if(useVerbose){
	    std::cout << "//////////////////////////////////////////////////////////////////////////////////////////\n";
	    std::cout << "TPC Track Start (x,y,z):\t(" << start.X() << "," << start.Y() << "," << start.Z() << ")\n";
	    std::cout << "TPC Track End (x,y,z):\t(" << end.X() << "," << end.Y() << "," << end.Z() << ")\n";
	    std::cout << "TPC Track length: " << tpcTrack.Length() << " cm\n";
	    std::cout << "(T0Min,T0Max)=(" << t0MinMax.first << "," << t0MinMax.second << ")\n";
	    std::cout << "Number of candidate CRT Hits: " << crtHits.size() << std::endl;
    }//end if(useVerbose)

    bool simple_cathode_crosscheck =( (std::abs(start.X()) < 210.215) != (std::abs(end.X()) < 210.215));
    // ====================== Matching Algorithm ========================== //
    //  std::vector<std::pair<sbn::crt::CRTHit, double>> t0Candidates;
    std::vector<matchCand_PCA> t0Candidates;

    TVector3 pca_start_pos, pca_end_pos, pca_start_dir, pca_end_dir, pca_start_crtcross, pca_end_crtcross;
    endpoint_PCA_ana(tpcTrack,true,pca_start_dir,pca_start_pos);
    endpoint_PCA_ana(tpcTrack,true,pca_end_dir,pca_end_pos);

    int count_crt_hits = 0;

    // Loop over all the CRT hits
    for(auto &crtHit : crtHits){

      TVector3 crthitposvec(crtHit.x_pos,crtHit.y_pos,crtHit.z_pos);

      if(useVerbose) std::cout << "CRT Hit #" << count_crt_hits << "position (x,y,z): (" << crtHit.x_pos << "," << crtHit.y_pos << "," << crtHit.z_pos << ")\n";

      // Check if hit is within the allowed t0 range
      double crtTime = GetCRTTime(crtHit,trigger_timestamp,IsData);  // units are us
      std::cout << "CRT Hit #" << count_crt_hits << " timestamp:\t" << crtTime << std::endl;

      // If track is stitched then try all hits
      if (!((crtTime >= t0MinMax.first - 10. && crtTime <= t0MinMax.second + 10.) 
            || t0MinMax.first == t0MinMax.second)) continue;

      // cut on CRT hit PE value
      if (crtHit.x_err>fMaxUncert) { 
		if(useVerbose) std::cout << "CRT Hit #" << count_crt_hits << " disqualified because X error is too large!\n";
		count_crt_hits++;
		continue; 
      }//end if (crtHit.x_err>fMaxUncert)
      if (crtHit.y_err>fMaxUncert) { 
		if(useVerbose) std::cout << "CRT Hit #" << count_crt_hits << " disqualified because Y error is too large!\n";
		count_crt_hits++;
		continue; 
      }//end if (crtHit.y_err>fMaxUncert)
      if (crtHit.z_err>fMaxUncert) { 
		if(useVerbose) std::cout << "CRT Hit #" << count_crt_hits << " disqualified because Z error is too large!\n";
		count_crt_hits++;
		continue; 
      }//end if (crtHit.z_err>fMaxUncert)
      if (tpcTrack.Length() < fMinTrackLength) { 
		if(useVerbose) std::cout << "CRT Hit #" << count_crt_hits << " disqualified because track length is too small!\n";
		count_crt_hits++;
		continue; 
      }//end if (tpcTrack.Length() < fMinTrackLength)

      geo::Point_t crtPoint(crtHit.x_pos, crtHit.y_pos, crtHit.z_pos);
    
      int crt_reg = AuxDetRegionNameToNum(crtHit.tagger);
      int thinrange;
      if(crt_reg==30) thinrange=1;
      else if(crt_reg==31||crt_reg==32||(crt_reg>=40&&crt_reg<=45)) thinrange = 0; 
      else if(crt_reg==33 || crt_reg==34 || crt_reg ==46 ||crt_reg==47) thinrange = 2; 

      std::pair<TVector3, TVector3> startEndDir;
      // dirmethod=2 is original algorithm, dirmethod=1 is simple algorithm for which SCE corrections are possible
      startEndDir = TrackDirection(detProp, tpcTrack, fTrackDirectionFrac, crtTime, driftDirection);
      TVector3 startDir = startEndDir.first;
      TVector3 endDir = startEndDir.second;
    
      // Calculate the distance of closest approach between TPC track vector and the CRT hit, SCE corrections are done inside but dropped
      double startDist_old = DistOfClosestApproach(detProp, start, startDir, crtHit, driftDirection, crtTime);
      double endDist_old = DistOfClosestApproach(detProp, end, endDir, crtHit, driftDirection, crtTime);

      // Calculate the distance between the crossing point and the CRT hit
      double xshift = driftDirection * crtTime * detProp.DriftVelocity();// std::cout << "drift velocity: " << detProp.DriftVelocity() << std::endl;
      TVector3 startposvec(pca_start_pos.X()+xshift,pca_start_pos.Y(),pca_start_pos.Z());
      TVector3 endposvec(pca_end_pos.X()+xshift,pca_end_pos.Y(),pca_end_pos.Z());

      double startDist = CRT_plane_dist(startposvec,pca_start_dir,crthitposvec,thinrange,pca_start_crtcross);
      double endDist = CRT_plane_dist(endposvec,pca_end_dir,crthitposvec,thinrange,pca_end_crtcross);
    
      if(useVerbose) std::cout << "CRT Hit #" << count_crt_hits << " (startDist,endDist)=(" << startDist << "," << endDist << ")\n";
      if(useVerbose) std::cout << "CRT Hit #" << count_crt_hits << " xshift for TPC Track: " << xshift << std::endl;

      count_crt_hits++;

      auto thisstart = start; 
      thisstart.SetX(start.X()+xshift);
      auto thisend = end; 
      thisend.SetX(end.X()+xshift);
      // repeat SCE correction for endpoints
      if (fSCE->EnableCalSpatialSCE() && fSCEposCorr) {
	geo::TPCID tpcid = fGeometryService->PositionToTPCID(thisstart);
	thisstart+= fSCE->GetCalPosOffsets(thisstart,tpcid.TPC);
	tpcid = fGeometryService->PositionToTPCID(thisend);
	thisend+= fSCE->GetCalPosOffsets(thisend,tpcid.TPC);

      }

      matchCand_PCA newmc;// = makeNULLmc();
      if (startDist<fDistanceLimit || endDist<fDistanceLimit) {
	double distS = startDist; //(crtPoint-thisstart).R();
	double distE = startDist; //(crtPoint-thisend).R();
	if (distS <= distE && startDist<fDistanceLimit){ 
	  newmc.dca = startDist;
	  newmc.extrapLen = distS;
	  newmc.best_DCA_pos=0;
	}//end if(distS < distE)
	else if(distE<=distS && endDist<fDistanceLimit ){
	  newmc.dca = endDist;
	  newmc.extrapLen = distE;
	  newmc.best_DCA_pos=1;
	}//end else if(distE<=distS && endDist<fDistanceLimit )
	else continue;
	newmc.thishit = crtHit;
	newmc.t0= crtTime;
	newmc.simple_cathodecrosser = simple_cathode_crosscheck;
	newmc.t0min = t0MinMax.first;
	newmc.t0max = t0MinMax.second;
	newmc.crtTime = crtTime;
	newmc.startDir = startDir;
	newmc.endDir = endDir;
	newmc.tpc_track_start.SetXYZ(thisstart.X(),thisstart.Y(),thisstart.Z());
	newmc.tpc_track_end.SetXYZ(thisend.X(),thisend.Y(),thisend.Z());
	newmc.PCA_start_pos.SetXYZ(startposvec.X(),startposvec.Y(),startposvec.Z());
	newmc.PCA_end_pos.SetXYZ(endposvec.X(),endposvec.Y(),endposvec.Z());
	newmc.PCA_start_dir.SetXYZ(pca_start_dir.X(),pca_start_dir.Y(),pca_start_dir.Z());
	newmc.PCA_end_dir.SetXYZ(pca_end_dir.X(),pca_end_dir.Y(),pca_end_dir.Z());
	newmc.PCA_start_crtplanecross.SetXYZ(pca_start_crtcross.X(),pca_start_crtcross.Y(),pca_start_crtcross.Z());
	newmc.PCA_end_crtplanecross.SetXYZ(pca_end_crtcross.X(),pca_end_crtcross.Y(),pca_end_crtcross.Z());
	newmc.driftdir = driftDirection;
	newmc.PCA_DCA_start = startDist_old;
	newmc.PCA_DCA_end = endDist_old;
	newmc.PCA_planedist_start = startDist;
	newmc.PCA_planedist_end = endDist;

	t0Candidates.push_back(newmc);

      }//end if (startDist<fDistanceLimit || endDist<fDistanceLimit) 
    }//end loop over CRT Hits


      //std::cout << " found " << t0Candidates.size() << " candidates" << std::endl;
    matchCand_PCA bestmatch;// = makeNULLmc();
    if(t0Candidates.size() > 0){
      // Find candidate with shortest DCA or DCA/L value
      bestmatch=t0Candidates[0];
      double sin_angle = bestmatch.dca/bestmatch.extrapLen;
      if (fDCAoverLength) { // Use dca/extrapLen to judge best
	for(auto &thisCand : t0Candidates){
	  double this_sin_angle = thisCand.dca/thisCand.extrapLen;
	  if (bestmatch.dca<0 )bestmatch=thisCand;
	  else if (this_sin_angle<sin_angle && thisCand.dca>=0)bestmatch=thisCand;
	}//end for(auto &thisCand : t0Candidates)
      }//end if (fDCAoverLength)
      else { // use Dca to judge best
	for(auto &thisCand : t0Candidates){
	  if (bestmatch.dca<0 )bestmatch=thisCand;
	  else if (thisCand.dca<bestmatch.dca && thisCand.dca>=0)bestmatch=thisCand;
	}//end for(auto &thisCand : t0Candidates)
      }//end else [use DCA for best match method]
    }//end if(t0Candidates.size() > 0)

    if(useVerbose) std::cout << "best match has dca of " << bestmatch.dca << std::endl;
    if(useVerbose) std::cout <<"########################################################";

    return bestmatch;

  }//end function defn

/*
 *
 *
   matchCand_PCA CRTT0MatchAlg::GetClosestCRTHit_PCA(detinfo::DetectorPropertiesData const& detProp,
					    recob::Track const& tpcTrack, std::vector<art::Ptr<recob::Hit>> const& hits, 
					    std::vector<sbn::crt::CRTHit> const& crtHits, uint64_t trigger_timestamp, bool IsData) {

    auto start = tpcTrack.Vertex();
    auto end   = tpcTrack.End();



    // Get the drift direction from the TPC
    int driftDirection = TPCGeoUtil::DriftDirectionFromHits(fGeometryService, hits);
    //std::cout << "size of hit in a track: " << hits.size() << ", driftDirection: "<< driftDirection 
    //	      << " , tpc: "<< hits[0]->WireID().TPC << std::endl; //<< " , intpc: "<< icarus::TPCGeoUtil::DetectedInTPC(hits) << std::endl;
    std::pair<double, double> xLimits = TPCGeoUtil::XLimitsFromHits(fGeometryService, hits);
    // Get the allowed t0 range
    std::pair<double, double> t0MinMax = TrackT0Range(detProp, start.X(), end.X(), driftDirection, xLimits);

    return GetClosestCRTHit_PCA(detProp, tpcTrack, t0MinMax, crtHits, driftDirection, trigger_timestamp, IsData);

  }

 *
 *
 * */


//Attempt at inverting the search algorithm to find the best TPC track match based on the CRT Hit position
//This will likely be rough to start but will be worked on as new ideas come to me -TB 04/03/23
tpc_cand CRTT0MatchAlg::GetClosestTPCTrack(detinfo::DetectorPropertiesData const& detProp, sbn::crt::CRTHit const& this_crthit, 
				std::vector<recob::Track> const& tpcTracks, std::vector<std::vector<art::Ptr<recob::Hit>>> const& allHits, //should contain a vector of tracks as well as a vector of vectors of Hits
				uint64_t& trigger_timestamp, bool IsData) const {

	size_t numtrks = tpcTracks.size();
	for(size_t thistrk = 0; thistrk < numtrks; thistrk++){

		std::vector<art::Ptr<recob::Hit>> thesehits = allHits[thistrk];

		int driftDirection = TPCGeoUtil::DriftDirectionFromHits(fGeometryService,thesehits);	
		std::pair<double,double> xLimits = TPCGeoUtil::XLimitsFromHits(fGeometryService, thesehits);
		std::pair<double,double> t0MinMax = TrackT0Range(detProp, tpcTracks[thistrk].Vertex().X(), tpcTracks[thistrk].End().X(), driftDirection, xLimits);
		std::cout << t0MinMax.first << t0MinMax.second;
	}//end for(size_t thistrk = 0; thistrk < numtrks; thistrk++)

	tpc_cand thiscand; 
	return thiscand;



}//end definition of GetClosestTPCTrack function



//First new added function is to do the PCA analysis on a given endpoint
void CRTT0MatchAlg::endpoint_PCA_ana(recob::Track trk, bool usestartpt, TVector3 &bestfit_dir, TVector3 &bestfit_pos){


//	size_t ntrk = trk.NPoints(); size_t nvalidtrk = trk.CountValidPoints();
	Eigen::Vector3d meanPos(Eigen::Vector3d::Zero());

	//create container for all hit point locations as well as calculate average positions along all axes
	float xavg = 0, yavg = 0, zavg = 0;
	std::vector<TVector3> hit_points;

	size_t firstpt;
	geo::Point_t thispt; geo::Point_t lastpt;
	TVector3 pos_of_endpoint;

	if(usestartpt) {
		firstpt = trk.FirstValidPoint(); 
		thispt = trk.LocationAtPoint(firstpt);
		xavg+=thispt.X(); yavg+=thispt.Y(); zavg+=thispt.Z();
		TVector3 temppt(thispt.X(),thispt.Y(),thispt.Z());
		pos_of_endpoint.SetXYZ(thispt.X(),thispt.Y(),thispt.Z());
		hit_points.push_back(temppt);
	}//end if(usestartpt)
	else if(!usestartpt) { 
		firstpt = trk.LastValidPoint(); 
		thispt = trk.LocationAtPoint(firstpt);
		xavg+=thispt.X(); yavg+=thispt.Y(); zavg+=thispt.Z();
		TVector3 temppt(thispt.X(),thispt.Y(),thispt.Z());
		pos_of_endpoint.SetXYZ(thispt.X(),thispt.Y(),thispt.Z());
		hit_points.push_back(temppt);
	}//end else if(!usestartpt)

	float distance_from_endpoint = 0;
	float distlimit = 10;
	int num_points_counter = 1;

	bool forcedbreak = false;


	while(distance_from_endpoint<distlimit || (size_t)num_points_counter<40){

		size_t nextpt;
		if(usestartpt) {
			nextpt = trk.NextValidPoint(firstpt+1);
		}//end if(usestartpt)
		else if(!usestartpt) { 
			nextpt = trk.PreviousValidPoint(firstpt-1);
		}//end else if(!usestartpt)

		if(nextpt==firstpt) {

			forcedbreak = true;
			break;

		}//end if(nextpt==firstpt)

		geo::Point_t thispt = trk.LocationAtPoint(nextpt);


		xavg+=thispt.X(); yavg+=thispt.Y(); zavg+=thispt.Z();
		TVector3 temppt(thispt.X(),thispt.Y(),thispt.Z());

		distance_from_endpoint = vect_dist(pos_of_endpoint,temppt);

		hit_points.push_back(temppt);

		firstpt = nextpt;


		num_points_counter++;
	}//end loop over points in TPC track, adding to what will eventually go in PCA analysis

	if(!forcedbreak) {



	xavg = xavg/num_points_counter; yavg = yavg/num_points_counter; zavg = zavg/num_points_counter;
	meanPos(0) = xavg; meanPos(1) = yavg; meanPos(2) = zavg;
	bestfit_pos.SetX(xavg); bestfit_pos.SetY(yavg); bestfit_pos.SetZ(zavg);

	//Make covariance matrix now that we have averages
	float xi2 = 0, yi2 = 0, zi2 = 0, xiy = 0, xiz = 0, yiz = 0;
	for(size_t i=0; i<hit_points.size(); i++){

		float tempx = hit_points[i].X() - xavg;
		float tempy = hit_points[i].Y() - yavg;
		float tempz = hit_points[i].Z() - zavg;

		xi2+=(tempx*tempx); yi2+=(tempy*tempy); zi2+=(tempz*tempz);
		xiy+=(tempx*tempy); xiz+=(tempx*tempz); yiz+=(tempy*tempz);

	}//end loop over points in data
	xi2 = xi2/hit_points.size(); yi2 = yi2/hit_points.size(); zi2 = zi2/hit_points.size();
	xiy = xiy/hit_points.size(); xiz = xiz/hit_points.size(); yiz = yiz/hit_points.size();

	Eigen::Matrix3d covmat; 
	covmat << xi2, xiy, xiz, xiy, yi2, yiz, xiz, yiz, zi2;

	Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigenMat(covmat);

	if(eigenMat.info() == Eigen::ComputationInfo::Success){

	      reco::PrincipalComponents::EigenValues recobEigenVals = eigenMat.eigenvalues().cast<float>();
	      reco::PrincipalComponents::EigenVectors recobEigenVecs = eigenMat.eigenvectors().cast<float>();

	      float maxeval = FLT_MIN; int maxvalpos = -1;

	      for(int tempindex=0; tempindex<3; tempindex++){

			if(maxeval<(float)recobEigenVals(tempindex)) {

				maxeval = (float)recobEigenVals(tempindex);
				maxvalpos = tempindex;

			}//end if(maxeval<(float)recobEigenVals)

	      }//end for(int tempindex=0; tempindex<3; tempindex++)

	      bestfit_dir.SetXYZ(recobEigenVecs(0,maxvalpos), recobEigenVecs(1,maxvalpos), recobEigenVecs(2,maxvalpos));

	}//end if(eigenMat.info() == Eigen::ComputationInfo::Success)

	}//end if(!forcedbreak)

}//end definition of void endpoint_PCA_ana(recob::Track trk, int startorend, TVector3 &bestfit_dir, TVector3 &bestfit_pos)

//Second function, a quick tool to calculate distance between two TVector3 objects
float CRTT0MatchAlg::vect_dist(TVector3 vec1, TVector3 vec2){

//	std::cout << "finding distance between two vectors.\n";
//	std::cout << "vector 1 (x,y,z):\t(" << vec1.X() << "," << vec1.Y() << "," << vec1.Z() << ")\n"; 
//	std::cout << "vector 2 (x,y,z):\t(" << vec2.X() << "," << vec2.Y() << "," << vec2.Z() << ")\n"; 

	float xdiff = vec1.X()-vec2.X();
	float ydiff = vec1.Y()-vec2.Y();
	float zdiff = vec1.Z()-vec2.Z();

//	std::cout << "Difference in X:\t" << xdiff << std::endl;
//	std::cout << "Difference in Y:\t" << ydiff << std::endl;
//	std::cout << "Difference in Z:\t" << zdiff << std::endl;

	xdiff=xdiff*xdiff; ydiff=ydiff*ydiff; zdiff=zdiff*zdiff;

//	std::cout << "Difference^2 in X:\t" << xdiff << std::endl;
//	std::cout << "Difference^2 in Y:\t" << ydiff << std::endl;
//	std::cout << "Difference^2 in Z:\t" << zdiff << std::endl;

	float returnval = std::sqrt(xdiff+ydiff+zdiff);

//	std::cout << "distance betwen vectors found to be:\t" << returnval << std::endl;

	return returnval;


}//end definition of double CRTT0MatchAlg::vect_dist(TVector3 vec1, TVector3 vec2)

//Third addition, a function to calculate the distance between a CRT Hit and where
//a given TPC track will intersect the plane of the CRT, which *should*
//be the "real" hit position
double CRTT0MatchAlg::CRT_plane_dist(TVector3 TPC_track_pos, TVector3 TPC_track_dir, TVector3 CRT_hit_pos, int plane_axis/*0 for x, 1 for y, 2 for z*/, TVector3 &TPC_track_crosspoint){

	double tfactor, returnval = -1;

	if(plane_axis==0){

		tfactor = (CRT_hit_pos.X() - TPC_track_pos.X())/TPC_track_dir.X();
		
		TPC_track_crosspoint.SetX(CRT_hit_pos.X());
		TPC_track_crosspoint.SetY(TPC_track_pos.Y()+TPC_track_dir.Y()*tfactor);
		TPC_track_crosspoint.SetZ(TPC_track_pos.Z()+TPC_track_dir.Z()*tfactor);

		double ydiff, zdiff;
		ydiff = TPC_track_crosspoint.Y() - CRT_hit_pos.Y();
		zdiff = TPC_track_crosspoint.Z() - CRT_hit_pos.Z();

		returnval = sqrt(ydiff*ydiff + zdiff*zdiff);

	}//end if(plane_axis==0)
	else if(plane_axis==1){

		tfactor = (CRT_hit_pos.Y() - TPC_track_pos.Y())/TPC_track_dir.Y();
		
		TPC_track_crosspoint.SetX(TPC_track_pos.X()+TPC_track_dir.X()*tfactor);
		TPC_track_crosspoint.SetY(CRT_hit_pos.Y());
		TPC_track_crosspoint.SetZ(TPC_track_pos.Z()+TPC_track_dir.Z()*tfactor);

		double xdiff, zdiff;
		xdiff = TPC_track_crosspoint.X() - CRT_hit_pos.X();
		zdiff = TPC_track_crosspoint.Z() - CRT_hit_pos.Z();

		returnval = sqrt(xdiff*xdiff + zdiff*zdiff);

	}//end else if (plane_axis==1)
	else if(plane_axis==2){

		tfactor = (CRT_hit_pos.Z() - TPC_track_pos.Z())/TPC_track_dir.Z();
		
		TPC_track_crosspoint.SetX(TPC_track_pos.X()+TPC_track_dir.X()*tfactor);
		TPC_track_crosspoint.SetY(TPC_track_pos.Y()+TPC_track_dir.Y()*tfactor);
		TPC_track_crosspoint.SetZ(CRT_hit_pos.Z());

		double ydiff, xdiff;
		ydiff = TPC_track_crosspoint.Y() - CRT_hit_pos.Y();
		xdiff = TPC_track_crosspoint.X() - CRT_hit_pos.X();

		returnval = sqrt(ydiff*ydiff + xdiff*xdiff);

	}//end else if(plane_axis==2)

	return returnval;


}//end double CRTT0MatchAlg::CRT_plane_dist(TVector3 TPC_track_pos, TVector3 TPC_track_dir, TVector3 CRT_hit_pos, int plane_axis/*0 for x, 1 for y, 2 for z*/)

int CRTT0MatchAlg::AuxDetRegionNameToNum(string reg)
{
    if(reg == "Top")        return 30;
    if(reg == "RimWest")    return 31;
    if(reg == "RimEast")    return 32;
    if(reg == "RimSouth")   return 33;
    if(reg == "RimNorth")   return 34;
    if(reg == "WestSouth")  return 40;
    if(reg == "WestCenter") return 41;
    if(reg == "WestNorth")  return 42;
    if(reg == "EastSouth")  return 43;
    if(reg == "EastCenter") return 44;
    if(reg == "EastNorth")  return 45;
    if(reg == "South")      return 46;
    if(reg == "North")      return 47;
    if(reg == "Bottom")     return 50;
    mf::LogError("CRT") << "region not found!" << '\n';
    return INT_MAX;
}//GetAuxDetRegionNum()















}
