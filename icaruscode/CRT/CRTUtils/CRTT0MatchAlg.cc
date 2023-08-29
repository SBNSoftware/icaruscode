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
    fTSMode             = pset.get<int>("TSMode", 1);
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
		crtTime = crthit.ts0_ns/1e3;
	}//end else if(!IsData)

	return crtTime;

    }//end definition of double CRTT0MatchAlg::GetCRTTime(sbn::crt:CRTHit const& crthit, uint64_t trigger_timestamp, bool isdata) const


}
