#include "CRTT0MatchAlg_evtdisp.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()

namespace icarus{

  CRTT0MatchAlg_evtdisp::CRTT0MatchAlg_evtdisp(const fhicl::ParameterSet& pset)
  {
    this->reconfigure(pset);
    return;
  }

  CRTT0MatchAlg_evtdisp::CRTT0MatchAlg_evtdisp() = default;
  
  void CRTT0MatchAlg_evtdisp::reconfigure(const fhicl::ParameterSet& pset){

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

    fGeometryService    = lar::providerFrom<geo::Geometry>();//GeometryService;
    fSCE                = lar::providerFrom<spacecharge::SpaceChargeService>();

    return;

  }
 
  match_geometry makeNULLmc_evtdisp (){
    sbn::crt::CRTHit hit;
    match_geometry null;
    null.thishit = hit;
    null.t0 = -99999;
    null.dca = -99999;
    null.extrapLen = -99999;

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

    null.is_best_DCA_startDir = false;
    null.is_best_DCA_endDir = false;

    return null;
  }

  // Utility function that determines the possible t0 range of a track
  std::pair<double, double> CRTT0MatchAlg_evtdisp::TrackT0Range(detinfo::DetectorPropertiesData const& detProp,
							double startX, double endX, int driftDirection, 
							std::pair<double, double> xLimits){

    // If track is stitched return zeros
    if(driftDirection == 0) return std::make_pair(0, 0);

    //std::pair<double, double> result; // unused
    double Vd = driftDirection * detProp.DriftVelocity();

//  std::cout << " [ driftdirn, vd ] = [ " << driftDirection << " , " << Vd << " ]" << std::endl;

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

    //  if (t0min>2500)  std::cout << " t0 min " << t0min << " t0max " << t0max << std::endl;
    return std::make_pair(std::min(t0min, t0max), std::max(t0min, t0max));


  } // CRTT0MatchAlg_evtdisp::TrackT0Range()


  double CRTT0MatchAlg_evtdisp::DistOfClosestApproach(detinfo::DetectorPropertiesData const& detProp,
					      TVector3 trackPos, TVector3 trackDir, 
					      sbn::crt::CRTHit crtHit, int driftDirection, double t0){

    //double minDist = 99999;

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

    // calculate distance of closest approach (DCA)
    //  default is the distance to the point specified by the CRT hit (Simple DCA)
    //    useBox is the distance to the closest edge of the rectangle with the CRT hit at the center and the sides defined
    //   the position uncertainties on the CRT hits.
    double thisdca;

    if (fDCAuseBox) thisdca =   DistToCrtHit(crtHit, trackPos, end);
    else thisdca =  SimpleDCA(crtHit, trackPos, trackDir);
    return thisdca;

  } // CRTT0MatchAlg_evtdisp::DistToOfClosestApproach()


  std::pair<TVector3, TVector3> CRTT0MatchAlg_evtdisp::TrackDirectionAverage(recob::Track track, double frac)
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

  } // CRTT0MatchAlg_evtdisp::TrackDirectionAverage()


  std::pair<TVector3, TVector3> CRTT0MatchAlg_evtdisp::TrackDirection(detinfo::DetectorPropertiesData const& detProp,
							      recob::Track track, double frac, 
							      double CRTtime, int driftDirection){
          
    size_t nTrackPoints = track.NPoints();
    int midPt = (int)floor(nTrackPoints*frac);
    geo::Point_t startP = track.Start();
    geo::Point_t endP = track.End();
    geo::Point_t midP = track.LocationAtPoint(midPt);

    double xshift = driftDirection * CRTtime * detProp.DriftVelocity();
    TVector3  startPoint = {startP.X()+xshift,startP.Y(),startP.Z()};
    TVector3  endPoint = {endP.X()+xshift,endP.Y(),endP.Z()};
    TVector3  midPoint = {midP.X()+xshift,midP.Y(),midP.Z()};
    if (fSCE->EnableCalSpatialSCE() && fSCEposCorr) {

      // Apply the shift depending on which TPC the track is in                                 
      geo::Point_t fTrackPos = startP;
      fTrackPos.SetX(startPoint.X());
      geo::TPCID tpcid = fGeometryService->PositionToTPCID(fTrackPos);                        
      geo::Vector_t fPosOffsets = fSCE->GetCalPosOffsets(geo::Point_t{fTrackPos.X(),fTrackPos.Y(),fTrackPos.Z()},tpcid.TPC);
      startPoint.SetX(fTrackPos.X() + fPosOffsets.X());                                       
      startPoint.SetY(fTrackPos.Y() + fPosOffsets.Y());                                       
      startPoint.SetZ(fTrackPos.Z() + fPosOffsets.Z());                                       
      fTrackPos = endP;
      fTrackPos.SetX(endPoint.X());
      tpcid = fGeometryService->PositionToTPCID(fTrackPos);
      //      fPosOffsets = fSCE->GetCalPosOffsets(fTrackPos,tpcid.TPC);
      fPosOffsets = fSCE->GetCalPosOffsets(geo::Point_t{fTrackPos.X(),fTrackPos.Y(),fTrackPos.Z()},tpcid.TPC);
      endPoint.SetX(fTrackPos.X() + fPosOffsets.X());
      endPoint.SetY(fTrackPos.Y() + fPosOffsets.Y());
      endPoint.SetZ(fTrackPos.Z() + fPosOffsets.Z());
      fTrackPos = midP;
      fTrackPos.SetX(midPoint.X());
      tpcid = fGeometryService->PositionToTPCID(fTrackPos);
      //fPosOffsets = fSCE->GetCalPosOffsets(fTrackPos,tpcid.TPC);
      fPosOffsets = fSCE->GetCalPosOffsets(geo::Point_t{fTrackPos.X(),fTrackPos.Y(),fTrackPos.Z()},tpcid.TPC);
      midPoint.SetX(fTrackPos.X() + fPosOffsets.X());
      midPoint.SetY(fTrackPos.Y() + fPosOffsets.Y());
      midPoint.SetZ(fTrackPos.Z() + fPosOffsets.Z());
    }
    
    TVector3 startDir = {midPoint.X()-startPoint.X(),midPoint.Y()-startPoint.Y(),midPoint.Z()-startPoint.Z()};
    float norm = startDir.Mag();
    if (norm>0)  startDir *=(1.0/norm);
    TVector3 endDir = {midPoint.X()-endPoint.X(),midPoint.Y()-endPoint.Y(),midPoint.Z()-endPoint.Z()};    
    norm = endDir.Mag();
    if (norm>0)  endDir *=(1.0/norm);
    
    return std::make_pair(startDir, endDir);
    
  } // CRTT0MatchAlg_evtdisp::TrackDirection()                                                                  

  // Keeping ClosestCRTHit function for backward compatibility only
  std::vector<match_geometry> CRTT0MatchAlg_evtdisp::GetClosestCRTHit(detinfo::DetectorPropertiesData const& detProp,
					    recob::Track tpcTrack, std::vector<art::Ptr<recob::Hit>> hits, 
					    std::vector<sbn::crt::CRTHit> crtHits, uint64_t trigger_timestamp, const art::Event& event, bool IsData) {

    auto start = tpcTrack.Vertex<TVector3>();
    auto end   = tpcTrack.End<TVector3>();

    // Get the drift direction from the TPC
    int driftDirection = TPCGeoUtil::DriftDirectionFromHits(fGeometryService, hits);
//  std::cout << "size of hit in a track: " << hits.size() << ", driftDirection: "<< driftDirection 
//	      << " , tpc: "<< hits[0]->WireID().TPC << std::endl; //<< " , intpc: "<< icarus::TPCGeoUtil::DetectedInTPC(hits) << std::endl;
    std::pair<double, double> xLimits = TPCGeoUtil::XLimitsFromHits(fGeometryService, hits);
    // Get the allowed t0 range
    std::pair<double, double> t0MinMax = TrackT0Range(detProp, start.X(), end.X(), driftDirection, xLimits);

    return GetClosestCRTHit(detProp, tpcTrack, t0MinMax, crtHits, driftDirection, trigger_timestamp, event, IsData);

  }

  std::vector<match_geometry> CRTT0MatchAlg_evtdisp::GetClosestCRTHit(detinfo::DetectorPropertiesData const& detProp,
								      recob::Track tpcTrack, std::vector<sbn::crt::CRTHit> crtHits, 	
								      const art::Event& event, uint64_t trigger_timestamp, bool IsData) {
    std::vector<match_geometry> matchcanvec;
    for(const auto& trackLabel : fTPCTrackLabel){
      auto tpcTrackHandle = event.getValidHandle<std::vector<recob::Track>>(trackLabel);
      if (!tpcTrackHandle.isValid()) continue;

      art::FindManyP<recob::Hit> findManyHits(tpcTrackHandle, event, trackLabel);
      for (auto const& tpcTrack : (*tpcTrackHandle)){
	std::vector<match_geometry> tempgeovec;
	std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(tpcTrack.ID());
        tempgeovec = GetClosestCRTHit(detProp, tpcTrack, hits, crtHits, trigger_timestamp,event, IsData);
	matchcanvec.insert(std::end(matchcanvec),std::begin(tempgeovec),std::end(tempgeovec));
	tempgeovec.clear();

      }
    }
    return matchcanvec;
  }


  std::vector<match_geometry> CRTT0MatchAlg_evtdisp::GetClosestCRTHit(detinfo::DetectorPropertiesData const& detProp,
					    recob::Track tpcTrack, std::pair<double, double> t0MinMax, 
					    std::vector<sbn::crt::CRTHit> crtHits, int driftDirection, uint64_t& trigger_timestamp, const art::Event& event, bool IsData) {

    auto start = tpcTrack.Vertex<TVector3>();
    auto end   = tpcTrack.End<TVector3>();

    float temp_start_x, temp_end_x; temp_start_x = start.X(); temp_end_x = end.X();

    // ====================== Matching Algorithm ========================== //
    //  std::vector<std::pair<sbn::crt::CRTHit, double>> t0Candidates;
    std::vector<match_geometry> t0Candidates;
    icarus::match_geometry this_candidate = makeNULLmc_evtdisp();

    int hit_id = 0;

    // Loop over all the CRT hits
    for(auto &crtHit : crtHits){
      // Check if hit is within the allowed t0 range
      double crtTime = -99999.;  // units are us
      if (fTSMode == 1) {
	crtTime = ((double)(int)crtHit.ts1_ns) * 1e-3 + fTimeCorrection;
      }
      else {

	if(IsData){
		crtTime = double(trigger_timestamp - crtHit.ts0_ns)/1e3;
		crtTime = -crtTime+1e6;
	}//end if(IsData)
	else if(!IsData){
		crtTime = 1600 - crtHit.ts0_ns/1e3;
	}//end else if (!IsData)
      }
      // If track is stitched then try all hits
      if (!((crtTime >= t0MinMax.first - 10. && crtTime <= t0MinMax.second + 10.) 
            || t0MinMax.first == t0MinMax.second)) { /*std::cout <<"hit not within T0 range\n";*/ continue; }

      std::cout << "passed ....................... " << std::endl;
      // cut on CRT hit PE value
      if (crtHit.peshit<fPEcut) continue;
      if (crtHit.x_err>fMaxUncert) continue;
      if (crtHit.y_err>fMaxUncert) continue;
      if (crtHit.z_err>fMaxUncert) continue;

      TVector3 crtPoint(crtHit.x_pos, crtHit.y_pos, crtHit.z_pos);
  
      //Calculate Track direction
      std::pair<TVector3, TVector3> startEndDir;
      // dirmethod=2 is original algorithm, dirmethod=1 is simple algorithm for which SCE corrections are possible
      if (fDirMethod==2)  startEndDir = TrackDirectionAverage(tpcTrack, fTrackDirectionFrac);
      else startEndDir = TrackDirection(detProp, tpcTrack, fTrackDirectionFrac, crtTime, driftDirection);
      TVector3 startDir = startEndDir.first;
      TVector3 endDir = startEndDir.second;
 
	bool check_east_1 = temp_start_x < -210.29 && temp_end_x > -210.14;
	bool check_east_2 = temp_end_x < -210.29 && temp_start_x > -210.14;
	bool check_west_1 = temp_start_x < 210.14 && temp_end_x > 210.29;
	bool check_west_2 = temp_end_x < 210.14 && temp_start_x > 210.29;
	
	bool check_east = check_east_1 || check_east_2;
	bool check_west = check_west_1 || check_west_2;

	bool simple_cathode_crosscheck = check_east || check_west;


	this_candidate.simple_cathodecrosser = simple_cathode_crosscheck;
	this_candidate.t0min = t0MinMax.first;
	this_candidate.t0max = t0MinMax.second;
	this_candidate.crtTime = crtTime;

	this_candidate.hit_id = hit_id; hit_id++;
	this_candidate.track_id = tpcTrack.ID();

	this_candidate.startDir.SetXYZ(startDir.X(),startDir.Y(),startDir.Z());
	this_candidate.endDir.SetXYZ(endDir.X(),endDir.Y(),endDir.Z());

	this_candidate.crt_hit_pos.SetXYZ(crtPoint.X(), crtPoint.Y(), crtPoint.Z());

   
      // Calculate the distance between the crossing point and the CRT hit, SCE corrections are done inside but dropped
      double startDist = DistOfClosestApproach(detProp, start, startDir, crtHit, driftDirection, crtTime);
      this_candidate.simpleDCA_startDir = startDist;
      double endDist = DistOfClosestApproach(detProp, end, endDir, crtHit, driftDirection, crtTime);
      this_candidate.simpleDCA_endDir = endDist;

    
      double xshift = driftDirection * crtTime * detProp.DriftVelocity();
      auto thisstart = start; 
      thisstart.SetX(start.X()+xshift);
      auto thisend = end; 
      thisend.SetX(end.X()+xshift);

      // repeat SCE correction for endpoints
      if (fSCE->EnableCalSpatialSCE() && fSCEposCorr) {
	geo::Point_t temppt = {thisstart.X(),thisstart.Y(),thisstart.Z()};
	geo::TPCID tpcid = fGeometryService->PositionToTPCID(temppt);
	geo::Vector_t  fPosOffsets = fSCE->GetCalPosOffsets(temppt,tpcid.TPC);
	thisstart[0] += fPosOffsets.X();
	thisstart[1] += fPosOffsets.Y();
	thisstart[2] += fPosOffsets.Z();
	temppt.SetX(thisend.X());
	temppt.SetY(thisend.Y());
	temppt.SetZ(thisend.Z());
	tpcid = fGeometryService->PositionToTPCID(temppt);
	fPosOffsets = fSCE->GetCalPosOffsets(temppt,tpcid.TPC);
	thisend[0] += fPosOffsets.X();
	thisend[1] += fPosOffsets.Y();
	thisend[2] += fPosOffsets.Z();
      }

	this_candidate.tpc_track_start.SetXYZ(thisstart.X(),thisstart.Y(),thisstart.Z());
	this_candidate.tpc_track_end.SetXYZ(thisend.X(),thisend.Y(),thisend.Z());

//      match_geometry newmc = makeNULLmc_evtdisp();
      if (startDist<fDistanceLimit || endDist<fDistanceLimit) {
	double distS = (crtPoint-thisstart).Mag();
	double distE =  (crtPoint-thisend).Mag();

	if (distS < distE){ 
	  this_candidate.thishit = crtHit;
	  this_candidate.t0= crtTime;
	  this_candidate.dca = startDist;
	  this_candidate.extrapLen = distS;
//	  t0Candidates.push_back(this_candidate);
	}
	else{
	  this_candidate.thishit = crtHit;
	  this_candidate.t0= crtTime;
	  this_candidate.dca = endDist;
	  this_candidate.extrapLen = distE;
//	  t0Candidates.push_back(this_candidate);
	}
	t0Candidates.push_back(this_candidate);
      }
    }


    //  std::cout << " found " << t0Candidates.size() << " candidates" << std::endl;
    icarus::match_geometry bestmatch = makeNULLmc_evtdisp();
    icarus::match_geometry thismatch = makeNULLmc_evtdisp();
    double min_DCA = DBL_MAX;
    int best_candidate; int start0_or_end1; 
    if(t0Candidates.size() > 0){
      // Find candidate with shortest DCA value
	for(int i=0; i <(int)t0Candidates.size(); i++){

      	  thismatch=t0Candidates[i]; 

	  if (thismatch.simpleDCA_startDir<0 || thismatch.simpleDCA_endDir<0 ){
		if(thismatch.simpleDCA_startDir < thismatch.simpleDCA_endDir && thismatch.simpleDCA_startDir < min_DCA) {
			best_candidate = i; start0_or_end1 = 0;
			min_DCA = thismatch.simpleDCA_startDir;
			bestmatch = thismatch;
		}//end if startDCA<endDCA
		else if(thismatch.simpleDCA_startDir > thismatch.simpleDCA_endDir && thismatch.simpleDCA_endDir < min_DCA){
			best_candidate = i; start0_or_end1 = 1;
			min_DCA = thismatch.simpleDCA_endDir;
			bestmatch = thismatch;
		}//end if startDCA>endDCA
		else if(thismatch.simpleDCA_startDir == thismatch.simpleDCA_endDir && thismatch.simpleDCA_startDir < min_DCA && thismatch.simpleDCA_endDir < min_DCA){
			best_candidate = i; start0_or_end1 = 2;
			min_DCA = thismatch.simpleDCA_endDir;
			bestmatch = thismatch;
		}//end if startDCA==endDCA
	  }//end check for either to be <0  
	else if ((thismatch.simpleDCA_startDir<min_DCA && thismatch.simpleDCA_startDir >=0) || (thismatch.simpleDCA_endDir<min_DCA && thismatch.simpleDCA_endDir>=0)){
		if(thismatch.simpleDCA_startDir < thismatch.simpleDCA_endDir){
			best_candidate = i; start0_or_end1 = 0;
			min_DCA = thismatch.simpleDCA_startDir;
			bestmatch = thismatch;
		}//end if startDCA<endDCA
		else if(thismatch.simpleDCA_endDir < thismatch.simpleDCA_startDir){
			best_candidate = i; start0_or_end1 = 1;
			min_DCA = thismatch.simpleDCA_endDir;
			bestmatch = thismatch;
		}//end if endDCA<startDCA
		else if(thismatch.simpleDCA_startDir==thismatch.simpleDCA_endDir){
			best_candidate = i; start0_or_end1 = 2;	
			min_DCA = thismatch.simpleDCA_endDir;
			bestmatch = thismatch;
		}//end if startDCA==endDCA
	}//end check for either to be <min_DCA but >0
      }//end loop over candidates

	if(start0_or_end1 == 0){

		t0Candidates[best_candidate].is_best_DCA_startDir = true;

	}//end if(start0_or_end1==0)
	else if(start0_or_end1 == 1){

		t0Candidates[best_candidate].is_best_DCA_endDir = true;

	}//end if(start0_or_end1==1)
	else if(start0_or_end1 == 2){

		t0Candidates[best_candidate].is_best_DCA_startDir = true;
		t0Candidates[best_candidate].is_best_DCA_endDir = true;
	
	}//end if(start0_or_end1==2)

    }//end if(t0Candidates.size()>0)


    /*
    f(t0Candidates.size() > 0){
      // Find candidate with shortest DCA or DCA/L value
      bestmatch=t0Candidates[0];
      double sin_angle = bestmatch.dca/bestmatch.extrapLen;
      if (fDCAoverLength) { // Use dca/extrapLen to judge best
	for(auto &thisCand : t0Candidates){
	  double this_sin_angle = thisCand.dca/thisCand.extrapLen;
	  if (bestmatch.dca<0 )bestmatch=thisCand;
	  else if (this_sin_angle<sin_angle && thisCand.dca>=0)bestmatch=thisCand;
	}
      }
      else { // use Dca to judge best
	for(auto &thisCand : t0Candidates){
	  if (bestmatch.dca<0 )bestmatch=thisCand;
	  else if (thisCand.dca<bestmatch.dca && thisCand.dca>=0)bestmatch=thisCand;
	}
      }
    }*/

    //  std::cout << "best match has dca of " << bestmatch.dca << std::endl;
//    return bestmatch;
	return t0Candidates;
  }

  // Simple distance of closest approach between infinite track and centre of hit
  double CRTT0MatchAlg_evtdisp::SimpleDCA(sbn::crt::CRTHit hit, TVector3 start, TVector3 direction){

    TVector3 pos (hit.x_pos, hit.y_pos, hit.z_pos);
    TVector3 end = start + direction;
    double denominator = direction.Mag();
    double numerator = (pos - start).Cross(pos - end).Mag();
    return numerator/denominator;

  }

  // Minimum distance from infinite track to CRT hit assuming that hit is a 2D square
  double CRTT0MatchAlg_evtdisp::DistToCrtHit(sbn::crt::CRTHit hit, TVector3 start, TVector3 end){

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
  double CRTT0MatchAlg_evtdisp::LineSegmentDistance(TVector3 start1, TVector3 end1, TVector3 start2, TVector3 end2){

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

  bool IsCathodeCrosser(const art::Event& event, recob::Track track, double &timestamp){

//  auto this_tpcTrackHandle = event.getValidHandle<std::vector<recob::Track>>(track);

	/*PFParticles*/
/*	art::Handle< std::vector<recob::PFParticle> > pfpsHandleW;
  	std::vector< art::Ptr<recob::PFParticle> > pfpsW;
	  if ( event.getByLabel("pandoraGausCryoW",pfpsHandleW) ) {
	     art::fill_ptr_vector(pfpsW,pfpsHandleW); 
	  }	
	  else{
//	    mf::LogWarning("TrkShwDuality") << "Event failed to find recob::PFParticle of a label.";
	    return false;
	  }
	
	art::Handle< std::vector<recob::PFParticle> > pfpsHandleE;
	std::vector< art::Ptr<recob::PFParticle> > pfpsE;
	  if ( event.getByLabel("pandoraGausCryoE",pfpsHandleE) ) {
	    art::fill_ptr_vector(pfpsE,pfpsHandleE);
	  }
	  else{
//	    mf::LogWarning("TrkShwDuality") << "Event failed to find recob::PFParticle of a label.";
	    return false;
	  }
*/	
	/*Tracks*/
/*	art::Handle< std::vector<recob::Track> > trksHandleW;
	std::vector< art::Ptr<recob::Track> > trksW;
	  if ( event.getByLabel("pandoraTrackGausCryoW",trksHandleW) ) {
	    art::fill_ptr_vector(trksW,trksHandleW);
	  }
	  else{
//	    mf::LogWarning("TrkShwDuality") << "Event failed to find recob::Track of a label.";
	    return false;
	  }
	
	art::Handle< std::vector<recob::Track> > trksHandleE;
	std::vector< art::Ptr<recob::Track> > trksE;
	  if ( event.getByLabel("pandoraTrackGausCryoE",trksHandleE) ) {
	    art::fill_ptr_vector(trksE,trksHandleE);
	  }
	  else{
//	    mf::LogWarning("TrkShwDuality") << "Event failed to find recob::Track of a label.";
	    return false;
	  }


	//T0 <-> PFParticle
	art::FindManyP<anab::T0> fmt0W(pfpsHandleW, event, "pandoraGausCryoW");
  	  if( !fmt0W.isValid() ){
//   	    mf::LogWarning("TrkShwDuality") << "Event failed to find many for anab::T0 <-> recob::PFParticle (West).";
	    return false;
  	  }

    	  art::FindManyP<anab::T0> fmt0E(pfpsHandleE, event, "pandoraGausCryoE");
  	  if( !fmt0E.isValid() ){
//    	    mf::LogWarning("TrkShwDuality") << "Event failed to find many for anab::T0 <-> recob::PFParticle (East).";
    	    return false;
  	  }

 	//Track <-> PFParticle
	  art::FindManyP<recob::Track> fmtrkW(pfpsHandleW, event, "pandoraTrackGausCryoW");
  	  if( !fmtrkW.isValid() ){
//    	    mf::LogWarning("TrkShwDuality") << "Event failed to find many for recob::Track <-> recob::PFParticle (West).";
    	    return false;
  	  }

  	  art::FindManyP<recob::Track> fmtrkE(pfpsHandleE, event, "pandoraTrackGausCryoE");
  	  if( !fmtrkE.isValid() ){
//    	    mf::LogWarning("TrkShwDuality") << "Event failed to find many for recob::Track <-> recob::PFParticle (East).";
    	    return false;
  	  }
*/
return false;

  }//end definition of IsCathodeCrosser

  // Intersection between axis-aligned cube and infinite line
  // (https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-box-intersection)
  std::pair<TVector3, TVector3> CRTT0MatchAlg_evtdisp::CubeIntersection(TVector3 min, TVector3 max, TVector3 start, TVector3 end){

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

    // Find final intersection points if(tzmin > tmin) tmin = tzmin;
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

}
