////////////////////////////////////////////////////////////////////////
// Class:       TrackCaloSkimmer
// Plugin Type: analyzer (art v3_06_03)
// File:        TrackCaloSkimmer_module.cc
//
// Generated at Mon May 17 09:46:34 2021 by Gray Putnam using cetskelgen
// from cetlib version v3_11_01.
////////////////////////////////////////////////////////////////////////

#include "TrackCaloSkimmer.h"

#include "art/Utilities/make_tool.h"

// Global functions / data for fitting
const size_t MAX_N_FIT_DATA = 30;

static int N_FIT_DATA = -1;
static double FIT_RR[MAX_N_FIT_DATA];
static double FIT_DQDX[MAX_N_FIT_DATA];

void ConstResiduals(int &npar, double *g, double &result, double *par, int flag) {
  double ret = 0;

  double C = *par;

  for (int i = 0; i < N_FIT_DATA; i++) {
    double diff = FIT_DQDX[i] - C;
    ret += diff*diff;
  }

  result = sqrt(ret);
}

void ExpResiduals(int &npar, double *g, double &result, double *par, int flag) {
  double ret = 0;

  double A = par[0];
  double R = par[1];

  for (int i = 0; i < N_FIT_DATA; i++) {
    double diff = FIT_DQDX[i] - A*exp(-FIT_RR[i]/R);
    ret += diff*diff;
  }

  result = sqrt(ret);
}

sbn::TrackCaloSkimmer::~TrackCaloSkimmer() {
  delete fTrack;
}

sbn::TrackCaloSkimmer::TrackCaloSkimmer(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    fFitExp(2),
    fFitConst(1)
{
  // Grab config
  fPFPproducer  = p.get< art::InputTag > ("PFPproducer","pandoraGausCryo0");
  fT0Producer   = p.get< art::InputTag > ("T0producer", "pandoraGausCryo0");
  fCALOproducer = p.get< art::InputTag > ("CALOproducer");
  fTRKproducer  = p.get< art::InputTag > ("TRKproducer" );
  fRequireT0 = p.get<bool>("RequireT0", false);

  std::vector<fhicl::ParameterSet> selection_tool_configs(p.get<std::vector<fhicl::ParameterSet>>("SelectionTools"), {});
  for (const fhicl::ParameterSet &p: selection_tool_configs) {
    fSelectionTools.push_back(art::make_tool<sbn::ITCSSelectionTool>(p));
  }

  // Setup meta info
  fMeta.iproc = -1;
  fMeta.ifile = -1;
  const char *process_str = std::getenv("PROCESS");
  if (process_str) {
    try {
      fMeta.iproc = std::stoi(process_str);
    }
    catch (...) {}
  }

  // Output stuff
  fTrack = new sbn::TrackInfo();

  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("TrackCaloSkim", "Calo Tree");
  fTree->Branch("trk", &fTrack);
}

void sbn::TrackCaloSkimmer::analyze(art::Event const& e)
{
  unsigned evt = e.event();
  unsigned sub = e.subRun();
  unsigned run = e.run();
  std::cout << "[TrackCaloSkimmer::analyzeEvent] Run: " << run << ", SubRun: " << sub << ", Event: "<< evt << ", Is Data: " << e.isRealData() << std::endl;

  fMeta.evt = evt;
  fMeta.subrun = sub;
  fMeta.run = run;

  std::vector<art::Ptr<recob::PFParticle>> PFParticleList;
  try {
    art::ValidHandle<std::vector<recob::PFParticle>> pfparticles = e.getValidHandle<std::vector<recob::PFParticle>>(fPFPproducer);
    art::fill_ptr_vector(PFParticleList, pfparticles);
  }
  catch(...) {
      std::cout << "PFP's with tag: " << fPFPproducer << " not present.\n";
      // Real data may have missing products -- just ignore the event
      if (e.isRealData()) return;
      else throw;
  }

  art::FindManyP<anab::T0> fmT0(PFParticleList, e, fT0Producer);

  // Now we don't need to guard access to further data. If this is an empty event it should be caught by PFP's or Hit's
  art::ValidHandle<std::vector<recob::Track>> tracks = e.getValidHandle<std::vector<recob::Track>>(fTRKproducer); 

  // Track - associated data
  art::FindManyP<recob::Track> fmTracks(PFParticleList, e, fTRKproducer);
  art::FindManyP<recob::Hit, recob::TrackHitMeta> fmtrkHits(tracks, e, fTRKproducer);
  art::FindManyP<anab::Calorimetry> fmCalo(tracks, e, fCALOproducer);

  // service data

  // Build global track info
  std::vector<GlobalTrackInfo> track_infos;
  for (const recob::Track &t: *tracks) {
    track_infos.push_back({
      t.Start(), t.End(), t.StartDirection(), t.EndDirection()
    });
  }

  for (art::Ptr<recob::PFParticle> p_pfp: PFParticleList) {
    const recob::PFParticle &pfp = *p_pfp;

    const std::vector<art::Ptr<recob::Track>> thisTrack = fmTracks.at(p_pfp.key());
    if (thisTrack.size() != 1)
      continue;

    art::Ptr<recob::Track> trkPtr = thisTrack.at(0);

    std::vector<art::Ptr<anab::Calorimetry>> emptyCaloVector;
    const std::vector<art::Ptr<anab::Calorimetry>> &calo = fmCalo.isValid() ? fmCalo.at(trkPtr.key()) : emptyCaloVector;

    std::vector<art::Ptr<recob::Hit>> emptyHitVector;
    const std::vector<art::Ptr<recob::Hit>> &trkHits  = fmtrkHits.isValid() ? fmtrkHits.at(trkPtr.key()) : emptyHitVector;

    std::vector<const recob::TrackHitMeta*> emptyTHMVector;
    const std::vector<const recob::TrackHitMeta*> &trkHitMetas = fmtrkHits.isValid() ? fmtrkHits.data(trkPtr.key()) : emptyTHMVector;

    float t0 = std::numeric_limits<float>::signaling_NaN();
    if (fmT0.isValid() && fmT0.at(p_pfp.key()).size()) t0 = fmT0.at(p_pfp.key()).at(0)->Time();

    if (fRequireT0 && fmT0.at(p_pfp.key()).size() == 0) {
      continue;
    }

    std::cout << "Processing new track! ID: " << trkPtr->ID() << " time: " << t0 << std::endl;

    FillTrack(*trkPtr, pfp, t0, trkHits, trkHitMetas, calo, track_infos);
  }
}

// helpers
float HitMinTime(const std::vector<sbn::HitInfo> &hits, bool TPCE) {
  double min = -1;

  for (const sbn::HitInfo &h: hits) {
    bool hit_is_TPCE = h.tpc <= 1;
    if (h.oncalo && hit_is_TPCE == TPCE) {
      if (min < 0. || h.time < min) min = h.time;
    } 
  }

  return min;
}

float HitMaxTime(const std::vector<sbn::HitInfo> &hits, bool TPCE) {
  double max = -1;

  for (const sbn::HitInfo &h: hits) {
    bool hit_is_TPCE = h.tpc <= 1;
    if (h.oncalo && hit_is_TPCE == TPCE) {
      if (max < 0. || h.time > max) max = h.time;
    } 
  }

  return max;
}

void sbn::TrackCaloSkimmer::FillTrack(const recob::Track &track, 
    const recob::PFParticle &pfp, float t0, 
    const std::vector<art::Ptr<recob::Hit>> &hits,
    const std::vector<const recob::TrackHitMeta*> &thms,
    const std::vector<art::Ptr<anab::Calorimetry>> &calo,
    const std::vector<GlobalTrackInfo> &tracks) {

  // Reset the track object
  *fTrack = sbn::TrackInfo();

  // Fill top level stuff
  fTrack->meta = fMeta;
  fTrack->t0 = t0;
  fTrack->id = track.ID();
  fTrack->clear_cosmic_muon = pfp.Parent() == recob::PFParticle::kPFParticlePrimary;
  fTrack->ndaughters = pfp.Daughters().size();

  fTrack->length = track.Length();
  fTrack->start_x = track.Start().X();
  fTrack->start_y = track.Start().Y();
  fTrack->start_z = track.Start().Z();
  fTrack->end_x = track.End().X();
  fTrack->end_y = track.End().Y();
  fTrack->end_z = track.End().Z();
  fTrack->dir_x = track.StartDirection().X();
  fTrack->dir_y = track.StartDirection().Y();
  fTrack->dir_z = track.StartDirection().Z();

  // Fill each hit
  for (unsigned i_hit = 0; i_hit < hits.size(); i_hit++) {
    sbn::HitInfo hinfo = MakeHit(*hits[i_hit], hits[i_hit].key(), *thms[i_hit], track, calo);
    if (hinfo.plane == 0) {
      fTrack->hits0.push_back(hinfo);
    }
    else if (hinfo.plane == 1) {
      fTrack->hits1.push_back(hinfo);
    } 
    else if (hinfo.plane == 2) {
      fTrack->hits2.push_back(hinfo);
    }
  }

  // Hit summary info
  fTrack->hit_min_time_p0_tpcE = HitMinTime(fTrack->hits0, true);
  fTrack->hit_max_time_p0_tpcE = HitMaxTime(fTrack->hits0, true);
  fTrack->hit_min_time_p0_tpcW = HitMinTime(fTrack->hits0, false);
  fTrack->hit_max_time_p0_tpcW = HitMaxTime(fTrack->hits0, false);
  fTrack->hit_min_time_p1_tpcE = HitMinTime(fTrack->hits1, true);
  fTrack->hit_max_time_p1_tpcE = HitMaxTime(fTrack->hits1, true);
  fTrack->hit_min_time_p1_tpcW = HitMinTime(fTrack->hits1, false);
  fTrack->hit_max_time_p1_tpcW = HitMaxTime(fTrack->hits1, false);
  fTrack->hit_min_time_p2_tpcE = HitMinTime(fTrack->hits2, true);
  fTrack->hit_max_time_p2_tpcE = HitMaxTime(fTrack->hits2, true);
  fTrack->hit_min_time_p2_tpcW = HitMinTime(fTrack->hits2, false);
  fTrack->hit_max_time_p2_tpcW = HitMaxTime(fTrack->hits2, false);

  // Try fitting the constant and exponentials to the tail of dQ/dx v. RR on the collection plane
  std::vector<double> fit_rr;
  std::vector<double> fit_dqdx;

  for (const HitInfo &h: fTrack->hits2) {
    if (h.oncalo && h.rr > 0. && h.rr < 5.) {
      fit_rr.push_back(h.rr);
      fit_dqdx.push_back(h.dqdx);
    }
  }

  // Copy the fit data to the global array
  for (unsigned i = 0; i < fit_rr.size() && i < MAX_N_FIT_DATA; i++) {
    FIT_RR[i] = fit_rr[i];
    FIT_DQDX[i] = fit_dqdx[i];

    std::cout << "FIT: " << FIT_RR[i] << " " << FIT_DQDX[i] << std::endl;
  }
  N_FIT_DATA = std::min(fit_rr.size(), MAX_N_FIT_DATA);

  fTrack->n_fit_point = N_FIT_DATA;
  if (fTrack->n_fit_point > 2) { // need more points than params
    // Fit the Exponential
    fFitExp.SetFCN(ExpResiduals);
    fFitExp.SetParameter(0, "A", *std::max_element(fit_dqdx.begin(), fit_dqdx.end()), 200, 0, 5000);
    fFitExp.SetParameter(1, "R", 10., 0.5, 0, 1000);
    fFitExp.ExecuteCommand("MIGRAD", 0, 0);

    double A = fFitExp.GetParameter(0);
    double R = fFitExp.GetParameter(1);

    int nparam;
    double param[2] {A, R};
    double residuals = -1;
    ExpResiduals(nparam, NULL, residuals, param, 0);

    fTrack->exp_fit_A = A; 
    fTrack->exp_fit_R = R;
    fTrack->exp_fit_residuals = residuals;

    // Fit the Constant
    fFitConst.SetFCN(ConstResiduals);
    fFitConst.SetParameter(0, "C", std::accumulate(fit_dqdx.begin(), fit_dqdx.end(), 0.), 200, 0, 5000);
    fFitConst.ExecuteCommand("MIGRAD", 0, 0);

    double C = fFitConst.GetParameter(0);

    double cresiduals = -1;
    ConstResiduals(nparam, NULL, cresiduals, &C, 0);

    fTrack->const_fit_C = C;    
    fTrack->const_fit_residuals = cresiduals;

    std::cout << "EXP FIT. A: " << A << " R: " << R << " resid: " << residuals << std::endl;
    std::cout << "CONST FIT. C: " << C << " resid: " << cresiduals << std::endl;
  }

  // get information on nearby tracks
  for (const GlobalTrackInfo &othr: tracks) {
    if ((track.End() - othr.start).r() < 50. || (track.End() - othr.end).r() < 50.) {
      fTrack->tracks_near_end_dist.push_back(std::min((track.End() - othr.start).r(), (track.End() - othr.end).r())); 
      fTrack->tracks_near_end_costh.push_back(
        (track.End() - othr.start).r() < (track.End() - othr.end).r() ? 
          track.EndDirection().Dot(othr.dir) : track.EndDirection().Dot(othr.enddir)); 
    }
  }

  // Save?
  bool select = false;
  if (!fSelectionTools.size()) select = true;

  // Take the OR of each selection tool
  for (const std::unique_ptr<sbn::ITCSSelectionTool> &t: fSelectionTools) {
    if (t->Select(*fTrack)) {
      select = true;
      break;
    }
  }
  
  // Save!
  if (select) {
    std::cout << "Track Selected!\n";
    fTree->Fill();
  }
}

sbn::HitInfo sbn::TrackCaloSkimmer::MakeHit(const recob::Hit &hit,
    unsigned hkey,
    const recob::TrackHitMeta &thm,
    const recob::Track &trk,
    const std::vector<art::Ptr<anab::Calorimetry>> &calo) {
  // HitInfo to save
  sbn::HitInfo hinfo;

  // information from the hit object
  hinfo.integral = hit.Integral();
  hinfo.sumadc = hit.SummedADC();
  hinfo.width = hit.RMS();
  hinfo.time = hit.PeakTime();
  hinfo.mult = hit.Multiplicity();
  hinfo.wire = hit.WireID().Wire;
  hinfo.plane = hit.WireID().Plane;
  hinfo.tpc = hit.WireID().TPC;

  // Information from the TrackHitMeta
  bool badhit = (thm.Index() == std::numeric_limits<unsigned int>::max()) ||
                    (!trk.HasValidPoint(thm.Index()));

  hinfo.ontraj = !badhit;

  // Save trajectory information if we can
  if (!badhit) {
    geo::Point_t loc = trk.LocationAtPoint(thm.Index());
    hinfo.x = loc.X();
    hinfo.y = loc.Y();
    hinfo.z = loc.Z();

    // And determine if the Hit is on a Calorimetry object
    for (const art::Ptr<anab::Calorimetry> &c: calo) {
      if (c->PlaneID().Plane != hinfo.plane) continue;

      // Found the plane! Now find the hit:
      for (unsigned i_calo = 0; i_calo < c->dQdx().size(); i_calo++) {
        if (c->TpIndices()[i_calo] == hkey) { // "TpIndices" match to the hit key
          // Fill the calo information associated with the hit
          hinfo.oncalo = true;
          hinfo.pitch = c->TrkPitchVec()[i_calo];
          hinfo.dqdx = c->dQdx()[i_calo];
          hinfo.rr = c->ResidualRange()[i_calo];
          break;
        }
      }
      break;
    }
  }

  return hinfo;
}

DEFINE_ART_MODULE(sbn::TrackCaloSkimmer)
