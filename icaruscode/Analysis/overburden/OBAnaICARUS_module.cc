////////////////////////////////////////////////////////////////////////
// Class:       OBAnaICARUS
// Plugin Type: analyzer (art v3_05_01)
// File:        OBAnaICARUS_module.cc
//
// Generated at Mon Nov 30 14:16:07 2020 by Biswaranjan Behera using cetskelgen
// from cetlib version v3_10_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art_root_io/TFileService.h"

#include <memory>
#include <map>

#include "TTree.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/MCBase/MCShower.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "larcoreobj/SummaryData/POTSummary.h"

#include "larcorealg/CoreUtils/ParticleFilters.h" // util::PositionInVolumeFilter
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()

#include <boost/uuid/uuid.hpp>            // uuid class
#include <boost/uuid/uuid_generators.hpp> // generators
#include <boost/uuid/uuid_io.hpp>         // streaming operators etc.

namespace obana {
  class OBAnaICARUS;
}


class obana::OBAnaICARUS : public art::EDAnalyzer {
public:
  explicit OBAnaICARUS(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  OBAnaICARUS(OBAnaICARUS const&) = delete;
  OBAnaICARUS(OBAnaICARUS&&) = delete;
  OBAnaICARUS& operator=(OBAnaICARUS const&) = delete;
  OBAnaICARUS& operator=(OBAnaICARUS&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;
  void beginSubRun(art::SubRun const& sr) override;

  using Point_t = std::array<double, 3>;

private:

  // Declare member data here.
  /// Clear vectors
  void clear_vectors();

  /// Finds the ancestors that was created in the Ovrburden, if any
  int FindMotherInOverburden(simb::MCParticle);

  /// Check if the point is in the detector
  bool InDetector(const double& x, const double& y, const double& z) const;

  /// Check if an MCP passes through the detector, sets step to the step index when particle is in detector
  bool InDetector(const simb::MCParticle&, int & step);

  /// Saves the pi0 info to a separate tree, give the pi0 track id
  void SavePi0ShowerInfo(int pi0_track_id);

  /// Saves the muon info to a separate tree, give the muon track id
  void SaveMuonShowerInfo(int muon_track_id);


  /// Configures and returns a particle filter
  std::unique_ptr<util::PositionInVolumeFilter> CreateParticleVolumeFilter
  (std::set<std::string> const& vol_names) const;

  std::map<unsigned int, simb::MCParticle> _trackid_to_mcparticle;

  std::unique_ptr<util::PositionInVolumeFilter> _part_filter;


  std::string _mctruth_producer    = "generator"; // For storing POT information
  std::string _mcparticle_producer = "largeant";
  std::string _mctrack_producer    = "mcreco";
  std::string _mcshower_producer   = "mcreco";

  bool _save_pi0_tree = true;
  bool _save_muon_tree = true;
  bool _simulating_dirt = false;

  std::vector<unsigned int> _pi0_ids;
  std::vector<unsigned int> _muon_ids;


  std::vector<std::string> _overburden_volumes = {"volShieldingLid", "volShieldingTop", "volMezzanineLid"};

  double xminc0 = -358.49;
  double xmaxc0 = -61.94;
  double xminc1 = 61.94;
  double xmaxc1 = 358.49;
  double ymin = -181.86;
  double ymax = 134.96;
  double zmin = -894.951;
  double zmax = 894.951;

  double _x_max{std::numeric_limits<double>::min()}; //!< x-max of volume box used to determine whether to save track information
  double _x_min{std::numeric_limits<double>::max()}; //!< x-min of volume box used to determine whether to save track information
  double _y_max{std::numeric_limits<double>::min()}; //!< y-max of volume box used to determine whether to save track information
  double _y_min{std::numeric_limits<double>::max()}; //!< y-min of volume box used to determine whether to save track information
  double _z_max{std::numeric_limits<double>::min()}; //!< z-max of volume box used to determine whether to save track information
  double _z_min{std::numeric_limits<double>::max()}; //!< z-max of volume box used to determine whether to save track information

  boost::uuids::uuid _uuid; ///< A unique ID to identify different events in files with same event number
  std::string _uuid_str; ///< Same as uuid, but converted to string

  TTree* _tree;

  int _run, _subrun, _event;

  float _nu_e;
  int _nu_pdg;
  int _nu_ccnc;
  int _nu_mode;
  int _nu_int_type;
  float _nu_vtx_x;
  float _nu_vtx_y;
  float _nu_vtx_z;
  float _nu_px;
  float _nu_py;
  float _nu_pz;
  int _nu_pip_mult; ///< Pi0 multiplicity
  int _nu_pi0_mult; ///< Pi plus multiplicity
  int _nu_p_mult; ///< Proton multiplicity
  std::vector<int> _pars_pdg; ///< All other particles produced - pdg code
  std::vector<float> _pars_e; ///< All other particles produced - energy


  std::vector<float> _mcp_px;
  std::vector<float> _mcp_py;
  std::vector<float> _mcp_pz;
  std::vector<float> _mcp_e;
  std::vector<float> _mcp_vx;
  std::vector<float> _mcp_vy;
  std::vector<float> _mcp_vz;
  std::vector<float> _mcp_endx;
  std::vector<float> _mcp_endy;
  std::vector<float> _mcp_endz;
  std::vector<float> _mcp_pdg;
  std::vector<float> _mcp_mother;
  std::vector<float> _mcp_status_code;
  std::vector<std::string> _mcp_process;
  std::vector<std::string> _mcp_end_process;
  std::vector<bool> _mcp_intpc;
  std::vector<float> _mcp_intpc_e;
  std::vector<float> _mcp_trackid;
  std::vector<float> _mcp_intpc_nu_e;
  //std::string _neut_par_uuid;
  std::vector<int> _neut_daughters_pdg; ///< All the neutron daughters 
  std::vector<float> _neut_daughters_e; ///< All the neutron daughters energy
  std::vector<float> _neut_daughters_px; ///< All the neutron daughters momentrum in x direction 
  std::vector<float> _neut_daughters_py; ///< All the neutron daughters momentrum in y direction
  std::vector<float> _neut_daughters_pz; ///< All the neutron daughters momentrum in z direction
  std::vector<std::string> _neut_daughters_startprocess; ///< All the neutron daughters start process
  std::vector<std::string> _neut_daughters_endprocess; ////< All the neutron daughters end process
  std::vector<float> _neut_daughters_start_x; ///< All the neutron daughters start x
  std::vector<float> _neut_daughters_start_y; ///< All the neutron daughters start y
  std::vector<float> _neut_daughters_start_z; ///< All the neutron daughters start z
  std::vector<float> _neut_daughters_end_x; ///< All the neutron daughters end x
  std::vector<float> _neut_daughters_end_y; ///< All the neutron daughters end y
  std::vector<float> _neut_daughters_end_z; ///< All the neutron daughters end z
  std::vector<float> _neut_daughters_trackid; ///< All the neutron daughters trackid


  std::vector<int> _neut_grand_daughters_pdg; ///< All the neutron daughters
  std::vector<float> _neut_grand_daughters_e; ///< All the neutrons grand daughters energy
  std::vector<float> _neut_grand_daughters_px; ///< All the neutrons grand daughters momentrum in x direction
  std::vector<float> _neut_grand_daughters_py; ///< All the neutrons grand daughters momentrum in y direction
  std::vector<float> _neut_grand_daughters_pz; ///< All the neutrons grand daughters momentrum in z direction
  std::vector<std::string> _neut_grand_daughters_startprocess; ///< All the neutrons grand daughters start process
  std::vector<std::string> _neut_grand_daughters_endprocess; ////< All the neutrons grand daughters end process
  std::vector<float> _neut_grand_daughters_start_x; ///< All the neutrons grand daughters start x
  std::vector<float> _neut_grand_daughters_start_y; ///< All the neutrons grand daughters start y
  std::vector<float> _neut_grand_daughters_start_z; ///< All the neutrons grand daughters start z
  std::vector<float> _neut_grand_daughters_end_x; ///< All the neutrons grand daughters end x
  std::vector<float> _neut_grand_daughters_end_y; ///< All the neutrons grand daughters end y
  std::vector<float> _neut_grand_daughters_end_z; ///< All the neutrons grand daughters end z
  std::vector<float> _neut_grand_daughters_trackid; ///< All the neutrons grand daughters trackid


  std::vector<float> _mcs_pdg;
  std::vector<std::string> _mcs_process;
  std::vector<float> _mcs_start_x;
  std::vector<float> _mcs_start_y;
  std::vector<float> _mcs_start_z;
  std::vector<float> _mcs_end_x;
  std::vector<float> _mcs_end_y;
  std::vector<float> _mcs_end_z;
  std::vector<float> _mcs_start_px;
  std::vector<float> _mcs_start_py;
  std::vector<float> _mcs_start_pz;
  std::vector<float> _mcs_start_e;
  std::vector<double> _mcs_charge_col;
  std::vector<double> _mcs_charge_ind2;
  std::vector<double> _mcs_charge_ind1;
  std::vector<float> _mcs_mother_pdg;
  std::vector<float> _mcs_mother_trackid;
  std::vector<float> _mcs_mother_start_x;
  std::vector<float> _mcs_mother_start_y;
  std::vector<float> _mcs_mother_start_z;
  std::vector<float> _mcs_mother_start_e;
  std::vector<float> _mcs_mother_end_x;
  std::vector<float> _mcs_mother_end_y;
  std::vector<float> _mcs_mother_end_z;
  std::vector<float> _mcs_mother_end_e;
  std::vector<std::string> _mcs_mother_process;
  std::vector<float> _mcs_ancestor_pdg;
  std::vector<std::string> _mcs_ancestor_process;
  std::vector<float> _mcs_ancestor_start_e;
  std::vector<float> _mcs_ancestor_end_e;
  std::vector<float> _mcs_mother_in_ob_trackid;
  int _n_mcs_lt1; ///< Number of MC showers with energy less than 10 MeV
  int _n_mcs_lt1_from_ob; ///< Number of MC showers with energy less than 10 MeV and coming from OB

  std::vector<float> _mct_pdg;
  std::vector<std::string> _mct_process;
  std::vector<float> _mct_start_x;
  std::vector<float> _mct_start_y;
  std::vector<float> _mct_start_z;
  std::vector<float> _mct_end_x;
  std::vector<float> _mct_end_y;
  std::vector<float> _mct_end_z;
  std::vector<float> _mct_start_px;
  std::vector<float> _mct_start_py;
  std::vector<float> _mct_start_pz;
  std::vector<float> _mct_start_e;
  std::vector<float> _mct_mother_pdg;
  std::vector<std::string> _mct_mother_process;
  std::vector<float> _mct_ancestor_pdg;
  std::vector<std::string> _mct_ancestor_process;
  std::vector<float> _mct_mother_in_ob_trackid;
  int _n_mct_lt1; ///< Number of MC tracks with energy less than 10 MeV
  int _n_mct_lt1_from_ob; ///< Number of MC tracks with energy less than 10 MeV and coming from OB

  TTree* _pi0_tree; ///< A pi0 TTree, one entry per pi0

  float _pi0_par_e; ///< The pi0 energy
  float _pi0_par_start_x; ///< The pi0 start x
  float _pi0_par_start_y; ///< The pi0 start y
  float _pi0_par_start_z; ///< The pi0 start z
  float _pi0_par_end_x; ///< The pi0 end x
  float _pi0_par_end_y; ///< The pi0 end y
  float _pi0_par_end_z; ///< The pi0 end z
  int _pi0_par_mother_pdg; ///< The pi0 mother pdg
  float _pi0_par_mother_e; ///< The pi0 mother energy
  float _pi0_par_mother_end_e; ///< The pi0 mother energy at end
  int _pi0_par_ancestor_trackid; ///< The pi0 primary particle ancestor track id (allows easy pi0 clustering by cosmic interaction)
  std::string _pi0_par_ancestor_uuid; ///< The pi0 unique ID for the event

  std::vector<int> _pi0_daughters_pdg; ///< All the pi0 daughters (usually two photons) (pdg)
  std::vector<float> _pi0_daughters_e; ///< All the pi0 daughters (usually two photons) (energy)
  std::vector<std::string> _pi0_daughters_startprocess; ///< All the pi0 daughters (usually two photons) (energy) (start process)
  std::vector<std::string> _pi0_daughters_endprocess; ////< All the pi0 daughters (usually two photons) (energy) (end process)
  std::vector<float> _pi0_daughters_start_x; ///< All the pi0 daughters (usually two photons) (start x)
  std::vector<float> _pi0_daughters_start_y; ///< All the pi0 daughters (usually two photons) (start y)
  std::vector<float> _pi0_daughters_start_z; ///< All the pi0 daughters (usually two photons) (start z)
  std::vector<float> _pi0_daughters_end_x; ///< All the pi0 daughters (usually two photons) (end x)
  std::vector<float> _pi0_daughters_end_y; ///< All the pi0 daughters (usually two photons) (end y)
  std::vector<float> _pi0_daughters_end_z; ///< All the pi0 daughters (usually two photons) (end z)

  std::vector<int> _pi0_event_particles_pdg; ///< All the particles produced together with the pi0
  std::vector<float> _pi0_event_particles_e; ///< All the particles produced together with the pi0

  std::vector<int> _pi0_genealogy_pdg; ///< The full pi0 genealogy, from its mother all the way to the ancestor (pdg)
  std::vector<std::string> _pi0_genealogy_startprocess; ///< The full pi0 genealogy, from its mother all the way to the ancestor (start process)
  std::vector<std::string> _pi0_genealogy_endprocess; ///< The full pi0 genealogy, from its mother all the way to the ancestor (end process)
  std::vector<int> _pi0_genealogy_mother; ///< The full pi0 genealogy, from its mother all the way to the ancestor (mother track id)
  std::vector<float> _pi0_genealogy_e; ///< The full pi0 genealogy, from its mother all the way to the ancestor (energy)
  std::vector<float> _pi0_genealogy_start_x; ///< The full pi0 genealogy, from its mother all the way to the ancestor (start x)
  std::vector<float> _pi0_genealogy_start_y; ///< The full pi0 genealogy, from its mother all the way to the ancestor (start y)
  std::vector<float> _pi0_genealogy_start_z; ///< The full pi0 genealogy, from its mother all the way to the ancestor (start z)
  std::vector<float> _pi0_genealogy_end_x; ///< The full pi0 genealogy, from its mother all the way to the ancestor (end x)
  std::vector<float> _pi0_genealogy_end_y; ///< The full pi0 genealogy, from its mother all the way to the ancestor (end y)
  std::vector<float> _pi0_genealogy_end_z; ///< The full pi0 genealogy, from its mother all the way to the ancestor (end z)
  std::vector<float> _pi0_genealogy_trackid; ///< The full pi0 genealogy, from its mother all the way to the ancestor (trackid)

  TTree* _muon_tree; ///< A muon TTree, one entry per muon

  float _muon_par_e; ///< The muon energy
  float _muon_par_start_x; ///< The muon start x
  float _muon_par_start_y; ///< The muon start y
  float _muon_par_start_z; ///< The muon start z
  float _muon_par_end_x; ///< The muon end x
  float _muon_par_end_y; ///< The muon end y
  float _muon_par_end_z; ///< The muon end z
  int _muon_par_mother_pdg; ///< The muon mother pdg
  float _muon_par_mother_e; ///< The muon mother energy
  float _muon_par_mother_end_e; ///< The muon mother energy at end
  int _muon_par_ancestor_trackid; ///< The muon primary particle ancestor track id (allows easy muon clustering by cosmic interaction)
  std::string _muon_par_ancestor_uuid; ///< The muon unique ID for the event

  std::vector<int> _muon_daughters_pdg; ///< All the muon daughters (usually two photons) (pdg)
  std::vector<float> _muon_daughters_e; ///< All the muon daughters (usually two photons) (energy)
  std::vector<std::string> _muon_daughters_startprocess; ///< All the muon daughters (usually two photons) (energy) (start process)
  std::vector<std::string> _muon_daughters_endprocess; ////< All the muon daughters (usually two photons) (energy) (end process)
  std::vector<float> _muon_daughters_start_x; ///< All the muon daughters (usually two photons) (start x)
  std::vector<float> _muon_daughters_start_y; ///< All the muon daughters (usually two photons) (start y)
  std::vector<float> _muon_daughters_start_z; ///< All the muon daughters (usually two photons) (start z)
  std::vector<float> _muon_daughters_end_x; ///< All the muon daughters (usually two photons) (end x)
  std::vector<float> _muon_daughters_end_y; ///< All the muon daughters (usually two photons) (end y)
  std::vector<float> _muon_daughters_end_z; ///< All the muon daughters (usually two photons) (end z)

  std::vector<int> _muon_event_particles_pdg; ///< All the particles produced together with the muon
  std::vector<float> _muon_event_particles_e; ///< All the particles produced together with the muon

  std::vector<int> _muon_genealogy_pdg; ///< The full muon genealogy, from its mother all the way to the ancestor (pdg)
  std::vector<std::string> _muon_genealogy_startprocess; ///< The full muon genealogy, from its mother all the way to the ancestor (start process)
  std::vector<std::string> _muon_genealogy_endprocess; ///< The full muon genealogy, from its mother all the way to the ancestor (end process)
  std::vector<int> _muon_genealogy_mother; ///< The full muon genealogy, from its mother all the way to the ancestor (mother track id)
  std::vector<float> _muon_genealogy_e; ///< The full muon genealogy, from its mother all the way to the ancestor (energy)
  std::vector<float> _muon_genealogy_start_x; ///< The full muon genealogy, from its mother all the way to the ancestor (start x)
  std::vector<float> _muon_genealogy_start_y; ///< The full muon genealogy, from its mother all the way to the ancestor (start y)
  std::vector<float> _muon_genealogy_start_z; ///< The full muon genealogy, from its mother all the way to the ancestor (start z)
  std::vector<float> _muon_genealogy_end_x; ///< The full muon genealogy, from its mother all the way to the ancestor (end x)
  std::vector<float> _muon_genealogy_end_y; ///< The full muon genealogy, from its mother all the way to the ancestor (end y)
  std::vector<float> _muon_genealogy_end_z; ///< The full muon genealogy, from its mother all the way to the ancestor (end z)
  std::vector<float> _muon_genealogy_trackid; ///< The full muon genealogy, from its mother all the way to the ancestor (trackid)



  TTree* _sr_tree; ///< TTree filled per subrun
  int _sr_run; ///< Subrun Run number
  int _sr_subrun; ///< Subrun Subrun number
  double _sr_begintime; ///< Subrun start time
  double _sr_endtime; ///< Subrun end time
  double _sr_pot; ///< Subrun POT
};


obana::OBAnaICARUS::OBAnaICARUS(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  _save_pi0_tree = p.get<bool>("SavePi0Tree", true);
  _simulating_dirt = p.get<bool>("SimulatingDirt", false);

  art::ServiceHandle<art::TFileService> fs;
  _tree = fs->make<TTree>("OBTree","");

  _tree->Branch("run",     &_run,     "run/I");
  _tree->Branch("subrun",     &_subrun,     "subrun/I");
  _tree->Branch("event",     &_event,     "event/I");


  _tree->Branch("nu_e", &_nu_e, "nu_e/F");
  _tree->Branch("nu_pdg", &_nu_pdg, "nu_pdg/I");
  _tree->Branch("nu_ccnc", &_nu_ccnc, "nu_ccnc/I");
  _tree->Branch("nu_mode", &_nu_mode, "nu_mode/I");
  _tree->Branch("nu_int_type", &_nu_int_type, "nu_int_type/I");
  _tree->Branch("nu_vtx_x", &_nu_vtx_x, "nu_vtx_x/F");
  _tree->Branch("nu_vtx_y", &_nu_vtx_y, "nu_vtx_y/F");
  _tree->Branch("nu_vtx_z", &_nu_vtx_z, "nu_vtx_z/F");
  _tree->Branch("nu_px", &_nu_px, "nu_px/F");
  _tree->Branch("nu_py", &_nu_py, "nu_px/F");
  _tree->Branch("nu_pz", &_nu_pz, "nu_px/F");
  _tree->Branch("nu_pip_mult", &_nu_pip_mult, "nu_pip_mult/I");
  _tree->Branch("nu_pi0_mult", &_nu_pi0_mult, "nu_pi0_mult/I");
  _tree->Branch("nu_p_mult", &_nu_p_mult, "nu_p_mult/I");
  _tree->Branch("pars_pdg", "std::vector<int>", &_pars_pdg);
  _tree->Branch("pars_e", "std::vector<float>", &_pars_e);


  _tree->Branch("mcp_px", "std::vector<float>", &_mcp_px);
  _tree->Branch("mcp_py", "std::vector<float>", &_mcp_py);
  _tree->Branch("mcp_pz", "std::vector<float>", &_mcp_pz);
  _tree->Branch("mcp_e", "std::vector<float>", &_mcp_e);
  _tree->Branch("mcp_vx", "std::vector<float>", &_mcp_vx);
  _tree->Branch("mcp_vy", "std::vector<float>", &_mcp_vy);
  _tree->Branch("mcp_vz", "std::vector<float>", &_mcp_vz);
  _tree->Branch("mcp_endx", "std::vector<float>", &_mcp_endx);
  _tree->Branch("mcp_endy", "std::vector<float>", &_mcp_endy);
  _tree->Branch("mcp_endz", "std::vector<float>", &_mcp_endz);
  _tree->Branch("mcp_pdg", "std::vector<float>", &_mcp_pdg);
  _tree->Branch("mcp_mother", "std::vector<float>", &_mcp_mother);
  _tree->Branch("mcp_status_code", "std::vector<float>", &_mcp_status_code);
  _tree->Branch("mcp_process", "std::vector<std::string>", &_mcp_process);
  _tree->Branch("mcp_end_process", "std::vector<std::string>", &_mcp_end_process);
  _tree->Branch("mcp_intpc", "std::vector<bool>", &_mcp_intpc);
  _tree->Branch("mcp_intpc_e", "std::vector<float>", &_mcp_intpc_e);
  _tree->Branch("mcp_trackid", "std::vector<float>", &_mcp_trackid);
  _tree->Branch("mcp_intpc_nu_e", "std::vector<float>", &_mcp_intpc_nu_e);

  //_tree->Branch("neut_par_uuid", "std::string", &_neut_par_uuid);
  _tree->Branch("neut_daughters_pdg", "std::vector<int>", &_neut_daughters_pdg);
  _tree->Branch("neut_daughters_px", "std::vector<float>", &_neut_daughters_px);
  _tree->Branch("neut_daughters_py", "std::vector<float>", &_neut_daughters_py);
  _tree->Branch("neut_daughters_pz", "std::vector<float>", &_neut_daughters_pz);
  _tree->Branch("neut_daughters_e", "std::vector<float>", &_neut_daughters_e);
  _tree->Branch("neut_daughters_startprocess", "std::vector<std::string>", &_neut_daughters_startprocess);
  _tree->Branch("neut_daughters_endprocess", "std::vector<std::string>", &_neut_daughters_endprocess);
  _tree->Branch("neut_daughters_start_x", "std::vector<float>", &_neut_daughters_start_x);
  _tree->Branch("neut_daughters_start_y", "std::vector<float>", &_neut_daughters_start_y);
  _tree->Branch("neut_daughters_start_z", "std::vector<float>", &_neut_daughters_start_z);
  _tree->Branch("neut_daughters_end_x", "std::vector<float>", &_neut_daughters_end_x);
  _tree->Branch("neut_daughters_end_y", "std::vector<float>", &_neut_daughters_end_y);
  _tree->Branch("neut_daughters_end_z", "std::vector<float>", &_neut_daughters_end_z);
  _tree->Branch("neut_daughters_trackid", "std::vector<float>", &_neut_daughters_trackid);


  _tree->Branch("neut_grand_daughters_pdg", "std::vector<int>", &_neut_grand_daughters_pdg);
  _tree->Branch("neut_grand_daughters_px", "std::vector<float>", &_neut_grand_daughters_px);
  _tree->Branch("neut_grand_daughters_py", "std::vector<float>", &_neut_grand_daughters_py);
  _tree->Branch("neut_grand_daughters_pz", "std::vector<float>", &_neut_grand_daughters_pz);
  _tree->Branch("neut_grand_daughters_e", "std::vector<float>", &_neut_grand_daughters_e);
  _tree->Branch("neut_grand_daughters_startprocess", "std::vector<std::string>", &_neut_grand_daughters_startprocess);
  _tree->Branch("neut_grand_daughters_endprocess", "std::vector<std::string>", &_neut_grand_daughters_endprocess);
  _tree->Branch("neut_grand_daughters_start_x", "std::vector<float>", &_neut_grand_daughters_start_x);
  _tree->Branch("neut_grand_daughters_start_y", "std::vector<float>", &_neut_grand_daughters_start_y);
  _tree->Branch("neut_grand_daughters_start_z", "std::vector<float>", &_neut_grand_daughters_start_z);
  _tree->Branch("neut_grand_daughters_end_x", "std::vector<float>", &_neut_grand_daughters_end_x);
  _tree->Branch("neut_grand_daughters_end_y", "std::vector<float>", &_neut_grand_daughters_end_y);
  _tree->Branch("neut_grand_daughters_end_z", "std::vector<float>", &_neut_grand_daughters_end_z);
  _tree->Branch("neut_grand_daughters_trackid", "std::vector<float>", &_neut_grand_daughters_trackid);

  _tree->Branch("mct_pdg", "std::vector<float>", &_mct_pdg);
  _tree->Branch("mct_process", "std::vector<std::string>", &_mct_process);
  _tree->Branch("mct_start_x", "std::vector<float>", &_mct_start_x);
  _tree->Branch("mct_start_y", "std::vector<float>", &_mct_start_y);
  _tree->Branch("mct_start_z", "std::vector<float>", &_mct_start_z);
  _tree->Branch("mct_end_x", "std::vector<float>", &_mct_end_x);
  _tree->Branch("mct_end_y", "std::vector<float>", &_mct_end_y);
  _tree->Branch("mct_end_z", "std::vector<float>", &_mct_end_z);
  _tree->Branch("mct_start_px", "std::vector<float>", &_mct_start_px);
  _tree->Branch("mct_start_py", "std::vector<float>", &_mct_start_py);
  _tree->Branch("mct_start_pz", "std::vector<float>", &_mct_start_pz);
  _tree->Branch("mct_start_e", "std::vector<float>", &_mct_start_e);
  _tree->Branch("mct_mother_pdg", "std::vector<float>", &_mct_mother_pdg);
  _tree->Branch("mct_ancestor_pdg", "std::vector<float>", &_mct_ancestor_pdg);
  _tree->Branch("mct_mother_process", "std::vector<std::string>>", &_mct_mother_process);
  _tree->Branch("mct_ancestor_process", "std::vector<std::string>>", &_mct_ancestor_process);
  _tree->Branch("mct_mother_in_ob_trackid", "std::vector<float>", &_mct_mother_in_ob_trackid);
  _tree->Branch("mct_n_lt1",     &_n_mct_lt1,     "mct_n_lt1/I");
  _tree->Branch("mct_n_lt1_from_ob",     &_n_mct_lt1_from_ob,     "mct_n_lt1_from_ob/I");

  _tree->Branch("mcs_pdg", "std::vector<float>", &_mcs_pdg);
  _tree->Branch("mcs_process", "std::vector<std::string>", &_mcs_process);
  _tree->Branch("mcs_start_x", "std::vector<float>", &_mcs_start_x);
  _tree->Branch("mcs_start_y", "std::vector<float>", &_mcs_start_y);
  _tree->Branch("mcs_start_z", "std::vector<float>", &_mcs_start_z);
  _tree->Branch("mcs_end_x", "std::vector<float>", &_mcs_end_x);
  _tree->Branch("mcs_end_y", "std::vector<float>", &_mcs_end_y);
  _tree->Branch("mcs_end_z", "std::vector<float>", &_mcs_end_z);
  _tree->Branch("mcs_start_px", "std::vector<float>", &_mcs_start_px);
  _tree->Branch("mcs_start_py", "std::vector<float>", &_mcs_start_py);
  _tree->Branch("mcs_start_pz", "std::vector<float>", &_mcs_start_pz);
  _tree->Branch("mcs_start_e", "std::vector<float>", &_mcs_start_e);
  _tree->Branch("mcs_charge_col", "std::vector<double>", &_mcs_charge_col);
  _tree->Branch("mcs_charge_ind2", "std::vector<double>", &_mcs_charge_ind2);
  _tree->Branch("mcs_charge_ind1", "std::vector<double>", &_mcs_charge_ind1);
  _tree->Branch("mcs_mother_pdg", "std::vector<float>", &_mcs_mother_pdg);
  _tree->Branch("mcs_mother_trackid", "std::vector<float>", &_mcs_mother_trackid);
  _tree->Branch("mcs_mother_start_x", "std::vector<float>", &_mcs_mother_start_x);
  _tree->Branch("mcs_mother_start_y", "std::vector<float>", &_mcs_mother_start_y);
  _tree->Branch("mcs_mother_start_z", "std::vector<float>", &_mcs_mother_start_z);
  _tree->Branch("mcs_mother_start_e", "std::vector<float>", &_mcs_mother_start_e);
  _tree->Branch("mcs_mother_end_x", "std::vector<float>", &_mcs_mother_end_x);
  _tree->Branch("mcs_mother_end_y", "std::vector<float>", &_mcs_mother_end_y);
  _tree->Branch("mcs_mother_end_z", "std::vector<float>", &_mcs_mother_end_z);
  _tree->Branch("mcs_mother_end_e", "std::vector<float>", &_mcs_mother_end_e);
  _tree->Branch("mcs_ancestor_pdg", "std::vector<float>", &_mcs_ancestor_pdg);
  _tree->Branch("mcs_mother_process", "std::vector<std::string>>", &_mcs_mother_process);
  _tree->Branch("mcs_ancestor_process", "std::vector<std::string>>", &_mcs_ancestor_process);
  _tree->Branch("mcs_ancestor_start_e", "std::vector<float>", &_mcs_ancestor_start_e);
  _tree->Branch("mcs_ancestor_end_e", "std::vector<float>", &_mcs_ancestor_end_e);
  _tree->Branch("mcs_mother_in_ob_trackid", "std::vector<float>", &_mcs_mother_in_ob_trackid);
  _tree->Branch("mcs_n_lt1",     &_n_mcs_lt1,     "mcs_n_lt1/I");
  _tree->Branch("mcs_n_lt1_from_ob",     &_n_mcs_lt1_from_ob,     "mcs_n_lt1_from_ob/I");


  if (_save_pi0_tree) {
    _pi0_tree = fs->make<TTree>("Pi0Tree","");

    _pi0_tree->Branch("run", &_run, "run/I");
    _pi0_tree->Branch("subrun", &_subrun, "subrun/I");
    _pi0_tree->Branch("event", &_event, "event/I");

    _pi0_tree->Branch("pi0_par_e", &_pi0_par_e, "pi0_par_e/F");
    _pi0_tree->Branch("pi0_par_start_x", &_pi0_par_start_x, "pi0_par_start_x/F");
    _pi0_tree->Branch("pi0_par_start_y", &_pi0_par_start_y, "pi0_par_start_y/F");
    _pi0_tree->Branch("pi0_par_start_z", &_pi0_par_start_z, "pi0_par_start_z/F");
    _pi0_tree->Branch("pi0_par_end_x", &_pi0_par_end_x, "pi0_par_end_x/F");
    _pi0_tree->Branch("pi0_par_end_y", &_pi0_par_end_y, "pi0_par_end_y/F");
    _pi0_tree->Branch("pi0_par_end_z", &_pi0_par_end_z, "pi0_par_end_z/F");
    _pi0_tree->Branch("pi0_par_mother_pdg", &_pi0_par_mother_pdg, "pi0_par_mother_pdg/I");
    _pi0_tree->Branch("pi0_par_mother_e", &_pi0_par_mother_e, "pi0_par_mother_e/F");
    _pi0_tree->Branch("pi0_par_mother_end_e", &_pi0_par_mother_end_e, "pi0_par_mother_end_e/F");
    _pi0_tree->Branch("pi0_par_ancestor_trackid", &_pi0_par_ancestor_trackid, "pi0_par_ancestor_trackid/I");
    _pi0_tree->Branch("pi0_par_ancestor_uuid", "std::string", &_pi0_par_ancestor_uuid);

    _pi0_tree->Branch("pi0_event_particles_pdg", "std::vector<int>", &_pi0_event_particles_pdg);
    _pi0_tree->Branch("pi0_event_particles_e", "std::vector<float>", &_pi0_event_particles_e);

    _pi0_tree->Branch("pi0_daughters_pdg", "std::vector<int>", &_pi0_daughters_pdg);
    _pi0_tree->Branch("pi0_daughters_e", "std::vector<float>", &_pi0_daughters_e);
    _pi0_tree->Branch("pi0_daughters_startprocess", "std::vector<std::string>", &_pi0_daughters_startprocess);
    _pi0_tree->Branch("pi0_daughters_endprocess", "std::vector<std::string>", &_pi0_daughters_endprocess);
    _pi0_tree->Branch("pi0_daughters_start_x", "std::vector<float>", &_pi0_daughters_start_x);
    _pi0_tree->Branch("pi0_daughters_start_y", "std::vector<float>", &_pi0_daughters_start_y);
    _pi0_tree->Branch("pi0_daughters_start_z", "std::vector<float>", &_pi0_daughters_start_z);
    _pi0_tree->Branch("pi0_daughters_end_x", "std::vector<float>", &_pi0_daughters_end_x);
    _pi0_tree->Branch("pi0_daughters_end_y", "std::vector<float>", &_pi0_daughters_end_y);
    _pi0_tree->Branch("pi0_daughters_end_z", "std::vector<float>", &_pi0_daughters_end_z);

    _pi0_tree->Branch("pi0_genealogy_pdg", "std::vector<int>", &_pi0_genealogy_pdg);
    _pi0_tree->Branch("pi0_genealogy_startprocess", "std::vector<std::string>", &_pi0_genealogy_startprocess);
    _pi0_tree->Branch("pi0_genealogy_endprocess", "std::vector<std::string>", &_pi0_genealogy_endprocess);
    _pi0_tree->Branch("pi0_genealogy_mother", "std::vector<int>", &_pi0_genealogy_mother);
    _pi0_tree->Branch("pi0_genealogy_e", "std::vector<float>", &_pi0_genealogy_e);
    _pi0_tree->Branch("pi0_genealogy_start_x", "std::vector<float>", &_pi0_genealogy_start_x);
    _pi0_tree->Branch("pi0_genealogy_start_y", "std::vector<float>", &_pi0_genealogy_start_y);
    _pi0_tree->Branch("pi0_genealogy_start_z", "std::vector<float>", &_pi0_genealogy_start_z);
    _pi0_tree->Branch("pi0_genealogy_end_x", "std::vector<float>", &_pi0_genealogy_end_x);
    _pi0_tree->Branch("pi0_genealogy_end_y", "std::vector<float>", &_pi0_genealogy_end_y);
    _pi0_tree->Branch("pi0_genealogy_end_z", "std::vector<float>", &_pi0_genealogy_end_z);
    _pi0_tree->Branch("pi0_genealogy_trackid", "std::vector<float>", &_pi0_genealogy_trackid);


    // Add the neutrino information also to this tree so we can better study pi0s
    // when looking at beam neutrinos
    _pi0_tree->Branch("nu_e", &_nu_e, "nu_e/F");
    _pi0_tree->Branch("nu_pdg", &_nu_pdg, "nu_pdg/I");
    _pi0_tree->Branch("nu_ccnc", &_nu_ccnc, "nu_ccnc/I");
    _pi0_tree->Branch("nu_mode", &_nu_mode, "nu_mode/I");
    _pi0_tree->Branch("nu_vtx_x", &_nu_vtx_x, "nu_vtx_x/F");
    _pi0_tree->Branch("nu_vtx_y", &_nu_vtx_y, "nu_vtx_y/F");
    _pi0_tree->Branch("nu_vtx_z", &_nu_vtx_z, "nu_vtx_z/F");
    _pi0_tree->Branch("nu_pip_mult", &_nu_pip_mult, "nu_pip_mult/I");
    _pi0_tree->Branch("nu_pi0_mult", &_nu_pi0_mult, "nu_pi0_mult/I");
    _pi0_tree->Branch("nu_p_mult", &_nu_p_mult, "nu_p_mult/I");
    _pi0_tree->Branch("pars_pdg", "std::vector<int>", &_pars_pdg);
    _pi0_tree->Branch("pars_e", "std::vector<float>", &_pars_e);

  }

 if (_save_muon_tree) {
    _muon_tree = fs->make<TTree>("MuonTree","");

    _muon_tree->Branch("run", &_run, "run/I");
    _muon_tree->Branch("subrun", &_subrun, "subrun/I");
    _muon_tree->Branch("event", &_event, "event/I");

    _muon_tree->Branch("muon_par_e", &_muon_par_e, "muon_par_e/F");
    _muon_tree->Branch("muon_par_start_x", &_muon_par_start_x, "muon_par_start_x/F");
    _muon_tree->Branch("muon_par_start_y", &_muon_par_start_y, "muon_par_start_y/F");
    _muon_tree->Branch("muon_par_start_z", &_muon_par_start_z, "muon_par_start_z/F");
    _muon_tree->Branch("muon_par_end_x", &_muon_par_end_x, "muon_par_end_x/F");
    _muon_tree->Branch("muon_par_end_y", &_muon_par_end_y, "muon_par_end_y/F");
    _muon_tree->Branch("muon_par_end_z", &_muon_par_end_z, "muon_par_end_z/F");
    _muon_tree->Branch("muon_par_mother_pdg", &_muon_par_mother_pdg, "muon_par_mother_pdg/I");
    _muon_tree->Branch("muon_par_mother_e", &_muon_par_mother_e, "muon_par_mother_e/F");
    _muon_tree->Branch("muon_par_mother_end_e", &_muon_par_mother_end_e, "muon_par_mother_end_e/F");
    _muon_tree->Branch("muon_par_ancestor_trackid", &_muon_par_ancestor_trackid, "muon_par_ancestor_trackid/I");
    _muon_tree->Branch("muon_par_ancestor_uuid", "std::string", &_muon_par_ancestor_uuid);

    _muon_tree->Branch("muon_event_particles_pdg", "std::vector<int>", &_muon_event_particles_pdg);
    _muon_tree->Branch("muon_event_particles_e", "std::vector<float>", &_muon_event_particles_e);

    _muon_tree->Branch("muon_daughters_pdg", "std::vector<int>", &_muon_daughters_pdg);
    _muon_tree->Branch("muon_daughters_e", "std::vector<float>", &_muon_daughters_e);
    _muon_tree->Branch("muon_daughters_startprocess", "std::vector<std::string>", &_muon_daughters_startprocess);
    _muon_tree->Branch("muon_daughters_endprocess", "std::vector<std::string>", &_muon_daughters_endprocess);
    _muon_tree->Branch("muon_daughters_start_x", "std::vector<float>", &_muon_daughters_start_x);
    _muon_tree->Branch("muon_daughters_start_y", "std::vector<float>", &_muon_daughters_start_y);
    _muon_tree->Branch("muon_daughters_start_z", "std::vector<float>", &_muon_daughters_start_z);
    _muon_tree->Branch("muon_daughters_end_x", "std::vector<float>", &_muon_daughters_end_x);
    _muon_tree->Branch("muon_daughters_end_y", "std::vector<float>", &_muon_daughters_end_y);
    _muon_tree->Branch("muon_daughters_end_z", "std::vector<float>", &_muon_daughters_end_z);

    _muon_tree->Branch("muon_genealogy_pdg", "std::vector<int>", &_muon_genealogy_pdg);
    _muon_tree->Branch("muon_genealogy_startprocess", "std::vector<std::string>", &_muon_genealogy_startprocess);
    _muon_tree->Branch("muon_genealogy_endprocess", "std::vector<std::string>", &_muon_genealogy_endprocess);
    _muon_tree->Branch("muon_genealogy_mother", "std::vector<int>", &_muon_genealogy_mother);
    _muon_tree->Branch("muon_genealogy_e", "std::vector<float>", &_muon_genealogy_e);
    _muon_tree->Branch("muon_genealogy_start_x", "std::vector<float>", &_muon_genealogy_start_x);
    _muon_tree->Branch("muon_genealogy_start_y", "std::vector<float>", &_muon_genealogy_start_y);
    _muon_tree->Branch("muon_genealogy_start_z", "std::vector<float>", &_muon_genealogy_start_z);
    _muon_tree->Branch("muon_genealogy_end_x", "std::vector<float>", &_muon_genealogy_end_x);
    _muon_tree->Branch("muon_genealogy_end_y", "std::vector<float>", &_muon_genealogy_end_y);
    _muon_tree->Branch("muon_genealogy_end_z", "std::vector<float>", &_muon_genealogy_end_z);
    _muon_tree->Branch("muon_genealogy_trackid", "std::vector<float>", &_muon_genealogy_trackid);

    // Add the neutrino information also to this tree so we can better study pi0s
    // when looking at beam neutrinos
    _muon_tree->Branch("nu_e", &_nu_e, "nu_e/F");
    _muon_tree->Branch("nu_pdg", &_nu_pdg, "nu_pdg/I");
    _muon_tree->Branch("nu_ccnc", &_nu_ccnc, "nu_ccnc/I");
    _muon_tree->Branch("nu_mode", &_nu_mode, "nu_mode/I");
    _muon_tree->Branch("nu_vtx_x", &_nu_vtx_x, "nu_vtx_x/F");
    _muon_tree->Branch("nu_vtx_y", &_nu_vtx_y, "nu_vtx_y/F");
    _muon_tree->Branch("nu_vtx_z", &_nu_vtx_z, "nu_vtx_z/F");
    _muon_tree->Branch("nu_pip_mult", &_nu_pip_mult, "nu_pip_mult/I");
    _muon_tree->Branch("nu_pi0_mult", &_nu_pi0_mult, "nu_pi0_mult/I");
    _muon_tree->Branch("nu_p_mult", &_nu_p_mult, "nu_p_mult/I");
    _muon_tree->Branch("pars_pdg", "std::vector<int>", &_pars_pdg);
    _muon_tree->Branch("pars_e", "std::vector<float>", &_pars_e);

  }

  _sr_tree = fs->make<TTree>("pottree","");
  _sr_tree->Branch("run", &_sr_run, "run/I");
  _sr_tree->Branch("subrun", &_sr_subrun, "subrun/I");
  _sr_tree->Branch("begintime", &_sr_begintime, "begintime/D");
  _sr_tree->Branch("endtime", &_sr_endtime, "endtime/D");
  _sr_tree->Branch("pot", &_sr_pot, "pot/D");

  // Iterate over all TPC's to get bounding box that covers volumes of each individual TPC in the detector
  art::ServiceHandle<geo::Geometry const> geo;
  for (auto const& tpc : geo->Iterate<geo::TPCGeo>()) {
    _x_min = std::min(_x_min, tpc.BoundingBox().MinX());
    _y_min = std::min(_y_min, tpc.BoundingBox().MinY());
    _z_min = std::min(_z_min, tpc.BoundingBox().MinZ());
    _x_max = std::max(_x_max, tpc.BoundingBox().MaxX());
    _y_max = std::max(_y_max, tpc.BoundingBox().MaxY());
    _z_max = std::max(_z_max, tpc.BoundingBox().MaxZ());
  }

  std::cout << "TPC limits: " << std::endl;
  std::cout << "\tx_max\t" << _x_max << std::endl;
  std::cout << "\tx_min\t" << _x_min << std::endl;
  std::cout << "\ty_max\t" << _y_max << std::endl;
  std::cout << "\ty_min\t" << _y_min << std::endl;
  std::cout << "\tz_max\t" << _z_max << std::endl;
  std::cout << "\tz_min\t" << _z_min << std::endl;
}

void obana::OBAnaICARUS::analyze(art::Event const& e)
{
  // Implementation of required member function here.
  _run = e.id().run();
  _subrun = e.id().subRun();
  _event =  e.id().event();

  //  std::cout << "Run: " << _run << ",  subrun: " << _subrun <<",  event: " << _event << std::endl;
  _uuid = boost::uuids::random_generator()();
  _uuid_str = boost::lexical_cast<std::string>(_uuid) + "-" + _event;
  //std::cout << "Event uuid " << _uuid_str << std::endl;


  std::set<std::string> volnameset(_overburden_volumes.begin(), _overburden_volumes.end());
  _part_filter = CreateParticleVolumeFilter(volnameset);

  const geo::GeometryCore*           fGeometry;             ///< pointer to Geometry service
  fGeometry           = lar::providerFrom<geo::Geometry>();


  //
  // MCTruth
  //
  auto mct_h = e.getHandle<std::vector<simb::MCTruth>>(_mctruth_producer);
  if(mct_h){
    //
    // Loop over the neutrino interactions in this event
    //
    for (simb::MCTruth const& mct : *mct_h) {
      // if (mct_v.at(i)->Origin() != simb::kBeamNeutrino) {
      //   std::cout << "[OverburdenAna] MCTruth from generator does not have neutrino origin?!" << std::endl;
      // }

      if(!mct.NeutrinoSet()) {
        break;
      }

      _nu_e = mct.GetNeutrino().Nu().E();
      _nu_pdg = mct.GetNeutrino().Nu().PdgCode();
      _nu_ccnc = mct.GetNeutrino().CCNC();
      _nu_mode = mct.GetNeutrino().Mode();
      _nu_int_type = mct.GetNeutrino().InteractionType();
      _nu_vtx_x = mct.GetNeutrino().Nu().Vx();
      _nu_vtx_y = mct.GetNeutrino().Nu().Vy();
      _nu_vtx_z = mct.GetNeutrino().Nu().Vz();
      _nu_px = mct.GetNeutrino().Nu().Px();
      _nu_py = mct.GetNeutrino().Nu().Py();
      _nu_pz = mct.GetNeutrino().Nu().Pz();


      //      std::cout<< _nu_vtx_x << "\t"<<_nu_vtx_y <<"\t"<< _nu_vtx_z << std::endl;

      // Do not save neutrinos interacting in the detector if we are using a dirt sample
      if (_simulating_dirt && InDetector(_nu_vtx_x, _nu_vtx_y, _nu_vtx_z)) {
	return;
	//break; 
	//continue;
      }

      // std::cout << " 2nd time ----- Run: " << _run << ",  subrun: " << _subrun <<",  event: " << _event <<  std::endl;

      _pars_pdg.clear();
      _pars_e.clear();
      _nu_pip_mult = 0;
      _nu_pi0_mult = 0;
      _nu_p_mult = 0; 

      for (int p = 0; p < mct.NParticles(); p++) {
        auto const & mcp = mct.GetParticle(p);

	if (mcp.StatusCode() != 1) continue;

	_pars_pdg.push_back(mcp.PdgCode());
	_pars_e.push_back(mcp.E());

	if (mcp.PdgCode() == 111) {
	  _nu_pi0_mult++;
	  //	  std::cout << "There is a GENIE pi0 with energy " << mcp.E() << std::endl;
	} else if (std::abs(mcp.PdgCode()) == 211) {
	  _nu_pip_mult++;
	}
	else if (std::abs(mcp.PdgCode()) == 2112) {
	  _nu_p_mult++;
	}
      }
      // std::cout << "_nu_e " << _nu_e << std::endl;
    }
  } else {
    std::cout << "MCTruth product " << _mctruth_producer << " not found..." << std::endl;
  }

  //std::cout << " 3rd time ----------------------- Run: " << _run << ",  subrun: " << _subrun <<",  event: " << _event <<  std::endl;

  auto mcp_h = e.getHandle<std::vector<simb::MCParticle>>(_mcparticle_producer);
  if(!mcp_h){
    std::cout << "MCParticle product " << _mcparticle_producer << " not found..." << std::endl;
    _tree->Fill();
    return;
    //throw std::exception();
  }

  clear_vectors();

  for (simb::MCParticle const& mcp : *mcp_h) {
    _trackid_to_mcparticle[mcp.TrackId()] = mcp;
  }

  for (simb::MCParticle const& mcp : *mcp_h) {
    //      bool in_det = InDetector(mcp);

    // Only save the MCP if it's a primary, or if it crosses the det
    //    if (mcp.Process() == "primary" || in_det) {
    if (mcp.Process() == "primary") {

      _mcp_px.push_back(mcp.Px());
      _mcp_py.push_back(mcp.Py());
      _mcp_pz.push_back(mcp.Pz());
      _mcp_e.push_back(mcp.E());

      _mcp_vx.push_back(mcp.Vx());
      _mcp_vy.push_back(mcp.Vy());
      _mcp_vz.push_back(mcp.Vz());
      _mcp_endx.push_back(mcp.EndX());
      _mcp_endy.push_back(mcp.EndY());
      _mcp_endz.push_back(mcp.EndZ());

      _mcp_pdg.push_back(mcp.PdgCode());
      _mcp_mother.push_back(mcp.Mother());
      _mcp_status_code.push_back(mcp.StatusCode());
      _mcp_process.push_back(mcp.Process());
      _mcp_end_process.push_back(mcp.EndProcess());
      _mcp_trackid.push_back(mcp.TrackId());

      //      _mcp_intpc.push_back(InDetector(mcp));
      // _neut_par_uuid = _uuid_str + "-" + std::to_string(mcp.TrackId());

      int step = 0;
      bool in_det = InDetector(mcp, step);
      _mcp_intpc.push_back(in_det);      
      
      if (mcp.PdgCode() == 2112 && in_det){ /// save the daughter of neutron entring to TPCs
        unsigned int nSec = mcp.NumberDaughters();
	for (size_t d = 0; d < nSec; ++d) {
          auto d_search = _trackid_to_mcparticle.find(mcp.Daughter(d));
	  if (d_search != _trackid_to_mcparticle.end()) {
	    auto const& daughter = d_search->second;
	    _neut_daughters_pdg.push_back(daughter.PdgCode());
	    _neut_daughters_px.push_back(daughter.Px());
	    _neut_daughters_py.push_back(daughter.Py());
	    _neut_daughters_pz.push_back(daughter.Pz());	  
	    _neut_daughters_e.push_back(daughter.E());
	    _neut_daughters_startprocess.push_back(daughter.Process());
	    _neut_daughters_endprocess.push_back(daughter.EndProcess());
	    _neut_daughters_start_x.push_back(daughter.Vx());
	    _neut_daughters_start_y.push_back(daughter.Vy());
	    _neut_daughters_start_z.push_back(daughter.Vz());
	    _neut_daughters_end_x.push_back(daughter.EndX());
	    _neut_daughters_end_y.push_back(daughter.EndY());
	    _neut_daughters_end_z.push_back(daughter.EndZ());
	    _neut_daughters_trackid.push_back(daughter.TrackId());
	    
	    // std::cout << "no. of daughter: \t" << nSec <<", pdg of daughter: \t" << daughter.PdgCode()
	    //	      << ", daughter process: \t" << daughter.Process() << ", trackid: \t" << daughter.TrackId()<< std::endl; 
	    /// save the grand daughter of neutron if daughters are from pions
	    if (std::abs(daughter.PdgCode())==211 || daughter.PdgCode()==111 || daughter.PdgCode()==2112){
	      unsigned int ngd = daughter.NumberDaughters();
	      for (size_t gd = 0; gd < ngd ; ++gd) {
		auto gd_search = _trackid_to_mcparticle.find(daughter.Daughter(gd));
		if (gd_search != _trackid_to_mcparticle.end()) {
		  auto const& granddaughter = gd_search->second;
		  //std::cout << "no. of grand daughter: \t" << ngd <<", pdg of grand daughter: \t" << granddaughter.PdgCode()
		  //	    << ", grand daughter process: \t" << granddaughter.Process() << ", trackid: \t" << granddaughter.TrackId()<< std::endl; 

		  _neut_grand_daughters_pdg.push_back(granddaughter.PdgCode());
		  _neut_grand_daughters_px.push_back(granddaughter.Px());
		  _neut_grand_daughters_py.push_back(granddaughter.Py());
		  _neut_grand_daughters_pz.push_back(granddaughter.Pz());
		  _neut_grand_daughters_e.push_back(granddaughter.E());
		  _neut_grand_daughters_startprocess.push_back(granddaughter.Process());
		  _neut_grand_daughters_endprocess.push_back(granddaughter.EndProcess());
		  _neut_grand_daughters_start_x.push_back(granddaughter.Vx());
		  _neut_grand_daughters_start_y.push_back(granddaughter.Vy());
		  _neut_grand_daughters_start_z.push_back(granddaughter.Vz());
		  _neut_grand_daughters_end_x.push_back(granddaughter.EndX());
		  _neut_grand_daughters_end_y.push_back(granddaughter.EndY());
		  _neut_grand_daughters_end_z.push_back(granddaughter.EndZ());
		  _neut_grand_daughters_trackid.push_back(granddaughter.TrackId());
                } // end of grand daughter trackid
              } // end of loop over grand daughter
            } // end of searching grand daughter
          } // end of  daughter trackid
        }// end of loop over daughter
      } //end of searching daughter
      
      if (in_det) {
	/*
          if (mcp.PdgCode()==13 && mcp.E() >= 0.01 ){
          std::cout << "photon: " << mcp.PdgCode()  << " , has energy : " << mcp.E()
          << " , process: " << mcp.Process() << " , and mother is : "
          << mcp.Mother() << std::endl;
	  
          std::cout << "step: " << step << " , has energy : " << mcp.E(step) << std::endl;
          if ((mcp.E(step)-0.105658) > 0.0) std::cout << "step: " << step << " , has energy : " << mcp.E(step)-0.105658 << std::endl;
	  
	  }*/
        _mcp_intpc_nu_e.push_back(mcp.E(step));
        _mcp_intpc_e.push_back(mcp.E(step-1));
        //_mcp_intpc_e_previousstep.push_back(mcp.E(step-1));
	
        //	std::cout << "step: " << step << " , has energy : " << mcp.E(step) << "\t"<< mcp.E(step)-0.105658 << std::endl;
        //std::cout << "previous step: " << step-1 << " , has energy : " << mcp.E(step-1) << "\t"<< mcp.E(step-1)-0.105658 << std::endl;
      } else {
        _mcp_intpc_e.push_back(-9999.);
        _mcp_intpc_nu_e.push_back(-9999.);
	//_mcp_intpc_e_previousstep.push_back(-9999.);
      }// end of in detector
    } // end of primary selection 
  } // end of mcp loop

/*
  for (unsigned int i =0; i < _mcp_intpc_e.size(); i++){
  if (_mcp_pdg[i]==13)
  std::cout << "event: " <<_event << " , has energy : " << _mcp_intpc_e[i] << std::endl;
  }
*/

  //
  // MCTrack
  //
  auto mc_track_h = e.getHandle<std::vector<sim::MCTrack>>(_mctrack_producer);
  if(!mc_track_h){
    std::cout << "MCTrack product " << _mctrack_producer << " not found..." << std::endl;
    throw std::exception();
  }

  for (sim::MCTrack const& mc_track : *mc_track_h) {
    // std::cout << "MCTrack " << i << ": ancestor pdg " << mc_track.AncestorPdgCode()
    //                               << ", ancestor process " << mc_track.AncestorProcess()
    //                               << ", mother pdg " << mc_track.MotherPdgCode()
    //                               << ", mother process " << mc_track.MotherProcess()
    //                              << "| PDG " << mc_track.PdgCode()
    //                               << ", process " << mc_track.Process()
    //                               << std::endl;

    auto iter = _trackid_to_mcparticle.find(mc_track.TrackID());
    int mother_in_ob = -1;
    if (iter != _trackid_to_mcparticle.end()) {
      mother_in_ob = FindMotherInOverburden(iter->second);
    }

    // Don't save MCT with energy less than 1 MeV
    if (mc_track.Start().E() < 1) { // MeV
      _n_mct_lt1 ++;
      if (mother_in_ob != -1){
        _n_mct_lt1_from_ob ++;
      }
      continue;
    }

    // Don't save MCS that are not in the TPCs
    if (mc_track.size() == 0) {
      continue;
    }

    geo::Point_t mctrackstartPoint(mc_track.Start().X(),mc_track.Start().Y(), mc_track.Start().Z());
    geo::Point_t mctrackendPoint(mc_track.End().X(),mc_track.End().Y(), mc_track.End().Z());

    const geo::TPCGeo* mctpcstartGeo = fGeometry->PositionToTPCptr(mctrackstartPoint);
    const geo::TPCGeo* mctpcendGeo = fGeometry->PositionToTPCptr(mctrackendPoint);

    if (!mctpcstartGeo || !mctpcendGeo) continue;

    //    if (tpcGeo->ID() != C0id) continue; // point not in cryostat 0
    if (!mctpcstartGeo->ActiveBoundingBox().ContainsPosition(mctrackstartPoint) ||
        !mctpcendGeo->ActiveBoundingBox().ContainsPosition(mctrackendPoint)) continue; // out of active volume


    _mct_pdg.push_back(mc_track.PdgCode());
    _mct_process.push_back(mc_track.Process());

    _mct_start_x.push_back(mc_track.Start().X());
    _mct_start_y.push_back(mc_track.Start().Y());
    _mct_start_z.push_back(mc_track.Start().Z());

    _mct_end_x.push_back(mc_track.End().X());
    _mct_end_y.push_back(mc_track.End().Y());
    _mct_end_z.push_back(mc_track.End().Z());

    _mct_start_px.push_back(mc_track.Start().Px());
    _mct_start_py.push_back(mc_track.Start().Py());
    _mct_start_pz.push_back(mc_track.Start().Pz());
    _mct_start_e.push_back(mc_track.Start().E());

    _mct_mother_pdg.push_back(mc_track.MotherPdgCode());
    _mct_mother_process.push_back(mc_track.MotherProcess());
    _mct_ancestor_pdg.push_back(mc_track.AncestorPdgCode());
    _mct_ancestor_process.push_back(mc_track.AncestorProcess());

    _mct_mother_in_ob_trackid.push_back(mother_in_ob);

  }

  //
  // MCShower
  //
  auto mc_shower_h = e.getHandle<std::vector<sim::MCShower>>(_mcshower_producer);
  if(!mc_shower_h){
    std::cout << "MCShower product " << _mcshower_producer << " not found..." << std::endl;
    throw std::exception();
  }

  for (sim::MCShower const& mc_shower : *mc_shower_h) {
    // std::cout << "MCShower " << i << ": ancestor pdg " << mc_shower.AncestorPdgCode()
    //                               << ", ancestor process " << mc_shower.AncestorProcess()
    //                               << ", mother pdg " << mc_shower.MotherPdgCode()
    //                               << ", mother process " << mc_shower.MotherProcess()
    //                               << "| PDG " << mc_shower.PdgCode()
    //                               << ", process " << mc_shower.Process()
    //                               << std::endl;

    auto iter = _trackid_to_mcparticle.find(mc_shower.TrackID());
    int mother_in_ob = -1;
    if (iter != _trackid_to_mcparticle.end()) {
      mother_in_ob = FindMotherInOverburden(iter->second);
    }

    // Don't save MCS with energy less than 1 MeV
    if (mc_shower.Start().E() < 1) { // MeV
      _n_mcs_lt1 ++;
      if (mother_in_ob != -1){
        _n_mcs_lt1_from_ob ++;
      }
      continue;
    }

    /// Don't save MCS that are not in the TPCs
    // Special case for photon showers, which can start outside
   
    bool end_in_det = InDetector(mc_shower.End().X(), mc_shower.End().Y(), mc_shower.End().Z());

    if (!end_in_det) {
      continue;
    }

    _mcs_pdg.push_back(mc_shower.PdgCode());
    _mcs_process.push_back(mc_shower.Process());

    _mcs_start_x.push_back(mc_shower.Start().X());
    _mcs_start_y.push_back(mc_shower.Start().Y());
    _mcs_start_z.push_back(mc_shower.Start().Z());

    _mcs_end_x.push_back(mc_shower.End().X());
    _mcs_end_y.push_back(mc_shower.End().Y());
    _mcs_end_z.push_back(mc_shower.End().Z());

    _mcs_start_px.push_back(mc_shower.Start().Px());
    _mcs_start_py.push_back(mc_shower.Start().Py());
    _mcs_start_pz.push_back(mc_shower.Start().Pz());
    _mcs_start_e.push_back(mc_shower.Start().E());
    _mcs_charge_col.push_back(mc_shower.Charge(2));
    _mcs_charge_ind2.push_back(mc_shower.Charge(1));
    _mcs_charge_ind1.push_back(mc_shower.Charge(0));

    _mcs_mother_pdg.push_back(mc_shower.MotherPdgCode());
    _mcs_mother_trackid.push_back(mc_shower.MotherTrackID());
    _mcs_mother_start_x.push_back(mc_shower.MotherStart().X());
    _mcs_mother_start_y.push_back(mc_shower.MotherStart().Y());
    _mcs_mother_start_z.push_back(mc_shower.MotherStart().Z());
    _mcs_mother_start_e.push_back(mc_shower.MotherStart().E());
    _mcs_mother_end_x.push_back(mc_shower.MotherEnd().X());
    _mcs_mother_end_y.push_back(mc_shower.MotherEnd().Y());
    _mcs_mother_end_z.push_back(mc_shower.MotherEnd().Z());
    _mcs_mother_end_e.push_back(mc_shower.MotherEnd().E());

    _mcs_mother_process.push_back(mc_shower.MotherProcess());
    _mcs_ancestor_pdg.push_back(mc_shower.AncestorPdgCode());
    _mcs_ancestor_process.push_back(mc_shower.AncestorProcess());
    _mcs_ancestor_start_e.push_back(mc_shower.AncestorStart().E());
    _mcs_ancestor_end_e.push_back(mc_shower.AncestorEnd().E());
    _mcs_mother_in_ob_trackid.push_back(mother_in_ob);

    if (mc_shower.MotherPdgCode() == 111 && _save_pi0_tree) {
      SavePi0ShowerInfo(mc_shower.MotherTrackID());
    }

    if (std::abs(mc_shower.MotherPdgCode()) == 13 && _save_muon_tree) {
      SaveMuonShowerInfo(mc_shower.MotherTrackID());
    }
    
  }


  _tree->Fill();
}

void obana::OBAnaICARUS::SavePi0ShowerInfo(int pi0_track_id) {
  //std::cout << "SavePi0ShowerInfo***, pi0_track_id: " << pi0_track_id << std::endl;
  auto it = std::find(_pi0_ids.begin(), _pi0_ids.end(), pi0_track_id);
  if (it != _pi0_ids.end()) {
    return;
  }
  //std::cout << "SavePi0ShowerInfo***, got it" << std::endl;
  _pi0_ids.push_back(pi0_track_id);

  // Get the pi0 MCParticle
  auto iter = _trackid_to_mcparticle.find(pi0_track_id);
  if (iter == _trackid_to_mcparticle.end()) {
   return;
   }

 
  simb::MCParticle pi0_mcp = iter->second;
  // std::cout << "Pi0 MCP, PDG = " << pi0_mcp.PdgCode() << ", E = " << pi0_mcp.E() << ", n daughters " << pi0_mcp.NumberDaughters() << std::endl;

  // Save the information on the pi0 itself
  _pi0_par_e = pi0_mcp.E();
  _pi0_par_start_x = pi0_mcp.Vx();
  _pi0_par_start_y = pi0_mcp.Vy();
  _pi0_par_start_z = pi0_mcp.Vz();
  _pi0_par_end_x = pi0_mcp.EndX();
  _pi0_par_end_y = pi0_mcp.EndY();
  _pi0_par_end_z = pi0_mcp.EndZ();

  // Get the pi0 mother MCParticle
  //  iter = _trackid_to_mcparticle.find(pi0_mcp.Mother());
  // if (iter == _trackid_to_mcparticle.end()) {
  //  return;
  // }

  simb::MCParticle pi0_mother_mcp;
  if (pi0_mcp.Mother() == 0) {
    // Use itself if this pi0 is a primary
    pi0_mother_mcp = pi0_mcp;
  } else {
    iter = _trackid_to_mcparticle.find(pi0_mcp.Mother());
    if (iter == _trackid_to_mcparticle.end()) {
      return;
    }
    pi0_mother_mcp = iter->second;
  }

  // simb::MCParticle pi0_mother_mcp = iter->second;
  _pi0_par_mother_pdg = pi0_mother_mcp.PdgCode();
  _pi0_par_mother_e = pi0_mother_mcp.E();
  _pi0_par_mother_end_e = pi0_mother_mcp.EndE();

  // std::cout << "Pi0 MCP Mother, PDG = " << pi0_mother_mcp.PdgCode() << ", E = " << pi0_mother_mcp.E() << std::endl;

  // Get the daughters of the pi0 mother
  for (int d = 0; d < pi0_mother_mcp.NumberDaughters(); d++) {
    iter = _trackid_to_mcparticle.find(pi0_mother_mcp.Daughter(d));
    if (iter != _trackid_to_mcparticle.end()) {
      simb::MCParticle daughter = iter->second;
      _pi0_event_particles_pdg.push_back(daughter.PdgCode());
      _pi0_event_particles_e.push_back(daughter.E());
    }
  }

  // Save the mother first...
  _pi0_genealogy_pdg.push_back(pi0_mother_mcp.PdgCode());
  _pi0_genealogy_startprocess.push_back(pi0_mother_mcp.Process());
  _pi0_genealogy_endprocess.push_back(pi0_mother_mcp.EndProcess());
  _pi0_genealogy_mother.push_back(pi0_mother_mcp.Mother());
  _pi0_genealogy_e.push_back(pi0_mother_mcp.E());
  _pi0_genealogy_start_x.push_back(pi0_mother_mcp.Vx());
  _pi0_genealogy_start_y.push_back(pi0_mother_mcp.Vy());
  _pi0_genealogy_start_z.push_back(pi0_mother_mcp.Vy());
  _pi0_genealogy_end_x.push_back(pi0_mother_mcp.EndX());
  _pi0_genealogy_end_y.push_back(pi0_mother_mcp.EndY());
  _pi0_genealogy_end_z.push_back(pi0_mother_mcp.EndZ());
  _pi0_genealogy_trackid.push_back(pi0_mother_mcp.TrackId());


  // ... then save all the other ancestors
  simb::MCParticle mcp = pi0_mother_mcp;
  while(true) {
    iter = _trackid_to_mcparticle.find(mcp.Mother());
    if (iter == _trackid_to_mcparticle.end()) {
      break;
    }
    mcp = iter->second;
    _pi0_genealogy_pdg.push_back(mcp.PdgCode());
    _pi0_genealogy_startprocess.push_back(mcp.Process());
    _pi0_genealogy_endprocess.push_back(mcp.EndProcess());
    _pi0_genealogy_mother.push_back(mcp.Mother());
    _pi0_genealogy_e.push_back(mcp.E());
    _pi0_genealogy_start_x.push_back(mcp.Vx());
    _pi0_genealogy_start_y.push_back(mcp.Vy());
    _pi0_genealogy_start_z.push_back(mcp.Vy());
    _pi0_genealogy_end_x.push_back(mcp.EndX());
    _pi0_genealogy_end_y.push_back(mcp.EndY());
    _pi0_genealogy_end_z.push_back(mcp.EndZ());
    _pi0_genealogy_trackid.push_back(mcp.TrackId());
  }

  _pi0_par_ancestor_trackid = _pi0_genealogy_trackid.back();
  _pi0_par_ancestor_uuid = _uuid_str + "-" + std::to_string(_pi0_par_ancestor_trackid);
  //std::cout << "Pi0 uuid " << _pi0_par_ancestor_uuid << std::endl;

  // Get the dautghers of this pi0
  for (int d = 0; d < pi0_mcp.NumberDaughters(); d++) {
    iter = _trackid_to_mcparticle.find(pi0_mcp.Daughter(d));
    if (iter != _trackid_to_mcparticle.end()) {
      simb::MCParticle daughter = iter->second;
      _pi0_daughters_pdg.push_back(daughter.PdgCode());
      _pi0_daughters_e.push_back(daughter.E());
      _pi0_daughters_startprocess.push_back(daughter.Process());
      _pi0_daughters_endprocess.push_back(daughter.EndProcess());
      _pi0_daughters_start_x.push_back(daughter.Vx());
      _pi0_daughters_start_y.push_back(daughter.Vy());
      _pi0_daughters_start_z.push_back(daughter.Vz());
      _pi0_daughters_end_x.push_back(daughter.EndX());
      _pi0_daughters_end_y.push_back(daughter.EndY());
      _pi0_daughters_end_z.push_back(daughter.EndZ());
    }
  }

  // Fill the tree and reset the variables
  _pi0_tree->Fill();


  _pi0_event_particles_pdg.clear();
  _pi0_event_particles_e.clear();

  _pi0_daughters_pdg.clear();
  _pi0_daughters_e.clear();
  _pi0_daughters_startprocess.clear();
  _pi0_daughters_endprocess.clear();
  _pi0_daughters_start_x.clear();
  _pi0_daughters_start_y.clear();
  _pi0_daughters_start_z.clear();
  _pi0_daughters_end_x.clear();
  _pi0_daughters_end_y.clear();
  _pi0_daughters_end_z.clear();

  _pi0_genealogy_pdg.clear();
  _pi0_genealogy_startprocess.clear();
  _pi0_genealogy_endprocess.clear();
  _pi0_genealogy_mother.clear();
  _pi0_genealogy_e.clear();
  _pi0_genealogy_start_x.clear();
  _pi0_genealogy_start_y.clear();
  _pi0_genealogy_start_z.clear();
  _pi0_genealogy_end_x.clear();
  _pi0_genealogy_end_y.clear();
  _pi0_genealogy_end_z.clear();
  _pi0_genealogy_trackid.clear();
}


void obana::OBAnaICARUS::SaveMuonShowerInfo(int muon_track_id) {
  //std::cout << "SaveMuonShowerInfo***, muon_track_id: " << muon_track_id << std::endl;
  auto it = std::find(_muon_ids.begin(), _muon_ids.end(), muon_track_id);
  if (it != _muon_ids.end()) {
    return;
  }
  //std::cout << "SaveMuonShowerInfo***, got it" << std::endl;
  _muon_ids.push_back(muon_track_id);

  // Get the muon MCParticle
  auto iter = _trackid_to_mcparticle.find(muon_track_id);
  if (iter == _trackid_to_mcparticle.end()) {
    return;
  }
  simb::MCParticle muon_mcp = iter->second;
  // std::cout << "Muon MCP, PDG = " << muon_mcp.PdgCode() << ", E = " << muon_mcp.E() << ", n daughters " << muon_mcp.NumberDaughters() << std::endl;

  // Save the information on the muon itself
  _muon_par_e = muon_mcp.E();
  _muon_par_start_x = muon_mcp.Vx();
  _muon_par_start_y = muon_mcp.Vy();
  _muon_par_start_z = muon_mcp.Vz();
  _muon_par_end_x = muon_mcp.EndX();
  _muon_par_end_y = muon_mcp.EndY();
  _muon_par_end_z = muon_mcp.EndZ();


  simb::MCParticle muon_mother_mcp;
  if (muon_mcp.Mother() == 0) {
    // Use itself if this pi0 is a primary
    muon_mother_mcp = muon_mcp;
  } else {
    iter = _trackid_to_mcparticle.find(muon_mcp.Mother());
    if (iter == _trackid_to_mcparticle.end()) {
      return;
    }
    muon_mother_mcp = iter->second;
  }

  // Get the muon mother MCParticle
  //iter = _trackid_to_mcparticle.find(muon_mcp.Mother());
  //if (iter == _trackid_to_mcparticle.end()) {
  // return;
  //}

  //simb::MCParticle muon_mother_mcp = iter->second;
  _muon_par_mother_pdg = muon_mother_mcp.PdgCode();
  _muon_par_mother_e = muon_mother_mcp.E();
  _muon_par_mother_end_e = muon_mother_mcp.EndE();

  // std::cout << "Muon MCP Mother, PDG = " << muon_mother_mcp.PdgCode() << ", E = " << muon_mother_mcp.E() << std::endl;

  // Get the daughters of the muon mother
  for (int d = 0; d < muon_mother_mcp.NumberDaughters(); d++) {
    iter = _trackid_to_mcparticle.find(muon_mother_mcp.Daughter(d));
    if (iter != _trackid_to_mcparticle.end()) {
      simb::MCParticle daughter = iter->second;
      _muon_event_particles_pdg.push_back(daughter.PdgCode());
      _muon_event_particles_e.push_back(daughter.E());
    }
  }

  // Save the mother first...
  _muon_genealogy_pdg.push_back(muon_mother_mcp.PdgCode());
  _muon_genealogy_startprocess.push_back(muon_mother_mcp.Process());
  _muon_genealogy_endprocess.push_back(muon_mother_mcp.EndProcess());
  _muon_genealogy_mother.push_back(muon_mother_mcp.Mother());
  _muon_genealogy_e.push_back(muon_mother_mcp.E());
  _muon_genealogy_start_x.push_back(muon_mother_mcp.Vx());
  _muon_genealogy_start_y.push_back(muon_mother_mcp.Vy());
  _muon_genealogy_start_z.push_back(muon_mother_mcp.Vy());
  _muon_genealogy_end_x.push_back(muon_mother_mcp.EndX());
  _muon_genealogy_end_y.push_back(muon_mother_mcp.EndY());
  _muon_genealogy_end_z.push_back(muon_mother_mcp.EndZ());
  _muon_genealogy_trackid.push_back(muon_mother_mcp.TrackId());


  // ... then save all the other ancestors
  simb::MCParticle mcp = muon_mother_mcp;
  while(true) {
    iter = _trackid_to_mcparticle.find(mcp.Mother());
    if (iter == _trackid_to_mcparticle.end()) {
      break;
    }
    mcp = iter->second;
    _muon_genealogy_pdg.push_back(mcp.PdgCode());
    _muon_genealogy_startprocess.push_back(mcp.Process());
    _muon_genealogy_endprocess.push_back(mcp.EndProcess());
    _muon_genealogy_mother.push_back(mcp.Mother());
    _muon_genealogy_e.push_back(mcp.E());
    _muon_genealogy_start_x.push_back(mcp.Vx());
    _muon_genealogy_start_y.push_back(mcp.Vy());
    _muon_genealogy_start_z.push_back(mcp.Vy());
    _muon_genealogy_end_x.push_back(mcp.EndX());
    _muon_genealogy_end_y.push_back(mcp.EndY());
    _muon_genealogy_end_z.push_back(mcp.EndZ());
    _muon_genealogy_trackid.push_back(mcp.TrackId());
  }

  _muon_par_ancestor_trackid = _muon_genealogy_trackid.back();
  _muon_par_ancestor_uuid = _uuid_str + "-" + std::to_string(_muon_par_ancestor_trackid);
  //std::cout << "Muon uuid " << _muon_par_ancestor_uuid << std::endl;

  // Get the dautghers of this muon
  for (int d = 0; d < muon_mcp.NumberDaughters(); d++) {
    iter = _trackid_to_mcparticle.find(muon_mcp.Daughter(d));
    if (iter != _trackid_to_mcparticle.end()) {
      simb::MCParticle daughter = iter->second;
      _muon_daughters_pdg.push_back(daughter.PdgCode());
      _muon_daughters_e.push_back(daughter.E());
      _muon_daughters_startprocess.push_back(daughter.Process());
      _muon_daughters_endprocess.push_back(daughter.EndProcess());
      _muon_daughters_start_x.push_back(daughter.Vx());
      _muon_daughters_start_y.push_back(daughter.Vy());
      _muon_daughters_start_z.push_back(daughter.Vz());
      _muon_daughters_end_x.push_back(daughter.EndX());
      _muon_daughters_end_y.push_back(daughter.EndY());
      _muon_daughters_end_z.push_back(daughter.EndZ());
    }
  }

  // Fill the tree and reset the variables
  _muon_tree->Fill();


  _muon_event_particles_pdg.clear();
  _muon_event_particles_e.clear();

  _muon_daughters_pdg.clear();
  _muon_daughters_e.clear();
  _muon_daughters_startprocess.clear();
  _muon_daughters_endprocess.clear();
  _muon_daughters_start_x.clear();
  _muon_daughters_start_y.clear();
  _muon_daughters_start_z.clear();
  _muon_daughters_end_x.clear();
  _muon_daughters_end_y.clear();
  _muon_daughters_end_z.clear();

  _muon_genealogy_pdg.clear();
  _muon_genealogy_startprocess.clear();
  _muon_genealogy_endprocess.clear();
  _muon_genealogy_mother.clear();
  _muon_genealogy_e.clear();
  _muon_genealogy_start_x.clear();
  _muon_genealogy_start_y.clear();
  _muon_genealogy_start_z.clear();
  _muon_genealogy_end_x.clear();
  _muon_genealogy_end_y.clear();
  _muon_genealogy_end_z.clear();
  _muon_genealogy_trackid.clear();

}


int obana::OBAnaICARUS::FindMotherInOverburden(simb::MCParticle mcp) {

  if (mcp.Process() == "primary") {
    return -1;
  }

  // We use this filter not to actualy filter, but to check if
  // the vertex of this particle is in the OB
  bool vtx_in_ob = _part_filter->mustKeep(Point_t{{ mcp.Vx(),
          mcp.Vy(),
          mcp.Vz() }});

  auto iter = _trackid_to_mcparticle.find(mcp.Mother());
  if (iter == _trackid_to_mcparticle.end()) {
    return -1;
  }
  auto mother = iter->second;


  if (vtx_in_ob                          // If this particle has a vertex in the OB
      && mcp.Process() != "primary"        // and this particle is not a primary one
      // && mother.Process() == "primary"     // and the mother of it is a primary
      ) {
    return mcp.TrackId();                // Then return it, as is something created by a primary in the OB
  }

  return FindMotherInOverburden(mother);

}

std::unique_ptr<util::PositionInVolumeFilter> obana::OBAnaICARUS::CreateParticleVolumeFilter
(std::set<std::string> const& vol_names) const
{

  // if we don't have favourite volumes, don't even bother creating a filter
  if (vol_names.empty()) return {};

  auto const& geom = *art::ServiceHandle<geo::Geometry const>();

  std::vector<std::vector<TGeoNode const*>> node_paths
    = geom.FindAllVolumePaths(vol_names);
  //std::cout << "Found " << node_paths.size() << " node paths." << std::endl;

  // collection of interesting volumes
  util::PositionInVolumeFilter::AllVolumeInfo_t GeoVolumePairs;
  GeoVolumePairs.reserve(node_paths.size()); // because we are obsessed

  //for each interesting volume, follow the node path and collect
  //total rotations and translations
  for (size_t iVolume = 0; iVolume < node_paths.size(); ++iVolume) {
    std::vector<TGeoNode const*> path = node_paths[iVolume];

    TGeoTranslation* pTransl = new TGeoTranslation(0.,0.,0.);
    TGeoRotation* pRot = new TGeoRotation();
    for (TGeoNode const* node: path) {
      TGeoTranslation thistranslate(*node->GetMatrix());
      TGeoRotation thisrotate(*node->GetMatrix());
      pTransl->Add(&thistranslate);
      *pRot=*pRot * thisrotate;
    }

    //for some reason, pRot and pTransl don't have tr and rot bits set correctly
    //make new translations and rotations so bits are set correctly
    TGeoTranslation* pTransl2 = new TGeoTranslation(pTransl->GetTranslation()[0],
                                                    pTransl->GetTranslation()[1],
                                                    pTransl->GetTranslation()[2]);
    double phi=0.,theta=0.,psi=0.;
    pRot->GetAngles(phi,theta,psi);
    TGeoRotation* pRot2 = new TGeoRotation();
    pRot2->SetAngles(phi,theta,psi);

    TGeoCombiTrans* pTransf = new TGeoCombiTrans(*pTransl2,*pRot2);
    GeoVolumePairs.emplace_back(node_paths[iVolume].back()->GetVolume(), pTransf);

  }

  return std::make_unique<util::PositionInVolumeFilter>(std::move(GeoVolumePairs));

} // CreateParticleVolumeFilter()


void obana::OBAnaICARUS::clear_vectors() {

  _n_mct_lt1 = 0;
  _n_mct_lt1_from_ob = 0;
  _n_mcs_lt1 = 0;
  _n_mcs_lt1_from_ob = 0;


  _mcp_px.clear();
  _mcp_py.clear();
  _mcp_pz.clear();
  _mcp_e.clear();
  _mcp_vx.clear();
  _mcp_vy.clear();
  _mcp_vz.clear();
  _mcp_endx.clear();
  _mcp_endy.clear();
  _mcp_endz.clear();
  _mcp_pdg.clear();
  _mcp_mother.clear();
  _mcp_status_code.clear();
  _mcp_process.clear();
  _mcp_end_process.clear();
  _mcp_intpc.clear();
  _mcp_intpc_e.clear();
  _mcp_trackid.clear();
  _mcp_intpc_nu_e.clear();

  _neut_daughters_pdg.clear();
  _neut_daughters_e.clear();
  _neut_daughters_px.clear();
  _neut_daughters_py.clear();
  _neut_daughters_pz.clear();
  _neut_daughters_startprocess.clear();
  _neut_daughters_endprocess.clear();
  _neut_daughters_start_x.clear();
  _neut_daughters_start_y.clear();
  _neut_daughters_start_z.clear();
  _neut_daughters_end_x.clear();
  _neut_daughters_end_y.clear();
  _neut_daughters_end_z.clear();
  _neut_daughters_trackid.clear();

  _neut_grand_daughters_pdg.clear();
  _neut_grand_daughters_e.clear();
  _neut_grand_daughters_px.clear();
  _neut_grand_daughters_py.clear();
  _neut_grand_daughters_pz.clear();
  _neut_grand_daughters_startprocess.clear();
  _neut_grand_daughters_endprocess.clear();
  _neut_grand_daughters_start_x.clear();
  _neut_grand_daughters_start_y.clear();
  _neut_grand_daughters_start_z.clear();
  _neut_grand_daughters_end_x.clear();
  _neut_grand_daughters_end_y.clear();
  _neut_grand_daughters_end_z.clear();
  _neut_grand_daughters_trackid.clear();

  _mcs_pdg.clear();
  _mcs_process.clear();
  _mcs_start_x.clear();
  _mcs_start_y.clear();
  _mcs_start_z.clear();
  _mcs_end_x.clear();
  _mcs_end_y.clear();
  _mcs_end_z.clear();
  _mcs_start_px.clear();
  _mcs_start_py.clear();
  _mcs_start_pz.clear();
  _mcs_start_e.clear();
  _mcs_charge_col.clear();
  _mcs_charge_ind2.clear();
  _mcs_charge_ind1.clear();
  _mcs_mother_pdg.clear();
  _mcs_mother_trackid.clear();
  _mcs_mother_process.clear();
  _mcs_mother_start_x.clear();
  _mcs_mother_start_y.clear();
  _mcs_mother_start_z.clear();
  _mcs_mother_start_e.clear();
  _mcs_mother_end_x.clear();
  _mcs_mother_end_y.clear();
  _mcs_mother_end_z.clear();
  _mcs_mother_end_e.clear();
  _mcs_ancestor_pdg.clear();
  _mcs_ancestor_process.clear();
  _mcs_ancestor_start_e.clear();
  _mcs_ancestor_end_e.clear();
  _mcs_mother_in_ob_trackid.clear();


  _mct_pdg.clear();
  _mct_process.clear();
  _mct_start_x.clear();
  _mct_start_y.clear();
  _mct_start_z.clear();
  _mct_end_x.clear();
  _mct_end_y.clear();
  _mct_end_z.clear();
  _mct_start_px.clear();
  _mct_start_py.clear();
  _mct_start_pz.clear();
  _mct_start_e.clear();
  _mct_mother_pdg.clear();
  _mct_mother_process.clear();
  _mct_ancestor_pdg.clear();
  _mct_ancestor_process.clear();
  _mct_mother_in_ob_trackid.clear();


  _pi0_ids.clear();
  _muon_ids.clear();
}

bool obana::OBAnaICARUS::InDetector(const double& x,
				    const double& y,
				    const double& z) const {

  return (((x > xminc0 && x < xmaxc0)
	   || (x > xminc1 && x < xmaxc1))
	  && y > ymin   && y < ymax
	  && z > zmin   && z < zmax);
}

/*
bool obana::OBAnaICARUS::InDetector(art::Ptr<simb::MCTruth> mctruth){
  const geo::GeometryCore*           fGeometry;             ///< pointer to Geometry service
  fGeometry           = lar::providerFrom<geo::Geometry>();
  for (size_t i = 0; i < mctruth.size(); i++) {
  geo::Point_t trackPoint(pos.X(),pos.Y(),pos.Z());

    const geo::TPCGeo* tpcGeo = fGeometry->PositionToTPCptr(trackPoint);
}
*/

bool obana::OBAnaICARUS::InDetector(simb::MCParticle const& mcp, int & step) {
  auto t = mcp.Trajectory();
  //  std::cout<< "size of trajectory: "<< mcp.NumberTrajectoryPoints() << "\t" << t.size() << std::endl;
  const geo::GeometryCore*           fGeometry;             ///< pointer to Geometry service
  fGeometry           = lar::providerFrom<geo::Geometry>();

  for (size_t i = 0; i < t.size(); i++) {
    // std::cout << "step: " << i << " , has energy : " << mcp.E(i) << "\t" <<mcp.PdgCode()<< std::endl;
    const TLorentzVector& pos = mcp.Position(i);
    geo::Point_t trackPoint(pos.X(),pos.Y(),pos.Z());

    const geo::TPCGeo* tpcGeo = fGeometry->PositionToTPCptr(trackPoint);

    if (!tpcGeo) continue;
    //    if (tpcGeo->ID() != C0id) continue; // point not in cryostat 0
    if (tpcGeo->ActiveBoundingBox().ContainsPosition(trackPoint)){
      step = i;
      //      std::cout << "Indide the boundary, step: " << step << " , has energy : " << mcp.E(step) << std::endl;
      return true; 
    }
    // if (InDetector(t.X(i), t.Y(i), t.Z(i))) return true;
  }
  return false;
}


void obana::OBAnaICARUS::beginSubRun(art::SubRun const& sr) {

  _sr_run       = sr.run();
  _sr_subrun    = sr.subRun();
  _sr_begintime = sr.beginTime().value();
  _sr_endtime   = sr.endTime().value();

  art::Handle<sumdata::POTSummary> pot_handle;
  sr.getByLabel(_mctruth_producer, pot_handle);

  if (pot_handle.isValid()) {
    _sr_pot = pot_handle->totpot;
  } else {
    _sr_pot = 0.;
  }
  //std::cout << "POT for this subrun: " << _sr_pot << std::endl;

  _sr_tree->Fill();
  //std::cout << "POT for this subrun, after filling: " << _sr_pot << std::endl;
}


DEFINE_ART_MODULE(obana::OBAnaICARUS)
