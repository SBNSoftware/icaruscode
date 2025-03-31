#include "larrecodnn/NuGraph/Tools/LoaderToolBase.h"

#include "art/Utilities/ToolMacros.h"

#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include <torch/torch.h>

class IcarusNuGraphLoader : public LoaderToolBase {

public:
  /**
   *  @brief  Constructor
   *
   *  @param  pset
   */
  IcarusNuGraphLoader(const fhicl::ParameterSet& pset);

  /**
   *  @brief  Virtual Destructor
   */
  virtual ~IcarusNuGraphLoader() noexcept = default;

  /**
   * @brief loadData function
   *
   * @param art::Event event record, list of input, idsmap
   */
  void loadData(art::Event& e,
                vector<art::Ptr<recob::Hit>>& hitlist,
                vector<NuGraphInput>& inputs,
                vector<vector<size_t>>& idsmap) override;

  int stitchedPlane(const geo::WireID wid) const {
    int plane = wid.Plane;
    if(wid.TPC==2 || wid.TPC==3) {
      if(wid.Plane==1) plane=2;
      else if(wid.Plane==2) plane=1;
    }
    return plane;
  }
  float stitchedTime(const geo::WireID wid, float timein) const {
    float time = timein;
    if(wid.TPC==2 || wid.TPC==3) {
      //correction = 2*(tpcgeo.DriftDistance()/detProp.DriftVelocity()-clockData.TriggerOffsetTPC())/clockData.TPCClock().TickPeriod() = 6442.15
      time = 6442.15 - timein;
    }
    return time;
  }
  size_t stitchedWire(const geo::WireID wid) const {
    size_t wire = wid.Wire;
    int plane = stitchedPlane(wid);
    if(wid.TPC==1 || wid.TPC==3) {
      if(plane==1 || plane == 2) {
	wire = wid.Wire + 2535; //2535 is the last part of the wires in cryos 0 an 2 before the cut in z=0
      }
    }
    return wire;
  }

private:
  art::InputTag hitInput;
  art::InputTag spsInput;
};

IcarusNuGraphLoader::IcarusNuGraphLoader(const fhicl::ParameterSet& p)
  : hitInput{p.get<art::InputTag>("hitInput")}, spsInput{p.get<art::InputTag>("spsInput")}
{}

void IcarusNuGraphLoader::loadData(art::Event& e,
                              vector<art::Ptr<recob::Hit>>& hitlist,
                              vector<NuGraphInput>& inputs,
                              vector<vector<size_t>>& idsmap)
{
  //
  art::Handle<vector<recob::Hit>> hitListHandle;
  if (e.getByLabel(hitInput, hitListHandle)) { art::fill_ptr_vector(hitlist, hitListHandle); }
  //
  idsmap = std::vector<std::vector<size_t>>(planes.size(), std::vector<size_t>());
  for (auto h : hitlist) {
    //
    size_t plane = stitchedPlane(h->WireID());
    idsmap[plane].push_back(h.key());
  }

  vector<int32_t> hit_table_hit_id_data;
  vector<int32_t> hit_table_local_plane_data;
  vector<float> hit_table_local_time_data;
  vector<int32_t> hit_table_local_wire_data;
  vector<float> hit_table_integral_data;
  vector<float> hit_table_rms_data;
  vector<int32_t> spacepoint_table_spacepoint_id_data;
  vector<int32_t> spacepoint_table_hit_id_u_data;
  vector<int32_t> spacepoint_table_hit_id_v_data;
  vector<int32_t> spacepoint_table_hit_id_y_data;

  // hit table
  for (auto h : hitlist) {

    geo::WireID wireid = h->WireID();
    size_t plane = stitchedPlane(wireid);
    double time = stitchedTime(wireid,h->PeakTime());
    size_t wire = stitchedWire(wireid);

    hit_table_hit_id_data.push_back(h.key());
    hit_table_local_plane_data.push_back(plane);
    hit_table_local_time_data.push_back(time);
    hit_table_local_wire_data.push_back(wire);
    hit_table_integral_data.push_back(h->Integral());
    hit_table_rms_data.push_back(h->RMS());
  }
  std::cout << "loader has nhits=" << hit_table_hit_id_data.size() << std::endl;

  // Get spacepoints from the event record
  art::Handle<vector<recob::SpacePoint>> spListHandle;
  vector<art::Ptr<recob::SpacePoint>> splist;
  if (e.getByLabel(spsInput, spListHandle)) { art::fill_ptr_vector(splist, spListHandle); }
  // Get assocations from spacepoints to hits
  vector<vector<art::Ptr<recob::Hit>>> sp2Hit(splist.size());
  if (splist.size() > 0) {
    art::FindManyP<recob::Hit> fmp(spListHandle, e, spsInput);
    for (size_t spIdx = 0; spIdx < sp2Hit.size(); ++spIdx) {
      if (splist[spIdx]->Chisq()>0.5) continue;
      if (fmp.at(spIdx).size()!=3) continue;
      sp2Hit[spIdx] = fmp.at(spIdx);
    }
  }

  // space point table
  for (size_t i = 0; i < splist.size(); ++i) {
    if (sp2Hit[i].size()==0) continue;
    spacepoint_table_spacepoint_id_data.push_back(i);
    spacepoint_table_hit_id_u_data.push_back(-1);
    spacepoint_table_hit_id_v_data.push_back(-1);
    spacepoint_table_hit_id_y_data.push_back(-1);
    for (size_t j = 0; j < sp2Hit[i].size(); ++j) {
      //
      size_t plane = stitchedPlane(sp2Hit[i][j]->WireID());
      if (plane == 0) spacepoint_table_hit_id_u_data.back() = sp2Hit[i][j].key();
      if (plane == 1) spacepoint_table_hit_id_v_data.back() = sp2Hit[i][j].key();
      if (plane == 2) spacepoint_table_hit_id_y_data.back() = sp2Hit[i][j].key();
    }
  }
  std::cout << "loader has nsps=" << spacepoint_table_hit_id_u_data.size() << std::endl;

  inputs.emplace_back("hit_table_hit_id", hit_table_hit_id_data);
  inputs.emplace_back("hit_table_local_plane", hit_table_local_plane_data);
  inputs.emplace_back("hit_table_local_time", hit_table_local_time_data);
  inputs.emplace_back("hit_table_local_wire", hit_table_local_wire_data);
  inputs.emplace_back("hit_table_integral", hit_table_integral_data);
  inputs.emplace_back("hit_table_rms", hit_table_rms_data);

  inputs.emplace_back("spacepoint_table_spacepoint_id", spacepoint_table_spacepoint_id_data);
  inputs.emplace_back("spacepoint_table_hit_id_u", spacepoint_table_hit_id_u_data);
  inputs.emplace_back("spacepoint_table_hit_id_v", spacepoint_table_hit_id_v_data);
  inputs.emplace_back("spacepoint_table_hit_id_y", spacepoint_table_hit_id_y_data);
}
DEFINE_ART_CLASS_TOOL(IcarusNuGraphLoader)
