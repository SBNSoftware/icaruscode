/**
 * @file   icaruscode/TPC/NuGraph/ICARUSNuGraphLoader_tool.cc
 * @author Giuseppe Cerati (cerati@fnal.gov)
 */

#include "larrecodnn/NuGraph/Tools/LoaderToolBase.h"

#include "art/Utilities/ToolMacros.h"

#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h" // geo::WireID
#include <utility> // std::move()
#include <torch/torch.h>

#include "StitchingUtils.h"

class ICARUSNuGraphLoader : public LoaderToolBase {

/**
 * @brief Tool to collect the inputs needed by NuGraph from ICARUS data products.
 *
 * Read hit and space points from the event record, and package them for usage in NuGraph. This tool is called from larrecodnn/NuGraph/NuGraphInference_module.cc.
 * Hits are stitched so that the 4 ICARUS TPCs in a cryostat are viewed in a single time vs wire plane.
 * Only space points with chi2<MinChiSq (currently set to 0.5) are considered.
 *
 */

public:
  /**
   *  @brief  Constructor
   *
   *  @param  pset
   */
  ICARUSNuGraphLoader(const fhicl::ParameterSet& pset);

  /**
   * @brief loadData function
   *
   * @param art::Event event record, list of input, idsmap
   */
  void loadData(art::Event& e,
                vector<art::Ptr<recob::Hit>>& hitlist,
                vector<NuGraphInput>& inputs,
                vector<vector<size_t>>& idsmap) override;

private:

  art::InputTag const fHitInput;
  art::InputTag const fSpsInput;
  static constexpr double MinChiSq = 0.5; ///< Threshold to consider a space point good.
};

ICARUSNuGraphLoader::ICARUSNuGraphLoader(const fhicl::ParameterSet& p)
  : fHitInput{p.get<art::InputTag>("hitInput")}, fSpsInput{p.get<art::InputTag>("spsInput")}
{}

void ICARUSNuGraphLoader::loadData(art::Event& e,
                              vector<art::Ptr<recob::Hit>>& hitlist,
                              vector<NuGraphInput>& inputs,
                              vector<vector<size_t>>& idsmap)
{
  //
  auto hitListHandle = e.getValidHandle<vector<recob::Hit>>(fHitInput);
  art::fill_ptr_vector(hitlist, hitListHandle);
  //
  idsmap = std::vector<std::vector<size_t>>(planes.size(), std::vector<size_t>());
  for (auto h : hitlist) {
    //
    size_t plane = util::stitchedPlane(h->WireID());
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
    size_t plane = util::stitchedPlane(wireid);
    double time = util::stitchedTime(wireid,h->PeakTime());
    size_t wire = util::stitchedWire(wireid);

    hit_table_hit_id_data.push_back(h.key());
    hit_table_local_plane_data.push_back(plane);
    hit_table_local_time_data.push_back(time);
    hit_table_local_wire_data.push_back(wire);
    hit_table_integral_data.push_back(h->Integral());
    hit_table_rms_data.push_back(h->RMS());
  }
  mf::LogDebug{ "ICARUSNuGraphLoader" } << "loader has nhits=" << hit_table_hit_id_data.size();

  // Get spacepoints from the event record
  art::Handle<vector<recob::SpacePoint>> spListHandle;
  vector<art::Ptr<recob::SpacePoint>> splist;
  if (e.getByLabel(fSpsInput, spListHandle)) { art::fill_ptr_vector(splist, spListHandle); }
  // Get assocations from spacepoints to hits
  vector<vector<art::Ptr<recob::Hit>>> sp2Hit(splist.size());
  if (!splist.empty()) {
    art::FindManyP<recob::Hit> fmp(spListHandle, e, fSpsInput);
    for (size_t spIdx = 0; spIdx < sp2Hit.size(); ++spIdx) {
      if (splist[spIdx]->Chisq()>MinChiSq) continue;
      // only space points with hits on all planes are enough for NuGraph
      if (fmp.at(spIdx).size()!=3) continue;
      sp2Hit[spIdx] = fmp.at(spIdx);
    }
  }

  // space point table
  size_t spidx = 0;
  for (size_t i = 0; i < splist.size(); ++i) {
    if (sp2Hit[i].empty()) continue;
    spacepoint_table_spacepoint_id_data.push_back(spidx++);
    spacepoint_table_hit_id_u_data.push_back(-1);
    spacepoint_table_hit_id_v_data.push_back(-1);
    spacepoint_table_hit_id_y_data.push_back(-1);
    for (size_t j = 0; j < sp2Hit[i].size(); ++j) {
      //
      size_t plane = util::stitchedPlane(sp2Hit[i][j]->WireID());
      if (plane == 0) spacepoint_table_hit_id_u_data.back() = sp2Hit[i][j].key();
      if (plane == 1) spacepoint_table_hit_id_v_data.back() = sp2Hit[i][j].key();
      if (plane == 2) spacepoint_table_hit_id_y_data.back() = sp2Hit[i][j].key();
    }
  }
  mf::LogDebug{ "ICARUSNuGraphLoader" }  << "loader has nsps=" << spacepoint_table_hit_id_u_data.size();

  inputs.emplace_back("hit_table_hit_id", std::move(hit_table_hit_id_data));
  inputs.emplace_back("hit_table_local_plane", std::move(hit_table_local_plane_data));
  inputs.emplace_back("hit_table_local_time", std::move(hit_table_local_time_data));
  inputs.emplace_back("hit_table_local_wire", std::move(hit_table_local_wire_data));
  inputs.emplace_back("hit_table_integral", std::move(hit_table_integral_data));
  inputs.emplace_back("hit_table_rms", std::move(hit_table_rms_data));

  inputs.emplace_back("spacepoint_table_spacepoint_id", std::move(spacepoint_table_spacepoint_id_data));
  inputs.emplace_back("spacepoint_table_hit_id_u", std::move(spacepoint_table_hit_id_u_data));
  inputs.emplace_back("spacepoint_table_hit_id_v", std::move(spacepoint_table_hit_id_v_data));
  inputs.emplace_back("spacepoint_table_hit_id_y", std::move(spacepoint_table_hit_id_y_data));
}
DEFINE_ART_CLASS_TOOL(ICARUSNuGraphLoader)
