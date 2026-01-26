/**
 * @file   icaruscode/TPC/NuGraph/ICARUSNuGraphMultiLoader_tool.cc
 * @author Leonardo Lena (https://github.com/leonardo-lena) based on G. Cerati (cerati@fnal.gov) work.
 */

#include "larrecodnn/NuGraph/Tools/LoaderToolBase.h"

#include "art/Utilities/ToolMacros.h"

#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Slice.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h" // geo::WireID
#include <utility> // std::move()
#include <torch/torch.h>
#include <vector>

#include "StitchingUtils.h"

class ICARUSNuGraphMultiLoader : public LoaderToolBase {

/**
 * @brief Tool to collect the inputs needed by NuGraph from ICARUS data products.
 *
 * Read hit and space points from the event record, and package them for usage in NuGraph.
 * This tool is called from icaruscode/TPC/NuGraph/ICARUSNuGraphInference_module.cc.
 * Hits are stitched so that the 4 ICARUS TPCs in a cryostat are viewed in a single time vs wire plane.
 * Only space points with chi2<MinChiSq (currently set to 0.5) and with 3 associated hits are considered.
 * The module assumes the input hits to all be from the same slice.
 */

public:
  /**
   *  @brief  Constructor
   *
   *  @param  pset
   */
  ICARUSNuGraphMultiLoader(const fhicl::ParameterSet& params);

  /**
   * @brief loadData function
   *
   * @param art::Event event record, list of input, idsmap
   */
  void loadData(art::Event& event,
                vector<art::Ptr<recob::Hit>>& hitsInSlice,
                vector<NuGraphInput>& graphinputs,
                vector<vector<size_t>>& singleIdsmap) override;

private:
  static constexpr double minChiSq = 0.5; ///< Threshold to consider a space point good.
  art::InputTag fPandoraLabel;
  art::InputTag fSlicesLabel;
  art::InputTag fClusterLabel;
};

ICARUSNuGraphMultiLoader::ICARUSNuGraphMultiLoader(const fhicl::ParameterSet& params)
  : fPandoraLabel(params.get<art::InputTag>("PandoraLabel", "pandoraGausCryoE"))
  , fSlicesLabel(params.get<art::InputTag>("SlicesLabel", "NCCSlicesCryoE"))
  , fClusterLabel(params.get<art::InputTag>("ClusterLabel", "cluster3DCryoE")) {}

void ICARUSNuGraphMultiLoader::loadData(art::Event& event,
                              std::vector<art::Ptr<recob::Hit>>& hitsInSlice,
                              std::vector<NuGraphInput>& graphinputs,
                              std::vector<std::vector<size_t>>& singleIdsmap)
{
  std::vector<int32_t> hit_table_hit_id_data;
  std::vector<int32_t> hit_table_local_plane_data;
  std::vector<float> hit_table_local_time_data;
  std::vector<int32_t> hit_table_local_wire_data;
  std::vector<float> hit_table_integral_data;
  std::vector<float> hit_table_rms_data;
  std::vector<int32_t> spacepoint_table_spacepoint_id_data;
  std::vector<int32_t> spacepoint_table_hit_id_u_data;
  std::vector<int32_t> spacepoint_table_hit_id_v_data;
  std::vector<int32_t> spacepoint_table_hit_id_y_data;

  art::ValidHandle<std::vector<recob::SpacePoint>> spacePointsHandle = event.getValidHandle<std::vector<recob::SpacePoint>>(fClusterLabel);

  art::FindOneP<recob::Slice> find1SliceFromHit(hitsInSlice, event, fPandoraLabel);
  art::Ptr<recob::Slice> slice = find1SliceFromHit.at(0); // we assume the hits are all from the same slice.

  art::FindManyP<recob::SpacePoint> findMSpsFromSlice(std::vector<art::Ptr<recob::Slice>>{slice}, event, fSlicesLabel);
  std::vector<art::Ptr<recob::SpacePoint>> spacePointsInSlice = findMSpsFromSlice.at(0);

  art::FindManyP<recob::Hit> findMHitsFromSps(spacePointsHandle, event, fClusterLabel);

  for (const art::Ptr<recob::SpacePoint>&  spacePoint : spacePointsInSlice) {
    std::vector<art::Ptr<recob::Hit>> hitsOfSpacePoint = findMHitsFromSps.at(spacePoint.key());

    if (spacePoint->Chisq() > minChiSq || hitsOfSpacePoint.size() != 3) continue;
    spacepoint_table_spacepoint_id_data.push_back(spacePoint.key());

    for (const art::Ptr<recob::Hit>& hit : hitsOfSpacePoint) {
      geo::WireID wireId = hit->WireID();
      size_t plane = util::stitchedPlane(wireId);
      if (plane == 0) spacepoint_table_hit_id_u_data.push_back(hit.key());
      if (plane == 1) spacepoint_table_hit_id_v_data.push_back(hit.key());
      if (plane == 2) spacepoint_table_hit_id_y_data.push_back(hit.key());
    }
  }

  std::set<std::tuple<size_t, size_t, double>> coordsSet;

  for (const art::Ptr<recob::Hit>& hit : hitsInSlice) {
    geo::WireID wireId = hit->WireID();
    size_t plane = util::stitchedPlane(wireId);
    double time = util::stitchedTime(wireId, hit->PeakTime());
    size_t wire = util::stitchedWire(wireId);

    auto insertReturn = coordsSet.insert(std::make_tuple(plane, wire, time));

    if (insertReturn.second) {
      hit_table_hit_id_data.push_back(hit.key());
      hit_table_local_plane_data.push_back(plane);
      hit_table_local_time_data.push_back(time);
      hit_table_local_wire_data.push_back(wire);
      hit_table_integral_data.push_back(hit->Integral());
      hit_table_rms_data.push_back(hit->RMS());

      singleIdsmap[plane].push_back(hit.key());
    } else if (debug) {
      std::cout << "Overlapping hit found at (plane, wire, time): (" << plane << ", " << wire << ", " << time << ").\n";
    }
    
  }

  graphinputs.emplace_back("hit_table_hit_id", std::move(hit_table_hit_id_data));
  graphinputs.emplace_back("hit_table_local_plane", std::move(hit_table_local_plane_data));
  graphinputs.emplace_back("hit_table_local_time", std::move(hit_table_local_time_data));
  graphinputs.emplace_back("hit_table_local_wire", std::move(hit_table_local_wire_data));
  graphinputs.emplace_back("hit_table_integral", std::move(hit_table_integral_data));
  graphinputs.emplace_back("hit_table_rms", std::move(hit_table_rms_data));

  graphinputs.emplace_back("spacepoint_table_spacepoint_id", std::move(spacepoint_table_spacepoint_id_data));
  graphinputs.emplace_back("spacepoint_table_hit_id_u", std::move(spacepoint_table_hit_id_u_data));
  graphinputs.emplace_back("spacepoint_table_hit_id_v", std::move(spacepoint_table_hit_id_v_data));
  graphinputs.emplace_back("spacepoint_table_hit_id_y", std::move(spacepoint_table_hit_id_y_data));
}

DEFINE_ART_CLASS_TOOL(ICARUSNuGraphMultiLoader)
