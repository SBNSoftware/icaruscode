/**
 * @file icaruscode/TPC/NuGraph/ICARUSNuGraphInference_module.cc
 * @brief Implementation of `ICARUSNuGraphInference` _art_ module.
 * @author Leonardo Lena (https://github.com/leonardo-lena) based on previous work.
 * @date October 1, 2025
 */

#include <torch/script.h> // this is to be loaded first else it conflicts with... something and does not compile.

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "lardataobj/RecoBase/Slice.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <canvas/Persistency/Common/Ptr.h>
#include <memory>
#include <map>

#include "delaunator-header-only.hpp"
#include <vector>

#include "lardataobj/RecoBase/Hit.h"

#include "larrecodnn/NuGraph/Tools/DecoderToolBase.h"
#include "larrecodnn/NuGraph/Tools/LoaderToolBase.h"

using recob::Slice;
using recob::Hit;
using std::vector;

/**
 * @class ICARUSNuGraphInference
 * Takes a std::vector<art::Ptr<recob::Slice>> from the Event on which to run NuGraph inference.
 * The data is pre-processed using a separate LoaderTool and then Delauney graphs are constructed in this module.
 * NuGraph is then run on every single slice one at the time and the outputs are then flattened together to keep backward compatibility.
 * Two separate DecoderTool(s) - one for filter, one for semantic - are responsible for taking the NuGraph output and writing it to the Event.
 */
class ICARUSNuGraphInference : public art::EDProducer {
  public:
    explicit ICARUSNuGraphInference(const fhicl::ParameterSet& params);

    // Plugins should not be copied or assigned.
    ICARUSNuGraphInference(ICARUSNuGraphInference const&) = delete;
    ICARUSNuGraphInference(ICARUSNuGraphInference&&) = delete;
    ICARUSNuGraphInference& operator=(ICARUSNuGraphInference const&) = delete;
    ICARUSNuGraphInference& operator=(ICARUSNuGraphInference&&) = delete;

    // Required functions.
    void produce(art::Event& event) override;

  private:
  art::InputTag fSlicesLabel;
  art::InputTag fHitLabel;
  static constexpr double minChiSq = 0.5;
  vector<std::string> planes;
  size_t minHits;
  bool debug;
  vector<vector<float>> avgs;
  vector<vector<float>> devs;
  vector<float> pos_norm;
  torch::jit::script::Module model;
  std::unique_ptr<LoaderToolBase> _loaderTool;
  // decoder tools
  std::vector<std::unique_ptr<DecoderToolBase>> _decoderToolsVec;
};

ICARUSNuGraphInference::ICARUSNuGraphInference(const fhicl::ParameterSet& params) : art::EDProducer{params},
  fSlicesLabel(params.get<art::InputTag>("SlicesLabel", "NCCSlicesCryoE")),
  fHitLabel(params.get<art::InputTag>("HitLabel", "pandoraGausCryoE")),
  planes(params.get<vector<std::string>>("planes")),
  minHits(params.get<size_t>("minHits")),
  debug(params.get<bool>("debug", false)),
  pos_norm(params.get<vector<float>>("pos_norm")) {
    for (size_t ip = 0; ip < planes.size(); ++ip) {
      avgs.push_back(params.get<vector<float>>("avgs_" + planes[ip]));
      devs.push_back(params.get<vector<float>>("devs_" + planes[ip]));
    }

    _loaderTool = art::make_tool<LoaderToolBase>(params.get<fhicl::ParameterSet>("LoaderTool"));
    _loaderTool->setDebugAndPlanes(debug, planes);

    // configure and construct Decoder Tools
    auto const tool_psets = params.get<fhicl::ParameterSet>("DecoderTools");
    for (auto const& tool_pset_labels : tool_psets.get_pset_names()) {
      std::cout << "decoder label: " << tool_pset_labels << std::endl;
      auto const tool_pset = tool_psets.get<fhicl::ParameterSet>(tool_pset_labels);
      _decoderToolsVec.push_back(art::make_tool<DecoderToolBase>(tool_pset));
      _decoderToolsVec.back()->setDebugAndPlanes(debug, planes);
      _decoderToolsVec.back()->declareProducts(producesCollector());
    }

    cet::search_path sp("FW_SEARCH_PATH");
    model = torch::jit::load(sp.find_file(params.get<std::string>("modelFileName")));
  }

void ICARUSNuGraphInference::produce(art::Event& event) {
  const std::vector<art::Ptr<recob::Slice>> slices = event.getProduct<vector<art::Ptr<recob::Slice>>>(fSlicesLabel);

  art::FindManyP<Hit> findMHitsFromSlice(slices, event, fHitLabel);

  vector<vector<size_t>> idsmap = std::vector<std::vector<size_t>>(planes.size(), std::vector<size_t>());
  std::map<std::string, std::vector<float>> outputsMap;
  vector<NuGraphOutput> infer_output;

  for (size_t sliceIdx = 0; sliceIdx < slices.size(); ++sliceIdx) {
    if (debug) std::cout << "===== BEGIN SLICE " << sliceIdx+1 << "/" << slices.size() << " =====" << '\n';

    vector<NuGraphInput> graphinputs;
    vector<vector<size_t>> singleIdsmap = std::vector<std::vector<size_t>>(planes.size(), std::vector<size_t>());

    vector<art::Ptr<Hit>> hitsInSlice = findMHitsFromSlice.at(sliceIdx);

    _loaderTool->loadData(event, hitsInSlice, graphinputs, singleIdsmap);

    bool emptyPlane = false;

    for (size_t plane = 0; plane < singleIdsmap.size(); plane++) {
      if (singleIdsmap[plane].size() <= 0) {
        emptyPlane = true;
      }
    }

    if (emptyPlane) {
      continue;
      if (debug) std::cout << "Slice has atleast one empty plane. Skipping." << std::endl;
    } else {
      for (size_t plane = 0; plane < singleIdsmap.size(); plane++) {
        idsmap[plane].insert(idsmap[plane].end(), singleIdsmap[plane].begin(), singleIdsmap[plane].end());
      }
    }

    const vector<int32_t>* spids = nullptr;
    const vector<int32_t>* hitids_u = nullptr;
    const vector<int32_t>* hitids_v = nullptr;
    const vector<int32_t>* hitids_y = nullptr;
    const vector<int32_t>* hit_plane = nullptr;
    const vector<float>* hit_time = nullptr;
    const vector<int32_t>* hit_wire = nullptr;
    const vector<float>* hit_integral = nullptr;
    const vector<float>* hit_rms = nullptr;
    for (const auto& gi : graphinputs) {
      if (gi.input_name == "spacepoint_table_spacepoint_id")
        spids = &gi.input_int32_vec;
      else if (gi.input_name == "spacepoint_table_hit_id_u")
        hitids_u = &gi.input_int32_vec;
      else if (gi.input_name == "spacepoint_table_hit_id_v")
        hitids_v = &gi.input_int32_vec;
      else if (gi.input_name == "spacepoint_table_hit_id_y")
        hitids_y = &gi.input_int32_vec;
      else if (gi.input_name == "hit_table_local_plane")
        hit_plane = &gi.input_int32_vec;
      else if (gi.input_name == "hit_table_local_time")
        hit_time = &gi.input_float_vec;
      else if (gi.input_name == "hit_table_local_wire")
        hit_wire = &gi.input_int32_vec;
      else if (gi.input_name == "hit_table_integral")
        hit_integral = &gi.input_float_vec;
      else if (gi.input_name == "hit_table_rms")
        hit_rms = &gi.input_float_vec;
    }

    if (debug) std::cout << "begin idsmapRev creation" << '\n';
    // Reverse lookup from key to index in plane index
    std::map<size_t, size_t> idsmapRev;
    for (const auto& ipv : singleIdsmap) {
      for (size_t ih = 0; ih < ipv.size(); ih++) {
        idsmapRev.insert(std::make_pair(ipv.at(ih), ih));
      }
    }
    if (debug) std::cout << "end idsmapRev creation" << '\n';


    struct Edge {
      size_t n1;
      size_t n2;
      bool operator==(const Edge& other) const
      {
        if (this->n1 == other.n1 && this->n2 == other.n2)
          return true;
        else
          return false;
      };
    };

    // Delauney graph construction
    auto start_preprocess1 = std::chrono::high_resolution_clock::now();
    vector<vector<Edge>> edge2d(planes.size(), vector<Edge>());
    for (size_t p = 0; p < planes.size(); p++) {
      vector<double> coords;
      for (size_t i = 0; i < hit_plane->size(); ++i) {
        if (size_t(hit_plane->at(i)) != p) continue;
        coords.push_back(hit_time->at(i) * pos_norm.at(1));
        coords.push_back(hit_wire->at(i) * pos_norm.at(0));
        if (debug) {
          std::cout 
          << "Plane: " << p
          << " Hit: " << i 
          << " HitIdxInEdge: " << (coords.size() / 2) - 1
          << " time_g: " << hit_time->at(i)
          << " wire_g: " << hit_wire->at(i)
          << '\n';
        }
      }
      if (debug) std::cout << "Plane " << p << " has N hits=" << coords.size() / 2 << std::endl;
      if (coords.size() / 2 < 3) { continue; }
      delaunator::Delaunator d(coords);
      if (debug) std::cout << "Found N triangles=" << d.triangles.size() / 3 << std::endl;
      for (std::size_t i = 0; i < d.triangles.size(); i += 3) {
        //create edges in both directions
        Edge e;
        e.n1 = d.triangles.at(i);
        e.n2 = d.triangles.at(i+1);
        edge2d[p].push_back(e);
        e.n1 = d.triangles.at(i+1);
        e.n2 = d.triangles.at(i);
        edge2d[p].push_back(e);
        //
        e.n1 = d.triangles.at(i);
        e.n2 = d.triangles.at(i+2);
        edge2d[p].push_back(e);
        e.n1 = d.triangles.at(i+1);
        e.n2 = d.triangles.at(i);
        edge2d[p].push_back(e);
        //
        e.n1 = d.triangles.at(i+1);
        e.n2 = d.triangles.at(i+2);
        edge2d[p].push_back(e);
        e.n1 = d.triangles.at(i+2);
        e.n2 = d.triangles.at(i+1);
        edge2d[p].push_back(e);
        //
      }
      //sort and cleanup duplicate edges
      std::sort(edge2d[p].begin(), edge2d[p].end(), [](const auto& i, const auto& j) {
        return (i.n1 != j.n1 ? i.n1 < j.n1 : i.n2 < j.n2);
      });
      if (debug) {
        for (auto& e : edge2d[p]) {
          std::cout << "sorted plane=" << p << " e1=" << e.n1 << " e2=" << e.n2 << std::endl;
        }
      }
      edge2d.at(p).erase(std::unique(edge2d.at(p).begin(), edge2d.at(p).end()), edge2d.at(p).end());
    }

    if (debug) {
      for (size_t p = 0; p < planes.size(); p++) {
        for (auto& e : edge2d[p]) {
          std::cout << " plane=" << p << " e1=" << e.n1 << " e2=" << e.n2 << std::endl;
        }
      }
    }
    auto end_preprocess1 = std::chrono::high_resolution_clock::now();
    if (debug) std::cout << "end preprocess1" << '\n';
    std::chrono::duration<double> elapsed_preprocess1 = end_preprocess1 - start_preprocess1;

    // Nexus edges
    auto start_preprocess2 = std::chrono::high_resolution_clock::now();
    vector<vector<Edge>> edge3d(planes.size(), vector<Edge>());
    if (debug) std::cout << "start 3d edges creation" << '\n';
    for (size_t i = 0; i < spids->size(); ++i) {
      if (hitids_u->at(i) >= 0) {
        Edge e;
        e.n1 = idsmapRev[hitids_u->at(i)];
        e.n2 = i;
        edge3d[0].push_back(e);
      }
      if (hitids_v->at(i) >= 0) {
        Edge e;
        e.n1 = idsmapRev[hitids_v->at(i)];
        e.n2 = i;
        edge3d[1].push_back(e);
      }
      if (hitids_y->at(i) >= 0) {
        Edge e;
        e.n1 = idsmapRev[hitids_y->at(i)];
        e.n2 = i;
        edge3d[2].push_back(e);
      }
    }
    if (debug) std::cout << "end 3d edges creation" << '\n';

    // Prepare inputs
    auto x = torch::Dict<std::string, torch::Tensor>();
    auto batch = torch::Dict<std::string, torch::Tensor>();
    for (size_t p = 0; p < planes.size(); p++) {
      vector<float> nodeft;
      for (size_t i = 0; i < hit_plane->size(); ++i) {
        if (size_t(hit_plane->at(i)) != p) continue;
        nodeft.push_back((hit_wire->at(i) * pos_norm[0] - avgs[hit_plane->at(i)][0]) /
                         devs[hit_plane->at(i)][0]);
        nodeft.push_back((hit_time->at(i) * pos_norm[1] - avgs[hit_plane->at(i)][1]) /
                         devs[hit_plane->at(i)][1]);
        nodeft.push_back((hit_integral->at(i) - avgs[hit_plane->at(i)][2]) /
                         devs[hit_plane->at(i)][2]);
        nodeft.push_back((hit_rms->at(i) - avgs[hit_plane->at(i)][3]) / devs[hit_plane->at(i)][3]);
      }
      long int dim = nodeft.size() / 4;
      torch::Tensor ix = torch::zeros({dim, 4}, torch::dtype(torch::kFloat32));
      if (debug) {
        std::cout << "plane=" << p << std::endl;
        std::cout << std::scientific;
        for (size_t n = 0; n < nodeft.size(); n = n + 4) {
          std::cout << nodeft[n] << " " << nodeft[n + 1] << " " << nodeft[n + 2] << " "
                    << nodeft[n + 3] << " " << std::endl;
        }
      }
      for (size_t n = 0; n < nodeft.size(); n = n + 4) {
        ix[n / 4][0] = nodeft[n];
        ix[n / 4][1] = nodeft[n + 1];
        ix[n / 4][2] = nodeft[n + 2];
        ix[n / 4][3] = nodeft[n + 3];
      }
      x.insert(planes[p], ix);
      torch::Tensor ib = torch::zeros({dim}, torch::dtype(torch::kInt64));
      batch.insert(planes[p], ib);
    }

    auto edge_index_plane = torch::Dict<std::string, torch::Tensor>();
    for (size_t p = 0; p < planes.size(); p++) {
      long int dim = edge2d[p].size();
      torch::Tensor ix = torch::zeros({2, dim}, torch::dtype(torch::kInt64));
      for (size_t n = 0; n < edge2d[p].size(); n++) {
        ix[0][n] = int(edge2d[p][n].n1);
        ix[1][n] = int(edge2d[p][n].n2);
      }
      edge_index_plane.insert(planes[p], ix);
      if (debug) {
        std::cout << "plane=" << p << std::endl;
        std::cout << "2d edge size=" << edge2d[p].size() << std::endl;
        for (size_t n = 0; n < edge2d[p].size(); n++) {
          std::cout << edge2d[p][n].n1 << " ";
        }
        std::cout << std::endl;
        for (size_t n = 0; n < edge2d[p].size(); n++) {
          std::cout << edge2d[p][n].n2 << " ";
        }
        std::cout << std::endl;
      }
    }

    auto edge_index_nexus = torch::Dict<std::string, torch::Tensor>();
    for (size_t p = 0; p < planes.size(); p++) {
      long int dim = edge3d[p].size();
      torch::Tensor ix = torch::zeros({2, dim}, torch::dtype(torch::kInt64));
      for (size_t n = 0; n < edge3d[p].size(); n++) {
        ix[0][n] = int(edge3d[p][n].n1);
        ix[1][n] = int(edge3d[p][n].n2);
      }
      edge_index_nexus.insert(planes[p], ix);
      if (debug) {
        std::cout << "plane=" << p << std::endl;
        std::cout << "3d edge size=" << edge3d[p].size() << std::endl;
        for (size_t n = 0; n < edge3d[p].size(); n++) {
          std::cout << edge3d[p][n].n1 << " ";
        }
        std::cout << std::endl;
        for (size_t n = 0; n < edge3d[p].size(); n++) {
          std::cout << edge3d[p][n].n2 << " ";
        }
        std::cout << std::endl;
      }
    }

    long int spdim = spids->size();
    auto nexus = torch::empty({spdim, 0}, torch::dtype(torch::kFloat32));

    std::vector<torch::jit::IValue> inputs;
    inputs.push_back(x);
    inputs.push_back(edge_index_plane);
    inputs.push_back(edge_index_nexus);
    inputs.push_back(nexus);
    inputs.push_back(batch);

    // Run inference
    auto end_preprocess2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_preprocess2 = end_preprocess2 - start_preprocess2;
    if (debug) std::cout << "FORWARD!" << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    c10::IValue raw_outputs;

    {
      // Use NoGradGuard to disable gradient tracking.
      // Gradients are not needed for inference, so save some memory and computation
      // by disabling them.
      //
      // This is equivalent to `with torch.no_grad():` in Python.
      // Uses RAII, so is disabled once the guard goes out of scope.
      torch::NoGradGuard guard;
      raw_outputs = model.forward(inputs);
    }

    auto outputs = raw_outputs.toGenericDict();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    if (debug) {
      std::cout << "Time taken for inference: "
                << elapsed_preprocess1.count() + elapsed_preprocess2.count() + elapsed.count()
                << " seconds" << std::endl;
      std::cout << "output =" << outputs << std::endl;
    }

    for (const auto& elem1 : outputs) {
      if (elem1.value().isTensor()) {
        torch::Tensor tensor = elem1.value().toTensor();
        std::vector<float> vec(tensor.data_ptr<float>(), tensor.data_ptr<float>() + tensor.numel());
        std::string key = elem1.key().to<std::string>();
        auto it = std::find_if(infer_output.begin(), infer_output.end(), [key](const NuGraphOutput& output){return output.output_name == key;});
        if (it == infer_output.end()) {
          infer_output.emplace_back(NuGraphOutput(key, vec));
        }
        else {
          it->output_vec.insert(it->output_vec.end(), vec.begin(), vec.end());
        }
      }
      else if (elem1.value().isGenericDict()) {
        for (const auto& elem2 : elem1.value().toGenericDict()) {
          torch::Tensor tensor = elem2.value().toTensor();
          std::vector<float> vec(tensor.data_ptr<float>(), tensor.data_ptr<float>() + tensor.numel());
          std::string key = elem1.key().to<std::string>() + "_" + elem2.key().to<std::string>();
          auto it = std::find_if(infer_output.begin(), infer_output.end(), [key](const NuGraphOutput& output){return output.output_name == key;});
          if (it == infer_output.end()) {
            infer_output.emplace_back(NuGraphOutput(key, vec));
          }
          else {
            it->output_vec.insert(it->output_vec.end(), vec.begin(), vec.end());
          }
        }
      }
    }
  }

  size_t idsmapEntries = 0;
  for (std::vector<size_t> idsVec : idsmap) {
    for (size_t idCounter = 0; idCounter < idsVec.size(); ++idCounter) {
      idsmapEntries++;
    }
  }
  if (debug) {
    std::cout << "idsmap size: " << idsmapEntries;
    for (const NuGraphOutput& output : infer_output) {
      std::cout << output.output_name << " size: " << output.output_vec.size() << ' ';
    }
    std::cout << std::endl;
  }

  if (idsmapEntries == 0) {
    for (size_t i = 0; i < _decoderToolsVec.size(); i++) {
      _decoderToolsVec[i]->writeEmptyToEvent(event, idsmap);
    }
  } else {
    for (size_t i = 0; i < _decoderToolsVec.size(); i++) {
      _decoderToolsVec[i]->writeToEvent(event, idsmap, infer_output);
    }
  }
  if (debug) std::cout << "INFERENCE COMPLETE" << '\n';
}

DEFINE_ART_MODULE(ICARUSNuGraphInference)
