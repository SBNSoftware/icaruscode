////////////////////////////////////////////////////////////////////////
// \file    LArCVNEvaluator_module.cc
// \brief   Producer module creating CVN neural net results
// \author  Alexander Radovic - a.radovic@gmail.com
//          Saul Alonso Monsalve - saul.alonso.monsalve@cern.ch
////////////////////////////////////////////////////////////////////////

// C/C++ includes
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <zlib.h>

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Utilities/make_tool.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/Assns.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "icaruscode/ICARUSCVN/module_helpers/ICARUSITFNetHandler.h"
#include "icaruscode/ICARUSCVN/module_helpers/ICARUSPixelMapProducer.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Slice.h"
#include "larrecodnn/CVN/func/AssignLabels.h"
#include "larrecodnn/CVN/func/InteractionType.h"
#include "larrecodnn/CVN/func/PixelMap.h"
#include "larrecodnn/CVN/func/Result.h"

namespace lcvn {

class ICARUSCVNEvaluator : public art::EDProducer {
public:
  explicit ICARUSCVNEvaluator(fhicl::ParameterSet const &pset);

  void produce(art::Event &evt);
  std::vector<std::string> resolve_input_files(const std::string &pixelMapInpu);
  std::vector<std::string> load_from_folder(const std::string &folder_path);
  std::vector<unsigned char> decompress_pixelmap_gz(const std::string &filename);
  void load_and_predict(const int event,
                        std::unique_ptr<lcvn::ICARUSITFNetHandler> &handler,
                        const std::string &pixelMapInput);

private:
  // /// Module label for input pixel maps
  std::string fModule;
  std::string fPixelMapInput;
  std::string fPixelMapModuleLabel;
  art::InputTag fSliceLabel;
  art::InputTag fPFParticleModuleLabel;
  art::InputTag fT0Label;
  std::string fCVNType;
  bool fVerbose;
  std::unique_ptr<lcvn::ICARUSITFNetHandler> fTFHandler;
  ICARUSPixelMapHitProducer fPMProducer;
  // std::vector<EvalConfig> fEvaluators;
};

//.......................................................................
// ICARUSCVNEvaluator::ICARUSCVNEvaluator(fhicl::ParameterSet const &pset)
//     : EDProducer{pset} {
//   std::vector<fhicl::ParameterSet> evalPsets =
//       pset.get<std::vector<fhicl::ParameterSet>>("Evaluators");

//   for (auto const &epset : evalPsets) {
//     EvalConfig cfg{epset.get<std::string>("Module"),
//                    epset.get<std::string>("PixelMapInput"),
//                    epset.get<std::string>("PixelMapModuleLabel"),
//                    epset.get<art::InputTag>("SliceLabel", ""),
//                    epset.get<art::InputTag>("PFParticleModuleLabel", ""),
//                    epset.get<art::InputTag>("T0Label", ""),
//                    epset.get<std::string>("CVNType"),
//                    epset.get<bool>("verbose", false),
//                    art::make_tool<lcvn::ICARUSITFNetHandler>(
//                        epset.get<fhicl::ParameterSet>("ICARUSTFHandler")),
//                    ICARUSPixelMapHitProducer(
//                        epset.get<fhicl::ParameterSet>("PixelMapProducer"))};
//     fEvaluators.push_back(std::move(cfg));
//   }
//   produces<std::vector<lcvn::Result>>();
//   if (fEvaluators[0].fSliceLabel != "") {
//     produces<art::Assns<recob::Slice, lcvn::Result>>();
//   } else {
//     produces<art::Assns<lcvn::PixelMap, lcvn::Result>>();
//   }
// }
ICARUSCVNEvaluator::ICARUSCVNEvaluator(fhicl::ParameterSet const &pset)
    : EDProducer{pset},
      fModule(pset.get<std::string>("Module")),
      fPixelMapInput(pset.get<std::string>("PixelMapInput")),
      fPixelMapModuleLabel(pset.get<std::string>("PixelMapModuleLabel")),
      fSliceLabel(pset.get<art::InputTag>("SliceLabel", "")),
      fPFParticleModuleLabel(pset.get<art::InputTag>("PFParticleModuleLabel", "")),
      fT0Label(pset.get<art::InputTag>("T0Label", "")),
      fCVNType(pset.get<std::string>("CVNType")),
      fVerbose(pset.get<bool>("verbose", false)),
      fTFHandler{art::make_tool<lcvn::ICARUSITFNetHandler>(pset.get<fhicl::ParameterSet>("ICARUSTFHandler"))},
      fPMProducer(pset.get<fhicl::ParameterSet>("PixelMapProducer"))
  {
  produces<std::vector<lcvn::Result>>();
  if (fSliceLabel != "") {
    produces<art::Assns<recob::Slice, lcvn::Result>>();
  } else {
    produces<art::Assns<lcvn::PixelMap, lcvn::Result>>();
  }
}
//......................................................................
void ICARUSCVNEvaluator::produce(art::Event &evt) {

  /// Define containers for the things we're going to produce
  std::unique_ptr<std::vector<Result>> resultCol(new std::vector<Result>);
  auto assn1 = std::make_unique<art::Assns<recob::Slice, lcvn::Result>>();
  auto assn2 = std::make_unique<art::Assns<lcvn::PixelMap, lcvn::Result>>();
  bool isSlice = false;
  // for (auto &cfg : fEvaluators) {
    if (fCVNType == "TF" || fCVNType == "Tensorflow" || fCVNType == "TensorFlow") {
      if (fSliceLabel != "") { // by default use slice information
        if (fVerbose) std::cout << "Using slice information" << std::endl;
        isSlice = true;
        std::vector<art::Ptr<recob::Slice>> slcList;
        auto slcHandle = evt.getHandle<std::vector<recob::Slice>>(fSliceLabel);
        if (!slcHandle.isValid()) {
          throw cet::exception("ICARUSCVNEvaluator")
              << "Unable to get slices using label " << fSliceLabel;
        } else {
          if (fVerbose) std::cout << "Filling slice" << std::endl;
          art::fill_ptr_vector(slcList, slcHandle);
        }
        art::Handle<std::vector<recob::PFParticle>> PFPListHandle;
        std::vector<art::Ptr<recob::PFParticle>> PFPList;
        if (evt.getByLabel(fPFParticleModuleLabel, PFPListHandle)) {
          if (fVerbose) std::cout << "Filling pfp" << std::endl;
          art::fill_ptr_vector(PFPList, PFPListHandle);
        }

        art::FindManyP<recob::Hit> findManyHits(slcHandle, evt,
                                                fSliceLabel);
        art::FindManyP<recob::PFParticle> findManyPFPs(
            slcHandle, evt, fPFParticleModuleLabel);
        art::FindManyP<anab::T0> findManyT0s(PFPListHandle, evt, fT0Label);
        auto const detProp =
            art::ServiceHandle<detinfo::DetectorPropertiesService const>()
                ->DataFor(evt);
        for (auto const &slice : slcList) {
          if (fVerbose)
            std::cout << "********* " << fModule << " ********* "
                      << evt.run() << "  " << evt.subRun() << "  "
                      << evt.id().event() << "  " << slice->ID()
                      << "  **************\n";
          std::vector<float> pfp_T0_vec;
          if (findManyPFPs.isValid()) {
            std::vector<art::Ptr<recob::PFParticle>> slicePFPs = findManyPFPs.at(slice.key());
            if (slicePFPs.size()) {
              for (auto const &pfp : slicePFPs) {
                if (findManyT0s.isValid()) {
                  std::vector<art::Ptr<anab::T0>> T0_vec =
                      findManyT0s.at(pfp.key());
                  if (T0_vec.size()) {
                    for (auto const &T0 : T0_vec) {
                      pfp_T0_vec.push_back(T0->Time());
                    }
                  }
                }
              }
            }
          }

          float min_T0 = 0.;
          if (pfp_T0_vec.size()) {
            min_T0 = *min_element(pfp_T0_vec.begin(), pfp_T0_vec.end());
            if (fVerbose) std::cout << "min t0: " << min_T0 << std::endl;
          }

          if (findManyHits.isValid()) {
            std::vector<art::Ptr<recob::Hit>> slicehits = findManyHits.at(slice.key());
            if (fVerbose) std::cout << "slice key: " << slice.key() << " slice id: " << slice.id() << std::endl;
            fPMProducer.Set_fT0_value(min_T0);
            PixelMap pm = fPMProducer.ICARUSCreateMap(detProp, slicehits);
            auto nhits = fPMProducer.TotHits();
            pm.SetTotHits(nhits);
            // pm.fSliceID = slice->ID();
            std::string moduleSlcId = std::string(fModule) + " - slice id: " + std::to_string(slice.id()) + " - slice key: " + std::to_string(slice.key() + 1);
            std::vector<std::vector<float>> output = fTFHandler->Predict(pm, evt.id().event(), moduleSlcId);
            resultCol->emplace_back(output);
            util::CreateAssn(*this, evt, *resultCol, slice, *assn1);
	  } else {
            if (fVerbose) std::cout << "findManyHits: " << findManyHits.isValid() << std::endl;
          }
        }
      } else { // Try to read pixel maps saved in the file
        /// Load in the pixel maps
        load_and_predict(evt.id().event(), fTFHandler, fPixelMapInput);
      }
    } else {
      mf::LogError("ICARUSCVNEvaluator::produce")
          << "CVN Type not in the allowed list: Tensorflow" << std::endl;
      mf::LogError("ICARUSCVNEvaluator::produce")
          << "Exiting without processing events" << std::endl;
      return;
    }
  // }
  for (const auto &r : *resultCol) {
	for (const auto &row : r.fOutput) {
        	for (const auto &val : row) {
                  std::cout << val << " ";
                }
                std::cout << std::endl;
        }
  }

  evt.put(std::move(resultCol));
  if (isSlice) {
    evt.put(std::move(assn1));
  } else {
    evt.put(std::move(assn2));
  }
}

// --- Main prediction logic ---
void ICARUSCVNEvaluator::load_and_predict(
    const int event, std::unique_ptr<lcvn::ICARUSITFNetHandler> &handler,
    const std::string &pixelMapInput) {
  auto files = resolve_input_files(pixelMapInput);

  int fileCount = 0;
  for (const auto &f : files) {
    ++fileCount;
    std::cout << "Processing [" << fileCount << "/" << files.size() << "]"
              << "\n";

    try {
      size_t lastSlash = f.find_last_of('/');
      std::string fileName =
          (lastSlash != std::string::npos) ? f.substr(lastSlash + 1) : f;

      if (fileName.size() > 3 &&
          fileName.substr(fileName.size() - 3) == ".gz") {
        fileName = fileName.substr(0, fileName.size() - 3); // Remove ".gz"
      }
      std::vector<unsigned char> pixel_array = decompress_pixelmap_gz(f);
      std::vector<std::vector<float>> output =
          handler->PredictFromArray(pixel_array, fileName);
    } catch (const std::exception &e) {
      std::cerr << "Failed on " << f << ": " << e.what() << "\n";
    }
  }
}

std::vector<unsigned char>
ICARUSCVNEvaluator::decompress_pixelmap_gz(const std::string &filename) {
  constexpr size_t kPixelArraySize =
      3 * 500 * 500; // Set this to your actual image shape
  std::vector<unsigned char> pixel_array(kPixelArraySize);

  std::ifstream file(filename, std::ios::binary | std::ios::ate);
  if (!file)
    throw std::runtime_error("Cannot open file: " + filename);
  auto size = file.tellg();
  file.seekg(0, std::ios::beg);

  std::vector<char> compressed(size);
  file.read(compressed.data(), size);
  file.close();

  uLongf dest_len = kPixelArraySize;
  int res = uncompress(pixel_array.data(), &dest_len,
                       reinterpret_cast<Bytef *>(compressed.data()), size);

  if (res != Z_OK || dest_len != kPixelArraySize) {
    throw std::runtime_error(
        "Decompression failed or pixel array size mismatch in " + filename);
  }

  return pixel_array;
}

// --- Load all .gz files from a folder ---
std::vector<std::string>
ICARUSCVNEvaluator::load_from_folder(const std::string &folder_path) {
  std::vector<std::string> files;
  for (const auto &entry :
       std::filesystem::recursive_directory_iterator(folder_path)) {
    if (entry.is_regular_file() && entry.path().extension() == ".gz") {
      files.push_back(entry.path().string());
    }
  }
  if (files.empty())
    throw std::runtime_error("No .gz files found in folder: " + folder_path);
  return files;
}

// --- Load file list or single file ---
std::vector<std::string>
ICARUSCVNEvaluator::resolve_input_files(const std::string &pixelMapInput) {
  std::filesystem::path input(pixelMapInput);

  if (!std::filesystem::exists(input)) {
    throw std::runtime_error("Input path does not exist: " + pixelMapInput);
  }

  if (std::filesystem::is_regular_file(input)) {
    if (input.extension() == ".gz") {
      std::cout << "Single File" << std::endl;
      return {pixelMapInput}; // Single .gz file
    } else {
      // Text file containing list of .gz files
      std::ifstream listfile(pixelMapInput);
      if (!listfile)
        throw std::runtime_error("Cannot open list file: " + pixelMapInput);
      std::vector<std::string> paths;
      std::string line;
      while (std::getline(listfile, line)) {
        if (!line.empty())
          paths.push_back(line);
      }
      return paths;
    }
  } else if (std::filesystem::is_directory(input)) {
    return load_from_folder(pixelMapInput); // Folder
  }

  throw std::runtime_error("Unsupported input type: " + pixelMapInput);
}

DEFINE_ART_MODULE(lcvn::ICARUSCVNEvaluator)
} // namespace lcvn
////////////////////////////////////////////////////////////////////////
