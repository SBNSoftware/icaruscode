////////////////////////////////////////////////////////////////////////////////////
/// \file    ICARUSTFNetHandler.cxx
/// \brief   ICARUSTFNetHandler for CVN
/// \author  Varuna Meddage (copied from
/// larrecodnn/CVN/tools/TFNetHandler_tool.cc) modified by Felipe Wieler for
/// Icarus CVN
///////////////////////////////////////////////////////////////////////////////////

#include "cetlib/getenv.h"
#include <fstream>
#include <iostream>
#include <cstdlib>   // for std::getenv
#include <string>
#include <regex>

#include "art/Utilities/ToolMacros.h"
#include "canvas/Utilities/Exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "icaruscode/ICARUSCVN/module_helpers/ICARUSITFNetHandler.h"
#include "icaruscode/ICARUSCVN/tf/tf_graph.h"
#include "larrecodnn/CVN/func/CVNImageUtils.h"

namespace lcvn {

/// Wrapper for caffe::Net which handles construction and prediction
class ICARUSTFNetHandler : public ICARUSITFNetHandler {
public:
  /// Constructor which takes a pset with DeployProto and ModelFile fields
  explicit ICARUSTFNetHandler(const fhicl::ParameterSet &pset);

  /// Return prediction arrays for PixelMap
  std::vector<std::vector<float>> Predict(const PixelMap &pm, const int event, const std::string cryo) const override;
  std::vector<std::vector<float>> PredictFromArray(const std::vector<unsigned char> &pa, const std::string event) const override;

private:
  std::string fLibPath; ///< Library path (typically dune_pardata...)
  std::string
      fTFProtoBuf; ///< location of the tf .pb file in the above path or the
                   ///< directory containing model files in SavedModel format
                   ///< (set UseBundle = true in this case)
  bool fUseLogChargeScale;  ///< Is the charge using a log scale?
  unsigned int fImageWires; ///< Number of wires for the network to classify
  unsigned int fImageTDCs;  ///< Number of tdcs for the network to classify
  std::vector<bool> fReverseViews; ///< Do we need to reverse any views?
  bool fUseBundle; ///< Use a bundled model saved in the SavedModel format from
                   ///< Tensorflow
  std::vector<std::string> fInputs;
  std::vector<std::string> fOutputs;
  int fNInputs;
  int fNOutputs;
  std::unique_ptr<tf::Graph> fTFGraph; ///< Tensorflow graph
  bool fverbose;
};

std::string expandEnvVars(const std::string& text) {
    static std::regex env("\\$\\{([^}]+)\\}");
    std::smatch match;
    std::string result = text;
    while (std::regex_search(result, match, env)) {
        const char* val = std::getenv(match[1].str().c_str());
        result.replace(match.position(0), match.length(0),
                       val ? val : ("${" + match[1].str() + "}"));
    }
    std::cout << result << std::endl;
    return result;
}

ICARUSTFNetHandler::ICARUSTFNetHandler(const fhicl::ParameterSet &pset)
    // : fLibPath(cet::getenv(pset.get<std::string>("LibPath", "")))
    : fTFProtoBuf(expandEnvVars(pset.get<std::string>("TFProtoBuf"))),
      fUseLogChargeScale(pset.get<bool>("ChargeLogScale")),
      fImageWires(pset.get<unsigned int>("NImageWires")),
      fImageTDCs(pset.get<unsigned int>("NImageTDCs")),
      fReverseViews(pset.get<std::vector<bool>>("ReverseViews")),
      fUseBundle(pset.get<bool>("UseBundle")),
      fInputs(pset.get<std::vector<std::string>>("Inputs")),
      fOutputs(pset.get<std::vector<std::string>>("Outputs")),
      fNInputs(pset.get<int>("NInputs")), fNOutputs(pset.get<int>("NOutputs")),
      fverbose(pset.get<bool>("verbose")) {

  // Construct the TF Graph object. The empty vector {} is used since the
  // protobuf file gives the names of the output layer nodes
  mf::LogInfo("TFNetHandler")
      << "Loading network: " << fTFProtoBuf << std::endl;
  fTFGraph = tf::Graph::create(fTFProtoBuf.c_str(), fInputs, fOutputs,
                               fUseBundle, fNInputs, fNOutputs);
  if (!fTFGraph) {
    art::Exception(art::errors::Unknown)
        << "Tensorflow model not found or incorrect";
  }
}

// Check the network outputs
bool check(const std::vector<std::vector<float>> &outputs) {
  if (outputs.size() == 1)
    return true;
  size_t aux = 0;
  for (size_t o = 0; o < outputs.size(); ++o) {
    size_t aux2 = 0;

    for (size_t i = 0; i < outputs[o].size(); ++i)
      if (outputs[o][i] == 0.0 || outputs[o][i] == 1.0)
        aux2++;
    if (aux2 == outputs[o].size())
      aux++;
  }
  return aux == outputs.size() ? false : true;
}

// Fill outputs with value -3
void fillEmpty(std::vector<std::vector<float>> &outputs) {
  for (auto &output : outputs) {
    output.assign(output.size(), -3.0);
  }

  return;
}

std::vector<std::vector<float>>
ICARUSTFNetHandler::Predict(const PixelMap &pm, const int event, const std::string cryo) const {

  CVNImageUtils imageUtils(fImageWires, fImageTDCs, 3);
  // Configure the image utility
  imageUtils.SetViewReversal(fReverseViews);
  imageUtils.SetImageSize(fImageWires, fImageTDCs, 3);
  imageUtils.SetLogScale(fUseLogChargeScale);
  imageUtils.SetPixelMapSize(pm.NWire(), pm.NTdc());

  ImageVectorF thisImage;
  imageUtils.ConvertPixelMapToImageVectorF(pm, thisImage);
  std::vector<ImageVectorF> vecForTF;

  vecForTF.push_back(thisImage);
  std::vector<std::vector<std::string>> target_names;
  if (fNOutputs == 1){
    target_names = {
        {"CC Numu", "CC Nue", "Cosmic", "NC"}};
  }
  else if (fNOutputs == 3){
    target_names = {
        {"neutrino", "antineutrino"},
        {"CC Numu", "CC Nue", "NC"},
        {"CC QE", "CC Res", "CC DIS", "CC Other", "NC"}};
  }
  else if(fNOutputs == 7){
    target_names = {
        {"neutrino", "antineutrino", "NULL"},
        {"CC Numu", "CC Nue", "CC Nutau", "NC"},
        {"CC QE", "CC Res", "CC DIS", "CC Other", "NULL"},
        {"0", "1", "2", ">2"},
        {"0", "1", "2", ">2"},
        {"0", "1", "2", ">2"},
        {"0", "1", "2", ">2"}};
  }

  std::vector<std::vector<std::vector<float>>>
      cvnResults; // shape(samples, #outputs, output_size)
  bool status = false;

  int counter = 0;

  while (status == false) { // do until it gets a correct result
    cvnResults = fTFGraph->run(vecForTF);
    status = check(cvnResults[0]);

    counter++;
    if (counter == 10) {
      std::cout << "Error, CVN never outputing a correct result. Filling "
                   "result with zeros.";
      std::cout << std::endl;
      fillEmpty(cvnResults[0]);
      break;
    }
  };

  //TO MAINTAIN THE SBND FORMAT FOR THE CAFMAKER FILES.
  for (auto &outer : cvnResults) {
    for (auto &inner : outer) {
        inner.insert(inner.begin() + 2, 0.0f);
    }
  }

  if (fverbose) {
    
    std::string s = "Classifier summary: e" + std::to_string(event);
    int output_index = 0;
    for (auto const &output : cvnResults[0]) {
      s += "\nOutput " + std::to_string(output_index) + ": ";
      for (auto const v : output)
        s += std::to_string(v) + ", ";
      s += cryo;
      std::cout << s << std::endl;
      output_index++;
    }

    std::ofstream outfile("classifier_summary.txt", std::ios::app); // Save to file

    outfile << s;
    outfile << std::endl;

    outfile.close();
  }

  return cvnResults[0];
}

std::vector<std::vector<float>>
ICARUSTFNetHandler::PredictFromArray(const std::vector<unsigned char> &pa, const std::string event) const {

  CVNImageUtils imageUtils(fImageWires, fImageTDCs, 3);
  // Configure the image utility
  imageUtils.SetViewReversal(fReverseViews);
  imageUtils.SetImageSize(fImageWires, fImageTDCs, 3);
  imageUtils.SetLogScale(fUseLogChargeScale);

  ImageVectorF thisImage;
  imageUtils.ConvertPixelArrayToImageVectorF(pa, thisImage);
  std::vector<ImageVectorF> vecForTF;

  vecForTF.push_back(thisImage);
  std::vector<std::vector<std::string>> target_names;
  if (fNOutputs == 1){
    target_names = {
        {"CC Numu", "CC Nue", "NC"}};
  }
  else if (fNOutputs == 3){
    target_names = {
        {"neutrino", "antineutrino"},
        {"CC Numu", "CC Nue", "NC"},
        {"CC QE", "CC Res", "CC DIS", "CC Other", "NC"}};
  }
  else if(fNOutputs == 7){
    target_names = {
        {"neutrino", "antineutrino", "NULL"},
        {"CC Numu", "CC Nue", "CC Nutau", "NC"},
        {"CC QE", "CC Res", "CC DIS", "CC Other", "NULL"},
        {"0", "1", "2", ">2"},
        {"0", "1", "2", ">2"},
        {"0", "1", "2", ">2"},
        {"0", "1", "2", ">2"}};
  }

  std::vector<std::vector<std::vector<float>>>
      cvnResults; // shape(samples, #outputs, output_size)
  bool status = false;

  int counter = 0;

  while (status == false) { // do until it gets a correct result
    cvnResults = fTFGraph->run(vecForTF);
    status = check(cvnResults[0]);

    counter++;
    if (counter == 10) {
      std::cout << "Error, CVN never outputing a correct result. Filling "
                   "result with zeros.";
      std::cout << std::endl;
      fillEmpty(cvnResults[0]);
      break;
    }
  };

  if (fverbose) {
    
    std::string s = "Classifier summary: " + event;
    int output_index = 0;
    for (auto const &output : cvnResults[0]) {
      s += "\nOutput " + std::to_string(output_index) + ": ";
      auto max_it = std::max_element(output.begin(), output.end());
      int argmax = std::distance(output.begin(), max_it);
      if (output_index == 0)
        argmax = std::round(output[0]);
      for (auto const v : output)
        s += std::to_string(v) + ", ";
      s += " " + target_names[output_index][argmax]  + " - " + std::to_string(*max_it);
      std::cout << event << target_names[output_index][argmax] << " - " << *max_it
                << std::endl;
      output_index++;
    }

    std::ofstream outfile("classifier_summary.txt", std::ios::app); // Save to file

    outfile << s;
    outfile << std::endl;

    outfile.close();
  }

  return cvnResults[0];
}

} // namespace lcvn
DEFINE_ART_CLASS_TOOL(lcvn::ICARUSTFNetHandler)
