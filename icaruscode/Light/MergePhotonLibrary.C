/**
 * @file    MergePhotonLibrary.C
 * @brief   ROOT macro for `MergePhotonLibrary()` function.
 * @author  Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date    August 1, 2018
 * @version 1.0
 */
/// 
/// Real life example:
///     
///     root -l -q 'MergePhotonLibrary.C++(
///       "PhotonLibraryTest.root",
///       "root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/persistent/users/cfarnese/pl_01082018/test_hist_icarus_00*.root",
///       "v06_85_00",
///       "20180801",
///       "PhotonLibraryData",
///       "pmtresponse"
///       )'
///     
/// 
///

// ROOT libraries
#include "Compression.h" // ROOT::kLZ4
#include "TFile.h"
#include "TDirectory.h"
#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"
#include "TNamed.h"

// C/C++ libraries
#include <iostream>
#include <string>
#include <memory> // std::make_unique()


class TgDirectoryGuard {
  TDirectory* pOldDir = nullptr;
    public:
  TgDirectoryGuard(): pOldDir(gDirectory) {}
  ~TgDirectoryGuard() { if (pOldDir) pOldDir->cd(); gDirectory = pOldDir; }
}; // classTgDirectoryGuard


/**
 * @brief Collects data trees and puts together a photon library tree.
 * @param outputFile name of the file to create
 * @param inputFilePattern pattern to find the input files
 * @param version version of the software used to create the input trees
 * @param date string representing the date the input trees were created
 * @param inputTreeName (default: `"PhotonLibrary"`) input tree name
 * @param inputTreeDir (default: `"pmtresponse"`) directory path to input trees
 * @return an error code, `0` on success
 * 
 * 
 * The output file must not exist yet.
 * 
 * @note The support of wildcards is the same as in `TChain`.
 * 
 * Return codes
 * -------------
 * 
 * * `0`: success
 * * `1`: no file matched the input file pattern
 * * `2`: failed to create destination ROOT file
 * * `3`: failed to create destination directory in the destination ROOT file
 * 
 */
int MergePhotonLibrary(
  std::string const& outputFile,
  std::string const& inputFilePattern,
  std::string const& version,
  std::string const& date,
  std::string const& inputTreeName = "PhotonLibraryData",
  std::string const& inputTreeDir = "pmtresponse"
) {
  
  // protect the global state that ROOT so much loves
  TgDirectoryGuard gDirectoryGuard;
  
  //
  // parse the name
  //
  
  //
  // open source
  //
  TChain* pSourceTree = new TChain(inputTreeName.c_str(), "Photon library");
  
  std::string pattern = inputFilePattern;
  if (!inputTreeDir.empty()) pattern += '/' + inputTreeDir;
  if (!inputTreeName.empty()) pattern += '/' + inputTreeName;
  auto const nFiles = pSourceTree->Add(pattern.c_str());
  if (nFiles == 0) {
    std::cerr << "No file matched '" << pattern << "'" << std::endl;
    return 1;
  }
  std::cout << nFiles << " files matched." << std::endl;
  
  //
  // create destination
  //
  auto pFDest = std::make_unique<TFile>(
    outputFile.c_str(), "CREATE",
    ("Photon library version " + version + " (" + date + ")").c_str(),
    4 // compression level (ignored with fast tree cloning)
    );
  if (!pFDest->IsOpen()) return 2; // assume ROOT was already noisy about this
  pFDest->cd();
  
  // create the output directory
  TDirectory* pDestDir
    = inputTreeDir.empty()? pFDest.get(): pFDest->mkdir(inputTreeDir.c_str());
  if (!pDestDir) {
    std::cerr << "Creation of '" << inputTreeDir << "' in file '"
      << pDestDir->GetPath() << "' failed!" << std::endl;
    return 3;
  }
  pDestDir->cd();
  
  //
  // write tags in the destination
  //
  auto pVersion = new TNamed("Version", version.c_str());
  pVersion->Write();
  auto pDate = new TNamed("Date", date.c_str());
  pDate->Write();
  auto pSource = new TNamed("Source", pattern.c_str());
  pSource->Write();
  
  //
  // clone the tree
  //
  
  // prepare a new tree for cloning, but do not copy any data just jet
  std::cout << "Setting up the new tree..." << std::endl;
  TTree* pDestTree = pSourceTree->CloneTree(0);
  
  // ROOT 6.12 recommended compression: LZ4 level 4
  int const compressionLevel = ROOT::kLZ4 + 4;
  
  // set the compression level to all branches
  for (TObject* pBranch: *(pDestTree->GetListOfBranches()))
    static_cast<TBranch*>(pBranch)->SetCompressionSettings(compressionLevel);
  
  // now do copy the data
  std::cout << "Copying tree data..." << std::endl;
  pDestTree->CopyEntries(pSourceTree);
  
  std::cout << pDestTree->GetEntries() << " entries transferred." << std::endl;
  pDestTree->SetDirectory(pDestDir);
  
  // make sure that the tree is written...
  pDestTree->Write();
  
  //
  // close and go
  //
  pFDest->Close();
  
  return 0;
  
} // MergePhotonLibrary()
