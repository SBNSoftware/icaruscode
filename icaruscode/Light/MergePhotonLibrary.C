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
#include "TStopwatch.h"

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


/// Returns a string with the name of the specified compression algorithm.
std::string CompressionAlgorithmName(ROOT::ECompressionAlgorithm algo);


/**
 * @brief Collects data trees and puts together a photon library tree.
 * @param outputFile name of the file to create
 * @param inputFilePattern pattern to find the input files
 * @param version version of the software used to create the input trees
 * @param date string representing the date the input trees were created
 * @param inputTreeName (default: `"PhotonLibrary"`) input tree name
 * @param inputTreeDir (default: `"pmtresponse"`) directory path to input trees
 * @param compressionFactor (default: `8`) compression level (1 to 9, higher
 *                          takes longer and compresses more)
 * @param compressionAlgo (default: `ROOT::kLZMA`) compression algorithm
 *                        (see `ROOT::ECompressionAlgorithm`)
 * @return an error code, `0` on success
 * 
 * This script reads the trees with name `inputTreeName` from ROOT directory
 * `inputTreeDir` from all the ROOT files matching the `inputFilePattern`,
 * and clones it into a single tree. The tree is stored in a new ROOT file,
 * `outputFile`, in a ROOT directory also named as `inputTreeDir`.
 * The specified compression level is applied to all the branches of the tree
 * (the cloning is _not_ fast).
 *
 * The script works as a glorified `hadd`, and *will not combine multiple data
 * from the same voxel and channel*.
 * 
 * The output file must not exist yet.
 * 
 * @note The support of wildcards is the same as in `TChain::Add()`.
 * 
 * Compression level
 * ------------------
 * 
 * From ROOT 6.12 documentation, it appears that the LZ4 compression algorithm
 * (`ROOT::kLZ4`) is the most effective, but the _decompression_ time also
 * depends on the compression level, while ZLib and LZMA do only very mildly.
 * The recommendation is to use compression level `8` with LZMA (`ROOT::kLZMA`),
 * but only `4` with LZ4, and `1` with ZLib (`ROOT::ZLIB`). My observations
 * in the case of a photon library tree with 250M entries (1 GB of data for each
 * of three branches) is in the following table.
 * 
 * ------  -------  ------  -------  ---------  ------------ ---------  ------
 *  algo    level    time    voxel    channel    visibility   total      read 
 * ------  -------  ------  -------  ---------  ------------ ---------  ------
 *  ZLib    `4`        5'    20 MB     110 MB        370 MB    500 MB     34"
 *  ZLib    `8`       16'    10 MB      60 MB        280 MB    350 MB     37"
 *  ZLib    `9`       20'    10 MB      55 MB        280 MB    345 MB     34"
 *  LZ4     `4`        8'    20 MB      80 MB        370 MB    470 MB     40" 
 *  LZ4     `8`       14'    20 MB      65 MB        335 MB    420 MB     41"
 *  LZMA    `5`       22'     5 MB      45 MB        200 MB    250 MB     49"
 *  LZMA    `7`       30'     5 MB      45 MB        195 MB    245 MB     48"
 * ------  -------  ------  -------  ---------  ------------ ---------  ------
 * 
 * For each compression setting, the "time" (wall time) to run the whole merging
 * process is reported (in minutes). The overhead for reading the data from the
 * network with XRootD is about 2-3 minutes, and is not dependent on the
 * settings (it does depend on contingent state as caching, though). The "read"
 * time is the CPU time reported by `ReadPhotonLibrary()` macro, so the actual
 * time spent in reading the library will present an overhead, especially with
 * remote data source, proportional to the compressed size.
 *
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
  std::string const& inputTreeDir = "pmtresponse",
  int compressionLevel = 8,
  ROOT::ECompressionAlgorithm compressionAlgo = ROOT::kLZMA
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
  
  // set the compression level to all branches
  for (TObject* pObj: *(pDestTree->GetListOfBranches())) {
    TBranch* pBranch = static_cast<TBranch*>(pObj);
    pBranch->SetCompressionAlgorithm(compressionAlgo);
    if (compressionLevel >= 0) pBranch->SetCompressionLevel(compressionLevel);
  }
  
  // now do copy the data
  std::cout << "Copying tree data..." << std::endl;
  TStopwatch timer;
  pDestTree->CopyEntries(pSourceTree);
  timer.Stop();
  
  std::cout << pDestTree->GetEntries() << " entries transferred in " << timer.RealTime()
    << " seconds; compression factors:" << std::endl;
  for (TObject* pObj: *(pDestTree->GetListOfBranches())) {
    TBranch* pBranch = static_cast<TBranch*>(pObj);
    auto const branchSize = pBranch->GetTotBytes();
    auto const compressedSize = pBranch->GetZipBytes();
    std::cout << " * branch '" << pBranch->GetName() << "': "
      << branchSize << " ==> " << compressedSize << " bytes";
    if (branchSize > 0) {
      std::cout << " (compression factor: "
        << (static_cast<float>(compressedSize)/branchSize)
        << ", " << CompressionAlgorithmName(ROOT::ECompressionAlgorithm(pBranch->GetCompressionAlgorithm()))
        << " level " << pBranch->GetCompressionLevel()
        << ")";
    }
    std::cout << std::endl;
  }
  
  pDestTree->SetDirectory(pDestDir);
  
  // make sure that the tree is written...
  pDestTree->Write();
  
  //
  // close and go
  //
  pFDest->Close();
  
  return 0;
  
} // MergePhotonLibrary()


// -----------------------------------------------------------------------------
std::string CompressionAlgorithmName(ROOT::ECompressionAlgorithm algo) {
  switch (algo) {
    case ROOT::kUseGlobalCompressionSetting: return "global setting from rootrc";
    case ROOT::kZLIB: return "ZLib";
    case ROOT::kLZMA: return "LZMA";
    case ROOT::kOldCompressionAlgo: return "ZLib (legacy)";
    case ROOT::kLZ4: return "LZ4";
    case ROOT::kUndefinedCompressionAlgorithm: return "undefined";
    default: return "unknown?!?";
  } // switch
} // CompressionAlgorithmName()


// -----------------------------------------------------------------------------

