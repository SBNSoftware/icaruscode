/**
 * @file    ReadPhotonLibrary.C
 * @brief   ROOT macro for `ReadPhotonLibrary()` function.
 * @author  Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date    August 2, 2018
 * @version 1.0
 * 
 * Example:
 *     
 *     root -l -q 'ReadPhotonLibrary.C+("PhotonLibrary.root")'
 *
 */

// ROOT libraries
#include "Compression.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TBranch.h"
#include "TNamed.h"
#include "TStopwatch.h"

// C/C++ libraries
#include <iostream>
#include <string>


class TgDirectoryGuard {
  TDirectory* pOldDir = nullptr;
    public:
  TgDirectoryGuard(): pOldDir(gDirectory) {}
  ~TgDirectoryGuard() { if (pOldDir) pOldDir->cd(); gDirectory = pOldDir; }
}; // classTgDirectoryGuard


std::string CompressionAlgorithmName(ROOT::ECompressionAlgorithm algo);


/**
 * @brief Fully read a photon library tree.
 * @param treeFile ROOT file containing the tree
 * @param treePath (default: `pmtresponse/PhotonLibraryData`) full path to the ROOT tree
 * @return an error code, `0` on success
 * 
 * This script reads all the entries from the specified tree.
 * 
 * Return codes
 * -------------
 * 
 * * `0`: success
 * * `1`: no ROOT object available with path `treeDir/treePath`
 * * `2`: failed to create destination ROOT file
 * * `3`: no `treeDir` ROOT directory found
 * * `4`: the object in `treePath` is not a `TTree`
 * 
 */
int ReadPhotonLibrary(
  std::string const& treeFile,
  std::string const& treeName = "PhotonLibraryData",
  std::string const& treeDir = "pmtresponse"
) {
  
  // protect the global state that ROOT so much loves
  TgDirectoryGuard gDirectoryGuard;
  
  //
  // open source
  //
  TFile FSource(treeFile.c_str(), "READ");
  if (!FSource.IsOpen()) return 2; // assume ROOT was already noisy about this
  
  TDirectory* pDir
    = treeDir.empty()? (&FSource): FSource.GetDirectory(treeDir.c_str());
  if (!pDir) {
    std::cerr
      << "Can't find ROOT directory '" << treeDir << "' in '" << FSource.GetPath() << "'."
      << std::endl;
    return 3;
  }
  
  // read the tree
  TObject* pObj = nullptr;
  pObj = pDir->Get(treeName.c_str());
  if (!pObj) {
    std::cerr << "No object '" << treeName << "' available in '" << pDir->GetPath() << "'"
      << std::endl;
    return 1;
  }

  TTree* pTree = dynamic_cast<TTree*>(pObj);
  if (!pTree) {
    std::cerr << "Object '" << pObj->GetName() << " from '" << pDir->GetPath() << "' is a "
      << pObj->IsA()->GetName() << ", not compatible with TTree." << std::endl;
    return 4;
  }
 
  //
  // write some metadata
  //
  std::cout << "Input tree: " << pTree->GetTitle() << std::endl;
  pObj = pDir->Get("Version");
  std::cout << "Version:   '" << (pObj? pObj->GetTitle(): "unknown") << "'" << std::endl;
  pObj = pDir->Get("Date");
  std::cout << "Date:      '" << (pObj? pObj->GetTitle(): "unknown") << "'" << std::endl;
  pObj = pDir->Get("Source");
  std::cout << "Source:    '" << (pObj? pObj->GetTitle(): "unknown") << "'" << std::endl;
  
  //
  // read the tree
  //
  std::cout << "Reading all entries in the tree..." <<std::endl;
  TStopwatch timer;
  unsigned int iEntry = 0U;
  while (pTree->GetEntry(iEntry++, 1) > 0);
  timer.Stop();
  std::cout << iEntry << " entries read in " << timer.RealTime() << " seconds ("
    << timer.CpuTime() << "\" CPU time)." << std::endl;
  
  std::cout << "Tree has " << pTree->GetListOfBranches()->GetEntries()
    << " branches:" << std::endl;
  std::size_t totalSize = 0U, totalCompressedSize = 0U;
  for (TObject* pObj: *(pTree->GetListOfBranches())) {
    TBranch* pBranch = static_cast<TBranch*>(pObj);
    auto const branchSize = pBranch->GetTotBytes("*");
    auto const compressedSize = pBranch->GetZipBytes("*");
    totalSize += branchSize;
    totalCompressedSize += compressedSize;
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
  } // for branches 
  std::cout << " * total: "
    << totalSize << " ==> " << totalCompressedSize << " bytes";
  if (totalSize > 0) {
    std::cout << " (compression factor: "
      << (static_cast<float>(totalCompressedSize)/totalSize) << ")";
  }
  std::cout << std::endl;
  
  
  // close everything, everything gets deleted, goodbye
  return 0;
  
} // ReadPhotonLibrary()



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



