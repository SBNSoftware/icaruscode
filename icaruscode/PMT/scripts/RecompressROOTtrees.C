/**
 * @file    RecompressROOTtrees.C
 * @brief   ROOT macro creating a new ROOT file with differently compressed
 *          trees.
 * @author  Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date    October 10, 2020
 * 
 * Changes
 * --------
 * 
 * * version 1:
 *     original version
 * 
 * Examples
 * =========
 * 
 * Real life example:
 *     
 *     root -l -q 'RecompressROOTtrees.C(
 *       "PhotonLibraryTest.root",
 *       "PhotonLibraryTestNew.root",
 *       8, ROOT::kLZMA
 *       )'
 * 
 */

// ROOT libraries
#include "Compression.h" // ROOT::kLZ4
#include "TFile.h"
#include "TDirectory.h"
#include "TObjArray.h"
#include "TTree.h"
#include "TNamed.h"
// #include "TStopwatch.h"

// C/C++ libraries
#include <iostream>
#include <string>
#include <set>
#include <utility> // std::pair
#include <cassert>


// -----------------------------------------------------------------------------
static std::string const ScriptName { "RecompressROOTtrees.C" };
static unsigned int const ScriptVersion { 1 };


// -----------------------------------------------------------------------------
// --- forward declarations
// -----------------------------------------------------------------------------

/// Helper not to perturb the current ROOT directory.
class TgDirectoryGuard {
  TDirectory* pOldDir = nullptr;
    public:
  TgDirectoryGuard(): pOldDir(gDirectory) {}
  TgDirectoryGuard(TDirectory& newDir): TgDirectoryGuard() { newDir.cd(); }
  ~TgDirectoryGuard()
    {
      if (pOldDir) pOldDir->cd();
      else { gDirectory = nullptr; gFile = nullptr; }
    }
}; // classTgDirectoryGuard


/// Returns a string with the name of the specified compression algorithm.
std::string CompressionAlgorithmName(ROOT::ECompressionAlgorithm algo);

// -----------------------------------------------------------------------------
class FileCopy {
  
  struct Context_t {
    int level = 0;
  };
  
  int fCompressionLevel;
  ROOT::ECompressionAlgorithm fCompressionAlgo;
  
    public:
  
  FileCopy(
    int compressionLevel,
    ROOT::ECompressionAlgorithm compressionAlgo
    )
    : fCompressionLevel(compressionLevel)
    , fCompressionAlgo(compressionAlgo)
    {}
  
  
  /**
   * @brief Copies the whole content of the specified ROOT file into another.
   * @param destPath path of the file to be created (must not exist)
   * @param sourcePath path of the file to be copied
   * @return 0 on success, an error code otherwise
   * 
   * 
   * Error codes
   * ------------
   * 
   * * `2`: could not open the source file
   * * `3`: could not create the destination file
   * 
   */
  int copyFile
    (std::string const& destPath, std::string const& sourcePath) const;
  
  int copyDirectoryContent
    (TDirectory& dest, TDirectory const& source, Context_t context) const;
  
  int recompressTree
    (TDirectory& destDir, TTree& srcTree, Context_t context) const;
  
}; // class FileCopy


// -----------------------------------------------------------------------------

int RecompressROOTtrees(
  std::string const& sourcePath, std::string const& destPath,
  int compressionLevel,
  ROOT::ECompressionAlgorithm compressionAlgo
  )
{
  
  FileCopy copier(compressionLevel, compressionAlgo);
  
  int res = copier.copyFile(destPath, sourcePath);
  
  return res;
} // RecompressROOTtrees()


// -----------------------------------------------------------------------------
// ---  implementations follow
// -----------------------------------------------------------------------------
int FileCopy::copyFile
  (std::string const& destPath, std::string const& sourcePath) const
{
  
  TgDirectoryGuard dg;
  
  TFile srcFile { sourcePath.c_str(), "READ" };
  if (!srcFile.IsOpen()) {
    std::cerr << "Could not read file '" << sourcePath << "'!" << std::endl;
    return 2;
  }
  
  int const compression
    = ROOT::CompressionSettings(fCompressionAlgo, fCompressionLevel);
  
  TFile destFile { destPath.c_str(), "CREATE", "", compression };
  if (!destFile.IsOpen()) {
    std::cerr << "Could not create file '" << destPath
      << "'! (does it exist already?)" << std::endl;
    return 3;
  }
  
  std::cout << "'" << srcFile.GetPath() << "' => '" << destFile.GetPath() << "'"
    << std::endl;
  return copyDirectoryContent(destFile, srcFile, {});
  
} // FileCopy::copyFile()


int FileCopy::copyDirectoryContent
  (TDirectory& dest, TDirectory const& source, Context_t context) const
{
  
  TgDirectoryGuard dg(dest);
  std::string const indentStr = std::string(context.level, ' ');
  
  std::cout << indentStr << "'" << source.GetName() << "' => '"
    << dest.GetName() << "'" << std::endl;
  
  Context_t subContext = context;
  ++subContext.level;
  
  // it may be useless precaution, but I'll try not to load the same tree twice
  std::set<std::string> trees;
  unsigned int nSkipped = 0U;
  unsigned int nTreeErrors = 0U, nDirErrors = 0U;
  for (TObject* keyObj: *(source.GetListOfKeys())) {
    
    TKey* key = static_cast<TKey*>(keyObj);
    
    TClass const* objClass = TClass::GetClass(key->GetClassName());
    if (!objClass || !(objClass->InheritsFrom(TObject::Class()))) {
      std::cerr << indentStr << " ERROR: object '" << key->GetName()
        << "' of unknown class '" << key->GetClassName() << "' skipped!"
        << std::endl;
      ++nSkipped;
      continue;
    }
    
    std::cerr << indentStr << " '" << key->GetName()
      << "' ('" << key->GetClassName() << "')" << std::endl;
    
    if (objClass->InheritsFrom(TTree::Class())) {
      if (trees.insert(std::string{ key->GetName() }).second) {
        TTree* tree = dynamic_cast<TTree*>(key->ReadObj());
        assert(tree);
        int res = recompressTree(dest, *tree, subContext);
        if (res != 0) ++nTreeErrors;
        delete tree;
      }
    }
    else if (objClass->InheritsFrom(TDirectory::Class())) {
      TgDirectoryGuard dg; // just in case
      TDirectory* srcDir = dynamic_cast<TDirectory*>(key->ReadObj());
      assert(srcDir);
      TDirectory* destDir = dest.mkdir(srcDir->GetName());
      int res = copyDirectoryContent(*destDir, *srcDir, subContext);
      if (res != 0) ++nDirErrors;
      destDir->Write();
      delete destDir;
    }
    else {
      TObject* obj = key->ReadObj();
      dest.WriteTObject(obj);
      delete obj;
    }
    
  } // for
  
  return ((nSkipped + nTreeErrors + nDirErrors) == 0)? 0: 1;
  
} // FileCopy::copyDirectory()


int FileCopy::recompressTree
  (TDirectory& destDir, TTree& srcTree, Context_t context) const
{
  
  std::string const indentStr = std::string(context.level, ' ');
  
  std::cout << indentStr << " compression "
    << CompressionAlgorithmName(fCompressionAlgo)
    << " level " << fCompressionLevel
    << std::endl;
  
  TgDirectoryGuard dg(destDir);
  
  // the compression settings are the same as in the destination file
  TTree* destTree = srcTree.CloneTree();
  if (!destTree) return 1;
  std::cout << indentStr << "  (cloned)" << std::flush;
  
  destTree->SetDirectory(&destDir); // we like to be very very safe
  destTree->Write();
  std::cout << "  (written)" << std::flush;
  
  delete destTree;
  std::cout << "  (deleted)" << std::endl;
  
  return 0;
} // FileCopy::recompressTree()


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

