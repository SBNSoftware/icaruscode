/**
 * @file    MergePhotonLibrary.C
 * @brief   ROOT macro for `MergePhotonLibrary()` function.
 * @author  Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date    August 1, 2018
 * 
 * Changes
 * --------
 * 
 * * version 5:
 *     * missing voxel check reports as errors only blocks larger than
 *       `MinMissingVoxelBlockReport`
 * 
 * * version 4:
 *     * support file lists
 * 
 * * version 3:
 *     * added check on the voxel list
 * 
 * * version 2:
 *     * added support for metadata
 * 
 * 
 * Examples
 * =========
 * 
 * Real life example:
 *     
 *     root -l -q 'MergePhotonLibrary.C++(
 *       "PhotonLibraryTest.root",
 *       "root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/persistent/users/cfarnese/pl_01082018/test_hist_icarus_00*.root",
 *       "v06_85_00",
 *       "20180801",
 *       "PhotonLibraryData",
 *       "pmtresponse"
 *       )'
 *     
 * The patterns are expanded by ROOT (by `TChain::Add()`), which is not that
 * flexible. For example, it can't expand wildcards in the paths. So another
 * approach is to use file lists:
 *     
 *     ls /pnfs/icarus/scratch/user/petrillo/jobOutput/photonlibrary_builder_icarus/20200806/output/job*[0-9]/out/[0-9]*_0/Supplemental-*.root | pnfsToXRootD > photonlibrary_builder_icarus-20200806
 *     root -l -q 'MergePhotonLibrary.C++(
 *       "PhotonLibraryTest.root",
 *       "photonli
 * brary_builder_icarus-test20200806.filelist",
 *       "v08_61_00",
 *       "2020test1",
 *       "PhotonLibraryData",
 *       "pmtresponse"
 *       )'
 *     
 * Needless to say: reading the tree chain across hundreds of files takes some
 * time (especially in the latter case where they are not local), and `TChain`
 * does not offer any visual feedback on console when copy is ongoing
 * (for example, the latter command reported
 * `85779522 entries transferred in 481.944 seconds`).
 * 
 */

// ROOT libraries
#include "Compression.h" // ROOT::kLZ4
#include "TFile.h"
#include "TDirectory.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TObjArray.h"
#include "TTree.h"
#include "TBranch.h"
#include "TNamed.h"
#include "TStopwatch.h"
#include "RooInt.h"
#include "RooDouble.h"

// C/C++ libraries
#include <iostream>
#include <fstream>
#include <set>
#include <string>
#include <vector>
#include <memory> // std::make_unique()
#include <utility> // std::pair
#include <cassert>


// -----------------------------------------------------------------------------
static std::string const ScriptName { "MergePhotonLibrary.C" };
static unsigned int const ScriptVersion { 5 };

static unsigned int MinMissingVoxelBlockReport = 10U;


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


// -----------------------------------------------------------------------------

/// Structure containing all supported metadata.
struct MetadataSet_t {
  std::map<std::string, std::unique_ptr<TNamed>> index;
}; // MetadataSet_t


/// Returns whether the string `s` ends with the specified `suffix`.
bool endsWith(std::string const& s, std::string const& suffix);


/// Writes the content of `metadata` into the console output stream `out`.
template <typename Stream>
void dumpMetadata(Stream& out, MetadataSet_t const& metadata);


/// Returns a string with the name of the specified compression algorithm.
std::string CompressionAlgorithmName(ROOT::ECompressionAlgorithm algo);


/// Removes the last element from the specified path.
std::string ROOTdirectoryOf(std::string const& ROOTobjectPath);


/**
 * @brief Returns the path of the tree object for all trees in the `chain`.
 * @param chain the ROOT tree chain to extract source paths from
 * @return a list of pairs: file name `first`, tree path in ROOT `second`
 */
std::vector<std::pair<std::string, std::string>> extractSourceFilePathsFromChain
  (TChain const& chain);


/// Returns the list of metadata contained in the specified ROOT directory.
MetadataSet_t extractMetadata(TDirectory& dir);


/// Returns a pair of `TFile` and `TDirectory` from the specified path.
std::pair<std::unique_ptr<TFile>, TDirectory*> openROOTdirectory(
  std::string const& filePath, std::string const& dirPath,
  std::string const& mode = "READ"
  );


/**
 * @brief Merges the metadata from `src` into `dest`
 * @return the number of encountered errors (`0` means success)
 */
unsigned int MergeMetadata
  (MetadataSet_t& dest, MetadataSet_t const& src, std::string const& srcName);


/// Writes the content of the specified `metadata` into ROOT output directory.
void writeMetadata(MetadataSet_t const& metadata, TDirectory& outDir);


/**
 * @brief Extracts and copies metadata from the input files of `tree` chain.
 * @param tree the chain of trees that made input of the photon library
 * @param outFile the ROOT directory where to write the metadata information
 * @return the collected metadata
 */
MetadataSet_t collectMetadata(TChain& tree);


/// Extracts the number of expected voxels from the metadata.
unsigned int extractNVoxels(MetadataSet_t const& metadata);


/// Extracts the list of voxels from the library and reports any missing ones.
/// @return the number of voxels known to be missing
unsigned int voxelCheck
  (TTree& tree, MetadataSet_t const& metadata, unsigned int minBlockSize = 1U);


/// Adds to a `tree` chain all the files matching a pattern (matching by ROOT).
unsigned int populateTreeFromPattern(
  TChain& tree, std::string const& inputPattern,
  std::string const& treeDir, std::string const& treeName
  );


/// Adds to a `tree` chain all the files from the specified file list.
/// Lines starting with `#` are ignored.
unsigned int populateTreeFromFileList(
  TChain& tree, std::string const& fileListPath,
  std::string const& treeDir, std::string const& treeName
  );



/** ****************************************************************************
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
 * * `4`: missing voxels (may be physical... check if it's in regular blocks!)
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
  
  int ErrorCode = 0; // tracks non-fatal errors that will be reported at the end
  
  //
  // parse the name
  //
  
  //
  // open source
  //
  TChain* pSourceTree = new TChain(inputTreeName.c_str(), "Photon library");
  
  bool const bInputFromList = endsWith(inputFilePattern, ".list")
    || endsWith(inputFilePattern, ".filelist");
  
  unsigned int const nFiles =
    (bInputFromList? populateTreeFromFileList: populateTreeFromPattern)
      (*pSourceTree, inputFilePattern, inputTreeDir, inputTreeName)
    ;
  if (nFiles == 0) {
    std::cerr << "No file matched '" << inputFilePattern << "'" << std::endl;
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
  auto pSource = new TNamed("Source", inputFilePattern.c_str());
  pSource->Write();

  TNamed{ "MergerProgram", ScriptName.c_str() }.Write();
  TNamed{ "MergerVersion", std::to_string(ScriptVersion).c_str() }.Write();
  
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
  
  std::cout << pDestTree->GetEntries() << " entries transferred in "
    << timer.RealTime() << " seconds; compression factors:" << std::endl;
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
  // deal with metadata; errors are not "fatal"
  //
  timer.Start();
  MetadataSet_t metadata = collectMetadata(*pSourceTree);
  writeMetadata(metadata, *pDestDir);
  timer.Stop();
  
  std::cout << "Metadata extracted in " << timer.RealTime() << " seconds."
    << std::endl;
  
  //
  // collect the list of voxels;
  // report missing blocks only if 10-voxel or larger
  //
  unsigned int const nLargeMissingBlocks
    = voxelCheck(*pDestTree, metadata, MinMissingVoxelBlockReport);
  if (nLargeMissingBlocks > 0U) ErrorCode = 4;
  
  //
  // close and go
  //
  pFDest->Close();
  
  return ErrorCode;
  
} // MergePhotonLibrary()


// -----------------------------------------------------------------------------
// ---  implementations follow
// -----------------------------------------------------------------------------
bool endsWith(std::string const& s, std::string const& suffix) {
  
  return (s.length() >= suffix.length())
    && (s.rfind(suffix) == (s.length() - suffix.length()));
  
} // endsWith()


// -----------------------------------------------------------------------------
template <typename Stream>
void printMetadataValue(Stream& out, TNamed const& obj) {
  
  if (dynamic_cast<RooInt const*>(&obj)) {
    out << Int_t(static_cast<RooInt const&>(obj));
  }
  else if (dynamic_cast<RooDouble const*>(&obj)) {
    out << Double_t(static_cast<RooDouble const&>(obj));
  }
  else {
    out << obj.GetTitle();
  }
  
} // printMetadataValue()


template <typename Stream>
void dumpMetadata(Stream& out, MetadataSet_t const& metadata) {
  
  out << "Metadata contains " << metadata.index.size() << " values:";
  for (auto const& [ name, obj ]: metadata.index) {
    
    out << "\n  '" << obj->GetName() << "': ";
    printMetadataValue(out, *obj);
    if (name != obj->GetName()) out << " (as \"" << name << "\")";
    
  }
  out << std::endl;
  
} // dumpMetadata()


// -----------------------------------------------------------------------------
std::string CompressionAlgorithmName(ROOT::ECompressionAlgorithm algo);


/// Removes the last element from the specified path.
std::string ROOTdirectoryOf(std::string const& ROOTobjectPath) {
  
  auto const iSep = ROOTobjectPath.rfind('/');
  return (iSep == std::string::npos)
    ? std::string{}: ROOTobjectPath.substr(0, iSep);
  
} // ROOTdirectoryOf()


// -----------------------------------------------------------------------------
std::vector<std::pair<std::string, std::string>> extractSourceFilePathsFromChain
  (TChain const& chain)
{
  std::vector<std::pair<std::string, std::string>> paths;
  paths.reserve(chain.GetListOfFiles()->GetEntries());
  for (TObject const* pObj: *(chain.GetListOfFiles())) {
    auto const& elem = dynamic_cast<TChainElement const&>(*pObj);
    
    // this is documented in `TChain::AddFile()` (ROOT 6.20)
    std::string treePath = elem.GetName();
    std::string filePath = elem.GetTitle();
    
    paths.emplace_back(std::move(filePath), std::move(treePath));
  } // for files
  
  return paths;
} // extractSourceFilePathsFromChain()


// -----------------------------------------------------------------------------
MetadataSet_t extractMetadata(TDirectory& dir) {
  MetadataSet_t metadata;
  
  for (auto obj: *(dir.GetListOfKeys())) {
    
    auto& key = dynamic_cast<TKey&>(*obj);
    
    std::string const className = key.GetClassName();
    
    TClass const* ROOTclass = TClass::GetClass(className.c_str());
    if (!ROOTclass) {
      // this is not considered an error
      std::cerr
        << "WARNING: in '" << dir.GetPath() << "': object '" << key.GetName()
        << "' is of type '" << className << "' unknown to ROOT."
        << std::endl;
      continue;
    }
    if (!ROOTclass->InheritsFrom(TNamed::Class())) {
      std::cerr << "WARNING: in '" << dir.GetPath() << "': object '"
        << key.GetName() << "' of type '" << className
        << "' is not supported metadata."
        << std::endl;
      continue;
    }
    
    if (ROOTclass->InheritsFrom(TTree::Class())) continue; // photon library?
    
    metadata.index[key.GetName()] = std::unique_ptr<TNamed>
      { static_cast<TNamed*>(key.ReadObject<TNamed>()->Clone()) };
    
  } // for
  
//   std::cout << "Collected metadata:\n";
//   dumpMetadata(std::cout, metadata);
  
  return metadata;
} // extractMetadata()


// -----------------------------------------------------------------------------
std::pair<std::unique_ptr<TFile>, TDirectory*> openROOTdirectory(
  std::string const& filePath, std::string const& dirPath,
  std::string const& mode /* = "READ" */
) {
  TgDirectoryGuard dirGuard;
  
  auto file
    = std::unique_ptr<TFile>(TFile::Open(filePath.c_str(), mode.c_str()));
  if (!file || !file->IsOpen()) return { nullptr, nullptr };
  
  TDirectory* const dir
    = dirPath.empty()? file.get(): file->GetDirectory(dirPath.c_str());
  
  return { std::move(file), dir };
  
} // openROOTdirectory()


// -----------------------------------------------------------------------------
unsigned int MergeMetadata(
  MetadataSet_t& dest, MetadataSet_t const& src, std::string const& srcName
) {
  
  unsigned int nErrors = 0U;
  
  for (auto const& [ key, obj ]: src.index) {
    
    auto iDest = dest.index.find(key);
    if (iDest == dest.index.end()) {
      /*
      std::cout
        << "New metadata element '" << key << "' (integral) found in '"
        << srcName << "'." << std::endl;
      */
      dest.index.emplace
        (key, std::unique_ptr<TNamed>{ static_cast<TNamed*>(obj->Clone()) });
    }
    else {
      
      if (obj->InheritsFrom(RooInt::Class())) {
        
        auto const value = Int_t(static_cast<RooInt const&>(*obj));
        auto const metaValue
          = Int_t(static_cast<RooInt const&>(*(iDest->second)));
        if (metaValue != value) {
          std::cerr << "ERROR: metadata '" << key << "' (RooInt) from '"
            << srcName << "' has value " << value << " incompatible with "
            << metaValue << " from the previous input."
            << std::endl;
          ++nErrors;
          continue;
        }
      }
      else if (obj->InheritsFrom(RooDouble::Class())) {
        
        auto const value = Double_t(static_cast<RooDouble const&>(*obj));
        auto const metaValue
          = Double_t(static_cast<RooDouble const&>(*(iDest->second)));
        if (metaValue != value) {
          std::cerr << "ERROR: metadata '" << key << "' (RooDouble) from '"
            << srcName << "' has value " << value << " incompatible with "
            << metaValue << " from the previous input."
            << std::endl;
        }
      }
      else {
        
        std::string const value = obj->GetTitle();
        std::string const metaValue = iDest->second->GetTitle();
        if (metaValue != value) {
          std::cerr << "ERROR: metadata '" << key << "' from '"
            << srcName << "' has value '" << value << "' incompatible with '"
            << metaValue << "' from the previous input."
            << std::endl;
        }
      }
      
    }
    
  } // for all metadata
  
  return nErrors;
} // MergeMetadata()


// -----------------------------------------------------------------------------
void writeMetadata(MetadataSet_t const& metadata, TDirectory& outDir) {
  
  TgDirectoryGuard dirChanger(outDir);
  for (auto const& [ key, value ]: metadata.index) {
    value->Write();
  } // for integral metadata
  std::cout << metadata.index.size() << " metadata entries written into '"
    << outDir.GetPath() << "'." << std::endl;
  
  dumpMetadata(std::cout, metadata);
  
} // writeMetadata()


// -----------------------------------------------------------------------------
MetadataSet_t collectMetadata(TChain& tree) {
  
  std::cout << "Parsing the tree chain for metadata..." << std::endl;
  
  auto const& sourceList = extractSourceFilePathsFromChain(tree);
  std::cout << "   ... extracted " << sourceList.size() << " sources"
    << std::endl;
  
  unsigned int nErrors = 0U;
  
  //
  // collect and merge
  //
  
  // do not try to be overly generic...
  MetadataSet_t globalMetadata;
  
  for (auto const& [ filePath, objPath ]: sourceList) {
    
    auto&& [ sourceFile, sourceDir ]
      = openROOTdirectory(filePath, ROOTdirectoryOf(objPath));
    if (!sourceDir) {
      std::cerr
        << "Failed to open '" << filePath << "/" << ROOTdirectoryOf(objPath)
        << "'" << std::endl;
      ++nErrors;
      continue;
    }
    
    MetadataSet_t sourceMetadata = extractMetadata(*sourceDir);
    
    nErrors += MergeMetadata
      (globalMetadata, sourceMetadata, sourceDir->GetPath());
    
  } // for
  
  return globalMetadata;
  
} // collectMetadata()


// -----------------------------------------------------------------------------
unsigned int extractNVoxels(MetadataSet_t const& metadata) {
  
  auto readInt = [&md = metadata.index](std::string const& key)
    {
      auto const iMeta = md.find(key);
      if (iMeta == md.end()) return 0;
      auto metaObj = dynamic_cast<RooInt const*>(iMeta->second.get());
      return metaObj? int(Int_t(*metaObj)): 0;
    };
  
  int NVoxels = readInt("NVoxels");
  if (NVoxels <= 0) {
    NVoxels = readInt("NDivX") * readInt("NDivY") * readInt("NDivZ");
  }
  return NVoxels;
  
} // extractNVoxels()



// -----------------------------------------------------------------------------
unsigned int voxelCheck(
  TTree& tree, MetadataSet_t const& metadata,
  unsigned int minBlockSize /* = 1U */
) {
  
  /*
   * first collect all the observed voxel numbers (one by one!),
   * then sort them and find, count and report any gap
   */
  TStopwatch timer;
  
  //
  // collect the number of all voxels in the library (automatically sorted)
  //
  auto const nEntries = tree.GetEntriesFast();
  std::cout << "Extracting the list of voxels from " << nEntries
    << " entries in the merged photon library..." << std::endl;
  std::set<int> voxelsFound;
  
  std::string const VoxelBranchName = "Voxel";
  
  Int_t branchVoxel = -1;
  tree.SetBranchStatus("*", false);
  tree.SetBranchStatus(VoxelBranchName.c_str(), true);
  tree.SetBranchAddress(VoxelBranchName.c_str(), &branchVoxel);
  Long64_t iEntry = 0;
  timer.Start();
  while (tree.GetEntry(iEntry++) > 0) voxelsFound.insert(branchVoxel);
  timer.Stop();
  if (--iEntry != nEntries) {
    std::cerr << "ERROR: " << nEntries << " entries expected in the tree, but "
      << iEntry << " were read." << std::endl;
  }
  else {
    std::cout << nEntries << " entries read from the tree in "
      << timer.RealTime() << " seconds." << std::endl;
  }
  
  //
  // find the gaps; if no metadata is provided, missing voxels at the end might
  // pass undetected
  //
  int NExpectedVoxels = extractNVoxels(metadata);
  if (NExpectedVoxels > 0)
    std::cout << NExpectedVoxels << " voxels are expected." << std::endl;
  else
    std::cout << "The total number of voxels is not known." << std::endl;    
  
  unsigned int nMissingVoxels = 0U;
  unsigned int nMissingBlocks = 0U;
  unsigned int nLargeMissingBlocks = 0U;
  int firstMissing = 0;
  for (int const voxel: voxelsFound) {
    
    if (voxel == firstMissing) { // not actually missing...
      ++firstMissing;
      continue;
    }
    else {
      auto const missingInBlock
        = static_cast<unsigned int>(voxel - firstMissing);
      if (missingInBlock >= minBlockSize) {
        std::cerr << "Missing voxels: " << firstMissing;
        if (missingInBlock > 1)
          std::cerr << " to " << (voxel - 1) << " (" << missingInBlock << ")";
        std::cerr << std::endl;
        ++nLargeMissingBlocks;
      }
      nMissingVoxels += missingInBlock;
      ++nMissingBlocks;
    }
    firstMissing = voxel + 1;
  } // for
  if (firstMissing < NExpectedVoxels) {
    auto const missingInBlock
      = static_cast<unsigned int>(NExpectedVoxels - firstMissing);
    if (missingInBlock >= minBlockSize) {
      std::cerr << "Missing voxels: " << firstMissing;
      if (missingInBlock > 1) {
        std::cerr << " to " << (NExpectedVoxels - 1)
          << " (" << missingInBlock << ")";
      }
      std::cerr << std::endl;
      ++nLargeMissingBlocks;
    }
    nMissingVoxels += missingInBlock;
    ++nMissingBlocks;
  }
  
  if (nMissingVoxels > 0U) {
    std::cerr << " => " << nMissingVoxels << " voxels missing in "
      << nMissingBlocks << " blocks (" << nLargeMissingBlocks
      << " blocks at least " << minBlockSize << " voxel big)!" << std::endl;
  }
  
  return nLargeMissingBlocks;
} // voxelCheck()



// -----------------------------------------------------------------------------
unsigned int populateTreeFromPattern(
  TChain& tree, std::string const& inputPattern,
  std::string const& treeDir, std::string const& treeName
) {
  
  std::string pattern = inputPattern;
  if (!treeDir.empty()) pattern += '/' + treeDir;
  if (!treeName.empty()) pattern += '/' + treeName;
  
  std::cout << "Input files matching pattern: '" << pattern << "'" << std::endl;
  return tree.Add(pattern.c_str());
  
} // populateTreeFromPattern()


// -----------------------------------------------------------------------------
unsigned int populateTreeFromFileList(
  TChain& tree, std::string const& fileListPath,
  std::string const& treeDir, std::string const& treeName
) {
  /*
   * File list format:
   *  * the entire content of each line is the file name, including any space
   *    anywhere
   *  * empty lines and lines whose first character is a '#' are skipped
   * White space treatment could be more flexible. But it is not.
   */
  
  std::string localTreePath;
  if (!treeDir.empty()) localTreePath += '/' + treeDir;
  if (!treeName.empty()) localTreePath += '/' + treeName;

  std::ifstream fileList(fileListPath);
  if (!fileList.is_open()) {
    std::cerr << "ERROR: can't open file list '" << fileListPath << "'"
      << std::endl;
    return 0;
  }
  
  unsigned int nFiles = 0;
  unsigned int LineNo = 0;
  std::string fileName;
  while (std::getline(fileList, fileName)) {
    ++LineNo;
    if (fileName.empty()) continue;
    if (fileName.front() == '#') continue;
    
    if (!localTreePath.empty()) fileName += localTreePath;
    
    // std::cout << "Adding [#" << nFiles << "] '" << fileName
    //   << "' (at line " << LineNo << ")" << std::endl;
    int const added = tree.AddFile(fileName.c_str());
    if (added == 0) {
      std::cerr << "Failed to add tree from file '" << fileName
        << "' (filelist line " << LineNo << ")" << std::endl;
      continue;
    }
    else nFiles += added;
    
  } // while
  
  return nFiles;
  
} // populateTreeFromFileList()



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

