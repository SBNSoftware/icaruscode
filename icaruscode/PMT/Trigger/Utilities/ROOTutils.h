/**
 * @file   icaruscode/PMT/Trigger/Utilities/ROOTutils.h
 * @brief  A bunch of diverse utilities and futilities related to ROOT.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   April 1, 2019
 * @see    `icaruscode/PMT/Trigger/Utilities/ROOTutils.tcc`
 */

#ifndef ICARUSCODE_PMT_TRIGGER_UTILITIES_ROOTUTILS_H
#define ICARUSCODE_PMT_TRIGGER_UTILITIES_ROOTUTILS_H

// LArSoft libraries
#include "larcorealg/CoreUtils/enumerate.h"

// ROOT libraries
#include "TStyle.h"
#include "TDirectory.h"
#include "TAxis.h"

// C/C++ libraries
#include <utility> // std::pair<>
#include <vector>
#include <string>
#include <cassert>


namespace util::ROOT {
  
  // --- BEGIN -- Global state changers ----------------------------------------
  /// @name Global state changers
  /// @{
  // ---------------------------------------------------------------------------
  /**
   * @brief A class restoring the previous `TDirectory` on destruction.
   * 
   * When an instance of this object is created, the existing current directory
   * is saved, and it is then restored on destruction.
   * 
   * Additional methods allow finer control on the restoration feature.
   * 
   * Example:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * {
   *   TDirectoryChanger DirGuard(gDirectory->GetDirectory("subdir"));
   *   
   *   // everything here happens in the (existing) `subdir` directory.
   * }
   * 
   * // whatever follows happens under the previous directory
   * 
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * 
   */
  class TDirectoryChanger {
    
    TDirectory* pSaved = nullptr;
    TDirectory* pNew = nullptr;
    
      public:
    
    TDirectoryChanger() { save(); }
    TDirectoryChanger(TDirectory* dir) { save(); cd(dir); }
    TDirectoryChanger(std::string const& dir, std::string const& title = "")
      { save(); cd(dir, title); }
    
    ~TDirectoryChanger()      { restore(); }
    
    /// Stores the current directory as the one to be saved.
    void save()               { pSaved = currentDir(); }
    
    /// Immediately restores the old directory.
    /// It will still restored on destruction too.
    void restore() const      { if (pSaved) pSaved->cd(); }
    
    /// Do not restore the old directory on destruction.
    void forget()             { pSaved = nullptr; }
    
    /// Make the stored new directory as current again.
    void cd() const           { if (pNew) pNew->cd(); }
    
    /// Make the specified directory as current.
    void cd(TDirectory* dir)  { pNew = dir; cd(); }
    
    /// Make the specified directory as current, possibly creating it.
    void cd(std::string const& name, std::string const& title = "")
      {
        assert(currentDir());
        pNew = currentDir()->GetDirectory(name.c_str());
        if (!pNew) pNew = currentDir()->mkdir(name.c_str(), title.c_str());
        cd(); 
      }
    
    /// Returns a pointer to the directory that will be restored on destruction.
    TDirectory* saved() const { return pSaved; }
    
    /// Returns whether there is a directory to be restored on destruction.
    bool hasSaved() const     { return saved() != nullptr; }
    
    static TDirectory* currentDir() { return gDirectory; }
    
  }; // class TDirectoryChanger
  
  
  // ---------------------------------------------------------------------------
  /**
   * @brief A class restoring the previous `TStyle` on destruction.
   * 
   * When an instance of this object is created, the existing current style is
   * saved, and it is then restored on destruction.
   * 
   * Additional methods allow finer control on the restoration feature.
   * 
   * Example:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * {
   *   TStyleChanger StyleGuard
   *     (static_cast<TStyle*>(gDirectory->Get("newStyle")));
   *   
   *   // everything here happens in the new style
   * }
   * 
   * // whatever follows happens under the previous style
   * 
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * 
   */
  class TStyleChanger {
    
    TStyle* pSaved = nullptr;
    TStyle* pNew = nullptr;
    
      public:
    
    TStyleChanger() { save(); }
    TStyleChanger(TStyle* newStyle) { save(); cd(newStyle); }
    
    ~TStyleChanger() { restore(); }
    
    /// Stores the current style as the one to be saved.
    void save()               { pSaved = gStyle; }
    
    /// Immediately restores the old style.
    /// It will still restored on destruction too.
    void restore() const      { if (pSaved) pSaved->cd(); }
    
    /// Do not restore the old style on destruction.
    void forget()             { pSaved = nullptr; }
    
    /// Make the stored new style as current again.
    void cd() const           { if (pNew) pNew->cd(); }
    
    /// Make the specified style as current.
    void cd(TStyle* newStyle) { pNew = newStyle; cd(); }
    
    /// Returns a pointer to the style that will be restored on destruction.
    TStyle* saved() const     { return pSaved; }
    
    /// Returns whether there is a style to be restored on destruction.
    bool hasSaved() const     { return saved() != nullptr; }
    
  }; // class TStyleChanger
  
  /// @}
  // --- END -- Global state changers ------------------------------------------
  
  
  // ---------------------------------------------------------------------------
  /**
   * @brief Sets all the labels starting with the bin `first` (`1` by default).
   * @param pAxis the axis to label
   * @param labels the labels of each bin
   * @param first (default: `1`) start from this bin with the first label
   */
  inline void applyAxisLabels
    (TAxis* pAxis, std::vector<std::string> const& labels, int first = 1)
    {
      for (auto&& [iLabel, label]: util::enumerate(labels))
        pAxis->SetBinLabel(first + iLabel, label.c_str());
    } // applyAxisLabels()
  
  
  // ---------------------------------------------------------------------------
  /**
   * @brief Returns a variable size binning for the points.
   * @tparam Coll type of collection of the points
   * @param centralPoints set of points to build the bins around
   * @return a pair of bin edge and label collections
   * 
   * This algorithm attempts to create a variable width binning so that all
   * specified points fall in the middle of their respective bin.
   * 
   * @note This algorithm is known to be prone to failure for some point
   *       distributions.
   * 
   * The return value is a pair of collections. The first one is a sequence of
   * bin boundaries, starting with the lower edge of the bin for the first
   * point and ending with the upper edge of the last point. Therefore this
   * first collection has a number of element larger than the number of points
   * by one unit.
   * The second collection is a string representing the name of the bin.
   * The first element of this collection is the label of the first bin, that
   * is the label of the bin whose lower edge is the first element of the other
   * collection. The second collection has the same number of elements as the
   * number of points, and smaller by one unit than the edge collection.
   * 
   * Requirements
   * -------------
   * 
   * * the collection `Coll` must be 
   * * the type (`T`) of values in `Coll` needs to support the following
   *   operations:
   *     * convertible to double: `operator double (T)`
   *     * conversion to string (`to_string(T)`)
   * * `centralPoint` must contain at least two points
   * 
   */
  template <typename Coll>
  std::pair<std::vector<double>, std::vector<std::string>>
  makeVariableBinningAndLabels(Coll const& centralPoints);
  
  
  //----------------------------------------------------------------------------
  
} // namespace util::ROOT


//------------------------------------------------------------------------------
//--- template implementation
//------------------------------------------------------------------------------
#include "icaruscode/PMT/Trigger/Utilities/ROOTutils.tcc"

//------------------------------------------------------------------------------


#endif // ICARUSCODE_PMT_TRIGGER_UTILITIES_ROOTUTILS_H
