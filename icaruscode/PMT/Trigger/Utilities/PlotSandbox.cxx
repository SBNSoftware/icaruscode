/**
 * @file   useranalysis/Trigger/Design/Utilities/PlotSandbox.cxx
 * @brief  A helper to manage ROOT objects in a `art::TFileDirectory`.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   August 8, 2019
 * @see    `useranalysis/Trigger/Design/Utilities/PlotSandbox.h`
 */


// library header
#include "icaruscode/PMT/Trigger/Utilities/PlotSandbox.h"

// LArSoft libraries
#include "larcorealg/CoreUtils/values.h"

// framework libraries
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h" // MF_XXX() macros

// ROOT libraries

// C/C++ standard libraries
#include <string_view>
#include <utility> // std::forward()
#include <type_traits> // std::add_const_t<>


//------------------------------------------------------------------------------
//--- icarus::trigger::PlotSandbox::TFileDirectoryHelper
//------------------------------------------------------------------------------
auto icarus::trigger::PlotSandbox::TFileDirectoryHelper::create
  (art::TFileDirectory dir) -> TFileDirectoryHelper
{
  
  /*
   * Finding the ROOT directory is going to be tricky, since the interface
   * includes: `mkdir()`, `make()` and `makeAndRegister()`. Woah.
   * 
   * So the plan is:
   * 1. have `art::TFileDirectory` `make()` a new ROOT object in its directory;
   *    let that ROOT object be a `TDirectory`
   * 2. `TDirectory` is placed by default in the current directory (as usual),
   *    and it knows which its parent directory is
   * 3. after learning which that parent directory is, we delete the directory
   *    we just created; this also updates the mother directory
   */
  
  static constexpr const char* TestDirName = " PlotSandbox invalid name! ";
  
  // even if another directory with this name existed, it would not be replaced,
  // but rather another one would be created with this name:
  TDirectory* testDir = dir.make<TDirectory>(TestDirName, TestDirName);
  if (!testDir) {
    throw cet::exception("PlotSandbox") << "TFileDirectoryHelper::create() "
      "failed to figure out the ROOT directory!\n";
  }
  TDirectory* pROOTdir = testDir->GetMotherDir();
  MF_LOG_DEBUG("TFileDirectoryHelper")
    << "icarus::trigger::PlotSandbox::TFileDirectoryHelper::create(): "
    << "found parent directory: '" << pROOTdir->GetName() << "'";
  // ... and even if there are multiple directories with the same name, using
  // the pointer to the object makes this deletion affect only the object itself
  delete testDir;
  
  return { dir, pROOTdir };
} // icarus::trigger::PlotSandbox::TFileDirectoryHelper::create()


//------------------------------------------------------------------------------
//--- icarus::trigger::PlotSandbox
//------------------------------------------------------------------------------
void icarus::trigger::PlotSandbox::Data_t::resetSubboxParents
  (PlotSandbox const* newParent)
{
  for (auto& subbox: util::values(subBoxes)) subbox->setParent(newParent);
} // icarus::trigger::PlotSandbox::Data_t::resetSubboxParents()


//------------------------------------------------------------------------------
icarus::trigger::PlotSandbox::PlotSandbox(
  art::TFileDirectory parentDir,
  std::string name, std::string desc
  )
  : fData(std::move(name), std::move(desc), 
    name.empty()
      ? TFileDirectoryHelper::create(parentDir)
      : TFileDirectoryHelper::create(parentDir, name, desc)
    )
{}


//------------------------------------------------------------------------------
icarus::trigger::PlotSandbox::PlotSandbox(PlotSandbox&& from)
  : fData(std::move(from.fData))
{
  fData.resetSubboxParents(this); // adjust the parent pointers of the subboxes
} // PlotSandbox::PlotSandbox(PlotSandbox&&)


//------------------------------------------------------------------------------
std::string icarus::trigger::PlotSandbox::ID() const
  { return fData.parent? (fData.parent->ID() + '/' + name()): name(); }


//------------------------------------------------------------------------------
std::string icarus::trigger::PlotSandbox::processName
  (std::string const& name) const
{
  std::string const sandboxName = processedSandboxName();
  return sandboxName.empty()? name: name + '_' + sandboxName;
} // icarus::trigger::PlotSandbox::processName()


//------------------------------------------------------------------------------
std::string icarus::trigger::PlotSandbox::processTitle
  (std::string const& title) const
{
  std::string const sandboxDesc = processedSandboxDesc();
  return sandboxDesc.empty()? title: title + ' ' + sandboxDesc;
} // icarus::trigger::PlotSandbox::processTitle()


//------------------------------------------------------------------------------
auto icarus::trigger::PlotSandbox::findSandbox(std::string const& name)
  -> PlotSandbox*
  { return findSandbox(*this, name); }

auto icarus::trigger::PlotSandbox::findSandbox(std::string const& name) const
  -> PlotSandbox const*
  { return findSandbox(*this, name); }


//------------------------------------------------------------------------------
auto icarus::trigger::PlotSandbox::demandSandbox(std::string const& name)
  -> PlotSandbox&
  { return demandSandbox(*this, name); }

auto icarus::trigger::PlotSandbox::demandSandbox(std::string const& name) const
  -> PlotSandbox const&
  { return demandSandbox(*this, name); }


//------------------------------------------------------------------------------
icarus::trigger::PlotSandbox::PlotSandbox
  (PlotSandbox const& parent, std::string name, std::string desc)
  : PlotSandbox(parent.fData.outputDir.fDir, std::move(name), std::move(desc))
{
  setParent(&parent);
}

//------------------------------------------------------------------------------
std::string icarus::trigger::PlotSandbox::processPlotTitle
  (std::string const& title) const
{
  /*
   * We want to process the title of the plot, but this "title" string might
   * contain also axis labels (see e.g. `TH1::SetTitle()`).
   * 
   * We need to first identify the actual title portion of the "title",
   * then process it alone and merge it back (one Python line...).
   * 
   */
  
  // eat title until an unescaped ';' is found
  auto const tbegin = title.begin();
  auto const tend = title.end();
  auto atend = tbegin; // actual title end
  while (atend != tend) {
    if ((*atend == ';') && ((atend == tbegin) || (*std::prev(atend) != '\\')))
      break;
    ++atend;
  } // while
  
  return processTitle({ tbegin, atend }).append(atend, tend);
  
} // icarus::trigger::PlotSandbox::processedPlotTitle()


//------------------------------------------------------------------------------
std::string icarus::trigger::PlotSandbox::processedSandboxName() const
{
  std::string processed;
  
  // if there is no name, the processed name is empty (we don't recurse parents)
  if (!hasName()) return processed;
  
  processed = name();
  if (fData.parent) processed += '_' + fData.parent->processedSandboxName();
  
  return processed;
} // icarus::trigger::PlotSandbox::processedSandboxName()


//------------------------------------------------------------------------------
std::string icarus::trigger::PlotSandbox::processedSandboxDesc() const
{
  std::string processed;
  
  // if there is no name, the processed name is empty (we don't recurse parents)
  if (hasDescription()) processed += description();
  
  if (fData.parent) processed += ' ' + fData.parent->processedSandboxDesc();
  
  return processed;
} // icarus::trigger::PlotSandbox::processedSandboxDesc()


//------------------------------------------------------------------------------
std::pair<std::string, std::string> icarus::trigger::PlotSandbox::splitPath
  (std::string const& path, char sep /* = '/' */)
{
  
  auto const iSep = path.rfind(sep);
  if (iSep == std::string::npos) return { {}, path };
  else return { path.substr(0U, iSep), path.substr(iSep + 1) };
  
} // icarus::trigger::PlotSandbox::splitPath()


//------------------------------------------------------------------------------
std::string icarus::trigger::PlotSandbox::joinPath
  (std::initializer_list<std::string> pathElements, char sep /* = '/' */)
{
  if (std::empty(pathElements)) return {};
  
  auto stripSep = [sep](std::string const& s) -> std::string_view {
    return {
      s.data(),
      ((s.length() > 1) && (s.back() == sep))? s.length() - 1: s.length()
      };
    };
  
  auto iElem = pathElements.begin();
  auto const eend = pathElements.end();
  auto const& first = stripSep(*iElem);
  std::string s { first.begin(), first.end() };
  while (++iElem != eend) {
    auto const& elem = stripSep(*iElem);
    if (elem.empty()) continue;
    if (elem.front() != sep) s.push_back(sep);
    s += elem;
  } // while
  return s;
} // icarus::trigger::PlotSandbox::joinPath()


//------------------------------------------------------------------------------
