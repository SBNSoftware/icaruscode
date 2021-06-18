/**
 * @file   icaruscode/PMT/Trigger/Utilities/PlotSandbox.h
 * @brief  A helper to manage ROOT objects in a `art::TFileDirectory`.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   August 8, 2019
 * @see    `icaruscode/PMT/Trigger/Utilities/PlotSandbox.cxx`
 */

#ifndef ICARUSCODE_PMT_TRIGGER_UTILITIES_PLOTSANDBOX_H
#define ICARUSCODE_PMT_TRIGGER_UTILITIES_PLOTSANDBOX_H

// LArSoft libraries
#include "larcorealg/CoreUtils/span.h" // util::make_transformed_span(), ...

// framework libraries
#include "art_root_io/TFileDirectory.h"
#include "cetlib_except/exception.h"

// ROOT libraries
#include "TDirectory.h"
#include "TDirectoryFile.h"
#include "TList.h"
#include "TKey.h"
#include "TClass.h"
#include "TObject.h"

// C/C++ standard libraries
#include <string>
#include <map>
#include <iterator> // std::prev()
#include <utility> // std::pair<>
#include <functional> // std::hash<>
#include <initializer_list>
#include <memory> // std::unique_ptr<>


//------------------------------------------------------------------------------
namespace icarus::trigger::details {
  
  // TODO it seems make_transformed_span() does not support nesting properly;
  //      this should be fixed
  template <typename Map>
  struct map_dereferenced_values_impl {
    
    static constexpr std::size_t NElement = 1U; // this is the value
    
    template <typename T>
    static constexpr decltype(auto) iterate(T&& coll) noexcept
      {
        auto extractor = [](auto&& value) -> decltype(auto)
          { return *std::get<NElement>(value); };
        return util::make_transformed_span(coll, extractor);
      }
    
  }; // map_values_impl
  
  template <typename Map>
  decltype(auto) map_dereferenced_values(Map&& map)
    { return map_dereferenced_values_impl<std::decay_t<Map>>::iterate(std::forward<Map>(map)); }
  
} // namespace util::details


//------------------------------------------------------------------------------
namespace icarus::trigger { class PlotSandbox; }
/**
 * @brief A helper to manage ROOT objects with consistent naming.
 *
 * A sandbox includes a ROOT directory where all the objects are written.
 * It also provides a name pattern to modify a generic object name in one
 * specific to this sandbox, e.g. from `"HEnergy"` to `"HMuonNeutrinoEnergy"`.
 * Object descriptions (usually titles) may also be processed in the same way.
 * 
 * In addition, the sandbox can point to a parent sandbox, in which case the
 * names and descriptions are first processed by this sandbox and then passed to
 * that parent sandbox for further processing.
 * 
 * The sandbox has two characterizing strings:
 * * the _name_ is a key used in ROOT object names and it is also by default
 *   the name of the output ROOT directory; for example: `"MuonNeutrino"` or
 *   `"CC"`. This string becoming part of the ROOT object name, it is better to
 *   keep it short and simply alphanumeric.
 * * the _description_ is a short wording for the content of the sandbox, used
 *   in processing ROOT object titles; for example, `"#nu_{#mu}"` or
 *   `"charged current"`.
 * 
 * @note By convention the subdirectory names are not processed.
 * 
 */
class icarus::trigger::PlotSandbox {
  
  /*
   * `art::TFileDirectory` does not give access to the underlying ROOT directory
   * and to be able to read from it we need the following trick:
   */
  
  /// Contains both a `art::TFileDirectory` and the `TDirectory` it manages.
  struct TFileDirectoryHelper {
    art::TFileDirectory fDir;
    TDirectory* fROOTdir;
    
    TFileDirectoryHelper(art::TFileDirectory dir, TDirectory* ROOTdir)
      : fDir(dir), fROOTdir(ROOTdir) {}
    
    /// Creates a helper managing a subdirectory of `parentDir`.
    template <typename RootDir = TDirectoryFile>
    static TFileDirectoryHelper create(
      art::TFileDirectory parentDir,
      std::string const& subdir, std::string const& dirTitle = ""
      );
    
    static TFileDirectoryHelper create(art::TFileDirectory dir);
    
  }; // struct TFileDirectoryHelper
  
  
  /// The whole data in a convenient package!
  struct Data_t {
    
    std::string name; ///< The name/key representing this sandbox.
    std::string desc; ///< The description characterizing the sandbox.
    
    PlotSandbox const* parent = nullptr; ///< Optional parent sandbox.
    
    /// Contained sand boxes.
    std::map<std::string, std::unique_ptr<PlotSandbox>> subBoxes;
    
    TFileDirectoryHelper outputDir; ///< Output ROOT directory of the sandbox.
    
    Data_t() = default;
    Data_t(Data_t const&) = delete;
    Data_t(Data_t&&) = default;
    Data_t& operator= (Data_t const&) = delete;
    Data_t& operator= (Data_t&&) = default;
    
    Data_t
      (std::string&& name, std::string&& desc, TFileDirectoryHelper outputDir)
      : name(std::move(name)), desc(std::move(desc))
      , outputDir(std::move(outputDir))
      {}
    
    void resetSubboxParents(PlotSandbox const* newParent);
    
  } fData;
  
  
  /// Helper function for `findSandbox()` implementations.
  template <typename SandboxType>
  static auto* findSandbox(SandboxType& sandbox, std::string const& name);
  
  /// Helper function for `demandSandbox()` implementations.
  template <typename SandboxType>
  static auto& demandSandbox(SandboxType& sandbox, std::string const& name);
  
  
    public:
  
  /**
   * @brief Constructor: specifies all sandbox characteristics.
   * @param parentDir ROOT directory under which the sandbox is created
   * @param name the name of the sandbox
   * @param desc description of the sandbox
   * 
   * If the `name` is empty, as a special case, the sandbox is unnamed and it is
   * created directly in `parentDir`. In that case, the title of the output
   * directory (`parentDir`) is also not changed.
   * If `name` or `desc` are empty, the pertaining processing is not performed:
   * if `name` is empty, ROOT object names will be unchanged, and if `desc` is
   * empty their descriptions/titles will be.
   * 
   * To create a sandbox with a parent, call `addSubSandbox()` of that parent.
   */
  PlotSandbox
    (art::TFileDirectory parentDir, std::string name, std::string desc);
  
  PlotSandbox(PlotSandbox const&) = delete;
  PlotSandbox(PlotSandbox&& from);
  
  // no assignment supported in `art::TFileDirectory`
  PlotSandbox& operator=(PlotSandbox const&) = delete;
  PlotSandbox& operator=(PlotSandbox&& from) = delete;
  
  
  /// Virtual destructor. Default, but C++ wants it.
  virtual ~PlotSandbox() = default;
  
  /// Returns whether we have a non-empty name.
  bool hasName() const { return !fData.name.empty(); }
  
  /// Returns the sandbox name.
  std::string const& name() const { return fData.name; }
  
  /// Returns whether we have a non-empty description.
  bool hasDescription() const { return !fData.desc.empty(); }
  
  /// Returns the sandbox description.
  std::string const& description() const { return fData.desc; }
  
  /// Returns the sandbox description for use at the beginning of a sentence.
  /// @todo Not implemented yet.
  virtual std::string const& Description() const { return description(); }
  
  /// Returns a string ID for this sandbox.
  std::string ID() const;
  
  /// Processes the specified string as it were a name.
  virtual std::string processName(std::string const& name) const;
  
  /// Processes the specified string as it were a description or title.
  virtual std::string processTitle(std::string const& title) const;
  
  
  // --- BEGIN -- ROOT object management ---------------------------------------
  /// @name ROOT object management
  /// @{
  
  /**
   * @brief Fetches the object with the specified name from the sandbox.
   * @tparam Obj (default: `TObject`) type of the object to fetch
   * @param name unprocessed name and path of the object to fetch
   * @return a constant pointer to the requested object , or `nullptr` if not
   *         available or wrong type
   * @see `use()`
   * 
   * The `name` specification may contain a ROOT path. The directory component
   * of the path, defined by everything preceding a `/` character, is _not_
   * processed.
   * 
   * The fetched object is converted to the desired type via `dynamic_cast`.
   * If conversion fails, a null pointer is returned.
   */
  template <typename Obj = TObject>
  Obj const* get(std::string const& name) const;
  
  /**
   * @brief Fetches an object with the specified name to be modified.
   * @tparam Obj (default: `TObject`) type of the object to fetch
   * @param name unprocessed name and path of the object to fetch
   * @return a pointer to the requested object , or `nullptr` if not
   *         available or wrong type
   * @see `get()`
   * 
   * This method is fully equivalent to `get()`, with the difference that the
   * returned object may be modified.
   */
  template <typename Obj = TObject>
  Obj* use(std::string const& name) const;
  
  /**
   * @brief Fetches an object with the specified name to be modified.
   * @tparam Obj (default: `TObject`) type of the object to fetch
   * @param name unprocessed name and path of the object to fetch
   * @return the requested object
   * @throw cet::exception (category: `"PlotSandbox"`) if no object with `name`
   *        exists in the box
   * @see `get()`, `use()`
   * 
   * This method is equivalent to `use()`, with the difference that the
   * returned object must exist.
   */
  template <typename Obj = TObject>
  Obj& demand(std::string const& name) const;
  
  /**
   * @brief Fetches the directory with the specified name from the sandbox.
   * @tparam DirObj (default: `TDirectory`) type of ROOT directory object to get
   * @param path path of the directory to fetch within this sandbox
   * @return a constant pointer to the requested directory, or `nullptr` if not
   *         available or wrong type
   * 
   * The fetched object is converted to the desired type via `dynamic_cast`.
   * If conversion fails, a null pointer is returned.
   */
  template <typename DirObj = TDirectory>
  DirObj* getDirectory(std::string const& path) const;
  
  /**
   * @brief Creates a new ROOT object with the specified name and title.
   * @tparam Obj type of ROOT object to be created
   * @tparam Args types of the argumenst to be forwarded to the constructor
   * @param name unprocessed name of the new object
   * @param title unprocessed title of the new object
   * @param args additional arguments forwarded to the constructor
   * @return a pointer to the newly created object
   * 
   * The name and title are processed with `processName()` and `processTitle()`
   * method respectively, before the object is created.
   * 
   */
  template <typename Obj, typename... Args>
  Obj* make(std::string const& name, std::string const& title, Args&&... args);
  
  /// @}
  // --- END -- ROOT object management -----------------------------------------
  
  
  // --- BEGIN -- Contained sandboxes ------------------------------------------
  /// @name Contained sandboxes
  /// @{
  
  /**
   * @brief Creates a new sandbox contained in this one.
   * @tparam SandboxType (default: `PlotSandbox`) type of sandbox to be created
   * @tparam Args types of the additional arguments for the sandbox constructor
   * @param name name of the new contained sand box
   * @param desc description of the new contained sand box
   * @param args additional arguments for the sandbox constructor
   * @return a reference to the created sandbox
   * @throw cet::exception (category: `"PlotSandbox"`) if a sandbox with this
   *        name already exists.
   * 
   * The arguments of this method are equivalent to the ones of the constructor.
   * 
   * The new sand box parent is set to point to this sand box.
   */
  template <typename SandboxType = PlotSandbox, typename... Args>
  SandboxType& addSubSandbox(
    std::string const& name, std::string const& desc,
    Args&&... args
    );
  
  /// Returns the number of contained sand boxes.
  std::size_t nSubSandboxes() const { return fData.subBoxes.size(); }
  
  // @{
  /**
   * @brief Returns the first contained sandbox with the specified name.
   * @param name unprocessed name of the box to be retrieved
   * @return the requested contained sandbox, or `nullptr` if not found
   * @see `demandSandbox()`
   */
  PlotSandbox const* findSandbox(std::string const& name) const;
  PlotSandbox* findSandbox(std::string const& name);
  // @}
  
  // @{
  /**
   * @brief Returns the first contained sandbox with the specified name.
   * @param name unprocessed name of the box to be retrieved
   * @return the requested contained sandbox
   * @throw cet::exception (category: `"PlotSandbox"`) if no a sandbox with this
   *        `name` exists
   * @see `findSandbox()`
   */
  PlotSandbox const& demandSandbox(std::string const& name) const;
  PlotSandbox& demandSandbox(std::string const& name);
  // @}
  
  // @{
  /// Returns an object proper to iterate through all contained sand boxes.
  decltype(auto) subSandboxes() const
    { return details::map_dereferenced_values(fData.subBoxes); }
  decltype(auto) subSandboxes()
    { return details::map_dereferenced_values(fData.subBoxes); }
  // @}
  
  /// @}
  // --- END -- Contained sandboxes --------------------------------------------
  
  
  /// Dumps the hierarchy of sandboxes into the specified stream.
  template <typename Stream>
  void dump(Stream&& out, std::string indent, std::string firstIndent) const;
  
  /// Dumps the hierarchy of sandboxes into the specified stream.
  template <typename Stream>
  void dump(Stream&& out, std::string indent = "") const
    { dump(std::forward<Stream>(out), indent, indent); }
  
    protected:
  
  /**
   * @brief Constructor: specifies all sandbox characteristics.
   * @param parent the sandbox used as parent
   * @param name the name of the sandbox
   * @param desc description of the sandbox
   * 
   * Compared to the public constructor, this one picks the parent directory
   * to be the output directory of the `parent` sandbox, and it registers
   * that `parent` box as parent.
   * Note that this constructor does *not* update the list of subboxes of the
   * parent.
   */
  PlotSandbox(PlotSandbox const& parent, std::string name, std::string desc);
  
  
  /// Sets the parent of this box.
  virtual void setParent(PlotSandbox const* parent) { fData.parent = parent; }
  
  
  /// Applies title processing only at the title part of the string.
  std::string processPlotTitle(std::string const& title) const;


  /// Returns a processed version of the name of this sandbox.
  /// 
  /// This may be used when integrating the sandbox name into the processed
  /// object names.
  virtual std::string processedSandboxName() const;

  /// Returns a processed version of the description of this sandbox.
  /// 
  /// This may be used when integrating the sandbox description into the
  /// processed object descriptions.
  virtual std::string processedSandboxDesc() const;

  /// Dumps the content of this box (nosubboxes) into `out` stream.
  template <typename Stream>
  void dumpContent
    (Stream&& out, std::string indent, std::string firstIndent) const;
  
  
  /// Retrieves or, if not present, creates a ROOT subdirectory in the sandbox.
  /// Returns a pair with the directory and the name part of `path`.
  /// The directory part may be empty.
  static std::pair<std::string, std::string> splitPath
    (std::string const& path, char sep = '/');
  
  /// Merges the pieces of path that are not empty into a path.
  /// One separator at the end of each piece is ignored.
  static std::string joinPath
    (std::initializer_list<std::string> pathElements, char sep = '/');
  
}; // icarus::trigger::PlotSandbox


//------------------------------------------------------------------------------
//---  Standard library support
//------------------------------------------------------------------------------
namespace std {
  
  template <>
  struct hash<icarus::trigger::PlotSandbox> {
    auto operator() (icarus::trigger::PlotSandbox const& key) const
      { return std::hash<std::string>()(key.ID()); }
  }; // hash<PlotSandbox>
  
} // namespace std


//------------------------------------------------------------------------------
//---  template implementation
//------------------------------------------------------------------------------
template <typename RootDir /* = TDirectoryFile */>
auto icarus::trigger::PlotSandbox::TFileDirectoryHelper::create(
  art::TFileDirectory parentDir,
  std::string const& subdir, std::string const& dirTitle /* = "" */
  )
  -> TFileDirectoryHelper
{
  // NOTE: we only support a direct subdirectory of parentDir.
  // NOTE: if the directory already exists the results are undefined;
  //       we can't figure out if the direct
  
  // first we create the directory directly via ROOT,
  // but starting from the directory stored in `parentDir`
  TDirectory* pROOTdir
    = parentDir.make<RootDir>(subdir.c_str(), dirTitle.c_str());
  
  // then we create a `art::TFileDirectory` for the same directory;
  // this is the only way we can create a new `art::TFileDirectory`,
  // and it actually does not create the ROOT directory because it's lazy,
  // and it will create it only when it is needed, if not present yet.
  return { parentDir.mkdir(subdir, dirTitle), pROOTdir };
  
} // icarus::trigger::PlotSandbox::TFileDirectoryHelper::create()


//------------------------------------------------------------------------------
template <typename Obj /* = TObject */>
Obj const* icarus::trigger::PlotSandbox::get(std::string const& name) const
  { return use<std::add_const_t<Obj>>(name); }


//------------------------------------------------------------------------------
template <typename Obj /* = TObject */>
Obj* icarus::trigger::PlotSandbox::use(std::string const& name) const {
  
  auto [ objDir, objName ] = splitPath(name);
  
  TDirectory* dir = getDirectory(objDir);
  if (!dir) return nullptr;
  
  std::string const processedName = processName(objName);
  return dir->Get<Obj>(processedName.c_str());
  
} // icarus::trigger::PlotSandbox::use()


//------------------------------------------------------------------------------
template <typename Obj /* = TObject */>
Obj& icarus::trigger::PlotSandbox::demand(std::string const& name) const {
  
  auto* obj = use<Obj>(name);
  if (obj) return *obj;
  cet::exception e { "PlotSandbox" };
  e << "PlotSandbox::demand(): object '" << name
    << "' not available in the sandbox '" << ID() << "'"
    << "\nBox content: ";
  dumpContent(e, "", ""); // no indent
  
  throw e << "\n";
} // icarus::trigger::PlotSandbox::demand()


//------------------------------------------------------------------------------
template <typename DirObj /* = TDirectory */>
DirObj* icarus::trigger::PlotSandbox::getDirectory
  (std::string const& path) const
{
  TDirectory* pBaseDir = fData.outputDir.fROOTdir;
  return dynamic_cast<DirObj*>
    (path.empty()? pBaseDir: pBaseDir->GetDirectory(path.c_str()));
} // icarus::trigger::PlotSandbox::getDirectory()


//------------------------------------------------------------------------------
template <typename Obj, typename... Args>
Obj* icarus::trigger::PlotSandbox::make
  (std::string const& name, std::string const& title, Args&&... args)
{
  auto [ objDir, objName ] = splitPath(name);
  
  std::string const processedName = processName(objName);
  std::string const processedTitle = processPlotTitle(title);
  
  art::TFileDirectory destDir
    = objDir.empty()? fData.outputDir.fDir: fData.outputDir.fDir.mkdir(objDir);
  
  // see Redmine issue #23075
  return destDir.makeAndRegister<Obj>(
    processedName.c_str(), processedTitle.c_str(),
    processedName.c_str(), processedTitle.c_str(), std::forward<Args>(args)...
    );
} // icarus::trigger::PlotSandbox::make()


//------------------------------------------------------------------------------
template
  <typename SandboxType /* = icarus::trigger::PlotSandbox */, typename... Args>
SandboxType& icarus::trigger::PlotSandbox::addSubSandbox
  (std::string const& name, std::string const& desc, Args&&... args)
{
  // we can't use make_unique() because the constructor it needs is protected:
  auto [ it, bInserted ] = fData.subBoxes.try_emplace
    (name, new SandboxType(*this, name, desc, std::forward<Args>(args)...));
  if (!bInserted) {
    throw cet::exception("PlotSandbox")
      << "PlotSandbox::addSubSandbox(): a subbox with name '" << name
      << "' already exists in  box '" << ID() << "'.\n";
  }
  return *(it->second); // it iterator to the inserted element
} // icarus::trigger::PlotSandbox::addSubSandbox()


//------------------------------------------------------------------------------
template <typename SandboxType>
auto* icarus::trigger::PlotSandbox::findSandbox
  (SandboxType& sandbox, std::string const& name)
{
  auto const it = sandbox.fData.subBoxes.find(name);
  return (it == sandbox.fData.subBoxes.end())? nullptr: it->second.get();
}


//------------------------------------------------------------------------------
template <typename SandboxType>
auto& icarus::trigger::PlotSandbox::demandSandbox
  (SandboxType& sandbox, std::string const& name)
{
  auto* box = findSandbox(sandbox, name);
  if (box) return *box;
  
  cet::exception e { "PlotSandbox" };
  e << "PlotSandbox::demandSandbox(): box '" << name
    << "' not available in the sandbox '" << sandbox.ID() << "'";
  if (sandbox.nSubSandboxes()) {
    e << "\n" << "Available nested boxes (" << sandbox.nSubSandboxes() << "):";
    for (auto const& subbox: sandbox.subSandboxes())
      e << "\n * '" << subbox.ID() << "'";
  } // if has subboxes
  else {
    e << "  (no contained box!)";
  }
  throw e << "\n";

} // icarus::trigger::PlotSandbox::demandSandbox(SandboxType)



//------------------------------------------------------------------------------
template <typename Stream>
void icarus::trigger::PlotSandbox::dump
  (Stream&& out, std::string indent, std::string firstIndent) const
{
  out << firstIndent;
  if (hasName()) out << "Box '" << name() << "'";
  else           out << "Unnamed box";
  if (hasDescription()) out << " (\"" << description() << "\")";
  out << " [ID=" << ID() << "] with ";
  dumpContent(std::forward<Stream>(out), indent, "");
  
  if (nSubSandboxes()) {
    out << "\n" << indent << "Nested boxes (" << nSubSandboxes() << "):";
    for (auto const& subbox: subSandboxes()) {
      out << "\n";
      subbox.dump(std::forward<Stream>(out), indent + "  ");
    }
  } // if has subboxes
} // icarus::trigger::PlotSandbox::dump()


//------------------------------------------------------------------------------
template <typename Stream>
void icarus::trigger::PlotSandbox::dumpContent
  (Stream&& out, std::string indent, std::string firstIndent) const
{
  out << firstIndent;
  
  TDirectory const* pDir = fData.outputDir.fROOTdir;
  if (!pDir) {
    out << "no content available";
    return;
  }
  
  TList const* objects = pDir->GetList();
  TList const* keys = pDir->GetListOfKeys();
  if (objects && !objects->IsEmpty()) {
    out << objects->GetSize() << " direct entries:";
    for (TObject const* obj: *objects) {
      out << "\n" << indent << "  '" << obj->GetName() << "'  ["
        << obj->IsA()->GetName() << "]";
    } // for objects
  }
  else out << "no direct entries;";
  if (keys) {
    for (TObject const* keyObj: *keys) {
      auto key = dynamic_cast<TKey const*>(keyObj);
      if (!key) continue;
      if (objects->Contains(key->GetName())) continue; // already in objects
      out << "\n" << indent
        << "[KEY]  '" << key->GetName() << "'  ["
        << key->GetClassName() << "]"
        ;
    } // for
  } // if has keys
  
} // icarus::trigger::PlotSandbox::dumpContent()


//------------------------------------------------------------------------------


#endif // ICARUSCODE_PMT_TRIGGER_UTILITIES_PLOTSANDBOX_H
