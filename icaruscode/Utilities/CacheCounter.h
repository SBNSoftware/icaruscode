/**
 * @file   icaruscode/Utilities/CacheCounter.h
 * @brief  Simple utility to track cache changes.
 * @author Gianluca Petrillo (petrillo@slac.standord.edu)
 * @date   March 11, 2024
 * 
 * This library is header only.
 */

#ifndef ICARUSCODE_UTILITIES_CACHECOUNTER_H
#define ICARUSCODE_UTILITIES_CACHECOUNTER_H

// C/C++ standard libraries
#include <map>
#include <string>
#include <initializer_list>
#include <atomic>
#include <stdexcept> // std::logic_error


// -----------------------------------------------------------------------------
namespace util {
  class CacheCounter;
  class CacheGuard;
}

// -----------------------------------------------------------------------------
/**
 * @brief Counter of cache versions.
 * @see `util::CacheGuard`
 * 
 * This simple class contains an index of cache IDs.
 * A cache ID is a number that represents the status of a cache.
 * 
 * Multiple caches are supported by a single object instance.
 * Each cache is represented by a "tag", or name. The cache with empty name is
 * special in that it represent the whole set of caches: when updating it, all
 * the other cache IDs are also updated.
 * 
 * Note that the special ID value `NoCache` represents the idea that the cache
 * is never up to date and always needs to be refreshed (functionally equivalent
 * to not having any cache).
 * 
 * See `util::CacheGuard` for a usage example.
 */
class util::CacheCounter {
  
    public:
  
  using CacheID_t = unsigned long; ///< Type of cache identifier.
  using CacheTag_t = std::string; ///< Type of cache name.
  
  /// ID specifying that cache is not present yet.
  static constexpr CacheID_t NoCache = 0;
  
  /// Constructor: prepares the ID of the default cache tag.
  CacheCounter() { fCacheID.emplace("", NoCache); }
  
  /// Constructor: prepares the ID of the caches with the specified tags.
  CacheCounter(std::initializer_list<CacheTag_t> tags)
    : CacheCounter{}
    { addCacheTags(std::move(tags)); }
  
  /// Adds tags. The ID of existing ones are not changed.
  void addCacheTags(std::initializer_list<CacheTag_t> tags)
    { for (CacheTag_t const& tag: tags) fCacheID.emplace(tag, NoCache); }
  
  /// Returns the current ID of the cache with the specified tag.
  CacheID_t cacheID(CacheTag_t const& tag) const { return fCacheID.at(tag); }
  
  /// Returns whether the specified tag is tracked.
  bool hasTag(CacheTag_t const& tag) const { return fCacheID.count(tag) > 0; }
  
  /// Returns whether the cache `tag` was updated since it has `sinceID` ID.
  bool wasCacheUpdatedSince(CacheTag_t const& tag, CacheID_t sinceID) const
    { return cacheID(tag) > sinceID; }
  
  /// Returns whether the default cache was updated since it has `sinceID` ID.
  bool wasCacheUpdatedSince(CacheID_t sinceID) const
    { return wasCacheUpdatedSince({}, sinceID); }
  
  
    protected:
  
  /// ID of the current caches.
  mutable std::map<CacheTag_t, std::atomic<CacheID_t>> fCacheID;
  
  /// Increments the counter of the specified cache tag.
  /// The default (empty) tag increments all.
  CacheID_t updateCacheID(CacheTag_t const& cacheTag = "") const
    {
      if (cacheTag.empty()) {
        for (auto& ID: fCacheID) if (!ID.first.empty()) ++(ID.second);
      }
      else ++fCacheID.at(""); // always update the default cache
      return ++fCacheID.at(cacheTag); // return only the ID from the argument
    }
  
}; // util::CacheCounter


// -----------------------------------------------------------------------------
/**
 * @brief Tracks the status of a cache.
 * @see `util::CacheCounter`
 * 
 * This class tracks the status of a cache, represented by an integer (ID).
 * The cache is supposed to update (increment) that ID value each time it
 * updates the cache content itself. This object tracks that ID and reports
 * whether it has changed since the last time.
 * 
 * There are two elements of the interface:
 *  * `update()` unconditionally acquires the current ID of the cache; it also
 *      returns whether it's a new ID (like `wasUpdated()`).
 *  * `wasUpdated()` compares the current cache ID to the one acquired by the
 *      last `update()` call, and returns whether the current cache is newer.
 * 
 * While the cache ID can be any number, and the interface allows to directly
 * provide the cache ID, nonetheless this class has a special relationship with
 * `util::CacheCounter`: it can point to a specific tag name of any cache,
 * or it can point to a specific cache and tag, with a simpler interface.
 * If a `CacheGuard` object needs to track a specific cache or cache tag, that
 * must be specified at  construction time.
 * 
 * Note that the object always starts with no update.
 * 
 * Also note that a cache ID of `util::CacheCounter::NoCache` is specially
 * treated and a cache with that ID will never be considered up to date.
 * 
 * 
 * Example
 * --------
 * 
 * A database access class performs standard queries reading time-dependent
 * data from a database. It caches the result of the query, and it updates
 * the cache when a time is requested that is not the one present in the cache.
 * This class also derives from `util::CacheCounter` to offer the cache updating
 * track capability.
 * 
 * An analysis class uses the data from the database to produce results.
 * It stores the small amount of data derived by the database that is needed
 * for the computation, and it needs to refetch or recompute that data when
 * its local cache is outdated. To do this it has a `util::CacheGuard` data
 * member bound to the database access class above.
 * 
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * class DatabaseQuery: public util::CacheCounter {
 *   
 *   bool fetchData()
 *     {
 *       // ...
 *       return true;
 *     }
 *   
 *   void updateCache() { if (fetchData()) updateCacheID(); }
 *   
 *     public:
 *   
 *   using Data_t = ... ;
 *   
 *   Data_t query();
 *   
 * };
 * 
 * DatabaseQuery DBquery;
 * 
 * 
 * class DataUser {
 *   
 *   util::CacheGuard fDatabaseQueryCacheGuard;
 *   DatabaseQuery::Data_t fCachedQueryResult;
 *   
 *   DatabaseQuery::Data_t const& queryResults()
 *     {
 *       if (fDatabaseQueryCacheGuard.update())
 *         fCachedQueryResult = DBquery.query();
 *       return fCachedQueryResult;
 *     }
 *   
 *   
 *     public:
 *   
 *   DataUser(): fDatabaseQueryCache{ DBquery } {}
 *   
 *   void compute()
 *     {
 *       DatabaseQuery::Data_t const& data = queryResults();
 *       // ...
 *     }
 *   
 * }; // DataUser
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * 
 * 
 */
class util::CacheGuard {
  
  /// Pointer to the cache to monitor by default.
  CacheCounter const* fDefaultCache = nullptr;
  
  /// Tag of the cache to monitor.
  CacheCounter::CacheTag_t fCacheTag;
  
  /// The ID of the cache at the last check.
  CacheCounter::CacheID_t fLastCacheID = CacheCounter::NoCache;
  
  /// Returns the default cache object, throwing `std::logic_error` if not set.
  CacheCounter const& defaultCache() const
    {
      if (!fDefaultCache)
        throw std::logic_error{ "util::CacheGuard default cache was not set." };
      return *fDefaultCache;
    }
  
  /// Changes the default cache.
  void chooseCache(CacheCounter const* cache)
    { fDefaultCache = cache; fLastCacheID = CacheCounter::NoCache; }
  
    public:
  
  /// Constructor: prepares monitoring of cache `tag`.
  CacheGuard(CacheCounter::CacheTag_t tag = ""): fCacheTag{ std::move(tag) } {}
  
  /// Constructor: prepares monitoring of cache `tag` of the specified object.
  CacheGuard(CacheCounter const& cache, CacheCounter::CacheTag_t tag = "")
    : fDefaultCache{ &cache }, fCacheTag{ std::move(tag) } {}
  
  /// Changes the default cache.
  void setCache(CacheCounter const& cache) { chooseCache(&cache); }
  
  /// Removes the default cache.
  void unsetCache() { chooseCache(nullptr); }
  
  
  /// Returns the ID of the cache at the last update.
  CacheCounter::CacheID_t lastUpdateID() const { return fLastCacheID; }
  
  
  /// @{
  
  /// Returns whether the cache, which has `cacheID` ID, was updated since the
  /// last `update()` call.
  bool wasUpdated(CacheCounter::CacheID_t cacheID) const
    { return (cacheID != CacheCounter::NoCache) && (cacheID > fLastCacheID); }
  
  /// Returns whether the `cache` was updated since the last `update()` call.
  bool wasUpdated(CacheCounter const& cache) const
    { return wasUpdated(cache.cacheID(fCacheTag)); }
  
  /// Returns whether the default cache was updated since the last `update()`
  /// call.
  bool wasUpdated() const { return wasUpdated(defaultCache()); }
  
  /// @}
  
  
  /// @{
  
  /// Returns if the cache was outdated, and at the same time updates its ID.
  bool update(CacheCounter::CacheID_t cacheID)
    {
      if (cacheID == CacheCounter::NoCache) return true; // no tracking
      if (fLastCacheID == cacheID) return false;
      fLastCacheID = cacheID;
      return true;
    }
  
  /// Returns if the cache was outdated, and at the same time updates its ID.
  bool update(CacheCounter const& cache)
    { return update(cache.cacheID(fCacheTag)); }
  
  /// Returns if the default cache was outdated, and at the same time updates
  /// its ID.
  bool update() { return update(defaultCache()); }
  
  /// @}
  
}; // util::CacheGuard


// -----------------------------------------------------------------------------

#endif // ICARUSCODE_UTILITIES_CACHECOUNTER_H
