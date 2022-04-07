/**
 * @file   icaruscode/Decode/DecoderTools/details/KeyValuesData.cxx
 * @brief  Simple parsed data format.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   May 9, 2021
 * @see    icaruscode/Decode/DecoderTools/details/KeyValuesData.h
 */

// ICARUS libraries
#include "icaruscode/Decode/DecoderTools/details/KeyValuesData.h"

// C++ standard libraries
#include <ostream>
#include <utility> // std::move(), std::as_const()


// -----------------------------------------------------------------------------
// ---  icarus::KeyValuesData
// -----------------------------------------------------------------------------
auto icarus::KeyValuesData::makeItem(std::string key) -> Item& {
  if (hasItem(key)) throw DuplicateKey{ std::move(key) };
  return makeItemImpl(std::move(key));
} // icarus::KeyValuesData<>::makeItem()


// -----------------------------------------------------------------------------
auto icarus::KeyValuesData::makeOrFetchItem(std::string const& key) -> Item& {
  
  Item* item = findItem(key);
  return item? *item: makeItemImpl(key);
  
} // icarus::KeyValuesData<>::makeOrFetchItem()


// -----------------------------------------------------------------------------
auto icarus::KeyValuesData::findItem
  (std::string const& key) const noexcept -> Item const*
{
  for (auto const& item: fItems) if (key == item.key) return &item;
  return nullptr;
} // icarus::KeyValuesData<>::findItem() const


// -----------------------------------------------------------------------------
auto icarus::KeyValuesData::findItem
  (std::string const& key) noexcept -> Item*
{
  // no violations here: this is a non-const method, with the right to modify
  // object data; and this avoids code duplication.
  return const_cast<Item*>(std::as_const(*this).findItem(key));
} // icarus::KeyValuesData<>::findItem()


// -----------------------------------------------------------------------------
auto icarus::KeyValuesData::getItem(std::string const& key) const
  -> Item const&
{
  if (auto item = findItem(key); item) return *item;
  throw ItemNotFound(key);
} // icarus::KeyValuesData<>::getItem()


// -----------------------------------------------------------------------------
bool icarus::KeyValuesData::hasItem
  (std::string const& key) const noexcept
{
  return findItem(key);
} // icarus::KeyValuesData<>::hasItem()


// -----------------------------------------------------------------------------
bool icarus::KeyValuesData::empty() const noexcept
  { return fItems.empty(); }


// -----------------------------------------------------------------------------
std::size_t icarus::KeyValuesData::size() const noexcept
  { return fItems.size(); }


// -----------------------------------------------------------------------------
decltype(auto) icarus::KeyValuesData::items() const noexcept
  { return fItems; }


// -----------------------------------------------------------------------------
auto icarus::KeyValuesData::makeItemImpl(std::string key) -> Item& {
  fItems.emplace_back(std::move(key));
  return fItems.back();
} // icarus::KeyValuesData<>::makeItemImpl()


// -----------------------------------------------------------------------------
std::ostream& icarus::operator<< (std::ostream& out, KeyValuesData const& data)
{
  out << data.size() << " items";
  if (!data.empty()) {
    out << ':';
    for (auto const& item: data.items()) {
      out << "\n  '" << item.key << "' (" << item.values.size() << ")";
      if (item.values.empty()) continue;
      out << ':';
      for (auto const& value: item.values) out << " '" << value << '\'';
    } // for
  } // if has data
  return out << "\n";
} // icarus::operator<< (KeyValuesData)


// -----------------------------------------------------------------------------
