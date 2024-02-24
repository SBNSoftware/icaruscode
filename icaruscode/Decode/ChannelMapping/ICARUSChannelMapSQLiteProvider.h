/**
 * @file   icaruscode/Decode/ChannelMapping/ICARUSChannelMapSQLiteProvider.h
 * @author T. Usher (factorised by G. Petrillo, petrillo@slac.stanford.edu)
 * @see    icaruscode/Decode/ChannelMapping/ICARUSChannelMapSQLiteProvider.cxx
 */

#ifndef ICARUSCODE_DECODE_CHANNELMAPPING_ICARUSCHANNELMAPSQLITEPROVIDER_H
#define ICARUSCODE_DECODE_CHANNELMAPPING_ICARUSCHANNELMAPSQLITEPROVIDER_H

// ICARUS libraries
#include "icaruscode/Decode/ChannelMapping/ICARUSChannelMapProviderBase.h"
#include "icaruscode/Decode/ChannelMapping/ChannelMapSQLite.h"


// -----------------------------------------------------------------------------
namespace icarusDB { class ICARUSChannelMapSQLiteProvider; }
/**
 * @brief Interface to the SQLite ICARUS channel mapping database.
 * 
 * The database is normally distributed as a file in `icarus_data`.
 * 
 * 
 * The implementation is fully delegated to
 * `icarusDB::ICARUSChannelMapProviderBase`.
 * 
 * 
 * Configuration parameters
 * =========================
 * 
 * See `icarusDB::ICARUSChannelMapProviderBase`, except for:
 * 
 * * `ChannelMappingTool` (algorithm configuration): see
 *     `icarusDB::ChannelMapSQLite` configuration.
 * 
 */
class icarusDB::ICARUSChannelMapSQLiteProvider
  : public icarusDB::ICARUSChannelMapProviderBase<icarusDB::ChannelMapSQLite>
{
  using Base_t = ICARUSChannelMapProviderBase<icarusDB::ChannelMapSQLite>;
  using Base_t::Base_t; 
};


// -----------------------------------------------------------------------------

#endif // ICARUSCODE_DECODE_CHANNELMAPPING_ICARUSCHANNELMAPSQLITEPROVIDER_H
