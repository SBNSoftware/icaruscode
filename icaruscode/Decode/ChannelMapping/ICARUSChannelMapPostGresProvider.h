/**
 * @file   icaruscode/Decode/ChannelMapping/ICARUSChannelMapPostGresProvider.h
 * @author T. Usher (factorised by G. Petrillo, petrillo@slac.stanford.edu)
 * @see    icaruscode/Decode/ChannelMapping/ICARUSChannelMapPostGresProvider.cxx
 */

#ifndef ICARUSCODE_DECODE_CHANNELMAPPING_ICARUSCHANNELMAPPOSTGRESPROVIDER_H
#define ICARUSCODE_DECODE_CHANNELMAPPING_ICARUSCHANNELMAPPOSTGRESPROVIDER_H

// ICARUS libraries
#include "icaruscode/Decode/ChannelMapping/ICARUSChannelMapProviderBase.h"
#include "icaruscode/Decode/ChannelMapping/ChannelMapPostGres.h"


// -----------------------------------------------------------------------------
namespace icarusDB { class ICARUSChannelMapPostGresProvider; }
/**
 * @brief Interface to the PostgreSQL ICARUS channel mapping database.
 * 
 * The database is normally deployed on the network.
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
 *     `icarusDB::ChannelMapPostGres` configuration.
 * 
 */
class icarusDB::ICARUSChannelMapPostGresProvider
  : public icarusDB::ICARUSChannelMapProviderBase<icarusDB::ChannelMapPostGres>
{
  using Base_t = ICARUSChannelMapProviderBase<icarusDB::ChannelMapPostGres>;
  using Base_t::Base_t; 
};


// -----------------------------------------------------------------------------

#endif // ICARUSCODE_DECODE_CHANNELMAPPING_ICARUSCHANNELMAPPOSTGRESPROVIDER_H
