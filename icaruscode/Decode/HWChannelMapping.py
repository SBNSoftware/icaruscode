#-------------------------------------------------------------------
#
# This script is to generate lot of functions which will help      /
# on mapping between different important piece of channel mapping  /
#
# This is a draft script for testing
# Author: Biswaranjan.Behera@colostate.edu @ 13 May 2020
#
#-------------------------------------------------------------------

from DataLoader3 import DataLoader, DataQuery
import sys, os


queryUrl = "https://dbdata0vm.fnal.gov:9443/QE/hw/app/SQ/query" 
print("\n Hello! Welcome to DataQuery form ICARUS Hardware Database: \n")


#Data query
dataQuery = DataQuery(queryUrl)

# Mapping between ChannelID ----> Physical location
def physical_loc (ChannelID):
    
    geo_loc = dataQuery.query("icarus_hardware_dev", "daq_channels", 'readout_board_id',
                              'channel_id:eq:'+str(ChannelID))
    
    Flange_ID = dataQuery.query("icarus_hardware_dev", "readout_boards",
                                'flange_id',
                                'readout_board_id:eq:'+str(geo_loc)[2:-6])
    
    geo_loc = dataQuery.query("icarus_hardware_dev", "flanges", 'flange_pos_at_chimney',
                              'flange_id:eq:'+str(Flange_ID)[2:-6])
    
    return str(geo_loc)[2:-6]
    
    
# for example..    
print (physical_loc(45295))


# Mapping between PMTChannelID (Larsoft) From Physical location of PMT ID
def PMTChannelIDFromPhysicalPMTPositionID (pmtid):
    geo_loc = dataQuery.query("icarus_hardware_dev", "pmt_placements",
                              'channel_id',
                              'pmt_id:eq:'+str(pmtid))
    
    return str(geo_loc)[2:-6]
    
# Mapping between PMTChannelID (Larsoft) From Digitizer Label and Digitizer Channel
def PMTChannelIDFromDigitizerLabel (digitizerlabel, digitizerchannel):
    geo_loc = dataQuery.query("icarus_hardware_dev", "pmt_placements",
                              'channel_id',
                              ['digitizer_label:eq:'+str(digitizerlabel), 'digitizer_ch_number:eq:'+str(digitizerchannel)])
    
    
    return str(geo_loc)[2:-6]

print (PMTChannelIDFromPhysicalPMTPositionID (333))
print (PMTChannelIDFromDigitizerLabel ("EE-BOT-B", 10))
