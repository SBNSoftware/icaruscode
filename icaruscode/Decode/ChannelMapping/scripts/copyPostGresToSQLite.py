#! /usr/bin/env python
#---------------------------------------------------------------
#
# Name: copyPostGresToSQLite.py
#
# Purpose: Accessing the channel mapping database
#
# Created: 23-Oct-2020 T. Usher
#
#---------------------------------------------------------------

from __future__ import print_function
import sys, os, datetime
import psycopg2
import sqlite3
from   DataLoader3 import DataLoader, DataQuery

# Set up the columns for the three databases
readoutBoardsColumns = ['readout_board_id integer',
                        'flange_id integer',
                        'chimney_number integer',
                        'tpc_id text',
                        'create_time text',
                        'create_user text',
                        'update_time text',
                        'update_user text',
                        'fragement_id text']

daqChannelsColumns = ['channel_id integer',
                      'wire_number integer',
                      'readout_board_id integer',
                      'chimney_number integer',
                      'readout_board_slot integer',
                      'channel_number integer',
                      'create_time text',
                      'create_user text',
                      'update_time text',
                      'update_user text',
                      'plane text',
                      'cable_label_number text',
                      'channel_type text']

flangesColumns = ['flange_id integer',
                  'tpc_id text',
                  'chimney_number integer',
                  'flange_pos_at_chimney text',
                  'create_time text',
                  'update_user text',
                  'update_time text',
                  'create_user text',
                  'power_supply_ip_address text',
                  'power_supply_id integer']

pmtPlacementColumns = ['pmt_id integer',
                       'pmt_sn text',
                       'sector_label text',
                       'ch_number integer',
                       'pmt_position_code integer',
                       'hv_cable_label text',
                       'signal_cable_label text',
                       'light_fiber_label text',
                       'digitizer_label text',
                       'digitizer_ch_number integer',
                       'hv_supply_label text',
                       'hv_supply_ch_number integer',
                       'create_time text',
                       'update_user text',
                       'update_time text',
                       'create_user text',
                       'pmt_in_tpc_plane text',
                       'channel_id integer',
                       'fragment_id integer']

crtFEBChannelsColumns = ['feb_id text',
                         'feb_channel integer',
                         'pedestal float',
                         'threshold_adjust integer',
                         'bias integer',
                         'hg integer',
                         'create_time text',
                         'update_user text',
                         'update_time text',
                         'create_user text',
                         'channel_id integer',
                         'feb_index integer',
                         'mac_address integer']

topcrtFEBColumns = ['feb_barcode text',
                    'serialnum integer',
                    'mac_add8b text',
                    'mac_add integer',
                    'voltage float',
                    'ch0 integer',
                    'ch1 integer',
                    'ch2 integer',
                    'ch3 integer',
                    'ch4 integer',
                    'ch5 integer',
                    'ch6 integer',
                    'ch7 integer',
                    'ch8 integer',
                    'ch9 integer',
                    'ch10 integer',
                    'ch11 integer',
                    'ch12 integer',
                    'ch13 integer',
                    'ch14 integer',
                    'ch15 integer',
                    'ch16 integer',
                    'ch17 integer',
                    'ch18 integer',
                    'ch19 integer',
                    'ch20 integer',
                    'ch21 integer',
                    'ch22 integer',
                    'ch23 integer',
                    'ch24 integer',
                    'ch25 integer',
                    'ch26 integer',
                    'ch27 integer',
                    'ch28 integer',
                    'ch29 integer',
                    'ch30 integer',
                    'ch31 integer',
                    'create_time text',
                    'update_user text',
                    'update_time text',
                    'create_user text',
                    'feb_index integer']

# Define the function to create and fill each table
def copyTable(postGres, dbCurs, dbName, table, columns):
    createString = "CREATE TABLE " + table + " ("

    for columnIdx in range(len(columns)):
        if columnIdx > 0:
            createString += ", "
        createString += columns[columnIdx]

    createString += ")"

    print(createString)

    dbCurs.execute(createString)

    query = postGres.query(database=dbName, table=table, columns="*")

    for rowIdx in range(len(query)):
        rowList = query[rowIdx].split(',')
        if len(rowList) != len(columns):
            print("Length mismatch! Will skip this row (",rowIdx,")")
            continue
        insertString = "INSERT INTO " + table + " VALUES ("
        for idx in range(len(rowList)):
            if idx > 0:
                insertString += ", "
            if "text" in columns[idx]:
                fieldEntry = rowList[idx]
                if table == "daq_channels" and "plane" in columns[idx]:
                    fieldEntry = rowList[idx].upper()
                insertString += "\'" + fieldEntry + "\'"
            else:
                if "none" in rowList[idx]:
                    #print("Found none in column data, field:",columns[idx],", value: ",rowList[idx])
                    insertString += '0'
                else:
                    insertString += "\'" +rowList[idx]+ "\'"
        insertString += ")"
        if table == "daq_channels":
            print("idx:",rowIdx,"  -->",insertString)
        dbCurs.execute(insertString)

# Ok get down to extracting the info and copying to SQLite
queryurl = "https://dbdata0vm.fnal.gov:9443/QE/hw/app/SQ/query"

dataQuery = DataQuery(queryurl)

print(dataQuery)

dbName             = "icarus_hardware_prd"
readoutBoardsTable = "readout_boards"
daqChannelsTable   = "daq_channels"
flangesTable       = "flanges"
pmtPlacementTable  = "pmt_placements"
crtfebchannelsTable= "feb_channels"
topcrtfebTable     = "crtfeb"
##################################################################################

# Start by creating the new database and setting up the first table
sqliteDB = sqlite3.connect("ChannelMapICARUS.db")
dbCurs   = sqliteDB.cursor()

###################################################################################

copyTable(dataQuery, dbCurs, dbName, readoutBoardsTable, readoutBoardsColumns)

copyTable(dataQuery, dbCurs, dbName, daqChannelsTable, daqChannelsColumns)

copyTable(dataQuery, dbCurs, dbName, flangesTable, flangesColumns)

copyTable(dataQuery, dbCurs, dbName, pmtPlacementTable, pmtPlacementColumns)

copyTable(dataQuery, dbCurs, dbName, crtfebchannelsTable, crtFEBChannelsColumns)

copyTable(dataQuery, dbCurs, dbName, topcrtfebTable, topcrtFEBColumns)

###################################################################################

sqliteDB.commit()

sqliteDB.close()



