#! /usr/bin/env python
#---------------------------------------------------------------
#
# Name: siov_extracter.py
#
# Purpose: Extract data from conditions database into sqlite databases.
#
# Created: 9-Dec-2019,  H. Greenlee
#
# Usage:
#
# $ iov_extracter.py [-d|--devel]
#
# Options:
#
# -d - Use development database.
#
# Notes:
#
# 1.  Source ~uboonepro/.sqlaccess/prod_access.sh to define environment variables
#     for database access.
#
#---------------------------------------------------------------

from __future__ import print_function
import sys, os, datetime
import psycopg2
import sqlite3

# Calibration database names.

calibs = ['crt_channelstatus_data',
          'detpedestals_data',
          'detpmtgains_data',
          'dqdx_xnorm_data',
          'dqdx_xshape_plane0_data',
          'dqdx_xshape_plane0_mc',
          'dqdx_xshape_plane1_data',
          'dqdx_xshape_plane1_mc',
          'dqdx_xshape_plane2_data',
          'dqdx_xshape_plane2_mc',
          'dqdx_yz_plane0_data',
          'dqdx_yz_plane0_mc',
          'dqdx_yz_plane1_data',
          'dqdx_yz_plane1_mc',
          'dqdx_yz_plane2_data',
          'dqdx_yz_plane2_mc',
          'electronicscalib_data',
          'elifetime_data',
          'lightyieldscale_data']
         
calibs_dev = ['crt_channelstatus_data',
              'detpedestals_data',
              'detpmtgains_data',
              'dqdx_xnorm_data',
              'dqdx_xshape_plane0_data',
              'dqdx_xshape_plane1_data',
              'dqdx_xshape_plane2_data',
              'dqdx_yz_plane0_data',
              'dqdx_yz_plane0_mc',
              'dqdx_yz_plane1_data',
              'dqdx_yz_plane1_mc',
              'dqdx_yz_plane2_data',
              'dqdx_yz_plane2_mc',
              'electronicscalib_data',
              'elifetime_data',
              'lightyieldscale_data']

# Table suffixes.

suffixes = ['tags', 'iovs', 'tag_iovs', 'data']

# Parse command line.

devel = 0
dev = ''
args = sys.argv[1:]
while len(args) > 0:
    if args[0] == '-d' or args[0] == '--devel' :
        devel = 1
        dev = '_dev'
        del args[0]
    else:
        print('Unknown option %s' % args[0])
        del args[0]
        sys.exit(1)


# Open connection to conditions database (postgres).

try:
    host = os.environ['PUB_PSQL_READER_HOST']
    port = int(os.environ['PUB_PSQL_READER_PORT'])
    user = os.environ['PUB_PSQL_READER_USER']
    db = os.environ['PUB_PSQL_READER_DB']
    pw = os.environ['PUB_PSQL_READER_PASS']
except:
    print("Source ~/.sqlaccess/prod_access.sh first.")


if devel:
    host = 'ifdb04.fnal.gov'
    db = 'microboone_dev'
    calibs = calibs_dev

conn = psycopg2.connect(host=host, port=port, dbname=db, user=user, password=pw)
cur = conn.cursor()

# Loop over calibration databases.

for calib in calibs:

    print('Calibration database %s.' % calib)

    # Create sqlite3 database for thie calibration.

    sqlite_database = '%s.db' % calib
    if os.path.exists(sqlite_database):
        os.remove(sqlite_database)
    sqlite_conn = sqlite3.connect(sqlite_database)
    sqlite_cur = sqlite_conn.cursor()

    # Dictionary of columns.

    schema = {}    # schema[table] = list of columns.

    # Loop over tables.

    for suffix in suffixes:

        # Handle data table separately.

        table_name = '%s_%s' % (calib, suffix)
        schema[table_name] = []
        print('Processing table %s.' % table_name)

        # Construct sqlite query to create corresponding table.

        qtbl = 'CREATE TABLE IF NOT EXISTS %s (' % table_name

        # Get schema of this table.

        q = 'select column_name, data_type, is_nullable from information_schema.columns where table_name=%s'
        cur.execute(q, (table_name,))
        rows = cur.fetchall()
        first = True
        for row in rows:
            column_name = row[0]
            data_type = row[1]
            null = row[2]

            schema[table_name].append(column_name)

            if not first:
                qtbl += ', '

            # Convert postgres data type to sqlite data type.

            sqlite_type = ''
            if data_type == 'integer':
                sqlite_type = 'integer'
            elif data_type == 'bigint':
                sqlite_type = 'integer'
            elif data_type == 'text':
                sqlite_type = 'text'
            elif data_type == 'timestamp without time zone':
                sqlite_type = 'integer'
            elif data_type == 'boolean':
                sqlite_type = 'integer'
            elif data_type == 'real':
                sqlite_type = 'real'
            if sqlite_type == '':
                print('Unknown type %s' % data_type)
                sys.exit(1)

            #print column_name, data_type, sqlite_type

            qtbl += '%s %s' % (column_name, sqlite_type)

            # Add primary keys and constraints.

            if first and (suffix == 'tags' or suffix == 'iovs'):
                qtbl += ' PRIMARY KEY'
            elif null == 'NO':
                qtbl += ' NOT NULL'
            first = False

        # Add foreign keys.

        if suffix == 'tag_iovs':
            qtbl += ', FOREIGN KEY (iov_id) REFERENCES %s_iovs(iov_id)' % calib
            qtbl += ', FOREIGN KEY (tag) REFERENCES %s_tags(tag)' % calib
        elif suffix == 'data':
            qtbl += ', FOREIGN KEY (__iov_id) REFERENCES %s_iovs(iov_id)' % calib

        # Complete query.

        qtbl += ');'
        #print qtbl
        sqlite_cur.execute(qtbl)

        # Done creating table.

        # For data table, query max __iov_id.

        max_iov_id = 0
        if suffix == 'data':
            q = 'SELECT MAX(__iov_id) FROM %s' % table_name
            cur.execute(q)
            row = cur.fetchone()
            max_iov_id = row[0]
        print('Maximum iov_id = %d' % max_iov_id)    

        # Loop over iov ids.

        for iov_id in range(max_iov_id + 1):

            # Query contents of table from postgres database.

            first = True
            q = 'SELECT '
            for column in schema[table_name]:
                if not first:
                    q += ','
                q += column
                first = False

            q += ' FROM %s' % table_name

            # For data table, append iov_id constraint.

            if suffix == 'data':
                print('iov id = %d' % iov_id)
                q += ' WHERE __iov_id = %d' % iov_id
            q += ';'

            #print q
            cur.execute(q)
            rows = cur.fetchall()
            print('%d rows fetched.' % len(rows))
            now = datetime.datetime.now()
            for row in rows:

                # Insert row into sqlite database.

                q = 'INSERT INTO %s (' % table_name
                qval = 'VALUES('
                values = []
                n = 0
                for column in schema[table_name]:
                    element = row[n]
                    if n > 0:
                        q += ','
                        qval += ','
                    q += column
                    qval += '?'
                    if type(element) == type(now):
                        values.append(int(element.strftime('%s')))
                    else:
                        values.append(element)
                    n += 1
                qval += ')'
                q += ') %s;' % qval
                #print q
                #print values
                sqlite_cur.execute(q, tuple(values))


    # Done looping over tables for thie calibration database.

    # Close sqlite database.

    sqlite_conn.commit()
    sqlite_conn.close()

# Done looping over calibration databases.

sys.exit(0)

