#----------------------------------------------------------
#
# This was developed by SCD of Fermilab.
# This file was copied from here
# https://cdcvs.fnal.gov/redmine/projects/hardwaredb/files
#
#-----------------------------------------------------------


import base64
import hashlib
import json
import os
import random
import ssl
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from io import IOBase


class DataLoader(object):
    """Collects data intended for a hardware database table.  On command, the data
       is sent to a server for loading."""

    def __init__(self, password, url, group, table):
        """ Class constructor.

            Args:
                 password - Agreed upon password - for the group.
                 url - Http URL to the server used for loading.
                 group - Group the table is part of.
                 table - Postgresql database table data will be loaded into.
        """
        self.password = password
        self.group = group
        self.url = url
        self.args = "table=%s" % table
        self.urlWithArgs = "%s" % self.url
        self.data = {'table': table.lower(),
                     'rows':  []
                     }

    def addRow(self, row, mode='insert'):
        """ Adds a single row of data to the instance.  This row will be
            inserted or updated in the database.

            Args:
                 row - a dictionary containing a name/value pair
                       for each required table column.
                 mode - insert or update
        """
        if isinstance(row, dict) is False:
            raise Exception("row must be a dictionary")
        if mode not in ('insert', 'update'):
            raise Exception("mode must be insert or update")
        data = {k: (v.name, (lambda f: f.seek(0) or base64.b64encode(f.read()).decode())(v)) if isinstance(v, IOBase) else v for (k, v) in list(row.items())}
        (self.data['rows']).append((mode, data))

    def send(self, echoUrl=False):
        """Sends the data to the server for loading.

           Returns:
              Boolean indicating success and failure of the call.
              A code indicating Html return status.
              Text describing any error which returned.
        """
        # Repeats if there is a collision on the salt.
        while 1:
            jdata = json.dumps(self.data)
            random.seed(time.time())
            salt = '%s' % (random.random(),)
            sig = self.__signature(jdata, salt)
            # The Request call is sending as a POST, not as a GET.
            req = urllib.request.Request(self.urlWithArgs, jdata.encode(),
                                  {'X-Salt':      salt,
                                   'X-Signature': sig,
                                   'X-Group':     self.group,
                                   'X-Table':     self.data['table']
                                   })
            if echoUrl:
                print("URL: %s\n  %s" % (req.get_full_url(), req.header_items()))
            ssl_cert_file = os.environ.get("SSL_CERT_FILE")
            try:
                response = urllib.request.urlopen(req) if ssl_cert_file else urllib.request.urlopen(req, context=ssl.SSLContext())
            except urllib.error.HTTPError as val:
                if self.urlWithArgs.lower().startswith("https") and ssl_cert_file:
                    print("\n*** Please verify CA certificate provided with SSL_CERT_FILE environment variable!\n")
                retValue = False
                code = "%s %s" % (val.code, val.msg)
                text = val.read()
            else:
                retValue = True
                code = "%s %s" % (response.getcode(), response.msg)
                text = response.read()
            if text != "Signature Error":
                break
        return retValue, code, text

    def clearRows(self):
        """ Deletes all rows from the instance, readying it for
            the next set of data.
        """
        self.data['rows'] = []

    def __buildReq(self):
        return req

    def __signature(self, data, salt):
        m = hashlib.md5()
        m.update(self.password.encode())
        m.update(salt.encode())
        m.update(data.encode())
        return m.hexdigest()

    def __str__(self):
        retVal = "URL: %s\nURL with Args: %s\nTable:%s\nPassword: XXXXX\nGroup:%s\n" % (
            self.url, self.urlWithArgs, self.data['table'], self.group
        )
        rowCnt = 0
        rows = self.data['rows']
        if len(rows) == 0:
            retVal += "Rows: None\n"
        else:
            for row in rows:
                retVal += "Row %s:\n" % rowCnt
                for column in list(row.keys()):
                    retVal += "    %s: %s\n" % (column, str(row.get(column)))
                rowCnt += 1
        return retVal


class DataQuery:
    """ Supports simple user queries through the use of QueryEngine.
        (https://cdcvs.fnal.gov/redmine/projects/qengine/wiki)
    """

    def __init__(self, url):
        """ Class constructor.

            Args:
                 url - Http URL to QueryEngine.
        """
        self.url = url

    def query(self, database, table, columns, where=None, order=None, limit=None, echoUrl=False):
        """ Executes a simple query and returns the results in a list.  List data will
            be in the same order as listed in the columns attribute.

            Args:
                 database - The name of the database to be queried.  (This database must
                            be in QueryEngine's configuration file.)
                 table - The name of the table to query on.
                 columns - A comma separated string of the table columns to be returned.
                 where - (optional) <column>:<op>:<value> - can be repeated; seperated by ampersand (&)
                         op can be: lt, le, eq, ne, ge, gt
                 order - (optional) A comma separated string of columns designating row order in the returned list.
                         Start the string with a minus (-) for descending order.
                 limit - (optional) - A integer designating the maximum number of rows to be returned.
        """

        parameters = [
            ('dbname', database),
            ('t', table),
            ('c', columns),
            ('x', 'no'),
        ]
        if where is not None:
            parameters.append(('w', where))
        if order is not None:
            parameters.append(('o', order))
        if limit is not None:
            parameters.append(('l', limit))
        fullUrl = self.url + '?' + urllib.parse.urlencode(parameters, doseq=True)
        if echoUrl:
            print(("Url: %s" % fullUrl))
        req = urllib.request.Request(fullUrl)
        ssl_cert_file = os.environ.get("SSL_CERT_FILE")
        try:
            resp = urllib.request.urlopen(req) if ssl_cert_file else urllib.request.urlopen(req, context=ssl.SSLContext())
        except urllib.error.HTTPError as val:
            retValue = False
            code = "%s %s" % (val.code, val.msg)
            text = val.read().decode()
            print(f"\n*** ", code, text)
        except urllib.error.URLError:
            if fullUrl.lower().startswith("https") and ssl_cert_file:
                print("\n*** Please verify CA certificate provided with SSL_CERT_FILE environment variable!\n")
            pass
        else:
            text = resp.read().decode()
            data = text.split('\n')
            return data[1:]

        return []
