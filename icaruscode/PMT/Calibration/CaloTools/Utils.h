////////////////////////////////////////////////////////////////////////////////
// This is the space where spurios functions can be defined
//
// maito: ascarpell@bnl.gov
//
////////////////////////////////////////////////////////////////////////////////

#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <iterator>
#include <string>
#include <algorithm>

using namespace std;

namespace pmtcalo
{

  //////////////////////////////////////////////////////////////////////////////

  class CSVReader
  {
	  string m_filename;
	  char m_delimeter;

    public:
	     CSVReader(string filename, char delm = ',')
			 : m_filename(filename), m_delimeter(delm)
	     { }

	     // Function to fetch data from a CSV File
	     vector<vector<string> > getData()
       {
          ifstream file(m_filename);
          vector<vector<string> > dataList;

          string line = "";
          // Iterate through each line and split the content using delimeter
          while (getline(file, line))
          {
            istringstream s_stream(line);
            string token;
            vector<string> vec;
            while (getline(s_stream, token, m_delimeter))
            {
              vec.push_back(token);
            }
            dataList.push_back(vec);
          }
          // Close the File
          file.close();

          return dataList;
        };
  };


  //----------------------------------------------------------------------------


  int fileProgNumber( std::string filename )
  {

    // This function returns the progressive filenumber given  the structure
    // *_*_run1067_10_*_*.root. It assumes the progressive filenumber follows
    // the run number and that run is explicitly written

    int number=0;

    char delm('_');
    istringstream s_stream(filename);
    string token;
    vector<string> vec;
    int pos=0, count = 0;

    while ( std::getline(s_stream, token, delm) )
    {
      if( token.find("run") < token.size() ){ pos = count; }
      vec.push_back(token);
      count++;
    }

    number = stoi( vec[pos+1] );

    // Trow Exception if stoi not good

    return number;
  }

  //----------------------------------------------------------------------------


  int getBoardNumber( string board_name )
  {

    int board_num = -1;

    std::vector<string> board_names = { "WW-Top-A", "WW-Top-B", "WW-Top-C",
                                        "WW-Bot-A", "WW-Bot-B", "WW-Bot-C",
                                        "WE-Top-A", "WE-Top-B", "WE-Top-C",
                                        "WE-Bot-A", "WE-Bot-B", "WE-Bot-C",
                                        "EW-Top-A", "EW-Top-B", "EW-Top-C",
                                        "EW-Bot-A", "EW-Bot-B", "EW-Bot-C",
                                        "EE-Top-A", "EE-Top-B", "EE-Top-C",
                                        "EE-Bot-A", "EE-Bot-B", "EE-Bot-C" };

    std::vector<string>::iterator it
              = std::find( board_names.begin(), board_names.end(), board_name );
    if( it != board_names.end() )
    {
      board_num = std::distance( board_names.begin(), it );
    }

    return board_num;

  }


  //----------------------------------------------------------------------------


  bool isIlluminated( int optical_channel, int pmt_number )
  {

    bool m_isIlluminated=false;

    int laser_groups=10;

    // We compare the optical channel with the absolute pmt number.
    // NB: This assumption works only in the case all boards are taken and optical
    // channel number is changed in progressive increasing order.
    if( (pmt_number>=optical_channel*laser_groups)
                              & (pmt_number<(optical_channel+1)*laser_groups) ){
        m_isIlluminated=true;
    }

    return m_isIlluminated;
  }

}

#endif //UTILS_H
