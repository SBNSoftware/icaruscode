// -------------------------------------------------
// TPCChannelmapping.cxx 
// This is the first pass to this script
// It will improve gradually as we move on 
// Also coding style will gradually improve.
//
//
// Author : Biswaranjan Behera @ 02 April 2020
// email  : bishu@colostate.edu
//-------------------------------------------------

#include <time.h>
#include "wda.h"
#include <vector> 
#include <string>
#include <iostream> 
#include <stdio.h>
#include <cstring>

// -----------------------------------------------------
// This Function connects to Hardware Datbase TPC tables
// inputs : 1. name (for convinient name of table)
//          2. URL 
//          3. passing int error for helping printing 
//             on error messages 
//-----------------------------------------------------

Dataset ConnectDB(const char* name, const char* url, int* error)
{
  Dataset dbtable;
  dbtable =  getData(url, name, error);
  return dbtable;
}

// -------------------------------------------------
// This Function returns Fragment_ID.
// input : readout_board_id
//------------------------------------------------- 

char* Fragment_ID(int readout_board_ID)
{
  int error;
  Dataset ds = ConnectDB("icarus_hw_readoutboard", "https://dbdata0vm.fnal.gov:9443/QE/hw/app/SQ/query?dbname=icarus_hardware_dev&t=readout_boards", &error);

  if (error) {                                // Check for curl library errors
    fprintf(stderr, "error code=%d\n", error);    perror("error message");
  }
  if (getHTTPstatus(ds) != 200) {                        // Check for HTTP error
    fprintf(stderr, "HTTP code=%ld, message: '%s'\n", getHTTPstatus(ds), getHTTPmessage(ds));
  }

  char ss[10];
  int err;
  int fragmentid;
  int readoutboard_id;
  int nrows =  getNtuples(ds);
  Tuple tup;
  for (int rows = 0; rows < nrows; rows++ ){
    tup = getTuple(ds, rows);                                           // Get the row with double array
    
    if (tup != NULL) {
      //int nc = getNfields(tup);
      //for (i = 0; i < nc; i++) { 
      //len = getStringValue(tup, i ,8, sizeof (ss), &err);
      //fprintf(stderr, "[%d]: l=%d, s='%s'\n", i, len, ss); 
      //fprintf(stderr, "[4]: v=%i\n", getStringValue(tup, i ,ss, sizeof (ss), &err));
      // }
      readoutboard_id =  (int)getDoubleValue(tup, 0, &err);
      if (readoutboard_id == readout_board_ID){
	fragmentid = getStringValue(tup, 8 ,ss, sizeof (ss), &err);
      }
      releaseTuple(tup);
    }
  }
  std::cout << "length of string fragmentid: " << fragmentid << std::endl;
  char* fragment_id = new char[4];
  std::strcpy(fragment_id, ss);
  return fragment_id;
}

// -------------------------------------------------
// This Function returns vector of Channel_ID,
// As One read out board contains 64 channels.
// input : readout_board_id
//-------------------------------------------------

std::vector<int> Channel_ID(int readout_board_ID)
{
  Tuple tu;
  int error;
  int err;
  Dataset ds = ConnectDB("icarus_hardware_dev", "https://dbdata0vm.fnal.gov:9443/QE/hw/app/SQ/query?dbname=icarus_hardware_dev&t=daq_channels", &error);
 
  if (error) {                                // Check for curl library errors
    fprintf(stderr, "error code=%d\n", error);    perror("error message");
  }
  if (getHTTPstatus(ds) != 200) {                        // Check for HTTP error
    fprintf(stderr, "HTTP code=%ld, message: '%s'\n", getHTTPstatus(ds), getHTTPmessage(ds));
  }

  int channel_id;
  std::vector<int> channelidvec;
  int nrows =  getNtuples(ds);
  for (int rows = 0; rows < nrows; rows++ ){
    tu = getTuple(ds, rows);                                           // Get the row with double array 

    // If we get the names print them
    if (tu != NULL) {                                               // If everything is OK
      int readoutboardid = (int)getDoubleValue(tu, 2, &err);
      
      if (readoutboardid == readout_board_ID){
	channel_id = (int)getDoubleValue(tu, 0, &err);
	channelidvec.push_back(channel_id);
      }
      releaseTuple(tu);
    } 
  }
  return channelidvec;
}

// -------------------------------------------------
// This is the main Function
//-------------------------------------------------

int main()
{
  std::vector<int> ch_id_vec =  Channel_ID(854);
  for(std::vector<int>::iterator it = ch_id_vec.begin(); it != ch_id_vec.end(); it++) {
    std::cout << "channel_id: \t" << *it << std::endl;
  }
 
  std ::cout << Fragment_ID(89)  << std::endl; 

  std ::cout << Fragment_ID(854) << std::endl;

  return 0;
}
