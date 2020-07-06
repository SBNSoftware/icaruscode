// -------------------------------------------------
// This is the first pass to this script
// It will improve gradually as we move on 
// Also coding style will gradually improve.
//
//
// Author : Biswaranjan Behera @ 02 April 2020
// email  : bishu@colostate.edu
//-------------------------------------------------

#include "wda.h"

#include <time.h>
#include <vector> 
#include <string>
#include <iostream> 
#include <stdio.h>

namespace database
{

    // -----------------------------------------------------
    // This Function does the basic information retrieval 
    // One passes in the type of data requested and a reference to the data holder
    // The function connects via the libwda function to recover the
    // data and checks to insure there were not connection errors. 
    // The function returns the error status, if null then success
    //-----------------------------------------------------

    inline int GetDataset(const std::string name, const std::string url, const std::string& dataType, Dataset& dataSet)
    {
        const int   timeout(200);
        std::string dburl = url + "&t=" + dataType;

        int error(0);

        dataSet = getDataWithTimeout(dburl.data(),name.data(),timeout,&error);

        if (error)
        {
            std::string errorMsg = "Database access GetDataset failed with error " + std::to_string(error) + "\nDB url: " 
                                 + dburl + ", name: " + name + ", type: " + dataType;
            std::cout << "****> Database retrieval error, code: " << error << std::endl;
            throw std::runtime_error(errorMsg);
        }

        return error;
    }

    // -----------------------------------------------------
    // The aim of this function is to build a map between the
    // TPC Fragment IDs and the readout board IDs. Here we 
    // expect there will be some number of boards per Fragment
    //-----------------------------------------------------

    using ReadoutIDVec                = std::vector<unsigned int>;
    using TPCFragmentIDToReadoutIDMap = std::map<unsigned int, ReadoutIDVec>;

    inline int BuildTPCFragmentIDToReadoutIDMap(TPCFragmentIDToReadoutIDMap& fragmentBoardMap)
    {
        const unsigned int tpcIdentifier(0x00001000);
        const std::string  name("icarus_hw_readoutboard");
        const std::string  dburl("https://dbdata0vm.fnal.gov:9443/QE/hw/app/SQ/query?dbname=icarus_hardware_dev");
        const std::string  dataType("readout_boards");
        Dataset            dataset;

        // Recover the data from the database
        int error = GetDataset(name,dburl,dataType,dataset);

        // If there was an error the function above would have printed a message so bail out
        if (error) throw(std::exception());

        // Loop through the data to recover the channels
        // NOTE that we skip the first row because that is just the labels
        for(int row = 1; row < getNtuples(dataset); row++)
        {
            // Recover the row
            Tuple tuple = getTuple(dataset, row);

            if (tuple != NULL)
            {
                // Note that the fragment ID is stored in the database as a string which reads as a hex number
                // Meaning we have to read back as a string and decode to get the numerical value. 
                char fragmentBuffer[16];

                getStringValue(tuple, 8, fragmentBuffer, sizeof(fragmentBuffer), &error);

                if (error) throw std::runtime_error("Encountered error in trying to recover FragmentID from database");

                std::string fragmentIDString(fragmentBuffer,4);

                unsigned int fragmentID = std::stol(fragmentIDString,nullptr,16);

                if (!(fragmentID & tpcIdentifier)) continue;

                unsigned int readoutID = getLongValue(tuple, 0, &error);

                if (error) throw std::runtime_error("Encountered error in trying to recover Board ReadoutID from database");

                fragmentBoardMap[fragmentID].emplace_back(readoutID);

                releaseTuple(tuple);
            }
        }

        return error;
    }

    // -----------------------------------------------------
    // The aim of this function is to build a map between the
    // TPC readout board IDs and the associated channels. So for
    // each readout board ID we expect a number of channels back
    // from the data base. So the returned data structure will 
    // be a map of readout ID to a vector of channels
    //-----------------------------------------------------
    
    using ChannelVec                  = std::vector<unsigned int>;
    using TPCReadoutBoardToChannelMap = std::map<unsigned int, ChannelVec>;

    inline int BuildTPCReadoutBoardToChannelMap(TPCReadoutBoardToChannelMap& rbChanMap)
    {
        const std::string  name("icarus_hardware_dev");
        const std::string  dburl("https://dbdata0vm.fnal.gov:9443/QE/hw/app/SQ/query?dbname=icarus_hardware_dev");
        const std::string  dataType("daq_channels");
        Dataset            dataset;

        // Recover the data from the database
        int error = GetDataset(name,dburl,dataType,dataset);

        // If there was an error the function above would have printed a message so bail out
        if (error) throw std::runtime_error("Encountered error accessing the database with GetDataset");

        // Loop through the data to recover the channels
        for(int row = 0; row < getNtuples(dataset); row++)
        {
            // Recover the row
            Tuple tuple = getTuple(dataset, row);

            if (tuple != NULL)
            {
                unsigned int readoutBoardID = getLongValue(tuple, 2, &error);

                if (error) throw std::runtime_error("Encountered error when trying to read Board ReadoutID");

                unsigned int channelID = getLongValue(tuple, 0, &error);

                if (error) throw std::runtime_error("Encountered error when recovering the channel ID list");

                rbChanMap[readoutBoardID].emplace_back(channelID);
            }
        }

        return error;
    }
    
    // -------------------------------------------------
    // This Function returns Fragment_ID.
    // input : readout_board_id
    //------------------------------------------------- 
    
    inline char* Fragment_ID(int readout_board_ID)
    {
        Dataset ds;
        const std::string name("icarus_hw_readoutboard");
        const std::string dburl("https://dbdata0vm.fnal.gov:9443/QE/hw/app/SQ/query?dbname=icarus_hardware_dev");
        const std::string dataType = "readout_boards";

        int     error = GetDataset(name,dburl,dataType, ds);
    
        if (error) 
        {                                // Check for curl library errors
            fprintf(stderr, "error code=%d\n", error);    perror("error message");
        }
        if (getHTTPstatus(ds) != 200) 
        {                        // Check for HTTP error
            fprintf(stderr, "HTTP code=%ld, message: '%s'\n", getHTTPstatus(ds), getHTTPmessage(ds));
        }
    
        char ss[10];
        int err;
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
    	              getStringValue(tup, 8 ,ss, sizeof (ss), &err);
                }
                releaseTuple(tup);
            }
        }
        char* fragment_id = new char[4];
        std::strcpy(fragment_id, ss);
        return fragment_id;
    }
    
    // -------------------------------------------------
    // This Function returns vector of Channel_ID,
    // As One read out board contains 64 channels.
    // input : readout_board_id
    //-------------------------------------------------
    
    inline std::vector<int> Channel_ID(int readout_board_ID)
    {
        Dataset ds;
        const std::string name("icarus_hardware_dev");
        const std::string dburl("https://dbdata0vm.fnal.gov:9443/QE/hw/app/SQ/query?dbname=icarus_hardware_dev");
        const std::string dataType = "daq_channels";

        int     error = GetDataset(name,dburl,dataType, ds);
     
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
            Tuple tu = getTuple(ds, rows);                                           // Get the row with double array 
    
            // If we get the names print them
            if (tu != NULL) { 
                int err(0); 
                // If everything is OK
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


    //******************* PMT Channel Mapping ***********************

    //---------------------------------------------------------------
    // The aim of this function is to build a map between the
    // Larsoft PMT Channel ID and the physical location of PMT ID
    //---------------------------------------------------------------

    inline int PMTChannelIDFromPhysicalPMTPositionID(unsigned int pmtid)
    {
      const std::string  name("icarus_hw_readoutboard");
      const std::string  dburl("https://dbdata0vm.fnal.gov:9443/QE/hw/app/SQ/query?dbname=icarus_hardware_dev");
      const std::string  dataType("pmt_placements");
      Dataset            dataset;

      // Recover the data from the database
      int error = GetDataset(name,dburl,dataType,dataset);

      // If there was an error the function above would have printed a message so bail out
      if (error) throw(std::exception());

      unsigned int PhysicalPMTPositionID;
      unsigned int PMTChannelID;

      // Loop through the data to recover the channels
      // NOTE that we skip the first row because that is just the labels
      for(int row = 1; row < getNtuples(dataset); row++)
        {
          // Recover the row
          Tuple tuple = getTuple(dataset, row);

          if (tuple != NULL)
            {

              PhysicalPMTPositionID = getLongValue(tuple, 0, &error);

	      if (error) throw std::runtime_error("Encountered error when trying to read Physical PMT Position ID from database");

              if (PhysicalPMTPositionID == pmtid) {

                PMTChannelID = getLongValue(tuple, 17, &error);

                if (error) throw std::runtime_error("Encountered error in trying to recover Larsoft PMT Channel ID from database");

              }else continue;

              releaseTuple(tuple);
            }
        }

      return PMTChannelID;
    }

    //-----------------------------------------------------------------------------
    // The aim of this function is to build a map between the
    // Larsoft PMT Channel ID and the Digitizer channel number and Digitizer label
    //-----------------------------------------------------------------------------

    inline int PMTChannelIDFromDigitizer(std::string digitizerlabel, int digitizerchannelnumber)
    {
      const std::string  name("icarus_hw_readoutboard");
      const std::string  dburl("https://dbdata0vm.fnal.gov:9443/QE/hw/app/SQ/query?dbname=icarus_hardware_dev");
      const std::string  dataType("pmt_placements");
      Dataset            dataset;

      // Recover the data from the database
      int error = GetDataset(name,dburl,dataType,dataset);

      // If there was an error the function above would have printed a message so bail out
      if (error) throw(std::exception());

      unsigned int PMTChannelID;

      // Loop through the data to recover the channels
      // NOTE that we skip the first row because that is just the labels
      for(int row = 1; row < getNtuples(dataset); row++)
        {
          // Recover the row
          Tuple tuple = getTuple(dataset, row);

          if (tuple != NULL)
            {
              char digitizerBuffer[10];

	      getStringValue(tuple, 8, digitizerBuffer, sizeof(digitizerBuffer), &error);
	      std::string digitizerLabel(digitizerBuffer, sizeof(digitizerBuffer));

              int digitizerChannelNumber = getDoubleValue(tuple, 9, &error);


              int matchedcharacter = digitizerLabel.compare(digitizerlabel);

              // Not sure what I am doing wrong why 2 is coming for matchedcharacter,
              // but I have tried a simple compare string with same name it is 0.
              // std::string s1("EE-BOT-B");
              // std::string s2("EE-BOT-B");
              // int compare = s1.compare(s2);
              // if (compare == 0) std::cout << "strings are equal" << std::endl;

              if(matchedcharacter == 2 && digitizerChannelNumber == digitizerchannelnumber){

                PMTChannelID = getLongValue(tuple, 17, &error);

                if (error) throw std::runtime_error("Encountered error in trying to recover Larsoft PMT Channel ID from database");

              }else  continue;

	      releaseTuple(tuple);
            }
        }

      return PMTChannelID;

    }
    
    // -------------------------------------------------
    // This is the main Function
    //-------------------------------------------------
    
    //int main()
    //{
    //  std::vector<int> ch_id_vec =  Channel_ID(854);
    //  for(std::vector<int>::iterator it = ch_id_vec.begin(); it != ch_id_vec.end(); it++) {
    //    std::cout << "channel_id: \t" << *it << std::endl;
    //  }
    // 
    //  std ::cout << Fragment_ID(89)  << std::endl; 
    //
    //  std ::cout << Fragment_ID(854) << std::endl;
    //
    //  return 0;
    //}

} // end of namespace
