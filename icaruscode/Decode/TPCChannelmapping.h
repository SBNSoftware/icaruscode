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
    using CrateNameReadoutIDPair      = std::pair<std::string,ReadoutIDVec>;
    using TPCFragmentIDToReadoutIDMap = std::map<unsigned int, CrateNameReadoutIDPair>;

    inline int BuildTPCFragmentIDToReadoutIDMap(TPCFragmentIDToReadoutIDMap& fragmentBoardMap)
    {
        const unsigned int tpcIdentifier(0x00001000);
        const std::string  name("icarus_hw_readoutboard");
        const std::string  dburl("https://dbdata0vm.fnal.gov:9443/QE/hw/app/SQ/query?dbname=icarus_hardware_dev");
        const std::string  dataType("readout_boards");
        Dataset            dataset;

        // Recover the data from the database
        int error = GetDataset(name,dburl,dataType,dataset);

        // Include a by hand mapping of fragement ID to crate
        using FlangeIDToCrateMap = std::map<size_t,std::string>;

        FlangeIDToCrateMap flangeIDToCrateMap;

        flangeIDToCrateMap[19]  = "WW01T";
        flangeIDToCrateMap[68]  = "WW01M";
        flangeIDToCrateMap[41]  = "WW01B";
        flangeIDToCrateMap[11]  = "WW02";
        flangeIDToCrateMap[17]  = "WW03";
        flangeIDToCrateMap[36]  = "WW04";
        flangeIDToCrateMap[18]  = "WW05";
        flangeIDToCrateMap[58]  = "WW06";
        flangeIDToCrateMap[71]  = "WW07";
        flangeIDToCrateMap[14]  = "WW08";
        flangeIDToCrateMap[25]  = "WW09";
        flangeIDToCrateMap[34]  = "WW10";
        flangeIDToCrateMap[67]  = "WW11";
        flangeIDToCrateMap[33]  = "WW12";
        flangeIDToCrateMap[87]  = "WW13";
        flangeIDToCrateMap[10]  = "WW14";
        flangeIDToCrateMap[59]  = "WW15";
        flangeIDToCrateMap[95]  = "WW16";
        flangeIDToCrateMap[22]  = "WW17";
        flangeIDToCrateMap[91]  = "WW18";
        flangeIDToCrateMap[61]  = "WW19";
        flangeIDToCrateMap[55]  = "WW20T";
        flangeIDToCrateMap[97]  = "WW20M";
        flangeIDToCrateMap[100] = "WW20B";
        flangeIDToCrateMap[83]  = "WE01T";
        flangeIDToCrateMap[85]  = "WE01M";
        flangeIDToCrateMap[7]   = "WE01B";
        flangeIDToCrateMap[80]  = "WE02";
        flangeIDToCrateMap[52]  = "WE03";
        flangeIDToCrateMap[32]  = "WE04";
        flangeIDToCrateMap[70]  = "WE05";
        flangeIDToCrateMap[74]  = "WE06";
        flangeIDToCrateMap[46]  = "WE07";
        flangeIDToCrateMap[81]  = "WE08";
        flangeIDToCrateMap[63]  = "WE09";
        flangeIDToCrateMap[30]  = "WE10";
        flangeIDToCrateMap[51]  = "WE11";
        flangeIDToCrateMap[90]  = "WE12";
        flangeIDToCrateMap[23]  = "WE13";
        flangeIDToCrateMap[93]  = "WE14";
        flangeIDToCrateMap[92]  = "WE15";
        flangeIDToCrateMap[88]  = "WE16";
        flangeIDToCrateMap[73]  = "WE17";
        flangeIDToCrateMap[1]   = "WE18";
        flangeIDToCrateMap[66]  = "WE19";
        flangeIDToCrateMap[48]  = "WE20T";
        flangeIDToCrateMap[13]  = "WE20M";
        flangeIDToCrateMap[56]  = "WE20B";
        flangeIDToCrateMap[94]  = "EW01T";
        flangeIDToCrateMap[77]  = "EW01M";
        flangeIDToCrateMap[72]  = "EW01B";
        flangeIDToCrateMap[65]  = "EW02";
        flangeIDToCrateMap[4]   = "EW03";
        flangeIDToCrateMap[89]  = "EW04";
        flangeIDToCrateMap[37]  = "EW05";
        flangeIDToCrateMap[76]  = "EW06";
        flangeIDToCrateMap[49]  = "EW07";
        flangeIDToCrateMap[60]  = "EW08";
        flangeIDToCrateMap[21]  = "EW09";
        flangeIDToCrateMap[6]   = "EW10";
        flangeIDToCrateMap[62]  = "EW11";
        flangeIDToCrateMap[2]   = "EW12";
        flangeIDToCrateMap[29]  = "EW13";
        flangeIDToCrateMap[44]  = "EW14";
        flangeIDToCrateMap[9]   = "EW15";
        flangeIDToCrateMap[31]  = "EW16";
        flangeIDToCrateMap[98]  = "EW17";
        flangeIDToCrateMap[38]  = "EW18";
        flangeIDToCrateMap[99]  = "EW19";
        flangeIDToCrateMap[53]  = "EW20T";
        flangeIDToCrateMap[82]  = "EW20M";
        flangeIDToCrateMap[35]  = "EW20B";
        flangeIDToCrateMap[96]  = "EE01T";
        flangeIDToCrateMap[28]  = "EE01M";
        flangeIDToCrateMap[16]  = "EE01T";
        flangeIDToCrateMap[69]  = "EE02";
        flangeIDToCrateMap[20]  = "EE02";
        flangeIDToCrateMap[79]  = "EE02";
        flangeIDToCrateMap[50]  = "EE02";
        flangeIDToCrateMap[45]  = "EE02";
        flangeIDToCrateMap[84]  = "EE02";
        flangeIDToCrateMap[42]  = "EE02";
        flangeIDToCrateMap[39]  = "EE02";
        flangeIDToCrateMap[26]  = "EE02";
        flangeIDToCrateMap[64]  = "EE02";
        flangeIDToCrateMap[43]  = "EE02";
        flangeIDToCrateMap[47]  = "EE02";
        flangeIDToCrateMap[15]  = "EE02";
        flangeIDToCrateMap[3]   = "EE02";
        flangeIDToCrateMap[27]  = "EE02";
        flangeIDToCrateMap[24]  = "EE02";
        flangeIDToCrateMap[40]  = "EE02";
        flangeIDToCrateMap[75]  = "EE02";
        flangeIDToCrateMap[86]  = "EE20T";
        flangeIDToCrateMap[54]  = "EE20M";
        flangeIDToCrateMap[8]   = "EE20B";

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

                if (fragmentBoardMap.find(fragmentID) == fragmentBoardMap.end())
                {
                    unsigned int flangeID = getLongValue(tuple, 1, &error);

                    if (error) throw std::runtime_error("Encountered error in trying to recover Board Flange ID from database");

                    fragmentBoardMap[fragmentID].first = flangeIDToCrateMap[flangeID];
                }

                unsigned int readoutID = getLongValue(tuple, 0, &error);

                if (error) throw std::runtime_error("Encountered error in trying to recover Board ReadoutID from database");

                fragmentBoardMap[fragmentID].second.emplace_back(readoutID);

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
    
    using ChannelPlanePair            = std::pair<unsigned int, unsigned int>;
    using ChannelPlanePairVec         = std::vector<ChannelPlanePair>;
    using SlotChannelVecPair          = std::pair<unsigned int, ChannelPlanePairVec>;
    using TPCReadoutBoardToChannelMap = std::map<unsigned int, SlotChannelVecPair>;

    const unsigned int CHANNELSPERBOARD = 64;

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

        // Loop through the data to recover the channels, making sure to skip the first (header) row
        for(int row = 1; row < getNtuples(dataset); row++)
        {
            // Recover the row
            Tuple tuple = getTuple(dataset, row);

            if (tuple != NULL)
            {
                unsigned int readoutBoardID   = getLongValue(tuple, 2, &error);

                if (error) throw std::runtime_error("Encountered error when trying to read Board ReadoutID");

                if (rbChanMap.find(readoutBoardID) == rbChanMap.end())
                {
                    unsigned int readoutBoardSlot = getLongValue(tuple, 4, &error);

                    if (error) throw std::runtime_error("Encountered error when trying to read Board Readout slot");

                    rbChanMap[readoutBoardID].first = readoutBoardSlot;
                    rbChanMap[readoutBoardID].second.resize(CHANNELSPERBOARD);
                }

                unsigned int channelNum = getLongValue(tuple, 5, &error);

                if (error) throw std::runtime_error("Encountered error when trying to read channel number");

                unsigned int channelID = getLongValue(tuple, 0, &error);

                if (error) throw std::runtime_error("Encountered error when recovering the channel ID list");

                // Recover the plane identifier 
                char fragmentBuffer[16];

                getStringValue(tuple, 10, fragmentBuffer, sizeof(fragmentBuffer), &error);

                if (error) throw std::runtime_error("Encountered error when trying to read plane type");

                // Make sure lower case... (sigh...)
                for(size_t charIdx = 0; charIdx < sizeof(fragmentBuffer); charIdx++) fragmentBuffer[charIdx] = tolower(fragmentBuffer[charIdx]);

                unsigned int plane(3);

                if      (strstr(fragmentBuffer,"collection"))  plane = 2;
                else if (strstr(fragmentBuffer,"induction 2")) plane = 1;
                else if (strstr(fragmentBuffer,"induction 1")) plane = 0;

                if (plane > 2) std::cout << "YIKES!!! Plane is " << plane << " for channel " << channelID << " with type " << std::string(fragmentBuffer) << std::endl;

                rbChanMap[readoutBoardID].second[channelNum] = ChannelPlanePair(channelID,plane);
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

    // Ultimately, we want to map between Fragment IDs and a pair of numbers - the channel in the fragment and the LArSoft channel ID
    using DigitizerChannelChannelIDPair    = std::pair<size_t,size_t>;
    using DigitizerChannelChannelIDPairVec = std::vector<DigitizerChannelChannelIDPair>;
    using FragmentToDigitizerChannelMap    = std::map<size_t,DigitizerChannelChannelIDPairVec>;

    inline int BuildFragmentToDigitizerChannelMap(FragmentToDigitizerChannelMap& fragmentToDigitizerChannelMap)
    {
        // clearing is cleansing
        fragmentToDigitizerChannelMap.clear();

        // We need a mapping between a fragmentID and the digitizer label which will be available in the database
        using DigitizerToFragmentMap = std::map<std::string,size_t>;

        DigitizerToFragmentMap digitizerToFragmentMap;

        // At the time of writing the database does not contain the fragment IDs, it contains the digitizer name... sigh... so we build a map here
        // Build the fragment ID to digitizer map by hand
        digitizerToFragmentMap["WW-TOP-A"] = 0;
        digitizerToFragmentMap["WW-TOP-B"] = 1;
        digitizerToFragmentMap["WW-TOP-C"] = 2;
        digitizerToFragmentMap["WW-BOT-A"] = 3;
        digitizerToFragmentMap["WW-BOT-B"] = 4;
        digitizerToFragmentMap["WW-BOT-C"] = 5;
        digitizerToFragmentMap["WE-TOP-A"] = 6;
        digitizerToFragmentMap["WE-TOP-B"] = 7;
        digitizerToFragmentMap["WE-TOP-C"] = 8;
        digitizerToFragmentMap["WE-BOT-A"] = 9;
        digitizerToFragmentMap["WE-BOT-B"] = 10;
        digitizerToFragmentMap["WE-BOT-C"] = 11;
        digitizerToFragmentMap["EW-TOP-A"] = 12;
        digitizerToFragmentMap["EW-TOP-B"] = 13;
        digitizerToFragmentMap["EW-TOP-C"] = 14;
        digitizerToFragmentMap["EW-BOT-A"] = 15;
        digitizerToFragmentMap["EW-BOT-B"] = 16;
        digitizerToFragmentMap["EW-BOT-C"] = 17;
        digitizerToFragmentMap["EE-TOP-A"] = 18;
        digitizerToFragmentMap["EE-TOP-B"] = 19;
        digitizerToFragmentMap["EE-TOP-C"] = 20;
        digitizerToFragmentMap["EE-BOT-A"] = 21;
        digitizerToFragmentMap["EE-BOT-B"] = 22;
        digitizerToFragmentMap["EE-BOT-C"] = 23;

//        std::cout << "PMT local map has " << digitizerToFragmentMap.size() << " rows" << std::endl;
//        for(const auto& mapPair : digitizerToFragmentMap) std::cout << "  - label: " << mapPair.first << ", fragment: " << mapPair.second << ", size: " << mapPair.first.size() << std::endl;

        // Recover the information from the database on the mapping 
        const std::string  name("icarus_hw_readoutboard");
        const std::string  dburl("https://dbdata0vm.fnal.gov:9443/QE/hw/app/SQ/query?dbname=icarus_hardware_dev");
        const std::string  dataType("pmt_placements");
        Dataset            dataset;

        // Recover the data from the database
        int error = GetDataset(name,dburl,dataType,dataset);

        // If there was an error the function above would have printed a message so bail out
        if (error) throw(std::exception());

        // Ok, now we can start extracting the information
        // We do this by looping through the database and building the map from that
        for(int row = 1; row < getNtuples(dataset); row++)
        {
            // Recover the row
            Tuple tuple = getTuple(dataset, row);

            if (tuple != NULL)
            {
                char digitizerBuffer[10];

                // Recover the digitizer label first 
	            getStringValue(tuple, 8, digitizerBuffer, sizeof(digitizerBuffer), &error);

                if (error) throw std::runtime_error("Encountered error when trying to recover the PMT digitizer label");

	            std::string digitizerLabel(digitizerBuffer, 8); //sizeof(digitizerBuffer));

                // Recover the fragment id
                unsigned fragmentID = getLongValue(tuple, 16, &error);

                if (error) throw std::runtime_error("Encountered error when trying to recover the PMT fragment id");

                // Now recover the digitizer channel number
                unsigned int digitizerChannelNo = getLongValue(tuple, 9, &error);

                if (error) throw std::runtime_error("Encountered error when trying to recover the PMT digitizer channel number");

                // Finally, get the LArsoft channel ID
                unsigned int channelID = getLongValue(tuple, 17, &error);

                if (error) throw std::runtime_error("Encountered error when trying to recover the PMT channel ID");

//                std::cout << "  >> Searching " << digitizerToFragmentMap.size() << " rows for " << digitizerLabel << ", size; " << digitizerLabel.size() << std::endl;

                // Do we have corresponence?
                if (digitizerToFragmentMap.find(digitizerLabel) == digitizerToFragmentMap.end())
                {
                    std::cout << "No match in map for label: " << digitizerLabel << std::endl;
                    throw std::runtime_error("Could not find fragment ID corresponding to digitizer label");
                }

                // Fill the map
                fragmentToDigitizerChannelMap[fragmentID].emplace_back(digitizerChannelNo,channelID);

                releaseTuple(tuple);
            }
        }

        return error;
    }

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
