#include "icaruscode/CRT/CRTUtils/CRTBackTracker.h"

namespace icarus{
 namespace crt{

    CRTBackTracker::CRTBackTracker(const Config& config){
        this->reconfigure(config);
    }
    
    CRTBackTracker::CRTBackTracker(){}
    CRTBackTracker::~CRTBackTracker(){}

    //------------------------------------------------------------------------------    
    void CRTBackTracker::reconfigure(const Config& config){
    
        fCRTTrueHitLabel = config.CRTTrueHitLabel();
        fCRTDataLabel = config.CRTDataLabel();
        fCRTSimHitLabel = config.CRTSimHitLabel();
        fCRTTrackLabel = config.CRTTrackLabel();
        fRollupUnsavedIds = config.RollupUnsavedIds();
      
        return;
    }

    //-----------------------------------------------------------------------------    
    void CRTBackTracker::Initialize(const art::Event& event){
    
        // Clear those data structures!
        fTrueHitTrueIds.clear();
        fDataTrueIds.clear();
        fSimHitTrueIds.clear();
        fTrackTrueIds.clear();
   
        //CRTTrueHit
        art::Handle< std::vector<CRTHit>> crtTrueHitHandle;
        std::vector<art::Ptr<CRTHit> > crtTrueHitList;
        if (event.getByLabel(fCRTTrueHitLabel, crtTrueHitHandle)){
            art::fill_ptr_vector(crtTrueHitList, crtTrueHitHandle);

            art::FindManyP<sim::AuxDetIDE> findManyIdes(crtTrueHitHandle, event, fCRTTrueHitLabel);
            std::map<art::Ptr<CRTHit>, int> trueHitPtrMap;

            for(size_t hit_i = 0; hit_i < crtTrueHitList.size(); hit_i++){

                trueHitPtrMap[crtTrueHitList[hit_i]] = hit_i;
                std::vector<art::Ptr<sim::AuxDetIDE>> ides = findManyIdes.at(hit_i);

                for(size_t i = 0; i < ides.size(); i++){

                    int id = ides[i]->trackID; 
                    if(fRollupUnsavedIds) id = std::abs(id);
                    fTrueHitTrueIds[hit_i][id] += ides[i]->energyDeposited;
                }
            }
        }//if CRTTrueHit found
 
        //CRTData
        art::Handle< std::vector<CRTData>> crtDataHandle;
        std::vector<art::Ptr<CRTData> > crtDataList;
        std::map<art::Ptr<CRTData>, int> dataPtrMap;

        if (event.getByLabel(fCRTDataLabel, crtDataHandle)) {
            art::fill_ptr_vector(crtDataList, crtDataHandle);
            art::FindManyP<sim::AuxDetIDE> findManyIdes(crtDataHandle, event, fCRTDataLabel);

            for(size_t data_i = 0; data_i < crtDataList.size(); data_i++){
      
                dataPtrMap[crtDataList[data_i]] = data_i;
                // Get all the true IDs from all the IDEs in the hit
                std::vector<art::Ptr<sim::AuxDetIDE>> ides = findManyIdes.at(data_i);
 
                for(size_t i = 0; i < ides.size(); i++){
      
                    int id = ides[i]->trackID;
                    if(fRollupUnsavedIds) id = std::abs(id);
                    fDataTrueIds[data_i][id] += ides[i]->energyDeposited;
                }
            }
        }//if CRTData found

        //CRTSimHit
        art::Handle< std::vector<CRTHit>> crtSimHitHandle;
        std::vector<art::Ptr<CRTHit> > crtSimHitList;
        std::map<art::Ptr<CRTHit>, int> simHitPtrMap;
        if (event.getByLabel(fCRTSimHitLabel, crtSimHitHandle)){
            art::fill_ptr_vector(crtSimHitList, crtSimHitHandle);
            art::FindManyP<CRTData> findManyData(crtSimHitHandle, event, fCRTSimHitLabel);
      
            for(size_t hit_i = 0; hit_i < crtSimHitList.size(); hit_i++){
      
                simHitPtrMap[crtSimHitList[hit_i]] = hit_i;
                std::vector<art::Ptr<CRTData>> data = findManyData.at(hit_i);
      
                for(size_t data_i = 0; data_i < data.size(); data_i++){
      
                    int dataID = dataPtrMap[data[data_i]];
      
                    for(auto const& di : fDataTrueIds[dataID]){
                        fSimHitTrueIds[hit_i][di.first] += di.second;
                    }
                }
            }
        }//if CRTSimHit found
      
        //CRTTracks
        art::Handle< std::vector<CRTTrack>> crtTrackHandle;
        std::vector<art::Ptr<CRTTrack> > crtTrackList;
      
        if (event.getByLabel(fCRTTrackLabel, crtTrackHandle)) {
            art::fill_ptr_vector(crtTrackList, crtTrackHandle);
      
            art::FindManyP<CRTHit> findManyHits(crtTrackHandle, event, fCRTTrackLabel);
      
            for(size_t track_i = 0; track_i < crtTrackList.size(); track_i++){
      
                std::vector<art::Ptr<CRTHit>> hits = findManyHits.at(track_i);
      
                for(size_t hit_i = 0; hit_i < hits.size(); hit_i++){
      
                    int hitID = simHitPtrMap[hits[hit_i]];
      
                    for(auto const& hi : fSimHitTrueIds[hitID]){
      
                        fTrackTrueIds[track_i][hi.first] += hi.second;
      
                    }
                }
            }
        }//if CRTrack product found
      
    }//Initialize

    // ---------------------------------------------------------------------------------------    
    // Check that two CRT data products are the same
    bool CRTBackTracker::DataCompare(const CRTData& data1, const CRTData& data2){
    
      if(data1.Mac5() != data2.Mac5()) return false;
      if(data1.TTrig() != data2.TTrig()) return false;
      //if(data1.T1() != data2.T1()) return false;
      //if(data1.ADC() != data2.ADC()) return false;
    
      return true;
    }

    // ----------------------------------------------------------------------------------------    
    // Check that two CRT hits are the same
    bool CRTBackTracker::HitCompare(const CRTHit& hit1, const CRTHit& hit2){
    
      if(hit1.ts1_ns != hit2.ts1_ns) return false;
      if(hit1.plane != hit2.plane) return false;
      if(hit1.x_pos != hit2.x_pos) return false;
      if(hit1.y_pos != hit2.y_pos) return false;
      if(hit1.z_pos != hit2.z_pos) return false;
      if(hit1.x_err != hit2.x_err) return false;
      if(hit1.y_err != hit2.y_err) return false;
      if(hit1.z_err != hit2.z_err) return false;
      if(hit1.tagger != hit2.tagger) return false;
    
      return true;
    }

    // -----------------------------------------------------------------------------------------    
    // Check that two CRT tracks are the same
    bool CRTBackTracker::TrackCompare(const CRTTrack& track1, const CRTTrack& track2){
    
      if(track1.ts1_ns != track2.ts1_ns) return false;
      if(track1.plane1 != track2.plane1) return false;
      if(track1.x1_pos != track2.x1_pos) return false;
      if(track1.y1_pos != track2.y1_pos) return false;
      if(track1.z1_pos != track2.z1_pos) return false;
      if(track1.x1_err != track2.x1_err) return false;
      if(track1.y1_err != track2.y1_err) return false;
      if(track1.z1_err != track2.z1_err) return false;
      if(track1.plane2 != track2.plane2) return false;
      if(track1.x2_pos != track2.x2_pos) return false;
      if(track1.y2_pos != track2.y2_pos) return false;
      if(track1.z2_pos != track2.z2_pos) return false;
      if(track1.x2_err != track2.x2_err) return false;
      if(track1.y2_err != track2.y2_err) return false;
      if(track1.z2_err != track2.z2_err) return false;
    
      return true;
    }

    // --------------------------------------------------------------------------------------------    
    // Get all the true particle IDs that contributed to the CRT data product
    std::vector<int> CRTBackTracker::AllTrueIds(const art::Event& event, const CRTData& data){
    
      std::vector<int> ids;
    
      // Get a handle to the CRT data in the event
      auto crtDataHandle = event.getValidHandle<std::vector<CRTData>>(fCRTDataLabel);
      art::FindManyP<sim::AuxDetIDE> findManyIdes(crtDataHandle, event, fCRTDataLabel);
    
      // Find which one matches the data passed to the function
      int data_i = 0, index = 0;
   
      int nmatch=0; 
      for(auto const& crtData : (*crtDataHandle)){
    
        if(DataCompare(crtData, data)) {
            data_i = index;
            nmatch++;
        }
    
        index++;
      }
      if(nmatch==0) std::cout << "BACKTRACKER ERROR: no matches for procided data product found!" << std::endl;
      if(nmatch>1) std::cout << "BACKTRACKER ERROR: multiple matches for given data product found!" << std::endl;

      // Get all the true IDs from all the IDEs in the hit
      std::vector<art::Ptr<sim::AuxDetIDE>> ides = findManyIdes.at(data_i);
    
      for(size_t i = 0; i < ides.size(); i++){
    
        int id = ides[i]->trackID;
    
        if(fRollupUnsavedIds) id = std::abs(id);
    
        ids.push_back(id);
      }
    
      // Remove any repeated IDs
      std::sort(ids.begin(), ids.end());
      ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
    
      return ids;
    }

    // ---------------------------------------------------------------------------------------------    
    // Get all the true particle IDs that contributed to the sim CRT hit
    std::vector<int> CRTBackTracker::AllTrueIds(const art::Event& event, const CRTHit& hit){
    
        std::vector<int> ids;
      
        // Get a handle to the CRT hits in the event
        art::Handle< std::vector<CRTHit> > crtSimHitHandle;
        art::Handle< std::vector<CRTHit> > crtTrueHitHandle;
 
        // Find which one matches the hit passed to the function
        int hit_i = -1, tmp_hit_i = -1, index = 0;
        bool sim=false, findsim=false, findtrue=false;      
        if(event.getByLabel(fCRTSimHitLabel,crtSimHitHandle)) {
            findsim = true;
            for(auto const& crtHit : (*crtSimHitHandle)){
      
                if(HitCompare(crtHit, hit)) tmp_hit_i = index;
                index++;
            }

            if(tmp_hit_i!=-1) {
                hit_i = tmp_hit_i;
                tmp_hit_i=-1;
                sim=true;
            }
        }

        if(event.getByLabel(fCRTTrueHitLabel,crtTrueHitHandle)) {
            findtrue = true;
            for(auto const& crtHit : (*crtTrueHitHandle)){
                
                if(HitCompare(crtHit, hit)) tmp_hit_i = index;
                index++;
            }
            if(tmp_hit_i!=-1){
                if(hit_i!=-1) {
                    mf::LogError("CRTBackTracker") << "true-sim hit ID ambiguity!";
                    return ids;
                }
                hit_i=tmp_hit_i;
            }
        }
        if(hit_i==-1) {
            mf::LogError("CRTBackTracker") << "no match for passed hit found!";
            return ids;
        }
        if(!findsim && !findtrue) {
            mf::LogError("CRTBackTracker") << "no handle to CRTHits could be found!";
            return ids;
        }

        // Get the crt data associated to that hit and the IDEs associate to the data
        std::vector<art::Ptr<sim::AuxDetIDE>> ides;
      
        // Get all the true IDs from all the IDEs in the hit
        if(sim) {
            art::FindManyP<CRTData> findManyData(crtSimHitHandle, event, fCRTSimHitLabel);
            std::vector<art::Ptr<CRTData>> data;
            data = findManyData.at(hit_i);
            art::FindManyP<sim::AuxDetIDE> findManySimIdes(data, event, fCRTDataLabel);
            for(size_t i = 0; i < data.size(); i++){
      
                ides = findManySimIdes.at(i);
                for(size_t j = 0; j < ides.size(); j++){
      
                    int id = ides[j]->trackID;
                    if(fRollupUnsavedIds) id = std::abs(id);
                    ids.push_back(id);
                }
            }
        }//if sim
        else {
            art::FindManyP<sim::AuxDetIDE> findManyTrueIdes(crtTrueHitHandle, event, fCRTTrueHitLabel);
            ides = findManyTrueIdes.at(hit_i);
            for(auto const& ide : ides) {
                int id = ide->trackID;
                if(fRollupUnsavedIds) id = std::abs(id);
                ids.push_back(id);
            }
        }//else true

        // Remove any repeated IDs
        std::sort(ids.begin(), ids.end());
        ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
      
        return ids;
    }

    //-----------------------------------------------------------------------------------------------    
    // Get all the true particle IDs that contributed to the CRT track
    std::vector<int> CRTBackTracker::AllTrueIds(const art::Event& event, const CRTTrack& track){
    
      std::vector<int> ids;
    
      // Get a handle to the CRT tracks in the event
      auto crtTrackHandle = event.getValidHandle<std::vector<CRTTrack>>(fCRTTrackLabel);
      art::FindManyP<CRTHit> findManyHits(crtTrackHandle, event, fCRTTrackLabel);
    
      // Find which one matches the track passed to the function
      int track_i = 0, index = 0;
    
      for(auto const& crtTrack : (*crtTrackHandle)){
    
        if(TrackCompare(crtTrack, track)) track_i = index;
    
        index++;
    
      }
    
      // Get the crt hits associated to that hit and the data associate to the hits
      std::vector<art::Ptr<CRTHit>> hits = findManyHits.at(track_i);
      art::FindManyP<CRTData> findManyData(hits, event, fCRTSimHitLabel);
    
      // Get all the true IDs from all the IDEs in the track
      for(size_t i = 0; i < hits.size(); i++){
    
        std::vector<art::Ptr<CRTData>> data = findManyData.at(i);
        art::FindManyP<sim::AuxDetIDE> findManyIdes(data, event, fCRTDataLabel);
    
        for(size_t j = 0; j < data.size(); j++){
    
          std::vector<art::Ptr<sim::AuxDetIDE>> ides = findManyIdes.at(j);
    
          for(size_t k = 0; k < ides.size(); k++){
    
            int id = ides[k]->trackID;
    
            if(fRollupUnsavedIds) id = std::abs(id);
    
            ids.push_back(id);
    
          }
        }
      }
    
      // Remove any repeated IDs
      std::sort(ids.begin(), ids.end());
      ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
    
      return ids;
    }

    //------------------------------------------------------------------------------------------    
    // Get the true particle ID that contributed the most energy to the CRT data product
    int CRTBackTracker::TrueIdFromTotalEnergy(const art::Event& event, const CRTData& data){
    
      std::map<int, double> ids;
    
      // Get a handle to the CRT data in the event
      auto crtDataHandle = event.getValidHandle<std::vector<CRTData>>(fCRTDataLabel);
      art::FindManyP<sim::AuxDetIDE> findManyIdes(crtDataHandle, event, fCRTDataLabel);
    
      // Find which one matches the data passed to the function
      int data_i = 0, index = 0;
    
      for(auto const& crtData : (*crtDataHandle)){
    
        if(DataCompare(crtData, data)) data_i = index;
        index++;
      }
    
      // Get all the true IDs from all the IDEs in the hit
      std::vector<art::Ptr<sim::AuxDetIDE>> ides = findManyIdes.at(data_i);
    
      for(size_t i = 0; i < ides.size(); i++){
    
        int id = ides[i]->trackID;
        if(fRollupUnsavedIds) id = std::abs(id);
        ids[id] += ides[i]->energyDeposited;
      }
    
      // Find the true ID that contributed the most energy
      double maxEnergy = -1;
      int trueId = -99999;
    
      for(auto &id : ids){
    
        if(id.second > maxEnergy){
    
          maxEnergy = id.second;
          trueId = id.first;
        }
      }
    
      return trueId;
    }

    //-------------------------------------------------------------------------------------------
    int CRTBackTracker::TrueIdFromDataId(const art::Event& event, int data_i){
    
      if(fDataTrueIds.find(data_i) != fDataTrueIds.end()){ 
    
        double maxEnergy = -1;
        int trueId = -99999;
    
        for(auto &id : fDataTrueIds[data_i]){
    
          if(id.second > maxEnergy){
    
            maxEnergy = id.second;
            trueId = id.first;
          }
        }
    
        return trueId;
      }
    
      return -99999;
    }

    //--------------------------------------------------------------------------------------------
    // Get the true particle ID that contributed the most energy to the CRT hit
    int CRTBackTracker::TrueIdFromTotalEnergy(const art::Event& event, const CRTHit& hit){
    
      std::map<int, double> ids;
    
      // Get a handle to the CRT hits in the event
      auto crtHitHandle = event.getValidHandle<std::vector<CRTHit>>(fCRTSimHitLabel);
      art::FindManyP<CRTData> findManyData(crtHitHandle, event, fCRTSimHitLabel);
    
      // Find which one matches the hit passed to the function
      int hit_i = 0, index = 0;
    
      for(auto const& crtHit : (*crtHitHandle)){
    
        if(HitCompare(crtHit, hit)) hit_i = index;
        index++;
      }
    
      // Get the crt data associated to that hit and the IDEs associate to the data
      std::vector<art::Ptr<CRTData>> data = findManyData.at(hit_i);
      art::FindManyP<sim::AuxDetIDE> findManyIdes(data, event, fCRTDataLabel);
    
      // Get all the true IDs from all the IDEs in the hit
      for(size_t i = 0; i < data.size(); i++){
    
        std::vector<art::Ptr<sim::AuxDetIDE>> ides = findManyIdes.at(i);
    
        for(size_t j = 0; j < ides.size(); j++){
    
          int id = ides[j]->trackID;
          if(fRollupUnsavedIds) id = std::abs(id);
          ids[id] += ides[j]->energyDeposited;
        }
      }
    
      // Find the true ID that contributed the most energy
      double maxEnergy = -1;
      int trueId = -99999;
    
      for(auto &id : ids){
    
        if(id.second > maxEnergy){
    
          maxEnergy = id.second;
          trueId = id.first;
    
        }
      }
    
      return trueId;
    
    }

    //-----------------------------------------------------------------------------------------
    int CRTBackTracker::TrueIdFromHitId(const art::Event& event, int hit_i){
    
      if(fSimHitTrueIds.find(hit_i) != fSimHitTrueIds.end()){ 
    
        double maxEnergy = -1;
        int trueId = -99999;
    
        for(auto &id : fSimHitTrueIds[hit_i]){
    
          if(id.second > maxEnergy){
    
            maxEnergy = id.second;
            trueId = id.first;
          }
        }
    
        return trueId;
    
      }
    
      return -99999;
    
    }

    //--------------------------------------------------------------------------------------------    
    // Get the true particle ID that contributed the most energy to the CRT track
    int CRTBackTracker::TrueIdFromTotalEnergy(const art::Event& event, const CRTTrack& track){
    
      std::map<int, double> ids;
    
      // Get a handle to the CRT tracks in the event
      auto crtTrackHandle = event.getValidHandle<std::vector<CRTTrack>>(fCRTTrackLabel);
      art::FindManyP<CRTHit> findManyHits(crtTrackHandle, event, fCRTTrackLabel);
    
      // Find which one matches the track passed to the function
      int track_i = 0, index = 0;
    
      for(auto const& crtTrack : (*crtTrackHandle)){
    
        if(TrackCompare(crtTrack, track)) track_i = index;
        index++;
      }
    
      // Get the crt hits associated to that hit and the data associate to the hits
      std::vector<art::Ptr<CRTHit>> hits = findManyHits.at(track_i);
      art::FindManyP<CRTData> findManyData(hits, event, fCRTSimHitLabel);
    
      // Get all the true IDs from all the IDEs in the track
      for(size_t i = 0; i < hits.size(); i++){
    
        std::vector<art::Ptr<CRTData>> data = findManyData.at(i);
        art::FindManyP<sim::AuxDetIDE> findManyIdes(data, event, fCRTDataLabel);
    
        for(size_t j = 0; j < data.size(); j++){
    
          std::vector<art::Ptr<sim::AuxDetIDE>> ides = findManyIdes.at(j);
    
          for(size_t k = 0; k < ides.size(); k++){
    
            int id = ides[k]->trackID;
            if(fRollupUnsavedIds) id = std::abs(id);
            ids[id] += ides[k]->energyDeposited;
          }
        }
    
      }
    
      // Find the true ID that contributed the most energy
      double maxEnergy = -1;
      int trueId = -99999;
    
      for(auto &id : ids){
    
        if(id.second > maxEnergy){
    
          maxEnergy = id.second;
          trueId = id.first;
        }
      }
    
      return trueId;
    }

    //--------------------------------------------------------------------------
    int CRTBackTracker::TrueIdFromTrackId(const art::Event& event, int track_i){
    
      if(fTrackTrueIds.find(track_i) != fTrackTrueIds.end()){ 
    
        double maxEnergy = -1;
    
        int trueId = -99999;
    
        for(auto &id : fTrackTrueIds[track_i]){
    
          if(id.second > maxEnergy){
    
            maxEnergy = id.second;
            trueId = id.first;
          }
        }
    
        return trueId;
      }
      return -99999;
    }

 }//namespace crt
}//namespace icarus

