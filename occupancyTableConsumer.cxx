// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
///
/// \brief  Service Work Code for Alice Data
/// \author Rahul Verma (rahul.verma@iitb.ac.in) :: Marian I Ivanov (mivanov@cern.ch)

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"

#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/EventSelection.h"

#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "DerivedPIDTables.h"
#include "Common/DataModel/OccupancyTables.h"

#include <typeinfo>

#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsFT0/Digit.h"
#include "DataFormatsParameters/GRPLHCIFData.h"
// #include "MetadataHelper.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

int32_t nBCsPerOrbit = o2::constants::lhc::LHCMaxBunches;

struct occupancyTableConsumer{

  template<typename T>
  int64_t findTrackInList(const int64_t& target, const T& sortedArray){
    auto it = std::lower_bound(sortedArray.begin(), sortedArray.end(), target);
    if (it != sortedArray.end() && *it == target) {
      int index = std::distance(sortedArray.begin(), it);
      return index;
    } else {
      return -1; 
    }
  }
  
  // using MyBCTable = soa::Join<aod::BCsWithTimestamps, aod::OccIndexTable, aod::BCTFinfoTable>;
  using MyBCTable    = soa::Join<aod::BCsWithTimestamps, aod::BCTFinfoTable >;
  using MyCollisions = aod::Collisions; 
  using MyTracks     = soa::Join<aod::Tracks, aod::TrackToTmo, aod::TrackToTracksQA, aod::TracksExtra, aod::TracksDCA>;
  using MyTracksQA   = soa::Join<aod::TracksQAVersion, aod::TrackQAToTmo>;
  // using MyOccsTable  = soa::Join<aod::OccsBCsList, aod::OccsDet, aod::OccsTrackMult, aod::OccsMultExtra, aod::OccsRobust,
  //                                                 aod::OccsMeanDet, aod::OccsMeanTrkMult, aod::OccsMnMultExtra, aod::OccsMeanRobust>;
  // using MyTrackMeanOccs = soa::Join<aod::TrackMeanOccs0, aod::TrackMeanOccs1, aod::TrackMeanOccs4,
  //                                                        aod::TrackMeanOccs5, aod::TrackMeanOccs8>;

  using MyTrackMeanOccs = soa::Join<aod::TmoTrackIds, aod::TmoToTrackQA, aod::TmoPrim,  aod::TmoT0V0,  aod::TmoRT0V0Prim, aod::TwmoRT0V0Prim> ; 
  // using MyTmoTable = aod::TmoTrackId;

  int dfCount = 0;
  int trackCount = 0;
  int bcCount = 0 ;
  int checkCounter = 0;

  void process(o2::aod::Origins const& Origins
               ,MyBCTable const& BCs
               ,MyCollisions const& collisions
               ,MyTrackMeanOccs const& tmoTable
               ,MyTracks const& tracks
               ,MyTracksQA const& tracksQA /*,
              //  ,MyOccsTable const& occTables,
              //  ,MyTrackMeanOccs const& trackMeanOccs*/
              )
  { 
    auto bc = collisions.begin().bc_as<MyBCTable>();
    LOG(info)<<"DEBUG :: df_"<<dfCount<<" :: bc = "<<bc.tfId();
    dfCount++;

    //TMO to other table dereferencing (aod::tracks, aod::TracksQAVersion)
    checkCounter=0;
    for( const auto& trackMeanOccs : tmoTable){
      checkCounter++;
      if(checkCounter < 10){
        LOG(info)<<"DEBUG :: trackMeanOccs.globalIndex() = "<<trackMeanOccs.globalIndex()<<" :: trackMeanOccs.trackId() = "<<trackMeanOccs.trackId()<<" :: trackMeanOccs.trackQAId() = "<<trackMeanOccs.trackQAId();
      }
      // if(trackMeanOccs.trackId() != -1){
        auto track_From_TrackMeanOccTable = trackMeanOccs.track_as<MyTracks>();//Dereferencing TrackMeanOccupancy ==> aod::Tracks
        if(trackMeanOccs.trackId() !=  track_From_TrackMeanOccTable.globalIndex()) {
          LOG(info)<<"DEBUG :: ERROR ERROR ERROR :: trackMeanOccs.trackId() !=  track_From_TrackMeanOccTable.globalIndex() :: "<<trackMeanOccs.trackId()<<" !=  "<<track_From_TrackMeanOccTable.globalIndex();
        }
      // }
      // if(trackMeanOccs.trackQAId() != -1){
        auto trackQA_From_TrackMeanOccTable = trackMeanOccs.trackQA_as<MyTracksQA>();//Dereferencing TrackMeanOccupancy ==> aod::TracksQAVersion
        if(trackMeanOccs.trackId() != trackQA_From_TrackMeanOccTable.trackId()){
          LOG(info)<<"DEBUG :: ERROR ERROR ERROR :: trackMeanOccs.trackId() !=  trackQA_From_TrackMeanOccTable.trackId() :: "<<trackMeanOccs.trackId()<<" !=  "<<trackQA_From_TrackMeanOccTable.trackId(); //this table is bad
        }
      // }
    }

    //track to other table dereferencing (aod::TracksQAVersion, aod::TmoTrackId)
    checkCounter = 0;
    for( const auto& track : tracks){
      checkCounter++;
      if(checkCounter < 10){
        LOG(info)<<"DEBUG :: track.globalIndex() = "<<track.globalIndex()<<" :: track.tmoId() = "<<track.tmoId()<<" :: trackMeanOccs.trackQAId() = "<<track.trackQAId();
      }
      if(track.tmoId() != -1){
        auto tmo_From_Tracks = track.tmo_as<MyTrackMeanOccs>();//Dereferencing aod::Tracks ==> TrackMeanOccupancy
        if(track.globalIndex() !=  tmo_From_Tracks.trackId()) {
          LOG(info)<<"DEBUG :: ERROR ERROR ERROR :: track.globalIndex() !=  tmo_From_Tracks.trackId() :: "<<track.globalIndex()<<" !=  "<<tmo_From_Tracks.trackId(); //this is also failing
        }
      }
      if(track.trackQAId() != -1){
        auto trackQA_From_Tracks = track.trackQA_as<MyTracksQA>();//Dereferencing aod::Tracks ==> aod::TracksQAVersion
        if(track.globalIndex() != trackQA_From_Tracks.trackId()){
          LOG(info)<<"DEBUG :: ERROR ERROR ERROR :: track.globalIndex() != trackQA_From_Tracks.trackId() :: "<<track.globalIndex()<<" !=  "<<trackQA_From_Tracks.trackId(); //this is also failing 
        }
      }
    }

    //trackQA to other table dereferencing (aod::Tracks, aod::TmoTrackId)
    checkCounter = 0;
    for( const auto& trackQA : tracksQA){
      checkCounter++;
      if(checkCounter < 10){
        LOG(info)<<"DEBUG :: trackQA.globalIndex = "<<trackQA.globalIndex()<<" :: trackQA.trackId = "<<trackQA.trackId()<<" :::: trackQA.tmoId = "<<trackQA.tmoId();
      }
        if(trackQA.tmoId() != -1){
          auto tmo_From_TracksQA = trackQA.tmo_as<MyTrackMeanOccs>();//Dereferencing aod::TracksQAVersion ==> TrackMeanOccupancy
          if(trackQA.trackId() !=  tmo_From_TracksQA.trackId()) {
            LOG(info)<<"DEBUG :: ERROR ERROR ERROR :: trackQA.trackId() !=  tmo_From_TracksQA.trackId() :: "<<trackQA.trackId()<<" !=  "<<tmo_From_TracksQA.trackId();
          }
        }
        if(trackQA.trackId() != -1){
          auto track_From_TracksQA = trackQA.track_as<MyTracks>();//Dereferencing aod::TracksQAVersion ==> aod::Tracks
          if(trackQA.trackId() != track_From_TracksQA.globalIndex()){
            LOG(info)<<"DEBUG :: ERROR ERROR ERROR :: trackQA.trackId() != track_From_TracksQA.globalIndex() :: "<<trackQA.trackId()<<" !=  "<<track_From_TracksQA.globalIndex();
          }
        }
    }
/*    // if(dfCount > 10) {return;}
    LOG(info)<<"DEBUG :: df_"<<dfCount<<" :: DF_"<<Origins.iteratorAt(0).dataframeID();
    LOG(info)<<"DEBUG :: df_"<<dfCount<<" :: BCs.size()           = "<<BCs.size();
    LOG(info)<<"DEBUG :: df_"<<dfCount<<" :: occTables.size()     = "<<occTables.size();
    LOG(info)<<"DEBUG :: df_"<<dfCount<<" :: collisions.size()    = "<<collisions.size();
    LOG(info)<<"DEBUG :: df_"<<dfCount<<" :: tracks.size()        = "<<tracks.size();
    LOG(info)<<"DEBUG :: df_"<<dfCount<<" :: tracksQA.size()      = "<<tracksQA.size();
    LOG(info)<<"DEBUG :: df_"<<dfCount<<" :: trackMeanOccs.size() = "<<trackMeanOccs.size();

    for(auto const& bc : BCs){
      if( bcCount % 10000 == 0){
      LOG(info)<<"DEBUG :: df_"<<dfCount<<" :: bc.gloablIndex() = "<<bc.globalIndex()
                                        <<" :: bc.tfId() = "<<bc.tfId()
                                        <<" :: bc.bcInTF() = "<<bc.bcInTF();
      }
      bcCount++;
    }

    for(auto const& myOccTrack : trackMeanOccs){
      trackCount++;
      auto track = myOccTrack.track_as<MyTracks>(); // de-reference to actual track
      float pt = track.pt();
      float eta = track.eta();
      float occ1 = myOccTrack.meanOccPrimUnfm80();
      float occ2 = myOccTrack.meanOccRobustT0V0PrimUnfm80();
      float occ3 = myOccTrack.weightMeanOccPrimUnfm80();
      float occ4 = myOccTrack.weightMeanOccRobustT0V0PrimUnfm80();

      if(trackCount < 10){
        // LOG(info)<<"DEBUG 1:: "<<trackCount<<" :: "<<"myOccTrack.trackId() = "<<myOccTrack.trackId()<<" :: track.globalIndex() = "<<track.globalIndex()
        //          <<" :: occ1 "<<occ1
        //          <<" :: occ2 "<<occ2
        //          <<" :: occ3 "<<occ3
        //          <<" :: occ4 "<<occ4;
      }
    }

    //Step 01-Build Index List ;
    std::vector<int64_t> indexList;
    for(auto const& myOccTrack : trackMeanOccs){
      indexList.push_back(myOccTrack.trackId());
    }

    int trackCount2 = 0;    
    //Step 02-create a reusable iterator object
    auto myOccTrack = trackMeanOccs.begin();
    
    float occ1 = -9999999;
    float occ2 = -9999999;
    float occ3 = -9999999;
    float occ4 = -9999999;

    int trackWithoutCollision = 0;
    int trackWithCollisionButNoOccs = 0;

    for(auto const& track : tracks){
      int occTableIndex = findTrackInList(track.globalIndex(), indexList); //Step 03- Check if track is in index list or not 
      if( occTableIndex == -1){
        occ1 = -9999999;
        occ2 = -9999999;
        occ3 = -9999999;
        occ4 = -9999999;                              
      }
      else { myOccTrack = trackMeanOccs.iteratorAt(occTableIndex); 
        occ1 = myOccTrack.meanOccPrimUnfm80();
        occ2 = myOccTrack.meanOccRobustT0V0PrimUnfm80();
        occ3 = myOccTrack.weightMeanOccPrimUnfm80();
        occ4 = myOccTrack.weightMeanOccRobustT0V0PrimUnfm80();   
      } //Step 04- Get object containing mean occupancies for this track

      float pt = track.pt();
      float eta = track.eta();

      if(trackCount2 < 10){
        trackCount2++;
        // LOG(info)<<"DEBUG 2:: "<<trackCount2<<" :: track.globalIndex() = "<<track.globalIndex()
        //                                    <<" :: myOccTrack.trackId() = "<<myOccTrack.trackId();
      }

      if(track.collisionId() < 0){
        trackWithoutCollision++;
        // LOG(info)<<"DEBUG 3:: occTableIndex = "<<occTableIndex<<" :: track.collisionId() = "<<track.collisionId()<<" :: track.globalIndex() = "<<track.globalIndex();
      }
      if(track.collisionId() >= 0 && occTableIndex == -1){trackWithCollisionButNoOccs++;}
        
    }
    LOG(info)<<"DEBUG :: trackWithoutCollision       = "<<trackWithoutCollision;
    LOG(info)<<"DEBUG :: trackWithCollisionButNoOccs = "<<trackWithCollisionButNoOccs;
    LOG(info)<<"DEBUG :: total Tracks = "<<tracks.size()<<" :: "<<(trackWithCollisionButNoOccs+trackWithoutCollision+trackMeanOccs.size());
    LOG(info)<<"DEBUG :: ";
*/
  }//Process function ends

};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc){
  return WorkflowSpec{ 
                      adaptAnalysisTask<occupancyTableConsumer>(cfgc)
  };
}
