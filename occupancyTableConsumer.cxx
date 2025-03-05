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
  
  using MyBCTable    = soa::Join<aod::BCsWithTimestamps,aod::OccIndexTable>;
  using MyCollisions = aod::Collisions; 
  using MyTracks     = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA>;
  using MyTracksQA = aod::TracksQA_002;
  using MyOccsTable  = soa::Join<aod::OccsBCsList, aod::OccsDet, aod::OccsTrackMult, aod::OccsMultExtra, aod::OccsRobust,
                                                  aod::OccsMeanDet, aod::OccsMeanTrkMult, aod::OccsMnMultExtra, aod::OccsMeanRobust>;
  using MyTrackMeanOccs = soa::Join<aod::TrackMeanOccs0, aod::TrackMeanOccs1, aod::TrackMeanOccs4,
                                                         aod::TrackMeanOccs5, aod::TrackMeanOccs8>;
  int dfCount = 0;
  int trackCount = 0;
  void process(o2::aod::Origins const& Origins,
               MyBCTable const& BCs,
               MyCollisions const& collisions,
               MyTracks const& tracks,
               MyTracksQA const& tracksQA,
               MyOccsTable const& occTables,
               MyTrackMeanOccs const& trackMeanOccs
              )
  {  
    dfCount++;
    // if(dfCount > 10) {return;}
    LOG(info)<<"DEBUG :: df_"<<dfCount<<" :: DF_"<<Origins.iteratorAt(0).dataframeID();
    LOG(info)<<"DEBUG :: df_"<<dfCount<<" :: BCs.size()           = "<<BCs.size();
    LOG(info)<<"DEBUG :: df_"<<dfCount<<" :: occTables.size()     = "<<occTables.size();
    LOG(info)<<"DEBUG :: df_"<<dfCount<<" :: collisions.size()    = "<<collisions.size();
    LOG(info)<<"DEBUG :: df_"<<dfCount<<" :: tracks.size()        = "<<tracks.size();
    LOG(info)<<"DEBUG :: df_"<<dfCount<<" :: tracksQA.size()      = "<<tracksQA.size();
    LOG(info)<<"DEBUG :: df_"<<dfCount<<" :: trackMeanOccs.size() = "<<trackMeanOccs.size();

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
        LOG(info)<<"DEBUG 1:: "<<trackCount<<" :: "<<"myOccTrack.trackId() = "<<myOccTrack.trackId()<<" :: track.globalIndex() = "<<track.globalIndex()
                 <<" :: occ1 "<<occ1
                 <<" :: occ2 "<<occ2
                 <<" :: occ3 "<<occ3
                 <<" :: occ4 "<<occ4;
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

    for(auto const& track : tracks){

      int occTableIndex = findTrackInList(track.globalIndex(), indexList); //Step 03- Check if track is in index list or not 
      if( occTableIndex == -1){ continue; }
      else { myOccTrack = trackMeanOccs.iteratorAt(occTableIndex); } //Step 04- Get object containing mean occupancies for this track

      float pt = track.pt();
      float eta = track.eta();
      float occ1 = myOccTrack.meanOccPrimUnfm80();
      float occ2 = myOccTrack.meanOccRobustT0V0PrimUnfm80();
      float occ3 = myOccTrack.weightMeanOccPrimUnfm80();
      float occ4 = myOccTrack.weightMeanOccRobustT0V0PrimUnfm80();

      if(trackCount2 < 10){
        trackCount2++;
        LOG(info)<<"DEBUG 2:: "<<trackCount2<<" :: track.globalIndex() = "<<track.globalIndex()
                                           <<" :: myOccTrack.trackId() = "<<myOccTrack.trackId();
      }
    }
  }//Process function ends

};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc){
  return WorkflowSpec{ 
                      adaptAnalysisTask<occupancyTableConsumer>(cfgc)
  };
}
