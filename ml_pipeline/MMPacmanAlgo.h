///
///  \file   MMPacmanAlgo.h
///
///  \author Christopher Rogan
///          (crogan@cern.ch)
///
///  \date   2016 Sept
///


#ifndef MMPacmanAlgo_H
#define MMPacmanAlgo_H

//#include "include/Constants.hh"
//#include "include/MMClusterAlgo.hh"
//#include "MMCluster.h"

class MMPacmanAlgo : public TObject {

public:
  // clus_size - 1 is the number of holes allowed
  MMPacmanAlgo(int clus_size = 1,
               int hole_size = 0,
               double seed_thresh = 10.,
               double hit_thresh  = 2.);

  ~MMPacmanAlgo() {}

  std::vector<MMCluster*> Clustering(std::vector<MMHit> & hits);

  void SetClusterSize(int clus_size);
  void SetSeedThreshold(double thresh);
  void SetHitThreshold(double thresh);
  int GetGoodHits();
  void SetGoodHits(int nhit);
    bool IsGoodHit( MMHit& hit);

private:
  int m_good_hits; // number of hits in the last clustering
  int m_clus_size;
  int m_hole_size;
  double m_seed_thresh;
  double m_hit_thresh;

    ClassDef(MMPacmanAlgo, 1)

};

 MMPacmanAlgo::MMPacmanAlgo(int clus_size, int hole_size,
                                  double seed_thresh,
                                  double hit_thresh){ 
  m_clus_size = clus_size;
  m_hole_size = hole_size;
  m_seed_thresh = seed_thresh;
  m_hit_thresh = hit_thresh;
  m_good_hits = 0;
}

 std::vector<MMCluster*> MMPacmanAlgo::Clustering(std::vector<MMHit> & hits){
     
//    cout << " IN Clustering  \n"<< endl;
  std::vector<MMCluster*> cluster_list;
  m_good_hits = 0;
    
    MMHit seed;
    // forward step
    int Nhit = hits.size();
    for(int i = 0; i < Nhit; i++){
      MMHit& hit = hits[i];
      m_good_hits++;
 
        // new cluster if seed above thresh
          seed = hit;
          auto cluster = new MMCluster();
          cluster->AddHit(hit);
          cluster->SetIndexOfFirstPacmanHit(i);
          int last_channel = hit.Channel();
        // look for additional hits forward
        for(int j = i+1; j < Nhit; j++){
          MMHit & hit2 = hits[j];
            if(hit2.Channel() <= last_channel+m_clus_size){   // here check if the distance is less than the allowed holes
                i = j;               // move index so we don't look for seeds in this channel
                cluster->AddHit(hit2);
                last_channel = hit2.Channel();
             }
            else {
            break;
            }
        }
        if(cluster->size()>0)
        //    if(cluster->size()>1)
        cluster_list.push_back(cluster);
      
    }
     
// backward step   Canceled !!!
   
     
     
  return cluster_list;
}

 void MMPacmanAlgo::SetClusterSize(int clus_size){
  if(clus_size >= 1)
    m_clus_size = clus_size;
}

 void MMPacmanAlgo::SetSeedThreshold(double thresh){
  m_seed_thresh = thresh;
}

 void MMPacmanAlgo::SetHitThreshold(double thresh){
  m_hit_thresh = thresh;
}

 int MMPacmanAlgo::GetGoodHits(){
  return m_good_hits;
}

 void MMPacmanAlgo::SetGoodHits(int nhit){
  m_good_hits = nhit;
}

bool MMPacmanAlgo::IsGoodHit( MMHit& hit) {
    // if(!hit.IsChargeCalib())   return false;
    // if(hit.Layer() < 0)      return false;
    //if(hit.Radius() < 0)    return false;
    if(hit.Charge() < 0)      return false;
    //  if((hit.BCID() - hit.TrigBCID()) > m_max_BCID_diff)     return false;
    //if((hit.BCID() - hit.TrigBCID()) < m_min_BCID_diff)     return false;

  return true;
}

#endif

