///  \file   MMCluster.hh
///
///  \author Christopher Rogan
///          (crogan@cern.ch)
///
///  \date   2016 Sept
///


#ifndef MMCluster_H
#define MMCluster_H

#include "TObject.h"
#include "MMHit.h"
#include <algorithm>
#include <vector>



class MMCluster : public TObject{
private:
    std::vector<MMHit*> m_hits;
    long s_id=0;
    long m_id;
    long m_track_id = -1;
    int m_pacman = 0;


    
public:
    MMCluster();
    ~MMCluster();
    
    long GetID() ;

    void AddHit( MMHit& hit);
    void AddHitFront( MMHit& hit);
    void SortHits();
    void Epmty();

    double Channel() ;          // Cluster Center of Gravity in "strips"
    double Pos(double pitch) ;     // Cluster Center of Gravity in "mm" (pitch in mm)
    double Charge() ;
    
    int NHoles() ;          // returns the number of empty strips
    int MaxConsecutiveHoles() ;
    int ClusterWidth() ;

    bool Contains( MMHit& hit) ;
     
    int GetNHits() ;
    int GetLeadingCluster(std::vector<MMCluster*> VClusters);

    size_t size() ;
    MMHit  & Get(int i) ;
    MMHit  & operator [] (int i) ;
    
    std::vector<MMHit*> VecHits();

    
    std::vector<MMHit*>::iterator begin();
    std::vector<MMHit*>::iterator end();

    void SetIndexOfFirstPacmanHit(int i);
    int IndexOfFirstPacmanHit() ;
     
    ClassDef(MMCluster, 1)
    
};


// initialize counter
//long MMCluster::s_id = 0;
 MMCluster::MMCluster() : m_id(s_id++) {
}

MMCluster::~MMCluster() {
 for (auto h: m_hits) delete h;
}

 void MMCluster::AddHit(  MMHit & hit) {
  m_hits.push_back(new MMHit(hit));
}

 void MMCluster::AddHitFront( MMHit & hit) {
  m_hits.insert(m_hits.begin(), new MMHit(hit));
}

 void MMCluster::SortHits() {
  std::sort(m_hits.begin(), m_hits.end(),
            [](MMHit* a, MMHit* b){ return a->Strip() < b->Strip(); });
}

void MMCluster::Epmty() {
    m_hits.clear();
//    int Nhit = GetNHits();
//    for (int i = 0; i < Nhit; i++)  delete Get(i);
}

double MMCluster::Channel()  {
 double ch = 0;
 int Nhit = GetNHits();
 for(int i = 0; i < Nhit; i++)
   ch += double(Get(i).Strip())*Get(i).Charge();
 ch /= Charge();
 return ch;
}
double MMCluster::Pos(double pitch)  {
 double position = 0;
 double chnl = Channel();
    position=chnl*pitch;
    
  return position;
}

 double MMCluster::Charge()  {
  double Q = 0;
  int Nhit = GetNHits();
  for(int i = 0; i < Nhit; i++)
    Q += double(Get(i).Charge());

  return Q;
}


 int MMCluster::NHoles()  {
  int Nhit = GetNHits();
  int clus_size = Get(Nhit-1).Channel()-Get(0).Channel()+1;

  return clus_size - Nhit;
}

  int MMCluster::MaxConsecutiveHoles()  {
  int max_hole = -1;
  int Nhit = GetNHits();
  for (int i = 0; i < (Nhit - 1); i++){
    int hole = Get(i+1).Channel() - Get(i).Channel() - 1;
    if (hole > max_hole)
      max_hole = hole;
  }
  return max_hole;
}

 int MMCluster::ClusterWidth()  {
  return Get(GetNHits() - 1).Channel() - Get(0).Channel() + 1;
}

int MMCluster::GetNHits()  {
 return int(m_hits.size());
}

size_t MMCluster::size()  {
 return m_hits.size();
}

 MMHit & MMCluster::Get(int i)  {
 return *m_hits[i];
}

 MMHit & MMCluster::operator [] (int i)  {
 return Get(i);
}

 void MMCluster::SetIndexOfFirstPacmanHit(int i) {
  m_pacman = i;
}

int MMCluster::IndexOfFirstPacmanHit()  {
  return m_pacman;
}


 long MMCluster::GetID()  {
  return m_id;
}

 bool MMCluster::Contains( MMHit& hit)  {
  for (auto test: m_hits)
    if (
  //    test->Layer()  == hit.Layer() &&
  //      test->Radius() == hit.Radius() &&
        test->Strip()  == hit.Strip()  )
      return true;
  return false;
}

std::vector<MMHit*> MMCluster::VecHits()    {
    
    std::vector<MMHit*> hit_list;
    
    for (int i=0; i< m_hits.size(); i++)
    hit_list.push_back(m_hits.at(i));
    
    return hit_list;
}



int GetLeadingCluster(std::vector<MMCluster*> VClusters)  {
          long int leadsumQ=0;
          int index;
        
        for (int i=0; i< VClusters.size(); i++){
                if(VClusters.at(i)->Charge()>leadsumQ) {
                    leadsumQ = VClusters.at(i)->Charge();
                    index=i;
                }
        }
      return index;
}


#endif
