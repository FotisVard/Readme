//A class for track reconstruction
//Based on chain algorithm
//Amperiadou Dimitra

#include "TObject.h"


#ifndef MMTrackingChain_H
#define MMTrackingChain_H

class MMTrackingChain: public TObject{

private:
  std::vector<MMCluster*>VClusters1;
  std::vector<MMCluster*>VClusters2;
  std::vector<MMCluster*>VClusters3;
  std::vector<MMCluster*>VClusters4;
  float corr_pos1;
  float corr_pos2;
  float corr_pos3;
  float corr_pos4;
  int clpev1;
  int clpev2;
  int clpev3;
  int clpev4;
  bool gettrack;
  float pos1Z= 210.4;
  float pos2Z= 139.677;
  float pos3Z= 70.326;
  float pos4Z = 0.5;
  float positions[4][20];
  float positions3[3][20];
  int cl_index[4][20];
  int nothitted;
  std::vector<float>chi2;
  std::vector<int>clpev={0,0,0,0};
  std::vector<std::vector<MMCluster*>>VClusters={{0},{0},{0},{0}};
  std::vector<float>corr_pos={0,0,0,0};
  std::vector<float>posZ={pos1Z,pos2Z,pos3Z,pos4Z};
  std::vector<float>layerz={pos1Z,pos2Z,pos3Z,pos4Z};

public:
  MMTrackingChain();
  ~MMTrackingChain();
  void SetClusters(std::vector<MMCluster*>VCluster1,
		   std::vector<MMCluster*>VCluster2,
		   std::vector<MMCluster*>VCluster3,
		   std::vector<MMCluster*>VCluster4);
  void SetClustersPerEvent(int clusters1,int clusters2,
			   int clusters3,int clusters4);
  void SetPositionCorrections( float corrpos1,
			       float corrpos2,
			       float corrpos3,
			       float corrpos4);
  std::vector<float> Tracking15(float dch);
  std::vector<float> Tracking3(float dch);
  void SetPositions(MMCluster* VClusters1,
		    MMCluster* VClusters2,
		    MMCluster* VClusters3,
		    MMCluster* VClusters4,
		    float chi2,
		    int tracknum);
  void SetClustersIndexTrack(int i,int j,
			     int k,int l,
			     int tracknum);
  void SetPositions3(MMCluster* VClusters1,
		     MMCluster* VClusters2,
		     MMCluster* VClusters3,
		     float chi2,
		     int tracknum);
  std::vector<float> GetPositionsBestTrack();
  std::vector<float> GetPositionsBestTrack3();
  std::vector<int> GetClustersIndexBestTrack();
  std::vector<float> GetZpositions3();
  int GetNotHitted();
  void SetTrack(bool track);
  bool GetTrack();
  void Clear();
  

  ClassDef(MMTrackingChain, 1);

};

MMTrackingChain:: MMTrackingChain(){
  corr_pos1=0;
  corr_pos2=0;
  corr_pos3=0;
  corr_pos4=0;
  clpev1=0;
  clpev2=0;
  clpev3=0;
  clpev4=0;
  VClusters1={0};
  VClusters2={0};
  VClusters3={0};
  VClusters4={0};
  
  gettrack=false;
  for(int i=0;i<4;i++){
    for(int j=0;j<20;j++){
      positions[i][j]=0;
      
    }

  }
   

}

MMTrackingChain::~MMTrackingChain(){

  std::cout<<"IN DESTRUCTOR"<<endl;
  chi2.clear();


}


void MMTrackingChain::SetClusters(std::vector<MMCluster*>VCluster1,std::vector<MMCluster*>VCluster2,
				  std::vector<MMCluster*>VCluster3,std::vector<MMCluster*>VCluster4){
  VClusters1=VCluster1;
  VClusters2=VCluster2;
  VClusters3=VCluster3;
  VClusters4=VCluster4;
  //for the tracking of triple signal in MM
  std::cout<<"Trying the thing with the clusters"<<std::endl;
  VClusters[0]=VCluster1;
  VClusters[1]=VCluster2;
  VClusters[2]=VCluster3;
  VClusters[3]=VCluster4;

  std::cout<<"Clusters ok"<<std::endl;


}

void MMTrackingChain::SetClustersPerEvent(int clusters1,int clusters2,
					  int clusters3,int clusters4){
  clpev1=clusters1;
  clpev2=clusters2;
  clpev3=clusters3;
  clpev4=clusters4;

  //for the tracking of triple signal in MM
  clpev[0]=clusters1;
  clpev[1]=clusters2;
  clpev[2]=clusters3;
  clpev[3]=clusters4;

   std::cout<<"Trying the thing with the clusters"<<std::endl;


}

void MMTrackingChain:: SetPositionCorrections(float corrpos1,float corrpos2,
					      float corrpos3,float corrpos4){
  corr_pos1=corrpos1;
  corr_pos2=corrpos2;
  corr_pos3=corrpos3;
  corr_pos4=corrpos4;
  //for the tracking of triple signal in MM
  corr_pos[0]=corrpos1;
  corr_pos[1]=corrpos2;
  corr_pos[2]=corrpos3;
  corr_pos[3]=corrpos4;
   std::cout<<"Trying the thing with the clusters"<<std::endl;
  


}
//method for the tracking of event type 15
std::vector<float> MMTrackingChain:: Tracking15(float dx){
  
  float slope12,intrc12,expectx12,slope123,intrc123,expectx123;
  int cluster_in_track=0;
  //float dx=4;
  std::vector<float> layerx12(2);
  std::vector<float> layerx123(3);
  std::vector<float> layerx12error(2,1);
  std::vector<float> layerx123error(3,1);
  std::vector<float> layerx(4);
  std::vector<float> layerxerror(4,1);
  std::vector<float> layerz={pos1Z,pos2Z,pos3Z,pos4Z};
  std::vector<float> layerzerror(4,0);
  std:: vector<float>vecpos(4);
  float chisqr;
  int numoftracks=0;

  cout<<"track boolean:"<<GetTrack()<<endl;

  //clusters in 1st chamber
  for(int i=0;i<clpev1;i++){
    layerx12[0]=VClusters1.at(i)->Channel();
    std::cout<<"channel 1:"<<layerx12[0]<<std::endl;
    cluster_in_track=1;
    
    //clusters in 2nd chamber
    for(int j=0;j<clpev2;j++){

      layerx12[1]=VClusters2.at(j)->Channel();
      //layerx.push_back(VClusters2.at(j)->Channel());
      std::cout<<"channel 2:"<<layerx12[1]<<std::endl;
      cluster_in_track=2;

      //linear fit for the two clusters of 1st and 2nd chamber
      std::cout<<"Fit from 1st and 2nd layer"<<endl;
      TGraph *track12=new TGraph(2,&layerx12[0],&layerz[0]);
      //TGraphErrors *track12=new TGraphErrors(2,&layerx12[0],&layerz[0],&layerx12error[0],&layerzerror[0]);
      TF1 *trackfit12=new TF1("trackfit","pol1",layerx12[0],layerx12[1]);
      track12->Fit(trackfit12);
      slope12=trackfit12->GetParameter(1);
      intrc12=trackfit12->GetParameter(0);
      expectx12=(pos3Z-intrc12)/slope12;
      std::cout<<"expected x12="<<expectx12<<endl;

      layerx123[0]=layerx12[0];
      layerx123[1]=layerx12[1];

      //clusters in 3rd chamber
      std::cout<<"clpev3="<<clpev3<<endl;
      for(int k=0;k<clpev3;k++){
	std::cout<<"channel 3:"<<VClusters3.at(k)->Channel()<<endl;

	//find a cluster on 3rd layer in the interval around the expected x
	if(expectx12-dx<=VClusters3.at(k)->Channel() && VClusters3.at(k)->Channel()<=expectx12+dx){
	  std::cout<<"channel 3:"<<VClusters3.at(k)->Channel()<<std::endl;
	  layerx123[2]=VClusters3.at(k)->Channel();
	  cluster_in_track=3;	  
	}else{
	  continue;
	}
	


	if(cluster_in_track==3){

	  //new fit with 1st,2nd and 3rd cluster of the 3 layers

	  std::cout<<"Fit from 1st,2nd and 3rd cluster"<<endl;
	 	
	  TGraph *track123=new TGraph(3,&layerx123[0],&layerz[0]);
	  //TGraphErrors *track123=new TGraphErrors(3,&layerx123[0],&layerz[0],&layerx123error[0],&layerzerror[0]);
	  TF1 *trackfit123=new TF1("trackfit123","pol1",layerx123[0],layerx123[2]);
	  track123->Fit(trackfit123);
	  slope123=trackfit123->GetParameter(1);
	  intrc123=trackfit123->GetParameter(0);
	  expectx123=(pos4Z-intrc123)/slope123;
	  std::cout<<"expected x123="<<expectx123<<endl;
	  layerx[0]=layerx123[0];
	  layerx[1]=layerx123[1];
	  layerx[2]=layerx123[2];
	  //4th chamber
	  for(int l=0;l<clpev4;l++){
	    std::cout<<"Channel 4:"<<VClusters4.at(l)->Channel()<<endl;
	  
	    if(expectx123-dx<=VClusters4.at(l)->Channel() && VClusters4.at(l)->Channel()<=expectx123+dx){
	      std::cout<<"channel 4:"<<VClusters4.at(l)->Channel()<<std::endl;
	      layerx[3]=VClusters4.at(l)->Channel();
	      cluster_in_track=4;
	    }else{
	      continue;
	    }

	    cout<<"cluster in track="<<cluster_in_track<<endl;
	    
	    if(cluster_in_track==4){//finds a track

	      cout<<"in cluster in track==4"<<endl;

	      numoftracks+=1;
	      cout<<"number of tracks"<<numoftracks<<endl;
	      TGraph* track1234=new TGraph(4,&layerx[0],&layerz[0]);
	      //TGraphErrors *track1234=new TGraphErrors(4,&layerx[0],&layerz[0],&layerxerror[0],&layerzerror[0]);
	      TF1 *trackfit = new TF1("trackfit","pol1",layerx[0],layerx[3]);
	      track1234->Fit(trackfit);
	      chisqr=trackfit->GetChisquare();
	      SetPositions(VClusters1.at(i),VClusters2.at(j),VClusters3.at(k),VClusters4.at(l),chisqr,numoftracks);
	      SetClustersIndexTrack(i,j,k,l,numoftracks);	      
	      SetTrack(true);

	    }//if clusters_in_track=4

	    cluster_in_track=3;
	  }//4th chamber

	}// if clusters_in_track=3

	cluster_in_track=2;
      }//3rd chamber

      cluster_in_track=1;
    }//2nd chamber


    cluster_in_track=0;
  }//1st chamber

  cout<<"track boolean:"<<GetTrack()<<endl;

  if(GetTrack()){

    cout<<"in get track"<<endl;

    vecpos=GetPositionsBestTrack();
    return vecpos;
  }

  return vecpos;

}

//tracking method for the 3 signal events
std::vector<float> MMTrackingChain::Tracking3(float dx){

  float slope12,intrc12,expectx12;
  int cluster_in_track=0;
  //float dx=4;
  std::vector<float> layerx12(2);
  std::vector<float> layerx(3);
  //std::vector<float> layerz={pos1Z,pos2Z,pos3Z,pos4Z};
  std:: vector<float>vecpos(3);
  float chisqr;
  int numoftracks=0;

  
  for(int it=0; it<clpev.size(); it++){

    if(clpev.at(it)==0) nothitted=it;

  }

  std::cout<<"nothitted="<<nothitted<<std::endl;

  //erase the element of the vectors that containts info about the not hitted MM
  //so as to run the for loops in the hitted chambers
  //the 4th position will be found with extrapolation
  clpev.erase(clpev.begin()+nothitted);
  VClusters.erase(VClusters.begin()+nothitted);
  corr_pos.erase(corr_pos.begin()+nothitted);
  layerz.erase(layerz.begin()+nothitted);

  //cluster in 1st hitted chamber
  for(int i=0;i<clpev[0];i++){
    layerx12[0]=VClusters[0].at(i)->Channel();
    std::cout<<"triple tracking Channel 1"<<layerx12[0]<<std::endl;
    cluster_in_track=1;

    //clusters in 2nd hitted chamber
    for(int j=0;j<clpev[1];j++){
      layerx12[1]=VClusters[1].at(j)->Channel();
      std::cout<<" triple tracking channel 2="<<layerx12[1]<<std::endl;
      cluster_in_track=2;

      //linear fit for the two clusters of the two hitted chambers
      std::cout<<"triple tracking Fit from 1st and 2nd layer"<<endl;
      TGraph *track12=new TGraph(2,&layerx12[0],&layerz[0]);
      TF1 *trackfit12=new TF1("trackfit","pol1",layerx[0],layerx[1]);
      track12->Fit(trackfit12);
      slope12=trackfit12->GetParameter(1);
      intrc12=trackfit12->GetParameter(0);
      expectx12=(layerz[2]-intrc12)/slope12;
      std::cout<<"expected x12="<<expectx12<<endl;


      layerx[0]=layerx12[0];
      layerx[1]=layerx12[1];

      //clusters in 3rd hitted chamber
      for(int k=0;k<clpev[2];k++){

	//find a cluster on the 3rd hitted chamber in the interval around dx
	if(expectx12-dx<=VClusters[2].at(k)->Channel() && VClusters[2].at(k)->Channel()<=expectx12+dx){
	  std::cout<<"channel 3:"<<VClusters[2].at(k)->Channel()<<std::endl;
	  layerx[2]=VClusters[2].at(k)->Channel();
	  cluster_in_track=3;	  
	}else{
	  continue;
	}

	std::cout<<"cluster in track="<<cluster_in_track<<std::endl;

	if(cluster_in_track==3){

	  std::cout<<"in cluster in track==3"<<std::endl;
	  numoftracks+=1;
	  TGraph *track123=new TGraph(3,&layerx[0],&layerz[0]);
	  TF1 *trackfit = new TF1("trackfit","pol1",layerx[0],layerx[2]);
	  track123->Fit(trackfit);
	  chisqr=trackfit->GetChisquare();
	  
	  SetPositions3(VClusters[0].at(i),VClusters[1].at(j),VClusters[2].at(k),chisqr,numoftracks);       
	  SetTrack(true);
	     


	}//if cluster in track==3

	cluster_in_track=2;

      }// 3rd chamber

      cluster_in_track=1;

    }// 2nd chamber

    cluster_in_track=0;
  }//1st chamber	
  
  
  cout<<"track boolean:"<<GetTrack()<<endl;

  if(GetTrack()){

    cout<<"in get track"<<endl;

    vecpos=GetPositionsBestTrack3();
    return vecpos;
  }

  return vecpos;

}


void MMTrackingChain::SetPositions(MMCluster* VClusters1,MMCluster* VClusters2,
				   MMCluster* VClusters3,MMCluster* VClusters4,
				   float chisqr,int tracknum){

  chi2.push_back(chisqr);

 

  positions[0][tracknum-1]=VClusters1->Pos(0.250)+corr_pos1;
  positions[1][tracknum-1]=VClusters2->Pos(0.250)+corr_pos2;
  positions[2][tracknum-1]=VClusters3->Pos(0.250)+corr_pos3;
  positions[3][tracknum-1]=VClusters4->Pos(0.250)+corr_pos4;
  

  
  

  // pos[0]=VClusters1->Pos(0.250)+corr_pos1;
  // cout<<"Cluster1 pos:"<<VClusters1->Pos(0.250)<<endl;
  // pos[1]=VClusters2->Pos(0.250)+corr_pos2;
  // cout<<"Cluster2 pos:"<<VClusters2->Pos(0.250)<<endl;
  // pos[2]=VClusters3->Pos(0.250)+corr_pos3;
  // cout<<"Cluster3 pos:"<<VClusters3->Pos(0.250)<<endl;
  // pos[3]=VClusters4->Pos(0.250)+corr_pos4;
  // cout<<"Cluster4 pos:"<<VClusters4->Pos(0.250)<<endl;
			       

}

void MMTrackingChain::SetClustersIndexTrack(int i,int j,int k,int l,int tracknum){
  
  cl_index[0][tracknum-1]=i;
  cl_index[1][tracknum-1]=j;
  cl_index[2][tracknum-1]=k;
  cl_index[3][tracknum-1]=l;



}

void MMTrackingChain::SetPositions3(MMCluster* VClusters1, MMCluster* VClusters2,
			       MMCluster* VClusters3,float chisqr,
			       int tracknum){
  chi2.push_back(chisqr);
  
  positions3[0][tracknum-1]=VClusters1->Pos(0.250)+corr_pos[0];
  positions3[1][tracknum-1]=VClusters2->Pos(0.250)+corr_pos[1];
  positions3[2][tracknum-1]=VClusters3->Pos(0.250)+corr_pos[2];
  
}

std::vector<float> MMTrackingChain::GetPositionsBestTrack(){
  //finds the minimum chi2 from all the tracks that corresponds to the best track in the event


  std::vector<float>pos(4);
  int p=0;
  float min=chi2[p];
  

  for(int i=0;i<chi2.size();i++){
    cout<<"chi2="<<chi2[i]<<endl;
    if(chi2[i]<min) {
      min=chi2[i];
      p=i;
     
    }
    
  }
  cout<<"p="<<p<<endl;

  for(int i=0;i<4;i++){
    for(int j=0;j<5;j++){
      cout<<"positions"<<positions[i][j]<<"\t";


    }
    cout<<"\n"<<endl;

  }
  
 
  cout<<"Cluster1 pos:"<<positions[0][p]<<endl;
  cout<<"Cluster2 pos:"<<positions[1][p]<<endl;
  cout<<"Cluster3 pos:"<<positions[2][p]<<endl;
  cout<<"Cluster4 pos:"<<positions[3][p]<<endl;

  for(int i=0;i<4;i++){
    pos[i]=positions[i][p];

  }

  return pos;

}

std::vector<int> MMTrackingChain::GetClustersIndexBestTrack(){

  std::vector<int>index(4);
  int p=0;
  float min=chi2[p];
  for(int i=0;i<chi2.size();i++){
    if(chi2[i]<min){
      min=chi2[i];
      p=i;
    }

  }
  
  for(int i=0;i<4;i++){
    index[i]=cl_index[i][p];

  }
  
 
  
  return index;
  

}

std::vector<float> MMTrackingChain::GetPositionsBestTrack3(){
  //finds the minimum chi2 from all the tracks that corresponds to the best track in the event


  std::cout<<"chi2 size="<<chi2.size()<<std::endl;

  std::vector<float>pos(3);
  int p=0;
  float min=chi2[p];
  

  for(int i=0;i<chi2.size();i++){
    cout<<"chi2="<<chi2[i]<<endl;
    if(chi2[i]<min) {
      min=chi2[i];
      p=i;
     
    }
    
  }
  cout<<"p="<<p<<endl;

  for(int i=0;i<3;i++){
    for(int j=0;j<5;j++){
      cout<<"positions"<<positions3[i][j]<<"\t";


    }
    cout<<"\n"<<endl;

  }
  
  pos[0]=positions3[0][p];
  cout<<"Cluster1 pos:"<<positions[0][p]<<endl;
  pos[1]=positions3[1][p];
  cout<<"Cluster2 pos:"<<positions[1][p]<<endl;
  pos[2]=positions3[2][p];
  cout<<"Cluster3 pos:"<<positions[2][p]<<endl;
  return pos;

}

std::vector<float> MMTrackingChain:: GetZpositions3(){

  return layerz;


}

int MMTrackingChain::GetNotHitted(){

  return nothitted;


}

void MMTrackingChain:: SetTrack( bool track){
  
  gettrack=track;
  
}
bool MMTrackingChain:: GetTrack(){
  
  return gettrack;

}




#endif

