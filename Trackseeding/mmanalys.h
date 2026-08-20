//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Dec 16 19:23:21 2016 by ROOT version 5.34/35
// from TTree apv_raw/apv_raw
// found on file: run268.root_FDCTplus.root
//////////////////////////////////////////////////////////

#ifndef mmanalys_h
#define mmanalys_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <iostream>
#include <map>


// Header file for the classes stored in the TTree if any.
#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.
using namespace std;


class mmanalys {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   ULong64_t       evt;
   UInt_t          error;
   Int_t           daqTimeSec;
   Int_t           daqTimeMicroSec;
   Int_t           srsTimeStamp;
   UInt_t          srsTrigger;
   vector<unsigned int> *srsFec;
   vector<unsigned int> *srsChip;
   vector<unsigned int> *srsChan;
   vector<string>  *mmChamber;
   vector<int>     *mmLayer;
   vector<char>    *mmReadout;
   vector<int>     *mmStrip;
//   vector<vector<short> > *raw_q;
   vector<short>   *max_q;
   vector<int>     *t_max_q;
   vector<float>   *FD_charge;
   vector<float>   *FD_time;
    

   // List of branches
   TBranch        *b_evt;   //!
   TBranch        *b_error;   //!
   TBranch        *b_daqTimeSec;   //!
   TBranch        *b_daqTimeMicroSec;   //!
   TBranch        *b_srsTimeStamp;   //!
   TBranch        *b_srsTrigger;   //!
   TBranch        *b_srsFec;   //!
   TBranch        *b_srsChip;   //!
   TBranch        *b_srsChan;   //!
   TBranch        *b_mmChamber;   //!
   TBranch        *b_mmLayer;   //!
   TBranch        *b_mmReadout;   //!
   TBranch        *b_mmStrip;   //!
//  TBranch        *b_raw_q;   //!
   TBranch        *b_max_q;   //!
   TBranch        *b_t_max_q;   //!
   TBranch        *b_FD_charge;   //!
   TBranch        *b_FD_time;   //!

   mmanalys(TTree *tree=0);
   virtual ~mmanalys();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
    ///////////////////// ======================================================
    virtual void plot_Raw();
    virtual void plot_Clusters();
    virtual void plot_Residuals();
    virtual void Setstyle();
    virtual void lsm( double *x, double *y, int n , double *a, double *b);


    
 /*
    virtual void sort(vector<float> &x, vector<float> &y, vector<float> &q, vector<float> &yer, vector<MMHit> &a1);
    virtual void remove_holes(vector<MMHit> &a1, map< int, vector<MMHit> > &m, int hole_size);
    virtual float get_cluster_charge(vector<MMHit> &hitlist);
    virtual float get_cluster_length(vector<MMHit> &hitlist);
    virtual float get_cluster_maxq(vector<MMHit> &hitlist);
    virtual void cluster_pos(vector<MMHit> &z1m, float *pos, int chamber , float &pitch );
    virtual void transform(vector<MMHit> &hitlist,int chamber,float pitch);
    virtual void cluster_pos_utpc(vector<MMHit> &z,float drift_gap,float drift_vel,float pitch,float *pos,float *a1gr,float *a0gr,float *siga1,float *siga0,float *chi2);
    virtual void find_longest(map<int, vector<MMHit> > &m, vector<MMHit> &hitlist );
    virtual void find_heaviest(map<int, vector<MMHit> > &m, vector<MMHit> &hitlist);
    virtual void plot();
    virtual void Setstyle();
*/

    /*    virtual void chisquare_chek(float *x, float *y, int z, float *a1gr, float *a0gr, float *chi2);
    virtual void track_fit(float *x, float *y, float *q, int z, float *a1gr, float *a0gr, float *a1, float *a0, float *chi2 );
    virtual void correct(vector<MMHit> &hitlist);
    virtual void circle(float *x,float *y,float *x0,float *y0,float *r);
    virtual void Setstyle();
    virtual void fill_histo_cluster(map<int, vector<MMHit> > &m, TH1F &histo );
    virtual void cluster_pos_utpc(vector<MMHit> &z,float drift_gap,float drift_vel,float pitch,float *pos,float *a1gr,float *a0gr,float *siga1,float *siga0,float *chi2);
    virtual void correct_lorenz(vector<MMHit> &hitlist,float  aL);
    virtual void print_cluster(vector<MMHit> &z );
    virtual void plot();
 */
/////////////==============================================
    
    
};

#endif

#ifdef mmanalys_cxx
mmanalys::mmanalys(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

       TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("run126_FDCTplus.root");

  if (!f || !f->IsOpen()) {
      f = new TFile("run126_FDCTplus.root");
}
      f->GetObject("apv_raw",tree);
   }
   Init(tree);
}

mmanalys::~mmanalys()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t mmanalys::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t mmanalys::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void mmanalys::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   srsFec = 0;
   srsChip = 0;
   srsChan = 0;
   mmChamber = 0;
   mmLayer = 0;
   mmReadout = 0;
   mmStrip = 0;
//   raw_q = 0;
   max_q = 0;
   t_max_q = 0;
   FD_charge = 0;
   FD_time = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("evt", &evt, &b_evt);
   fChain->SetBranchAddress("error", &error, &b_error);
   fChain->SetBranchAddress("daqTimeSec", &daqTimeSec, &b_daqTimeSec);
   fChain->SetBranchAddress("daqTimeMicroSec", &daqTimeMicroSec, &b_daqTimeMicroSec);
   fChain->SetBranchAddress("srsTimeStamp", &srsTimeStamp, &b_srsTimeStamp);
   fChain->SetBranchAddress("srsTrigger", &srsTrigger, &b_srsTrigger);
   fChain->SetBranchAddress("srsFec", &srsFec, &b_srsFec);
   fChain->SetBranchAddress("srsChip", &srsChip, &b_srsChip);
   fChain->SetBranchAddress("srsChan", &srsChan, &b_srsChan);
   fChain->SetBranchAddress("mmChamber", &mmChamber, &b_mmChamber);
   fChain->SetBranchAddress("mmLayer", &mmLayer, &b_mmLayer);
   fChain->SetBranchAddress("mmReadout", &mmReadout, &b_mmReadout);
   fChain->SetBranchAddress("mmStrip", &mmStrip, &b_mmStrip);
//   fChain->SetBranchAddress("raw_q", &raw_q, &b_raw_q);
   fChain->SetBranchAddress("max_q", &max_q, &b_max_q);
   fChain->SetBranchAddress("t_max_q", &t_max_q, &b_t_max_q);
   fChain->SetBranchAddress("FD_charge", &FD_charge, &b_FD_charge);
   fChain->SetBranchAddress("FD_time", &FD_time, &b_FD_time);
   Notify();
}

Bool_t mmanalys::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void mmanalys::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t mmanalys::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef mmanalys_cxx
