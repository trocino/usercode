#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"


const UInt_t maxln  = 50;
const UInt_t maxjn  = 50;
const UInt_t maxmet = 50;

// Tree variables
Int_t cat;
Int_t ngenITpu;
Float_t l1_px;
Float_t l1_py;
Float_t l1_pz;
Float_t l1_en;
Float_t l1_ptErr;
Int_t   l1_id;
Float_t l2_px;
Float_t l2_py;
Float_t l2_pz;
Float_t l2_en;
Float_t l2_ptErr;
Int_t   l2_id;
Int_t   ln;
Float_t ln_px[maxln];
Float_t ln_py[maxln];
Float_t ln_pz[maxln];
Float_t ln_en[maxln];
Float_t ln_ptErr[maxln];
Int_t   ln_id[maxln];
Int_t   jnum;
Float_t jn_px[maxjn];
Float_t jn_py[maxjn];
Float_t jn_pz[maxjn];
Float_t jn_en[maxjn];
Float_t htvec_px;
Float_t htvec_py;
Float_t met_pt[maxmet];
Float_t met_phi[maxmet];

// Initialize variables
void initializeTreeVariables() {

  cat = 0;
  ngenITpu = 0;
  l1_px = 0.;
  l1_py = 0.;
  l1_pz = 0.;
  l1_en = 0.;
  l1_ptErr = 0.;
  l1_id = 0;
  l2_px = 0.;
  l2_py = 0.;
  l2_pz = 0.;
  l2_en = 0.;
  l2_ptErr = 0.;
  l2_id = 0;
  ln = 0;
  for(unsigned int i=0; i<maxln; ++i) {
    ln_px[i] = 0.;
    ln_py[i] = 0.;
    ln_pz[i] = 0.;
    ln_en[i] = 0.;
    ln_ptErr[i] = 0.;
    ln_id[i] = 0;
  }
  jnum = 0;
  for(unsigned int i=0; i<maxln; ++i) {
    jn_px[i] = 0.;
    jn_py[i] = 0.;
    jn_pz[i] = 0.;
    jn_en[i] = 0.;
  }
  htvec_px = 0.;
  htvec_py = 0.;
  for(unsigned int i=0; i<maxmet; ++i) {
    met_pt[i] = 0.;
    met_phi[i] = 0.;
  }
}

// Set branch addresses
void attachToTree(TChain *t) {

  t->SetBranchAddress("cat", &cat);
  t->SetBranchAddress("ngenITpu", &ngenITpu);
  t->SetBranchAddress("l1_px", &l1_px);
  t->SetBranchAddress("l1_py", &l1_py);
  t->SetBranchAddress("l1_pz", &l1_pz);
  t->SetBranchAddress("l1_en", &l1_en);
  t->SetBranchAddress("l1_ptErr", &l1_ptErr);
  t->SetBranchAddress("l1_id", &l1_id);
  t->SetBranchAddress("l2_px", &l2_px);
  t->SetBranchAddress("l2_py", &l2_py);
  t->SetBranchAddress("l2_pz", &l2_pz);
  t->SetBranchAddress("l2_en", &l2_en);
  t->SetBranchAddress("l2_ptErr", &l2_ptErr);
  t->SetBranchAddress("l2_id", &l2_id);
  t->SetBranchAddress("ln", &ln);
  t->SetBranchAddress("ln_px", &ln_px);
  t->SetBranchAddress("ln_py", &ln_py);
  t->SetBranchAddress("ln_pz", &ln_pz);
  t->SetBranchAddress("ln_en", &ln_en);
  t->SetBranchAddress("ln_ptErr", &ln_ptErr);
  t->SetBranchAddress("ln_id", &ln_id);
  t->SetBranchAddress("jn", &jnum);
  t->SetBranchAddress("jn_px", &jn_px);
  t->SetBranchAddress("jn_py", &jn_py);
  t->SetBranchAddress("jn_pz", &jn_pz);
  t->SetBranchAddress("jn_en", &jn_en);
  t->SetBranchAddress("htvec_px", &htvec_px);
  t->SetBranchAddress("htvec_py", &htvec_py);
  t->SetBranchAddress("met_pt", &met_pt);
  t->SetBranchAddress("met_phi", &met_phi);

}

