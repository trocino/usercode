#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"


const UInt_t maxmc  = 50;
const UInt_t maxln  = 50;
const UInt_t maxen  = 50;
const UInt_t maxjn  = 50;
const UInt_t maxmet = 50;
const UInt_t maxg   = 50;

// Tree variables
Int_t cat;
Int_t mccat;
Int_t ngenITpu;
Int_t nvtx;
Int_t nmcparticles;
Float_t l1_px;
Float_t l1_py;
Float_t l1_pz;
Float_t l1_en;
Float_t l1_ptErr;
Int_t   l1_id;
Int_t   l1_pid;
Float_t l2_px;
Float_t l2_py;
Float_t l2_pz;
Float_t l2_en;
Float_t l2_ptErr;
Int_t   l2_id;
Int_t   l2_pid;
Int_t   ln;
Float_t ln_px[maxln];
Float_t ln_py[maxln];
Float_t ln_pz[maxln];
Float_t ln_en[maxln];
Float_t ln_ptErr[maxln];
Int_t   ln_id[maxln];
Float_t en_corren[maxen];
Int_t   jnum;
Float_t jn_px[maxjn];
Float_t jn_py[maxjn];
Float_t jn_pz[maxjn];
Float_t jn_en[maxjn];
Float_t jn_genpt[maxjn];
Float_t jn_btag1[maxjn];
Float_t jn_btag2[maxjn];
Bool_t  jn_tightId[maxjn];
Float_t htvec_px;
Float_t htvec_py;
Float_t met_pt[maxmet];
Float_t met_phi[maxmet];
Float_t mc_px[maxmc];
Float_t mc_py[maxmc];
Float_t mc_pz[maxmc];
Float_t mc_en[maxmc];
Int_t   mc_id[maxmc];
Int_t   gn;
Float_t g_px[maxg];
Float_t g_py[maxg];
Float_t g_pz[maxg];
Float_t g_en[maxg];
Float_t g_corren[maxg];
Float_t g_correnerr[maxg];
Float_t g_iso1[maxg];
Float_t g_iso2[maxg];
Float_t g_iso3[maxg];
Float_t g_sihih[maxg];
Float_t g_r9[maxg];
Float_t g_trkVeto[maxg];
Float_t g_conv[maxg];

// Initialize variables
void initializeTreeVariables() {

  cat = 0;
  mccat = 0;
  ngenITpu = 0;
  nvtx = 0;
  nmcparticles = 0;
  l1_px = 0.;
  l1_py = 0.;
  l1_pz = 0.;
  l1_en = 0.;
  l1_ptErr = 0.;
  l1_id = 0;
  l1_pid = -1;
  l2_px = 0.;
  l2_py = 0.;
  l2_pz = 0.;
  l2_en = 0.;
  l2_ptErr = 0.;
  l2_id = 0;
  l2_pid = -1;
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
    jn_genpt[i] = 0.;
    jn_btag1[i] = 0.;
    jn_btag2[i] = 0.;
    jn_tightId[i] = false;
  }
  for(unsigned int i=0; i<maxen; ++i) {
    en_corren[i] = 0.;
  }
  htvec_px = 0.;
  htvec_py = 0.;
  for(unsigned int i=0; i<maxmet; ++i) {
    met_pt[i] = 0.;
    met_phi[i] = 0.;
  }
  for(unsigned int i=0; i<maxmc; ++i) {
    mc_px[i] = 0.;
    mc_py[i] = 0.;
    mc_pz[i] = 0.;
    mc_en[i] = 0.;
    mc_id[i] = 0;
  }
  gn = 0;
  for(unsigned int i=0; i<maxg; ++i) {
    g_px[i] = 0.;
    g_py[i] = 0.;
    g_pz[i] = 0.;
    g_en[i] = 0.;
    g_corren[i] = 0.;
    g_correnerr[i] = 0.;
    g_iso1[i] = 0.;
    g_iso2[i] = 0.;
    g_iso3[i] = 0.;
    g_sihih[i] = 0.;
    g_r9[i] = 0.;
    g_trkVeto[i] = 0;
    g_conv[i] = 0;
  }
}

// Set branch addresses
void attachToTree(TChain *t) {

  t->SetBranchAddress("cat", &cat);
  t->SetBranchAddress("mccat", &mccat);
  t->SetBranchAddress("ngenITpu", &ngenITpu);
  t->SetBranchAddress("nvtx", &nvtx);
  t->SetBranchAddress("nmcparticles", &nmcparticles);
  t->SetBranchAddress("l1_px", &l1_px);
  t->SetBranchAddress("l1_py", &l1_py);
  t->SetBranchAddress("l1_pz", &l1_pz);
  t->SetBranchAddress("l1_en", &l1_en);
  t->SetBranchAddress("l1_ptErr", &l1_ptErr);
  t->SetBranchAddress("l1_id", &l1_id);
  t->SetBranchAddress("l1_pid", &l1_pid);
  t->SetBranchAddress("l2_px", &l2_px);
  t->SetBranchAddress("l2_py", &l2_py);
  t->SetBranchAddress("l2_pz", &l2_pz);
  t->SetBranchAddress("l2_en", &l2_en);
  t->SetBranchAddress("l2_ptErr", &l2_ptErr);
  t->SetBranchAddress("l2_id", &l2_id);
  t->SetBranchAddress("l2_pid", &l2_pid);
  t->SetBranchAddress("ln", &ln);
  t->SetBranchAddress("ln_px", &ln_px);
  t->SetBranchAddress("ln_py", &ln_py);
  t->SetBranchAddress("ln_pz", &ln_pz);
  t->SetBranchAddress("ln_en", &ln_en);
  t->SetBranchAddress("ln_ptErr", &ln_ptErr);
  t->SetBranchAddress("ln_id", &ln_id);
  t->SetBranchAddress("en_corren", &en_corren);
  t->SetBranchAddress("jn", &jnum);
  t->SetBranchAddress("jn_px", &jn_px);
  t->SetBranchAddress("jn_py", &jn_py);
  t->SetBranchAddress("jn_pz", &jn_pz);
  t->SetBranchAddress("jn_en", &jn_en);
  t->SetBranchAddress("jn_genpt", &jn_genpt);
  t->SetBranchAddress("jn_btag1", &jn_btag1);
  t->SetBranchAddress("jn_btag2", &jn_btag2);
  t->SetBranchAddress("jn_tightId", &jn_tightId);
  t->SetBranchAddress("htvec_px", &htvec_px);
  t->SetBranchAddress("htvec_py", &htvec_py);
  t->SetBranchAddress("met_pt", &met_pt);
  t->SetBranchAddress("met_phi", &met_phi);
  t->SetBranchAddress("mc_px", &mc_px);
  t->SetBranchAddress("mc_py", &mc_py);
  t->SetBranchAddress("mc_pz", &mc_pz);
  t->SetBranchAddress("mc_en", &mc_en);
  t->SetBranchAddress("mc_id", &mc_id);
  t->SetBranchAddress("gn", &gn);
  t->SetBranchAddress("g_px", &g_px);
  t->SetBranchAddress("g_py", &g_py);
  t->SetBranchAddress("g_pz", &g_pz);
  t->SetBranchAddress("g_en", &g_en);
  t->SetBranchAddress("g_corren", &g_corren);
  t->SetBranchAddress("g_correnerr", &g_correnerr);
  t->SetBranchAddress("g_iso1", &g_iso1);
  t->SetBranchAddress("g_iso2", &g_iso2);
  t->SetBranchAddress("g_iso3", &g_iso3);
  t->SetBranchAddress("g_sihih", &g_sihih);
  t->SetBranchAddress("g_r9", &g_r9);
  t->SetBranchAddress("g_trkVeto", &g_trkVeto);
  t->SetBranchAddress("g_conv", &g_conv);

}

