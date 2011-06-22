/** \class fitLLOut
 *  Root macro for the contruction of a likelihood ratio discriminant for the
 *  selection of ZZ->llvv events and WW background
 *  
 *  The macro reads the discrimination variables from a TTree (or TChain)
 *  both for signal and background events
 *
 *  To run the macro 
 *  root [0] .x fitLLOut.r
 * 
 *  $Date: 2008/03/14 15:41:47 $
 *  $Revision: 1.2 $
 *  \author G. Cerminara - NEU Boston & INFN Torino
 *
 * NOTE: need the following settings (ROOT v. 5.28 for RooFit classes)
 *  source /afs/cern.ch/sw/lcg/external/gcc/4.3.2/x86_64-slc5/setup.csh
 *  source /afs/cern.ch/sw/lcg/app/releases/ROOT/5.28.00c/x86_64-slc5-gcc43-opt/root/bin/thisroot.csh
 *
 */

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TNtuple.h"
#include "TString.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
// #include "TMultiLayerPerceptron.h"
#include "TCanvas.h"
#include "TMLPAnalyzer.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TGraph.h"
#include "macros.C"
//#include "plotSignBckg.r"
#include "TGaxis.h"

#include "THStack.h"
#include "TEventList.h"

#include "RooFit.h"

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooArgSet.h"
#include "RooGaussian.h"
#include "RooBreitWigner.h"
#include "RooFitResult.h"
#include "RooGlobalFunc.h"
#include "RooGaussModel.h"
#include "RooVoigtian.h"
#include "RooPolynomial.h"
#include "RooAddPdf.h"
#include "RooDataHist.h"
#include "RooExponential.h"
#include "RooNLLVar.h"
#include "RooExtendPdf.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooKeysPdf.h"
#include "RooHistPdf.h"

#include "TLimitDataSource.h"
#include "TConfidenceLevel.h"
#include "TLimit.h"

//#include "main_distrib_funct.C"

#include <iostream>
#include "extra_root_lib/EnsembleProducer.h"
#include "extra_root_lib/TreeSum.h"
#include "extra_root_lib/Utils.h"
//#include "extra_root_lib/LumiNormalization.h"
#include "extra_root_lib/LumiNorm.h"
#include "extra_root_lib/TreeCombiner.h"
#include "extra_root_lib/RndmGenerator.h"
#include "extra_root_lib/LikelihoodDisc.h"
#endif

using namespace std;
using namespace RooFit;

#define protected public 

// -------------------------------------------------------------------
// Some configuration variables
int form=2;
TString epoch="run2a";
bool dqApplied=true;


float massMinCut=76.1876;
float massMaxCut=106.1876;
double lumiFactor=1.;

TString inputDir="/tmp/trocino/AllData";
TString inputFileName="ZZllvvAnalyzer_whatsample.root";
bool doNormalizeLuminosity=false;
bool doNormalizeXSection=false;

// Produce ensemble plots
bool doEnsemble=false;
bool doSigPlusBck=true;
bool doBckOnly=false;

// Plot all the ensemble plots
bool drawEnsemblePlots=false;

bool doToyMc=false;
int nEvToyMc=1000;
bool extended=false;
bool cutOnLikelihood=false;
double lOutCut=0.5;

int nBins=10;
// -------------------------------------------------------------------
bool debug=false;

bool doData=false;
bool doL11=false;

double expo(double *x, double *par); // Implementation of exponential function for TF1 class
void fitLLOut();
void normalize(TH1F *histo);
void normalize(TH2F *histo);
//double poisson(double *x, double *par);
//double f_poisson(double n, double mu, double norm);
RooDataSet* buildDataSet(const TString &name, const TString &title,
			 TNtuple *ntuple, const LikelihoodDisc &disc, const TString &finalState);
//RooDataSet mergeDiMuAndDiEm(const RooDataSet& dimuDataSet, const RooDataSet& diemDataSet);
RooDataSet buildToyMCDataSample(const TString &finalState, const int nEv, const RooAbsPdf &pdf);

double function(double *x, double *par);


void fitLLOut() {
  // Load needed macros
  gROOT->LoadMacro("macros.C");
  //gROOT->LoadMacro("plotSignBckg.r");
  //gROOT->LoadMacro("main_distrib_funct.C");
//   gSystem->Load("libRooFit");
//   gSystem->Load("extra_root_lib/libEvent.so");
//   gSystem->AddLinkedLibs("extra_root_lib/libEvent.so");

  // Get the style
//   gROOT->LoadMacro("/afs/cern.ch/user/t/trocino/Utilities/danystyle.C");
  //setTDRStyle();
  //TStyle * style=getStyle("tdr");
  // Style options
  //style->SetOptStat("RMEN");
  //style->SetOptStat("E");
  //style->SetOptFit(101);
  //style->SetOptStat(1111111);
  //style->SetOptTitle(1);
  //style->cd();       

  // ------------------------------------------------------------------------------
  vector<TString> signalSample;
  signalSample.push_back("zz_llnn");

  /*
  vector<TString> bckSampleDiem;
  bckSampleDiem.push_back("ww");
  bckSampleDiem.push_back("tt");
  bckSampleDiem.push_back("wz");
  bckSampleDiem.push_back("zz_llll");
  bckSampleDiem.push_back("z_15-60");
  bckSampleDiem.push_back("z_60-130");
  bckSampleDiem.push_back("z_130-250");
  bckSampleDiem.push_back("ztt_15-60");
  bckSampleDiem.push_back("ztt_60-130");
  bckSampleDiem.push_back("ztt_130-250");
  bckSampleDiem.push_back("wgamma_rad");
  bckSampleDiem.push_back("wgamma_prod");
  bckSampleDiem.push_back("wjets_w0lp");
  bckSampleDiem.push_back("wjets_w1lp");
  bckSampleDiem.push_back("wjets_w2lp");
  bckSampleDiem.push_back("wjets_w3lp");
  bckSampleDiem.push_back("wjets_w4lp");
  bckSampleDiem.push_back("wjets_w5lp");
  */

  vector<TString> bckSampleDimu;
  //bckSampleDimu.push_back("wz");
  bckSampleDimu.push_back("ww");
  bckSampleDimu.push_back("tt");
  bckSampleDimu.push_back("t");
  bckSampleDimu.push_back("zz_x");
  bckSampleDimu.push_back("dy_mm");
  // FIXME: add W+jets

  /*
  TString XSectionFile=inputDir+"XsectionNew.txt";
  // Get the normalization for the two final states
  LumiNormalization lumiNormDiem(XSectionFile, "./Luminosity.txt", epoch, "diem", dqApplied, inputDir);
  lumiNormDiem.setAdditionalScale(lumiFactor);
  lumiNormDiem.addData();
  lumiNormDiem.addMC(signalSample);
  lumiNormDiem.addMC(bckSampleDiem);
  */

  /*
  LumiNormalization lumiNormDimu(XSectionFile, "./Luminosity.txt", epoch, "dimu", dqApplied, inputDir);
  lumiNormDimu.setAdditionalScale(lumiFactor);
  lumiNormDimu.addData();
  lumiNormDimu.addMC(signalSample);
  lumiNormDimu.addMC(bckSampleDimu);
  */
  LumiNorm lumiNormDimu(inputDir, inputFileName, doNormalizeLuminosity, doNormalizeXSection);
  lumiNormDimu.setAdditionalScale(lumiFactor);
  lumiNormDimu.addData();
  lumiNormDimu.addMC(signalSample);
  lumiNormDimu.addMC(bckSampleDimu);

  // ------------------------------------------------------------------------------
  // Build signal and background ntuples with correct weight
  // Get the signal ntuple (with correct weight)
  /*
  TreeCombiner comboSignalDiem(&lumiNormDiem, "diem", inputDir);
  comboSignalDiem.setMassCut(massMinCut, massMaxCut);
  TNtuple *sigNtuple_diem=comboSignalDiem.getCombinedTree(signalSample);
  */

  TreeCombiner comboSignalDimu(&lumiNormDimu, inputDir, inputFileName);
  comboSignalDimu.setMassCut(massMinCut, massMaxCut);
  comboSignalDimu.requireFullSection(true);
  TNtuple *sigNtuple_dimu=comboSignalDimu.getCombinedTree(signalSample);

  // Get the backgorund ntuple (with correct weight)
  /*
  TreeCombiner comboBckDiem(&lumiNormDiem, "diem", inputDir);
  comboBckDiem.setMassCut(massMinCut, massMaxCut);
  TNtuple *bckNtuple_diem=comboBckDiem.getCombinedTree(bckSampleDiem);
  */

  TreeCombiner comboBckDimu(&lumiNormDimu, inputDir, inputFileName);
  comboBckDimu.setMassCut(massMinCut, massMaxCut);
  comboBckDimu.requireFullSection(true);
  TNtuple *bckNtuple_dimu=comboBckDimu.getCombinedTree(bckSampleDimu);

  double massMin=massMinCut;
  double massMax=massMaxCut;
  int massBins=(int)(massMax-massMin)/3;

  // ------------------------------------------------------------------------------
  // Build the likelihood object and train it
  /*
  LikelihoodDisc discDiem(sigNtuple_diem,bckNtuple_diem);
  discDiem.separateSamples(true); // Use different samples for training and test
  discDiem.setMassCut(massMinCut,massMaxCut); // Set the mass cut
  discDiem.addVariable("L0", "diLeptInvMass", "M_{ll} (GeV)", massBins, massMin, massMax, "pol6");

  pair<TH1F *, TF1 *> lrL0=discDiem.getLR("L0");
  newCanvas(lrL0.first->GetName(),form);
  TF1 *myFunc=new TF1("L1Fit",function,50,130,13);
  TF1 *polA=new TF1("polA","pol4",50,82);
  TF1 *g1=new TF1("g1","gaus",82,102);
  TF1 *polB=new TF1("polB","pol4",102,130);
  //TF1 *polC=new TF1("polB","pol3",110,130);
  cout << "----- polB ------------------------------------------" << endl;
  lrL0.first->Fit(polB,"R");
  cout << "----- polA ------------------------------------------" << endl;
  lrL0.first->Fit(polA,"R+");
  cout << "----- g1 ------------------------------------------" << endl;
  lrL0.first->Fit(g1,"R+");

 
  Double_t par[13];
  polA->GetParameters(&par[0]);
  g1->GetParameters(&par[5]);
  polB->GetParameters(&par[8]);
  //polC->GetParameters(&par[11]);

  myFunc->SetParameters(par);


  discDiem.addVariable("L1", "diLeptInvMass", "M_{ll} (GeV)", massBins, massMin, massMax, myFunc);
  double LeadLeptPt_bins_diem[18]={20,29,38,47,56,65,74,83,92,101,110,119,128,137,146,155,173,200};
  discDiem.addVariable("L2", "LeadLeptPt", "Leading lepton pT (GeV)",
		       17,20,200, LeadLeptPt_bins_diem,"pol4");

  discDiem.addVariable("L3", "DeltaPhiLeadLeptZ","#Delta#phi lead l-Z", 24,0,1.4,"pol4");
  discDiem.addVariable("L4","leptMinusCmCosTheta","Cos(#theta^{*}) lept -", 25, -1, 1, "pol4");
  discDiem.addVariable("L5","TransMass4Body","4 body Transv Mass",25,40,400,"pol5");
  discDiem.addVariable("L6","MinMass4Body","Min 4 body Inv Mass",25,100,400,"pol6");
  if(doL11) {
    double MinvTot_bins_diem[25]={80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,330,400};
    discDiem.addVariable("L11","MinvTot","Minv Tot (dilept + MET)",24,80,400,MinvTot_bins_diem, "pol4");
  }
  */
  
  LikelihoodDisc discDimu(sigNtuple_dimu, bckNtuple_dimu);
  discDimu.separateSamples(true);
  discDimu.setMassCut(massMinCut, massMaxCut);
  discDimu.useMassCut(false);

  // ---- Build fitting function for dilepton mass ----
  discDimu.addVariable("L0", "diLeptInvMass", "M_{ll} [GeV/c^{2}]", massBins, massMin, massMax, "pol4");
  pair<TH1D *, TF1 *> lrL0=discDimu.getLR("L0");
  newCanvas(lrL0.first->GetName(), form);
  TF1 *myFunc=new TF1("L1Fit", function, massMin, massMax, 13);
  TF1 *polA=new TF1("polA", "pol2", massMin, 86.);
  TF1 *g1=new TF1("g1", "gaus", 86., 96.);
  TF1 *polB=new TF1("polB", "pol2", 96., massMax);
  cout << "----- polB ------------------------------------------" << endl;
  lrL0.first->Fit(polB, "R");
  cout << "----- polA ------------------------------------------" << endl;
  lrL0.first->Fit(polA, "R+");
  cout << "----- g1 ------------------------------------------" << endl;
  lrL0.first->Fit(g1, "R+"); 
  Double_t par[13]={0.};
  polA->GetParameters(&par[0]);
  g1->GetParameters(&par[5]);
  polB->GetParameters(&par[8]);
  myFunc->SetParameters(par);
  // --------------------------------------------------

  discDimu.addVariable("L1", "diLeptInvMass", "M_{ll} [GeV/c^{2}]", massBins, massMin, massMax, myFunc);
  discDimu.addVariable("L2", "leadPt", "Leading lepton p_{T} [GeV/c^{2}]", 10, 20, 200, "pol4");
  discDimu.addVariable("L3", "dileptLeadDeltaPhi", "#Delta#phi lead l-Z", 10, 0.2, 1., "pol4");
  discDimu.addVariable("L4", "leptMinusCmCosTheta", "Cos(#theta^{*}) lept -", 10, -1., 1., "pol4");

  /*
  discDimu.addVariable("L1", "diLeptInvMass", "M_{ll} [GeV/c^{2}]", massBins, massMin, massMax, "pol4");
  double LeadLeptPt_bins_dimu[18]={20, 29, 38, 47, 56, 65, 74, 83, 92, 101, 110, 119, 128, 137, 146, 155, 173, 200};
  discDimu.addVariable("L2", "leadPt", "Leading lepton p_{T} [GeV/c^{2}]",
		       17, 20, 200, LeadLeptPt_bins_dimu, "pol4");
  discDimu.addVariable("L3", "dileptLeadDeltaPhi", "#Delta#phi lead l-Z", 24, 0, 1.4, "pol4");
  discDimu.addVariable("L4", "leptMinusCmCosTheta", "Cos(#theta^{*}) lept -", 25, -1, 1, "pol4");
  */

  //discDimu.addVariable("L5", "TransMass4Body", "4 body Transv Mass", 25, 40, 400, "pol5");
  //discDimu.addVariable("L6", "MinMass4Body", "Min 4 body Inv Mass", 25, 100, 400, "pol6");
  //if(doL11) {
  //  double MinvTot_bins_dimu[24]={80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 330};
  //  discDimu.addVariable("L11", "MinvTot", "Minv Tot (dilept + MET)", 23, 80, 330, MinvTot_bins_dimu, "pol4");
  //}


  // ------------------------------------------------------------------------------
  // Plot the likelihood ratio for the different components
  /*
  pair<TH1F *, TF1 *> lrL1_diem=discDiem.getLR("L1");
  newCanvas(TString(lrL1_diem.first->GetName())+"_diem",form);
  lrL1_diem.first->Draw();
  
  pair<TH1F *, TF1 *> lrL2_diem=discDiem.getLR("L2");
  newCanvas(TString(lrL2_diem.first->GetName())+"_diem",form);
  lrL2_diem.first->Draw();

  pair<TH1F *, TF1 *> lrL3_diem=discDiem.getLR("L3");
  newCanvas(TString(lrL3_diem.first->GetName())+"_diem",form);
  lrL3_diem.first->Draw();

  pair<TH1F *, TF1 *> lrL4_diem=discDiem.getLR("L4");
  newCanvas(TString(lrL4_diem.first->GetName())+"_diem",form);
  lrL4_diem.first->Draw();

  pair<TH1F *, TF1 *> lrL5_diem=discDiem.getLR("L5");
  newCanvas(TString(lrL5_diem.first->GetName())+"_diem",form);
  lrL5_diem.first->Draw();

  pair<TH1F *, TF1 *> lrL6_diem=discDiem.getLR("L6");
  newCanvas(TString(lrL6_diem.first->GetName())+"_diem",form);
  lrL6_diem.first->Draw();

  if(doL11) {
    pair<TH1F *, TF1 *> lrL11_diem=discDiem.getLR("L11");
    newCanvas(TString(lrL11_diem.first->GetName())+"_diem",form);
    lrL11_diem.first->Draw();
  }
  */

  pair<TH1D *, TF1 *> lrL1_dimu=discDimu.getLR("L1");
  if( !(lrL1_dimu.first && lrL1_dimu.second) ) {
    cout << "[fitLLOut] *** Error: in lrL1_dimu, TH1D * = " << lrL1_dimu.first
	 << " and TF1 * = " << lrL1_dimu.second << endl;
    throw std::exception();
  }
  newCanvas(TString(lrL1_dimu.first->GetName())+"_dimu",form);
  //   TCanvas *ppp1=new TCanvas(TString(lrL1_dimu.first->GetName())+"_dimu", 
  // 			    TString(lrL1_dimu.first->GetName())+"_dimu", 
  // 			    600, 600);
  //   ppp1->cd();
  lrL1_dimu.first->Draw();

  pair<TH1D *, TF1 *> lrL2_dimu=discDimu.getLR("L2");
  if( !(lrL2_dimu.first && lrL2_dimu.second) ) {
    cout << "[fitLLOut] *** Error: in lrL2_dimu, TH1D * = " << lrL2_dimu.first
	 << " and TF1 * = " << lrL2_dimu.second << endl;
    throw std::exception();
  }
  newCanvas(TString(lrL2_dimu.first->GetName())+"_dimu",form);
  //   TCanvas *ppp2=new TCanvas(TString(lrL2_dimu.first->GetName())+"_dimu",
  // 			    TString(lrL2_dimu.first->GetName())+"_dimu",
  // 			    600, 600);
  //   ppp2->cd();
  lrL2_dimu.first->Draw();

  pair<TH1D *, TF1 *> lrL3_dimu=discDimu.getLR("L3");
  if( !(lrL3_dimu.first && lrL3_dimu.second) ) {
    cout << "[fitLLOut] *** Error: in lrL3_dimu, TH1D * = " << lrL3_dimu.first
	 << " and TF1 * = " << lrL3_dimu.second << endl;
    throw std::exception();
  }
  newCanvas(TString(lrL3_dimu.first->GetName())+"_dimu",form);
  //   TCanvas *ppp3=new TCanvas(TString(lrL3_dimu.first->GetName())+"_dimu",
  // 			    TString(lrL3_dimu.first->GetName())+"_dimu",
  // 			    600, 600);
  //   ppp3->cd();
  lrL3_dimu.first->Draw();

  pair<TH1D *, TF1 *> lrL4_dimu=discDimu.getLR("L4");
  if( !(lrL4_dimu.first && lrL4_dimu.second) ) {
    cout << "[fitLLOut] *** Error: in lrL4_dimu, TH1D * = " << lrL4_dimu.first
	 << " and TF1 * = " << lrL4_dimu.second << endl;
    throw std::exception();
  }
  newCanvas(TString(lrL4_dimu.first->GetName())+"_dimu",form);
  //   TCanvas *ppp4=new TCanvas(TString(lrL4_dimu.first->GetName())+"_dimu",
  // 			    TString(lrL4_dimu.first->GetName())+"_dimu",
  // 			    600, 600);
  //   ppp4->cd();
  lrL4_dimu.first->Draw();

  //pair<TH1F *, TF1 *> lrL5_dimu=discDimu.getLR("L5");
  //newCanvas(TString(lrL5_dimu.first->GetName())+"_dimu",form);
  //lrL5_dimu.first->Draw();

  //pair<TH1F *, TF1 *> lrL6_dimu=discDimu.getLR("L6");
  //newCanvas(TString(lrL6_dimu.first->GetName())+"_dimu",form);
  //lrL6_dimu.first->Draw();

  //if(doL11) {
  //  pair<TH1F *, TF1 *> lrL11_dimu=discDimu.getLR("L11");
  //  newCanvas(TString(lrL11_dimu.first->GetName())+"_dimu",form);
  //  lrL11_dimu.first->Draw();
  //}


  // ------------------------------------------------------------------------------
//   // Instantiate the variables to be read from the tree
//   float diLeptInvMass, LeadLeptPt, NoLeadLeptPt,
//     DeltaPhiLeadLeptZ, leptMinusCmCosTheta,
//     DiLeptDeltaEta, TransMass4Body, MinMass4Body, weight;
//   float run, event, type; // must be float to cope with TNtuple syntax


  // Open signal and background trees
  /*
  int sigEntries_diem=sigNtuple_diem->GetEntries();
  int bckEntries_diem=bckNtuple_diem->GetEntries();
  cout << "Signal TTree # of entries:" << endl
       << "    - diem: " << sigEntries_diem << " entries" << endl
       << "    - dimu: " << sigEntries_dimu << " entries" << endl;
  */
  int sigEntries_dimu=sigNtuple_dimu->GetEntries();
  int bckEntries_dimu=bckNtuple_dimu->GetEntries();
  cout << "Signal TTree # of entries:" << endl
    //<< "    - diem: " << sigEntries_diem << " entries" << endl
       << "    - dimu: " << sigEntries_dimu << " entries" << endl;
  cout << "Background TTree # of entries:" << endl
    //<< "    - diem: " << bckEntries_diem << " entries" << endl
       << "    - dimu: " << bckEntries_dimu << " entries" << endl;


  // ------------------------------------------------------------------------------
  // Build the histos of the likelihood output
  double binMin=0;
  double binMax=1;

  // Choose the variables to be used
  vector<TString> L1xL2xL3xL4;
  L1xL2xL3xL4.push_back("L1"); L1xL2xL3xL4.push_back("L2");
  L1xL2xL3xL4.push_back("L3"); L1xL2xL3xL4.push_back("L4");
  //if(doL11) L1xL2xL3xL4.push_back("L11");

  /*
  //// Diem
  // Signal:
  TH1F* hTestStatSigDiem=discDiem.evaluate(L1xL2xL3xL4,sigNtuple_diem, "SigDiem", nBins);
  // Background
  TH1F* hTestStatBckDiem=discDiem.evaluate(L1xL2xL3xL4,bckNtuple_diem, "BckDiem", nBins);
  */

  //// Dimu
  // Signal:
  TH1F* hTestStatSigDimu=discDimu.evaluate(L1xL2xL3xL4, sigNtuple_dimu, "SigDimu", nBins);
  // Background
  TH1F* hTestStatBckDimu=discDimu.evaluate(L1xL2xL3xL4, bckNtuple_dimu, "BckDimu", nBins);


  // ------------------------------------------------------------------------------
  
  // ------------------------------------------------------------------------------
  // Get the number of signal and background samples
  /*
  double nSigEvAfterPreselDiem=hTestStatSigDiem->Integral();
  double nBckEvAfterPreselDiem=hTestStatBckDiem->Integral();
  */
  double nSigEvAfterPreselDimu=hTestStatSigDimu->Integral();
  double nBckEvAfterPreselDimu=hTestStatBckDimu->Integral();
  // ------------------------------------------------------------------------------

  // Create the TNtuple of sig+bck with correct weight
  TFile *file=new TFile("/tmp/trocino/fitLLOut.root", "RECREATE");
  file->cd();
  /*
  vector<TString> sampleListDiem;
  copy(signalSample.begin(), signalSample.end(), back_inserter(sampleListDiem));
  copy(bckSampleDiem.begin(), bckSampleDiem.end(), back_inserter(sampleListDiem));

  TreeCombiner comboDiem(&lumiNormDiem, "diem", inputDir);
  comboDiem.setMassCut(massMinCut, massMaxCut);
  file->cd();
  TNtuple *allMcTreeDiem=comboDiem.getCombinedTree(sampleListDiem);

  file->cd();
  allMcTreeDiem->Write();

  // ------------------------------------------------------------------------------
  // Plot the original distributions
  plot1DHistos(hTestStatSigDiem,hTestStatBckDiem,1,false,true);
  */
  
  vector<TString> sampleListDimu;
  copy(signalSample.begin(), signalSample.end(), back_inserter(sampleListDimu));
  copy(bckSampleDimu.begin(), bckSampleDimu.end(), back_inserter(sampleListDimu));
  
  TreeCombiner comboDimu(&lumiNormDimu, inputDir, inputFileName);
  comboDimu.setMassCut(massMinCut, massMaxCut);
  comboDimu.requireFullSection(true);
  file->cd();
  TNtuple *allMcTreeDimu=comboDimu.getCombinedTree(sampleListDimu);

  file->cd();
  allMcTreeDimu->Write();
  
  // ------------------------------------------------------------------------------
  // Plot the original distributions
  //plot1DHistos(hTestStatSigDimu, hTestStatBckDimu, 1, false, true);

  // Plot the sum of the signal and background distributions
  /*
  newCanvas("hStack_diem", form);
  THStack *hsTestStatDiem=new THStack("hsTestStatDiem","Likelihood Out - diem -");
  hsTestStatDiem->SetMinimum(0.001);
  //TLegend* l=new TLegend(0.8,0.6,0.95,0.97);
  //l->SetMargin(0.5);
  hTestStatSigDiem->SetFillColor(8);
  hTestStatSigDiem->SetFillStyle(3002);
  hsTestStatDiem->Add(hTestStatSigDiem,"hist");
  hTestStatBckDiem->SetFillColor(50);
  hTestStatBckDiem->SetFillStyle(3002);
  hsTestStatDiem->Add(hTestStatBckDiem,"hist");
  hsTestStatDiem->Draw("E");
  */

  newCanvas("hStack_dimu",form);
  //   TCanvas *tmpCanv=new TCanvas("hStack_dimu", "hStack_dimu", 600, 600);
  //   tmpCanv->cd();
  THStack * hsTestStatDimu=new THStack("hsTestStatDimu", "Likelihood Out - dimu -");
  hsTestStatDimu->SetMinimum(0.001);
  //TLegend* l=new TLegend(0.8,0.6,0.95,0.97);
  //l->SetMargin(0.5);
  hTestStatSigDimu->SetFillColor(8);
  hTestStatSigDimu->SetFillStyle(3002);
  hsTestStatDimu->Add(hTestStatSigDimu,"hist");
  hTestStatBckDimu->SetFillColor(50);
  hTestStatBckDimu->SetFillStyle(3002);
  hsTestStatDimu->Add(hTestStatBckDimu,"hist");
  hsTestStatDimu->Draw("E");

  // ------------------------------------------------------------------------------
  // Define the variables to be used in the fit
  RooRealVar likeLihoodOut("LikeLihoodOut", "LL out", binMin, binMax);
  RooRealVar wgh("weight","weight", 0, 100); // event weight
  // Define the set of variables: mass and weight
  //RooArgSet setOfVar(likeLihoodOut,wgh);

  // -------------------------------------------------------------------------------------
  // Get the histos and build the binned samples
  file->cd();

  /*
  RooDataSet signalMC_diem=*buildDataSet("signalMC_diem", "Signal MC dataset", sigNtuple_diem, discDiem, "diem");
  // Un-binned sample (from back tree)
  RooDataSet bckMC_diem=*buildDataSet("bckMC_diem", "Bck MC dataset", bckNtuple_diem, discDiem, "diem");
  // Un-binned sample signal + background
  RooDataSet allMC_diem=*buildDataSet("allMC_diem", "Sig+Bck MC dataset", allMcTreeDiem, discDiem, "diem");
  */

  RooDataSet signalMC_dimu=*buildDataSet("signalMC_dimu", "Signal MC dataset", sigNtuple_dimu, discDimu, "dimu");
  // Un-binned sample (from back tree)
  RooDataSet bckMC_dimu=*buildDataSet("bckMC_dimu", "Bck MC dataset", bckNtuple_dimu, discDimu, "dimu");
  // Un-binned sample signal + background
  RooDataSet allMC_dimu=*buildDataSet("allMC_dimu", "Sig+Bck MC dataset", allMcTreeDimu, discDimu, "dimu");

  RooCategory tagCat("FinalStateCat","Final state category");
  /*tagCat.defineType("diem");*/
  tagCat.defineType("dimu");

  RooArgSet setOfVar(likeLihoodOut,wgh,tagCat);
  RooDataSet allMC_allFinalState("AllFinalState", "AllFinalState", setOfVar, WeightVar(wgh));
  //RooDataSet allMC_allFinalState("AllFinalState", "AllFinalState", setOfVar, "weight"); // OLD
  //allMC_allFinalState.setWeightVar(wgh);                                                // OLD
  /*allMC_allFinalState.append(allMC_diem);*/
  allMC_allFinalState.append(allMC_dimu);


  // -------------------------------------------------------------------------------------
  // Plot the datasets
  /*
  RooPlot *plotSig_diem=likeLihoodOut.frame(binMin, binMax, nBins);
  signalMC_diem.plotOn(plotSig_diem,DataError(RooAbsData::SumW2));
  signalMC_diem.statOn(plotSig_diem);
  
  RooDataHist dataHistSigDiem("dataHistSigDiem","dataHistSigDiem",likeLihoodOut,hTestStatSigDiem);
  RooHistPdf pdfSig_diem("hPdf_sigDiem","hPdf_sigDiem",likeLihoodOut,dataHistSigDiem);

  RooDataHist dataHistBckDiem("dataHistBckDiem","dataHistBckDiem",likeLihoodOut,hTestStatBckDiem);
  RooHistPdf pdfBck_diem("hPdf_bckDiem","hPdf_bckDiem",likeLihoodOut,dataHistBckDiem);
  */

  RooPlot *plotSig_dimu=likeLihoodOut.frame(binMin, binMax, nBins);
  signalMC_dimu.plotOn(plotSig_dimu,DataError(RooAbsData::SumW2));
  signalMC_dimu.statOn(plotSig_dimu);

  RooDataHist dataHistSigDimu("dataHistSigDimu","dataHistSigDimu",likeLihoodOut,hTestStatSigDimu);
  RooHistPdf pdfSig_dimu("pdfSig_dimu","pdfSig_dimu",likeLihoodOut,dataHistSigDimu);

  RooDataHist dataHistBckDimu("dataHistBckDimu","dataHistBckDimu",likeLihoodOut,hTestStatBckDimu);
  RooHistPdf pdfBck_dimu("pdfBck_dimu","pdfBck_dimu",likeLihoodOut,dataHistBckDimu);

  //RooKeysPdf pdfSig_diem("pdfSig_diem","pdfSig_diem",likeLihoodOut,signalMC_diem);
  //RooKeysPdf pdfSig_dimu("pdfSig_dimu","pdfSig_dimu",likeLihoodOut,signalMC_dimu);
  //RooKeysPdf pdfBck_diem("pdfBck_diem","pdfBck_diem",likeLihoodOut,bckMC_diem);
  //RooKeysPdf pdfBck_dimu("pdfBck_dimu","pdfBck_dimu",likeLihoodOut,bckMC_dimu);
  //return;

  // -------------------------------------------------------------------------------------
  // Fit the likelihood output - diem
  /*
  pdfSig_diem.plotOn(plotSig_diem,LineColor(kBlue));
  // Draw the results
  newCanvas("Signal_diem");  
  plotSig_diem->Draw();
  */

  pdfSig_dimu.plotOn(plotSig_dimu,LineColor(kBlue));
  // Draw the results
  newCanvas("Signal_dimu");
  //   TCanvas *nnn=new TCanvas("Signal_dimu", "Signal_dimu", 600, 600);
  //   nnn->cd();
  plotSig_dimu->Draw();

  // ---------------------------------------------------------------------------------
  // Fit the background shape
  /*
  RooPlot *plotBck_diem=likeLihoodOut.frame(binMin, binMax, nBins);
  bckMC_diem.plotOn(plotBck_diem,DataError(RooAbsData::SumW2));
  bckMC_diem.statOn(plotBck_diem);

  pdfBck_diem.plotOn(plotBck_diem,LineColor(kGreen));
  // Draw the results
  newCanvas("Background_diem");  
  plotBck_diem->Draw();
  */

  RooPlot *plotBck_dimu=likeLihoodOut.frame(binMin, binMax, nBins);
  bckMC_dimu.plotOn(plotBck_dimu, DataError(RooAbsData::SumW2));
  bckMC_dimu.statOn(plotBck_dimu);

  pdfBck_dimu.plotOn(plotBck_dimu,LineColor(kGreen));
  // Draw the results
  newCanvas("BackgroundMinv_dimu");  
  //   TCanvas *mmm=new TCanvas("BackgroundMinv_dimu", "BackgroundMinv_dimu", 600, 600);
  //   mmm->cd();
  plotBck_dimu->Draw();


  /*
  // ---------------------------------------------------------------------------------
  // Fit the signal+background shape - diem
  cout << endl;
  cout << "-------------------------------------------------------------------------" << endl;
  cout << "[SIG+BCK] Fitting signal + backgorund shape - diem" << endl;
  RooPlot *plotMassAll_diem=likeLihoodOut.frame(binMin, binMax, nBins);
  allMC_diem.plotOn(plotMassAll_diem,DataError(RooAbsData::SumW2));
  allMC_diem.statOn(plotMassAll_diem);

  // Define the two pdfs
  // The signal fraction
  RooRealVar fracSigDiem("fracSigDiem","fraction of signal",0.2,0,1);
  fracSigDiem.setConstant(kFALSE);

  // Build the sum of signal and background pdfs  
  RooArgList pdfListDiem(pdfSig_diem, pdfBck_diem);
  RooArgList coeffListDiem(fracSigDiem);
  RooAddPdf pdfSum_diem("pdfSum_diem","Sig + Back pdf",pdfListDiem,coeffListDiem);
  
  // Fit
  RooFitResult *fitResultSum_diem=pdfSum_diem.fitTo(allMC_diem,Save(kTRUE),Extended(extended) );
  fitResultSum_diem->Print("v");
  pdfSum_diem.plotOn(plotMassAll_diem,LineColor(kRed));
  //   pdfSum_diem.plotOn(plotMassAll_diem,LineColor(kRed),Components(RooArgList(polBgAll_diem)));
  pdfSum_diem.paramOn(plotMassAll_diem,Label("f x Sign + (1-f) x Bck"));

  newCanvas("AllMinv_diem");
  plotMassAll_diem->Draw();

  cout << "--------------------------------------------------------------" << endl;
  */

  // Fit the signal+background shape - dimu
  cout << endl;
  cout << "-------------------------------------------------------------------------" << endl;
  cout << "[SIG+BCK] Fitting signal + backgorund shape - dimu" << endl;
  RooPlot *plotMassAll_dimu=likeLihoodOut.frame(binMin, binMax, nBins);
  allMC_dimu.plotOn(plotMassAll_dimu,DataError(RooAbsData::SumW2));
  allMC_dimu.statOn(plotMassAll_dimu);

  // The signal fraction
  RooRealVar fracSigDimu("fracSigDimu", "fraction of signal", 0.2, 0, 1);
  fracSigDimu.setConstant(kFALSE);

  // Build the sum of signal and background pdfs  
  RooArgList pdfListDimu(pdfSig_dimu, pdfBck_dimu);
  RooArgList coeffListDimu(fracSigDimu);
  RooAddPdf pdfSum_dimu("pdfSum_dimu", "Sig + Back pdf", pdfListDimu, coeffListDimu);

  // Fit
  RooFitResult *fitResultSum_dimu=pdfSum_dimu.fitTo(allMC_dimu, Save(kTRUE), Extended(extended), SumW2Error(kTRUE));
  fitResultSum_dimu->Print("v");
  pdfSum_dimu.plotOn(plotMassAll_dimu,LineColor(kRed));
  pdfSum_dimu.paramOn(plotMassAll_dimu,Label("f x Sig + (1-f) x Bck"));
  

  newCanvas("AllMinv_dimu");
  //   TCanvas *tmpCanv2=new TCanvas("AllMinv_dimu", "AllMinv_dimu", 600, 600);
  //   tmpCanv2->cd();
  plotMassAll_dimu->Draw();
 
  cout << "--------------------------------------------------------------" << endl;

  //-------------------------------------------------------------------------------
  /*  
  // Perform a simultaneous fit 
  RooSimultaneous simPdf("SimPdf","SimPdf",tagCat) ;
  simPdf.addPdf(pdfSum_diem,"diem") ;
  simPdf.addPdf(pdfSum_dimu,"dimu") ;
  RooFitResult *fitResultSimultan=simPdf.fitTo(allMC_allFinalState,Save(kTRUE));
  fitResultSimultan->Print("v");

  RooPlot *frame=likeLihoodOut.frame(binMin, binMax, nBins);
  newCanvas("test");
  allMC_allFinalState.plotOn(frame, Cut("FinalStateCat==FinalStateCat::diem"), DataError(RooAbsData::SumW2));
  tagCat.setLabel("diem");
  simPdf.plotOn(frame, LineColor(kRed), Slice(tagCat), ProjWData(tagCat,allMC_allFinalState));  
  frame->Draw();

  newCanvas("test1");
  RooPlot *frame1=likeLihoodOut.frame(binMin, binMax, nBins);
  allMC_allFinalState.plotOn(frame1, Cut("FinalStateCat==FinalStateCat::dimu"), DataError(RooAbsData::SumW2));
  tagCat.setLabel("dimu");
  simPdf.plotOn(frame1, LineColor(kRed), Slice(tagCat), ProjWData(tagCat, allMC_allFinalState));
  frame1->Draw();
  */

  //-------------------------------------------------------------------------------

  int nBins2=20;

  // Toy MC 
  // Test the significance acievable with the fit on MC ensembles generated using
  // the fitted distributions
  if(doToyMc) {
    cout << "---------- PRODUCE TOY MC ENSEMBLES -----------------------------------" << endl;
    int printLevel=-1;

    if(debug)
      printLevel=1;

    // Get the f value fitted on the full sample
    /*double fPrevDiem=fracSigDiem.getVal();*/
    double fPrevDimu=fracSigDimu.getVal();
    
    // Book some histograms
    /*
    TH1F *hFrac_toyMcSigBckDiem=new TH1F("hFrac_toyMcSigBckDiem", "f for sig+bck toy mc diem",20,0,1);
    TH1F *hSign_toyMcSigBckDiem=new TH1F("hSign_toyMcSigBckDiem", "Signif Sig+bck toy mc diem",50,0,10);
    */
    TH1F *hFrac_toyMcSigBckDimu=new TH1F("hFrac_toyMcSigBckDimu", "f for sig+bck toy mc dimu",20,0,1);
    TH1F *hSign_toyMcSigBckDimu=new TH1F("hSign_toyMcSigBckDimu", "Signif Sig+bck toy mc dimu",50,0,10);
    //TH1F *hFrac_toyMcSigBckSimultFit=new TH1F("hFrac_toyMcSigBckSimultFit",
    //			"f for sig+bck toy mc SimultFit",20,0,1);
    /*
    TH1F *hSign_toyMcSigBckSimultFit=new TH1F("hSign_toyMcSigBckSimultFit", "Signif Sig+bck toy mc SimultFit",50,0,10);
    */  
    // Define the variable to be used in the fit
    RooRealVar likeLihoodOut2("LikeLihoodOut","LikeLihoodOut",binMin,binMax);
    // Define the set of variables (ll out + tag for simult. fit)
    RooArgSet setOfVarUnw(likeLihoodOut2,tagCat);
    // Instantiate a random generator
    RndmGenerator rndm(37);
    
    for(int u=0; u!=nEvToyMc; ++u) { // Loop ove Toy MC experiments
      // Se the initial value of the parameter to what obtained on MC
      /*
      fracSigDiem.setVal(fPrevDiem);
      fracSigDiem.setConstant(kFALSE);
      */

      fracSigDimu.setVal(fPrevDimu);      
      fracSigDimu.setConstant(kFALSE);
  
      /*    
      // Extract the # of signal events accordingly to poisson distr.
      int nEvSigDiem=rndm.extractPoisson(nSigEvAfterPreselDiem);
      // Extract the # of background events accordingly to poisson distr.
      int nEvBckDiem=rndm.extractPoisson(nBckEvAfterPreselDiem);
      */
      // Extract the # of signal events accordingly to poisson distr.
      int nEvBckDimu=rndm.extractPoisson(nBckEvAfterPreselDimu);
      // Extract the # of background events accordingly to poisson distr.
      int nEvSigDimu=rndm.extractPoisson(nSigEvAfterPreselDimu);


      // Build the toy data sets
      /*
      RooDataSet toyMCSignDiem=	buildToyMCDataSample("diem", nEvSigDiem, pdfSig_diem);
      RooDataSet toyMCBckDiem=	buildToyMCDataSample("diem", nEvBckDiem, pdfBck_diem);
      */

      RooDataSet toyMCSignDimu=buildToyMCDataSample("dimu", nEvSigDimu, pdfSig_dimu);
      RooDataSet toyMCBckDimu=buildToyMCDataSample("dimu", nEvBckDimu, pdfBck_dimu);
      
      // Create the signal + background ensembles
      /*
      RooDataSet toyEnsDiem=toyMCSignDiem;
      toyEnsDiem.append(toyMCBckDiem);

      RooPlot *toyMCPlotDiem=likeLihoodOut2.frame(binMin, binMax, nBins2);
      toyEnsDiem.plotOn(toyMCPlotDiem);
      toyEnsDiem.statOn(toyMCPlotDiem);

      // Fit the diem ensemble
      RooFitResult *fitToyEnsDiem
	=pdfSum_diem.fitTo(toyEnsDiem,Save(kTRUE),PrintLevel(printLevel));
      if(debug)
	fitToyEnsDiem->Print("v");
      pdfSum_diem.plotOn(toyMCPlotDiem,LineColor(8));
      pdfSum_diem.paramOn(toyMCPlotDiem,Label("f x BW + (1-f) x pol2"));
      hFrac_toyMcSigBckDiem->Fill(fracSigDiem.getVal());
      */

      RooDataSet toyEnsDimu=toyMCSignDimu;
      toyEnsDimu.append(toyMCBckDimu);

      RooPlot *toyMCPlotDimu=likeLihoodOut2.frame(binMin, binMax, nBins2);
      toyEnsDimu.plotOn(toyMCPlotDimu);
      toyEnsDimu.statOn(toyMCPlotDimu);

      // Fit the dimu ensemble
      RooFitResult *fitToyEnsDimu=pdfSum_dimu.fitTo(toyEnsDimu, Save(kTRUE), SumW2Error(kTRUE), PrintLevel(printLevel));
      if(debug) fitToyEnsDimu->Print("v");
      pdfSum_dimu.plotOn(toyMCPlotDimu,LineColor(9));
      pdfSum_dimu.paramOn(toyMCPlotDimu,Label("f x BW + (1-f) x pol2"));
      hFrac_toyMcSigBckDimu->Fill(fracSigDimu.getVal());

      /*
      // Fit the two ensebles simultaneously
      // Perform a simultaneous fit 
      RooDataSet allMC_allFinalState=toyEnsDiem;
      allMC_allFinalState.append(toyEnsDimu);

      // define the joint pdf
      RooSimultaneous simPdf("SimPdf","SimPdf",tagCat) ;
      simPdf.addPdf(pdfSum_diem,"diem") ;
      simPdf.addPdf(pdfSum_dimu,"dimu") ;

      RooFitResult *fitResultSimultan=simPdf.fitTo(allMC_allFinalState,Save(kTRUE),
						     PrintLevel(printLevel));
      if(debug)
	fitResultSimultan->Print("v");
      */

      // Refit it in the backgorund only hipotesis
      // Build the sum of signal and background pdfs
      /*
      fracSigDiem.setVal(0);
      fracSigDiem.setConstant(kTRUE);
      */

      fracSigDimu.setVal(0);
      fracSigDimu.setConstant(kTRUE);

      // Build the likelihood of the background only hypot.
      /*
      RooNLLVar nllBckOnlyDiem("nllBckOnlyDiem", "nllBckOnlyDiem", pdfSum_diem, toyEnsDiem);
      double sigDiem=Utils::getSignifFromMinNLL(fitToyEnsDiem->minNll(),nllBckOnlyDiem.getVal());
      hSign_toyMcSigBckDiem->Fill(sigDiem);
      */

      RooNLLVar nllBckOnlyDimu("nllBckOnlyDimu", "nllBckOnlyDimu", pdfSum_dimu, toyEnsDimu);
      double sigDimu=Utils::getSignifFromMinNLL(fitToyEnsDimu->minNll(),nllBckOnlyDimu.getVal());
      hSign_toyMcSigBckDimu->Fill(sigDimu);

      /*
      RooNLLVar nllBckOnlySimultFit("nllBckOnlySimultFit","nllBckOnlySimultFit", simPdf,allMC_allFinalState);
      double sigSimulFit=Utils::getSignifFromMinNLL(fitResultSimultan->minNll(), nllBckOnlySimultFit.getVal());
      hSign_toyMcSigBckSimultFit->Fill(sigSimulFit);
      */

      if(debug) {
	/*
	cout << "- diem: " << endl;
	cout << "  The minNLL sig+bck: " << fitToyEnsDiem->minNll() << endl;
	cout << "  The minNLL: bck only (f=0)" << nllBckOnlyDiem.getVal() << endl;
	cout << "  S: " << sigDiem << endl;
	*/

	cout << "- dimu: " << endl;
	cout << "  The minNLL sig+bck: " << fitToyEnsDimu->minNll() << endl;
	cout << "  The minNLL: bck only (f=0)" << nllBckOnlyDimu.getVal() << endl;
	cout << "  S: " <<  sigDimu << endl;

	/*
	cout << "- SimultFit: " << endl;
	cout << "  The minNLL sig+bck: " << fitResultSimultan->minNll() << endl;
	cout << "  The minNLL: bck only (f=0)" << nllBckOnlySimultFit.getVal() << endl;
	cout << "  S: " << sigSimulFit << endl;
	*/
      }

      if(drawEnsemblePlots) {
	/*
	new TCanvas();
	toyMCPlotDiem->Draw();     
	*/

	new TCanvas();
	toyMCPlotDimu->Draw();     
      } else {
	/*
	delete toyMCPlotDiem;
	delete fitToyEnsDiem;
	delete fitResultSimultan;
	*/
	delete toyMCPlotDimu;
	delete fitToyEnsDimu;
      }
    }

    //TCanvas *c1=newCanvas(hFrac_toyMcSigBckDiem->GetName(),2,1,1,-1);
    //c1->cd(0);
    /*
    newCanvas(hFrac_toyMcSigBckDiem->GetName(),form);
    hFrac_toyMcSigBckDiem->Draw();
    newCanvas(hSign_toyMcSigBckDiem->GetName(),form);
    hSign_toyMcSigBckDiem->Draw();
    */
    
    newCanvas(hFrac_toyMcSigBckDimu->GetName(),form);
    hFrac_toyMcSigBckDimu->Draw();
    newCanvas(hSign_toyMcSigBckDimu->GetName(),form);
    hSign_toyMcSigBckDimu->Draw();

    /*
    //newCanvas(hFrac_toyMcSigBckSimultFit->GetName(),form);
    //hFrac_toyMcSigBckSimultFit->Draw();
    newCanvas(hSign_toyMcSigBckSimultFit->GetName(),form);
    hSign_toyMcSigBckSimultFit->Draw();
    */


    // Background only toy mc
    /*
    TH1F *hFrac_toyMcBckDiem=new TH1F("hFrac_toyMcBckDiem", "f for bck only toy mc diem",20,0,1);
    TH1F *hSign_toyMcBckDiem=new TH1F("hSign_toyMcBckDiem", "Signif bck only toy mc diem",50,0,10);
    */

    TH1F *hFrac_toyMcBckDimu=new TH1F("hFrac_toyMcBckDimu", "f for bck only toy mc dimu",20,0,1);
    TH1F *hSign_toyMcBckDimu=new TH1F("hSign_toyMcBckDimu", "Signif bck only toy mc dimu",50,0,10);

    /*
    //TH1F *hFrac_toyMcBckSimultFit=new TH1F("hFrac_toyMcBckSimultFit", "f for bck only toy mc SimultFit",20,0,1);
    TH1F *hSign_toyMcBckSimultFit=new TH1F("hSign_toyMcBckSimultFit", "Signif bck only toy mc SimultFit",50,0,10);
    TH1F *hSign2_toyMcBckSimultFit=new TH1F("hSign2_toyMcBckSimultFit", "Signif^2 bck only toy mc SimultFit",50,0,10);
    */


    for(int u=0; u!=nEvToyMc; u++) {
      /*
      fracSigDiem.setVal(fPrevDiem);
      fracSigDiem.setConstant(kFALSE);
      int nEvBckDiem=rndm.extractPoisson(nBckEvAfterPreselDiem);
      RooDataSet toyMCBckDiem=buildToyMCDataSample("diem", nEvBckDiem, pdfBck_diem);
      RooPlot *toyMCPlotDiem=likeLihoodOut2.frame(binMin, binMax, nBins2);
      toyMCBckDiem.plotOn(toyMCPlotDiem);
      */

      fracSigDimu.setVal(fPrevDimu);
      fracSigDimu.setConstant(kFALSE);
      int nEvBckDimu=rndm.extractPoisson(nBckEvAfterPreselDimu);
      RooDataSet toyMCBckDimu=buildToyMCDataSample("dimu", nEvBckDimu, pdfBck_dimu);
      RooPlot *toyMCPlotDimu=likeLihoodOut2.frame(binMin, binMax, nBins2);
      toyMCBckDimu.plotOn(toyMCPlotDimu);

      /*
      RooDataSet toyMCBckAll=toyMCBckDiem;
      toyMCBckAll.append(toyMCBckDimu);
      */

      /*
      RooFitResult *fitToyEnsDiem=pdfSum_diem.fitTo(toyMCBckDiem,Save(kTRUE),PrintLevel(printLevel));
      if(debug)	fitToyEnsDiem->Print("v");
      pdfSum_diem.plotOn(toyMCPlotDiem,LineColor(6));
      pdfSum_diem.paramOn(toyMCPlotDiem,Label("f x BW + (1-f) x pol2"));
      hFrac_toyMcBckDiem->Fill(fracSigDiem.getVal());
      */

      RooFitResult *fitToyEnsDimu=pdfSum_dimu.fitTo(toyMCBckDimu, Save(kTRUE), SumW2Error(kTRUE), PrintLevel(printLevel));
      if(debug)	fitToyEnsDimu->Print("v");
      pdfSum_dimu.plotOn(toyMCPlotDimu,LineColor(49));
      pdfSum_dimu.paramOn(toyMCPlotDimu,Label("f x BW + (1-f) x pol2"));
      hFrac_toyMcBckDimu->Fill(fracSigDimu.getVal());

      /*
      // Define the joint pdf
      RooSimultaneous simPdf("SimPdf","SimPdf",tagCat) ;
      simPdf.addPdf(pdfSum_diem,"diem") ;
      simPdf.addPdf(pdfSum_dimu,"dimu") ;

      RooFitResult *fitResultSimultan=simPdf.fitTo(toyMCBckAll,Save(kTRUE),PrintLevel(printLevel));
      if(debug) fitResultSimultan->Print("v");
      */


      
      // Refit it in the backgorund only hipotesis
      // Build the sum of signal and background pdfs
      /*
      fracSigDiem.setVal(0);
      fracSigDiem.setConstant(kTRUE);
      */

      fracSigDimu.setVal(0);
      fracSigDimu.setConstant(kTRUE);


      // Build the likelihood of the background only hypot.
      // Build the likelihood of the background only hypot.
      /*
      RooNLLVar nllBckOnlyDiem("nllBckOnlyDiem", "nllBckOnlyDiem", pdfSum_diem, toyMCBckDiem);
      double sigDiem=Utils::getSignifFromMinNLL(fitToyEnsDiem->minNll(), nllBckOnlyDiem.getVal());
      */

      RooNLLVar nllBckOnlyDimu("nllBckOnlyDimu", "nllBckOnlyDimu", pdfSum_dimu, toyMCBckDimu);
      double sigDimu=Utils::getSignifFromMinNLL(fitToyEnsDimu->minNll(), nllBckOnlyDimu.getVal());

      /*
      RooNLLVar nllBckOnlySimultFit("nllBckOnlySimultFit", "nllBckOnlySimultFit", simPdf,toyMCBckAll);
      double sigSimulFit=Utils::getSignifFromMinNLL(fitResultSimultan->minNll(), nllBckOnlySimultFit.getVal());
      */

      if(debug) {
	/*
	cout << "- diem: " << endl;
	cout << "  The minNLL sig+bck: " << fitToyEnsDiem->minNll() << endl;
	cout << "  The minNLL: bck only (f=0)" << nllBckOnlyDiem.getVal() << endl;
	cout << "  S: " << sigDiem << endl;
	*/

	cout << "- dimu: " << endl;
	cout << "  The minNLL sig+bck: " << fitToyEnsDimu->minNll() << endl;
	cout << "  The minNLL: bck only (f=0)" << nllBckOnlyDimu.getVal() << endl;
	cout << "  S: " <<  sigDimu << endl;

	/*
	cout << "- SimultFit: " << endl;
	cout << "  The minNLL sig+bck: " << fitResultSimultan->minNll() << endl;
	cout << "  The minNLL: bck only (f=0)" << nllBckOnlySimultFit.getVal() << endl;
	cout << "  S: " << sigSimulFit << endl;
	*/
      }

      hSign_toyMcBckDimu->Fill(sigDimu);
      /*
      hSign_toyMcBckDiem->Fill(sigDiem);
      hSign_toyMcBckSimultFit->Fill(sigSimulFit);
      hSign2_toyMcBckSimultFit->Fill(sigSimulFit*sigSimulFit);
      */

      if(drawEnsemblePlots) {
	/*
	new TCanvas();
	toyMCPlotDiem->Draw();
	*/
	new TCanvas();
	toyMCPlotDimu->Draw();     
      } else {
	delete toyMCPlotDimu;
	delete fitToyEnsDimu;
	/*
	delete toyMCPlotDiem;
	delete fitToyEnsDiem;
	delete fitResultSimultan;
	*/
      }
    }

    /*
    newCanvas(hFrac_toyMcBckDiem->GetName(),form);
    hFrac_toyMcBckDiem->Draw();
    newCanvas(hSign_toyMcBckDiem->GetName(),form);
    //TF1 *gaus_diem=new TF1("gaus_diem","gaus",0.2,100);
    //hSign_toyMcBckDiem->Fit(gaus_diem,"","",0.2,100);
    hSign_toyMcBckDiem->Draw();
    */
    
    newCanvas(hFrac_toyMcBckDimu->GetName(),form);
    hFrac_toyMcBckDimu->Draw();
    newCanvas(hSign_toyMcBckDimu->GetName(),form);
    //TF1 *gaus_dimu=new TF1("gaus_dimu","gaus",0.2,100);
    //hSign_toyMcBckDimu->Fit(gaus_dimu,"","",0.2,100);
    hSign_toyMcBckDimu->Draw();

    /*
    //newCanvas(hFrac_toyMcBckSimultFit->GetName(),form);
    //hFrac_toyMcBckSimultFit->Draw();
    newCanvas(hSign_toyMcBckSimultFit->GetName(),form);
    hSign_toyMcBckSimultFit->Draw();
    newCanvas(hSign2_toyMcBckSimultFit->GetName(),form);
    //TF1 *gaus_simult=new TF1("gaus_simult","gaus",0.2,100);
    //hSign2_toyMcBckSimultFit->Fit(gaus_simult,"","",0.2,100);
    hSign2_toyMcBckSimultFit->Draw();
    */

    /*double medianS_sigPlusBckDiem=computeHistoMedian(hSign_toyMcSigBckDiem);*/
    double medianS_sigPlusBckDimu=computeHistoMedian(hSign_toyMcSigBckDimu);
    /*double medianS_sigPlusBckSimultFit=computeHistoMedian(hSign_toyMcSigBckSimultFit);*/

    //double diem_Pvalue=gaus_diem->Integral(medianS_sigPlusBckDiem,100)/gaus_diem->Integral(-100,100);
    //double dimu_Pvalue=gaus_dimu->Integral(medianS_sigPlusBckDimu,100)/gaus_dimu->Integral(-100,100);
    //double simult_Pvalue=gaus_simult->Integral(medianS_sigPlusBckSimultFit,100)/gaus_simult->Integral(-100,100);

    /*
    cout << "----- diem: " << endl;
    cout << "  Median Significance (Sig+Bck): " << medianS_sigPlusBckDiem << endl;
    cout << "  Median Significance (Bck only): " << computeHistoMedian(hSign_toyMcBckDiem)
	 << endl;
    double pValueDiem=computePValue(hSign_toyMcBckDiem,medianS_sigPlusBckDiem);
    cout << "  P-value: " << pValueDiem << endl;
    cout << "  # sigma: " << TMath::ErfInverse(1.-2.*pValueDiem)*sqrt(2.) << endl;
    */

    cout << "----- dimu: " << endl;
    cout << "  Median Significance (Sig+Bck): " << medianS_sigPlusBckDimu << endl;
    cout << "  Median Significance (Bck only): " << computeHistoMedian(hSign_toyMcBckDimu)
	 << endl;
    double pValueDimu=computePValue(hSign_toyMcBckDimu,medianS_sigPlusBckDimu);
    cout << "  P-value: " << pValueDimu << endl;
    cout << "  # sigma: " << TMath::ErfInverse(1.-2.*pValueDimu)*sqrt(2.) << endl;

    /*
    cout << "----- SimultFit: " << endl;
    cout << "  Median Significance (Sig+Bck): " << medianS_sigPlusBckSimultFit << endl;
    cout << "  Median Significance (Bck only): " << computeHistoMedian(hSign_toyMcBckSimultFit)
	 << endl;
    double pValueSimultFit=computePValue(hSign_toyMcBckSimultFit,medianS_sigPlusBckSimultFit);
    cout << "  P-value: " << pValueSimultFit << endl;
    cout << "  # sigma: " << TMath::ErfInverse(1.-2.*pValueSimultFit)*sqrt(2.) << endl;
    */
  }

  /*
    TEMPORARILY EXCLUDED
  if(doData) {
    cout << "------------------------------------------------------------------------------" << endl;
    cout << " ---------------- FITTING DATA: ----------------------------------------------" << endl;

    // Retrieve the data ntuple
    TString dataFileName=inputDir+"DiLeptDistr_data_diem.root";
    TFile *fileDataDiem=new TFile(dataFileName.Data());
    TNtuple *dataNtupleDiem=(TNtuple*) fileDataDiem->Get("DiLeptNtuple");

    // Build the dataset
    RooDataSet dataDiem=*buildDataSet("data_diem", "Data diem", dataNtupleDiem, discDiem, "diem");
    //RooRealVar likeLihoodOut("LikeLihoodOut","LikeLihoodOut",0,1);
    //RooRealVar wgh("weight","weight",0,100); // event weight
    //wgh.setVal(1);
    //likeLihoodOut.setVal(0.99);
    //RooCategory tagCat("FinalStateCat", "Final state category") ;
    //tagCat.defineType("diem");
    //tagCat.defineType("dimu");
    //tagCat.setLabel("diem");

    //RooArgSet setOfVar(likeLihoodOut, wgh, tagCat);

    //dataDiem.add(setOfVar, 1);
    //dataDiem.add(setOfVar, 1);

    RooPlot *dataPlotDiem=likeLihoodOut.frame(binMin, binMax, nBins);
    dataDiem.plotOn(dataPlotDiem);
    dataDiem.statOn(dataPlotDiem);
    
    fracSigDiem.setConstant(kFALSE);
    fracSigDimu.setConstant(kFALSE);

    RooFitResult *fitresData_Diem=pdfSum_diem.fitTo(dataDiem,Save(kTRUE));
    fitresData_Diem->Print("v");
    pdfSum_diem.plotOn(dataPlotDiem,LineColor(kBlack));
    pdfSum_diem.paramOn(dataPlotDiem,Label("Gauss (x) BW pdf"));
    
    new TCanvas();
    dataPlotDiem->Draw();

    fracSigDiem.setVal(0);
    fracSigDiem.setConstant(kTRUE);
    
    // Build the likelihood of the background only hypot.
    RooNLLVar nllBckOnlyDiem("nllBckOnlyDiem", "nllBckOnlyDiem", pdfSum_diem, dataDiem) ;
    
    cout << "- diem (data): " << endl;
    cout << "  The minNLL sig+bck: " << fitresData_Diem->minNll() << endl;
    cout << "  The minNLL: bck only (f=0)" << nllBckOnlyDiem.getVal() << endl;
    double sigDiem=Utils::getSignifFromMinNLL(fitresData_Diem->minNll(), nllBckOnlyDiem.getVal());
    cout << "  S: " << sigDiem  << endl;
    //double pValueData=computePValue(hSign_toyMcBckDiem,sigDiem);
    //cout << "  # sigma: " << TMath::ErfInverse(1.-2.*pValueData)*sqrt(2.) << endl;

    TH1* hSignal=hTestStatSigDiem;
    TH1* bh=hTestStatBckDiem;
    TH1* dh=discDiem.evaluateData(L1xL2xL3xL4,dataNtupleDiem,"data", -1,nBins);

    TLimitDataSource* mydatasource=new TLimitDataSource(hSignal, bh, dh);
    TConfidenceLevel *myconfidence=TLimit::ComputeLimit(mydatasource,50000);
    cout << "# S: " << myconfidence->GetStot() << endl;
    cout << "# B: " << myconfidence->GetBtot() << endl;
    cout << "# D: " << myconfidence->GetDtot() << endl;

    cout << " GetStatistic(): " << myconfidence->GetStatistic() << endl;
    //  Get the expected statistic value in the background only hypothesis
    cout << "  GetExpectedStatistic_b(Int_t sigma): " << myconfidence->GetExpectedStatistic_b() << endl;
    //  Get the expected statistic value in the signal plus background hypothesis 
    cout << " GetExpectedStatistic_sb(Int_t sigma): " << myconfidence->GetExpectedStatistic_sb() << endl;
    //  Get the Confidence Level for the background only
    cout << " CLb(bool use_sMC): " << myconfidence->CLb() << endl;
    //  Get the Confidence Level for the signal plus background hypothesis
    cout << " CLsb(bool use_sMC): " << myconfidence->CLsb() << endl;
    //  Get the Confidence Level defined by CLs=CLsb/CLb.
    //  This quantity is stable w.r.t. background fluctuations.
    cout << " CLs(bool use_sMC): " << myconfidence->CLs() << endl;
    //  Get the expected Confidence Level for the signal plus background hypothesis
    //  if there is only background.
    cout << " GetExpectedCLsb_b(Int_t sigma): " << myconfidence->GetExpectedCLsb_b() << endl;
    //  Get the expected Confidence Level for the signal plus background hypothesis
    //  if there is only background.
    cout << " GetExpectedCLs_b(Int_t sigma): " << myconfidence->GetExpectedCLs_b() << endl;

    //  Get the expected Confidence Level for the background only
    //  if there is signal and background.
    cout << " GetExpectedCLb_sb(Int_t sigma): " << myconfidence->GetExpectedCLb_sb() << endl; // Questo definisce expected significance P-value=1 - GetExpectedCLb_sb
    //  Get the expected Confidence Level for the background only
    //  if there is only background.
    cout << " GetExpectedCLb_b(Int_t sigma): " << myconfidence->GetExpectedCLb_b() << endl;

    //  Get average CLsb. 
    cout << " GetAverageCLsb(): " << myconfidence->GetAverageCLsb() << endl;
    //  Get average CLs. 
    cout << " GetAverageCLs(): " << myconfidence->GetAverageCLs() << endl;
    //  Get 3s probability. 
    cout << " Get3sProbability(): " << myconfidence->Get3sProbability() << endl;
    //  Get 5s probability. 
    cout << " Get5sProbability(): " << myconfidence->Get5sProbability() << endl;
 

    //myconfidence->Dump();

    cout << "  CLs    : " << myconfidence->CLs()  << endl;
    cout << "  CLsb   : " << myconfidence->CLsb() << endl;
    cout << "  CLb    : " << myconfidence->CLb()  << endl;
    cout << "< CLs >  : " << myconfidence->GetExpectedCLs_b()  << endl;
    cout << "< CLsb > : " << myconfidence->GetExpectedCLsb_b() << endl;
    cout << "< CLb >  : " << myconfidence->GetExpectedCLb_b()  << endl;
    TCanvas *c2=new TCanvas("c2");
    myconfidence->Draw();


    ////     Add stat uncertainty
    //cout << endl << "Computing limits with stat systematics... " << endl;
    //TConfidenceLevel *mystatconfidence=TLimit::ComputeLimit(mydatasource,50000,true);
    //cout << "CLs    : "   << mystatconfidence->CLs()  << endl;
    //cout << "CLsb   : "   << mystatconfidence->CLsb() << endl;
    //cout << "CLb    : "   << mystatconfidence->CLb()  << endl;
    //cout << "< CLs >  : " << mystatconfidence->GetExpectedCLs_b()  << endl;
    //cout << "< CLsb > : " << mystatconfidence->GetExpectedCLsb_b() << endl;
    //cout << "< CLb >  : " << mystatconfidence->GetExpectedCLb_b()  << endl;

    //// Add some systematics
    //cout << endl << "Computing limits with systematics... " << endl;
    //TVectorD errorb(2);
    //TVectorD errors(2);
    //TObjArray* names=new TObjArray();
    //TObjString name1("bg uncertainty");
    //TObjString name2("sig uncertainty");
    //names->AddLast(&name1);
    //names->AddLast(&name2);
    //errorb[0]=0.05; // error source 1: 5%
    //errorb[1]=0;    // error source 2: 0%
    //errors[0]=0;    // error source 1: 0%
    //errors[1]=0.01; // error source 2: 1%
    //TLimitDataSource* mynewdatasource=new TLimitDataSource();
    //mynewdatasource->AddChannel(signal,background,data,&errors,&errorb,names);
    //TConfidenceLevel *mynewconfidence=TLimit::ComputeLimit(mynewdatasource,50000,true);
    //cout << "CLs    : " << mynewconfidence->CLs()  << endl;
    //cout << "CLsb   : " << mynewconfidence->CLsb() << endl;
    //cout << "CLb    : " << mynewconfidence->CLb()  << endl;
    //cout << "< CLs >  : " << mynewconfidence->GetExpectedCLs_b()  << endl;
    //cout << "< CLsb > : " << mynewconfidence->GetExpectedCLsb_b() << endl;
    //cout << "< CLb >  : " << mynewconfidence->GetExpectedCLb_b()  << endl;

    //// show canonical -2lnQ plots in a new canvas
    //// - The histogram of -2lnQ for background hypothesis (full)
    //// - The histogram of -2lnQ for signal and background hypothesis (dashed)
    //TCanvas *c2=new TCanvas("c2");
    //myconfidence->Draw();
  
    //delete myconfidence;
    //delete mydatasource;
    //delete mystatconfidence;
    //delete mynewconfidence;
    //delete mynewdatasource;

    delete myconfidence;
    delete mydatasource;

  }
  */

  /*
    TEMPORARILY EXCLUDED
  if(doEnsemble) {
    // ---------------------------------------------------------------------------------
    EnsembleProducer prodBck("Bck", bckNtuple_diem, nBckEvAfterPreselDiem);
    prodBck.setSeed(3);
    TH1F *hEnsNevBck=new TH1F("hEnsNevBck","# events in the ensemble", 50, 0, 50);
    hEnsNevBck->Sumw2();
    for(int i=0; i!=prodBck.getNEnsembles(); i++) {
      hEnsNevBck->Fill(prodBck.getEnsembleN(i)->GetEntries());
    }

    newCanvas(hEnsNevBck->GetName(),form);
    hEnsNevBck->Draw();
  
    TH1F * poisBck=buildPoisson(prodBck.getExpectedAverage(), hEnsNevBck->Integral());
    poisBck->SetLineColor(kRed);
    poisBck->Draw("same");

    // ---------------------------------------------------------------------------------
    EnsembleProducer prodSig("Sig", sigNtuple_diem, nSigEvAfterPreselDiem);
    prodSig.setSeed(3);
    TH1F *hEnsNevSig=new TH1F("hEnsNevSig","# events in the ensemble", 50, 0, 50);
    hEnsNevSig->Sumw2();
    for(int i=0; i!=prodSig.getNEnsembles(); i++) {
      //const TNtuple *ensamble=prodSig.getEnsembleN(i);
      hEnsNevSig->Fill(prodSig.getEnsembleN(i)->GetEntries());
    }

    newCanvas(hEnsNevSig->GetName(),form);
    hEnsNevSig->Draw();
    TH1F * poisSig=buildPoisson(prodSig.getExpectedAverage(), hEnsNevSig->Integral());
    poisSig->SetLineColor(kRed);
    poisSig->Draw("same");


    TH1F *hMinNLL_SigBck=new TH1F("hMinNLL_SigBck", "Min -Log(L) Sig+Bck", 100,0,200);
    TH1F *hFracSig_SigBck=new TH1F("hFracSig_SigBck", "f Sig+Bck",20,0,1);
    TH1F *hSign_SigBck=new TH1F("hSign_SigBck", "Signif Sig+bck",50,0,5);


    //   doSigPlusBck=true;
    //   doBckOnly=true;

    if(doSigPlusBck) {
      for(int jj=0; jj!= 1; ++jj) {
	//for(int jj=0; jj!= prodBck.getNEnsembles(); jj++) {
	cout << "-------------------------------------------------------------------------" << endl;
	cout << "[SIG+BCK / ENS] Fitting signal + backgorund ensemble # " << jj << endl;

	cout << "    # sig: " <<  prodSig.getEnsembleN(jj)->GetEntries() << endl;
	cout << "    # bck: " <<  prodBck.getEnsembleN(jj)->GetEntries() << endl;
	  
	// Sum signal and background ensembles
	TreeSum sum(prodSig.getEnsembleN(jj), prodBck.getEnsembleN(jj));
	TNtuple *theNewNtpl=sum.getSumTree();
	cout << "    # tot entries: " << theNewNtpl->GetEntries() << endl;

	// Build the dataset
	RooDataSet ensemple_diem=*buildDataSet("ensemble_diem","Signal MC dataset",
						 theNewNtpl, discDiem, "diem");
	// Plot the dataset
	RooPlot *ensPlot_diem=likeLihoodOut.frame(binMin, binMax, nBins2);
	ensemple_diem.plotOn(ensPlot_diem);
	ensemple_diem.statOn(ensPlot_diem);
      
	// Fit the spectrum to get f
      
	// Build the sum of signal and background pdfs
	fracSigDiem.setVal(0.5);
	fracSigDiem.setConstant(kFALSE);
      
	RooFitResult *fitResultSumEns_diem
	=pdfSum_diem.fitTo(ensemple_diem,Save(kTRUE),Extended(kFALSE));
	fitResultSumEns_diem->Print();
	pdfSum_diem.plotOn(ensPlot_diem,LineColor(kBlue));
	pdfSum_diem.paramOn(ensPlot_diem,Label("f x Sig + (1-f) x Bckg"));

	hFracSig_SigBck->Fill(fracSigDiem.getVal());

	// Refit it in the backgorund only hipotesis
	// Build the sum of signal and background pdfs
	fracSigDiem.setVal(0);
	fracSigDiem.setConstant(kTRUE);

	TH1* hSignal=hTestStatSigDiem;
	TH1* bh=hTestStatBckDiem;
	
	TH1* dh=discDiem.evaluateData(L1xL2xL3xL4,theNewNtpl,"data", -1,nBins);
	new TCanvas();
	dh->Draw();

	TLimitDataSource* mydatasource=new TLimitDataSource(hSignal,bh,dh);
	TConfidenceLevel *myconfidence=TLimit::ComputeLimit(mydatasource,50000);
	cout << "----------------" << endl;
	cout << "# S: " << myconfidence->GetStot() << endl;
	cout << "# B: " << myconfidence->GetBtot() << endl;
	cout << "# D: " << myconfidence->GetDtot() << endl;

	cout << " GetStatistic(): " << myconfidence->GetStatistic() << endl;
	//  Get the expected statistic value in the background only hypothesis
	cout << "  GetExpectedStatistic_b(Int_t sigma): " << myconfidence->GetExpectedStatistic_b() << endl;
	//  Get the expected statistic value in the signal plus background hypothesis 
	cout << " GetExpectedStatistic_sb(Int_t sigma): " << myconfidence->GetExpectedStatistic_sb() << endl;
	//myconfidence->Dump();

	cout << "  CLs    : " << myconfidence->CLs()  << endl;
	cout << "  CLsb   : " << myconfidence->CLsb() << endl;
	cout << "  CLb    : " << myconfidence->CLb()  << endl;
	cout << "< CLs >  : " << myconfidence->GetExpectedCLs_b()  << endl;
	cout << "< CLsb > : " << myconfidence->GetExpectedCLsb_b() << endl;
	cout << "< CLb >  : " << myconfidence->GetExpectedCLb_b()  << endl;
	cout << "----------------" << endl;


	// Build the likelihood of the background only hypot.
	RooNLLVar nll("nll","nll",pdfSum_diem,ensemple_diem) ;
    
	cout << "The minNLL sig+bck: " << fitResultSumEns_diem->minNll() << endl;
	cout << "The minNLL: bck only (f=0)" << nll.getVal() << endl;
	cout << "S: " << Utils::getSignifFromMinNLL(fitResultSumEns_diem->minNll(),
						    nll.getVal())
	     << endl;
	hSign_SigBck->Fill(Utils::getSignifFromMinNLL(fitResultSumEns_diem->minNll(),
						      nll.getVal()));

	// Draw the plot
	if(drawEnsemblePlots) {
	  newCanvas(theNewNtpl->GetName(),form);
	  ensPlot_diem->Draw();
	}
	hMinNLL_SigBck->Fill(fitResultSumEns_diem->minNll());
	delete fitResultSumEns_diem;
      }
    }

    newCanvas(hSign_SigBck->GetName(),form);
    hSign_SigBck->Draw();


    TH1F *hMinNLLBck=new TH1F("hMinNLLBck", "Min -Log(L) Bck only", 100,0,200);
    TH1F *hFracSig_Bck=new TH1F("hFracSig_Bck", "f Bck Only",20,0,1);

    TH1F *hSign_Bck=new TH1F("hSign_Bck", "Signif bck only",50,0,5);

    if(doBckOnly) {
      //   for(int jj=0; jj!= 5; jj++) {
      for(int jj=0; jj!= prodBck.getNEnsembles(); jj++) {
	cout << "-------------------------------------------------------------------------" << endl;
	cout << "[BCK Only/ ENS] Fitting backgorund only ensemble # " << jj << endl;
      
	TNtuple *theEnsNtpl=prodBck.getEnsembleN(jj);
	cout << "#entries: " << theEnsNtpl->GetEntries() << endl;

	RooDataSet mass("mass","mass",theEnsNtpl,likeLihoodOut);


	RooPlot *massPlotEnsBck=likeLihoodOut.frame(binMin, binMax, nBins2);
	mass.plotOn(massPlotEnsBck);
	mass.statOn(massPlotEnsBck);


	// Build the sum of signal and background pdfs
	fracSigDiem.setVal(0.5);
	fracSigDiem.setConstant(kFALSE);
      
	RooFitResult *fitResultSumEnsBck_diem
	=pdfSum_diem.fitTo(mass,Save(kTRUE),Extended(kFALSE));
	fitResultSumEnsBck_diem->Print();
	pdfSum_diem.plotOn(massPlotEnsBck,LineColor(kBlue));
	pdfSum_diem.paramOn(massPlotEnsBck,Label("f x BW + (1-f) x pol2"));

	hFracSig_Bck->Fill(fracSigDiem.getVal());

	// Refit it in the backgorund only hipotesis
	// Build the sum of signal and background pdfs
	fracSigDiem.setVal(0);
	fracSigDiem.setConstant(kTRUE);
      
	// Build the likelihood of the background only hypot.
	RooNLLVar nll("nll","nll",pdfSum_diem,mass) ;
    
	cout << "The minNLL sig+bck: " << fitResultSumEnsBck_diem->minNll() << endl;
	cout << "The minNLL: bck only (f=0)" << nll.getVal() << endl;
	cout << "S: " << Utils::getSignifFromMinNLL(fitResultSumEnsBck_diem->minNll(),
						    nll.getVal())
	     << endl;

	hSign_Bck->Fill(Utils::getSignifFromMinNLL(fitResultSumEnsBck_diem->minNll(),
						   nll.getVal()));




	if(drawEnsemblePlots) {
	  newCanvas(theEnsNtpl->GetName(),form);
	  massPlotEnsBck->Draw();
	}


	hMinNLLBck->Fill(fitResultSumEnsBck_diem->minNll());
	delete fitResultSumEnsBck_diem;
      
      }
    }



    newCanvas(hSign_Bck->GetName(),form);
    hSign_Bck->Draw();


    newCanvas(hMinNLL_SigBck->GetName(),form);
    hMinNLL_SigBck->Draw();
    newCanvas(hMinNLLBck->GetName(),form);
    hMinNLLBck->Draw();
    newCanvas(hFracSig_SigBck->GetName(),form);
    hFracSig_SigBck->Draw();

    newCanvas(hFracSig_Bck->GetName(),form);
    hFracSig_Bck->Draw();

    double medianS_sigPlusBck=computeHistoMedian(hSign_SigBck);
    cout << "Median Sign sig+bck: " << medianS_sigPlusBck << endl;
    cout << "Median Sign bck only: " << computeHistoMedian(hSign_Bck) << endl;
    double pValue=computePValue(hSign_Bck,medianS_sigPlusBck);
    cout << "P-value: " << pValue << endl;
    cout << "# sigma: " << TMath::ErfInverse(1.-2.*pValue)*sqrt(2.) << endl;
  }
  */
  
  cout << "----------------------------------------------------------------------" << endl;
  /*
  cout << "----- diem channel:" << endl;
  cout << "      # entries in Signal tree (un-weighted events): "
       << sigNtuple_diem->GetEntries() << endl;
  cout << "      # entries in Background tree (un-weighted events): "
       << bckNtuple_diem->GetEntries() << endl;
  cout << "    ----" << endl;
  cout << "      # of un-weighted events between " << massMinCut << " GeV and " << massMaxCut
       << " GeV is " << hTestStatSigDiem->GetEntries() << " (sig) and "
       <<  hTestStatBckDiem->GetEntries() << " (bck)" << endl;
  cout << "      # of weighted events between " << massMinCut << " GeV and " << massMaxCut
       << " GeV is " << nSigEvAfterPreselDiem <<  " (sig) and "
       <<  nBckEvAfterPreselDiem << " (bck)" << endl << endl;
  */
  cout << "----- dimu channel:" << endl;
  cout << "      # entries in Signal tree (un-weighted events): "
       << sigNtuple_dimu->GetEntries() << endl;
  cout << "      # entries in Background tree (un-weighted events): "
       << bckNtuple_dimu->GetEntries() << endl;
  cout << "    ----" << endl;
  cout << "      # of un-weighted events between " << massMinCut << " GeV and " << massMaxCut
       << " GeV is " << hTestStatSigDimu->GetEntries() << " (sig) and "
       <<  hTestStatBckDimu->GetEntries() << " (bck)" << endl;
  cout << "      # of weighted events between " << massMinCut << " GeV and " << massMaxCut
       << " GeV is " << nSigEvAfterPreselDimu <<  " (sig) and "
       <<  nBckEvAfterPreselDimu << " (bck)" << endl;
  cout << "----------------------------------------------------------------------" << endl;


  // Write down all histograms in memory  
  //TString rootFileName="FitHistos.root";
  //TFile hfile(rootFileName.Data(),"RECREATE");
  //hfile.cd();


  file->cd();
  gROOT->GetList()->Write();
  gROOT->GetListOfCanvases()->Write();
  

}

void normalize(TH1F *histo) {
  histo->Scale(1./histo->Integral());
}

void normalize(TH2F *histo) {
  histo->Scale(1./histo->Integral());
}




double expo(double *x, double *par){
  double X=x[0];
  double tau=0.018;
  return par[0]*TMath::Exp(-tau*X);
  
  //   return par[0]*TMath::Exp(-par[1]*X);
}


RooDataSet* buildDataSet(const TString& name, const TString& title,
			 TNtuple *ntuple, const LikelihoodDisc& disc, const TString& finalState) {

  RooRealVar likeLihoodOut("LikeLihoodOut", "LikeLihoodOut", 0, 1);
  RooRealVar wgh("weight", "weight", 0, 100); // event weight

  RooCategory tagCat("FinalStateCat", "Final state category") ;
  /*tagCat.defineType("diem");*/
  tagCat.defineType("dimu");

   RooArgSet setOfVar(likeLihoodOut,wgh,tagCat);
   //RooArgSet setOfVar(likeLihoodOut,tagCat);
   //RooDataSet *dataSet=new RooDataSet(name.Data(), title.Data(),setOfVar);
   RooDataSet *dataSet=new RooDataSet(name.Data(), title.Data(), setOfVar, WeightVar(wgh));
   //RooDataSet *dataSet=new RooDataSet(name.Data(), title.Data(), setOfVar, "weight");
   //dataSet->setWeightVar(wgh);
  dataSet->weightError(RooAbsData::SumW2);

  float diLeptInvMass, leadPt, /*subleadPt,*/
    dileptLeadDeltaPhi, leptMinusCmCosTheta,
    /*DiLeptDeltaEta, TransMass4Body, MinMass4Body,*/ 
    MinvTot, weight;
  /*
  float run, event, type; // must be float to 
  */

  // Open signal and background trees
  //ntuple->SetBranchAddress("run", &run);
  //ntuple->SetBranchAddress("event", &event);
  ntuple->SetBranchAddress("weight", &weight);
  //ntuple->SetBranchAddress("type", &type);
  ntuple->SetBranchAddress("diLeptInvMass", &diLeptInvMass);
  ntuple->SetBranchAddress("leadPt", &leadPt);
  //ntuple->SetBranchAddress("subleadPt", &subleadPt);
  ntuple->SetBranchAddress("dileptLeadDeltaPhi", &dileptLeadDeltaPhi);
  ntuple->SetBranchAddress("leptMinusCmCosTheta", &leptMinusCmCosTheta);
  if(doL11) ntuple->SetBranchAddress("MinvTot", &MinvTot);
  //ntuple->SetBranchAddress("DiLeptDeltaEta",&DiLeptDeltaEta);
  //ntuple->SetBranchAddress("TransMass4Body", &TransMass4Body);
  //ntuple->SetBranchAddress("MinMass4Body", &MinMass4Body);

  vector<TString> L1xL2xL3xL4;
  L1xL2xL3xL4.push_back("L1"); L1xL2xL3xL4.push_back("L2");
  L1xL2xL3xL4.push_back("L3"); L1xL2xL3xL4.push_back("L4");
  if(doL11)  L1xL2xL3xL4.push_back("L11"); // FIXME
  for(int entry=0; entry!=ntuple->GetEntries(); entry++) {
    ntuple->GetEntry(entry);
    if(name=="ensemble_diem") {
      weight=1;
    }

    // Cut on the di-lepton invariant mass....
    if(diLeptInvMass<massMinCut || diLeptInvMass>massMaxCut) continue;
    vector<float> L1xL2xL3xL4Vars;
    L1xL2xL3xL4Vars.push_back(diLeptInvMass);
    L1xL2xL3xL4Vars.push_back(leadPt);
    L1xL2xL3xL4Vars.push_back(dileptLeadDeltaPhi);
    L1xL2xL3xL4Vars.push_back(leptMinusCmCosTheta);
    if(doL11) L1xL2xL3xL4Vars.push_back(MinvTot); // FIXME
    
    //cout << "LL: " << disc.evaluate(L1xL2xL3xL4,L1xL2xL3xL4Vars) << endl;
    //cout << "weight: " << weight << endl;
    likeLihoodOut.setVal(disc.evaluate(L1xL2xL3xL4,L1xL2xL3xL4Vars));
    wgh.setVal(weight);
    tagCat.setLabel(finalState.Data());

    dataSet->add(setOfVar,weight);
  }
  ntuple->ResetBranchAddresses();
  return dataSet;
}


// RooDataSet mergeDiMuAndDiEm(const RooDataSet& dimuDataSet, const RooDataSet& diemDataSet) {
//   RooRealVar dileptminv("diLeptInvMass","M_{ll}",massMinCut,massMaxCut);
//   RooRealVar wgh("weight","weight",0,100); // event weight
//   RooCategory tagCat("FinalStateCat","Final state category") ;
//   tagCat.defineType("diem") ;
//   tagCat.defineType("dimu") ;
//   RooArgSet setOfVar(dileptminv,wgh,tagCat);
//   RooDataSet dataSet("pippo", "pippo",setOfVar,"weight");
//   dataSet.setWeightVar(wgh);
//   dataSet.weightError(RooAbsData::SumW2);
//   for(int i=0; i!=dimuDataSet.numEntries(); i++) {
//     double mass=dimuDataSet.get(i)->getRealValue("diLeptInvMass");
//     double weight=dimuDataSet.get(i)->getRealValue("weight");
//     cout << " mass=" << mass << " weight= " << weight 
// 	 << " final state: " << dimuDataSet.get(i)->getCatLabel("FinalStateCat") << endl;    
//     dileptminv.setVal(mass);
//     wgh.setVal(1);
//     tagCat.setLabel("diem");
//   }
//   return dataSet;
// }


RooDataSet buildToyMCDataSample(const TString& finalState, const int nEv, const RooAbsPdf& pdf) {
  RooRealVar likeLihoodOut("LikeLihoodOut","LikeLihoodOut",0.,1.);
  RooCategory tagCat("FinalStateCat","Final state category") ;
  /*tagCat.defineType("diem");*/
  tagCat.defineType("dimu");
  if(debug) cout << "   - NEV: " << nEv << endl;
  RooArgSet setOfVarUnw(likeLihoodOut, tagCat);

  RooDataSet ret("ToyMcDataSet", "ToyMcDataSet", setOfVarUnw);
  RooDataSet * tmp=pdf.generate(likeLihoodOut, nEv);

  for(int i=0; i!=nEv; ++i) {
    tagCat.setLabel(finalState.Data());
    double value=tmp->get(i)->getRealValue("LikeLihoodOut");
    //cout << " - [" << i << "] value: " << value << endl;
    likeLihoodOut.setVal(value);
    ret.add(setOfVarUnw);
  }
  delete tmp;
  return ret;
}


double function(double *x, double *par) {
  TF1 intGaus("intGaus", "gaus", 86, 96);
  intGaus.SetParameter(0, par[5]);
  intGaus.SetParameter(1, par[6]);
  intGaus.SetParameter(2, par[7]);

  if(x[0]<86) {
    return par[0]+
      par[1]*x[0]+
      par[2]*x[0]*x[0]+
      par[3]*x[0]*x[0]*x[0]+
      par[4]*x[0]*x[0]*x[0]*x[0];
  } 
  else if(86<=x[0] && x[0]<=96) {
    return intGaus.Eval(x[0]);
  } 
  else {             // x[0]>=96
    return par[8]+
      par[9]*x[0]+
      par[10]*x[0]*x[0]+
      par[11]*x[0]*x[0]*x[0]+
      par[12]*x[0]*x[0]*x[0]*x[0];
  }
}

