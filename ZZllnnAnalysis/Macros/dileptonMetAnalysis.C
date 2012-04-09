/*
 *  $Date: 2011/07/11 19:28:25 $
 *  \author D. Trocino   - Northeastern University
 */

/*
    Usage: 
      root -l
      .L dileptonMetAnalysis.C+
      dileptonMetAnalysis()

 */

#include "TSystem.h"
#include "TString.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TDirectoryFile.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TChain.h"
#include "TLorentzVector.h"
#include "TLine.h"
#include <iostream>
#include <map>
#include <sys/stat.h>
#include <sys/types.h>

#include "treeVariables.C"
#include "usefulFunctions.C"

using namespace std;

// 
// Function to fill plots
// 
void fillPlots(std::map<TString, TH1F*> &, TString, double wght = 1.0);

//
// Main function
//
void dileptonMetAnalysis() {

  // int loaded = -1;
  // if( gSystem->IsFileInIncludePath("usefulFunctions_C.so") )
  //   loaded = gSystem->Load("usefulFunctions_C.so");
  // if(loaded < 0)
  //   gSystem->CompileMacro("usefulFunctions.C");

  initializeGlobalVariables();

  // Samples location
  //TString indir="/home/daniele/Documents/Work/samples/2012-03-26_cmgTrees";
  TString indir="/home/daniele/Documents/Work/samples/2012-04-03_cmgTrees";
  TString basename_mc="MC_";
  TString basename_dt="Data_";
  TString base_mc=indir+"/"+basename_mc;
  TString base_dt=indir+"/"+basename_dt;

  // Integrated luminosity
  float intLumi = 4653.;
  TString integrLumi = "4653";

  // Pile-up scenario in data
  std::vector<double> dataPuDistribution;
  double alldatapileups[] = {1.22825e+07, 5.33316e+07, 1.27819e+08, 2.21102e+08, 3.09325e+08, 3.74101e+08, 4.09049e+08, 4.17488e+08, 4.06878e+08, 3.84466e+08, 3.55412e+08, 3.22755e+08, 2.88141e+08, 2.52593e+08, 2.17008e+08, 1.82346e+08, 1.49605e+08, 1.19698e+08, 9.33203e+07, 7.08679e+07, 5.24176e+07, 3.77691e+07, 2.65207e+07, 1.81566e+07, 1.21265e+07, 7.90616e+06, 5.03513e+06, 3.13451e+06, 1.90872e+06, 1.1377e+06, 664248, 380148, 213407, 117609, 63687.5, 0.}; // pixel-based
  // double alldatapileups[] = {1.344651e+07, 5.90653e+07, 1.409027e+08, 2.413012e+08, 3.337449e+08, 3.98711e+08, 4.301064e+08, 4.32283e+08, 4.138202e+08, 3.82846e+08, 3.451637e+08, 3.043438e+08, 2.62555e+08, 2.213308e+08, 1.819826e+08, 1.456898e+08, 1.134134e+08, 8.577886e+07, 6.301239e+07, 4.495959e+07, 3.116904e+07, 2.100786e+07, 1.377588e+07, 8796407, 5474418, 3323776, 1970638, 1142040, 647538.6, 359547.2, 195673.2, 104459.9, 54745.15, 28185.57, 28005.55, 0.008}; // hf-based

  unsigned int nDataPu = sizeof(alldatapileups)/sizeof(double);
  for(unsigned int pui=0; pui<nDataPu; ++pui) {
    dataPuDistribution.push_back( alldatapileups[pui] );
  }

  // Labels for all standard processes
  TString alllabels[]={"zh", "zz", "wz", "ww", "tt", "t_s", "tbar_s", "t_t", "tbar_t", "t_tw", "tbar_tw", "w", "ztautau", "z"};
  unsigned int nSmps=sizeof(alllabels)/sizeof(TString);

  // Labels for legend
  std::map<TString, TString> alllegendlabs;
  alllegendlabs["zh"]="ZH #rightarrow 2l+MET";
  alllegendlabs["zz"]="ZZ"; // alllegendlabs["zz"]="ZZ #rightarrow 2l2#nu";
  alllegendlabs["wz"]="WZ #rightarrow 3l#nu";
  alllegendlabs["ww"]="WW #rightarrow 2l2#nu";
  alllegendlabs["tt"]="tt";
  alllegendlabs["t_s"]="t (s ch.)";
  alllegendlabs["tbar_s"]="#bar{t} (s ch.)";
  alllegendlabs["t_t"]="t (t ch.)";
  alllegendlabs["tbar_t"]="#bar{t} (t ch.)";
  alllegendlabs["t_tw"]="tW";
  alllegendlabs["tbar_tw"]="#bar{t}W";
  alllegendlabs["w"]="W #rightarrow l#nu";
  alllegendlabs["ztautau"]="Z #rightarrow #tau#tau";
  alllegendlabs["z"]="Z #rightarrow ll";

  // Cross sections, process by process
  //double xSect[] = {0.01171, 5.9, 18.2, 43., 165., 3.19, 1.44, 41.92, 22.65, 7.87, 7.87, 31314., 3048.}; // old (for my own trees)
  double xSect[] = {0.01171, 4.287, 0.856, 4.78, 165., 3.19, 1.44, 41.92, 22.65, 7.87, 7.87, 31314., 3048., 3048.};
  unsigned int nXsect = sizeof(xSect)/sizeof(double);
  if( nXsect != nSmps ) {
    std::cout << " *************************** ERROR ***************************" << endl;
    std::cout << "    Number of cross sections != number of samples! Check!"      << std::endl;
    return;
  }

  // Branching ratios, process by process
  double brFra[] = {1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.};
  unsigned int nBrFra = sizeof(brFra)/sizeof(double);
  if( nBrFra != nSmps ) {
    std::cout << " *************************** ERROR ***************************" << endl;
    std::cout << "    Number of branching ratios != number of samples! Check!"    << std::endl;
    return;
  }

  // Colors and marker styles, process by process
  Color_t cols[]={kBlack, kRed, kMagenta, kViolet-1, kBlue, kCyan, kCyan, kCyan, kCyan, kCyan, kCyan, kSpring+3, kYellow-3, kYellow-7};
  Style_t mrks[]={25, 21, 22, 23, 24, 25, 25, 25, 25, 25, 25, 26, 27, 28, 29, 30, 31, 32};
  Style_t lines[]={1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  // Backup colors: kGreen, kPink, kAzure+1, kTeal+1, kOrange+1
  unsigned int nCols=sizeof(cols)/sizeof(Color_t);
  unsigned int nMrks=sizeof(mrks)/sizeof(Style_t);
  unsigned int nLines=sizeof(lines)/sizeof(Style_t);

  if( nCols<nSmps ) {
    std::cout << " ******************** ERROR ********************" << endl;
    std::cout << "    Not enough colors available, add colors!"     << std::endl;
    return;
  }

  if( nMrks<nSmps ) {
    std::cout << " ******************** ERROR ********************" << endl;
    std::cout << "    Not enough markers available, add styles!"    << std::endl;
    return;
  }

  if( nLines<nSmps ) {
    std::cout << " ************************ ERROR ************************" << endl;
    std::cout << "    Not enough line styles available, add line styles!"   << std::endl;
    return;
  }

  // Labels to be actually used
  TString alllabelstouse[]={"zz", "wz", "ww", "tt", "t_s", "tbar_s", "t_t", "tbar_t", "t_tw", "tbar_tw", "w", "z"};
  unsigned int nSmpsToUse=sizeof(alllabelstouse)/sizeof(TString);

  //  TString allhiggslabelstouse[]={"zh"};
  TString allhiggslabelstouse[]={};
  unsigned int nHiggsSmpsToUse=sizeof(allhiggslabelstouse)/sizeof(TString);

  if( nSmps < (nSmpsToUse+nHiggsSmpsToUse) ) {
    std::cout << " ******************** ERROR ********************" << endl;
    std::cout << "    Not enough samples available, add samples!"   << std::endl;
    return;
  }

  // Maps associating each process to its tree, color, marker and line style
  std::map<TString, std::vector<TString> > allfiles;
  std::map<TString, TChain*> alltrees;
  std::map<TString, double>  allgenev;
  std::map<TString, std::vector<double> >  allpileups;
  std::map<TString, double>  allxsect;
  std::map<TString, double>  allbrfra;
  std::map<TString, Color_t> allcols;
  std::map<TString, Style_t> allmrks;
  std::map<TString, Style_t> alllines;

  // Fill maps with files
  allfiles["zh"].push_back(base_mc+"ZH.root");
  allfiles["zz"].push_back(base_mc+"ZZ_0.root"); allfiles["zz"].push_back(base_mc+"ZZ_1.root");
  allfiles["wz"].push_back(base_mc+"WZ_0.root"); allfiles["wz"].push_back(base_mc+"WZ_1.root");
  allfiles["ww"].push_back(base_mc+"WW.root");
  allfiles["tt"].push_back(base_mc+"TTJets.root");
  allfiles["t_tw"].push_back(base_mc+"SingleT_tW.root"); allfiles["tbar_tw"].push_back(base_mc+"SingleTbar_tW.root");
  allfiles["t_t"].push_back(base_mc+"SingleT_t.root");   allfiles["tbar_t"].push_back(base_mc+"SingleTbar_t.root");
  allfiles["t_s"].push_back(base_mc+"SingleT_s.root");   allfiles["tbar_s"].push_back(base_mc+"SingleTbar_s.root");
  allfiles["w"].push_back(base_mc+"WJetsToLNu.root"); 
  allfiles["ztautau"].push_back(base_mc+"DYJetsToTauTau_0.root"); allfiles["ztautau"].push_back(base_mc+"DYJetsToTauTau_1.root"); allfiles["ztautau"].push_back(base_mc+"DYJetsToTauTau_2.root"); 
  allfiles["ztautau"].push_back(base_mc+"DYJetsToTauTau_3.root"); allfiles["ztautau"].push_back(base_mc+"DYJetsToTauTau_4.root"); allfiles["ztautau"].push_back(base_mc+"DYJetsToTauTau_5.root"); 
  allfiles["ztautau"].push_back(base_mc+"DYJetsToTauTau_6.root"); allfiles["ztautau"].push_back(base_mc+"DYJetsToTauTau_7.root"); allfiles["ztautau"].push_back(base_mc+"DYJetsToTauTau_8.root"); 
  allfiles["ztautau"].push_back(base_mc+"DYJetsToTauTau_9.root"); 
  allfiles["z"].push_back(base_mc+"DYJetsToLL_0.root"); allfiles["z"].push_back(base_mc+"DYJetsToLL_1.root"); allfiles["z"].push_back(base_mc+"DYJetsToLL_2.root"); 
  allfiles["z"].push_back(base_mc+"DYJetsToLL_3.root"); allfiles["z"].push_back(base_mc+"DYJetsToLL_4.root"); allfiles["z"].push_back(base_mc+"DYJetsToLL_5.root"); 
  allfiles["z"].push_back(base_mc+"DYJetsToLL_6.root"); allfiles["z"].push_back(base_mc+"DYJetsToLL_7.root"); allfiles["z"].push_back(base_mc+"DYJetsToLL_8.root"); 
  allfiles["z"].push_back(base_mc+"DYJetsToLL_9.root"); 

  for(unsigned int k=0; k<nSmps; ++k) {

    TH1D *hmcPileup = 0;
    alltrees[alllabels[k]]=new TChain("evAnalyzer/data");
    allgenev[alllabels[k]]=0.;

    for( unsigned int kk=0; kk<allfiles[alllabels[k]].size(); ++kk ) {

      // Fill maps with trees
      alltrees[alllabels[k]]->Add( allfiles[alllabels[k]][kk] );

      // ---------------------------
      // Fill maps with numbers of generated events
      allgenev[alllabels[k]] += ((TH1D*)(TFile::Open( allfiles[alllabels[k]][kk] )->Get("evAnalyzer/h2zz/cutflow")))->GetBinContent(1);


      // ---------------------------
      // Fill maps with pile-up scenarios
      if(kk==0) 
	hmcPileup = (TH1D*)(TFile::Open( allfiles[alllabels[k]][kk] )->Get("evAnalyzer/h2zz/pileup")); 
      else 
	hmcPileup->Add( (TH1D*)(TFile::Open( allfiles[alllabels[k]][kk] )->Get("evAnalyzer/h2zz/pileup")) );
    }

    // Fill maps with pile-up scenarios (cont'd)
    allpileups[alllabels[k]] = std::vector<double>();
    makeWeightDistribution(hmcPileup, dataPuDistribution, allpileups[alllabels[k]]);

    delete hmcPileup;

    // ---------------------------
    // Fill maps with cross sections
    allxsect[alllabels[k]]=xSect[k];

    // ---------------------------
    // Fill maps with branching ratios
    allbrfra[alllabels[k]]=brFra[k];

    // ---------------------------
    // Fill maps with colors
    allcols[alllabels[k]]=cols[k];

    // ---------------------------
    // Fill maps with marker styles
    allmrks[alllabels[k]]=mrks[k];

    // ---------------------------
    // Fill maps with line styles
    alllines[alllabels[k]]=lines[k];
  }

  // Data samples
  TChain *t_data=new TChain("evAnalyzer/data");
  // DoubleMu
  t_data->Add(base_dt+"DoubleMu2011A_0.root"); t_data->Add(base_dt+"DoubleMu2011A_1.root"); 
  t_data->Add(base_dt+"DoubleMu2011B_0.root"); t_data->Add(base_dt+"DoubleMu2011B_1.root"); 
  // DoubleElectron
  t_data->Add(base_dt+"DoubleElectron2011A_0.root"); t_data->Add(base_dt+"DoubleElectron2011A_1.root"); 
  t_data->Add(base_dt+"DoubleElectron2011B_0.root"); t_data->Add(base_dt+"DoubleElectron2011B_1.root"); 
  // MuUG
  t_data->Add(base_dt+"MuEG2011A_0.root"); t_data->Add(base_dt+"MuEG2011A_1.root"); 
  t_data->Add(base_dt+"MuEG2011B_0.root"); t_data->Add(base_dt+"MuEG2011B_1.root"); 


  //
  // ::::::::::::::::::::::::::::::::::::::::::::::::::::
  // ::::::::::::::::::::::::::::::::::::::::::::::::::::
  // ::::::::::::::::::::::::::::::::::::::::::::::::::::
  //

  // 
  // Final states to plot
  // 
  std::map<int, TString> finalStates;
  finalStates[1] = "mm";
  finalStates[2] = "ee";

  // 
  // Variables to be plotted + binning and range limits (optional)
  // 
  //   Arguments (all TString): 
  //    - variable "label"
  //    - x axis title [unit]
  //    - y axis title; if empty: "Events/<bin+unit>" 
  //    - string containing: n. bins, first bin, last bin
  //
  addVariable("dileptMass", "Dilepton invariant mass [GeV/c^{2}]", "", 60, 60., 120.);
  addVariable("dileptPt", "Dilepton p_{T} [GeV/c]", "", 60, 0., 120.);
  addVariable("jetNumber", "PFJet number (p_{T} > 30 GeV/c)", "", 6, -0.5, 5.5);
  addVariable("d0RedMet", "D0 Reduced MET [GeV]", "", 100, 0., 100.); 
  addVariable("metPtBalance", "PF MET/p_{T}(Z)", "", 20, -5., 5.); 
  addVariable("deltaPhiJetMet", "#Delta#phi(jet,MET) [rad]", "", 6., 0., 3.141592654); 
  addVariable("jetCsv", "CSV discriminator (PFJet p_{T} > 20 GeV)", "", 10, -2., 8.);
  addVariable("extraLeptonNumber", "Additional lepton number (loose sel.)", "", 4, -0.5, 3.5);

  //
  // ::::::::::::::::::::::::::::::::::::::::::::::::::::
  // ::::::::::::::::::::::::::::::::::::::::::::::::::::
  // ::::::::::::::::::::::::::::::::::::::::::::::::::::
  //

  //
  // Preselection cuts 
  //
  //  Arguments (both TString): 
  //   - cut name (will appear in the cut-flow table and plot
  //   - cut expression (will be used in Draw() function)
  // 
  //addCut("Lepton flavor", "cat==2");  // 1: mm, 2: ee, 3: em
  // - // - // addCut("Lepton isolation", "getDetIso(l1_trkIso,l1_ecalIso,l1_hcalIso,rho)<0.15 && getDetIso(l2_trkIso,l2_ecalIso,l2_hcalIso,rho)<0.15"); 
  //addCut("Invariant mass", "abs(getMass(l1_en+l2_en, l1_px+l2_px, l1_py+l2_py, l1_pz+l2_pz)-91.1876)<15.");
  // - // - // addCut("Invariant mass", "abs(getMass(l1_en+l2_en, l1_px+l2_px, l1_py+l2_py, l1_pz+l2_pz) - 91.1876) < 10."); // for overall reweighting
  // - // - // addCut("Invariant mass", "dileptInvMass>46.1876 && dileptInvMass<136.1876"); // for non-resonant background measurement (for plots!)
  // - // - // addCut("Invariant mass", "(dileptInvMass>46.1876 && dileptInvMass<76.1876) || (dileptInvMass>106.1876 && dileptInvMass<136.1876)"); // for non-resonant background measurement (for numbers!)
  //addCut("Jet veto", "jn==0");
  // - // - // addCut("Jet max pt", "jet_pt_0<30");
  // - // - // addCut("DeltaPhi(jet, MET)", "jet_min_deltaPhiJetMET>0.5");
  // - // - // addCut("Anti b-tagging", "jet_max_TCHE<2.0");
  //addCut("Third lepton veto", "ln==0");
  // - // - // addCut("Third lepton anti-veto", "numberOfLeptons>0.5");            // for WZ measurement
  // - // - // addCut("Third lepton flavor", "abs(abs(thirdLept_flavor)-13)<0.1"); // for WZ measurement (only when 3rd lepton is a muon)
  // - // - // addCut("ThirdLept-MET mass", "thirdLept_METleptTransvMass>0");
  // - // - // -- Diagnostic...
  // - // - // // addCut("Third lepton veto", "abs(LeptonVeto-1)<0.1"); //addCut("abs(thirdLept_flavor)==13");
  // - // - // // addCut("", "abs(LeptonVeto-"+lept3_all+")>0.5");
  // - // - // -- end Third lepton
  // - // - // addCut("D0 Red-MET", "D0redMet>40");
  // - // - // addCut("D0 Red-MET", d0RedMetOppositeJets+">55");
  // - // - // addCut("D0 Red-MET", "getD0RedMet(Flavor, dilepPROJLong, uncertPROJLong, unclPROJLong, sumOppositeJetPROJLong, dilepPROJPerp, uncertPROJPerp, unclPROJPerp, sumOppositeJetPROJPerp, 1)>50");
  //addCut("PF MET", "met_pt[0]>70");
  // - // - // addCut("D0 Red-MET", d0RedMetOppositeJets_cut);
  // - // - // addCut("CMS Red-MET", cmsRedMetOppositeJets+">70");
  // - // - // addCut("DeltaPhi(jet, MET)", "jet_min_deltaPhiJetMET>0.5");
  // - // - // addCut("MET-lept trans. mass.", "leadMetTransMass+subleadMetTransMass>150");

  // 
  // Weights
  // 
  //  Usage:
  //   addWeight( "weight-expression", sample-label )                       // to re-weight single sample
  //     OR
  //   addWeight( "weight-expression", sample-labels, nunmber-of-samples )  // to re-weight groups of samples 
  // 
  //  Arguments: 
  //   - "weight-expression": weight (TString)
  //   - sample-label: label of the sample to re-weight (const char*)
  //   - sample-labels: array of labels of the samples to re-weight (TString*)
  //   - nunmber-of-samples: number of samples to re-weight, i.e. dimension of the array above (unsigned int)
  // 
  //addWeight("getPuWeights(ngenITpu)", alllabels, nSmps); // PU reweighting (N.B. the vector of PU weights will change from sample to sample)
  //addWeight("getOverallNorm(cat)", alllabels, nSmps);    // overall normalization (from Z peak)

  /*
  // N.B. 
  // If you pass "alllabels, nSmps" as parameters, Higgs samples will be reweighted too. 
  // If you pass "alllabelstouse, nSmpsToUse" as parameters, Higgs samples will NOT be reweighted!
  // In this case, need to reweight Higgs samples separately, by passing "allhiggslabelstouse, nHiggsSmpsToUse" as parameters
  */

  // 
  // Variables to be plotted + binning and range limits (optional)
  //
  //   Arguments (all TString): 
  //    - name of branch in tree
  //    - x axis title [unit]
  //    - y axis title; if empty: "Events/<bin+unit>" 
  //    - string containing: n. bins, first bin, last bin
  //
  //addVariable("getMass(l1_en+l2_en, l1_px+l2_px, l1_py+l2_py, l1_pz+l2_pz)", "Dilepton invariant mass [GeV/c^{2}]", "", "60,60,120");
  //addVariable("getMass(l1_en+l2_en, l1_px+l2_px, l1_py+l2_py, l1_pz+l2_pz)", "Dilepton invariant mass [GeV/c^{2}]", "", "12,75,105");
  //addVariable("getMass(l1_en+l2_en, l1_px+l2_px, l1_py+l2_py, l1_pz+l2_pz)", "Dilepton invariant mass [GeV/c^{2}]", "", "1,81.1876,101.1876"); // for overall reweighting (mZ +- 10 GeV/c^2)
  //addVariable("getMass(l1_en+l2_en, l1_px+l2_px, l1_py+l2_py, l1_pz+l2_pz)", "Dilepton invariant mass [GeV/c^{2}]", "", "3,46.1876,136.1876"); // for non-resonant background measurement
  //addVariable("jn", "PFJet number (p_{T} > 15 GeV/c)", "", "6,-0.5,5.5");
  //addVariable("ln", "Additional lepton number (loose sel.)", "", "4,-0.5,3.5");
  //addVariable("met_pt[0]", "PF MET [GeV]", "", "100,0.,200.");
  //addVariable("met_pt[1]", "PV-assoc. MET [GeV]", "", "100,0.,200.");
  //addVariable("met_pt[6]", "PV-assoc. + fwd MET [GeV]", "", "100,0.,200.");
  //addVariable("met_pt[12]", "PV-assoc. + fwd MET + #beta corr. [GeV]", "", "100,0.,200.");
  //addVariable("getD0RedMet(l1_px, l1_py, l1_ptErr, l2_px, l2_py, l2_ptErr, htvec_px, htvec_py, met_pt, met_phi, cat)", "D0 Reduced MET [GeV]", "", "100,0,100"); 
  //addVariable("getD0RedMet(l1_px, l1_py, l1_ptErr, l2_px, l2_py, l2_ptErr, met_pt, met_phi, cat)", "D0 Reduced MET [GeV]", "", "100,0,400"); // ONLY FOR EVENTS WITH NO JETS (FOR NOW!)

  //addVariable("nvtx", "Number of good reconstructed vertices", "", "50, 0.5, 50.5");

  //addVariable("patMet", "PF MET [geV]", "", "100,0.,400.");
  //addVariable("getD0RedMet(Flavor, dilepPROJLong, uncertPROJLong, unclPROJLong, sumOppositeJetPROJLong, dilepPROJPerp, uncertPROJPerp, unclPROJPerp, sumOppositeJetPROJPerp)", 
  //  	      "D0 Reduced MET [GeV]", "", "100,0.,400.");
  //addVariable("getCmsRedMet(Flavor, dilepPROJLong, uncertPROJLong, unclPROJLong, sumAllJetPROJLong, dilepPROJPerp, uncertPROJPerp, unclPROJPerp, sumAllJetPROJPerp)", 
  //  	      "CMS Ind. min. red. MET (all jets) [GeV]", "", "100,0.,400.");

  /// RedMET components
  //addVariable("dilepPROJLong", "Dilepton longitudinal p_{T} [GeV/c]", "", "40,0.,200.");
  //addVariable("dilepPROJPerp", "Dilepton perpendicular p_{T} [GeV/c]", "", "60,0.,300.");
  //addVariable("uncertPROJLong", "Lepton uncert. corr. longitudinal [GeV/c]", "", "50,-50.,0.");
  //addVariable("uncertPROJPerp", "Lepton uncert. corr. perpendicular [GeV/c]", "", "40,-2.,2.");
  //addVariable("sumOppositeJetPROJLong", "Sum of opposite jets longitudinal p_{T} [GeV/c]", "", "51,-250.,5.");
  //addVariable("sumOppositeJetPROJPerp", "Sum of opposite jets perpendicular p_{T} [GeV/c]", "", "61,-300.,5.");
  //addVariable("sumAllJetPROJLong", "Sum of jets longitudinal p_{T} [GeV/c]", "", "40,-250.,150.");
  //addVariable("sumAllJetPROJPerp", "Sum of jets perpendicular p_{T} [GeV/c]", "", "40,-250.,150.");
  //addVariable("METPROJLong", "Longitudinal pfMET [GeV/c]", "", "50,-150.,100.");
  //addVariable("METPROJPerp", "Perpendicular pfMET [GeV/c]", "", "60,-200.,100.");
  //addVariable("unclPROJLong", "Longitudinal pfMET+dilepton [GeV]", "", "30,-100.,200.");
  //addVariable("unclPROJPerp", "Perpendicular pfMET+dilepton [GeV]", "", "50,-100.,400.");
  //addVariable("unclPROJLong+sumOppositeJetPROJLong", "Longitudinal pfMET+dilepton+opp.jets [GeV]", "", "30,-100.,200.");
  //addVariable("unclPROJPerp+sumOppositeJetPROJPerp", "Perpendicular pfMET+dilepton+opp.jets [GeV]", "", "50,-100.,400.");
  //addVariable("unclPROJLong+sumAllJetPROJLong", "Longitudinal unclustered energy [GeV]", "", "30,-100.,200.");
  //addVariable("unclPROJPerp+sumAllJetPROJPerp", "Perpendicular unclustered energy [GeV]", "", "50,-100.,400.");
  //addVariable("recoilOppositeJetPROJLong", "Longitudinal recoil (opp.jets) [GeV]", "", "41,-200.,5.");
  //addVariable("recoilOppositeJetPROJPerp", "Perpendicular recoil (opp.jets) [GeV]", "", "41,-400.,10.");
  //addVariable("recoilAllJetPROJLong", "Longitudinal recoil [GeV]", "", "41,-200.,5.");
  //addVariable("recoilAllJetPROJPerp", "Perpendicular recoil [GeV]", "", "41,-400.,10.");

  //addVariable("getD0RedMet(Flavor, dilepPROJLong, uncertPROJLong, unclPROJLong, sumOppositeJetPROJLong, dilepPROJPerp, uncertPROJPerp, unclPROJPerp, sumOppositeJetPROJPerp)", 
  //  	      "D0 Reduced MET [GeV]", "", "50,0,150");
  //addVariable("getCmsRedMet(Flavor, dilepPROJLong, uncertPROJLong, unclPROJLong, sumAllJetPROJLong, dilepPROJPerp, uncertPROJPerp, unclPROJPerp, sumAllJetPROJPerp)", 
  //  	      "CMS Ind. min. red. MET (all jets) [GeV]", "", "50,0,150");
  //addVariable("getCmsRedMet(Flavor, dilepPROJLong, uncertPROJLong, unclPROJLong, sumOppositeJetPROJLong, dilepPROJPerp, uncertPROJPerp, unclPROJPerp, sumOppositeJetPROJPerp)", 
  //  	      "CMS Ind. min. red. MET (all jets) [GeV]", "", "50,0,150");
  //addVariable("getCmsRedMet(Flavor, dilepPROJLong, uncertPROJLong, sumOppositeJetPROJLong, sumOppositeJetPROJLong, dilepPROJPerp, uncertPROJPerp, sumOppositeJetPROJPerp, sumOppositeJetPROJPerp)", 
  //  	      "CMS Ind. min. red. MET (all jets) [GeV]", "", "50,0,150");
  //addVariable("getCmsRedMet(Flavor, dilepPROJLong, uncertPROJLong, unclPROJLong, unclPROJLong, dilepPROJPerp, uncertPROJPerp, unclPROJPerp, unclPROJPerp)", 
  //  	      "CMS Ind. min. red. MET (all jets) [GeV]", "", "50,0,150");
  //addVariable("getCmsRedMet(Flavor, dilepPROJLong, uncertPROJLong, unclPROJLong, sumAllJetPROJLong, dilepPROJPerp, uncertPROJPerp, unclPROJPerp, sumAllJetPROJPerp)", 
  //  	      "CMS Ind. min. red. MET (all jets) [GeV]", "", "50,0,150");
  //addVariable("getCmsRedMet(Flavor, dilepPROJLong, uncertPROJLong, sumAllJetPROJLong, sumAllJetPROJLong, dilepPROJPerp, uncertPROJPerp, sumAllJetPROJPerp, sumAllJetPROJPerp)", 
  //  	      "CMS Ind. min. red. MET (all jets) [GeV]", "", "50,0,150");
  //addVariable("getCmsRedMet(Flavor, dilepPROJLong, uncertPROJLong, unclPROJLong, unclPROJLong, dilepPROJPerp, uncertPROJPerp, unclPROJPerp, unclPROJPerp)", 
  //  	      "CMS Ind. min. red. MET (all jets) [GeV]", "", "50,0,150");

  // addVariable(d0RedMetOppositeJets, "D0 Reduced MET [GeV]", "", "50,0,150");
  //addVariable("jet_min_deltaPhiJetMET", "#Delta#phi(jet,MET) [rad]", "", "6,0,3.141592654");
  // addVariable("leadMetTransMass+subleadMetTransMass", "M_{T}(MET,l_{1}) + M_{T}(MET,l_{2}) [GeV/c^{2}]", "", "40,0,400");
  // addVariable("myDeltaPhi(leadPhi, patMetPhi)", "#Delta#phi(l_{1},MET) [rad]", "", "16,0,3.141592654");
  // addVariable("myDeltaPhi(subleadPhi, patMetPhi)", "#Delta#phi(l_{2},MET) [rad]", "", "16,0,3.141592654");
  //addVariable("myDeltaPhi(dileptPhi, patMetPhi)", "#Delta#phi(dilepton,MET) [rad]", "", "8,0.,3.141592654");
  //addVariable("myDeltaPhi(subleadPhi, patMetPhi):myDeltaPhi(leadPhi, patMetPhi)", "#Delta#phi(l_{1},MET) [rad]", "#Delta#phi(l_{1},MET) [rad]", "32,0,3.141592654,32,0,3.141592654");
  //addVariable("csCosThetaAbs(leadPt, leadEta, leadPhi, leadCharge, subleadPt, subleadEta, subleadPhi)", 
  //	      "cos^{2}#theta_{CS}", "", "8,0,1");
  //addVariable("dileptMetDeltaPhi(dilepPROJLong, dilepPROJPerp, METPROJLong, METPROJPerp)", "#Delta#phi(dilepton,MET)", "rad", "16,0,3.141592654");

  //addVariable("patMetPhi", "PF MET #phi [rad]", "", "16,2.,3.141592654");

  //addVariable("dileptPt", "Dilepton p_{T} [GeV/c]", "", "40,0.,400.");
  // addVariable("dileptMetTransMass", "Dilepton-MET transv. mass [GeV/c^{2}]", "", "50,0.,500.");
  // addVariable("dileptMetTransMassZ", "Dilepton-MET transv. mass [GeV/c^{2}]", "", "50,0.,500.");
  // addVariable("dileptMetTransMassZZ", "Dilepton-MET transv. mass [GeV/c^{2}]", "", "50,0.,500.");
  // addVariable("myDeltaPhi(leadPhi, subleadPhi)", "#Delta#phi(l_{1},l_{2} [rad]", "", "24,-3.141592654,3.141592654");


  //addVariable("dileptInvMass", "Dilepton invariant mass [GeV/c^{2}]", "", "15,76,106");
  //addVariable("leadMetTransMass+subleadMetTransMass", "M_{T}(MET,l_{1}) + M_{T}(MET,l_{2}) [GeV/c^{2}]", "", "25,100,350");
  // addVariable("abs(abs(leadPhi-subleadPhi)-(6.283185307*(abs(leadPhi-subleadPhi)>3.141592654)))", "#Delta#phi(l_{1},l_{2}) [rad]", "", "24,6-6,6");//"24,0,3.141592654");
  // addVariable("dileptLeadDeltaPhi", "#Delta#phi(dilept,l_{1}) [rad]", "", "12,0,1.2");
  // addVariable("leptMinusCmCosTheta", "cos#theta_{CM}(l^{-})", "", "12,-1,1");

  //addVariable("dileptMetTransMass", "Dilepton-MET transverse mass [GeV/c^{2}]", "", "");
  //addVariable("dileptMetTransMassZ", "Dilepton-MET transverse mass [GeV/c^{2}]", "", "");
  //addVariable("dileptMetTransMassZZ", "Dilepton-MET transverse mass [GeV/c^{2}]", "", "");
  //addVariable("jet_max_JBP", "max b-tag discriminator JBP", "", "");
  //addVariable("jet_max_SSVHE", "max b-tag discriminator SSVHE", "", "");

  //  - binning format: xbin, xmin, xmax, (ybin, ymin, ymax, (zbin, zmin, zmax))
  //addVariable("dileptInvMass", "Dilepton invariant mass [GeV/c^{2}]", "", "60,60,120");
  //addVariable(cmsRedMetAllJets, "CMS Ind. min. red. MET (all jets) [GeV]", "", "50,0,150");
  //addVariable("jet_min_deltaPhiJetMET", "36,0,3.141592654");
  //addVariable("jet_max_TCHE", "50,-2,8");
  //addVariable("jet_pt_0", "First PF Jet p_{T} [GeV/c]", "", "35,15,85");
  //addVariable("jet_pt_1", "50,15,115");
  //addVariable("jet_pt_2", "50,15,115");
  //addVariable("thirdLept_flavor");
  //addVariable("thirdLept_pixHits");
  //addVariable("thirdLept_trkHits");
  //addVariable("thirdLept_isGlobalMuonPT");
  //addVariable("thirdLept_isTrackerMuon");
  //addVariable("thirdLept_muMatches");
  //addVariable("thirdLept_trkChi2");
  //addVariable("thirdLept_pt");
  //addVariable("thirdLept_corrRelIso");
  //addVariable("thirdLept_METleptTransvMass", "3^{rd} lept-MET transverse mass [GeV/c^{2}]", "", "60, 0, 180");

  //addVariable(d0RedMetAllJets, "100,0,200");
  //addVariable(d0RedMetOppositeJets, "100,0,200");
  //addVariable(cmsRedMet, "100,0,200");
  //addVariable(cmsRedMetAllJets, "100,0,200");
  //addVariable(cmsRedMetOppositeJets, "100,0,200");
  //addVariable("CMSredMet", "100,0,200");
  //addVariable("totalNVertex", "Number of reconstructed vertices", "", "35, 0.5, 35.5");
  //addVariable("goodNVertex", "Number of good reconstructed vertices", "", "35, 0.5, 35.5");

  //
  // ::::::::::::::::::::::::::::::::::::::::::::::::::::
  // ::::::::::::::::::::::::::::::::::::::::::::::::::::
  // ::::::::::::::::::::::::::::::::::::::::::::::::::::
  //
  // Drawing options
  //  0: 1D, no stack;  1: 1D, stack;  2: 2D, scatter
  int drawOpts=1;
  if(drawOpts==2) {
    allcols["z"]=kOrange+1;
  }
  //
  // Stack plots? (yes: "";  no: "nostack")
  TString dostack="";
  //
  // Global drawing options (parameter to be passed to "Draw" function)
  TString gdopt="";
  //
  // Legend options
  TString legopt="lpf";
  //
  // Invert plot order?
  bool invertOrder=false;
  //
  // Draw data?
  bool drawData=true;

  //
  // Draw single-top samples altogether
  bool mergeTopSamples=true;

  //
  // Bin-by-bin comparison?
  bool binbybinComp=true;

  // 
  // Draw only MC cumulative plot? (Sum of all MC samples)
  bool drawOnlyCumulative = false;

  // 
  // Produce plots and/or cut-table?
  bool doPlot=true;
  bool doCutTable=false;

  double canvx(500.), canvy(500.);
  if(binbybinComp && drawData && (!doCutTable)) {
    //canvx=500.;
    canvy=625.;
  }

  // -- // -- // -- // -- // -- // -- // -- // -- // -- // -- // -- // 
  // -- // -- // -- // -- // -- // -- // -- // -- // -- // -- // -- // 
  // -- // -- // -- // -- // -- // -- // -- // -- // -- // -- // -- // 

  // 
  // Plots
  // 
  if(doPlot) {
    std::map<TString, THStack*> allstacks;
    std::map<TString, TCanvas*> allcanvas;
    std::map<TString, TPad*>    allpads1;
    //std::map<TString, TPad*>  allpads2;
    std::map<TString, TPad*>    allpads3;
    std::map<TString, TLegend*> alllegends;
    //std::map<TString, std::vector<TH1F*> > allhistos;
    //std::map<TString, std::vector<TH1F*> > allhiggshistos;
    std::map<TString, TH1F*> allhistos;
    std::map<TString, TH1F*> allhiggshistos;
    std::map<TString, TH1F*> allmchistos;
    //std::map<TString, TH1F*> allmcbkghistos;
    std::map<TString, TH1F*> alldatahistos;
    //std::map<TString, TH1F*> allpoisdiffhistos;
    std::map<TString, TH1F*> allfracdiffhistos;

    // if(binnings[i].Length()==0) {
    // 	cout << " *** New rule: you MUST select a binning!!! ***" << endl;
    // 	throw std::exception();
    // 	return;
    // }
    // TString thisbins=binnings[i]( binnings[i].Index("(")+1, binnings[i].Index(")")-binnings[i].Index("(")-1);
    // TObjArray *bintok=thisbins.Tokenize(","); bintok->SetOwner(kTRUE);
    // TString binsStr=((TObjString*)(*bintok)[0])->GetString(); 
    // TString lowBinStr=((TObjString*)(*bintok)[1])->GetString(); 
    // TString highBinStr=((TObjString*)(*bintok)[2])->GetString(); 
    // int bins=binsStr.Atoi();
    // double lowBin=lowBinStr.Atof();
    // double highBin=highBinStr.Atof();

    // // Some manipulation of the string, to avoid symbols in the names of canvases/stacks/etc. 
    // TString varlab=variables[i];
    // varlab.ReplaceAll(" ", ""); 
    // varlab.ReplaceAll("[", "_"); varlab.ReplaceAll("]", "_"); 
    // varlab.ReplaceAll("(", "_"); varlab.ReplaceAll(")", "_"); 
    // varlab.ReplaceAll("<", "LT"); varlab.ReplaceAll(">", "GT"); 
    // varlab.ReplaceAll("+", "_P_"); varlab.ReplaceAll("-", "_M_"); 
    // varlab.ReplaceAll("*", "_T_"); varlab.ReplaceAll("/", "_O_"); 
    // varlab.ReplaceAll(",", "_V_"); 
    /*
    allstacks.push_back(new THStack( ("stack_"+varlab).Data(), ("stack_"+varlab).Data() ));
    allcanvas.push_back(new TCanvas( ("canvas_"+varlab).Data(), ("canvas_"+varlab).Data(), canvx, canvy ));
    if(binbybinComp && drawData) {
	allcanvas.back()->cd();
	allpads1.push_back(new TPad( ("pad1_"+varlab).Data(), ("pad1_"+varlab).Data(), 0., 1.-canvx/canvy, 1., 1. ));
	//allpads2.push_back(new TPad( ("pad2_"+varlab).Data(), ("pad2_"+varlab).Data(), 0., (1.-canvx/canvy)/2., 1., 1.-canvx/canvy ));
	allpads3.push_back(new TPad( ("pad3_"+varlab).Data(), ("pad3_"+varlab).Data(), 0., 0., 1., 1.-canvx/canvy ));
	allpads1.back()->Draw();
	//allpads2.back()->Draw();
	allpads3.back()->Draw();
    }
    else {
	allcanvas.back()->cd();
	allpads1.push_back(new TPad( ("pad1_"+varlab).Data(), ("pad1_"+varlab).Data(), 0., 0., 1., 1. ));
	allpads1.back()->Draw();
    }
    alllegends.push_back(new TLegend(0.76, 0.64, 0.98, 0.94));
    alllegends.back()->SetFillColor(0);
    alllegends.back()->SetFillStyle(0);
    alllegends.back()->SetBorderSize(0);
    alllegends.back()->SetLineColor(0);
    allhistos.push_back(std::vector<TH1F*>());
    THStack *stk=allstacks.back();
    */

    // 
    // All histos and other stuff
    // 
    for(std::map<int, TString>::iterator finStat = finalStates.begin(); finStat!=finalStates.end(); ++finStat) {
      for(std::vector<TString>::iterator varsToPlot = variables.begin(); varsToPlot!=variables.end(); ++varsToPlot) {
	TString plotIdx = finStat->second+"_"+(*varsToPlot);
	allstacks[plotIdx] = new THStack( ("stack_"+plotIdx).Data(), ("stack_"+plotIdx).Data() ); 
	//allhistos[plotIdx] = std::vector<TH1F*>(); 
	//allhiggshistos[plotIdx] = std::vector<TH1F*>();
	alllegends[plotIdx] = new TLegend(0.76, 0.64, 0.98, 0.94);
	alllegends[plotIdx]->SetFillColor(0);
	alllegends[plotIdx]->SetFillStyle(0);
	alllegends[plotIdx]->SetBorderSize(0);
	alllegends[plotIdx]->SetLineColor(0);
      }
    }

    int j_ah = -1;          // index on allhistos elements
    unsigned int max_t = 0; // total number of top samples
    unsigned int idx_t = 0; // index (1-max_t) of top sample
    for(unsigned int t0=0; t0<nSmpsToUse; ++t0) {
	if( alllabelstouse[t0].BeginsWith("t_") || 
	    alllabelstouse[t0].BeginsWith("tbar_") ) 
	  max_t++;
    }

    // 
    // Loop over samples
    //
    for(unsigned int j00=0; j00<nSmpsToUse; ++j00) {
	unsigned int j=(invertOrder ? nSmpsToUse-1-j00 : j00);

	//if( alllabelstouse[j].CompareTo("zh") == 0 ) continue;

	TString plotname = "h_"+alllabelstouse[j]+"_";
	TString sumplot = "";
	bool addThisPlot = true;

	bool mergingTop = (alllabelstouse[j].BeginsWith("t_") || alllabelstouse[j].BeginsWith("tbar_")) && mergeTopSamples;

	if( mergingTop ) {
	  idx_t += 1;
	  plotname = "h_t_";
	  if( idx_t>1 ) {
	    sumplot = "+";
	    addThisPlot = false;
	  }
	}

	// Switch to pile-up scenarion for this MC sample
	puWeights = allpileups[alllabelstouse[j]];

	TString dopt="";
	if(j_ah > -1) dopt="same";
	//allpads1[i]->cd();
	if(addThisPlot) {
	  for(std::map<int, TString>::iterator finStat = finalStates.begin(); finStat!=finalStates.end(); ++finStat) {
	    for(std::vector<TString>::iterator varsToPlot = variables.begin(); varsToPlot!=variables.end(); ++varsToPlot) {
	      // TString plotIdx = (*varsToPlot)+"_"+finStat->second;
	      // allhistos[plotIdx].push_back( new TH1F( (plotname+plotIdx).Data(), (plotname+plotIdx).Data(), 
	      // 					      binsUInt[*varsToPlot], firstBinsDouble[*varsToPlot], lastBinsDouble[*varsToPlot] ) );
	      // allhistos[plotIdx].back()->Sumw2();
	      TString plotIdx = plotname+finStat->second+"_"+(*varsToPlot);
	      allhistos[plotIdx] = new TH1F( plotIdx.Data(), plotIdx.Data(), binsUInt[*varsToPlot], 
					     firstBinsDouble[*varsToPlot], lastBinsDouble[*varsToPlot] );
	      allhistos[plotIdx]->Sumw2();
	      j_ah += 1;
	    }
	  }
	}

	// alltrees[alllabelstouse[j]]->Draw( (variables[i]+">>"+sumplot+plotname).Data(),
	// 				   ("("+generalWeight+")*("+allweightspersample[alllabelstouse[j]]+")*("+allcuts+")").Data(),
	// 				   dopt.Data() );

	unsigned int nEvt = alltrees[alllabelstouse[j]]->GetEntries();
	std::cout << " - " << alllabelstouse[j].Data() << "  " << nEvt << std::endl;

	initializeTreeVariables();
	attachToTree(alltrees[alllabelstouse[j]]);

	// General weight
	float genWeight = intLumi * allxsect[alllabelstouse[j]] * allbrfra[alllabelstouse[j]] / allgenev[alllabelstouse[j]];

	// 
	// Loop over events
	// 
	for(unsigned int iEvt=0; iEvt<nEvt; ++iEvt) {
	  alltrees[alllabelstouse[j]]->GetEntry(iEvt);

	  if( finalStates.count(cat)>0 ) {

	    TString plotIdx = plotname+finalStates[cat]+"_";
	    // Event weight
	    float evtWeight = genWeight * getPuWeights(ngenITpu);

	    fillPlots(allhistos, plotIdx, evtWeight);

	    // float invmass = getMass(l1_en+l2_en, l1_px+l2_px, l1_py+l2_py, l1_pz+l2_pz);
	    // allhistos["mass"].back()->Fill(invmass, evtWeight);
	  }
	}
	// 
	// End loop over events
	// 

	if(mergingTop && idx_t<max_t) continue;  // stop here until you have all the top samples

	for(std::map<int, TString>::iterator finStat = finalStates.begin(); finStat!=finalStates.end(); ++finStat) {
	  for(std::vector<TString>::iterator varsToPlot = variables.begin(); varsToPlot!=variables.end(); ++varsToPlot) {
	    TString fsAndVar = finStat->second+"_"+(*varsToPlot);
	    TString plotIdx = plotname+fsAndVar;

	    if( allmchistos.count(fsAndVar)==0 ) {
	      allmchistos["h_sumAllMc_"+fsAndVar] = (TH1F*)allhistos[plotIdx]->Clone( ("h_sumAllMc_"+fsAndVar).Data() );
	    }
	    else {
	      allmchistos["h_sumAllMc_"+fsAndVar]->Add(allhistos[plotIdx]);
	    }

	    // Drawing options
	    if(drawOpts==0) {        // 1D, no stack
	      allhistos[plotIdx]->SetLineWidth(2);
	      allhistos[plotIdx]->SetLineColor(allcols[alllabelstouse[j]]);
	      dostack="nostack";
	      gdopt="HIST";
	    }
	    else if(drawOpts==1) {   // 1D, stack
	      allhistos[plotIdx]->SetLineWidth(1); 
	      allhistos[plotIdx]->SetFillColor(allcols[alllabelstouse[j]]);
	      gdopt="HIST";
	    }
	    else if(drawOpts==2) {   // 2D, scatter
	      allhistos[plotIdx]->SetMarkerColor(allcols[alllabelstouse[j]]);
	      allhistos[plotIdx]->SetFillColor(allcols[alllabelstouse[j]]);  // only for legend
	      gdopt="SCAT";
	      drawData=false;
	      legopt="f";
	    }

	    TString thislabel = alllegendlabs[alllabelstouse[j]];
	    if(mergingTop) thislabel = "t";
	    alllegends[fsAndVar]->AddEntry(allhistos[plotIdx], thislabel.Data(), legopt.Data());
	    allstacks[fsAndVar]->Add(allhistos[plotIdx]);
	  }
	}
    }
    // 
    // End loop over samples
    //


    // 
    // Loop over Higgs samples
    //
    //allhiggshistos["mass"] = std::vector<TH1F*>();

    for(unsigned int j=0; j<nHiggsSmpsToUse; ++j) {
	TString plotname = "h_"+allhiggslabelstouse[j]+"_";

	// Switch to pile-up scenarion for this MC sample
	puWeights = allpileups[allhiggslabelstouse[j]];

	for(std::map<int, TString>::iterator finStat = finalStates.begin(); finStat!=finalStates.end(); ++finStat) {
	  for(std::vector<TString>::iterator varsToPlot = variables.begin(); varsToPlot!=variables.end(); ++varsToPlot) {
	    TString plotIdx = plotname+finStat->second+"_"+(*varsToPlot);
	    allhiggshistos[plotIdx] = new TH1F( plotIdx.Data(), plotIdx.Data(), binsUInt[*varsToPlot], 
						firstBinsDouble[*varsToPlot], lastBinsDouble[*varsToPlot] );
	    allhiggshistos[plotIdx]->Sumw2();
	  }
	}

	unsigned int nEvt = alltrees[allhiggslabelstouse[j]]->GetEntries();
	std::cout << " - " << allhiggslabelstouse[j].Data() << "  " << nEvt << std::endl;

	initializeTreeVariables();
	attachToTree(alltrees[allhiggslabelstouse[j]]);

	// General weight
	float genWeight = intLumi * allxsect[allhiggslabelstouse[j]] * allbrfra[allhiggslabelstouse[j]] / allgenev[allhiggslabelstouse[j]];

	// 
	// Loop over events
	// 
	for(unsigned int iEvt=0; iEvt<nEvt; ++iEvt) {
	  alltrees[allhiggslabelstouse[j]]->GetEntry(iEvt);

	  if( finalStates.count(cat)>0 ) {

	    TString plotIdx = plotname+finalStates[cat]+"_";
	    // Event weight
	    float evtWeight = genWeight * getPuWeights(ngenITpu);

	    fillPlots(allhiggshistos, plotIdx, evtWeight);

	    // float invmass = getMass(l1_en+l2_en, l1_px+l2_px, l1_py+l2_py, l1_pz+l2_pz);
	    // allhistos["mass"].back()->Fill(invmass, evtWeight);
	  }
	}

	for(std::map<int, TString>::iterator finStat = finalStates.begin(); finStat!=finalStates.end(); ++finStat) {
	  for(std::vector<TString>::iterator varsToPlot = variables.begin(); varsToPlot!=variables.end(); ++varsToPlot) {
	    TString fsAndVar = finStat->second+"_"+(*varsToPlot);
	    TString plotIdx = plotname+fsAndVar;

	    // Drawing options for Higgs sample
	    if(drawOpts==0 || drawOpts==1) {        // 1D, no stack
	      allhiggshistos[plotIdx]->SetLineWidth(2);
	      allhiggshistos[plotIdx]->SetLineColor(allcols[allhiggslabelstouse[j]]);
	      allhiggshistos[plotIdx]->SetLineStyle(alllines[allhiggslabelstouse[j]]);
	    }
	    else if(drawOpts==2) {   // 2D, scatter
	      allhiggshistos[plotIdx]->SetMarkerColor(allcols[allhiggslabelstouse[j]]);
	      allhiggshistos[plotIdx]->SetFillColor(allcols[allhiggslabelstouse[j]]);  // only for legend
	    }

	    alllegends[fsAndVar]->AddEntry(allhiggshistos[plotIdx], alllegendlabs[allhiggslabelstouse[j]].Data(), legopt.Data());
	  }
	}
    }
    // 
    // End loop over Higgs samples
    //

    // 
    // Data plots
    //
    if(drawData) {
      //allpads1[i]->cd();
        //alldatahistos["mass"] = new TH1F( "h_data_mass", "h_data_mass", 60, 60., 120.);

	for(std::map<int, TString>::iterator finStat = finalStates.begin(); finStat!=finalStates.end(); ++finStat) {
	  for(std::vector<TString>::iterator varsToPlot = variables.begin(); varsToPlot!=variables.end(); ++varsToPlot) {
	    TString plotIdx = "h_data_"+finStat->second+"_"+(*varsToPlot);
	    alldatahistos[plotIdx] = new TH1F( plotIdx.Data(), plotIdx.Data(), binsUInt[*varsToPlot], 
					       firstBinsDouble[*varsToPlot], lastBinsDouble[*varsToPlot] );
	  }
	}

	// THIS IS THE STANDARD, UNCOMMENT IT!!!
	// t_data->Draw( (variables[i]+">>h_data_"+varlab).Data(),
	// 	      allcuts.Data(),
	// 	      "pesame" );
	// THIS IS ONLY FOR THIRD-LEPTON-REWEIGHTING FOR WZ MEASUREMENT!!!
	//std::cout << "WATCH OUT!!! You're using the third-lepton-reweighting for WZ measurement!!!" << std::endl;
	//t_data->Draw( (variables[i]+">>h_data_"+varlab).Data(),
	//	      ("(getWzInefficiency(Flavor, thirdLept_flavor, run, thirdLept_eta))*("+allcuts+")").Data(),
	//	      "pesame" );
	// alldatahistos.push_back( (TH1F*)gDirectory->Get( ("h_data_"+varlab).Data() ) );

	initializeTreeVariables();
	attachToTree(t_data);

	unsigned int nEvt = t_data->GetEntries();
	std::cout << " - data " << nEvt << std::endl;

	// 
	// Loop over events
	// 
	for(unsigned int iEvt=0; iEvt<nEvt; ++iEvt) {
	  t_data->GetEntry(iEvt);

	  if( finalStates.count(cat)>0 ) {
	    TString plotIdx = "h_data_"+finalStates[cat]+"_";
	    fillPlots(alldatahistos, plotIdx, 1.0);
	  }
	}

	for(std::map<int, TString>::iterator finStat = finalStates.begin(); finStat!=finalStates.end(); ++finStat) {
	  for(std::vector<TString>::iterator varsToPlot = variables.begin(); varsToPlot!=variables.end(); ++varsToPlot) {
	    TString fsAndVar = finStat->second+"_"+(*varsToPlot);
	    TString plotIdx = "h_data_"+fsAndVar;
	    alldatahistos[plotIdx]->SetMarkerStyle(20);
	    alllegends[fsAndVar]->AddEntry(alldatahistos[plotIdx], "data", legopt.Data());
	  }
	}
    }
    // 
    // End data plot
    //

    //
    // Final Draw
    //
    for(std::map<int, TString>::iterator finStat = finalStates.begin(); finStat!=finalStates.end(); ++finStat) {
      for(std::vector<TString>::iterator varsToPlot = variables.begin(); varsToPlot!=variables.end(); ++varsToPlot) {
	TString fsAndVar = finStat->second+"_"+(*varsToPlot);

	allcanvas[fsAndVar] = new TCanvas( ("canvas_"+fsAndVar).Data(), ("canvas_"+fsAndVar).Data(), canvx, canvy );
	if(binbybinComp && drawData) {
	  allcanvas[fsAndVar]->cd();
	  allpads1[fsAndVar] = new TPad( ("pad1_"+fsAndVar).Data(), ("pad1_"+fsAndVar).Data(), 0., 1.-canvx/canvy, 1., 1. );
	  //allpads2[fsAndVar] = new TPad( "pad2_mass", "pad2_mass", 0., (1.-canvx/canvy)/2., 1., 1.-canvx/canvy );
	  allpads3[fsAndVar] = new TPad( ("pad3_"+fsAndVar).Data(), ("pad3_"+fsAndVar).Data(), 0., 0., 1., 1.-canvx/canvy );
	  allpads1[fsAndVar]->Draw();
	  //allpads2[fsAndVar]->Draw();
	  allpads3[fsAndVar]->Draw();
	}
	else {
	  allcanvas[fsAndVar]->cd();
	  allpads1[fsAndVar] = new TPad( ("pad1_"+fsAndVar).Data(), ("pad1_"+fsAndVar).Data(), 0., 0., 1., 1. );
	  allpads1[fsAndVar]->Draw();
	}

	allpads1[fsAndVar]->cd();

	if(drawOpts!=2) 
	  allpads1[fsAndVar]->SetLogy(1);

	// Draw all MC samples
	if(drawOnlyCumulative==false) 
	  allstacks[fsAndVar]->Draw((dostack+gdopt).Data());
	else 
	  allmchistos["h_sumAllMc_"+fsAndVar]->Draw(gdopt.Data());

	// Draw all Higgs samples
	for(unsigned int j=0; j<nHiggsSmpsToUse; ++j) {
	  TString plotIdx = "h_"+allhiggslabelstouse[j]+"_"+fsAndVar;
	  allhiggshistos[plotIdx]->Draw((gdopt+"same").Data());
	}

	// Draw data
	if(drawData) 
	  alldatahistos["h_data_"+fsAndVar]->Draw((gdopt+"pesame").Data());

	if(drawOnlyCumulative==false) 
	  alllegends[fsAndVar]->Draw();


	// Axis, labels & Co.
	if(yTitles[*varsToPlot].EndsWith("Events/")) {
	  double onebin = ( lastBinsDouble[*varsToPlot] - firstBinsDouble[*varsToPlot] ) / binsUInt[*varsToPlot];
	  yTitles[*varsToPlot] += onebin;
	  TString unit = xTitles[*varsToPlot]( xTitles[*varsToPlot].Index("[")+1, ( xTitles[*varsToPlot].Index("]")-xTitles[*varsToPlot].Index("[")-1 ) );
	  if(unit.Length()>0) {
	    unit.Prepend(" ");
	    yTitles[*varsToPlot] += unit;
	  }
	  if(yTitles[*varsToPlot].EndsWith("/1")) yTitles[*varsToPlot].ReplaceAll("/1","");
	}
	if(drawOnlyCumulative==false) {
	  allstacks[fsAndVar]->GetXaxis()->SetTitle(xTitles[*varsToPlot].Data());
	  allstacks[fsAndVar]->GetXaxis()->SetTitleSize(0.05);
	  allstacks[fsAndVar]->GetXaxis()->SetLabelSize(0.04);
	  allstacks[fsAndVar]->GetYaxis()->SetTitle(yTitles[*varsToPlot].Data());
	  allstacks[fsAndVar]->GetYaxis()->SetTitleSize(0.05);
	  allstacks[fsAndVar]->GetYaxis()->SetLabelSize(0.04);
	  allstacks[fsAndVar]->GetYaxis()->SetTitleOffset(1.3);
	}
	else {
	  allmchistos["h_sumAllMc_"+fsAndVar]->GetXaxis()->SetTitle(xTitles[*varsToPlot].Data());
	  allmchistos["h_sumAllMc_"+fsAndVar]->GetXaxis()->SetTitleSize(0.05);
	  allmchistos["h_sumAllMc_"+fsAndVar]->GetXaxis()->SetLabelSize(0.04);
	  allmchistos["h_sumAllMc_"+fsAndVar]->GetYaxis()->SetTitle(yTitles[*varsToPlot].Data());
	  allmchistos["h_sumAllMc_"+fsAndVar]->GetYaxis()->SetTitleSize(0.05);
	  allmchistos["h_sumAllMc_"+fsAndVar]->GetYaxis()->SetLabelSize(0.04);
	  allmchistos["h_sumAllMc_"+fsAndVar]->GetYaxis()->SetTitleOffset(1.3);
	}

    /*
    // Axis, labels & Co.
    if(ytitles[i].EndsWith("Events/")) {
	int xbin=allmchistos[i]->GetXaxis()->GetNbins();
	double xmin=allmchistos[i]->GetXaxis()->GetXmin();
	double xmax=allmchistos[i]->GetXaxis()->GetXmax();
	double onebin=(xmax-xmin)/xbin;
	ytitles[i]+=onebin;
	TString unit=xtitles[i]( xtitles[i].Index("[")+1, ( xtitles[i].Index("]")-xtitles[i].Index("[")-1 ) );
	if(unit.Length()>0) {
	  unit.Prepend(" ");
	  ytitles[i]+=unit;
	}
	if(ytitles[i].EndsWith("/1")) ytitles[i].ReplaceAll("/1","");
    }
    if(drawOnlyCumulative==false) {
	stk->GetXaxis()->SetTitle(xtitles[i].Data());
	stk->GetXaxis()->SetTitleSize(0.05);
	stk->GetXaxis()->SetLabelSize(0.04);
	stk->GetYaxis()->SetTitle(ytitles[i].Data());
	stk->GetYaxis()->SetTitleSize(0.05);
	stk->GetYaxis()->SetLabelSize(0.04);
	stk->GetYaxis()->SetTitleOffset(1.3);
    }
    else {
	allmchistos[i]->GetXaxis()->SetTitle(xtitles[i].Data());
	allmchistos[i]->GetXaxis()->SetTitleSize(0.05);
	allmchistos[i]->GetXaxis()->SetLabelSize(0.04);
	allmchistos[i]->GetYaxis()->SetTitle(ytitles[i].Data());
	allmchistos[i]->GetYaxis()->SetTitleSize(0.05);
	allmchistos[i]->GetYaxis()->SetLabelSize(0.04);
	allmchistos[i]->GetYaxis()->SetTitleOffset(1.3);
    }
    */

    // "Bin-by-bin" data-MC comparison
	if(binbybinComp && drawData) {	// Start bin-by-bin data-MC comparison

	  double normdt = alldatahistos["h_data_"+fsAndVar]->Integral(0, -1);
	  double normmc = allmchistos["h_sumAllMc_"+fsAndVar]->Integral(0, -1);
	  if(false) cout << normdt << " - " << normmc << endl;  // just to avoid warnings at compilation time ("unused variable")
	  // allpoisdiffhistos.push_back( (TH1F*)allmchistos["h_sumAllMc_"+fsAndVar]->Clone( ("h_poisDiff_"+varlab).Data() ) );
	  // allpoisdiffhistos.back()->Reset();
	  // allpoisdiffhistos.back()->GetXaxis()->SetTitle("");
	  // allpoisdiffhistos.back()->GetYaxis()->SetTitle("No. Poisson #sigma");
	  // allpoisdiffhistos.back()->GetYaxis()->SetTitleSize(0.12);
	  // allpoisdiffhistos.back()->GetYaxis()->SetTitleOffset(0.5);
	  // allpoisdiffhistos.back()->GetXaxis()->SetLabelSize(0.10);
	  // allpoisdiffhistos.back()->GetYaxis()->SetLabelSize(0.10);
	  // allpoisdiffhistos.back()->SetMarkerStyle(20);
	  // allpoisdiffhistos.back()->SetMarkerSize(0.8);
	  // allpoisdiffhistos.back()->SetMarkerColor(kBlue);
	  allfracdiffhistos[fsAndVar] = (TH1F*)allmchistos["h_sumAllMc_"+fsAndVar]->Clone( "h_fracDiff_mass" );
	  allfracdiffhistos[fsAndVar]->Reset();
	  allfracdiffhistos[fsAndVar]->GetXaxis()->SetTitle("");
	  allfracdiffhistos[fsAndVar]->GetYaxis()->SetTitle("Rel. difference");
	  allfracdiffhistos[fsAndVar]->GetYaxis()->SetTitleSize(0.12);
	  allfracdiffhistos[fsAndVar]->GetYaxis()->SetTitleOffset(0.5);
	  allfracdiffhistos[fsAndVar]->GetXaxis()->SetLabelSize(0.10);
	  allfracdiffhistos[fsAndVar]->GetYaxis()->SetLabelSize(0.10);
	  allfracdiffhistos[fsAndVar]->GetYaxis()->SetRangeUser(-1., 1.);
	  allfracdiffhistos[fsAndVar]->SetMarkerStyle(20);
	  allfracdiffhistos[fsAndVar]->SetMarkerSize(0.8);
	  allfracdiffhistos[fsAndVar]->SetMarkerColor(kRed);
	  allfracdiffhistos[fsAndVar]->SetLineColor(kRed);

	  for(unsigned int ibin=1; ibin<=binsUInt[*varsToPlot]; ++ibin) {

	    double ndt = alldatahistos["h_data_"+fsAndVar]->GetBinContent(ibin);
	    double nmc = allmchistos["h_sumAllMc_"+fsAndVar]->GetBinContent(ibin);
	    double err_nmc = allmchistos["h_sumAllMc_"+fsAndVar]->GetBinError(ibin); 
	    //double err_tot = sqrt(err_nmc*err_nmc  + ndt); // correct: taking sqrt[sum(weights^2)] as MC error
	    double err_tot = sqrt(nmc + ndt); // not-so-correct: taking sqrt(N) as MC error
	    double chi2bin=0.;
	    double fracbin=0.;
	    double err_fracbin=0.;

	    //if(nmc>1e-6 || ndt>1e-6) {
	    //if(err_tot>1e-6) {
	    if(nmc>1e-6 && ndt>1e-6) {

	      chi2bin = (ndt-nmc)/err_tot;

	      if(nmc>1e-6) {
		fracbin = (ndt-nmc)/nmc;
		err_fracbin = fracbin*sqrt( 1./(chi2bin*chi2bin) + (err_nmc/nmc)*(err_nmc/nmc) );

		//std::cout << alldatahistos["h_data_"+fsAndVar]->GetBinCenter(ibin) << ", " << (ndt/normdt)/(nmc/normmc) << std::endl;
		std::cout << alldatahistos["h_data_"+fsAndVar]->GetBinCenter(ibin) << ", " << (ndt)/(nmc) << std::endl;

	      }

	    }

	    // allpoisdiffhistos[fsAndVar]->SetBinContent(ibin, chi2bin);
	    allfracdiffhistos[fsAndVar]->SetBinContent(ibin, fracbin);
	    allfracdiffhistos[fsAndVar]->SetBinError(ibin, err_fracbin);

	  } // end for(int ibin=0; ibin<=nbins+1; ++ibin)

	  TLine *ll=new TLine(firstBinsDouble[*varsToPlot], 0., lastBinsDouble[*varsToPlot], 0.);
	  ll->SetLineWidth(2);
	  ll->SetLineStyle(7);

	  //allpads2[fsAndVar]->cd();
	  //allpoisdiffhistos[fsAndVar]->Draw("axis");
	  // ll->Draw();
	  // allpoisdiffhistos[fsAndVar]->Draw("psame");
	  //allpads2[fsAndVar]->SetGridx(1);
	  //allpads2[fsAndVar]->SetGridy(1);

	  allpads3[fsAndVar]->cd();
	  allfracdiffhistos[fsAndVar]->Draw("axis");
	  ll->Draw();
	  allfracdiffhistos[fsAndVar]->Draw("pesame");
	  allpads3[fsAndVar]->SetGridx(1);
	  allpads3[fsAndVar]->SetGridy(1);


	} // end for(std::vector<TString>::iterator varsToPlot = variables.begin(); varsToPlot!=variables.end(); ++varsToPlot)
      } // end for(std::map<int, TString>::iterator finStat = finalStates.begin(); finStat!=finalStates.end(); ++finStat)
    } // end if(drawData)
  } // end if(doPlot)


  // -- // -- // -- // -- // -- // -- // -- // -- // -- // -- // -- // 
  // -- // -- // -- // -- // -- // -- // -- // -- // -- // -- // -- // 
  // -- // -- // -- // -- // -- // -- // -- // -- // -- // -- // -- // 

  return;

}


void fillPlots(std::map<TString, TH1F*> & histos, TString plotlab, double wght) {

  double thismass = getMass(l1_en+l2_en, l1_px+l2_px, l1_py+l2_py, l1_pz+l2_pz);
  histos[plotlab+"dileptMass"]->Fill(thismass, wght);
  if( fabs(thismass-91.1876)>10. ) return;

  histos[plotlab+"jetNumber"]->Fill(jnum, wght);
  if( jnum>0 ) return;

  // histos[plotlab+"extraLeptonNumber"]->Fill(ln, wght);
  // if( ln>0 ) return;

  // double thisD0redMet = getD0RedMet(l1_px, l1_py, l1_ptErr, l2_px, l2_py, l2_ptErr, htvec_px, htvec_py, met_pt[0], met_phi[0], cat);
  // histos[plotlab+"d0RedMet"]->Fill(thisD0redMet, wght);
  // if( thisD0redMet<60. ) return;

  return;

}

