/*
 *  $Date: 2012/06/25 08:02:09 $
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
#include "TDatime.h"
#include <iostream>
#include <iomanip>
#include <map>
#include <sys/stat.h>
#include <sys/types.h>

#include "treeVariables.C"
#include "usefulFunctions.C"
#include "efficiency.cpp"

using namespace std;

// 
// Function to fill plots
// 
void fillPlots(std::map<TString, TH1F*> &, TString, double wght = 1.0, double puWght = -1.0);
void printEvents(std::map<int, TString> &, 
		 std::map<TString, TH1F*> &, unsigned int, TString *,
		 std::map<TString, TH1F*> &, unsigned int, TString *, bool, 
		 std::map<TString, TH1F*> &, 
		 std::map<TString, TH1F*> &);


//
// Main function
//
void wz_dileptonMetAnalysis() {

  // int loaded = -1;
  // if( gSystem->IsFileInIncludePath("usefulFunctions_C.so") )
  //   loaded = gSystem->Load("usefulFunctions_C.so");
  // if(loaded < 0)
  //   gSystem->CompileMacro("usefulFunctions.C");

  // Set name of output folder 
  TDatime tdt;
  TString outputfolder("wz_plots_");
  outputfolder += tdt.GetYear();
  outputfolder += "-";
  TString tmpstr; tmpstr += tdt.GetMonth(); if( tmpstr.Length()==1 ) tmpstr.Prepend("0"); outputfolder += tmpstr; 
  outputfolder += "-";
  tmpstr.Clear(); tmpstr += tdt.GetDay(); if( tmpstr.Length()==1 ) tmpstr.Prepend("0"); outputfolder += tmpstr; 
  outputfolder += "_";
  tmpstr.Clear(); tmpstr += tdt.GetHour(); if( tmpstr.Length()==1 ) tmpstr.Prepend("0"); outputfolder += tmpstr; 
  outputfolder += "-";
  tmpstr.Clear(); tmpstr += tdt.GetMinute(); if( tmpstr.Length()==1 ) tmpstr.Prepend("0"); outputfolder += tmpstr; 
  outputfolder += "-";
  tmpstr.Clear(); tmpstr += tdt.GetSecond(); if( tmpstr.Length()==1 ) tmpstr.Prepend("0"); outputfolder += tmpstr; 
  std::cout << "Creating output folder \"" << outputfolder.Data() << "\"" << std::endl;
  gSystem->Exec( ("mkdir "+outputfolder).Data() );

  initializeGlobalVariables();

  // Samples location
  //TString indir="/home/daniele/Documents/Work/samples/2012-03-26_cmgTrees";
  //TString indir="/home/daniele/Documents/Work/samples/2012-04-03_cmgTrees";
  //TString indir="/home/daniele/Documents/Work/samples/2012-05-14_trees"; // last used
  TString indir="root://eoscms//eos/cms/store/cmst3/user/querten/12_04_14_HZZ2l2v_ntuples"; // for lxplus
  TString higgsindir="rfio:/castor/cern.ch/user/t/trocino/ZZllnn/InvisibleHiggs/2012-04-16_invisibleHiggs/trees"; // for lxplus
  TString basename_mc="MC_";
  TString basename_dt="Data_";
  TString base_mc=indir+"/"+basename_mc;
  TString higgsbase_mc=higgsindir+"/"+basename_mc;
  TString base_dt=indir+"/"+basename_dt;

  // Integrated luminosity
  float intLumi = 5035.;
  TString integrLumi = "5035";

  // Pile-up scenario in data
  std::vector<double> dataPuDistribution;
  double alldatapileups[] = {1.22825e+07, 5.33316e+07, 1.27819e+08, 2.21102e+08, 3.09325e+08, 3.74101e+08, 4.09049e+08, 4.17488e+08, 4.06878e+08, 3.84466e+08, 3.55412e+08, 3.22755e+08, 2.88141e+08, 2.52593e+08, 2.17008e+08, 1.82346e+08, 1.49605e+08, 1.19698e+08, 9.33203e+07, 7.08679e+07, 5.24176e+07, 3.77691e+07, 2.65207e+07, 1.81566e+07, 1.21265e+07, 7.90616e+06, 5.03513e+06, 3.13451e+06, 1.90872e+06, 1.1377e+06, 664248, 380148, 213407, 117609, 63687.5, 0.}; // pixel-based
  // double alldatapileups[] = {1.344651e+07, 5.90653e+07, 1.409027e+08, 2.413012e+08, 3.337449e+08, 3.98711e+08, 4.301064e+08, 4.32283e+08, 4.138202e+08, 3.82846e+08, 3.451637e+08, 3.043438e+08, 2.62555e+08, 2.213308e+08, 1.819826e+08, 1.456898e+08, 1.134134e+08, 8.577886e+07, 6.301239e+07, 4.495959e+07, 3.116904e+07, 2.100786e+07, 1.377588e+07, 8796407, 5474418, 3323776, 1970638, 1142040, 647538.6, 359547.2, 195673.2, 104459.9, 54745.15, 28185.57, 28005.55, 0.008}; // hf-based

  unsigned int nDataPu = sizeof(alldatapileups)/sizeof(double);
  for(unsigned int pui=0; pui<nDataPu; ++pui) {
    dataPuDistribution.push_back( alldatapileups[pui] );
  }

  // Labels for all standard processes
  TString alllabels[]={"zh105", "zh115", "zh125", "zh150", "zz", "wz", "ww", "tt", "t_s", "tbar_s", "t_t", "tbar_t", "t_tw", "tbar_tw", "zzx", "w", "z"};
  unsigned int nSmps=sizeof(alllabels)/sizeof(TString);

  // Labels for legend
  std::map<TString, TString> alllegendlabs;
  alllegendlabs["zh105"]="ZH(105) #rightarrow 2l+MET";
  alllegendlabs["zh115"]="ZH(115) #rightarrow 2l+MET";
  alllegendlabs["zh125"]="ZH(125) #rightarrow 2l+MET";
  alllegendlabs["zh150"]="ZH(150) #rightarrow 2l+MET";
  alllegendlabs["zz"]="ZZ #rightarrow 2l2#nu";
  alllegendlabs["zzx"]="ZZ #rightarrow X"; 
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
  alllegendlabs["z"]="Z #rightarrow ll";

  // Cross sections, process by process
  //double xSect[] = {0.01171, 5.9, 18.2, 43., 165., 3.19, 1.44, 41.92, 22.65, 7.87, 7.87, 31314., 3048.}; // old (for my own trees)
  //double xSect[] = {0.04001, 0.02981, 0.02231, 0.01171, 4.287, 0.856, 4.78, 165., 3.19, 1.44, 41.92, 22.65, 7.87, 7.87, 4.287, 31314., 3048.};
  double xSect[] = {0.04001, 0.02981, 0.02231, 0.01171, 5.9, 0.856, 4.78, 165., 3.19, 1.44, 41.92, 22.65, 7.87, 7.87, 5.9, 31314., 3048.};
  unsigned int nXsect = sizeof(xSect)/sizeof(double);
  if( nXsect != nSmps ) {
    std::cout << " *************************** ERROR ***************************" << endl;
    std::cout << "    Number of cross sections != number of samples! Check!"      << std::endl;
    return;
  }

  // Branching ratios, process by process
  double brFra[] = {1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.};
  unsigned int nBrFra = sizeof(brFra)/sizeof(double);
  if( nBrFra != nSmps ) {
    std::cout << " *************************** ERROR ***************************" << endl;
    std::cout << "    Number of branching ratios != number of samples! Check!"    << std::endl;
    return;
  }

  // Correction factors for PU re-weighting, process by process
  // zz	        CF = (77686 +- 278.722) / (78867.3 +- 304.917) = 0.985022 +- 0.00519546
  // wz	        CF = (205335 +- 453.139) / (208323 +- 494.938) = 0.985657 +- 0.00319611
  // ww	        CF = (172559 +- 415.402) / (176613 +- 458.705) = 0.977047 +- 0.00346001
  // tt	        CF = (98733 +- 314.218) / (101587 +- 348.506) = 0.971905 +- 0.00454799
  // t	        CF = (96185 +- 310.137) / (99457 +- 343.772) = 0.967101 +- 0.00457143
  // zzx	CF = (111301 +- 333.618) / (113072 +- 363.994) = 0.984336 +- 0.00432967
  // w	        CF = (449812 +- 670.68) / (483117 +- 763.878) = 0.931061 +- 0.00202346
  // z	        CF = (6.697e+06 +- 2587.86) / (6.90884e+06 +- 2847.21) = 0.969338 +- 0.000547617
  // zh105	CF = (15114 +- 122.939) / (15511.4 +- 136.691) = 0.974381 +- 0.0116853
  // zh115	CF = (15646 +- 125.084) / (15911.9 +- 135.785) = 0.98329 +- 0.011498
  // zh125	CF = (15403 +- 124.109) / (15713.3 +- 136.341) = 0.980251 +- 0.0116071
  // zh150	CF = (15917 +- 126.163) / (16212.4 +- 139.321) = 0.981777 +- 0.0114777
  // 
  // {"zh105", "zh115", "zh125", "zh150", "zz", "wz", "ww", "tt", "t_s", "tbar_s", "t_t", "tbar_t", "t_tw", "tbar_tw", "zzx", "w", "z"};
  double puCorrFact[] = {0.974381, 0.98329, 0.980251, 0.981777, 0.985022, 0.985657, 0.977047, 0.971905, 0.967101, 0.967101, 0.967101, 0.967101, 0.967101, 0.967101, 0.984336, 0.931061, 0.979652};
  unsigned int nPuCorrFact = sizeof(puCorrFact)/sizeof(double);
  if( nPuCorrFact != nSmps ) {
    std::cout << " *************************** ERROR ***************************" << endl;
    std::cout << "    Number of PU corr. factors != number of samples! Check!"    << std::endl;
    return;
  }


  // Colors and marker styles, process by process
  Color_t cols[]={kBlack, kBlack, kBlack, kBlack, kRed, kMagenta, kViolet-1, kBlue, kCyan, kCyan, kCyan, kCyan, kCyan, kCyan, kGreen, kSpring+3, kYellow-7};
  Style_t mrks[]={25, 25, 25, 25, 21, 22, 23, 24, 25, 25, 25, 25, 25, 25, 26, 27, 28, 29, 30, 31, 32};
  Style_t lines[]={1, 2, 3, 4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
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
  TString alllabelstouse[]={"zz", "wz", "ww", "tt", "t_s", "tbar_s", "t_t", "tbar_t", "t_tw", "tbar_tw", "zzx", "w", "z"};
  unsigned int nSmpsToUse=sizeof(alllabelstouse)/sizeof(TString);

  TString allhiggslabelstouse[]={};
  //TString allhiggslabelstouse[]={"zh105", "zh115", "zh125", "zh150"};
  unsigned int nHiggsSmpsToUse=sizeof(allhiggslabelstouse)/sizeof(TString);

  if( nSmps < (nSmpsToUse+nHiggsSmpsToUse) ) {
    std::cout << " ******************** ERROR ********************" << endl;
    std::cout << "    Not enough samples available, add samples!"   << std::endl;
    return;
  }

  // Maps associating each process to its tree, color, marker and line style
  std::map<TString, std::vector<TString> > allfiles;
  std::map<TString, TChain*> alltrees;
  std::map<TString, TChain*> alldatatrees;
  std::map<TString, double>  allgenev;
  std::map<TString, std::vector<double> >  allpileups;
  std::map<TString, double>  allxsect;
  std::map<TString, double>  allbrfra;
  std::map<TString, double>  allpucfs;
  std::map<TString, Color_t> allcols;
  std::map<TString, Style_t> allmrks;
  std::map<TString, Style_t> alllines;

  // Fill maps with files
  allfiles["zh105"].push_back(higgsbase_mc+"ZH105.root");
  allfiles["zh115"].push_back(higgsbase_mc+"ZH115.root");
  allfiles["zh125"].push_back(higgsbase_mc+"ZH125.root");
  allfiles["zh150"].push_back(higgsbase_mc+"ZH150.root");
  // allfiles["zh105"].push_back(base_mc+"ZH105.root");
  // allfiles["zh115"].push_back(base_mc+"ZH115.root");
  // allfiles["zh125"].push_back(base_mc+"ZH125.root");
  // allfiles["zh150"].push_back(base_mc+"ZH150.root");
  allfiles["zz"].push_back(base_mc+"ZZ_0.root"); allfiles["zz"].push_back(base_mc+"ZZ_1.root");
  allfiles["wz"].push_back(base_mc+"WZ_0.root"); allfiles["wz"].push_back(base_mc+"WZ_1.root");
  allfiles["ww"].push_back(base_mc+"WW.root");
  allfiles["tt"].push_back(base_mc+"TTJets.root");
  allfiles["t_tw"].push_back(base_mc+"SingleT_tW.root"); allfiles["tbar_tw"].push_back(base_mc+"SingleTbar_tW.root");
  allfiles["t_t"].push_back(base_mc+"SingleT_t.root");   allfiles["tbar_t"].push_back(base_mc+"SingleTbar_t.root");
  allfiles["t_s"].push_back(base_mc+"SingleT_s.root");   allfiles["tbar_s"].push_back(base_mc+"SingleTbar_s.root");
  allfiles["zzx"].push_back(base_mc+"ZZ_0.root"); allfiles["zzx"].push_back(base_mc+"ZZ_1.root");
  allfiles["w"].push_back(base_mc+"WJetsToLNu.root"); 
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
    // Fill maps with PU correction factors
    allpucfs[alllabels[k]]=puCorrFact[k];

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

  // // Data samples
  // TChain *t_data=new TChain("evAnalyzer/data");
  // // DoubleMu
  // t_data->Add(base_dt+"DoubleMu2011A_0.root"); t_data->Add(base_dt+"DoubleMu2011A_1.root"); 
  // t_data->Add(base_dt+"DoubleMu2011B_0.root"); t_data->Add(base_dt+"DoubleMu2011B_1.root"); 
  // // // DoubleElectron
  // t_data->Add(base_dt+"DoubleElectron2011A_0.root"); t_data->Add(base_dt+"DoubleElectron2011A_1.root"); 
  // t_data->Add(base_dt+"DoubleElectron2011B_0.root"); t_data->Add(base_dt+"DoubleElectron2011B_1.root"); 
  // // // MuUG
  // t_data->Add(base_dt+"MuEG2011A_0.root"); t_data->Add(base_dt+"MuEG2011A_1.root"); 
  // t_data->Add(base_dt+"MuEG2011B_0.root"); t_data->Add(base_dt+"MuEG2011B_1.root"); 

  // Data samples
  // DoubleMu
  alldatatrees["mm"] = new TChain("evAnalyzer/data");
  alldatatrees["mm"]->Add(base_dt+"DoubleMu2011A_0.root"); alldatatrees["mm"]->Add(base_dt+"DoubleMu2011A_1.root"); 
  alldatatrees["mm"]->Add(base_dt+"DoubleMu2011B_0.root"); alldatatrees["mm"]->Add(base_dt+"DoubleMu2011B_1.root"); 
  // DoubleElectron
  alldatatrees["ee"] = new TChain("evAnalyzer/data");
  alldatatrees["ee"]->Add(base_dt+"DoubleElectron2011A_0.root"); alldatatrees["ee"]->Add(base_dt+"DoubleElectron2011A_1.root"); 
  alldatatrees["ee"]->Add(base_dt+"DoubleElectron2011B_0.root"); alldatatrees["ee"]->Add(base_dt+"DoubleElectron2011B_1.root"); 
  // MuEG
  alldatatrees["em"] = new TChain("evAnalyzer/data");
  alldatatrees["em"]->Add(base_dt+"MuEG2011A_0.root"); alldatatrees["em"]->Add(base_dt+"MuEG2011A_1.root"); 
  alldatatrees["em"]->Add(base_dt+"MuEG2011B_0.root"); alldatatrees["em"]->Add(base_dt+"MuEG2011B_1.root"); 

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
  //  N.B.: variables used for selection should be added with "addPrintVariable(...)";
  //        "auxiliary" variables (not used for selection) should be added with "addVariable(...)"
  // 
  //   Arguments (all TString): 
  //    - variable "label"
  //    - x axis title [unit]
  //    - y axis title; if empty: "Events/<bin+unit>" 
  //    - string containing: n. bins, first bin, last bin
  //
  //           Variable             Title X                                     Y   N    x0     xN           Plot   Print
  addVariable( "dileptMass",        "Dilepton invariant mass [GeV/c^{2}]",      "", 60,  60.,   120.,        false, false );
  addVariable( "dileptPt",          "Dilepton p_{T} [GeV/c]",                   "", 60,  0.,    120.,        false, false );
  addVariable( "jetNumber",         "PFJet number (p_{T} > 30 GeV/c)",          "", 6,   -0.5,  5.5,         false, false );
  addVariable( "cmsIndMinRedMet",   "CMS reduced MET [GeV]",                    "", 100, 0.,    300.,        false, false ); 
  //addVariable( "d0RedMet",          "D0 reduced MET [GeV]",                     "", 100, 0.,    300.,        false, false ); 
  addVariable( "metPtBalance",      "PF MET/p_{T}(Z)",                          "", 30,  0.,    3.,          false, false ); 
  addVariable( "deltaPhiJetMet",    "#Delta#phi(jet,MET) [rad]",                "", 18., 0.,    3.141592654, false, false ); 
  addVariable( "jetCsv",            "CSV discriminator (PFJet p_{T} > 20 GeV)", "", 22,  -1.1,  1.1,         false, false );
  //addVariable( "extraLeptonNumber", "Additional lepton number (loose sel.)",    "", 4,   -0.5,  3.5,         false, true  );

  // WZ background
  // addVariable( "extraLeptonGenId",  "3^{rd} lepton gen. ID",              "",           7,   9.5,  16.5, false, false );
  // addVariable( "extraMuonGenPt",    "3^{rd} muon gen. p_{T} [GeV/c]",     "",           20,  0.,   100., false, false );
  // addVariable( "extraMuonEffPt",    "3^{rd} muon gen. p_{T} [GeV/c]",     "Efficiency", 20,  0.,   100., false, false );
  // addVariable( "extraMuonGenEta",   "3^{rd} muon gen. #eta",              "",           24,  -6.,  6.,   false, false );
  // addVariable( "extraMuonEffEta",   "3^{rd} muon gen. #eta",              "Efficiency", 24,  -6.,  6.,   false, false );
  // addVariable( "extraEleGenPt",     "3^{rd} electron gen. p_{T} [GeV/c]", "",           20,  0.,   100., false, false );
  // addVariable( "extraEleEffPt",     "3^{rd} electron gen. p_{T} [GeV/c]", "Efficiency", 20,  0.,   100., false, false );
  // addVariable( "extraEleGenEta",    "3^{rd} electron gen. #eta",          "",           24,  -6.,  6.,   false, false );
  // addVariable( "extraEleEffEta",    "3^{rd} electron gen. #eta",          "Efficiency", 24,  -6.,  6.,   false, false );

  // addVariable( "extraLeptonGenIdEM","3^{rd} lepton gen. ID e-#mu",        "",           1,   0.,   2.,   false, false );
  // addVariable( "extraLeptonGenIdAll","3^{rd} lepton gen. ID all",         "",           1,   0.,   2.,   false, false );
  // addVariable( "extraMuonGenAll",   "3^{rd} muon gen.",                   "",           1,   0.,   2.,   false, false );
  // addVariable( "extraMuonGenAcc",   "3^{rd} muon gen. accept.",           "",           1,   0.,   2.,   false, false );
  // addVariable( "extraEleGenAll",    "3^{rd} electron gen.",               "",           1,   0.,   2.,   false, false );
  // addVariable( "extraEleGenAcc",    "3^{rd} electron gen. accept.",       "",           1,   0.,   2.,   false, false );
  // // addVariable( "extraLepton0",      "0 extra leptons",                    "",           1,   0.,   2.,   false, true  );
  // // addVariable( "extraLepton1",      "1 extra lepton",                     "",           1,   0.,   2.,   false, true  );

  addVariable( "extraLeptonNumber",      "3^{rd} lepton number",               "",           4,   -0.5, 3.5,  true,  false );
  addVariable( "extraLeptonNumberAlpha", "3^{rd} lepton number",               "",           4,   -0.5, 3.5,  true,  false );

  addVariable( "thirdLeptMetTransMass", "3^{rd} lepton-MET trans. mass [GeV/c^{2}]", "", 25, 180., 280., true,  false );
  addVariable( "thirdLeptMetTransMassAlpha", "3^{rd} lepton-MET trans. mass [GeV/c^{2}]", "", 25, 180., 280., true,  false );

  addVariable( "finalSelection", "Final selection", "", 1, 0., 2., false, true );

  unsigned int numberCuts = allVarsToPrint.size();

  // One extra plots: the cut-flow
  if(numberCuts>0)
    addVariable( "cutFlow", "", "Events", numberCuts, 0.5, float(numberCuts)+0.5, true, false );

  //
  // ::::::::::::::::::::::::::::::::::::::::::::::::::::
  // ::::::::::::::::::::::::::::::::::::::::::::::::::::
  // ::::::::::::::::::::::::::::::::::::::::::::::::::::
  //

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
  bool doCutTable=true;

  double canvx(500.), canvy(500.);
  if(binbybinComp && drawData) {
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
    std::map<TString, TH1F*> allhistos;
    std::map<TString, TH1F*> allhiggshistos;
    std::map<TString, TH1F*> allmchistos;
    std::map<TString, TH1F*> alldatahistos;
    //std::map<TString, TH1F*> allpoisdiffhistos;
    std::map<TString, TH1F*> allfracdiffhistos;

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
	      TString plotIdx = plotname+finStat->second+"_"+(*varsToPlot);
	      allhistos[plotIdx] = new TH1F( plotIdx.Data(), plotIdx.Data(), binsUInt[*varsToPlot], 
					     firstBinsDouble[*varsToPlot], lastBinsDouble[*varsToPlot] );
	      allhistos[plotIdx]->Sumw2();
	      j_ah += 1;
	    }
	    // Only for cut-flow plot
	    for(unsigned int h=0; h<allVarsToPrint.size(); ++h) {
	      allhistos[plotname+finStat->second+"_cutFlow"]->GetXaxis()->SetBinLabel(h+1, allVarsToPrint[h]);
	    }
	  }
	}

	unsigned int nEvt = alltrees[alllabelstouse[j]]->GetEntries();
	std::cout << " - " << alllabelstouse[j].Data() << "  " << nEvt << std::endl;

	initializeTreeVariables();
	attachToTree(alltrees[alllabelstouse[j]]);

	// General weight
	float genWeight = intLumi * allxsect[alllabelstouse[j]] * allbrfra[alllabelstouse[j]] * allpucfs[alllabelstouse[j]] / allgenev[alllabelstouse[j]];

	// 
	// Loop over events
	// 
	for(unsigned int iEvt=0; iEvt<nEvt; ++iEvt) {
	  alltrees[alllabelstouse[j]]->GetEntry(iEvt);

	  // Split ZZ into ZZ->2l2n and ZZ->X
	  // 
	  // - Select only ZZ -> 2l2n
	  if( alllabelstouse[j].CompareTo("zz") == 0 ) {
	    if( isZZ2l2nu(mccat) == false ) {
	      continue;
	    }
	  }
	  else if( alllabelstouse[j].CompareTo("zzx") == 0 ) {
	    if( isZZ2l2nu(mccat) == true ) {
	      continue;
	    }
	  }

	  // if( finalStates.count(cat)>0 &&     // one of the chosen categories (1: mm; 2: ee; 3: em)
	  //     l1_id * l2_id < 0          ) {  // same charge
	  if( finalStates.count(cat)>0 ) {  // one of the chosen categories (1: mm; 2: ee; 3: em)

	    TString plotIdx = plotname+finalStates[cat]+"_";
	    // Event weight
	    float evtWeight = genWeight * getPuWeights(ngenITpu);

	    fillPlots(allhistos, plotIdx, evtWeight);
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

	    if( allmchistos.count("h_sumAllMc_"+fsAndVar)==0 ) {
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
	  // Only for cut-flow plot
	  for(unsigned int h=0; h<allVarsToPrint.size(); ++h) {
	    allhiggshistos[plotname+finStat->second+"_cutFlow"]->GetXaxis()->SetBinLabel(h+1, allVarsToPrint[h]);
	  }
	}

	unsigned int nEvt = alltrees[allhiggslabelstouse[j]]->GetEntries();
	std::cout << " - " << allhiggslabelstouse[j].Data() << "  " << nEvt << std::endl;

	initializeTreeVariables();
	attachToTree(alltrees[allhiggslabelstouse[j]]);

	// General weight
	float genWeight = intLumi * allxsect[allhiggslabelstouse[j]] * allbrfra[allhiggslabelstouse[j]] * allpucfs[allhiggslabelstouse[j]] / allgenev[allhiggslabelstouse[j]];

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
      for(std::map<int, TString>::iterator finStat = finalStates.begin(); finStat!=finalStates.end(); ++finStat) {
	for(std::vector<TString>::iterator varsToPlot = variables.begin(); varsToPlot!=variables.end(); ++varsToPlot) {
	  TString plotIdx = "h_data_"+finStat->second+"_"+(*varsToPlot);
	  alldatahistos[plotIdx] = new TH1F( plotIdx.Data(), plotIdx.Data(), binsUInt[*varsToPlot], 
					     firstBinsDouble[*varsToPlot], lastBinsDouble[*varsToPlot] );
	}
	// Only for cut-flow plot
	for(unsigned int h=0; h<allVarsToPrint.size(); ++h) {
	  alldatahistos["h_data_"+finStat->second+"_cutFlow"]->GetXaxis()->SetBinLabel(h+1, allVarsToPrint[h]);
	}
      }

      // initializeTreeVariables();
      // attachToTree(t_data);

      // unsigned int nEvt = t_data->GetEntries();
      // std::cout << " - data " << nEvt << std::endl;

      // 
      // Loop over events
      // 
      for(std::map<int, TString>::iterator finStat = finalStates.begin(); finStat!=finalStates.end(); ++finStat) {
	TChain *this_t_data = alldatatrees[(finStat->second).Data()];
	initializeTreeVariables();
	attachToTree( this_t_data );

	unsigned int nEvt = this_t_data->GetEntries();
	std::cout << " - data (" << (finStat->second).Data() << ") " << nEvt << std::endl;

	TString plotIdx = "h_data_"+finStat->second+"_";

	for(unsigned int iEvt=0; iEvt<nEvt; ++iEvt) {
	  this_t_data->GetEntry(iEvt);
	  if(cat != finStat->first) continue; 
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

    // // Only for WZ background
    // allhistos["h_wz_mm_extraMuonEffPt"]->Divide(  allhistos["h_wz_mm_extraMuonGenPt"] ); 
    // allhistos["h_wz_mm_extraMuonEffEta"]->Divide( allhistos["h_wz_mm_extraMuonGenEta"] ); 
    // allhistos["h_wz_ee_extraMuonEffPt"]->Divide(  allhistos["h_wz_ee_extraMuonGenPt"] ); 
    // allhistos["h_wz_ee_extraMuonEffEta"]->Divide( allhistos["h_wz_ee_extraMuonGenEta"] ); 
    // allhistos["h_wz_mm_extraEleEffPt"]->Divide(  allhistos["h_wz_mm_extraEleGenPt"] ); 
    // allhistos["h_wz_mm_extraEleEffEta"]->Divide( allhistos["h_wz_mm_extraEleGenEta"] ); 
    // allhistos["h_wz_ee_extraEleEffPt"]->Divide(  allhistos["h_wz_ee_extraEleGenPt"] ); 
    // allhistos["h_wz_ee_extraEleEffEta"]->Divide( allhistos["h_wz_ee_extraEleGenEta"] ); 

    // N.B. only "allVarsToPlot" is used, not all "variables"!!!
    // 
    for(std::map<int, TString>::iterator finStat = finalStates.begin(); finStat!=finalStates.end(); ++finStat) {
      for(std::vector<TString>::iterator varsToPlot = allVarsToPlot.begin(); varsToPlot!=allVarsToPlot.end(); ++varsToPlot) {
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
	  approxToN(onebin, 0, -2);
	  yTitles[*varsToPlot] += onebin;
	  if( yTitles[*varsToPlot].Length()>15 ) {
	    int lastFigIdx = yTitles[*varsToPlot].Length()-1;
	    char lastFig = (yTitles[*varsToPlot].Data())[lastFigIdx]; 
	    char lastBut1Fig = (yTitles[*varsToPlot].Data())[lastFigIdx-1]; 
	    TString lastBut1FigStr(lastBut1Fig);
	    yTitles[*varsToPlot] = yTitles[*varsToPlot].Strip(TString::kTrailing, lastFig);
	    while( yTitles[*varsToPlot].EndsWith(lastBut1FigStr.Data()) ) 
	      yTitles[*varsToPlot] = yTitles[*varsToPlot].Strip(TString::kTrailing, lastBut1Fig);
	  }
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
		//std::cout << alldatahistos["h_data_"+fsAndVar]->GetBinCenter(ibin) << ", " << (ndt)/(nmc) << std::endl;
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

	  // Save plots
	  
	  allcanvas[fsAndVar]->SaveAs( (outputfolder+"/"+allcanvas[fsAndVar]->GetName()+".png").Data() );
	  allcanvas[fsAndVar]->SaveAs( (outputfolder+"/"+allcanvas[fsAndVar]->GetName()+".root").Data() );

	} // end for(std::vector<TString>::iterator varsToPlot = allVarsToPlot.begin(); varsToPlot!=allVarsToPlot.end(); ++varsToPlot)
      } // end for(std::map<int, TString>::iterator finStat = finalStates.begin(); finStat!=finalStates.end(); ++finStat)
    } // end if(drawData)

    if(doCutTable) {
      printEvents(finalStates, 
		  allhiggshistos, nHiggsSmpsToUse, allhiggslabelstouse,
		  allhistos,      nSmpsToUse,      alllabelstouse,     mergeTopSamples, 
		  allmchistos, 
		  alldatahistos);
    }

    std::cout << std::endl << std::endl;

    // 
    // Print out result of WZ measurement
    // 
    /*
    double var_wz_mm_m_all(0.), err_wz_mm_m_all(0.), var_wz_mm_m_pss(0.), err_wz_mm_m_pss(0.), var_wz_mm_m_acc(0.), err_wz_mm_m_acc(0.); 
    double var_wz_mm_e_all(0.), err_wz_mm_e_all(0.), var_wz_mm_e_pss(0.), err_wz_mm_e_pss(0.), var_wz_mm_e_acc(0.), err_wz_mm_e_acc(0.); 
    double var_wz_ee_m_all(0.), err_wz_ee_m_all(0.), var_wz_ee_m_pss(0.), err_wz_ee_m_pss(0.), var_wz_ee_m_acc(0.), err_wz_ee_m_acc(0.); 
    double var_wz_ee_e_all(0.), err_wz_ee_e_all(0.), var_wz_ee_e_pss(0.), err_wz_ee_e_pss(0.), var_wz_ee_e_acc(0.), err_wz_ee_e_acc(0.); 

    var_wz_mm_m_all = allhistos["h_wz_mm_extraMuonGenAll"]->IntegralAndError(0, -1, err_wz_mm_m_all);
    var_wz_mm_m_pss = allhistos["h_wz_mm_extraMuonGenAcc"]->IntegralAndError(0, -1, err_wz_mm_m_pss);
    var_wz_mm_m_acc = var_wz_mm_m_pss/var_wz_mm_m_all;
    err_wz_mm_m_acc = var_wz_mm_m_acc * sqrt( pow( (err_wz_mm_m_pss/var_wz_mm_m_pss), 2) + pow( (err_wz_mm_m_all/var_wz_mm_m_all), 2) );

    std::cout << "Channel mm:" << std::endl;
    std::cout << " - acceptance third mu = ";
    std::cout <<  var_wz_mm_m_acc << " +- " << err_wz_mm_m_acc << " (" << err_wz_mm_m_acc/var_wz_mm_m_acc*100 << "%)" << std::endl;

    var_wz_mm_e_all = allhistos["h_wz_mm_extraEleGenAll"]->IntegralAndError(0, -1, err_wz_mm_e_all);
    var_wz_mm_e_pss = allhistos["h_wz_mm_extraEleGenAcc"]->IntegralAndError(0, -1, err_wz_mm_e_pss);
    var_wz_mm_e_acc = var_wz_mm_e_pss/var_wz_mm_e_all;
    err_wz_mm_e_acc = var_wz_mm_e_acc * sqrt( pow( (err_wz_mm_e_pss/var_wz_mm_e_pss), 2) + pow( (err_wz_mm_e_all/var_wz_mm_e_all), 2) );

    std::cout << " - acceptance third e  = ";
    std::cout <<  var_wz_mm_e_acc << " +- " << err_wz_mm_e_acc << " (" << err_wz_mm_e_acc/var_wz_mm_e_acc*100 << "%)" << std::endl;

    var_wz_ee_m_all = allhistos["h_wz_ee_extraMuonGenAll"]->IntegralAndError(0, -1, err_wz_ee_m_all);
    var_wz_ee_m_pss = allhistos["h_wz_ee_extraMuonGenAcc"]->IntegralAndError(0, -1, err_wz_ee_m_pss);
    var_wz_ee_m_acc = var_wz_ee_m_pss/var_wz_ee_m_all;
    err_wz_ee_m_acc = var_wz_ee_m_acc * sqrt( pow( (err_wz_ee_m_pss/var_wz_ee_m_pss), 2) + pow( (err_wz_ee_m_all/var_wz_ee_m_all), 2) );

    std::cout << "Channel ee:" << std::endl;
    std::cout << " - acceptance third mu = ";
    std::cout <<  var_wz_ee_m_acc << " +- " << err_wz_ee_m_acc << " (" << err_wz_ee_m_acc/var_wz_ee_m_acc*100 << "%)" << std::endl;

    var_wz_ee_e_all = allhistos["h_wz_ee_extraEleGenAll"]->IntegralAndError(0, -1, err_wz_ee_e_all);
    var_wz_ee_e_pss = allhistos["h_wz_ee_extraEleGenAcc"]->IntegralAndError(0, -1, err_wz_ee_e_pss);
    var_wz_ee_e_acc = var_wz_ee_e_pss/var_wz_ee_e_all;
    err_wz_ee_e_acc = var_wz_ee_e_acc * sqrt( pow( (err_wz_ee_e_pss/var_wz_ee_e_pss), 2) + pow( (err_wz_ee_e_all/var_wz_ee_e_all), 2) );

    std::cout << " - acceptance third e  = ";
    std::cout <<  var_wz_ee_e_acc << " +- " << err_wz_ee_e_acc << " (" << err_wz_ee_e_acc/var_wz_ee_e_acc*100 << "%)" << std::endl;

    double var_wz_mm_emu(0.), err_wz_mm_emu(0.), var_wz_mm_all(0.), err_wz_mm_all(0.), var_wz_mm_idemu(0.), err_wz_mm_idemu(0.); 
    double var_wz_ee_emu(0.), err_wz_ee_emu(0.), var_wz_ee_all(0.), err_wz_ee_all(0.), var_wz_ee_idemu(0.), err_wz_ee_idemu(0.); 
    var_wz_mm_emu = allhistos["h_wz_mm_extraLeptonGenIdEM"]->IntegralAndError(0, -1, err_wz_mm_emu);
    var_wz_mm_all = allhistos["h_wz_mm_extraLeptonGenIdAll"]->IntegralAndError(0, -1, err_wz_mm_all);
    var_wz_mm_idemu = var_wz_mm_emu/var_wz_mm_all;
    err_wz_mm_idemu = var_wz_mm_idemu * sqrt( pow( (err_wz_mm_emu/var_wz_mm_emu), 2) + pow( (err_wz_mm_all/var_wz_mm_all), 2) );

    std::cout << "Channel mm:" << std::endl;
    std::cout << " - fraction of e+mu = ";
    std::cout <<  var_wz_mm_idemu << " +- " << err_wz_mm_idemu << " (" << err_wz_mm_idemu/var_wz_mm_idemu*100 << "%)" << std::endl;

    var_wz_ee_emu = allhistos["h_wz_ee_extraLeptonGenIdEM"]->IntegralAndError(0, -1, err_wz_ee_emu);
    var_wz_ee_all = allhistos["h_wz_ee_extraLeptonGenIdAll"]->IntegralAndError(0, -1, err_wz_ee_all);
    var_wz_ee_idemu = var_wz_ee_emu/var_wz_ee_all;
    err_wz_ee_idemu = var_wz_ee_idemu * sqrt( pow( (err_wz_ee_emu/var_wz_ee_emu), 2) + pow( (err_wz_ee_all/var_wz_ee_all), 2) );

    std::cout << "Channel ee:" << std::endl;
    std::cout << " - fraction of e+mu = ";
    std::cout <<  var_wz_ee_idemu << " +- " << err_wz_ee_idemu << " (" << err_wz_ee_idemu/var_wz_ee_idemu*100 << "%)" << std::endl;
    */


    std::cout << "MC:" << std::endl;
    // --- MM ---
    std::cout << " - Channel mm:" << std::endl;
    double wz_mm_bin0 = allhistos["h_wz_mm_extraLeptonNumber"]->GetBinContent(1);
    double err_wz_mm_bin0 = allhistos["h_wz_mm_extraLeptonNumber"]->GetBinError(1);
    double wz_mm_bin1 = allhistos["h_wz_mm_extraLeptonNumber"]->GetBinContent(2);
    double err_wz_mm_bin1 = allhistos["h_wz_mm_extraLeptonNumber"]->GetBinError(2);
    double dt_mm_bin1 = alldatahistos["h_data_mm_extraLeptonNumber"]->GetBinContent(2);
    double err_dt_mm_bin1 = alldatahistos["h_data_mm_extraLeptonNumber"]->GetBinError(2);

    std::cout << " 0 lept.: " << wz_mm_bin0 << " +- "  << err_wz_mm_bin0 << std::endl;
    std::cout << " 1 lept.: " << wz_mm_bin1 << " +- "  << err_wz_mm_bin1 << std::endl;
    std::cout << std::endl;
    double wz_mm_sf = wz_mm_bin0/wz_mm_bin1;
    double relerr_wz_mm_sf = sqrt( pow(err_wz_mm_bin0/wz_mm_bin0, 2) + pow(err_wz_mm_bin1/wz_mm_bin1, 2) );
    double err_wz_mm_sf = wz_mm_sf * relerr_wz_mm_sf; 

    std::cout << " scale factor = (0 lept.)/(1 lept.) = " << wz_mm_sf << " +- "  
	      << err_wz_mm_sf << " (" << relerr_wz_mm_sf*100 << "%)" << std::endl;
    std::cout << std::endl;

    // --- EE ---
    std::cout << std::endl;
    std::cout << " - Channel ee:" << std::endl;
    double wz_ee_bin0 = allhistos["h_wz_ee_extraLeptonNumber"]->GetBinContent(1);
    double err_wz_ee_bin0 = allhistos["h_wz_ee_extraLeptonNumber"]->GetBinError(1);
    double wz_ee_bin1 = allhistos["h_wz_ee_extraLeptonNumber"]->GetBinContent(2);
    double err_wz_ee_bin1 = allhistos["h_wz_ee_extraLeptonNumber"]->GetBinError(2);
    double dt_ee_bin1 = alldatahistos["h_data_ee_extraLeptonNumber"]->GetBinContent(2);
    double err_dt_ee_bin1 = alldatahistos["h_data_ee_extraLeptonNumber"]->GetBinError(2);

    std::cout << " 0 lept.: " << wz_ee_bin0 << " +- "  << err_wz_ee_bin0 << std::endl;
    std::cout << " 1 lept.: " << wz_ee_bin1 << " +- "  << err_wz_ee_bin1 << std::endl;
    std::cout << std::endl;
    double wz_ee_sf = wz_ee_bin0/wz_ee_bin1;
    double relerr_wz_ee_sf = sqrt( pow(err_wz_ee_bin0/wz_ee_bin0, 2) + pow(err_wz_ee_bin1/wz_ee_bin1, 2) );
    double err_wz_ee_sf = wz_ee_sf * relerr_wz_ee_sf; 

    std::cout << " scale factor = (0 lept.)/(1 lept.) = " << wz_ee_sf << " +- "  
	      << err_wz_ee_sf << " (" << relerr_wz_ee_sf*100 << "%)" << std::endl;
    std::cout << std::endl;

    std::cout << "Data (full selection, full SF):" << std::endl;
    // --- MM ---
    std::cout << " - Channel mm:" << std::endl;

    double dt_mm_bin0 = dt_mm_bin1 * wz_mm_sf;
    double err_dt_mm_bin0 = err_dt_mm_bin1 * wz_mm_sf;
    double syserr_dt_mm_bin0 = dt_mm_bin0 * relerr_wz_mm_sf;

    std::cout << " * mm data in 1 lept. bin (with full selection): " 
	      << dt_mm_bin1 << " +- " << err_dt_mm_bin1 << std::endl;
    std::cout << " * Expected mm data in 0 lept. bin (with full SF): " << std::endl;
    std::cout << "    " << dt_mm_bin0 << " +- " << err_dt_mm_bin0 << " (stat) +- "
	      << syserr_dt_mm_bin0 << " (syst)" << std::endl;

    // --- EE ---
    std::cout << " - Channel ee:" << std::endl;

    double dt_ee_bin0 = dt_ee_bin1 * wz_ee_sf;
    double err_dt_ee_bin0 = err_dt_ee_bin1 * wz_ee_sf;
    double syserr_dt_ee_bin0 = dt_ee_bin0 * relerr_wz_ee_sf;

    std::cout << " * ee data in 1 lept. bin (with full selection): " 
	      << dt_ee_bin1 << " +- " << err_dt_ee_bin1 << std::endl;
    std::cout << " * Expected ee data in 0 lept. bin (with full SF): " << std::endl;
    std::cout << "    " << dt_ee_bin0 << " +- " << err_dt_ee_bin0 << " (stat) +- "
	      << syserr_dt_ee_bin0 << " (syst)" << std::endl;

    std::cout << std::endl;

    /*
    // 
    // NOOOOOO
    // WRONG!!!!!!!!!
    // 
    std::cout << std::endl;
    std::cout << "Scale factor data/MC (with looser selection, higher statistics)" << std::endl;
    // --- MM ---
    std::cout << " - Channel mm:" << std::endl;
    double mc_ls_mm_bin1 = allmchistos["h_sumAllMc_mm_extraLeptonNumberAlpha"]->GetBinContent(2);
    double err_mc_ls_mm_bin1 = allmchistos["h_sumAllMc_mm_extraLeptonNumberAlpha"]->GetBinError(2);
    double dt_ls_mm_bin1 = alldatahistos["h_data_mm_extraLeptonNumberAlpha"]->GetBinContent(2);
    double err_dt_ls_mm_bin1 = alldatahistos["h_data_mm_extraLeptonNumberAlpha"]->GetBinError(2);

    double dtToMc_mm_sf = dt_ls_mm_bin1/mc_ls_mm_bin1;
    double relerr_dtToMc_mm_sf = sqrt( pow(err_mc_ls_mm_bin1/mc_ls_mm_bin1, 2) + pow(err_dt_ls_mm_bin1/dt_ls_mm_bin1, 2) );
    double err_dtToMc_mm_sf = dtToMc_mm_sf * relerr_dtToMc_mm_sf; 

    std::cout << " Data/MC scale factor = (Data, 1 lept.)/(MC, 1 lept.) = " << dtToMc_mm_sf << " +- "  
	      << err_dtToMc_mm_sf << " (" << relerr_dtToMc_mm_sf*100 << "%)" << std::endl;
    std::cout << std::endl;

    // --- EE ---
    std::cout << " - Channel ee:" << std::endl;
    double mc_ls_ee_bin1 = allmchistos["h_sumAllMc_ee_extraLeptonNumberAlpha"]->GetBinContent(2);
    double err_mc_ls_ee_bin1 = allmchistos["h_sumAllMc_ee_extraLeptonNumberAlpha"]->GetBinError(2);
    double dt_ls_ee_bin1 = alldatahistos["h_data_ee_extraLeptonNumberAlpha"]->GetBinContent(2);
    double err_dt_ls_ee_bin1 = alldatahistos["h_data_ee_extraLeptonNumberAlpha"]->GetBinError(2);

    double dtToMc_ee_sf = dt_ls_ee_bin1/mc_ls_ee_bin1;
    double relerr_dtToMc_ee_sf = sqrt( pow(err_mc_ls_ee_bin1/mc_ls_ee_bin1, 2) + pow(err_dt_ls_ee_bin1/dt_ls_ee_bin1, 2) );
    double err_dtToMc_ee_sf = dtToMc_ee_sf * relerr_dtToMc_ee_sf; 

    std::cout << " Data/MC scale factor = (Data, 1 lept.)/(MC, 1 lept.) = " << dtToMc_ee_sf << " +- "  
	      << err_dtToMc_ee_sf << " (" << relerr_dtToMc_ee_sf*100 << "%)" << std::endl;
    std::cout << std::endl;

    // ------------------------------------

    std::cout << std::endl;
    std::cout << "Data (full selection):" << std::endl;
    // --- MM ---
    std::cout << " - Channel mm:" << std::endl;
    double dt_mm_bin1 = alldatahistos["h_data_mm_extraLeptonNumber"]->GetBinContent(2);
    double err_dt_mm_bin1 = alldatahistos["h_data_mm_extraLeptonNumber"]->GetBinError(2);

    std::cout << " 1 lept.: " << dt_mm_bin1 << " +- "  << err_dt_mm_bin1 << std::endl;

    double dt_mm_bin0 = dt_mm_bin1 * wz_mm_sf * dtToMc_mm_sf;
    double err_dt_mm_bin0 = err_dt_mm_bin1 * wz_mm_sf * dtToMc_mm_sf;
    double syserr_dt_mm_bin0 = dt_mm_bin0 * sqrt( pow( relerr_wz_mm_sf, 2 ) + pow( relerr_dtToMc_mm_sf, 2 ) );

    std::cout << " Expected mm data in 0 lept. bin: " << std::endl;
    std::cout << "  " << dt_mm_bin0 << " +- " << err_dt_mm_bin0 << " (stat) +- "
	      << syserr_dt_mm_bin0 << " (syst)" << std::endl;

    // --- EE ---
    std::cout << std::endl;
    std::cout << " - Channel ee:" << std::endl;
    double dt_ee_bin1 = alldatahistos["h_data_ee_extraLeptonNumber"]->GetBinContent(2);
    double err_dt_ee_bin1 = alldatahistos["h_data_ee_extraLeptonNumber"]->GetBinError(2);

    std::cout << " 1 lept.: " << dt_ee_bin1 << " +- "  << err_dt_ee_bin1 << std::endl;

    double dt_ee_bin0 = dt_ee_bin1 * wz_ee_sf * dtToMc_ee_sf;
    double err_dt_ee_bin0 = err_dt_ee_bin1 * wz_ee_sf * dtToMc_ee_sf;
    double syserr_dt_ee_bin0 = dt_ee_bin0 * sqrt( pow( relerr_wz_ee_sf, 2 ) + pow( relerr_dtToMc_ee_sf, 2 ) );

    std::cout << " Expected ee data in 0 lept. bin: " << std::endl;
    std::cout << "  " << dt_ee_bin0 << " +- " << err_dt_ee_bin0 << " (stat) +- "
	      << syserr_dt_ee_bin0 << " (syst)" << std::endl;
    */


    std::cout << std::endl;
    std::cout << "MC (LOOSE!):" << std::endl;
    // --- MM ---
    std::cout << " - Channel mm:" << std::endl;
    double wz_ls_mm_bin1 = allhistos["h_wz_mm_extraLeptonNumberAlpha"]->GetBinContent(2);
    double err_wz_ls_mm_bin1 = allhistos["h_wz_mm_extraLeptonNumberAlpha"]->GetBinError(2);
    double dt_ls_mm_bin1 = alldatahistos["h_data_mm_extraLeptonNumberAlpha"]->GetBinContent(2);
    double err_dt_ls_mm_bin1 = alldatahistos["h_data_mm_extraLeptonNumberAlpha"]->GetBinError(2);

    std::cout << " 1 lept. (loose): " << wz_ls_mm_bin1 << " +- "  << err_wz_ls_mm_bin1 << std::endl;
    std::cout << std::endl;
    double wz_ls_mm_sf = wz_mm_bin0/wz_ls_mm_bin1;
    double relerr_wz_ls_mm_sf = sqrt( pow(err_wz_mm_bin0/wz_mm_bin0, 2) + pow(err_wz_ls_mm_bin1/wz_ls_mm_bin1, 2) );
    double err_wz_ls_mm_sf = wz_ls_mm_sf * relerr_wz_ls_mm_sf; 

    std::cout << " scale factor = (0 lept.)/(1 lept. loose) = " << wz_ls_mm_sf << " +- "  
	      << err_wz_ls_mm_sf << " (" << relerr_wz_ls_mm_sf*100 << "%)" << std::endl;
    std::cout << std::endl;

    // --- EE ---
    std::cout << std::endl;
    std::cout << " - Channel ee:" << std::endl;
    double wz_ls_ee_bin1 = allhistos["h_wz_ee_extraLeptonNumberAlpha"]->GetBinContent(2);
    double err_wz_ls_ee_bin1 = allhistos["h_wz_ee_extraLeptonNumberAlpha"]->GetBinError(2);
    double dt_ls_ee_bin1 = alldatahistos["h_data_ee_extraLeptonNumberAlpha"]->GetBinContent(2);
    double err_dt_ls_ee_bin1 = alldatahistos["h_data_ee_extraLeptonNumberAlpha"]->GetBinError(2);

    std::cout << " 1 lept. (loose): " << wz_ls_ee_bin1 << " +- "  << err_wz_ls_ee_bin1 << std::endl;
    std::cout << std::endl;
    double wz_ls_ee_sf = wz_ee_bin0/wz_ls_ee_bin1;
    double relerr_wz_ls_ee_sf = sqrt( pow(err_wz_ee_bin0/wz_ee_bin0, 2) + pow(err_wz_ls_ee_bin1/wz_ls_ee_bin1, 2) );
    double err_wz_ls_ee_sf = wz_ee_sf * relerr_wz_ee_sf; 

    std::cout << " scale factor = (0 lept.)/(1 lept. loose) = " << wz_ls_ee_sf << " +- "  
	      << err_wz_ls_ee_sf << " (" << relerr_wz_ls_ee_sf*100 << "%)" << std::endl;
    std::cout << std::endl;


    // ------------------------------------

    std::cout << std::endl;
    std::cout << "Data (full selection, loose SF):" << std::endl;
    // --- MM ---
    std::cout << " - Channel mm:" << std::endl;

    double dt_ls_mm_bin0 = dt_ls_mm_bin1 * wz_ls_mm_sf;
    double err_dt_ls_mm_bin0 = err_dt_ls_mm_bin1 * wz_ls_mm_sf;
    double syserr_dt_ls_mm_bin0 = dt_ls_mm_bin0 * relerr_wz_ls_mm_sf;

    std::cout << " * mm data in 1 lept. bin (with loose selection): " 
	      << dt_ls_mm_bin0 << " +- " << err_dt_ls_mm_bin1 << std::endl;
    std::cout << " * Expected mm data in 0 lept. bin (with loose SF): " << std::endl;
    std::cout << "    " << dt_ls_mm_bin0 << " +- " << err_dt_ls_mm_bin0 << " (stat) +- "
	      << syserr_dt_ls_mm_bin0 << " (syst)" << std::endl;

    // --- EE ---
    std::cout << " - Channel ee:" << std::endl;

    double dt_ls_ee_bin0 = dt_ls_ee_bin1 * wz_ls_ee_sf;
    double err_dt_ls_ee_bin0 = err_dt_ls_ee_bin1 * wz_ls_ee_sf;
    double syserr_dt_ls_ee_bin0 = dt_ls_ee_bin0 * relerr_wz_ls_ee_sf;

    std::cout << " * ee data in 1 lept. bin (with loose selection): " 
	      << dt_ls_ee_bin0 << " +- " << err_dt_ls_ee_bin1 << std::endl;
    std::cout << " * Expected ee data in 0 lept. bin (with loose SF): " << std::endl;
    std::cout << "    " << dt_ls_ee_bin0 << " +- " << err_dt_ls_ee_bin0 << " (stat) +- "
	      << syserr_dt_ls_ee_bin0 << " (syst)" << std::endl;


  } // end if(doPlot)

  // -- // -- // -- // -- // -- // -- // -- // -- // -- // -- // -- // 
  // -- // -- // -- // -- // -- // -- // -- // -- // -- // -- // -- // 
  // -- // -- // -- // -- // -- // -- // -- // -- // -- // -- // -- // 

  return;
}


void fillPlots(std::map<TString, TH1F*> & histos, TString plotlab, double wght, double puWght) {

  bool isMC( nmcparticles>0 ); 
  if(puWght<0.) puWght = wght;

  // Energy corrections for electrons
  if( isMC==false && cat==2) {
    double enCorr1 = en_corren[l1_pid]/l1_en; if(enCorr1==0) enCorr1 = 1.0;
    double enCorr2 = en_corren[l2_pid]/l2_en; if(enCorr2==0) enCorr2 = 1.0;

    l1_px *= enCorr1; l1_py *= enCorr1; l1_pz *= enCorr1; l1_en *= enCorr1; 
    l2_px *= enCorr2; l2_py *= enCorr2; l2_pz *= enCorr2; l2_en *= enCorr2; 
  }

  double pt1 = getPt(l1_px, l1_py);
  double eta1 = getEta(l1_px, l1_py, l1_pz);
  double aeta1 = fabs(eta1);
  double pt2 = getPt(l2_px, l2_py);
  double eta2 = getEta(l2_px, l2_py, l2_pz);
  double aeta2 = fabs(eta2);

  // Trigger efficiencies, Reco+ID+Isolation efficiencies
  if( isMC ) {
    if( cat==1 ) {  // mm
      wght *= ( muonTriggerEfficiency(pt1, aeta1) * muonTriggerEfficiency(pt2, aeta2) );
      wght *= ( muonScaleFactor(pt1, aeta1)       * muonScaleFactor(pt2, aeta2) );
    }

    else if( cat==2 ) {  // ee
      wght *= ( electronTriggerEfficiency(pt1, aeta1) * electronTriggerEfficiency(pt2, aeta2) );
      wght *= ( electronScaleFactor(pt1, aeta1)       * electronScaleFactor(pt2, aeta2) );
    }

    else if( cat==3 ) {  // em, me
      if( abs(l1_id)==11 ) {
	wght *= ( electronTriggerEfficiency(pt1, aeta1) * muonTriggerEfficiency(pt2, aeta2) );
	wght *= ( electronScaleFactor(pt1, aeta1)       * muonScaleFactor(pt2, aeta2) );
      }
      else {
	wght *= ( muonTriggerEfficiency(pt1, aeta1) * electronTriggerEfficiency(pt2, aeta2) );
	wght *= ( muonScaleFactor(pt1, aeta1)       * electronScaleFactor(pt2, aeta2) );
      }
    }

    else {}
  }

  Float_t thismass = getMass(l1_en+l2_en, l1_px+l2_px, l1_py+l2_py, l1_pz+l2_pz);
  histos[plotlab+"dileptMass"]->Fill(thismass, wght);
  //bool peakBool(fabs(thismass-91.)<15.);
  //bool peakBool(fabs(thismass-91.)<10.);
  bool peakBool(fabs(thismass-91.)<7.);
  if(peakBool)
    if( allVarsToPrintBins.count("dileptMass")>0 ) histos[plotlab+"cutFlow"]->Fill( float(allVarsToPrintBins["dileptMass"]), wght);

  // Z pt > 30
  Float_t zPt = getPt(l1_px+l2_px, l1_py+l2_py);
  histos[plotlab+"dileptPt"]->Fill(zPt, wght);
  //if(zPt<30.) return;
  bool zPtBool(zPt>30.);
  if(peakBool && zPtBool)
    if( allVarsToPrintBins.count("dileptPt")>0 ) histos[plotlab+"cutFlow"]->Fill( float(allVarsToPrintBins["dileptPt"]), wght);

  /////
  //// WARNING!!!
  //// Use either "jn" or "ajn" jets!
  //// Make sure to change everything consistently!!!
  ////
  // // Using jn (PF charge-subtracted)
  // Int_t use_jn = 0;
  // Float_t use_jn_px[maxjn];
  // Float_t use_jn_py[maxjn];
  // Float_t use_jn_pz[maxjn];
  // Float_t use_jn_btag2[maxjn];
  // Bool_t use_jn_tightId[maxjn];
  // Float_t use_htvec_px = 0.;
  // Float_t use_htvec_py = 0.;

  // for(int jdx=0; jdx<jnum; ++jdx) {
  //   use_htvec_px += jn_px[jdx]; 
  //   use_htvec_py += jn_py[jdx]; 
  //   double j_pt = getPt(jn_px[jdx], jn_py[jdx]); 
  //   double j_eta = getEta(jn_px[jdx], jn_py[jdx], jn_pz[jdx]); 
  //   double jetCorrSF = 1.0; // jetSmearingFactor( jn_genpt[jdx], j_pt, j_eta ); 
  //   if( isMC==false ) jetCorrSF = 1.0; 
  //   if( jetCorrSF*j_pt<15. || jn_tightId[jdx]==0 ) continue;

  //   useV_jn_px[use_jn] = jetCorrSF*jn_px[jdx]; // with SF (N.B. also jn_en should be corrected)
  //   useV_jn_py[use_jn] = jetCorrSF*jn_py[jdx]; // with SF (N.B. also jn_en should be corrected)
  //   useV_jn_pz[use_jn] = jn_pz[jdx];
  //   useV_jn_btag2[use_jn] = jn_btag2[jdx];
  //   useV_jn_tightId[use_jn] = jn_tightId[jdx];
  //   use_jn += 1;
  // }

  // Using ajn (PF NON-charge-subtracted)
  Int_t use_jn = 0;
  Float_t use_jn_px[maxajn];
  Float_t use_jn_py[maxajn];
  Float_t use_jn_pz[maxajn];
  Float_t use_jn_btag2[maxajn];
  Bool_t use_jn_tightId[maxajn];
  Float_t use_htvec_px = 0.;
  Float_t use_htvec_py = 0.;

  for(int jdx=0; jdx<ajnum; ++jdx) {
    use_htvec_px += ajn_px[jdx]; 
    use_htvec_py += ajn_py[jdx]; 
    double j_pt = getPt(ajn_px[jdx], ajn_py[jdx]); 
    double j_eta = getEta(ajn_px[jdx], ajn_py[jdx], ajn_pz[jdx]); 
    double jetCorrSF = 1.0; // jetSmearingFactor( ajn_genpt[jdx], j_pt, j_eta ); 
    if( isMC==false ) jetCorrSF = 1.0; 
    if( jetCorrSF*j_pt<15. ) continue;

    use_jn_px[use_jn] = jetCorrSF*ajn_px[jdx]; // with SF (N.B. also jn_en should be corrected)
    use_jn_py[use_jn] = jetCorrSF*ajn_py[jdx]; // with SF (N.B. also jn_en should be corrected)
    use_jn_pz[use_jn] = ajn_pz[jdx];
    use_jn_btag2[use_jn] = ajn_btag2[jdx];
    use_jn_tightId[use_jn] = ajn_tightId[jdx];
    use_jn += 1;
  }

  // Anti-b-tag only jets>20
  Int_t idxHighestCsv = -1;
  //Float_t thisbtagval = getMaxValue(jnum, jn_btag2, idxHighestCsv, jn_px, jn_py, jn_pz, 20., 2.5);
  //Float_t thisbtagval = getMaxValue(jnum, jn_btag2, idxHighestCsv, jn_px, jn_py, jn_pz, 20., 2.5, jn_tightId);
  Float_t thisbtagval = getMaxValue(use_jn, use_jn_btag2, idxHighestCsv, use_jn_px, use_jn_py, use_jn_pz, 20., 2.5); // tightId not required
  histos[plotlab+"jetCsv"]->Fill(thisbtagval, wght);
  bool csvBool(thisbtagval<0.244);
  if(peakBool && zPtBool && csvBool)
    if( allVarsToPrintBins.count("jetCsv")>0 ) histos[plotlab+"cutFlow"]->Fill( float(allVarsToPrintBins["jetCsv"]), wght);
  // // Alternative method (used in non-res background macros)
  // std::vector<UInt_t> jets20Btagged = getListOfParticlesWithThreshold(use_jn, use_jn_btag2, 0.244, use_jn_px, use_jn_py, use_jn_pz, 20., 2.5); 
  // bool csvBool( jets20Btagged.size()==0 );
  // if(peakBool && zPtBool && csvBool)
  //   if( allVarsToPrintBins.count("jetCsv")>0 ) histos[plotlab+"cutFlow"]->Fill( float(allVarsToPrintBins["jetCsv"]), wght);

  // Jet veto only jets>30
  //std::vector<UInt_t> jets30 = getListOfParticlesWithPt(jnum, jn_px, jn_py, 30.);
  //std::vector<UInt_t> jets30 = getListOfParticlesWithPt(jnum, jn_px, jn_py, 30., 0, 0, 0, jn_tightId);
  std::vector<UInt_t> jets30 = getListOfParticlesWithPt(use_jn, use_jn_px, use_jn_py, 30.); // tightId already required
  unsigned int jets30N = jets30.size();
  histos[plotlab+"jetNumber"]->Fill(jets30N, wght);
  bool jets30NBool(jets30N==0);
  if(peakBool && zPtBool && csvBool && jets30NBool)
    if( allVarsToPrintBins.count("jetNumber")>0 ) histos[plotlab+"cutFlow"]->Fill( float(allVarsToPrintBins["jetNumber"]), wght);

  // Ind. minimized CMS RedMET > 50
  std::vector<double> redMetCompts(12, 0.); 
  //Float_t thisCMSredMet = getCMSRedMet(l1_px, l1_py, 0., l2_px, l2_py, 0., htvec_px, htvec_py, met_pt[0], met_phi[0], cat); // no lept. uncert.
  Float_t thisCMSredMet = getCMSRedMet(l1_px, l1_py, 0., l2_px, l2_py, 0., use_htvec_px, use_htvec_py, met_pt[0], met_phi[0], redMetCompts, cat); // no lept. uncert.
  histos[plotlab+"cmsIndMinRedMet"]->Fill(thisCMSredMet, wght);
  //bool redMetBool(thisCMSredMet>50.);
  //bool redMetBool(thisCMSredMet>55.);
  //bool redMetBool(thisCMSredMet>60.);
  bool redMetBool(thisCMSredMet>70.);
  if(peakBool && zPtBool && csvBool && jets30NBool && redMetBool)
    if( allVarsToPrintBins.count("cmsIndMinRedMet")>0 ) histos[plotlab+"cutFlow"]->Fill( float(allVarsToPrintBins["cmsIndMinRedMet"]), wght);

  // 0.4 < bal < 1.8
  Float_t metPtBal = ( zPt>0 ? met_pt[0]/zPt : -1. );
  histos[plotlab+"metPtBalance"]->Fill(metPtBal, wght);
  bool balBool(metPtBal>0.4 && metPtBal<1.8);
  if(peakBool && zPtBool && csvBool && jets30NBool && redMetBool && balBool)
    if( allVarsToPrintBins.count("metPtBalance")>0 ) histos[plotlab+"cutFlow"]->Fill( float(allVarsToPrintBins["metPtBalance"]), wght);

  Float_t dPhiJetMet = 4.0;
  Int_t idxClosestJet = -1;
  //dPhiJetMet = getParticleClosestInPhi(jnum, jn_px, jn_py, met_phi[0], idxClosestJet, 20.);  // loose cut
  //dPhiJetMet = getParticleClosestInPhi(jnum, jn_px, jn_py, met_phi[0], idxClosestJet, 20., jn_tightId);  // medium cut
  //dPhiJetMet = getParticleClosestInPhi(jnum, jn_px, jn_py, met_phi[0], idxClosestJet, 0., jn_tightId);  // tight cut
  dPhiJetMet = getParticleClosestInPhi(use_jn, use_jn_px, use_jn_py, met_phi[0], idxClosestJet, 20.); // tight cut (NB: tightId not required)
  histos[plotlab+"deltaPhiJetMet"]->Fill(dPhiJetMet, wght);
  bool jetMetPhiBool(dPhiJetMet>0.5);
  if(peakBool && zPtBool && csvBool && jets30NBool && redMetBool && balBool && jetMetPhiBool)
    if( allVarsToPrintBins.count("deltaPhiJetMet")>0 ) histos[plotlab+"cutFlow"]->Fill( float(allVarsToPrintBins["deltaPhiJetMet"]), wght);


  // pt>10 both ele and mu 
  std::vector<UInt_t> thirdlept10 = getListOfParticlesWithPt(ln, ln_px, ln_py, 10.);
  unsigned int thirdlept10N = thirdlept10.size();

  unsigned int thirdLeptId = 0;
  double thirdMetTransMass = 0.;

  if(thirdlept10N==1) {
    double tmpPt = 0.;
    if(ln>1) {
      for(int h=0; h<ln; ++h) {
	double tmpPt2 = getPt(ln_px[h], ln_py[h]);
	if(tmpPt2>tmpPt) {
	  thirdLeptId = h;
	  tmpPt = tmpPt2;
	}
      }
    }
    thirdMetTransMass = getTransMass(ln_px[thirdLeptId], ln_py[thirdLeptId], 91.18, 
				     met_pt[0], met_phi[0], 91.18 ); 
  }

  // Full selection
  if(peakBool && zPtBool && csvBool && jets30NBool && redMetBool && balBool && jetMetPhiBool) {
    histos[plotlab+"extraLeptonNumber"]->Fill(thirdlept10N, wght);
    if(thirdlept10N==1) {
      histos[plotlab+"thirdLeptMetTransMass"]->Fill(thirdMetTransMass, wght);
    }
    else if(thirdlept10N==0) {
      histos[plotlab+"finalSelection"]->Fill(1., wght);
      if( allVarsToPrintBins.count("finalSelection")>0 ) histos[plotlab+"cutFlow"]->Fill( float(allVarsToPrintBins["finalSelection"]), wght);
    }
    else {}
  }

  // Softer selection for data/MC normalization
  bool softRedMetBool(thisCMSredMet>40.);
  if(peakBool && zPtBool && csvBool && jets30NBool && softRedMetBool && balBool && jetMetPhiBool) {
    histos[plotlab+"extraLeptonNumberAlpha"]->Fill(thirdlept10N, wght);
    if(thirdlept10N==1) {
      histos[plotlab+"thirdLeptMetTransMassAlpha"]->Fill(thirdMetTransMass, wght);
    }
  }


  //if( thirdlept10N!=1 ) return; // 3rd lepton anti-veto: accept only if event has EXACTLY one extra lepton
  // if( thirdlept10N==0 ) {
  //   histos[plotlab+"extraLepton0"]->Fill(1., wght);
  // }
  // else if( thirdlept10N==1 ) {
  //   histos[plotlab+"extraLepton1"]->Fill(1., wght);
  // }

  //Float_t extraFact = getWzCorrection( cat, ln_id[thirdlept10[0]] );

  return;
}

void printEvents(std::map<int, TString> & finst, 
		 std::map<TString, TH1F*> & higgshistos, unsigned int nhiggslabs, TString *higgslabs,
		 std::map<TString, TH1F*> & histos,      unsigned int nlabs,      TString *labs,      bool mergetop, 
		 std::map<TString, TH1F*> & tothistos, 
		 std::map<TString, TH1F*> & datahistos) {

  std::cout << std::endl;

  // 
  // Loop over final states
  // 
  for(std::map<int, TString>::iterator fs=finst.begin(); fs!=finst.end(); ++fs) {

    std::cout << fs->second.Data() << std::endl;

    // Headers
    for(unsigned int h=0; h<allVarsToPrint.size(); ++h) {
      std::cout << "\t" << allVarsToPrint[h];
    }

    // 
    // Loop over Higgs samples
    // 
    for(unsigned int j=0; j<nhiggslabs; ++j) {
      std::cout << std::endl;
      std::cout << higgslabs[j];

      TString plotlab = "h_"+higgslabs[j]+"_"+fs->second+"_cutFlow";

      // 
      // Loop over variables
      // 
      for(unsigned int h=0; h<allVarsToPrint.size(); ++h) {
	double tmpVal = higgshistos[plotlab]->GetBinContent(h+1);
	double tmpErr = higgshistos[plotlab]->GetBinError(h+1);
	int errOrd = approxToN(tmpErr, 1);
	approxToN(tmpVal, 1, errOrd);
	std::cout << "\t" << tmpVal;
      }
    }

    // 
    // Loop over samples
    // 
    bool stoptop = false;
    for(unsigned int j=0; j<nlabs; ++j) {

      TString thislab = labs[j];
      bool mergingtop = (thislab.BeginsWith("t_") || thislab.BeginsWith("tbar_")) && mergetop;
      if(mergingtop && stoptop) continue;

      if(mergingtop) {
	thislab = "t";
	stoptop = true;
      }

      std::cout << std::endl;
      std::cout << thislab;

      TString plotlab = "h_"+thislab+"_"+fs->second+"_cutFlow";

      // 
      // Loop over variables
      // 
      for(unsigned int h=0; h<allVarsToPrint.size(); ++h) {
	double tmpVal = histos[plotlab]->GetBinContent(h+1);
	double tmpErr = histos[plotlab]->GetBinError(h+1);
	int errOrd = approxToN(tmpErr, 1);
	approxToN(tmpVal, 1, errOrd);
	std::cout << "\t" << tmpVal;
      }
    }

    // 
    // Total MC
    // 
    {
      std::cout << std::endl;
      std::cout << "Tot. MC";

      TString plotlab = "h_sumAllMc_"+fs->second+"_cutFlow";

      // 
      // Loop over variables
      // 
      for(unsigned int h=0; h<allVarsToPrint.size(); ++h) {
	double tmpVal = tothistos[plotlab]->GetBinContent(h+1);
	double tmpErr = tothistos[plotlab]->GetBinError(h+1);
	int errOrd = approxToN(tmpErr, 1);
	approxToN(tmpVal, 1, errOrd);
	std::cout << "\t" << tmpVal
		  << " +- " << tmpErr;
      }
    }

    // 
    // Data
    // 
    {
      std::cout << std::endl;
      std::cout << "data";

      TString plotlab = "h_data_"+fs->second+"_cutFlow";

      // 
      // Loop over variables
      // 
      for(unsigned int h=0; h<allVarsToPrint.size(); ++h) {
	double tmpVal = datahistos[plotlab]->GetBinContent(h+1);
	double tmpErr = datahistos[plotlab]->GetBinError(h+1);
	int errOrd = approxToN(tmpErr, 1);
	approxToN(tmpVal, 1, errOrd);
	std::cout << "\t" << tmpVal
		  << " +- " << tmpErr;
      }
    }

    std::cout << std::endl;
  }

  return;
}
