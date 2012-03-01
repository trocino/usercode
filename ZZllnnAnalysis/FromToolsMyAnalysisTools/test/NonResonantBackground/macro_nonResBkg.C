/*
 *  $Date: 2012/03/01 09:35:03 $
 *  \author D. Trocino   - Northeastern University
 */

/*
    Usage: 
      root -l
      .L macro_nonResBkg.C+
      macro_nonResBkg()

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

#include "usefulFunctions.C"

using namespace std;

//
// Main function
//
void macro_nonResBkg() {

  int loaded = -1;
  if( gSystem->IsFileInIncludePath("../usefulFunctions_C.so") )
    loaded = gSystem->Load("../usefulFunctions_C.so");
  if(loaded < 0)
    gSystem->CompileMacro("../usefulFunctions.C");

  initializeGlobalVariables();

  // Samples location
  TString indir="/home/daniele/Documents/Work/samples/analysis_2012-02-13";
  TString basename="ZZllvvAnalyzer_any_anyCharge_";
  TString base=indir+"/"+basename;

  // Labels for all processes
  TString alllabels[]={"zz", "wz", "ww", "tt", "t", "zzx", "w", "z"};
  unsigned int nSmps=sizeof(alllabels)/sizeof(TString);

  // Labels for legend
  std::map<TString, TString> alllegendlabs;
  alllegendlabs["zz"]="ZZ #rightarrow 2l2#nu";
  alllegendlabs["wz"]="WZ";
  alllegendlabs["ww"]="WW";
  alllegendlabs["tt"]="tt";
  alllegendlabs["t"]="t";
  alllegendlabs["w"]="W + jets";
  alllegendlabs["zzx"]="ZZ #rightarrow X";
  alllegendlabs["z"]="Z + jets";

  // Colors and marker styles, process by process
  Color_t cols[]={kRed, kMagenta, kViolet-1, kBlue, kCyan, kGreen, kSpring+3, kYellow-7};
  Style_t mrks[]={21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32};
  // Backup colors: kPink, kAzure+1, kTeal+1, kOrange+1
  unsigned int nCols=sizeof(cols)/sizeof(Color_t);
  unsigned int nMrks=sizeof(mrks)/sizeof(Color_t);
  if(nCols<nSmps) {
    std::cout << " ******************** ERROR ********************" << endl;
    std::cout << "    Not enough colors available, add colors!" << std::endl;
    return;
  }
  if(nMrks<nSmps) {
    std::cout << " ******************** ERROR ********************" << endl;
    std::cout << "    Not enough markers available, add styles!" << std::endl;
    return;
  }

  // Labels to be actually used
  TString alllabelstouse[]={"zz", "wz", "ww", "tt", "t", "zzx", "w", "z"};
  unsigned int nSmpsToUse=sizeof(alllabelstouse)/sizeof(TString);
  if(nSmps<nSmpsToUse) {
    std::cout << " ******************** ERROR ********************" << endl;
    std::cout << "    Not enough samples available, add samples!" << std::endl;
    return;
  }

  // Maps associating each process to its tree, color, marker style
  std::map<TString, TChain*> alltrees;
  std::map<TString, Color_t> allcols;
  std::map<TString, Style_t> allmrks;

  for(unsigned int k=0; k<nSmps; ++k) {
    alltrees[alllabels[k]]=new TChain("ntuple_nEvt");
    allcols[alllabels[k]]=cols[k];
    allmrks[alllabels[k]]=mrks[k];
  }

  // Fill maps with trees
  alltrees["zz"]->Add(base+"ZZ_only2l2n.root");
  alltrees["wz"]->Add(base+"WZ.root");
  alltrees["ww"]->Add(base+"WW.root");
  alltrees["tt"]->Add(base+"TTJets.root");
  alltrees["t"]->Add(base+"SingleT_tW.root"); alltrees["t"]->AddFile(base+"SingleTbar_tW.root");
  alltrees["t"]->Add(base+"SingleT_t.root"); alltrees["t"]->AddFile(base+"SingleTbar_t.root");
  alltrees["t"]->Add(base+"SingleT_s.root"); alltrees["t"]->AddFile(base+"SingleTbar_s.root");
  alltrees["w"]->Add(base+"WJetsToLNu_0.root"); alltrees["w"]->AddFile(base+"WJetsToLNu_1.root"); 
  alltrees["w"]->Add(base+"WJetsToLNu_2.root"); alltrees["w"]->AddFile(base+"WJetsToLNu_3.root");
  alltrees["zzx"]->Add(base+"ZZ_allBut2l2n.root");
  alltrees["z"]->Add(base+"DYJetsToLL_0.root"); alltrees["z"]->Add(base+"DYJetsToLL_1.root"); 

  // Data samples
  TChain *t_data=new TChain("ntuple");

  // May 10 rereco
  t_data->Add(base+"DoubleMuMay10ReReco_0.root"); t_data->AddFile(base+"DoubleMuMay10ReReco_1.root");
  t_data->Add(base+"DoubleElectronMay10ReReco_0.root"); t_data->AddFile(base+"DoubleElectronMay10ReReco_1.root"); 
  t_data->Add(base+"MuEGMay10ReReco.root"); 
  // Aug 05 rereco
  t_data->Add(base+"DoubleMuAug05ReReco.root");
  t_data->Add(base+"DoubleElectronAug05ReReco.root");
  t_data->Add(base+"MuEGAug05ReReco.root");
  // Prompt reco V4
  t_data->Add(base+"DoubleMuPromptRecoV4_0.root"); t_data->Add(base+"DoubleMuPromptRecoV4_1.root");
  t_data->Add(base+"DoubleElectronPromptRecoV4_0.root"); t_data->AddFile(base+"DoubleElectronPromptRecoV4_1.root");
  t_data->Add(base+"MuEGPromptRecoV4_0.root"); t_data->Add(base+"MuEGPromptRecoV4_1.root");
  // Prompt reco V6
  t_data->Add(base+"DoubleMuPromptRecoV6.root");
  t_data->Add(base+"DoubleElectronPromptRecoV6.root"); 
  t_data->Add(base+"MuEGPromptRecoV6.root");
  // Prompt reco 2011B
  t_data->Add(base+"DoubleMuPromptReco_2011B_0.root"); t_data->AddFile(base+"DoubleMuPromptReco_2011B_1.root"); t_data->AddFile(base+"DoubleMuPromptReco_2011B_2.root");
  t_data->Add(base+"DoubleElectronPromptReco_2011B_0.root"); t_data->AddFile(base+"DoubleElectronPromptReco_2011B_1.root"); 
  t_data->Add(base+"MuEGPromptReco_2011B_0.root"); t_data->AddFile(base+"MuEGPromptReco_2011B_1.root"); t_data->AddFile(base+"MuEGPromptReco_2011B_2.root");


  //
  // ::::::::::::::::::::::::::::::::::::::::::::::::::::
  // ::::::::::::::::::::::::::::::::::::::::::::::::::::
  // ::::::::::::::::::::::::::::::::::::::::::::::::::::

  //
  // Cuts for non-resonant background measurement
  // 
  TString finalStates[4] = {"1", "3", "2", "3"};  // 1: mm, 2: ee, 3: em
  TString pickAflavor[4] = {"1", "1", "2", "2"};  // Choose which version of reduced-MET to use
  // Peak or side-bands
  TString massDistr[2] = {"dileptInvMass>76.1876 && dileptInvMass<106.1876", "(dileptInvMass>46.1876 && dileptInvMass<76.1876) || (dileptInvMass>106.1876 && dileptInvMass<136.1876)"};

  //
  // ::::::::::::::::::::::::::::::::::::::::::::::::::::
  // ::::::::::::::::::::::::::::::::::::::::::::::::::::
  // ::::::::::::::::::::::::::::::::::::::::::::::::::::
  //

  //
  // Preselection cuts 
  //
  addCut("Cosmic rejection", "abs(CosmicMuon)<0.1");
  addCut("Lepton flavor", "abs(Flavor-LEPTFLAVTEMPL)<0.1");  // 1: mm, 2: ee, 3: em
  addCut("Lepton charge", "leadCharge*subleadCharge<0");
  addCut("Invariant mass", "INVMASSTEMPL"); 
  addCut("Jet veto", "Jet_number<0.5");
  addCut("Third lepton veto", "numberOfLeptons<0.5");
  addCut("D0 Red-MET", "getD0RedMet(Flavor, dilepPROJLong, uncertPROJLong, unclPROJLong, sumOppositeJetPROJLong, dilepPROJPerp, uncertPROJPerp, unclPROJPerp, sumOppositeJetPROJPerp, PICKAFLAVTEMPL)>50");

  // 
  // Weights 
  // 

  makeWeightDistribution("mc_observed_summer11_s4_35", "data_observed_prompt_all2011_35", puWeights); // PU scenario for MC and data

  addWeight("4653*crossSection*branchingRatio/GenEvents", alllabelstouse, nSmpsToUse); // standard normalization
  addWeight("getPuWeights(nObsPUinter0)", alllabelstouse, nSmpsToUse);                 // PU reweighting
  addWeight("getOverallNorm(Flavor)", alllabelstouse, nSmpsToUse);                     // overall normalization (from Z peak)

  // Variables to be plotted + binning and range limits (optional)
  //   arguments: 
  //    - name of branch in tree
  //    - x axis title [unit]
  //    - y axis title; if empty: "Events/<bin+unit>" 
  //    - string containing: n. bins, first bin, last bin
  addVariable("dileptInvMass", "Dilepton invariant mass [GeV/c^{2}]", "", "3,46.1876,136.1876"); 

  //
  // ::::::::::::::::::::::::::::::::::::::::::::::::::::
  // ::::::::::::::::::::::::::::::::::::::::::::::::::::
  // ::::::::::::::::::::::::::::::::::::::::::::::::::::
  //

  // 
  // Arrays for n. events + errors
  // 
  double mc_im_fs[2][4];
  double err_mc_im_fs[2][4];
  double dt_im_fs[2][4];
  double err_dt_im_fs[2][4];

  // Loop over all flavors, loop over peak/side-bands
  for(unsigned int a=0; a<2; ++a) {
    for(unsigned int b=0; b<4; ++b) {
      TString allcuts_new = allcuts;
      allcuts_new.ReplaceAll("LEPTFLAVTEMPL", finalStates[b].Data());
      allcuts_new.ReplaceAll("PICKAFLAVTEMPL", pickAflavor[b].Data());
      allcuts_new.ReplaceAll("INVMASSTEMPL", massDistr[a].Data());
      TH1D *hcount_mc=new TH1D("hcount_mc", "hcount_mc", 1, 0., 2.);
      TH1D *hcount_dt=new TH1D("hcount_dt", "hcount_dt", 1, 0., 2.);
      hcount_mc->Sumw2();
      for(unsigned int k=0; k<nSmpsToUse; ++k) {
	if(k==0)
	  alltrees[alllabelstouse[k]]->Draw( "1>>hcount_mc",
					     ("("+allweightspersample[alllabelstouse[k].Data()]+")*("+allcuts_new+")").Data(),
					     "goff" );
	else
	  alltrees[alllabelstouse[k]]->Draw( "1>>+hcount_mc",
					     ("("+allweightspersample[alllabelstouse[k].Data()]+")*("+allcuts_new+")").Data(),
					     "goff" );    
      }
      mc_im_fs[a][b] = hcount_mc->GetBinContent(1);
      err_mc_im_fs[a][b] = hcount_mc->GetBinError(1);

      t_data->Draw( "1>>hcount_dt", allcuts_new.Data(), "goff" );
      dt_im_fs[a][b] = hcount_dt->GetBinContent(1);
      err_dt_im_fs[a][b] = hcount_dt->GetBinError(1);

      delete hcount_mc;
      delete hcount_dt;
    }
  }

  // Print out results
  cout << "Monte Carlo" << endl;
  cout << "          \tmm                  "
       <<           "\tem                  "
       <<           "\tee                  "
       <<           "\tem" << std::endl;
  cout << "      peak\t" << mc_im_fs[0][0] << " +- " << err_mc_im_fs[0][0] 
       <<           "\t" << mc_im_fs[0][1] << " +- " << err_mc_im_fs[0][1] 
       <<           "\t" << mc_im_fs[0][2] << " +- " << err_mc_im_fs[0][2] 
       <<           "\t" << mc_im_fs[0][3] << " +- " << err_mc_im_fs[0][3] 
       << endl;
  cout << " sidebands\t" << mc_im_fs[1][0] << " +- " << err_mc_im_fs[1][0] 
       <<           "\t" << mc_im_fs[1][1] << " +- " << err_mc_im_fs[1][1] 
       <<           "\t" << mc_im_fs[1][2] << " +- " << err_mc_im_fs[1][2] 
       <<           "\t" << mc_im_fs[1][3] << " +- " << err_mc_im_fs[1][3] 
       << endl << endl;
  cout << "For mm:" << endl;
  double alpha_mc_mm = mc_im_fs[1][0]/mc_im_fs[1][1];
  double relErr_alpha_mc_mm =  sqrt( pow( (err_mc_im_fs[1][0]/mc_im_fs[1][0]), 2) + pow( (err_mc_im_fs[1][1]/mc_im_fs[1][1]), 2) );
  cout << " alpha = " << mc_im_fs[1][0] << " +- " << err_mc_im_fs[1][0] << " / " 
       << mc_im_fs[1][1] << " +- " << err_mc_im_fs[1][1] << " = " 
       << alpha_mc_mm << " +- " << alpha_mc_mm*relErr_alpha_mc_mm << " (" << relErr_alpha_mc_mm*100. << "%)" << endl;
  cout << " residual bkg = " << mc_im_fs[0][1]*alpha_mc_mm << " +- " << err_mc_im_fs[0][1]*alpha_mc_mm << " +- " << relErr_alpha_mc_mm*mc_im_fs[0][1] 
       << endl << endl;
  cout << "For ee:" << endl;
  double alpha_mc_ee = mc_im_fs[1][2]/mc_im_fs[1][3];
  double relErr_alpha_mc_ee =  sqrt( pow( (err_mc_im_fs[1][2]/mc_im_fs[1][2]), 2) + pow( (err_mc_im_fs[1][3]/mc_im_fs[1][3]), 2) );
  cout << " alpha = " << mc_im_fs[1][2] << " +- " << err_mc_im_fs[1][2] << " / " 
       << mc_im_fs[1][3] << " +- " << err_mc_im_fs[1][3] << " = " 
       << alpha_mc_ee << " +- " << alpha_mc_ee*relErr_alpha_mc_ee << " (" << relErr_alpha_mc_ee*100. << "%)" << endl;
  cout << " residual bkg = " << mc_im_fs[0][3]*alpha_mc_ee << " +- " << err_mc_im_fs[0][3]*alpha_mc_ee << " +- " << relErr_alpha_mc_ee*mc_im_fs[0][3] 
       << endl << endl;
  cout << "-------------------------------------------------------" << endl << endl;
  cout << "Data" << endl;
  cout << "          \tmm         "
       <<           "\tem         "
       <<           "\tee         "
       <<           "\tem" << std::endl;
  cout << "      peak\t" << dt_im_fs[0][0] << " +- " << err_dt_im_fs[0][0] 
       <<           "\t" << dt_im_fs[0][1] << " +- " << err_dt_im_fs[0][1] 
       <<           "\t" << dt_im_fs[0][2] << " +- " << err_dt_im_fs[0][2] 
       <<           "\t" << dt_im_fs[0][3] << " +- " << err_dt_im_fs[0][3] 
       << endl;
  cout << " sidebands\t" << dt_im_fs[1][0] << " +- " << err_dt_im_fs[1][0] 
       <<           "\t" << dt_im_fs[1][1] << " +- " << err_dt_im_fs[1][1] 
       <<           "\t" << dt_im_fs[1][2] << " +- " << err_dt_im_fs[1][2] 
       <<           "\t" << dt_im_fs[1][3] << " +- " << err_dt_im_fs[1][3] 
       << endl << endl;
  cout << "For mm:" << endl;
  double alpha_dt_mm = dt_im_fs[1][0]/dt_im_fs[1][1];
  double relErr_alpha_dt_mm =  sqrt( pow( (err_dt_im_fs[1][0]/dt_im_fs[1][0]), 2) + pow( (err_dt_im_fs[1][1]/dt_im_fs[1][1]), 2) );
  cout << " alpha = " << dt_im_fs[1][0] << " +- " << err_dt_im_fs[1][0] << " / " 
       << dt_im_fs[1][1] << " +- " << err_dt_im_fs[1][1] << " = " 
       << alpha_dt_mm << " +- " << alpha_dt_mm*relErr_alpha_dt_mm << " (" << relErr_alpha_dt_mm*100. << "%)" << endl;
  cout << " residual bkg = " << dt_im_fs[0][1]*alpha_dt_mm << " +- " << err_dt_im_fs[0][1]*alpha_dt_mm << " +- " << relErr_alpha_dt_mm*dt_im_fs[0][1] 
       << endl << endl;
  cout << "For ee:" << endl;
  double alpha_dt_ee = dt_im_fs[1][2]/dt_im_fs[1][3];
  double relErr_alpha_dt_ee =  sqrt( pow( (err_dt_im_fs[1][2]/dt_im_fs[1][2]), 2) + pow( (err_dt_im_fs[1][3]/dt_im_fs[1][3]), 2) );
  cout << " alpha = " << dt_im_fs[1][2] << " +- " << err_dt_im_fs[1][2] << " / " 
       << dt_im_fs[1][3] << " +- " << err_dt_im_fs[1][3] << " = " 
       << alpha_dt_ee << " +- " << alpha_dt_ee*relErr_alpha_dt_ee << " (" << relErr_alpha_dt_ee*100. << "%)" << endl;
  cout << " residual bkg = " << dt_im_fs[0][3]*alpha_dt_ee << " +- " << err_dt_im_fs[0][3]*alpha_dt_ee << " +- " << relErr_alpha_dt_ee*dt_im_fs[0][3] 
       << endl << endl;


  return;

}

