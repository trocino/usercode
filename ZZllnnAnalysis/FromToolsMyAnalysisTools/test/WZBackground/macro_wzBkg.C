/*
 *  $Date: 2011/07/11 19:28:25 $
 *  \author D. Trocino   - Northeastern University
 */

/*
    Usage: 
      root -l
      .L macro_wzBkg.C+
      macro_wzBkg()
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
void macro_wzBkg() {

  int loaded = -1;
  if( gSystem->IsFileInIncludePath("usefulFunctions_C.so") )
    loaded = gSystem->Load("usefulFunctions_C.so");
  if(loaded < 0)
    gSystem->CompileMacro("usefulFunctions.C");

  initializeGlobalVariables();

  // Samples location
  TString indir="/home/daniele/Documents/Work/samples/2012-01-09_analysis_eleCorrections_detIso_rhoCorr";
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

  //
  // Preselection cuts 
  // 
  addCut("Cosmic rejection", "abs(CosmicMuon)<0.1");
  addCut("Lepton flavor", "abs(Flavor-1)<0.1");  // 1: mm, 2: ee, 3: em
  addCut("Lepton charge", "leadCharge*subleadCharge<0");
  addCut("Invariant mass", "dileptInvMass>76.1876 && dileptInvMass<106.1876");
  addCut("Jet veto", "Jet_number<0.5");
  addCut("Third lepton anti-veto", "numberOfLeptons>0.5"); // anti-third lepton veto
  addCut("D0 Red-MET", "getD0RedMet(Flavor, dilepPROJLong, uncertPROJLong, unclPROJLong, sumOppositeJetPROJLong, dilepPROJPerp, uncertPROJPerp, unclPROJPerp, sumOppositeJetPROJPerp)>50");

  // 
  // Weights 
  // 

  makeWeightDistribution("mc_observed_summer11_s4_35", "data_observed_prompt_all2011_35", puWeights); // PU scenario for MC and data

  addWeight("4653*crossSection*branchingRatio/GenEvents", alllabelstouse, nSmpsToUse); // standard normalization
  addWeight("getPuWeights(nObsPUinter0)", alllabelstouse, nSmpsToUse);                 // PU reweighting
  addWeight("getOverallNorm(Flavor)", alllabelstouse, nSmpsToUse);                     // overall normalization (from Z peak)

  // For WZ measurement
  addWeight("getWzInefficiency(Flavor, thirdLept_flavor, run, thirdLept_eta)", alllabelstouse, nSmpsToUse);

  // Variables to be plotted + binning and range limits (optional)
  //   arguments: 
  //    - name of branch in tree
  //    - x axis title [unit]
  //    - y axis title; if empty: "Events/<bin+unit>" 
  //    - string containing: n. bins, first bin, last bin
  addVariable("dileptInvMass", "Dilepton invariant mass [GeV/c^{2}]", "", "6,70,110");

  //
  // ::::::::::::::::::::::::::::::::::::::::::::::::::::
  // ::::::::::::::::::::::::::::::::::::::::::::::::::::
  // ::::::::::::::::::::::::::::::::::::::::::::::::::::
  //

  // 
  // Options for plots and tables
  // 

  // Want to draw plots? (Here always false)
  bool doPlot=false;

  // Want to print a cut-flow table? (Here always true)
  bool doCutTable=true;

  // Use data, or only MC?
  bool drawData=true;

  // Invert order of processes (in table and plots)?
  bool invertOrder=false;

  //
  // Drawing options (if doPlot=true)
  //
  //  0: 1D, no stack;  1: 1D, stack;  2: 2D, scatter
  int drawOpts=1;
  if(drawOpts==2) {
    allcols["z"]=kOrange+1;
  }
  //
  TString dostack="";
  TString gdopt="";
  TString legopt="lpf";

  // Bin-by-bin comparison?
  bool binbybinComp=false;

  // Draw only MC cumulative plot?
  bool drawOnlyCumulative = false;

  double canvx(500.), canvy(500.);
  if(binbybinComp && drawData && (!doCutTable)) {
    //canvx=500.;
    canvy=750.;
  }

  // -- // -- // -- // -- // -- // -- // -- // -- // -- // -- // -- // 
  // -- // -- // -- // -- // -- // -- // -- // -- // -- // -- // -- // 
  // -- // -- // -- // -- // -- // -- // -- // -- // -- // -- // -- // 

  if(doPlot) {
    std::vector<THStack*> allstacks;
    std::vector<TCanvas*> allcanvas;
    std::vector<TPad*>    allpads1;
    std::vector<TPad*>    allpads2;
    std::vector<TPad*>    allpads3;
    std::vector<TLegend*> alllegends;
    std::vector<std::vector<TH1F*> > allhistos;
    std::vector<TH1F*> allmchistos;
    std::vector<TH1F*> allmcbkghistos;
    std::vector<TH1F*> alldatahistos;
    std::vector<TH1F*> allpoisdiffhistos;
    std::vector<TH1F*> allfracdiffhistos;

    for(unsigned int i=0; i<nVars; ++i) {
      if(binnings[i].Length()==0) {
	cout << " *** New rule: you MUST select a binning!!! ***" << endl;
	throw std::exception();
	return;
      }
      TString thisbins=binnings[i]( binnings[i].Index("(")+1, binnings[i].Index(")")-binnings[i].Index("(")-1);
      TObjArray *bintok=thisbins.Tokenize(","); bintok->SetOwner(kTRUE);
      TString binsStr=((TObjString*)(*bintok)[0])->GetString(); 
      TString lowBinStr=((TObjString*)(*bintok)[1])->GetString(); 
      TString highBinStr=((TObjString*)(*bintok)[2])->GetString(); 
      int bins=binsStr.Atoi();
      double lowBin=lowBinStr.Atof();
      double highBin=highBinStr.Atof();
      
      TString varlab=variables[i];
      varlab.ReplaceAll(" ", ""); 
      varlab.ReplaceAll("(", "_"); varlab.ReplaceAll(")", "_"); 
      varlab.ReplaceAll("<", "LT"); varlab.ReplaceAll(">", "GT"); 
      varlab.ReplaceAll("+", "_P_"); varlab.ReplaceAll("-", "_M_"); 
      varlab.ReplaceAll("*", "_T_"); varlab.ReplaceAll("/", "_O_"); 
      allstacks.push_back(new THStack( ("stack_"+varlab).Data(), ("stack_"+varlab).Data() ));
      allcanvas.push_back(new TCanvas( ("canvas_"+varlab).Data(), ("canvas_"+varlab).Data(), canvx, canvy ));
      if(binbybinComp && drawData) {
	allcanvas.back()->cd();
	allpads1.push_back(new TPad( ("pad1_"+varlab).Data(), ("pad1_"+varlab).Data(), 0., 1.-canvx/canvy, 1., 1. ));
	allpads2.push_back(new TPad( ("pad2_"+varlab).Data(), ("pad2_"+varlab).Data(), 0., (1.-canvx/canvy)/2., 1., 1.-canvx/canvy ));
	allpads3.push_back(new TPad( ("pad3_"+varlab).Data(), ("pad3_"+varlab).Data(), 0., 0., 1., (1.-canvx/canvy)/2. ));
	allpads1.back()->Draw();
	allpads2.back()->Draw();
	allpads3.back()->Draw();
      }
      else {
	allcanvas.back()->cd();
	allpads1.push_back(new TPad( ("pad1_"+varlab).Data(), ("pad1_"+varlab).Data(), 0., 0., 1., 1. ));
	allpads1.back()->Draw();
      }
      alllegends.push_back(new TLegend(0.76, 0.6, 0.98, 0.91));
      alllegends.back()->SetFillColor(0);
      alllegends.back()->SetFillStyle(0);
      alllegends.back()->SetBorderSize(0);
      alllegends.back()->SetLineColor(0);
      allhistos.push_back(std::vector<TH1F*>());
      THStack *stk=allstacks.back();

      for(unsigned int j0=0; j0<nSmpsToUse; ++j0) {
	unsigned int j=(invertOrder ? nSmpsToUse-1-j0 : j0);
	TString dopt="";
	if(j0!=0) dopt="same";
	allpads1[i]->cd();
	allhistos[i].push_back( new TH1F( ("h_"+alllabelstouse[j]+"_"+varlab).Data(), ("h_"+alllabelstouse[j]+"_"+varlab).Data(), bins, lowBin, highBin) );
	allhistos[i].back()->Sumw2();
	alltrees[alllabelstouse[j]]->Draw( (variables[i]+">>h_"+alllabelstouse[j]+"_"+varlab).Data(),
					   ("("+allweightspersample[alllabelstouse[j].Data()]+")*("+allcuts+")").Data(),
					   dopt.Data() );
	// alltrees[alllabelstouse[j]]->Draw( (variables[i]+">>h_"+alllabelstouse[j]+"_"+varlab).Data(),
	// 				   ("("+allweights+")*("+allcuts+")").Data(),
	// 				   dopt.Data() );
	if(allmchistos.size()==i) {
	  allmchistos.push_back( (TH1F*)allhistos[i][0]->Clone( ("h_sumAllMc_"+varlab).Data() ) );
	  //allmchistos.back()->Sumw2();
	}
	else {
	  allmchistos[i]->Add(allhistos[i][j0]);
	}
	if(allmcbkghistos.size()==i) {
	  if(!alllabelstouse[j].EndsWith("zz")) {
	    allmcbkghistos.push_back( (TH1F*)allhistos[i][0]->Clone( ("h_sumAllBkgMc_"+varlab).Data()) );
	    //allmcbkghistos.back()->Sumw2();
	  }
	}
	else {
	  allmcbkghistos[i]->Add(allhistos[i][j0]);
	}

	// Drawing options
	if(drawOpts==0) {        // 1D, no stack
	  allhistos[i][j0]->SetLineWidth(2);
	  allhistos[i][j0]->SetLineColor(allcols[alllabelstouse[j]]);
	  dostack="nostack";
	  gdopt="HIST";
	}
	else if(drawOpts==1) {   // 1D, stack
	  allhistos[i][j0]->SetLineWidth(2);
	  allhistos[i][j0]->SetFillColor(allcols[alllabelstouse[j]]);
	  gdopt="HIST";
	}
	else if(drawOpts==2) {   // 2D, scatter
	  allhistos[i][j0]->SetMarkerColor(allcols[alllabelstouse[j]]);
	  allhistos[i][j0]->SetFillColor(allcols[alllabelstouse[j]]);  // only for legend
	  gdopt="SCAT";
	  drawData=false;
	  legopt="f";
	}
	alllegends[i]->AddEntry(allhistos[i][j0], alllegendlabs[alllabelstouse[j]].Data(), legopt.Data());
	stk->Add(allhistos[i][j0]);
      }
      if(drawData) {
	allpads1[i]->cd();
	alldatahistos.push_back( new TH1F( ("h_data_"+varlab).Data(), ("h_data_"+varlab).Data(), bins, lowBin, highBin) );
	// THIS IS THE STANDARD, UNCOMMENT IT!!!
	t_data->Draw( (variables[i]+">>h_data_"+varlab).Data(),
		      allcuts.Data(),
		      "pesame" );
	// THIS IS ONLY FOR THIRD-LEPTON-REWEIGHTING FOR WZ MEASUREMENT!!!
	//std::cout << "WATCH OUT!!! You're using the third-lepton-reweighting for WZ measurement!!!" << std::endl;
	//t_data->Draw( (variables[i]+">>h_data_"+varlab).Data(),
	//	      ("(getWzInefficiency(Flavor, thirdLept_flavor, run, thirdLept_eta))*("+allcuts+")").Data(),
	//	      "pesame" );
	// alldatahistos.push_back( (TH1F*)gDirectory->Get( ("h_data_"+varlab).Data() ) );
	alldatahistos[i]->SetMarkerStyle(20);
	alllegends[i]->AddEntry(alldatahistos[i], "data", legopt.Data());
      }

      /// Final Draw
      allpads1[i]->cd();
      if(drawOpts!=2) 
      	allpads1[i]->SetLogy(1);
      if(drawOnlyCumulative==false) 
	stk->Draw((dostack+gdopt).Data());
      else 
	allmchistos[i]->Draw(gdopt.Data());
      if(drawData) alldatahistos[i]->Draw((gdopt+"pesame").Data());
      if(drawOnlyCumulative==false) 
	alllegends[i]->Draw();
      // Axis, labels & Co.
      //if(ytitles[i].EqualTo("Events/")) {  // doesn't exist in my version of ROOT... O_O
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
	//stk->GetXaxis()->SetTitleOffset(0.9);
	stk->GetYaxis()->SetTitle(ytitles[i].Data());
	stk->GetYaxis()->SetTitleSize(0.05);
	stk->GetYaxis()->SetLabelSize(0.04);
	stk->GetYaxis()->SetTitleOffset(1.3);
      }
      else {
	allmchistos[i]->GetXaxis()->SetTitle(xtitles[i].Data());
	allmchistos[i]->GetXaxis()->SetTitleSize(0.05);
	allmchistos[i]->GetXaxis()->SetLabelSize(0.04);
	//allmchistos[i]->GetXaxis()->SetTitleOffset(0.9);
	allmchistos[i]->GetYaxis()->SetTitle(ytitles[i].Data());
	allmchistos[i]->GetYaxis()->SetTitleSize(0.05);
	allmchistos[i]->GetYaxis()->SetLabelSize(0.04);
	allmchistos[i]->GetYaxis()->SetTitleOffset(1.3);
      }

      if(binbybinComp && drawData) {	// Start bin-by-bin data-MC comparison

	int nbins=alldatahistos[i]->GetNbinsX();
	double normdt = alldatahistos[i]->Integral(0, -1);
	double normmc = allmchistos[i]->Integral(0, -1);
	if(false) cout << normdt << " - " << normmc << endl;  // just to avoid warnings at compilation time ("unused variable")
	allpoisdiffhistos.push_back( (TH1F*)allmchistos[i]->Clone( ("h_poisDiff_"+varlab).Data() ) );
	allpoisdiffhistos.back()->Reset();
	allpoisdiffhistos.back()->GetXaxis()->SetTitle("");
	allpoisdiffhistos.back()->GetYaxis()->SetTitle("No. Poisson #sigma");
	allpoisdiffhistos.back()->GetYaxis()->SetTitleSize(0.12);
	allpoisdiffhistos.back()->GetYaxis()->SetTitleOffset(0.5);
	allpoisdiffhistos.back()->GetXaxis()->SetLabelSize(0.10);
	allpoisdiffhistos.back()->GetYaxis()->SetLabelSize(0.10);
	allpoisdiffhistos.back()->SetMarkerStyle(20);
	allpoisdiffhistos.back()->SetMarkerSize(0.8);
	allpoisdiffhistos.back()->SetMarkerColor(kBlue);
	allfracdiffhistos.push_back( (TH1F*)allmchistos[i]->Clone( ("h_fracDiff_"+varlab).Data() ) );
	allfracdiffhistos.back()->Reset();
	allfracdiffhistos.back()->GetXaxis()->SetTitle("");
	allfracdiffhistos.back()->GetYaxis()->SetTitle("Rel. difference");
	allfracdiffhistos.back()->GetYaxis()->SetTitleSize(0.12);
	allfracdiffhistos.back()->GetYaxis()->SetTitleOffset(0.5);
	allfracdiffhistos.back()->GetXaxis()->SetLabelSize(0.10);
	allfracdiffhistos.back()->GetYaxis()->SetLabelSize(0.10);
	allfracdiffhistos.back()->GetYaxis()->SetRangeUser(-1., 1.);
	allfracdiffhistos.back()->SetMarkerStyle(20);
	allfracdiffhistos.back()->SetMarkerSize(0.8);
	allfracdiffhistos.back()->SetMarkerColor(kRed);
	allfracdiffhistos.back()->SetLineColor(kRed);

	for(int ibin=1; ibin<=nbins; ++ibin) {

	  double ndt = alldatahistos[i]->GetBinContent(ibin);
	  double nmc = allmchistos[i]->GetBinContent(ibin);
	  double err_nmc = allmchistos[i]->GetBinError(ibin); 
	  //double err_tot = sqrt(err_nmc*err_nmc  + ndt); // correct: taking sqrt[sum(weights^2)] as MC error
	  double err_tot = sqrt(nmc  + ndt); // not-so-correct: taking sqrt(N) as MC error
	  double chi2bin=0.;
	  double fracbin=0.;
	  double err_fracbin=0.;

	  //if(nmc>1e-6 || ndt>1e-6) {
	  if(err_tot>1e-6) {

	    chi2bin = (ndt-nmc)/err_tot;

	    if(nmc>1e-6) {
	      fracbin = (ndt-nmc)/nmc;
	      err_fracbin = fracbin*sqrt( 1./(chi2bin*chi2bin) + (err_nmc/nmc)*(err_nmc/nmc) );

	      std::cout << alldatahistos[i]->GetBinCenter(ibin) << ", " << (ndt/normdt)/(nmc/normmc) << std::endl;
	      //std::cout << alldatahistos[i]->GetBinCenter(ibin) << ", " << (ndt)/(nmc) << std::endl;

	    }

	  }

	  allpoisdiffhistos[i]->SetBinContent(ibin, chi2bin);
	  allfracdiffhistos[i]->SetBinContent(ibin, fracbin);
	  allfracdiffhistos[i]->SetBinError(ibin, err_fracbin);

	} // end for(int ibin=0; ibin<=nbins+1; ++ibin)

	TLine *ll=new TLine(lowBin, 0., highBin, 0.);
	ll->SetLineWidth(2);
	ll->SetLineStyle(7);

	allpads2[i]->cd();
	allpoisdiffhistos[i]->Draw("axis");
	ll->Draw();
	allpoisdiffhistos[i]->Draw("psame");
	allpads2[i]->SetGridx(1);
	allpads2[i]->SetGridy(1);

	allpads3[i]->cd();
	allfracdiffhistos[i]->Draw("axis");
	ll->Draw();
	allfracdiffhistos[i]->Draw("pesame");
	allpads3[i]->SetGridx(1);
	allpads3[i]->SetGridy(1);

      } // end if(drawData)
    } // end for(unsigned int i=0; i<nVars; ++i)
  } // end if(doPlot)


  // -- // -- // -- // -- // -- // -- // -- // -- // -- // -- // -- // 
  // -- // -- // -- // -- // -- // -- // -- // -- // -- // -- // -- // 
  // -- // -- // -- // -- // -- // -- // -- // -- // -- // -- // -- // 

  if(doCutTable) {
    const unsigned int nCuts=cutCascade.size();
    double totMc[nCuts];
    TCanvas *ccount=new TCanvas("ccount", "ccount", 100, 100);
    TH1D *hcount=new TH1D("hcount", "hcount", 1, 0., 2.);
    hcount->Sumw2();
    TH1D *hcountTot[nCuts]; //=new TH1D("hcountTot", "hcountTot", 1, 0., 2.);
    TH1D *hWzEffFact=new TH1D("hWzEffFact", "hWzEffFact", 1, 0., 2.);
    TH1D *hWzEffFactData=new TH1D("hWzEffFactData", "hWzEffFactData", 1, 0., 2.);
    const int totSamples((int)nSmpsToUse);
    TH1F *hEventsMc[totSamples];
    THStack *sAllEvents=new THStack("AllEvents", "Number of events");
    TH1F *hEventsData;
    TLegend *nevtsleg=new TLegend(0.76, 0.6, 0.98, 0.91);
    nevtsleg->SetFillColor(0);
    nevtsleg->SetFillStyle(0);
    nevtsleg->SetBorderSize(0);
    nevtsleg->SetLineColor(0);
    TCanvas *cEvents=new TCanvas("cEvents", "canvas_nEvents", canvx, canvy);

    //
    // MC sample rows
    for(int j0=-1; j0<totSamples; ++j0) {
      unsigned int j=(invertOrder ? nSmpsToUse-1-j0 : j0);
      if(j0==-1) {  // first row
      }
      else {
	TString head=alllabelstouse[j];
	std::cout << head.Data();
	hEventsMc[j]=new TH1F(("EventsMc_"+alllabelstouse[j]).Data(), ("Number of "+alllabelstouse[j]+" events").Data(), nCuts, 0, nCuts);
	hEventsMc[j]->SetLineColor(allcols[alllabelstouse[j]]);
	hEventsMc[j]->SetMarkerColor(allcols[alllabelstouse[j]]);
	hEventsMc[j]->SetMarkerStyle(allmrks[alllabelstouse[j]]);
      }
      for(unsigned int k=0; k<nCuts; ++k) {
	TString hcountTotName = "hcountTot_";
	hcountTotName += k;
	if(j0==-1) {  // first row
	  std::cout << "\t" << cutSet[k].Data();
	  totMc[k]=0.;
	  hcountTot[k] = new TH1D(hcountTotName.Data(), hcountTotName.Data(), 1, 0., 2.);
	  hcountTot[k]->Sumw2();
	}
	else { 
	  alltrees[alllabelstouse[j]]->Draw( "1>>hcount", ("("+allweightspersample[alllabelstouse[j]]+")*("+cutCascade[k]+")").Data());
	  if(j0==0) alltrees[alllabelstouse[j]]->Draw( ("1>>"+hcountTotName).Data(), ("("+allweightspersample[alllabelstouse[j]]+")*("+cutCascade[k]+")").Data());
	  else      alltrees[alllabelstouse[j]]->Draw( ("1>>+"+hcountTotName).Data(), ("("+allweightspersample[alllabelstouse[j]]+")*("+cutCascade[k]+")").Data());
	  ///////// For WZ efficiency factor error
	  // if(k==nCuts-1) {
	  //   if(j0==0) alltrees[alllabelstouse[j]]->Draw( "1>>hWzEffFact", ("("+allweightspersample[alllabelstouse[j]]+")*("+cutCascade[k]+")*(getWzInefficiencyError(Flavor, thirdLept_flavor, run, thirdLept_eta))").Data());
	  //   else      alltrees[alllabelstouse[j]]->Draw( "1>>+hWzEffFact", ("("+allweightspersample[alllabelstouse[j]]+")*("+cutCascade[k]+")*(getWzInefficiencyError(Flavor, thirdLept_flavor, run, thirdLept_eta))").Data());
	  // }
	  //////////

	  std::cout << "\t" << hcount->GetBinContent(1) << " +- " << hcount->GetBinError(1);
	  totMc[k]+=hcount->GetBinContent(1);
	  hEventsMc[j]->SetBinContent(k+1, hcount->GetBinContent(1));
	  hEventsMc[j]->GetXaxis()->SetBinLabel(k+1, cutSet[k].Data());
	}
      }
      if(j0!=-1) sAllEvents->Add(hEventsMc[j]);
      if(j0!=-1) nevtsleg->AddEntry(hEventsMc[j], alllegendlabs[alllabelstouse[j]].Data(), "pe");
      std::cout << std::endl;
    }

    std::cout << "Tot. MC";
    for(unsigned int k=0; k<nCuts; ++k) {
      std::cout << "\t" << hcountTot[k]->GetBinContent(1) << " +- " << hcountTot[k]->GetBinError(1);
    }
    std::cout << std::endl;
    //
    // Data row
    if(drawData) {
      hEventsData=new TH1F("EventsMc_data", "Number of events", nCuts, 0, nCuts);
      hEventsData->SetMarkerStyle(20);
      std::cout << "Data";
      for(unsigned int k=0; k<nCuts; ++k) {
	// THIS IS THE STANDARD, UNCOMMENT IT!!!
	//t_data->Draw( "1>>hcount", cutCascade[k].Data() );
	// THIS IS ONLY FOR THIRD-LEPTON-REWEIGHTING FOR WZ MEASUREMENT!!!
	t_data->Draw( "1>>hcount", ("(getWzInefficiency(Flavor, thirdLept_flavor, run, thirdLept_eta))*("+cutCascade[k]+")").Data() );
	///////// For WZ efficiency factor error
	// if(k==nCuts-1) {
	//   t_data->Draw( "1>>hWzEffFactData", ("(getWzInefficiency(Flavor, thirdLept_flavor, run, thirdLept_eta))*("+cutCascade[k]+")*(getWzInefficiencyError(Flavor, thirdLept_flavor, run, thirdLept_eta))").Data() );
	// }
	///////// 
	std::cout << "\t" << hcount->GetBinContent(1) << " +- " << hcount->GetBinError(1);
	hEventsData->SetBinContent(k+1, hcount->GetBinContent(1));
      }
      nevtsleg->AddEntry(hEventsData, "data", "lp");
      std::cout << std::endl;
    }

    delete ccount;

    std::cout << std::endl;

    std::cout << "  WZ eff. factor cumulative error MC:   " << hWzEffFact->GetBinContent(1) << std::endl;
    std::cout << "  WZ eff. factor cumulative error Data: " << hWzEffFactData->GetBinContent(1) << std::endl;

    cEvents->cd();
    cEvents->SetLogy();
    sAllEvents->Draw("penostack");
    if(drawData) hEventsData->Draw("pesame");
    double yMin=( invertOrder ? hEventsMc[totSamples-1]->GetBinContent(nCuts) : hEventsMc[0]->GetBinContent(nCuts) );
    yMin/=10.;
    if(yMin>10e-6) sAllEvents->SetMinimum(yMin);
    nevtsleg->Draw();

  } // end if(doCutTable) 

  return;

}

