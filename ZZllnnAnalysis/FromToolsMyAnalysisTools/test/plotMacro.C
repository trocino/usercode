/*
 *  $Date: 2011/06/01 18:02:26 $
 *  \author D. Trocino   - Northeastern University
 */

#include "TString.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TDirectoryFile.h"
#include "TFile.h"
#include "TH1F.h"
#include <iostream>
#include <map>
#include <sys/stat.h>
#include <sys/types.h>

using namespace std;


TString inputdir("/tmp/trocino/AllData");
TString namefile("ZZllvvAnalyzer_0.root");

TString outname("ZZllvvPlots.root");
TString outdirname("ZZllvvPlots");

// -- For plotMacro -- //
//TString plotnametempl("DileptKin_whatever_hMass");
//TString xlab("Dilepton inv. mass [GeV/c^{2}]");
//TString plotnametempl("RedMetStd_whatever_hRedMET");
//TString xlab("Reduced MET [GeV]");
//TString plotnametempl("SelJetKin_whatever_hNObj");
//TString xlab("# selected jets");
//TString plotnametempl("JetKin_whatever_hNObj");  // ALL jets
//TString xlab("# jets");                          // ALL jets
//TString plotnametempl("LeptonLead_whatever_hNLept");
//TString xlab("# extra leptons");
//TString plotnametempl("SelectedTracks_whatever_hNTrack");
//TString xlab("# extra isolated tracks");
//TString plotnametempl("LeptonLead_whatever_hPt");
//TString xlab("Leading lepton p_{T} [GeV/c]");
//TString plotnametempl("DileptKin_whatever_hCmCosThetaN");
//TString xlab("cos#theta*(l^{-})");
//TString plotnametempl("DileptKin_whatever_hDeltaPhiDL");
//TString xlab("#Delta#phi(dilept, lead) [rad]");
TString plotnametempl("DileptKin_whatever_hPt");
TString xlab("Dilepton p_{T} [GeV/c]");


// -- For plotNumbers -- //
//TString plotnamefix("EventCounter");
//TString xlabN("");
TString plotnamefix("PreEventCounter");
TString xlabN("");
//TString plotnamefix("hNVertexAll");
//TString xlabN("# vertices");
//TString plotnamefix("hNVertexGood");
//TString xlabN("# good vertices");


// -- For all -- //
int rebinValue(1);
bool reversePlotsOrder(false);
float xbinMin(0);
float xbinMax(-1);
int logscale(1);
bool stacked(false);
bool doNormalize(false);
//if(!stacked) doNormalize=true;
bool plotData(false);
//if(!stacked) plotData=false;
bool drawAsHist(false);

// -- Only for 2mu channel -- // 
TString samplesMC[]={"ZZtoAnything_only2l2n_Spring11","WZtoAnything_Spring11","WWtoAnything_Spring11","TTJets_madgraph_Spring11"};
//TString samplesMC[]={"ZZtoAnything_only2l2n_Spring11","WZtoAnything_Spring11","WWtoAnything_Spring11","ZZtoAnything_allBut2l2n_Spring11","TTJets_madgraph_Spring11","TToBLNu_Spring11","WJetsToLNu_Spring11","DYToTauTau_M-20_Spring11","DYToMuMu_M-20_Spring11"};

int sizeMc=sizeof(samplesMC)/sizeof(TString);

Color_t cols[]={kRed, /*kPink,*/ kMagenta, kViolet-1, kBlue, kAzure+1, kCyan, kTeal+1, kGreen+1, /*kSpring,*/ kYellow-7, kOrange+1};


// ---- Initialize sample-axis title map ---- // 
void initiateMap(map<TString, TString> &process) {

  process["ZZtoAnything_only2l2n_Spring11"]="ZZ#rightarrow2l2#nu";
  process["ZZtoAnything_allBut2l2n_Spring11"]="ZZ incl.";
  process["WZtoAnything_Spring11"]="WZ";
  process["WWtoAnything_Spring11"]="WW";
  process["TTJets_madgraph_Spring11"]="t#bar{t}+jets";
  process["TToBLNu_Spring11"]="single t";
  process["WJetsToLNu_Spring11"]="W+jets";
  process["DYToTauTau_M-20_Spring11"]="DY#rightarrow2#tau";
  process["DYToMuMu_M-20_Spring11"]="DY#rightarrow2#mu";
  process["DYToEE_M-20_Spring11"]="DY#rightarrow2e";
  process["QCD_Pt-30to1000_Spring11"]="QCD";

  return;
}


// ---- Function to make a dir ---- // 
void checkAndMakeDir(TString namedir) {
  struct stat mystat;
  if(stat(namedir.Data(), &mystat)!=0) {  // dir <namedir> doesn't exist
    int status=mkdir(namedir.Data(), 0777);
    if(status==-1) {
      cout << " Cannot create directory " << namedir.Data() <<", exiting" << endl;
      return;
    }
  }
  else {
    if(stat(namedir.Data(), &mystat)==0 && ! S_ISDIR(mystat.st_mode)) {  // <namedir> exists, but is not a directory
      namedir+="_2";
      int status=mkdir(namedir.Data(), 0777);
      if(status==-1) {
	cout << " Cannot create directory " << namedir.Data() <<", exiting" << endl;
	return;
      }
    }
  }
}


// --------------------------------------------------------------- //
// --  Macro for plots in subdirectories  ------------------------ //
// --------------------------------------------------------------- //

void plotMacro() {

  /*
  TString samplesMC[]={"WWtoAnything_Spring11","WZtoAnything_Spring11","ZZtoAnything_only2l2n_Spring11","ZZtoAnything_allBut2l2n_Spring11","TToBLNu_tW-channel_Spring11","TToBLNu_t-channel_Spring11","TToBLNu_s-channel_Spring11","TTJets_madgraph_Spring11","WJetsToLNu_Spring11","DYToEE_M-20_Spring11","DYToMuMu_M-20_Spring11","DYToTauTau_M-20_Spring11","QCD_Pt-30to50_Spring11","QCD_Pt-50to80_Spring11","QCD_Pt-80to120_Spring11","QCD_Pt-120to170_Spring11","QCD_Pt-170to300_Spring11","QCD_Pt-300to470_Spring11","QCD_Pt-470to600_Spring11","QCD_Pt-600to800_Spring11","QCD_Pt-800to1000_Spring11","DoubleMuon-v3","DoubleMuon-v5","DoubleMuon-v6","DoubleMuon-v7","DoubleMuon-v8","DoubleMuon-v9","DoubleElectron-v3","DoubleElectron-v5","DoubleElectron-v6","DoubleElectron-v7","DoubleElectron-v8","DoubleElectron-v9","MuEG-v3","MuEG-v5","MuEG-v6","MuEG-v7","MuEG-v8","MuEG-v9"};
  */

  /*
  TString samplesMC[]={"ZZtoAnything_only2l2n_Spring11","WZtoAnything_Spring11","WWtoAnything_Spring11","ZZtoAnything_allBut2l2n_Spring11","TTJets_madgraph_Spring11","TToBLNu_Spring11","WJetsToLNu_Spring11","QCD_Pt-30to1000_Spring11","DYToEE_M-20_Spring11","DYToTauTau_M-20_Spring11","DYToMuMu_M-20_Spring11"};
  */

  map<TString, TString> proc;
  initiateMap(proc);

  TString sampleData="AllData";

  TString subdir[]={"cut1_Dilepton", "cut2_Mass-window", "cut3_RedMET-cut", "cut4_Jet-veto", "cut5_Lepton-veto"};
  int sizeSdir=sizeof(subdir)/sizeof(TString);

  vector<THStack *> stacks;
  vector<TH1F *>    dataplots;
  vector<TCanvas *> ccs;
  vector<TLegend *> lleg;
  vector<TFile *>   allFiles;

  // Old, not consistent with the output of "submission.py --copy"
  //  TFile *fdata=new TFile((inputdir+sampleData+"/"+namefile).Data(), "READ");

  // Assuming all data files (from DoubleMu, DoubleElectron, MuEG streams) 
  // have been merged into a single file named "ZZllvvAnalyzer_AllData.root"
  TString namefiletmp=namefile;
  namefiletmp.ReplaceAll("0", sampleData.Data());
  TFile *fdata=new TFile((inputdir+"/"+namefiletmp).Data(), "READ");

  vector<float> maxV(sizeSdir, 0.);
  vector<float> minV(sizeSdir, 1.e9);

  //  for(vector<TString>::iterator sdr=subdir.begin(); sdr!=subdir.end(); ++sdr) {
  for(int j1=0; j1<sizeSdir; ++j1) {
    TString plotname=plotnametempl;
    plotname.ReplaceAll("whatever", subdir[j1].Data());
    stacks.push_back(new THStack(("stack_"+plotname).Data(), plotname.Data()) );
    ccs.push_back( new TCanvas( ("canvas_"+plotname).Data(), ("canvas_"+subdir[j1]+"_"+plotname).Data(), 600, 600) );
    lleg.push_back(new TLegend(0.75, 0.6, 1., 1.));
    lleg.back()->SetFillColor(0);
    lleg.back()->SetFillStyle(0);
    lleg.back()->SetBorderSize(0);
    lleg.back()->SetLineColor(0);

    if(plotData) {
      TDirectoryFile *dir=(TDirectoryFile*)fdata->Get(subdir[j1].Data());
      TH1F *h=(TH1F*)dir->Get(plotname.Data());
      h->SetName( (TString("data_")+h->GetName()).Data() );
      h->SetMarkerStyle(20);
      h->Rebin(rebinValue);
      h->GetXaxis()->SetRangeUser(xbinMin, xbinMax);
      if(doNormalize) {
	h->Scale(1./h->Integral());
      }
      dataplots.push_back(h);
      lleg.back()->AddEntry(dataplots[0], "data", "pe");
      maxV[j1]=(h->GetBinContent(h->GetMaximumBin())>maxV[j1] ? h->GetBinContent(h->GetMaximumBin()) : maxV[j1]);
      minV[j1]=(h->GetBinContent(h->GetMaximumBin())<minV[j1] && h->GetBinContent(h->GetMaximumBin())>1.e-7 ? h->GetBinContent(h->GetMaximumBin()) : minV[j1]);
    }
  }

  int idx1(0);
  for(int idx0=0; idx0<sizeMc; ++idx0) {
    idx1=(reversePlotsOrder ? sizeMc-idx0-1 : idx0);
    // Old, not consistent with the output of "submission.py --copy"
    //    allFiles.push_back(new TFile((inputdir+samplesMC[idx1]+"/"+namefile).Data(), "READ"));
    TString namefilemctmp=namefile;
    namefilemctmp.ReplaceAll("0", samplesMC[idx1].Data());
    allFiles.push_back(new TFile((inputdir+"/"+namefilemctmp).Data(), "READ"));
    TFile *f=allFiles.back();
    //    for(vector<TString>::iterator sdr=subdir.begin(); sdr!=subdir.end(); ++sdr, ++idx2) {
    for(int idx2=0; idx2<sizeSdir; ++idx2) {
      TDirectoryFile *dir=(TDirectoryFile*)f->Get(subdir[idx2].Data());
      TString plotname=plotnametempl;
      plotname.ReplaceAll("whatever", subdir[idx2].Data());
      TH1F *h=(TH1F*)dir->Get(plotname.Data());
      //if(idx1==0) firstHisto[idx2]=h;
      h->SetName( ( samplesMC[idx1]+"_"+h->GetName() ).Data() );
      if(stacked) {
	h->SetFillColor(cols[idx1]); 
	h->SetLineColor(kBlack);
      }
      else {
	h->SetLineColor(cols[idx1]);
	h->SetLineWidth(2);
      }
      h->Rebin(rebinValue);
      h->GetXaxis()->SetRangeUser(xbinMin, xbinMax);
      if(doNormalize) {
	h->Scale(1./h->Integral());
      }
      stacks[idx2]->Add(h);
      lleg[idx2]->AddEntry(h, proc[samplesMC[idx1]].Data(), "f");
      maxV[idx2]=(h->GetBinContent(h->GetMaximumBin())>maxV[idx2] ? h->GetBinContent(h->GetMaximumBin()) : maxV[idx2]);
      minV[idx2]=(h->GetBinContent(h->GetMaximumBin())<minV[idx2] && h->GetBinContent(h->GetMaximumBin())>1.e-7 ? h->GetBinContent(h->GetMaximumBin()) : minV[idx2]);
    } // end for( ... subdir ...)
  } // end for( ... samples ...)


  checkAndMakeDir(outdirname);

  TFile *outfile=new TFile((outdirname+"/"+outname).Data(), "UPDATE");

  vector<TH1F *> frames;
  for(int idx3=0; idx3<sizeSdir; ++idx3) {

    TDirectory *dir=(TDirectory*)outfile->Get(subdir[idx3].Data());
    if(dir==0) dir=outfile->mkdir(subdir[idx3].Data());
    dir->cd();

    int xbin=((TH1F*)stacks[0]->GetHists()->At(0))->GetXaxis()->GetNbins();
    float xmin=((TH1F*)stacks[0]->GetHists()->At(0))->GetXaxis()->GetXmin();
    float xmax=((TH1F*)stacks[0]->GetHists()->At(0))->GetXaxis()->GetXmax();
    float onebin=(xmax-xmin)/xbin;
    TString ylab("Entries/");
    ylab+=onebin;
    TString unit=xlab( xlab.Index("[")+1, ( xlab.Index("]")-xlab.Index("[")-1 ) );
    if(unit.Length()!=0) unit.Prepend(" ");
    ylab+=unit;
    if(ylab.EndsWith("/1")) ylab.ReplaceAll("/1","");
    if(logscale) {
      minV[idx3]/=2.;
      maxV[idx3]*=2.;
    }
    else {
      minV[idx3]=0.;
      maxV[idx3]*=1.2;
    }
    TH1F *fr=new TH1F( (subdir[idx3]+"_frame").Data(), (subdir[idx3]+"_frame").Data(), xbin, xmin, xmax);
    fr->SetMinimum(minV[idx3]);
    fr->SetMaximum(maxV[idx3]);
    fr->GetXaxis()->SetRangeUser(xbinMin, xbinMax);
    fr->GetXaxis()->SetTitle(xlab.Data());
    fr->GetXaxis()->SetTitleSize(0.06);
    //fr->GetXaxis()->SetTitleOffset(0.9);
    fr->GetYaxis()->SetTitle(ylab.Data());
    fr->GetYaxis()->SetTitleSize(0.06);
    fr->GetYaxis()->SetTitleOffset(1.3);
    frames.push_back(fr);
    ccs[idx3]->cd();
    fr->Draw();
    TString opt("histsame");
    if(!stacked) opt+="nostack";
    stacks[idx3]->Draw(opt.Data());
    if(plotData) dataplots[idx3]->Draw("samepe");
    lleg[idx3]->Draw();
    ccs[idx3]->SetLogy(logscale);
    ccs[idx3]->Modified();
    ccs[idx3]->Update();
    //ccs[idx3]->SaveAs( (TString(ccs[idx3]->GetName())+".png").Data() );
    TString finalname=ccs[idx3]->GetName();
    finalname.ReplaceAll("canvas_", "");
    ccs[idx3]->Write(finalname.Data(), TObject::kOverwrite);
    checkAndMakeDir( (outdirname+"/"+subdir[idx3]).Data() );
    ccs[idx3]->SaveAs((outdirname+"/"+subdir[idx3]+"/"+finalname+".png").Data());

    outfile->cd();
  }

}



// --------------------------------------------------------------- //
// --  Macro for plots in main directory  ------------------------ //
// --------------------------------------------------------------- //

void plotNumbers() {

  /*
  TString inputdir="/tmp/trocino/analysis_try/";
  TString namefile="ZZllvvAnalyzer_0.root";
  */

  /*
  vector<TString> samples={"WWtoAnything_Spring11","WZtoAnything_Spring11","ZZtoAnything_only2l2n_Spring11","ZZtoAnything_allBut2l2n_Spring11","TToBLNu_tW-channel_Spring11","TToBLNu_t-channel_Spring11","TToBLNu_s-channel_Spring11","TTJets_madgraph_Spring11","WJetsToLNu_Spring11","DYToEE_M-20_Spring11","DYToMuMu_M-20_Spring11","DYToTauTau_M-20_Spring11","QCD_Pt-30to50_Spring11","QCD_Pt-50to80_Spring11","QCD_Pt-80to120_Spring11","QCD_Pt-120to170_Spring11","QCD_Pt-170to300_Spring11","QCD_Pt-300to470_Spring11","QCD_Pt-470to600_Spring11","QCD_Pt-600to800_Spring11","QCD_Pt-800to1000_Spring11","DoubleMuon-v3","DoubleMuon-v5","DoubleMuon-v6","DoubleMuon-v7","DoubleMuon-v8","DoubleMuon-v9","DoubleElectron-v3","DoubleElectron-v5","DoubleElectron-v6","DoubleElectron-v7","DoubleElectron-v8","DoubleElectron-v9","MuEG-v3","MuEG-v5","MuEG-v6","MuEG-v7","MuEG-v8","MuEG-v9"};
  */

  /*
  TString samplesMC[]={"ZZtoAnything_only2l2n_Spring11","WZtoAnything_Spring11","WWtoAnything_Spring11","ZZtoAnything_allBut2l2n_Spring11","TTJets_madgraph_Spring11","TToBLNu_Spring11","WJetsToLNu_Spring11","DYToMuMu_M-20_Spring11","DYToTauTau_M-20_Spring11","DYToEE_M-20_Spring11","QCD_Pt-30to1000_Spring11"};
  int sizeMc=sizeof(samplesMC)/sizeof(TString);
  */

  /*
  // Only for 2mu channel
  vector<TString> samplesMC={"ZZtoAnything_only2l2n_Spring11","WZtoAnything_Spring11","WWtoAnything_Spring11","ZZtoAnything_allBut2l2n_Spring11","TTJets_madgraph_Spring11","TToBLNu_Spring11","WJetsToLNu_Spring11","QCD_Pt-30to1000_Spring11","DYToTauTau_M-20_Spring11","DYToMuMu_M-20_Spring11"};
  */

  /*
  Color_t cols[]={kRed, kPink, kMagenta, kViolet, kBlue, kAzure, kCyan, kTeal, kGreen, kSpring, kYellow, kOrange};
  */

  map<TString, TString> proc;
  initiateMap(proc);

  Style_t mrks[]={21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32};

  TString sampleData="AllData";

  TCanvas *ccc=new TCanvas( "ccc", "ccc", 600, 600);

  vector<TFile *> allFiles;

  // Old, not consistent with the output of "submission.py --copy"
  //  TFile *fdata=new TFile((inputdir+sampleData+"/"+namefile).Data(), "READ");

  // Assuming all data files (from DoubleMu, DoubleElectron, MuEG streams) 
  // have been merged into a single file named "ZZllvvAnalyzer_AllData.root"
  TString namefiletmp=namefile;
  namefiletmp.ReplaceAll("0", sampleData.Data());
  TFile *fdata=new TFile((inputdir+"/"+namefiletmp).Data(), "READ");
  TH1F *hdata=(TH1F*)fdata->Get(plotnamefix.Data());
  if(plotData) {
    hdata->SetMarkerStyle(20);
    hdata->SetMarkerSize(1.4);
    hdata->Rebin(rebinValue);
    hdata->GetXaxis()->SetRangeUser(xbinMin, xbinMax);
    if(doNormalize) {
      hdata->Scale(1./hdata->Integral());
    }
  }

  TLegend *leg=new TLegend(0.75, 0.6, 1., 1.);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetLineColor(0);
  if(plotData) leg->AddEntry(hdata, "data", "pe");

  float maxC=0.;
  float minC=0.;

  //  TH1F *firstPlot=0;
  THStack *stack=new THStack(("stack_"+plotnamefix).Data(), plotnamefix.Data()) ;

  for(int idx1=0; idx1<sizeMc; ++idx1) {
    // Old, not consistent with the output of "submission.py --copy"
    //    allFiles.push_back(new TFile((inputdir+samplesMC[idx1]+"/"+namefile).Data(), "READ"));
    TString namefilemctmp=namefile;
    namefilemctmp.ReplaceAll("0", samplesMC[idx1].Data());
    allFiles.push_back(new TFile((inputdir+"/"+namefilemctmp).Data(), "READ"));
    TFile *f=allFiles.back();
    TH1F *h=(TH1F*)f->Get(plotnamefix.Data());
    //if(idx1==0) firstPlot=h;
    h->SetName( ( samplesMC[idx1]+"_"+h->GetName() ).Data() );
    h->Rebin(rebinValue);
    h->GetXaxis()->SetRangeUser(xbinMin, xbinMax);
    if(doNormalize) {
      h->Scale(1./h->Integral());
    }
    if(stacked) {
      h->SetFillColor(cols[idx1]); 
      h->SetLineColor(kBlack);
    }
    else {
      h->SetLineColor(cols[idx1]);
      h->SetMarkerStyle(mrks[idx1]);
      h->SetMarkerColor(cols[idx1]);
      h->SetMarkerSize(1.4);
    }
    stack->Add(h);
    //ccc->cd();
    //if(idx1==0) h->Draw();
    //else h->Draw("same");
    leg->AddEntry(h, proc[samplesMC[idx1]].Data(), "pe");

    if(samplesMC[idx1].Contains("DYToMuMu")) {
      h->SetMarkerColor(kOrange+1);
    }
    maxC=(h->GetBinContent(h->GetMaximumBin())>maxC ? h->GetBinContent(h->GetMaximumBin()) : maxC);
    if(samplesMC[idx1].Contains("ZZtoAnything_only2l2n")) {
      if(plotnamefix.Contains("EventCount") && !plotnamefix.Contains("Pre"))
	minC=h->GetBinContent(6);
      else if(plotnamefix.Contains("PreEventCount"))
	minC=h->GetBinContent(4);
      else if(plotnamefix.Contains("Vertex"))
	minC=h->GetBinContent(h->GetMaximumBin());
    }
  } // end for( ... samples ...)

  if(plotnamefix.Contains("EventCount") && !plotnamefix.Contains("Pre")) {
    if(logscale) {
      minC/=10.;
      maxC*=10.;
    }
    else {
      minC=0.;
      maxC*=1.2;
    }
  }
  else if(plotnamefix.Contains("Vertex")) {
    if(logscale) {
      minC/=100.;
      maxC*=2.;
    }
    else {
      minC=0.;
      maxC*=1.2;
    }
  }

  TString ylabN("Events");

  TH1F *frame=new TH1F( (plotnamefix+"_frame").Data(), (plotnamefix+"_frame").Data(), 
			((TH1F*)stack->GetHists()->At(0))->GetXaxis()->GetNbins(), 
			((TH1F*)stack->GetHists()->At(0))->GetXaxis()->GetXmin(), 
			((TH1F*)stack->GetHists()->At(0))->GetXaxis()->GetXmax());
  frame->SetMinimum(minC);
  frame->SetMaximum(maxC);
  frame->GetXaxis()->SetRangeUser(xbinMin, xbinMax);
  frame->GetXaxis()->SetTitle(xlabN.Data());
  frame->GetXaxis()->SetTitleSize(0.06);
  //frame->GetXaxis()->SetTitleOffset(0.9);
  frame->GetYaxis()->SetTitle(ylabN.Data());
  frame->GetYaxis()->SetTitleSize(0.06);
  frame->GetYaxis()->SetTitleOffset(1.3);

  ccc->cd();
  frame->Draw();
  TString opt("same");
  if(drawAsHist) opt+="hist";
  if(!stacked) opt+="nostack";
  stack->Draw(opt.Data());
  if(plotData) hdata->Draw("same");
  leg->Draw();
  ccc->SetLogy(logscale);
  ccc->Modified();
  ccc->Update();

  checkAndMakeDir(outdirname);

  TFile *outfile=new TFile((outdirname+"/"+outname).Data(), "UPDATE");
  outfile->cd();
  ccc->Write(plotnamefix.Data(), TObject::kOverwrite);
  ccc->SaveAs((outdirname+"/"+plotnamefix+".png").Data());

}

