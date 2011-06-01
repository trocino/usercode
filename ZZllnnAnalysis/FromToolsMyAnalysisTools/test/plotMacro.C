#include "TString.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TDirectoryFile.h"
#include "TFile.h"
#include "TH1F.h"
#include <iostream>

using namespace std;

void plotMacro() {

  TString inputdir="/tmp/trocino/analysis_try/";
  TString namefile="ZZllvvAnalyzer_0.root";
  /*
  vector<TString> samples={"WWtoAnything_Spring11","WZtoAnything_Spring11","ZZtoAnything_only2l2n_Spring11","ZZtoAnything_allBut2l2n_Spring11","TToBLNu_tW-channel_Spring11","TToBLNu_t-channel_Spring11","TToBLNu_s-channel_Spring11","TTJets_madgraph_Spring11","WJetsToLNu_Spring11","DYToEE_M-20_Spring11","DYToMuMu_M-20_Spring11","DYToTauTau_M-20_Spring11","QCD_Pt-30to50_Spring11","QCD_Pt-50to80_Spring11","QCD_Pt-80to120_Spring11","QCD_Pt-120to170_Spring11","QCD_Pt-170to300_Spring11","QCD_Pt-300to470_Spring11","QCD_Pt-470to600_Spring11","QCD_Pt-600to800_Spring11","QCD_Pt-800to1000_Spring11","DoubleMuon-v3","DoubleMuon-v5","DoubleMuon-v6","DoubleMuon-v7","DoubleMuon-v8","DoubleMuon-v9","DoubleElectron-v3","DoubleElectron-v5","DoubleElectron-v6","DoubleElectron-v7","DoubleElectron-v8","DoubleElectron-v9","MuEG-v3","MuEG-v5","MuEG-v6","MuEG-v7","MuEG-v8","MuEG-v9"};
  */

  TString samplesMC[]={"ZZtoAnything_only2l2n_Spring11","WZtoAnything_Spring11","WWtoAnything_Spring11","ZZtoAnything_allBut2l2n_Spring11","TTJets_madgraph_Spring11","TToBLNu_Spring11","WJetsToLNu_Spring11","DYToMuMu_M-20_Spring11","DYToTauTau_M-20_Spring11","DYToEE_M-20_Spring11","QCD_Pt-30to1000_Spring11"};
  int sizeMc=sizeof(samplesMC)/sizeof(TString);

  /*
  // Only for 2mu channel
  vector<TString> samplesMC={"ZZtoAnything_only2l2n_Spring11","WZtoAnything_Spring11","WWtoAnything_Spring11","ZZtoAnything_allBut2l2n_Spring11","TTJets_madgraph_Spring11","TToBLNu_Spring11","WJetsToLNu_Spring11","QCD_Pt-30to1000_Spring11","DYToTauTau_M-20_Spring11","DYToMuMu_M-20_Spring11"};
  */

  Color_t cols[]={kRed, kPink, kMagenta, kViolet, kBlue, kAzure, kCyan, kTeal, kGreen, kSpring, kYellow, kOrange};

  TString sampleData="AllData";

  TString subdir[]={"cut1_Dilepton", "cut2_Mass-window", "cut3_RedMET-cut", "cut4_Jet-veto", "cut5_Lepton-veto"};
  int sizeSdir=sizeof(subdir)/sizeof(TString);

  //TString plotcategory="LeptonLead_"; TString plottype="_hPt";
  TString plotcategory="LeptonLead_"; TString plottype="_hPt";

  vector<THStack *> stacks;
  vector<TH1F *>    dataplots;
  vector<TCanvas *> ccs;
  vector<TLegend *> lleg;
  vector<TFile *>   allFiles;

  TFile *fdata=new TFile((inputdir+sampleData+"/"+namefile).Data(), "READ");

  //  for(vector<TString>::iterator sdr=subdir.begin(); sdr!=subdir.end(); ++sdr) {
  for(int j1=0; j1<sizeSdir; ++j1) {
    TString plotname=plotcategory+subdir[j1]+plottype;
    stacks.push_back(new THStack(("stack_"+subdir[j1]+"_"+plotname).Data(), plotname.Data()) );
    ccs.push_back( new TCanvas( ("canvas_"+subdir[j1]+"_"+plotname).Data(), ("canvas_"+subdir[j1]+"_"+plotname).Data(), 600, 600) );
    lleg.push_back(new TLegend(0.6, 0.4, 1., 1.));

    TDirectoryFile *dir=(TDirectoryFile*)fdata->Get(subdir[j1].Data());
    TH1F *h=(TH1F*)dir->Get(plotname.Data());
    h->SetName( (TString("data_")+h->GetName()).Data() );
    h->SetMarkerStyle(21);
    //h->Rebin(4);
    dataplots.push_back(h);
    lleg.back()->SetFillColor(0);
    lleg.back()->SetFillStyle(0);
    lleg.back()->SetBorderSize(0);
    lleg.back()->SetLineColor(0);
    lleg.back()->AddEntry(dataplots[0], "data", "pe");
  }

  int idx0(0);
  //  int idx1(0);
  //  for(vector<TString>::iterator smpl=samplesMC.begin(); smpl!=samplesMC.end(); ++smpl, ++idx1) {
  for(int idx1=0; idx1<sizeMc; ++idx1) {
    allFiles.push_back(new TFile((inputdir+samplesMC[idx1]+"/"+namefile).Data(), "READ"));
    TFile *f=allFiles.back();
    //    for(vector<TString>::iterator sdr=subdir.begin(); sdr!=subdir.end(); ++sdr, ++idx2) {
    for(int idx2=0; idx2<sizeSdir; ++idx2, ++idx0) {
      TDirectoryFile *dir=(TDirectoryFile*)f->Get(subdir[idx2].Data());
      TString plotname=plotcategory+subdir[idx2]+plottype;
      TH1F *h=(TH1F*)dir->Get(plotname.Data());
      h->SetName( ( samplesMC[idx1]+"_"+h->GetName() ).Data() );
      h->SetLineColor(cols[idx1]);
      h->SetFillColor(cols[idx1]);
      //h->Rebin(4);
      stacks[idx2]->Add(h);
      TString leglab(samplesMC[idx1].Data());
      leglab.ReplaceAll("_Spring11", "");
      lleg[idx2]->AddEntry(h, leglab.Data(), "lf");
    } // end for( ... subdir ...)
    //    f.Close();
  } // end for( ... samples ...)

  for(int idx3=0; idx3<sizeSdir; ++idx3) {
    ccs[idx3]->cd();
    stacks[idx3]->Draw("hist");
    dataplots[idx3]->Draw("samepe");
    lleg[idx3]->Draw();
    ccs[idx3]->SetLogy();
    ccs[idx3]->Modified();
    ccs[idx3]->Update();
    //ccs[idx3]->SaveAs( (TString(ccs[idx3]->GetName())+".png").Data() );
  }

}


void plotNumbers() {

  TString inputdir="/tmp/trocino/analysis_try/";
  TString namefile="ZZllvvAnalyzer_0.root";
  /*
  vector<TString> samples={"WWtoAnything_Spring11","WZtoAnything_Spring11","ZZtoAnything_only2l2n_Spring11","ZZtoAnything_allBut2l2n_Spring11","TToBLNu_tW-channel_Spring11","TToBLNu_t-channel_Spring11","TToBLNu_s-channel_Spring11","TTJets_madgraph_Spring11","WJetsToLNu_Spring11","DYToEE_M-20_Spring11","DYToMuMu_M-20_Spring11","DYToTauTau_M-20_Spring11","QCD_Pt-30to50_Spring11","QCD_Pt-50to80_Spring11","QCD_Pt-80to120_Spring11","QCD_Pt-120to170_Spring11","QCD_Pt-170to300_Spring11","QCD_Pt-300to470_Spring11","QCD_Pt-470to600_Spring11","QCD_Pt-600to800_Spring11","QCD_Pt-800to1000_Spring11","DoubleMuon-v3","DoubleMuon-v5","DoubleMuon-v6","DoubleMuon-v7","DoubleMuon-v8","DoubleMuon-v9","DoubleElectron-v3","DoubleElectron-v5","DoubleElectron-v6","DoubleElectron-v7","DoubleElectron-v8","DoubleElectron-v9","MuEG-v3","MuEG-v5","MuEG-v6","MuEG-v7","MuEG-v8","MuEG-v9"};
  */

  TString samplesMC[]={"ZZtoAnything_only2l2n_Spring11","WZtoAnything_Spring11","WWtoAnything_Spring11","ZZtoAnything_allBut2l2n_Spring11","TTJets_madgraph_Spring11","TToBLNu_Spring11","WJetsToLNu_Spring11","DYToMuMu_M-20_Spring11","DYToTauTau_M-20_Spring11","DYToEE_M-20_Spring11","QCD_Pt-30to1000_Spring11"};
  int sizeMc=sizeof(samplesMC)/sizeof(TString);

  /*
  // Only for 2mu channel
  vector<TString> samplesMC={"ZZtoAnything_only2l2n_Spring11","WZtoAnything_Spring11","WWtoAnything_Spring11","ZZtoAnything_allBut2l2n_Spring11","TTJets_madgraph_Spring11","TToBLNu_Spring11","WJetsToLNu_Spring11","QCD_Pt-30to1000_Spring11","DYToTauTau_M-20_Spring11","DYToMuMu_M-20_Spring11"};
  */

  Color_t cols[]={kRed, kPink, kMagenta, kViolet, kBlue, kAzure, kCyan, kTeal, kGreen, kSpring, kYellow, kOrange};
  Style_t mrks[]={21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32};

  TString sampleData="AllData";

  TString plotname="EventCounter";

  TCanvas *ccc=new TCanvas( "ccc", "ccc", 600, 600);

  vector<TFile *> allFiles;

  TFile *fdata=new TFile((inputdir+sampleData+"/"+namefile).Data(), "READ");
  TH1F *hdata=(TH1F*)fdata->Get(plotname.Data());
  hdata->SetMarkerStyle(20);

  TLegend *leg=new TLegend(0.6, 0.4, 1., 1.);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetLineColor(0);
  leg->AddEntry(hdata, "data", "pe");

  float maxC=0.;
  float minC=0.;

  TH1F *firstPlot=0;

  for(int idx1=0; idx1<sizeMc; ++idx1) {
    allFiles.push_back(new TFile((inputdir+samplesMC[idx1]+"/"+namefile).Data(), "READ"));
    TFile *f=allFiles.back();
    TH1F *h=(TH1F*)f->Get(plotname.Data());
    if(idx1==0) firstPlot=h;
    h->SetName( ( samplesMC[idx1]+"_"+h->GetName() ).Data() );
    h->SetMarkerStyle(mrks[idx1]);
    h->SetMarkerColor(cols[idx1]);
    h->SetLineColor(cols[idx1]);
    //h->SetFillColor(cols[idx1]);
    ccc->cd();
    if(idx1==0) h->Draw();
    else h->Draw("same");
    TString leglab(samplesMC[idx1].Data());
    leglab.ReplaceAll("_Spring11", "");
    leg->AddEntry(h, leglab.Data(), "pe");

    if(samplesMC[idx1].Contains("DYToMuMu")) 
      maxC=h->GetBinContent(1);
    if(samplesMC[idx1].Contains("ZZtoAnything_only2l2n"))
      minC=h->GetBinContent(6);
  } // end for( ... samples ...)

  maxC*=2.5;
  minC/=0.1;

  firstPlot->GetYaxis()->SetLimits(minC, maxC);

  ccc->cd();
  hdata->Draw("same");
  leg->Draw();
  ccc->SetLogy();
  ccc->Modified();
  ccc->Update();

}
