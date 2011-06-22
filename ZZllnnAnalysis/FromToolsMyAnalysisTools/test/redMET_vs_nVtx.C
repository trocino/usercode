void redMET_vs_nVtx() {

  TFile *f=new TFile("/tmp/trocino/AllData/ZZllvvAnalyzer_AllData.root", "READ");
  //TFile *f=new TFile("/tmp/trocino/AllData/ZZllvvAnalyzer_DYToMuMu_M-20_Spring11.root", "READ");
  //TH2F *h2=(TH2F*)f->Get("cut2_MassWindow/RedMetStd_cut2_MassWindow_hRedMETOverNVtxVsNVtx");
  TH2F *h2=(TH2F*)f->Get("cut2_MassWindow/RedMetStd_cut2_MassWindow_hRedMETVsNVtx");
  //TH2F *h2=(TH2F*)f->Get("cut2_MassWindow/RedMetStd_cut2_MassWindow_hRedMETCompLongOverNVtxVsNVtx");
  //TH2F *h2=(TH2F*)f->Get("cut2_MassWindow/RedMetStd_cut2_MassWindow_hRedMETCompLongVsNVtx");
  TCanvas *ccc=new TCanvas("ccc", "ccc", 600, 600);
  ccc->SetLogy();
  int firstBin=2;
  int step=1;
  const int nPlots=10;
  const int maxN=firstBin+nPlots*step;
  double nVtx[nPlots];
  double mean[nPlots];
  double meanErr[nPlots];
  double sigma[nPlots];
  double sigmaErr[nPlots];
  int cnt=0;
  TLegend *leg=new TLegend(0.2, 0.15, 0.35, 0.4);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetLineColor(0);
  for(int i=firstBin; i<maxN; i+=step, ++cnt) {
    TString name("redMET_");
    name+=(i-1);
    name+="Vtx";
    TH1D *h=h2->ProjectionY(name.Data(), i, i+step-1);
    h->Rebin(2);
    h->GetXaxis()->SetRange(2, h->FindBin(30.));
    h->Scale(1./h->Integral(2, 1));
    //h->Scale(1./h->Integral());
    h->SetLineWidth(2);
    if(i!=10) h->SetLineColor(i);
    else h->SetLineColor(30);
    h->SetFillColor(0);
    h->SetFillStyle(0);
    //h->GetXaxis()->SetRangeUser(0., 50.);
    if(i==firstBin) h->Draw("hist");
    else h->Draw("histsame");

    TString label("n_{vtx} = ");
    label+=(i-1);
    if(step!=1) {label+="-"; label+=(i+step-2);}
    leg->AddEntry(h, label.Data(), "l");

    // --- Fit --- //
    myfunc();
    TF1 *fct=gROOT->GetFunction("myfunc");
    if(i!=10) fct->SetLineColor(i);
    else fct->SetLineColor(30);
    fct->SetLineWidth(2);
    h->Fit("myfunc", "0R+", "", h->GetXaxis()->GetBinCenter(2), 30.);
    fct->Draw("same");

    nVtx[cnt]=((i-1)+(i+step-2))/2.;
    mean[cnt]=fct->GetParameter(1);
    meanErr[cnt]=fct->GetParError(1);
    sigma[cnt]=fct->GetParameter(2);
    sigmaErr[cnt]=fct->GetParError(2);
  }
  leg->Draw();

  TCanvas *ccc2=new TCanvas("ccc2", "ccc2", 600, 600);
  TGraphErrors *gre=new TGraphErrors(nPlots, nVtx, sigma, 0, sigmaErr);
  gre->SetMarkerStyle(20);
  gre->GetXaxis()->SetTitle("Number of reconstructed vertices");
  gre->GetYaxis()->SetTitle("Gaussian #sigma of redMET");
  gre->GetXaxis()->SetTitleSize(0.052);
  gre->GetYaxis()->SetTitleSize(0.052);
  gre->Draw("ape");
  /*
  gre->Fit("pol0", "+");
  TF1* f0=gre->GetFunction("pol0");
  f0->SetLineWidth(2);
  f0->SetLineColor(kRed);
  gre->Fit("pol1", "+");
  TF1* f1=gre->GetFunction("pol1");
  f1->SetLineWidth(2);
  f1->SetLineColor(kBlue);
  */

  TCanvas *ccc3=new TCanvas("ccc3", "ccc3", 600, 600);
  TGraphErrors *grem=new TGraphErrors(nPlots, nVtx, mean, 0, meanErr);
  grem->SetMarkerStyle(20);
  grem->GetXaxis()->SetTitle("Number of reconstructed vertices");
  grem->GetYaxis()->SetTitle("Gaussian mean of redMET");
  grem->GetXaxis()->SetTitleSize(0.052);
  grem->GetYaxis()->SetTitleSize(0.052);
  grem->Draw("ape");

  TCanvas *ccc4=new TCanvas("ccc4", "ccc4", 600, 600);
  TGraphErrors *grems=new TGraphErrors(nPlots, mean, sigma, meanErr, sigmaErr);
  grems->SetMarkerStyle(20);
  grems->GetXaxis()->SetTitle("Gaussian mean of redMET");
  grems->GetYaxis()->SetTitle("Gaussian #sigma of redMET");
  grems->GetXaxis()->SetTitleSize(0.052);
  grems->GetYaxis()->SetTitleSize(0.052);
  grems->Draw("ape");
}

Double_t myfunction(Double_t *x, Double_t *par) {
  Float_t n=x[0];
  //Double_t f=par[0]*TMath::Gaus(n, par[1], par[2]) + par[3];
  //Double_t f = par[0]*TMath::Exp(-1.*par[1]*n);
  Double_t f=par[0]*TMath::Gaus(n, par[1], par[2]);
    /*+ par[3]*TMath::Exp(-1.*par[4]*n);*/
  return f;
}

void myfunc() {
  TF1 *f1 = new TF1("myfunc", myfunction, 0., 100., 3);
  f1->SetParameter(0, 5.6);
  f1->SetParameter(1, -15.);
  f1->SetParameter(2, 9.5);
  //f1->SetParameter(3, 0.002);
  //f1->SetParameter(4, 0.1);
  f1->SetParNames("N_{1}", "#mu", "#sigma");
}

