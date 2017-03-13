#include "TDatime.h" 
#include "TStyle.h" 
#include "TString.h" 
#include "TSystem.h" 
#include "TFile.h" 
#include "TH1F.h" 
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TLegend.h" 
#include "TEfficiency.h"
#include <stdlib.h>
#include <iostream>

//#include "/afs/cern.ch/user/t/trocino/Utilities/danystyle.C" 

void mergeRootFiles(TString folder = "") {

  if(folder.Length()==0) { 
    std::cout << " *** ERROR: missing input folder! *** " << std::endl; 
    return; 
  } 

  TDatime tdt; 
  //TString outpath; 
  TString outpath = "/afs/cern.ch/user/t/trocino/www/Temp/"; 
  outpath += tdt.GetDate(); 
  outpath += "_"; 
  outpath += tdt.GetTime(); 

  //  TSystem tus; 
  gSystem->MakeDirectory(outpath.Data()); 
  TString cmmd = "cp /afs/cern.ch/user/t/trocino/www/Temp/index.php "; 
  gSystem->Exec( (cmmd+outpath).Data()); 

  TString segmentcoll[] = { 
    //// Strip: 768; Eta partitions: 4, 6, 8, 16 
    "s768p004", 
    "s768p006", 
    "s768p008", 
    "s768p012", 
    "s768p016", 
    //"s768p032", 
    //"s768p064", 
    //// Strip: 512; Eta partitions: 4, 6, 8, 16 
    "s512p004", 
    "s512p006", 
    "s512p008", 
    "s512p012", 
    "s512p016", 
    //"s512p032", 
    //"s512p064", 
    //// Strip: 384; Eta partitions: 4, 6, 8, 16 
    "s384p004", 
    "s384p006", 
    "s384p008", 
    "s384p012", 
    "s384p016", 
    //"s384p032", 
    //"s384p064", 
    //// Strip: 256; Eta partitions: 4, 6, 8, 16 
    "s256p004", 
    "s256p006", 
    "s256p008", 
    "s256p012", 
    "s256p016" //, 
    //"s256p032", 
    //"s256p064" 
  }; 

  if(sizeof(segmentcoll)==0) { 
    std::cout << " *** ERROR: segmentcoll can't be empty! Add {""}! *** " << std::endl; 
    return; 
  } 
  // static const size_t nrows = sizeof(segmentcoll)/sizeof(segmentcoll[0]); 

  double params[][5] = {
    // {-1.0,  4.0, -1.0, -1.0, -1.0}, //  0 
    // { 3.0,  4.0, -1.0, -1.0, -1.0}, //  1 
    // { 2.0,  4.0, -1.0, -1.0, -1.0}, //  2 
    // { 1.0,  4.0, -1.0, -1.0, -1.0}, //  3 
    // { 0.5,  4.0, -1.0, -1.0, -1.0}, //  4 
    // { 0.1,  4.0, -1.0, -1.0, -1.0}, //  5 

    // { 3.0, -1.0, -1.0, -1.0, -1.0}, //  6 
    // { 3.0,  3.0, -1.0, -1.0, -1.0}, //  7 
    // { 3.0,  2.0, -1.0, -1.0, -1.0}, //  8 
    // { 3.0,  1.0, -1.0, -1.0, -1.0}, //  9 
    // { 3.0,  0.5, -1.0, -1.0, -1.0}, // 10 

    // { 3.0,  4.0, 10.0, -1.0, -1.0}, // 11 
    // { 3.0,  4.0,  5.0, -1.0, -1.0}, // 12 
    // { 3.0,  4.0,  4.0, -1.0, -1.0}, // 13 
    // { 3.0,  4.0,  3.0, -1.0, -1.0}, // 14 
    // { 3.0,  4.0,  2.0, -1.0, -1.0}, // 15 
    // { 3.0,  4.0,  1.0, -1.0, -1.0}, // 16 
    // { 3.0,  4.0,  0.5, -1.0, -1.0}, // 17 

    // { 3.0,  4.0, -1.0, 20.0, -1.0}, // 18 
    // { 3.0,  4.0, -1.0, 10.0, -1.0}, // 19 
    // { 3.0,  4.0, -1.0,  5.0, -1.0}, // 20 
    // { 3.0,  4.0, -1.0,  4.0, -1.0}, // 21 
    // { 3.0,  4.0, -1.0,  3.0, -1.0}, // 22 
    // { 3.0,  4.0, -1.0,  2.0, -1.0}, // 23 
    // { 3.0,  4.0, -1.0,  1.0, -1.0}, // 24 

    // { 3.0,  4.0, -1.0, -1.0,  1.6}, // 25 
    // { 3.0,  4.0, -1.0, -1.0,  0.8}, // 26 
    // { 3.0,  4.0, -1.0, -1.0,  0.4}, // 27 
    // { 3.0,  4.0, -1.0, -1.0,  0.2}, // 28 
    // { 3.0,  4.0, -1.0, -1.0,  0.1}  // 29 
  }; 

  static const size_t nsegmts = sizeof(segmentcoll)==0 ? 0 : sizeof(segmentcoll)/sizeof(segmentcoll[0]); 
  static const size_t nparams = sizeof(params)==0 ? 0 : sizeof(params)/sizeof(params[0]); 
  static const size_t nrows = nparams>0 ? (nsegmts>0 ? nsegmts*nparams : nparams) : nsegmts; 

  std::map<Int_t, std::vector<Int_t> > striproll; 
  std::vector<std::pair<Int_t, Int_t> > strstriproll; 

  TString type[] = {"eff", "fake", "bkg"}; 
  static const size_t ntypes = sizeof(type)/sizeof(type[0]); 
  TString var[] = {"aeta", "phi", "pt", "allpu"}; 
  static const size_t nvars = sizeof(var)/sizeof(var[0]); 

  TH1F *histos_den[nrows][ntypes][nvars]; // 2: efficiency/fake rate/background rate;  3: eta, aeta, phi, pt, allpu 
  TH1F *histos_num[nrows][ntypes][nvars]; // 2: efficiency/fake rate/background rate;  3: eta, aeta, phi, pt, allpu 

  // Event/segment counters 
  TH1F *h_evt_seg = 0; 

  // Avarage fakes/background per event/segment 
  TH1F *h_fake_avg_perevt[nrows][ntypes-1][nvars]; 

  // Number of muons per segment 
  TH1F *h_mu_perseg_den_aeta[nrows]; 
  TH1F *h_mu_perseg_den_dphi[nrows]; 
  //TH1F *h_mu_perseg_num_aeta[nrows]; 
  //TH1F *h_mu_perseg_num_dphi[nrows]; 
  TH1F *h_mu_perseg_aeta[nrows]; 
  TH1F *h_mu_perseg_dphi[nrows]; 

  TString outdirs[nrows]; 

  TString pref = ""; 
  TString path = ""; 
  path += (pref + folder + "/"); 

  std::vector<TString> filenames; 
  size_t nfiles = 0; 

  int ipos = 0; // Type has to be "Ssiz_t", a.k.a. "int" 
  TString cmmd2 = "ls "; 
  cmmd2 += path; 
  cmmd2 += " | grep root"; 
  TString filelist = gSystem->GetFromPipe(cmmd2.Data()); 
  TString thisfile; 

  while(filelist.Tokenize(thisfile, ipos, "\n")) {
    filenames.push_back(path+thisfile); 
    ++nfiles; 
  } 

  TFile *fin(0); 

  // Loop over files 
  printf("Reading files...\n"); 
  size_t ndashes = 100; 
  for(size_t idash=1; idash<=ndashes; ++idash) printf("-"); 
  printf("\n"); 

  for(size_t i=0; i<nfiles; ++i) {
    fin = 0; 
    fin = TFile::Open( (filenames[i]).Data() );
    if(fin==0 || fin->IsZombie()) {
      std::cout << "  File " << filenames[i].Data() << " not found!" << std::endl;
      return; 
    }
    else {
      //std::cout << "  File " << filenames[i].Data() << " found!" << std::endl;
    }

    // Loop over folders (i.e. selection paramters) 
    //for(size_t j=0; j<nrows; ++j) { 
    for(size_t j0=0; j0<nsegmts; ++j0) { 
    for(size_t j1=0; j1<nparams; ++j1) { 
      size_t j = j0*nparams + j1; 
      // Only for first file (i = 0) 
      if(i==0) { 
	TString dirname = ""; 
	if(nparams>0) { 
	  // dirname.Form("PullX%.1fDiffX%.1fPullY%.1fDiffY%.1fDiffPhi%.1f", params[j1][0], params[j1][1], params[j1][2], params[j1][3], params[j1][4]); 
	  dirname.Form("%.1fDiffX%.1fPullY%.1fDiffY%.1fDiffPhi%.1f", params[j1][0], params[j1][1], params[j1][2], params[j1][3], params[j1][4]); 
	  dirname.ReplaceAll(".", "p"); 
	  dirname.ReplaceAll("-", "m"); 
	  //if(dirname.EndsWith("0")==true && dirname.EndsWith("m1p0")==false) dirname.Remove(TString::kTrailing, '0'); 
	} 
	//outdirs[j] = dirname; 
	TString srstr = segmentcoll[j0](segmentcoll[j0].Index('s')+1, segmentcoll[j0].Index('p')-segmentcoll[j0].Index('s')-1); 
	int nstrips = srstr.Atoi(); 
	srstr = segmentcoll[j0](segmentcoll[j0].Index('p')+1, segmentcoll[j0].Length()-segmentcoll[j0].Index('p')-1); 
	int nrolls = srstr.Atoi(); 
	striproll[nstrips].push_back(nrolls); 
	strstriproll.push_back(std::pair<Int_t, Int_t>(nstrips, nrolls)); 

	outdirs[j] = TString("me0Ana")+segmentcoll[j0]+dirname; 

	gSystem->MakeDirectory( (outpath+"/"+outdirs[j]).Data() ); 
	gSystem->Exec( (cmmd+outpath+"/"+outdirs[j]).Data()); 
      } 

      // Only for first folder (j = 0), but for every file (i) 
      if(j==0) { 
	TH1F *htmp = (TH1F*)fin->Get( (outdirs[j]+"/evt_seg_counter").Data() ); 
	if(htmp==0) {
	  std::cout << "  Histo " << (outdirs[j]+"/evt_seg_counter").Data() << " (file " << i << ") not found!" << std::endl;
	  return; 
	} 
	if(i==0) {
	  h_evt_seg = (TH1F*)htmp->Clone("new_seg_evt"); 
	  h_evt_seg->SetDirectory(0); 
	} 
	else { 
	  h_evt_seg->Add(htmp); 
	} 
      } 

      // N. muons per event/segment 
      TH1F *htmp_aeta_den = (TH1F*)fin->Get( (outdirs[j]+"/nmu_perseg_den_aeta").Data() ); 
      TH1F *htmp_dphi_den = (TH1F*)fin->Get( (outdirs[j]+"/nmu_perseg_den_dphi").Data() ); 
      TH1F *htmp_aeta_num = (TH1F*)fin->Get( (outdirs[j]+"/nmu_perseg_num_aeta").Data() ); 
      TH1F *htmp_dphi_num = (TH1F*)fin->Get( (outdirs[j]+"/nmu_perseg_num_dphi").Data() ); 
      if(htmp_aeta_den==0 || htmp_dphi_den==0 || htmp_aeta_num==0 || htmp_dphi_num==0) {
	std::cout << "  Some histos among " 
		  << (outdirs[j]+"/[nmu_perseg_den_aeta|nmu_perseg_den_dphi|nmu_perseg_num_aeta|nmu_perseg_num_dphi]").Data() 
		  << " (file " << i << ") not found!" << std::endl 
		  << "  Files: [" << htmp_aeta_den << "|" << htmp_dphi_den << "|" << htmp_aeta_num << "|" << htmp_dphi_num << "]" 
		  << std::endl; 
	return; 
      } 
      if(i==0) {
	h_mu_perseg_den_aeta[j] = (TH1F*)htmp_aeta_den->Clone( (TString("mu_perseg_den_aeta_")+outdirs[j]).Data() ); 
	h_mu_perseg_den_dphi[j] = (TH1F*)htmp_dphi_den->Clone( (TString("mu_perseg_den_dphi_")+outdirs[j]).Data() ); 
	h_mu_perseg_aeta[j] = (TH1F*)htmp_aeta_num->Clone( (TString("mu_perseg_aeta_")+outdirs[j]).Data() ); 
	h_mu_perseg_dphi[j] = (TH1F*)htmp_dphi_num->Clone( (TString("mu_perseg_dphi_")+outdirs[j]).Data() ); 
	h_mu_perseg_den_aeta[j]->SetDirectory(0); 
	h_mu_perseg_den_dphi[j]->SetDirectory(0); 
	h_mu_perseg_aeta[j]->SetDirectory(0); 
	h_mu_perseg_dphi[j]->SetDirectory(0); 
	// h_mu_perseg_den_aeta[j]->Sumw2(); 
	// h_mu_perseg_den_dphi[j]->Sumw2(); 
	// h_mu_perseg_aeta[j]->Sumw2(); 
	// h_mu_perseg_dphi[j]->Sumw2(); 
      } 
      else { 
	h_mu_perseg_den_aeta[j]->Add(htmp_aeta_den); 
	h_mu_perseg_den_dphi[j]->Add(htmp_dphi_den); 
	h_mu_perseg_aeta[j]->Add(htmp_aeta_num); 
	h_mu_perseg_dphi[j]->Add(htmp_dphi_num); 
      } 

      // Loop over types (efficiency, fake rate, background) 
      for(size_t k=0; k<ntypes; ++k) { 

	// Loop over vars (eta, aeta, phi, pt) 
	for(size_t w=0; w<nvars; ++w) { 

	  TString dentype = type[k]; 
	  if(dentype.EqualTo("bkg")) dentype = "fake"; 

	  //TH1F *htmp = (TH1F*)fin->Get( (outdirs[j]+"/"+type[k]+"_den_"+var[w]).Data() ); 
	  TH1F *htmp = (TH1F*)fin->Get( (outdirs[j]+"/"+dentype+"_den_"+var[w]).Data() ); 
	  if(htmp==0) {
	    std::cout << "  Histo " << (outdirs[j]+"/"+dentype+"_den_"+var[w]).Data() << " (file " << i << ") not found!" << std::endl;
	    return; 
	  } 
	  else { 
	    //std::cout << "  Histo " << (dirname+"/"+type+"_den_"+var).Data() << " (" << i << ") found!" << std::endl;
	    // if(w==1) { 
	    //   printf("  ==== File %3d, plot %15s: n. bins: %2d\n", int(i), htmp->GetName(), htmp->GetNbinsX()); 
	    // } 
	  } 

	  if(i==0) {
	    //std::cout << " --- Cloning *den* ... "; 
	    histos_den[j][k][w] = (TH1F*)htmp->Clone( (outdirs[j]+"_"+dentype+"_den_"+var[w]).Data() );
	    histos_den[j][k][w]->SetDirectory(0); 
	    //std::cout << "cloned!" << std::endl; 
	  } 
	  else { 
	    //std::cout << " --- Adding *den* ... "; 
	    histos_den[j][k][w]->Add(htmp); 
	    //std::cout << "added!" << std::endl; 
	  } 
	  htmp = 0; 
    
	  htmp = (TH1F*)fin->Get( (outdirs[j]+"/"+type[k]+"_num_"+var[w]).Data() ); 
	  if(htmp==0) {
	    std::cout << "  Histo " << (outdirs[j]+"/"+type[k]+"_num_"+var[w]).Data() << " (file " << i << ") not found!" << std::endl;
	    return; 
	  } 
	  else {
	    //std::cout << "  Histo " << (dirname+"/"+type+"_num_"+var).Data() << " (" << i << ") found!" << std::endl;
	  } 


	  if(i==0) {
	    //std::cout << " --- Cloning *num* ... "; 
	    histos_num[j][k][w] = (TH1F*)htmp->Clone( (outdirs[j]+"_"+type[k]+"_num_"+var[w]).Data() ); 
	    histos_num[j][k][w]->SetDirectory(0); 
	    //std::cout << "cloned!" << std::endl; 

	    // For average fakes/background 
	    if(k>0) { // only for fakes and background 
	      h_fake_avg_perevt[j][k-1][w] = (TH1F*)htmp->Clone( (outdirs[j]+"_"+type[k]+"_avg_"+var[w]).Data() ); 
	      h_fake_avg_perevt[j][k-1][w]->SetDirectory(0); 
	    } 
	  }
	  else {
	    //std::cout << " --- Adding *num* ... "; 
	    histos_num[j][k][w]->Add(htmp); 
	    //std::cout << "added!" << std::endl; 

	    // For average fakes/background 
	    if(k>0) { // only for fakes and background 
	      h_fake_avg_perevt[j][k-1][w]->Add(htmp); 
	    } 
	  } 

	} // End loop over vars (eta, aeta, phi, pt): for(size_t w=0; w<nvars; ++w)
      } // End loop over types (efficiency, fake rate): for(size_t k=0; k<ntypes; ++k)    
    } // End loop over folders (i.e. selection paramters): for(size_t j=0; j<nrows; ++j)
    } 

    fin->Close(); 

    for(size_t idash=1; idash<=ndashes; ++idash) { 
      size_t eqmax = (i+1)*ndashes/nfiles; 
      if(idash<=eqmax) printf("="); 
      else            printf("-"); 
    } 
    printf("\n"); 

  } //  End loop over files: for(size_t i=0; i<nfiles; ++i)  


  //tdrStyle->SetPadTopMargin(0.07);
  //tdrStyle->SetPadBottomMargin(0.13);
  tdrStyle->SetPadLeftMargin(0.10);
  tdrStyle->SetPadRightMargin(0.30); // (0.04);

  TCanvas *cc = new TCanvas("cc", "cc", 1000, 700); 
  cc->Draw();

  TCanvas *cc2 = new TCanvas("cc2", "cc2", 1000, 700); 
  cc2->Draw();

  TCanvas *ccall = new TCanvas("ccall", "ccall", 1000, 700); 

  float legxmin(0.75), legxmax(0.95), legymin(0.14), legymax(0.34); 


  float xxx[nrows], yyy[nrows], zzz[nrows], exl[nrows], exh[nrows], eyl[nrows], eyh[nrows], ezl[nrows], ezh[nrows]; 
  //Color_t cols[9] = {kBlack, kBlue, kViolet, kMagenta-9, kRed, kOrange+1, kYellow+2, kGreen+2, kCyan+2}; 
  //Style_t stls[9] = {20, 21, 22, 23, 24, 25, 26, 27, 28}; 
  Color_t cols[5] = {kBlue, kMagenta-9, kRed, kOrange+1, kGreen+2}; 
  Style_t stls[4] = {20, 24, 33, 27}; 
  Style_t stll[4] = {1, 7, 2, 3}; 

  for(size_t k=0; k<ntypes; ++k) { // Loop over types (efficiency, fake rate, background) 
    for(size_t w=0; w<nvars; ++w) { // Loop over vars (eta, aeta, phi, pt) 
      TString finname = (type[k]+"_"+var[w]); 

      ccall->Clear(); 

      // k: type[] = {"eff", "fake", "bkg"} 
      // w:  var[] = {"eta", "aeta", "phi", "pt", "allpu"} 

      /*
      if( (w==0 && k==0) || (w==0 && k==1) ) { 
	legxmin = 0.42; 
	legxmax = 0.68; 
	legymin = 0.12; 
	legymax = 0.40; 
      } 
      else if( (w==1 && k==0) || (w==3 && k==1) ) { 
	legxmin = 0.15; 
	legxmax = 0.40; 
	legymin = 0.65; 
	legymax = 0.95; 
      } 
      else if( (w==1 && k==1) ) { 
	legxmin = 0.72; 
	legxmax = 0.98; 
	legymin = 0.65; 
	legymax = 0.95; 
      } 
      else { 
	legxmin = 0.72; 
	legxmax = 0.98; 
	legymin = 0.12; 
	legymax = 0.40; 
      } 
      */ 

      legxmin = 0.72; 
      legxmax = 1.00; 
      legymin = 0.12; 
      legymax = 0.95; 

      TLegend *legall = new TLegend(legxmin, legymin, legxmax, legymax); 
      legall->SetFillColor(kWhite);
      legall->SetFillStyle(1001);
      legall->SetLineStyle(0);
      legall->SetLineWidth(0);
      legall->SetBorderSize(0);
      //legall->SetTextFont(42);

      std::vector<TGraphAsymmErrors*> tgaes, tgaes2, tgaes3, tgaes4; 

      // MultiGraphs 
      TMultiGraph *mgall  = new TMultiGraph(); 
      TMultiGraph *mgall2 = new TMultiGraph(); 
      TMultiGraph *mgall3 = new TMultiGraph(); 
      TMultiGraph *mgall4 = new TMultiGraph(); 

      //for(size_t j=0; j<nrows; ++j) { // Loop over folders (i.e. selection paramters) 
      for(size_t j0=0; j0<nsegmts; ++j0) { 
      for(size_t j1=0; j1<nparams; ++j1) { 
	size_t j = j0*nparams + j1; 
	// Let's start with the number of muons per segment 
	if(k==0 && w==0) { // Only once... 
	  h_mu_perseg_aeta[j]->Divide(h_mu_perseg_den_aeta[j]); 
	  h_mu_perseg_aeta[j]->SetNameTitle( (TString("nmu_avg_perseg_aeta_")+outdirs[j]).Data(), (TString("nmu_avg_perseg_aeta_")+outdirs[j]).Data() ); 
	  cc->cd(); 
	  tgaes3.push_back(new TGraphAsymmErrors(h_mu_perseg_aeta[j])); 
	  TGraphAsymmErrors *tgae3 = tgaes3.back(); 
	  tgae3->SetLineWidth(2); 
	  tgae3->SetLineColor(cols[j%5]); 
	  tgae3->SetMarkerSize(1.0); 
	  tgae3->SetMarkerColor(cols[j%5]); 
	  tgae3->SetMarkerStyle(stls[(j/5)%4]); 
	  tgae3->SetLineStyle(stll[(j/5)%4]); 
	  tgae3->GetXaxis()->SetTitle("Segment |#it{#eta}|"); 
	  tgae3->GetXaxis()->SetTitleOffset(0.84); 
	  tgae3->GetYaxis()->SetTitle("Muons/segment"); 
	  tgae3->GetYaxis()->SetTitleOffset(0.60); 
	  tgae3->Draw("APL"); 
	  cc->SaveAs( (outpath+"/"+outdirs[j]+"/avgMuons_perSeg_aeta.png").Data() ); 

	  h_mu_perseg_dphi[j]->Divide(h_mu_perseg_den_dphi[j]); 
	  //h_mu_perseg_dphi[j]->GetXaxis()->SetRangeUser(0., 1.); 
	  h_mu_perseg_dphi[j]->SetNameTitle( (TString("nmu_avg_perseg_dphi_")+outdirs[j]).Data(), (TString("nmu_avg_perseg_dphi_")+outdirs[j]).Data() ); 
	  cc->cd(); 
	  tgaes4.push_back(new TGraphAsymmErrors(h_mu_perseg_dphi[j])); 
	  TGraphAsymmErrors *tgae4 = tgaes4.back(); 
	  tgae4->SetLineWidth(2); 
	  tgae4->SetLineColor(cols[j%5]); 
	  tgae4->SetMarkerSize(1.0); 
	  tgae4->SetMarkerColor(cols[j%5]); 
	  tgae4->SetMarkerStyle(stls[(j/5)%4]); 
	  tgae4->SetLineStyle(stll[(j/5)%4]); 
	  tgae4->GetXaxis()->SetTitle("Track-segment #Delta#phi"); 
	  tgae4->GetXaxis()->SetTitleOffset(0.84); 
	  tgae4->GetYaxis()->SetTitle("Muons/segment"); 
	  tgae4->GetYaxis()->SetTitleOffset(0.60); 
	  tgae4->Draw("APL"); 
	  cc->SaveAs( (outpath+"/"+outdirs[j]+"/avgMuons_perSeg_dphi.png").Data() ); 
	} 


	if(var[w].EqualTo("phi")) { 
	  //histos_den[j][k][w]->Rebin(10); 
	  //histos_num[j][k][w]->Rebin(10); 
	  //if(k>0) h_fake_avg_perevt[j][k-1][w]->Rebin(10); 
	} 
	else if(var[w].EqualTo("pt")) { 
	  //double pt_bins[] = {0., 5., 10., 15., 20., 25., 30., 40., 50., 60., 80., 200.}; 
	  double pt_bins[] = {0., 1., 2., 3., 4., 5., 7.5, 10., 15., 20., 35., 50., 100., 200.}; 
	  int n_pt_bins = sizeof(pt_bins)/sizeof(double)-1; 
	  TString newdenname = (outdirs[j]+"_"+type[k]+"_den_"+var[w]+"_rebin"); 
	  TString newnumname = (outdirs[j]+"_"+type[k]+"_num_"+var[w]+"_rebin"); 
	  TString newavgname; 
	  TH1F *hdentmp = new TH1F(newdenname.Data(), newdenname.Data(), n_pt_bins, pt_bins); 
	  TH1F *hnumtmp = new TH1F(newnumname.Data(), newnumname.Data(), n_pt_bins, pt_bins); 
	  TH1F *havgtmp = 0; 
	  if(k>0) { 
	    newavgname = (outdirs[j]+"_"+type[k]+"_avg_"+var[w]+"_rebin"); 
	    havgtmp = new TH1F(newavgname.Data(), newavgname.Data(), n_pt_bins, pt_bins); 
	  } 
	  for(int ibin=0; ibin<=histos_den[j][k][w]->GetNbinsX(); ++ibin) { 
	    // Den 
	    int thisbin = hdentmp->FindBin( histos_den[j][k][w]->GetBinCenter(ibin) ); 
	    float thiscont = hdentmp->GetBinContent(thisbin) + histos_den[j][k][w]->GetBinContent(ibin); 
	    float thiserr = hdentmp->GetBinError(thisbin); 
	    thiserr *= thiserr; 
	    thiserr += histos_den[j][k][w]->GetBinError(ibin)*histos_den[j][k][w]->GetBinError(ibin); 
	    thiserr = sqrt(thiserr); 
	    hdentmp->SetBinContent(thisbin, thiscont); 
	    hdentmp->SetBinError(thisbin, thiserr); 
	    // Num 
	    thiscont = hnumtmp->GetBinContent(thisbin) + histos_num[j][k][w]->GetBinContent(ibin); 
	    thiserr = hnumtmp->GetBinError(thisbin); 
	    thiserr *= thiserr; 
	    thiserr += histos_num[j][k][w]->GetBinError(ibin)*histos_num[j][k][w]->GetBinError(ibin); 
	    thiserr = sqrt(thiserr); 
	    hnumtmp->SetBinContent(thisbin, thiscont); 
	    hnumtmp->SetBinError(thisbin, thiserr); 
	    // Avg 
	    if(k>0) { 
	      thiscont = havgtmp->GetBinContent(thisbin) + h_fake_avg_perevt[j][k-1][w]->GetBinContent(ibin); 
	      thiserr = havgtmp->GetBinError(thisbin); 
	      thiserr *= thiserr; 
	      thiserr += h_fake_avg_perevt[j][k-1][w]->GetBinError(ibin)*h_fake_avg_perevt[j][k-1][w]->GetBinError(ibin); 
	      thiserr = sqrt(thiserr); 
	      havgtmp->SetBinContent(thisbin, thiscont); 
	      havgtmp->SetBinError(thisbin, thiserr); 	      
	    } 
	  } 
	  histos_den[j][k][w] = hdentmp;                  // Bye bye, old plot! 
	  histos_num[j][k][w] = hnumtmp;                  // Bye bye, old plot! 
	  if(k>0) h_fake_avg_perevt[j][k-1][w] = havgtmp; // Bye bye, old plot! 
	} 
	tgaes.push_back(new TGraphAsymmErrors()); 
	TGraphAsymmErrors *tgae = tgaes.back();  // this is for "efficiency-like" division (efficiency, fake or bkg rate) 
	tgae->SetNameTitle( (finname+"_"+outdirs[j]).Data(), (finname+"_"+outdirs[j]).Data()); 
	tgae->Divide(histos_num[j][k][w], histos_den[j][k][w], "cl=0.683 b(1,1) mode" ); 
	tgae->SetLineWidth(2); 
	tgae->SetLineColor(cols[j%5]); 
	tgae->SetMarkerSize(1.0); 
	tgae->SetMarkerColor(cols[j%5]); 
	tgae->SetMarkerStyle(stls[(j/5)%4]); 
	tgae->SetLineStyle(stll[(j/5)%4]); 
	cc->cd(); 
	tgae->Draw("APL"); 

	TGraphAsymmErrors *tgae2 = 0; 
	if(k>0) { // only for fakes and background 
	  //h_fake_avg_perevt[j][k-1][w]->Scale(1./h_evt_seg->GetBinContent(1)); 
	  TH1F *h_fake_avg_perevt_tmp_den = (TH1F*)h_fake_avg_perevt[j][k-1][w]->Clone((TString(h_fake_avg_perevt[j][k-1][w]->GetName())+"_tmp_den").Data()); 
	  h_fake_avg_perevt_tmp_den->Reset("ICES"); 
	  for(size_t ibin=1; ibin<=h_fake_avg_perevt_tmp_den->GetNbinsX(); ++ibin) { 
	    if(h_fake_avg_perevt[j][k-1][w]->GetBinContent(ibin)==0.) continue; 
	    h_fake_avg_perevt_tmp_den->SetBinContent(ibin, h_evt_seg->GetBinContent(1)); 
	    h_fake_avg_perevt_tmp_den->SetBinError(ibin, h_evt_seg->GetBinError(1)); 
	  } 
	  tgaes2.push_back(new TGraphAsymmErrors()); 
	  tgae2 = tgaes2.back(); 
	  tgae2->Divide(h_fake_avg_perevt[j][k-1][w], h_fake_avg_perevt_tmp_den, "pois"); 
	  tgae2->SetNameTitle( (finname+"_perevt_"+outdirs[j]).Data(), (finname+"_perevt_"+outdirs[j]).Data()); 
	  tgae2->SetLineWidth(2); 
	  tgae2->SetLineColor(cols[j%5]); 
	  tgae2->SetMarkerSize(1.0); 
	  tgae2->SetMarkerColor(cols[j%5]); 
	  tgae2->SetMarkerStyle(stls[(j/5)%4]); 
	  tgae2->SetLineStyle(stll[(j/5)%4]); 
	  cc2->cd(); 
	  tgae2->Draw("APL"); 
	} 

	if(var[w].EqualTo("pt")) { 
	  tgae->GetXaxis()->SetTitle("Muon #it{p}_{T} [GeV]"); 
	  tgae->GetXaxis()->SetLimits(0., 35.); 
	  if(k>0) { 
	    tgae->GetXaxis()->SetLimits(0., 200.); 
	    tgae2->GetXaxis()->SetLimits(0., 200.); 
	    cc->SetLogx(); 
	    cc2->SetLogx(); 
	    tgae2->GetXaxis()->SetTitle("Muon #it{p}_{T} [GeV]"); 
	    //tgae2->GetXaxis()->SetLimits(0.3, 200.); 
	  } 
	} 

	else if(var[w].EqualTo("eta")) { 
	  tgae->GetXaxis()->SetTitle("Muon #it{#eta}"); 
	  if(k>0) tgae2->GetXaxis()->SetTitle("Muon #it{#eta}"); 
	} 
	else if(var[w].EqualTo("aeta")) { 
	  tgae->GetXaxis()->SetTitle("Muon |#it{#eta}|"); 
	  if(k>0) tgae2->GetXaxis()->SetTitle("Muon |#it{#eta}|"); 
	} 
	else if(var[w].EqualTo("phi")) { 
	  tgae->GetXaxis()->SetTitle("Muon #it{#phi} [rad]"); 

	  // Compute integrated efficiency 
	  float central       = histos_num[j][k][w]->Integral()/histos_den[j][k][w]->Integral(); 
	  float minusonesigma = central - TEfficiency::ClopperPearson(histos_den[j][k][w]->Integral(), histos_num[j][k][w]->Integral(), 0.683, 0); 
	  float plusonesigma  = TEfficiency::ClopperPearson(histos_den[j][k][w]->Integral(), histos_num[j][k][w]->Integral(), 0.683, 1) - central; 

	  // Compute integrated efficiency 
	  Double_t bkg_ctr(0.), bkg_err(0.); 
	  if(k>0) { 
	    tgae2->GetXaxis()->SetTitle("Muon #it{#phi} [rad]"); 
	    bkg_ctr = h_fake_avg_perevt[j][k-1][w]->IntegralAndError(0, -1, bkg_err) / h_evt_seg->GetBinContent(1); 
	    bkg_err /= h_evt_seg->GetBinContent(1); 
	  } 

	  if(type[k].EqualTo("eff")) { 
	    //std::cout << "  *** Efficiency = "; 
	    xxx[j] = central; 
	    exl[j] = minusonesigma; 
	    exh[j] = plusonesigma; 
	  } 
	  else if(type[k].EqualTo("fake")) { 
	    //std::cout << "  *** Fake rate = "; 
	    //// Fake rate 
	    // yyy[j] = central; 
	    // eyl[j] = minusonesigma; 
	    // eyh[j] = plusonesigma; 
	    //// Fakes per event 
	    yyy[j] = bkg_ctr; 
	    eyl[j] = bkg_err; 
	    eyh[j] = bkg_err; 
	  } 
	  else if(type[k].EqualTo("bkg")) { 
	    //std::cout << "  *** Background rate = "; 
	    //// Background rate 
	    // zzz[j] = central; 
	    // ezl[j] = minusonesigma; 
	    // ezh[j] = plusonesigma; 
	    //// Background per event 
	    zzz[j] = bkg_ctr; 
	    ezl[j] = bkg_err; 
	    ezh[j] = bkg_err; 
	  } 
	  else {} 

	  //std::cout << central << " - " << minusonesigma << " + " << plusonesigma << std::endl; 
	} 
	else if(var[w].EqualTo("allpu")) { 
	  // tgae->GetXaxis()->SetTitle("Simulated vertices"); 
	  // if(k>0) tgae2->GetXaxis()->SetTitle("Simulated vertices"); 
	  tgae->GetXaxis()->SetTitle("Average"); 
	  if(k>0) tgae2->GetXaxis()->SetTitle("Average"); 
	} 
	else {} 

	tgae->GetXaxis()->SetTitleOffset(0.84); 
	tgae->GetYaxis()->SetTitleOffset(0.60); 
	if(k>0) { 
	  tgae2->GetXaxis()->SetTitleOffset(0.84); 
	  tgae2->GetYaxis()->SetTitleOffset(0.60); 
	} 

	if(type[k].EqualTo("eff")) { 
	  tgae->GetYaxis()->SetTitle("Efficiency"); 
	} 
	else if(type[k].EqualTo("fake")) { 
	  tgae->GetYaxis()->SetTitle("Fake rate"); 
	  tgae2->GetYaxis()->SetTitle("Fakes/event"); 
	} 
	else if(type[k].EqualTo("bkg")) { 
	  tgae->GetYaxis()->SetTitle("Background rate"); 
	  tgae2->GetYaxis()->SetTitle("Background/event"); 
	} 
	else {} 

	cc->cd(); 
	cc->SaveAs( (outpath+"/"+outdirs[j]+"/"+finname+".png").Data() ); 
	cc->SetLogx(0); 
	if(k>0) { 
	  cc2->cd(); 
	  cc2->SaveAs( (outpath+"/"+outdirs[j]+"/"+finname+"_perevt.png").Data() ); 
	  cc2->SetLogx(0); 
	} 

	//delete tgae; 
      } // End loop over folders (i.e. selection parameters) 
      } 

      TString xytitle = ""; 
      TString xytitle2 = ""; 
      TString xytitle3 = ";Segment |#it{#eta}|;Muons/segment"; 
      TString xytitle4 = ";Track-segment #Delta#phi;Muons/segment"; 

      if(var[w].EqualTo("pt")) { 
	xytitle  += ";Muon #it{p}_{T} [GeV]"; 
	xytitle2 += ";Muon #it{p}_{T} [GeV]"; 
      } 
      else if(var[w].EqualTo("eta")) { 
	xytitle  += ";Muon #it{#eta}"; 
	xytitle2 += ";Muon #it{#eta}"; 
      } 
      else if(var[w].EqualTo("aeta")) { 
	xytitle  += ";Muon |#it{#eta}|"; 
	xytitle2 += ";Muon |#it{#eta}|"; 
      } 
      else if(var[w].EqualTo("phi")) { 
	xytitle  += ";Muon #it{#phi} [rad]"; 
	xytitle2 += ";Muon #it{#phi} [rad]"; 
      }
      else if(var[w].EqualTo("allpu")) { 
	// xytitle  += ";Simulated vertices"; 
	// xytitle2 += ";Simulated vertices"; 
	xytitle  += ";Average"; 
	xytitle2 += ";Average"; 
      }

      if(type[k].EqualTo("eff")) { 
	xytitle += ";Efficiency"; 
      } 
      else if(type[k].EqualTo("fake")) { 
	xytitle  += ";Fake rate"; 
	xytitle2 += ";Fakes/event"; 
      } 
      else if(type[k].EqualTo("bkg")) { 
	xytitle  += ";Background rate"; 
	xytitle2 += ";Background/event"; 
      } 
      else {} 


      TString label; 

      // All 
      // mgall->Add(tgaes[0], "lep"); label.Form("#Deltax/#sigma_{x} < %.1f", params[0][0]); legall->AddEntry(tgaes[0], label.Data(), "LEP"); 
      // mgall->Add(tgaes[1], "lep"); label.Form("#Deltax/#sigma_{x} < %.1f", params[1][0]); legall->AddEntry(tgaes[1], label.Data(), "LEP"); 
      // mgall->Add(tgaes[2], "lep"); label.Form("#Deltax/#sigma_{x} < %.1f", params[2][0]); legall->AddEntry(tgaes[2], label.Data(), "LEP"); 
      // mgall->Add(tgaes[3], "lep"); label.Form("#Deltax/#sigma_{x} < %.1f", params[3][0]); legall->AddEntry(tgaes[3], label.Data(), "LEP"); 
      // mgall->Add(tgaes[4], "lep"); label.Form("#Deltax/#sigma_{x} < %.1f", params[4][0]); legall->AddEntry(tgaes[4], label.Data(), "LEP"); 
      // mgall->Add(tgaes[5], "lep"); label.Form("#Deltax/#sigma_{x} < %.1f", params[5][0]); legall->AddEntry(tgaes[5], label.Data(), "LEP"); 
      for(size_t j=0; j<nrows; ++j) { // Loop over folders (i.e. selection paramters) 
	mgall->Add(tgaes[j], "lep"); 
	label.Form("%d strips, %d rolls", strstriproll[j].first, strstriproll[j].second); 
	legall->AddEntry(tgaes[j], label.Data(), "LEP"); 
      } 
      ccall->cd(); 
      //ccall->SetLogy(); 
      ccall->SetGridx(); 
      ccall->SetGridy(); 
      mgall->SetTitle(xytitle.Data()); 
      mgall->Draw("a"); 
      if(var[w].EqualTo("pt")) { 
	if(k==0) 
	  mgall->GetXaxis()->SetLimits(0.0,  35.0); 
	else { 
	  ccall->SetLogx(); 
	  mgall->GetXaxis()->SetLimits(0.3, 200.0); 
	} 
      } 
      mgall->GetXaxis()->SetTitleOffset(0.84); 
      mgall->GetYaxis()->SetTitleOffset(0.60); 
      legall->Draw(); 
      ccall->SaveAs( (outpath+"/"  +finname+".png").Data() ); 
      ccall->SaveAs( (outpath+"/"  +finname+".root").Data() ); 
      ccall->SetLogx(0); 
      ccall->SetGridx(0); 
      ccall->SetGridy(0); 

      if(k>0) { 
	// mgall2->Add(tgaes2[0], "lep"); 
	// mgall2->Add(tgaes2[1], "lep"); 
	// mgall2->Add(tgaes2[2], "lep"); 
	// mgall2->Add(tgaes2[3], "lep"); 
	// mgall2->Add(tgaes2[4], "lep"); 
	// mgall2->Add(tgaes2[5], "lep"); 
	for(size_t j=0; j<nrows; ++j) { // Loop over folders (i.e. selection paramters) 
	  mgall2->Add(tgaes2[j], "lep"); 
	} 
	ccall->cd(); 
	ccall->Clear(); 
	ccall->SetGridx(); 
	ccall->SetGridy(); 
	mgall2->SetTitle(xytitle2.Data()); 
	mgall2->Draw("a"); 
	// mgall2->GetXaxis()->SetTitleOffset(0.84); 
	// mgall2->GetYaxis()->SetTitleOffset(0.60); 
	legall->Draw(); 
	//ccall->SetLogy(); 
	if(var[w].EqualTo("pt")) { 
	  mgall2->GetYaxis()->SetRangeUser(10.0, 400.0); 
	  ccall->SetLogy(); 
	  mgall2->GetYaxis()->UnZoom(); 
	  ccall->Modified(); 
	  ccall->Update(); 
	  ccall->SetLogx(); 
	  mgall2->GetXaxis()->SetLimits(0.3, 200.0); 
	} 
	mgall2->GetXaxis()->SetTitleOffset(0.84); 
	mgall2->GetYaxis()->SetTitleOffset(0.60); 
	ccall->SaveAs( (outpath+"/"  +finname+"_perevt.png").Data() ); 
	ccall->SaveAs( (outpath+"/"  +finname+"_perevt.root").Data() ); 
	ccall->SetLogx(0); 
	ccall->SetLogy(0); 
	ccall->SetGridx(0); 
	ccall->SetGridy(0); 
      } 

      if(k==0 && w==0) { 
	// mgall3->Add(tgaes3[0], "lep"); 
	// mgall3->Add(tgaes3[1], "lep"); 
	// mgall3->Add(tgaes3[2], "lep"); 
	// mgall3->Add(tgaes3[3], "lep"); 
	// mgall3->Add(tgaes3[4], "lep"); 
	// mgall3->Add(tgaes3[5], "lep"); 
	for(size_t j=0; j<nrows; ++j) { // Loop over folders (i.e. selection paramters) 
	  mgall3->Add(tgaes3[j], "lep"); 
	} 
	ccall->cd(); 
	ccall->Clear(); 
	mgall3->SetTitle(xytitle3.Data()); 
	mgall3->Draw("a"); 
	mgall3->GetXaxis()->SetTitleOffset(0.84); 
	mgall3->GetYaxis()->SetTitleOffset(0.60); 
	legall->Draw(); 
	//ccall->SetLogy(); 
	ccall->SetGridx(); 
	ccall->SetGridy(); 
	ccall->SaveAs( (outpath+"/nmu_perseg_aeta.png").Data() ); 
	ccall->SaveAs( (outpath+"/nmu_perseg_aeta.root").Data() ); 
	ccall->SetGridx(0); 
	ccall->SetGridy(0); 

	// mgall4->Add(tgaes4[0], "lep"); 
	// mgall4->Add(tgaes4[1], "lep"); 
	// mgall4->Add(tgaes4[2], "lep"); 
	// mgall4->Add(tgaes4[3], "lep"); 
	// mgall4->Add(tgaes4[4], "lep"); 
	// mgall4->Add(tgaes4[5], "lep"); 
	for(size_t j=0; j<nrows; ++j) { // Loop over folders (i.e. selection paramters) 
	  mgall4->Add(tgaes4[j], "lep"); 
	} 
	ccall->cd(); 
	ccall->Clear(); 
	mgall4->SetTitle(xytitle4.Data()); 
	mgall4->Draw("a"); 
	mgall4->GetXaxis()->SetTitleOffset(0.84); 
	mgall4->GetYaxis()->SetTitleOffset(0.60); 
	legall->Draw(); 
	//ccall->SetLogy(); 
	ccall->SetGridx(); 
	ccall->SetGridy(); 
	ccall->SaveAs( (outpath+"/nmu_perseg_dphi.png").Data() ); 
	ccall->SaveAs( (outpath+"/nmu_perseg_dphi.root").Data() ); 
	ccall->SetGridx(0); 
	ccall->SetGridy(0); 
      } 

    } // End loop over vars (eta, aeta, phi, pt) 
  } // End loop over types (efficiency, fake rate) 


  //// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= //// 
  //// ROCs 
  //// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= //// 
  // ROC efficiency vs fake rate 
  cc->cd(); 
  //TGraphAsymmErrors *tgroc = new TGraphAsymmErrors(nrows, xxx, yyy, exl, exh, eyl, eyh); 
  TGraphAsymmErrors *tgroc = new TGraphAsymmErrors(nrows, xxx, yyy, 0, 0, 0, 0); 
  tgroc->SetNameTitle("ROC", "ROC"); 
  tgroc->Draw("AP"); // Just to get the right size of the frame 
  //tgroc->SetMarkerStyle(20); 
  tgroc->SetMarkerSize(0); 
  tgroc->GetXaxis()->SetTitle("Efficiency"); 
  tgroc->GetXaxis()->SetTitleOffset(0.84); 
  //tgroc->GetYaxis()->SetTitle("Fake rate"); 
  tgroc->GetYaxis()->SetTitle("Avg. fakes/event"); 
  tgroc->GetYaxis()->SetTitleOffset(0.60); 
  //cc->SaveAs( (outpath+"/roc_eff_fake.png").Data() ); 

  //TLegend *leg = new TLegend(0.17, 0.60, 0.55, 0.95); 
  TLegend *leg = new TLegend(legxmin, legymin, legxmax, legymax); 
  leg->SetFillColor(kWhite);
  leg->SetFillStyle(1001);
  leg->SetLineStyle(0);
  leg->SetLineWidth(0);
  leg->SetBorderSize(0);
  //leg->SetTextFont(42);


  size_t idx1(0), idx2(0); 
  std::vector<TGraphAsymmErrors*> tgalls; 
  std::map<Int_t, std::vector<Int_t> >::const_iterator itsr = striproll.begin(); 
  for(; itsr!=striproll.end(); ++itsr) { 
    tgalls.push_back( new TGraphAsymmErrors(itsr->second.size()) ); 
    TGraphAsymmErrors *tgall = tgalls.back(); 
    for(size_t jdx=0; jdx<itsr->second.size(); ++jdx, ++idx1) { 
      tgall->SetPoint(jdx, xxx[idx1], yyy[idx1]); tgall->SetPointError(jdx, exl[idx1], exh[idx1], eyl[idx1], eyh[idx1]); 
    } 
    tgall->SetLineColor(cols[idx2%5]); 
    tgall->SetLineWidth(2); 
    tgall->SetMarkerStyle(stls[(idx2/5)%4]); 
    tgall->SetLineStyle(stll[(idx2/5)%4]); 
    tgall->SetMarkerColor(cols[idx2%5]); idx2++; 
    tgall->SetMarkerSize(1.6); 
    tgall->Draw("PL"); 
    char buff[99]; 
    sprintf(buff, "%d strips", itsr->first); 
    leg->AddEntry(tgall, buff, "PL"); 
  } 

  leg->Draw("same"); 
  //cc->SetLogy(); 
  cc->SetGridx(); 
  cc->SetGridy(); 
  cc->SaveAs( (outpath+"/roc_eff_fake.png").Data() ); 
  cc->SaveAs( (outpath+"/roc_eff_fake.root").Data() ); 
  cc->SetGridx(0); 
  cc->SetGridy(0); 


  //// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= //// 
  // ROC efficiency vs background rate 
  cc->Clear(); cc->cd(); 
  //TGraphAsymmErrors *tgroc = new TGraphAsymmErrors(nrows, xxx, zzz, exl, exh, ezl, ezh); 
  tgroc = new TGraphAsymmErrors(nrows, xxx, zzz, 0, 0, 0, 0); 
  tgroc->SetNameTitle("ROC2", "ROC2"); 
  tgroc->Draw("AP"); // Just to get the right size of the frame 
  //tgroc->SetMarkerStyle(20); 
  tgroc->SetMarkerSize(0); 
  tgroc->GetXaxis()->SetTitle("Efficiency"); 
  tgroc->GetXaxis()->SetTitleOffset(0.84); 
  //tgroc->GetYaxis()->SetTitle("Background rate"); 
  tgroc->GetYaxis()->SetTitle("Avg. background/event"); 
  tgroc->GetYaxis()->SetTitleOffset(0.60); 
  //cc->SaveAs( (outpath+"/roc_eff_fake.png").Data() ); 

  //leg = new TLegend(0.17, 0.60, 0.55, 0.95); 
  leg = new TLegend(legxmin, legymin, legxmax, legymax); 
  leg->SetFillColor(kWhite);
  leg->SetFillStyle(1001);
  leg->SetLineStyle(0);
  leg->SetLineWidth(0);
  leg->SetBorderSize(0);
  //leg->SetTextFont(42);


  idx1 = 0; 
  idx2 = 0; 
  itsr = striproll.begin(); 
  for(; itsr!=striproll.end(); ++itsr) { 
    tgalls.push_back( new TGraphAsymmErrors(itsr->second.size()) ); 
    TGraphAsymmErrors *tgall = tgalls.back(); 
    for(size_t jdx=0; jdx<itsr->second.size(); ++jdx, ++idx1) { 
      tgall->SetPoint(jdx, xxx[idx1], zzz[idx1]); tgall->SetPointError(jdx, exl[idx1], exh[idx1], ezl[idx1], ezh[idx1]); 
    } 
    tgall->SetLineColor(cols[idx2%5]); 
    tgall->SetLineWidth(2); 
    tgall->SetMarkerStyle(stls[(idx2/5)%4]); 
    tgall->SetLineStyle(stll[(idx2/5)%4]); 
    tgall->SetMarkerColor(cols[idx2%5]); idx2++; 
    tgall->SetMarkerSize(1.6); 
    tgall->Draw("PL"); 
    char buff[99]; 
    sprintf(buff, "%d strips", itsr->first); 
    leg->AddEntry(tgall, buff, "PL"); 
  } 

  leg->Draw("same"); 
  // tgall->GetXaxis()->SetTitle("Efficiency"); 
  // tgall->GetXaxis()->SetTitleOffset(0.84); 
  // tgall->GetYaxis()->SetTitle("Fake rate"); 
  // tgall->GetYaxis()->SetTitleOffset(1.10); 
  //cc->SetLogy(); 
  cc->SetGridx(); 
  cc->SetGridy(); 
  cc->SaveAs( (outpath+"/roc_eff_bkg.png").Data() ); 
  cc->SaveAs( (outpath+"/roc_eff_bkg.root").Data() ); 
  cc->SetGridx(0); 
  cc->SetGridy(0); 

  return; 
} 

