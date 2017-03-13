#include "TDatime.h" 
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

void mergeRootFiles_optim(TString folder = "") {

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
    // //// Strip: 768; Eta partitions: 4, 6, 8, 16 
    // "s768p004", 
    // "s768p006", 
    // "s768p008", 
    // "s768p012", 
    // "s768p016", 
    // //"s768p032", 
    // //"s768p064", 
    // //// Strip: 512; Eta partitions: 4, 6, 8, 16 
    // "s512p004", 
    "s512p006" //, 
    // "s512p008", 
    // "s512p012", 
    // "s512p016", 
    // //"s512p032", 
    // //"s512p064", 
    // //// Strip: 384; Eta partitions: 4, 6, 8, 16 
    // "s384p004", 
    // "s384p006", 
    // "s384p008", 
    // "s384p012", 
    // "s384p016", 
    // //"s384p032", 
    // //"s384p064", 
    // //// Strip: 256; Eta partitions: 4, 6, 8, 16 
    // "s256p004", 
    // "s256p006", 
    // "s256p008", 
    // "s256p012", 
    // "s256p016" //, 
    // //"s256p032", 
    // //"s256p064" 
  }; 

  if(sizeof(segmentcoll)==0) { 
    std::cout << " *** ERROR: segmentcoll can't be empty! Add {""}! *** " << std::endl; 
    return; 
  } 
  // static const size_t nrows = sizeof(segmentcoll)/sizeof(segmentcoll[0]); 

  double params[][5] = {
    {-1.0,  4.0, -1.0, -1.0, -1.0}, //  0 
    { 3.0,  4.0, -1.0, -1.0, -1.0}, //  1 (default) 
    { 2.0,  4.0, -1.0, -1.0, -1.0}, //  2 
    { 1.0,  4.0, -1.0, -1.0, -1.0}, //  3 
    { 0.5,  4.0, -1.0, -1.0, -1.0}, //  4 
    { 0.1,  4.0, -1.0, -1.0, -1.0}, //  5 

    { 3.0, -1.0, -1.0, -1.0, -1.0}, //  6 
    { 3.0,  3.0, -1.0, -1.0, -1.0}, //  7 
    { 3.0,  2.0, -1.0, -1.0, -1.0}, //  8 
    { 3.0,  1.0, -1.0, -1.0, -1.0}, //  9 
    { 3.0,  0.5, -1.0, -1.0, -1.0}, // 10 

    { 3.0,  4.0, 10.0, -1.0, -1.0}, // 11 
    { 3.0,  4.0,  5.0, -1.0, -1.0}, // 12 
    { 3.0,  4.0,  4.0, -1.0, -1.0}, // 13 
    { 3.0,  4.0,  3.0, -1.0, -1.0}, // 14 
    { 3.0,  4.0,  2.0, -1.0, -1.0}, // 15 
    { 3.0,  4.0,  1.0, -1.0, -1.0}, // 16 
    { 3.0,  4.0,  0.5, -1.0, -1.0}, // 17 

    { 3.0,  4.0, -1.0, 20.0, -1.0}, // 18 
    { 3.0,  4.0, -1.0, 10.0, -1.0}, // 19 
    { 3.0,  4.0, -1.0,  5.0, -1.0}, // 20 
    { 3.0,  4.0, -1.0,  4.0, -1.0}, // 21 
    { 3.0,  4.0, -1.0,  3.0, -1.0}, // 22 
    { 3.0,  4.0, -1.0,  2.0, -1.0}, // 23 
    { 3.0,  4.0, -1.0,  1.0, -1.0}, // 24 

    { 3.0,  4.0, -1.0, -1.0,  1.6}, // 25 
    { 3.0,  4.0, -1.0, -1.0,  0.8}, // 26 
    { 3.0,  4.0, -1.0, -1.0,  0.4}, // 27 
    { 3.0,  4.0, -1.0, -1.0,  0.2}, // 28 
    { 3.0,  4.0, -1.0, -1.0,  0.1}  // 29 
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
	  dirname.Form("PullX%.1fDiffX%.1fPullY%.1fDiffY%.1fDiffPhi%.1f", params[j1][0], params[j1][1], params[j1][2], params[j1][3], params[j1][4]); 
	  //dirname.Form("%.1fDiffX%.1fPullY%.1fDiffY%.1fDiffPhi%.1f", params[j1][0], params[j1][1], params[j1][2], params[j1][3], params[j1][4]); 
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


  TCanvas *cpullx   = new TCanvas("cpullx"  , "cpullx"  , 1000, 700); 
  TCanvas *cdiffx   = new TCanvas("cdiffx"  , "cdiffx"  , 1000, 700); 
  TCanvas *cpully   = new TCanvas("cpully"  , "cpully"  , 1000, 700); 
  TCanvas *cdiffy   = new TCanvas("cdiffy"  , "cdiffy"  , 1000, 700); 
  TCanvas *cdiffphi = new TCanvas("cdiffphi", "cdiffphi", 1000, 700); 

  float legxmin(0.75), legxmax(0.95), legymin(0.14), legymax(0.34); 


  float xxx[nrows], yyy[nrows], zzz[nrows], exl[nrows], exh[nrows], eyl[nrows], eyh[nrows], ezl[nrows], ezh[nrows]; 
  Color_t cols[9] = {kBlack, kBlue, kViolet, kMagenta-9, kRed, kOrange+1, kYellow+2, kGreen+2, kCyan+2}; 
  //Style_t stls[9] = {20, 21, 22, 23, 24, 25, 26, 27, 28}; 
  Style_t stls[4] = {20, 24, 33, 27}; 
  Style_t stll[4] = {1, 7, 2, 3}; 

  size_t ncols = sizeof(cols)/sizeof(cols[0]); 
  size_t nstls = sizeof(stls)/sizeof(stls[0]); 
  size_t nstll = sizeof(stll)/sizeof(stll[0]); 

  for(size_t k=0; k<ntypes; ++k) { // Loop over types (efficiency, fake rate, background) 
    for(size_t w=0; w<nvars; ++w) { // Loop over vars (eta, aeta, phi, pt) 
      TString finname = (type[k]+"_"+var[w]); 

      cpullx  ->Clear(); 
      cdiffx  ->Clear(); 
      cpully  ->Clear(); 
      cdiffy  ->Clear(); 
      cdiffphi->Clear(); 

      // k: type[] = {"eff", "fake", "bkg"} 
      // w:  var[] = {"eta", "aeta", "phi", "pt", "allpu"} 

      /* 
      if( (w==0 && k==0) || (w==0 && k==1) ) { 
	legxmin = 0.45; 
	legxmax = 0.65; 
	legymin = 0.14; 
	legymax = 0.34; 
      } 
      else if( (w==1 && k==0) || (w==3 && k==1) ) { 
	legxmin = 0.17; 
	legxmax = 0.37; 
	legymin = 0.71; 
	legymax = 0.91; 
      } 
      else if( (w==1 && k==1) ) { 
	legxmin = 0.75; 
	legxmax = 0.95; 
	legymin = 0.71; 
	legymax = 0.91; 
      } 
      else { 
	legxmin = 0.75; 
	legxmax = 0.95; 
	legymin = 0.14; 
	legymax = 0.34; 
      } 
      */ 

      legxmin = 0.72; 
      legxmax = 1.00; 
      legymin = 0.12; 
      legymax = 0.95; 

      TLegend *legpullx = new TLegend(legxmin, legymin, legxmax, legymax); 
      legpullx->SetFillColor(kWhite);
      legpullx->SetFillStyle(1001);
      legpullx->SetLineStyle(0);
      legpullx->SetLineWidth(0);
      legpullx->SetBorderSize(0);
      //legpullx->SetTextFont(42);

      TLegend *legdiffx = new TLegend(legxmin, legymin, legxmax, legymax); 
      legdiffx->SetFillColor(kWhite);
      legdiffx->SetFillStyle(1001);
      legdiffx->SetLineStyle(0);
      legdiffx->SetLineWidth(0);
      legdiffx->SetBorderSize(0);
      //legdiffx->SetTextFont(42);

      TLegend *legpully = new TLegend(legxmin, legymin, legxmax, legymax); 
      legpully->SetFillColor(kWhite);
      legpully->SetFillStyle(1001);
      legpully->SetLineStyle(0);
      legpully->SetLineWidth(0);
      legpully->SetBorderSize(0);
      //legpully->SetTextFont(42);

      TLegend *legdiffy = new TLegend(legxmin, legymin, legxmax, legymax); 
      legdiffy->SetFillColor(kWhite);
      legdiffy->SetFillStyle(1001);
      legdiffy->SetLineStyle(0);
      legdiffy->SetLineWidth(0);
      legdiffy->SetBorderSize(0);
      //legdiffy->SetTextFont(42);

      TLegend *legdiffphi = new TLegend(legxmin, legymin, legxmax, legymax); 
      legdiffphi->SetFillColor(kWhite);
      legdiffphi->SetFillStyle(1001);
      legdiffphi->SetLineStyle(0);
      legdiffphi->SetLineWidth(0);
      legdiffphi->SetBorderSize(0);
      //legdiffphi->SetTextFont(42);

      std::vector<TGraphAsymmErrors*> tgaes, tgaes2, tgaes3, tgaes4; 

      // MultiGraphs 
      TMultiGraph *mgpullx   = new TMultiGraph(); 
      TMultiGraph *mgdiffx   = new TMultiGraph(); 
      TMultiGraph *mgpully   = new TMultiGraph(); 
      TMultiGraph *mgdiffy   = new TMultiGraph(); 
      TMultiGraph *mgdiffphi = new TMultiGraph(); 

      TMultiGraph *mgpullx2   = new TMultiGraph(); 
      TMultiGraph *mgdiffx2   = new TMultiGraph(); 
      TMultiGraph *mgpully2   = new TMultiGraph(); 
      TMultiGraph *mgdiffy2   = new TMultiGraph(); 
      TMultiGraph *mgdiffphi2 = new TMultiGraph(); 

      TMultiGraph *mgpullx3   = new TMultiGraph(); 
      TMultiGraph *mgdiffx3   = new TMultiGraph(); 
      TMultiGraph *mgpully3   = new TMultiGraph(); 
      TMultiGraph *mgdiffy3   = new TMultiGraph(); 
      TMultiGraph *mgdiffphi3 = new TMultiGraph(); 

      TMultiGraph *mgpullx4   = new TMultiGraph(); 
      TMultiGraph *mgdiffx4   = new TMultiGraph(); 
      TMultiGraph *mgpully4   = new TMultiGraph(); 
      TMultiGraph *mgdiffy4   = new TMultiGraph(); 
      TMultiGraph *mgdiffphi4 = new TMultiGraph(); 


      //for(size_t j=0; j<nrows; ++j) { // Loop over folders (i.e. selection paramters) 
      for(size_t j0=0; j0<nsegmts; ++j0) { // Loop over segment collections  
      for(size_t j1=0; j1<nparams; ++j1) { // Loop over selection paramters 
	size_t j = j0*nparams + j1; 
	// Let's start with the number of muons per segment 
	if(k==0 && w==0) { // Only once... 
	  h_mu_perseg_aeta[j]->Divide(h_mu_perseg_den_aeta[j]); 
	  h_mu_perseg_aeta[j]->SetNameTitle( (TString("nmu_avg_perseg_aeta_")+outdirs[j]).Data(), (TString("nmu_avg_perseg_aeta_")+outdirs[j]).Data() ); 
	  cc->cd(); 
	  tgaes3.push_back(new TGraphAsymmErrors(h_mu_perseg_aeta[j])); 
	  TGraphAsymmErrors *tgae3 = tgaes3.back(); 
	  tgae3->SetLineWidth(2); 
	  tgae3->SetLineColor(cols[j%ncols]); 
	  tgae3->SetLineStyle(stll[j%nstll]); 
	  tgae3->SetMarkerSize(1.0); 
	  tgae3->SetMarkerColor(cols[j%ncols]); 
	  tgae3->SetMarkerStyle(stls[j%nstls]); 
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
	  tgae4->SetLineColor(cols[j%ncols]); 
	  tgae4->SetLineStyle(stll[j%nstll]); 
	  tgae4->SetMarkerSize(1.0); 
	  tgae4->SetMarkerColor(cols[j%ncols]); 
	  tgae4->SetMarkerStyle(stls[j%nstls]); 
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
	TGraphAsymmErrors *tgae = tgaes.back(); 
	tgae->SetNameTitle( (finname+"_"+outdirs[j]).Data(), (finname+"_"+outdirs[j]).Data()); 
	tgae->Divide(histos_num[j][k][w], histos_den[j][k][w], "cl=0.683 b(1,1) mode" ); 
	tgae->SetLineWidth(2); 
	tgae->SetLineColor(cols[j%ncols]);   if(j==5) tgae->SetLineColor(cols[0]); 
	tgae->SetLineStyle(stll[j%nstll]); 
	tgae->SetMarkerSize(1.0); 
	tgae->SetMarkerColor(cols[j%ncols]); if(j==5) tgae->SetMarkerColor(cols[0]); 
	tgae->SetMarkerStyle(stls[j%nstls]); if(j==5) tgae->SetMarkerStyle(stls[0]); 
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
	  tgae2->SetLineColor(cols[j%ncols]); 
	  tgae2->SetLineStyle(stls[j%nstll]); 
	  tgae2->SetMarkerSize(1.0); 
	  tgae2->SetMarkerColor(cols[j%ncols]); 
	  tgae2->SetMarkerStyle(stls[j%nstls]); 
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

      // Pull x: 0, 1, 2, 3, 4, 5  
      mgpullx->Add(tgaes[0], "lep"); label.Form("#Deltax/#sigma_{x} < %.1f", params[0][0]); legpullx->AddEntry(tgaes[0], label.Data(), "LEP"); 
      mgpullx->Add(tgaes[1], "lep"); label.Form("#Deltax/#sigma_{x} < %.1f", params[1][0]); legpullx->AddEntry(tgaes[1], label.Data(), "LEP"); 
      mgpullx->Add(tgaes[2], "lep"); label.Form("#Deltax/#sigma_{x} < %.1f", params[2][0]); legpullx->AddEntry(tgaes[2], label.Data(), "LEP"); 
      mgpullx->Add(tgaes[3], "lep"); label.Form("#Deltax/#sigma_{x} < %.1f", params[3][0]); legpullx->AddEntry(tgaes[3], label.Data(), "LEP"); 
      mgpullx->Add(tgaes[4], "lep"); label.Form("#Deltax/#sigma_{x} < %.1f", params[4][0]); legpullx->AddEntry(tgaes[4], label.Data(), "LEP"); 
      mgpullx->Add(tgaes[5], "lep"); label.Form("#Deltax/#sigma_{x} < %.1f", params[5][0]); legpullx->AddEntry(tgaes[5], label.Data(), "LEP"); 
      cpullx->cd(); 
      //cpullx->SetLogy(); 
      cpullx->SetGridx(); 
      cpullx->SetGridy(); 
      mgpullx->SetTitle(xytitle.Data()); 
      mgpullx->Draw("a"); 
      if(var[w].EqualTo("pt")) { 
	if(k==0) 
	  mgpullx->GetXaxis()->SetLimits(0.0,  35.0); 
	else { 
	  cpullx->SetLogx(); 
	  mgpullx->GetXaxis()->SetLimits(0.3, 200.0); 
	} 
      } 
      mgpullx->GetXaxis()->SetTitleOffset(0.84); 
      mgpullx->GetYaxis()->SetTitleOffset(0.60); 
      legpullx->Draw(); 
      //cpullx->SetLogy(); 
      cpullx->SaveAs( (outpath+"/pullx_"  +finname+".png").Data() ); 
      cpullx->SaveAs( (outpath+"/pullx_"  +finname+".root").Data() ); 
      cpullx->SetLogx(0); 
      cpullx->SetGridx(0); 
      cpullx->SetGridy(0); 

      if(k>0) { 
	mgpullx2->Add(tgaes2[0], "lep"); 
	mgpullx2->Add(tgaes2[1], "lep"); 
	mgpullx2->Add(tgaes2[2], "lep"); 
	mgpullx2->Add(tgaes2[3], "lep"); 
	mgpullx2->Add(tgaes2[4], "lep"); 
	mgpullx2->Add(tgaes2[5], "lep"); 
	cpullx->cd(); 
	cpullx->Clear(); 
	//cpullx->SetLogy(); 
	cpullx->SetGridx(); 
	cpullx->SetGridy(); 
	mgpullx2->SetTitle(xytitle2.Data()); 
	mgpullx2->Draw("a"); 
	if(var[w].EqualTo("pt")) { 
	  mgpullx2->GetYaxis()->SetRangeUser(10.0, 400.0); 
	  cpullx->SetLogy(); 
	  mgpullx2->GetYaxis()->UnZoom(); 
	  cpullx->Modified(); 
	  cpullx->Update(); 
	  cpullx->SetLogx(); 
	  mgpullx2->GetXaxis()->SetLimits(0.3, 200.0); 
	} 
	mgpullx2->GetXaxis()->SetTitleOffset(0.84); 
	mgpullx2->GetYaxis()->SetTitleOffset(0.60); 
	legpullx->Draw(); 
	//cpullx->SetLogy(); 
	cpullx->SaveAs( (outpath+"/pullx_"  +finname+"_perevt.png").Data() ); 
	cpullx->SaveAs( (outpath+"/pullx_"  +finname+"_perevt.root").Data() ); 
	cpullx->SetLogx(0); 
	cpullx->SetLogy(0); 
	cpullx->SetGridx(0); 
	cpullx->SetGridy(0); 
      } 

      if(k==0 && w==0) { 
	mgpullx3->Add(tgaes3[0], "lep"); 
	mgpullx3->Add(tgaes3[1], "lep"); 
	mgpullx3->Add(tgaes3[2], "lep"); 
	mgpullx3->Add(tgaes3[3], "lep"); 
	mgpullx3->Add(tgaes3[4], "lep"); 
	mgpullx3->Add(tgaes3[5], "lep"); 
	cpullx->cd(); 
	cpullx->Clear(); 
	mgpullx3->SetTitle(xytitle3.Data()); 
	mgpullx3->Draw("a"); 
	mgpullx3->GetXaxis()->SetTitleOffset(0.84); 
	mgpullx3->GetYaxis()->SetTitleOffset(0.60); 
	legpullx->Draw(); 
	//cpullx->SetLogy(); 
	cpullx->SaveAs( (outpath+"/pullx_nmu_perseg_aeta.png").Data() ); 
	cpullx->SaveAs( (outpath+"/pullx_nmu_perseg_aeta.root").Data() ); 
	cpullx->SetGridx(0); 
	cpullx->SetGridy(0); 

	mgpullx4->Add(tgaes4[0], "lep"); 
	mgpullx4->Add(tgaes4[1], "lep"); 
	mgpullx4->Add(tgaes4[2], "lep"); 
	mgpullx4->Add(tgaes4[3], "lep"); 
	mgpullx4->Add(tgaes4[4], "lep"); 
	mgpullx4->Add(tgaes4[5], "lep"); 
	cpullx->cd(); 
	cpullx->Clear(); 
	mgpullx4->SetTitle(xytitle4.Data()); 
	mgpullx4->Draw("a"); 
	mgpullx4->GetXaxis()->SetTitleOffset(0.84); 
	mgpullx4->GetYaxis()->SetTitleOffset(0.60); 
	legpullx->Draw(); 
	//cpullx->SetLogy(); 
	cpullx->SetGridx(); 
	cpullx->SetGridy(); 
	cpullx->SaveAs( (outpath+"/pullx_nmu_perseg_dphi.png").Data() ); 
	cpullx->SaveAs( (outpath+"/pullx_nmu_perseg_dphi.root").Data() ); 
	cpullx->SetGridx(0); 
	cpullx->SetGridy(0); 
      } 

      // Diff. x: 6, 1, 7, 8, 9, 10 
      mgdiffx->Add(tgaes[ 6], "lep"); label.Form("#Deltax < %.1f", params[ 6][1]); legdiffx->AddEntry(tgaes[ 6], label.Data(), "LEP"); 
      mgdiffx->Add(tgaes[ 1], "lep"); label.Form("#Deltax < %.1f", params[ 1][1]); legdiffx->AddEntry(tgaes[ 1], label.Data(), "LEP"); 
      mgdiffx->Add(tgaes[ 7], "lep"); label.Form("#Deltax < %.1f", params[ 7][1]); legdiffx->AddEntry(tgaes[ 7], label.Data(), "LEP"); 
      mgdiffx->Add(tgaes[ 8], "lep"); label.Form("#Deltax < %.1f", params[ 8][1]); legdiffx->AddEntry(tgaes[ 8], label.Data(), "LEP"); 
      mgdiffx->Add(tgaes[ 9], "lep"); label.Form("#Deltax < %.1f", params[ 9][1]); legdiffx->AddEntry(tgaes[ 9], label.Data(), "LEP"); 
      mgdiffx->Add(tgaes[10], "lep"); label.Form("#Deltax < %.1f", params[10][1]); legdiffx->AddEntry(tgaes[10], label.Data(), "LEP"); 
      cdiffx->cd(); 
      //cdiffx->SetLogy(); 
      cdiffx->SetGridx(); 
      cdiffx->SetGridy(); 
      mgdiffx->SetTitle(xytitle.Data()); 
      mgdiffx->Draw("a"); 
      if(var[w].EqualTo("pt")) { 
	if(k==0) 
	  mgdiffx->GetXaxis()->SetLimits(0.0,  35.0); 
	else { 
	  cdiffx->SetLogx(); 
	  mgdiffx->GetXaxis()->SetLimits(0.3, 200.0); 
	} 
      } 
      mgdiffx->GetXaxis()->SetTitleOffset(0.84); 
      mgdiffx->GetYaxis()->SetTitleOffset(0.60); 
      legdiffx->Draw(); 
      //cdiffx->SetLogy(); 
      cdiffx->SaveAs( (outpath+"/diffx_"  +finname+".png").Data() ); 
      cdiffx->SaveAs( (outpath+"/diffx_"  +finname+".root").Data() ); 
      cdiffx->SetLogx(0); 
      cdiffx->SetGridx(0); 
      cdiffx->SetGridy(0); 

      if(k>0) { 
	mgdiffx2->Add(tgaes2[ 6], "lep"); 
	mgdiffx2->Add(tgaes2[ 1], "lep"); 
	mgdiffx2->Add(tgaes2[ 7], "lep"); 
	mgdiffx2->Add(tgaes2[ 8], "lep"); 
	mgdiffx2->Add(tgaes2[ 9], "lep"); 
	mgdiffx2->Add(tgaes2[10], "lep"); 
	cdiffx->cd(); 
	cdiffx->Clear(); 
	//cdiffx->SetLogy(); 
	cdiffx->SetGridx(); 
	cdiffx->SetGridy(); 
	mgdiffx2->SetTitle(xytitle2.Data()); 
	mgdiffx2->Draw("a"); 
	if(var[w].EqualTo("pt")) { 
	  mgdiffx2->GetYaxis()->SetRangeUser(10.0, 400.0); 
	  cdiffx->SetLogy(); 
	  mgdiffx2->GetYaxis()->UnZoom(); 
	  cdiffx->Modified(); 
	  cdiffx->Update(); 
	  cdiffx->SetLogx(); 
	  mgdiffx2->GetXaxis()->SetLimits(0.3, 200.0); 
	} 
	mgdiffx2->GetXaxis()->SetTitleOffset(0.84); 
	mgdiffx2->GetYaxis()->SetTitleOffset(0.60); 
	legdiffx->Draw(); 
	//cdiffx->SetLogy(); 
	cdiffx->SaveAs( (outpath+"/diffx_"  +finname+"_perevt.png").Data() ); 
	cdiffx->SaveAs( (outpath+"/diffx_"  +finname+"_perevt.root").Data() ); 
	cdiffx->SetLogx(0); 
	cdiffx->SetLogy(0); 
	cdiffx->SetGridx(0); 
	cdiffx->SetGridy(0); 
      } 

      if(k==0 && w==0) { 
	mgdiffx3->Add(tgaes3[ 6], "lep"); 
	mgdiffx3->Add(tgaes3[ 1], "lep"); 
	mgdiffx3->Add(tgaes3[ 7], "lep"); 
	mgdiffx3->Add(tgaes3[ 8], "lep"); 
	mgdiffx3->Add(tgaes3[ 9], "lep"); 
	mgdiffx3->Add(tgaes3[10], "lep"); 
	cdiffx->cd(); 
	cdiffx->Clear(); 
	mgdiffx3->SetTitle(xytitle3.Data()); 
	mgdiffx3->Draw("a"); 
	mgdiffx3->GetXaxis()->SetTitleOffset(0.84); 
	mgdiffx3->GetYaxis()->SetTitleOffset(0.60); 
	legdiffx->Draw(); 
	//cdiffx->SetLogy(); 
	cdiffx->SaveAs( (outpath+"/diffx_nmu_perseg_aeta.png").Data() ); 
	cdiffx->SaveAs( (outpath+"/diffx_nmu_perseg_aeta.root").Data() ); 
	cdiffx->SetGridx(0); 
	cdiffx->SetGridy(0); 

	mgdiffx4->Add(tgaes4[ 6], "lep"); 
	mgdiffx4->Add(tgaes4[ 1], "lep"); 
	mgdiffx4->Add(tgaes4[ 7], "lep"); 
	mgdiffx4->Add(tgaes4[ 8], "lep"); 
	mgdiffx4->Add(tgaes4[ 9], "lep"); 
	mgdiffx4->Add(tgaes4[10], "lep"); 
	cdiffx->cd(); 
	cdiffx->Clear(); 
	mgdiffx4->SetTitle(xytitle4.Data()); 
	mgdiffx4->Draw("a"); 
	mgdiffx4->GetXaxis()->SetTitleOffset(0.84); 
	mgdiffx4->GetYaxis()->SetTitleOffset(0.60); 
	legdiffx->Draw(); 
	//cdiffx->SetLogy(); 
	cdiffx->SetGridx(); 
	cdiffx->SetGridy(); 
	//cdiffx->SetLogy(); 
	cdiffx->SaveAs( (outpath+"/diffx_nmu_perseg_dphi.png").Data() ); 
	cdiffx->SaveAs( (outpath+"/diffx_nmu_perseg_dphi.root").Data() ); 
	cdiffx->SetGridx(0); 
	cdiffx->SetGridy(0); 
      } 

      // ====================================================================== 

      //
      // Pull y:  1, 11, 12, 13, 14, 15, 16, 17 
      mgpully->Add(tgaes[ 1], "lep"); label.Form("#Deltay/#sigma_{y} < %.1f", params[ 1][2]); legpully->AddEntry(tgaes[ 1], label.Data(), "LEP"); 
      mgpully->Add(tgaes[11], "lep"); label.Form("#Deltay/#sigma_{y} < %.1f", params[11][2]); legpully->AddEntry(tgaes[11], label.Data(), "LEP"); 
      mgpully->Add(tgaes[12], "lep"); label.Form("#Deltay/#sigma_{y} < %.1f", params[12][2]); legpully->AddEntry(tgaes[12], label.Data(), "LEP"); 
      mgpully->Add(tgaes[13], "lep"); label.Form("#Deltay/#sigma_{y} < %.1f", params[13][2]); legpully->AddEntry(tgaes[13], label.Data(), "LEP"); 
      mgpully->Add(tgaes[14], "lep"); label.Form("#Deltay/#sigma_{y} < %.1f", params[14][2]); legpully->AddEntry(tgaes[14], label.Data(), "LEP"); 
      mgpully->Add(tgaes[15], "lep"); label.Form("#Deltay/#sigma_{y} < %.1f", params[15][2]); legpully->AddEntry(tgaes[15], label.Data(), "LEP"); 
      mgpully->Add(tgaes[16], "lep"); label.Form("#Deltay/#sigma_{y} < %.1f", params[16][2]); legpully->AddEntry(tgaes[16], label.Data(), "LEP"); 
      mgpully->Add(tgaes[17], "lep"); label.Form("#Deltay/#sigma_{y} < %.1f", params[17][2]); legpully->AddEntry(tgaes[17], label.Data(), "LEP"); 
      cpully->cd(); 
      //cpully->SetLogy(); 
      cpully->SetGridx(); 
      cpully->SetGridy(); 
      mgpully->SetTitle(xytitle.Data()); 
      mgpully->Draw("a"); 
      if(var[w].EqualTo("pt")) { 
	if(k==0) 
	  mgpully->GetXaxis()->SetLimits(0.0,  35.0); 
	else { 
	  cpully->SetLogx(); 
	  mgpully->GetXaxis()->SetLimits(0.3, 200.0); 
	} 
      } 
      mgpully->GetXaxis()->SetTitleOffset(0.84); 
      mgpully->GetYaxis()->SetTitleOffset(0.60); 
      legpully->Draw(); 
      //cpully->SetLogy(); 
      cpully->SaveAs( (outpath+"/pully_"  +finname+".png").Data() ); 
      cpully->SaveAs( (outpath+"/pully_"  +finname+".root").Data() ); 
      cpully->SetLogx(0); 
      cpully->SetGridx(0); 
      cpully->SetGridy(0); 

      if(k>0) { 
	mgpully2->Add(tgaes2[ 1], "lep"); 
	mgpully2->Add(tgaes2[11], "lep"); 
	mgpully2->Add(tgaes2[12], "lep"); 
	mgpully2->Add(tgaes2[13], "lep"); 
	mgpully2->Add(tgaes2[14], "lep"); 
	mgpully2->Add(tgaes2[15], "lep"); 
	mgpully2->Add(tgaes2[16], "lep"); 
	mgpully2->Add(tgaes2[17], "lep"); 
	cpully->cd(); 
	cpully->Clear(); 
	//cpully->SetLogy(); 
	cpully->SetGridx(); 
	cpully->SetGridy(); 
	mgpully2->SetTitle(xytitle2.Data()); 
	mgpully2->Draw("a"); 
	if(var[w].EqualTo("pt")) { 
	  mgpully2->GetYaxis()->SetRangeUser(10.0, 400.0); 
	  cpully->SetLogy(); 
	  mgpully2->GetYaxis()->UnZoom(); 
	  cpully->Modified(); 
	  cpully->Update(); 
	  cpully->SetLogx(); 
	  mgpully2->GetXaxis()->SetLimits(0.3, 200.0); 
	} 
	mgpully2->GetXaxis()->SetTitleOffset(0.84); 
	mgpully2->GetYaxis()->SetTitleOffset(0.60); 
	legpully->Draw(); 
	//cpully->SetLogy(); 
	cpully->SaveAs( (outpath+"/pully_"  +finname+"_perevt.png").Data() ); 
	cpully->SaveAs( (outpath+"/pully_"  +finname+"_perevt.root").Data() ); 
	cpully->SetLogx(0); 
	cpully->SetLogy(0); 
	cpully->SetGridx(0); 
	cpully->SetGridy(0); 
      } 

      if(k==0 && w==0) { 
	mgpully3->Add(tgaes3[ 1], "lep"); 
	mgpully3->Add(tgaes3[11], "lep"); 
	mgpully3->Add(tgaes3[12], "lep"); 
	mgpully3->Add(tgaes3[13], "lep"); 
	mgpully3->Add(tgaes3[14], "lep"); 
	mgpully3->Add(tgaes3[15], "lep"); 
	mgpully3->Add(tgaes3[16], "lep"); 
	mgpully3->Add(tgaes3[17], "lep"); 
	cpully->cd(); 
	cpully->Clear(); 
	mgpully3->SetTitle(xytitle3.Data()); 
	mgpully3->Draw("a"); 
	mgpully3->GetXaxis()->SetTitleOffset(0.84); 
	mgpully3->GetYaxis()->SetTitleOffset(0.60); 
	legpully->Draw(); 
	//cpully->SetLogy(); 
	cpully->SaveAs( (outpath+"/pully_nmu_perseg_aeta.png").Data() ); 
	cpully->SaveAs( (outpath+"/pully_nmu_perseg_aeta.root").Data() ); 
	cpully->SetGridx(0); 
	cpully->SetGridy(0); 

	mgpully4->Add(tgaes4[ 1], "lep"); 
	mgpully4->Add(tgaes4[11], "lep"); 
	mgpully4->Add(tgaes4[12], "lep"); 
	mgpully4->Add(tgaes4[13], "lep"); 
	mgpully4->Add(tgaes4[14], "lep"); 
	mgpully4->Add(tgaes4[15], "lep"); 
	mgpully4->Add(tgaes4[16], "lep"); 
	mgpully4->Add(tgaes4[17], "lep"); 
	cpully->cd(); 
	cpully->Clear(); 
	mgpully4->SetTitle(xytitle4.Data()); 
	mgpully4->Draw("a"); 
	mgpully4->GetXaxis()->SetTitleOffset(0.84); 
	mgpully4->GetYaxis()->SetTitleOffset(0.60); 
	legpully->Draw(); 
	//cpully->SetLogy(); 
	cpully->SetGridx(); 
	cpully->SetGridy(); 
	cpully->SaveAs( (outpath+"/pully_nmu_perseg_dphi.png").Data() ); 
	cpully->SaveAs( (outpath+"/pully_nmu_perseg_dphi.root").Data() ); 
	cpully->SetGridx(0); 
      } 

      // ======================================================================

      //
      // Diff. y:  1, 18, 19, 20, 21, 22, 23, 24 
      mgdiffy->Add(tgaes[ 1], "lep"); label.Form("#Deltay < %.1f", params[ 1][3]); legdiffy->AddEntry(tgaes[ 1], label.Data(), "LEP"); 
      mgdiffy->Add(tgaes[18], "lep"); label.Form("#Deltay < %.1f", params[18][3]); legdiffy->AddEntry(tgaes[18], label.Data(), "LEP"); 
      mgdiffy->Add(tgaes[19], "lep"); label.Form("#Deltay < %.1f", params[19][3]); legdiffy->AddEntry(tgaes[19], label.Data(), "LEP"); 
      mgdiffy->Add(tgaes[20], "lep"); label.Form("#Deltay < %.1f", params[20][3]); legdiffy->AddEntry(tgaes[20], label.Data(), "LEP"); 
      mgdiffy->Add(tgaes[21], "lep"); label.Form("#Deltay < %.1f", params[21][3]); legdiffy->AddEntry(tgaes[21], label.Data(), "LEP"); 
      mgdiffy->Add(tgaes[22], "lep"); label.Form("#Deltay < %.1f", params[22][3]); legdiffy->AddEntry(tgaes[22], label.Data(), "LEP"); 
      mgdiffy->Add(tgaes[23], "lep"); label.Form("#Deltay < %.1f", params[23][3]); legdiffy->AddEntry(tgaes[23], label.Data(), "LEP"); 
      mgdiffy->Add(tgaes[24], "lep"); label.Form("#Deltay < %.1f", params[24][3]); legdiffy->AddEntry(tgaes[24], label.Data(), "LEP"); 
      cdiffy->cd(); 
      //cdiffy->SetLogy(); 
      cdiffy->SetGridx(); 
      cdiffy->SetGridy(); 
      mgdiffy->SetTitle(xytitle.Data()); 
      mgdiffy->Draw("a"); 
      if(var[w].EqualTo("pt")) { 
	if(k==0) 
	  mgdiffy->GetXaxis()->SetLimits(0.0,  35.0); 
	else { 
	  cdiffy->SetLogx(); 
	  mgdiffy->GetXaxis()->SetLimits(0.3, 200.0); 
	} 
      } 
      mgdiffy->GetXaxis()->SetTitleOffset(0.84); 
      mgdiffy->GetYaxis()->SetTitleOffset(0.60); 
      legdiffy->Draw(); 
      //cdiffy->SetLogy(); 
      cdiffy->SaveAs( (outpath+"/diffy_"  +finname+".png").Data() ); 
      cdiffy->SaveAs( (outpath+"/diffy_"  +finname+".root").Data() ); 
      cdiffy->SetLogx(0); 
      cdiffy->SetGridx(0); 
      cdiffy->SetGridy(0); 

      if(k>0) { 
	mgdiffy2->Add(tgaes2[ 1], "lep"); 
	mgdiffy2->Add(tgaes2[18], "lep"); 
	mgdiffy2->Add(tgaes2[19], "lep"); 
	mgdiffy2->Add(tgaes2[20], "lep"); 
	mgdiffy2->Add(tgaes2[21], "lep"); 
	mgdiffy2->Add(tgaes2[22], "lep"); 
	mgdiffy2->Add(tgaes2[23], "lep"); 
	mgdiffy2->Add(tgaes2[24], "lep"); 
	cdiffy->cd(); 
	cdiffy->Clear(); 
	//cdiffy->SetLogy(); 
	cdiffy->SetGridx(); 
	cdiffy->SetGridy(); 
	mgdiffy2->SetTitle(xytitle2.Data()); 
	mgdiffy2->Draw("a"); 
	if(var[w].EqualTo("pt")) { 
	  mgdiffy2->GetYaxis()->SetRangeUser(10.0, 400.0); 
	  cdiffy->SetLogy(); 
	  mgdiffy2->GetYaxis()->UnZoom(); 
	  cdiffy->Modified(); 
	  cdiffy->Update(); 
	  cdiffy->SetLogx(); 
	  mgdiffy2->GetXaxis()->SetLimits(0.3, 200.0); 
	} 
	mgdiffy2->GetXaxis()->SetTitleOffset(0.84); 
	mgdiffy2->GetYaxis()->SetTitleOffset(0.60); 
	legdiffy->Draw(); 
	//cdiffy->SetLogy(); 
	cdiffy->SaveAs( (outpath+"/diffy_"  +finname+"_perevt.png").Data() ); 
	cdiffy->SaveAs( (outpath+"/diffy_"  +finname+"_perevt.root").Data() ); 
	cdiffy->SetLogx(0); 
	cdiffy->SetLogy(0); 
	cdiffy->SetGridx(0); 
	cdiffy->SetGridy(0); 
      } 

      if(k==0 && w==0) { 
	mgdiffy3->Add(tgaes3[ 1], "lep"); 
	mgdiffy3->Add(tgaes3[18], "lep"); 
	mgdiffy3->Add(tgaes3[19], "lep"); 
	mgdiffy3->Add(tgaes3[20], "lep"); 
	mgdiffy3->Add(tgaes3[21], "lep"); 
	mgdiffy3->Add(tgaes3[22], "lep"); 
	mgdiffy3->Add(tgaes3[23], "lep"); 
	mgdiffy3->Add(tgaes3[24], "lep"); 
	cdiffy->cd(); 
	cdiffy->Clear(); 
	mgdiffy3->SetTitle(xytitle3.Data()); 
	mgdiffy3->Draw("a"); 
	mgdiffy3->GetXaxis()->SetTitleOffset(0.84); 
	mgdiffy3->GetYaxis()->SetTitleOffset(0.60); 
	legdiffy->Draw(); 
	//cdiffy->SetLogy(); 
	cdiffy->SaveAs( (outpath+"/diffy_nmu_perseg_aeta.png").Data() ); 
	cdiffy->SaveAs( (outpath+"/diffy_nmu_perseg_aeta.root").Data() ); 
	cdiffy->SetGridx(0); 
	cdiffy->SetGridy(0); 

	mgdiffy4->Add(tgaes4[ 1], "lep"); 
	mgdiffy4->Add(tgaes4[18], "lep"); 
	mgdiffy4->Add(tgaes4[19], "lep"); 
	mgdiffy4->Add(tgaes4[20], "lep"); 
	mgdiffy4->Add(tgaes4[21], "lep"); 
	mgdiffy4->Add(tgaes4[22], "lep"); 
	mgdiffy4->Add(tgaes4[23], "lep"); 
	mgdiffy4->Add(tgaes4[24], "lep"); 
	cdiffy->cd(); 
	cdiffy->Clear(); 
	mgdiffy4->SetTitle(xytitle4.Data()); 
	mgdiffy4->Draw("a"); 
	mgdiffy4->GetXaxis()->SetTitleOffset(0.84); 
	mgdiffy4->GetYaxis()->SetTitleOffset(0.60); 
	legdiffy->Draw(); 
	//cdiffy->SetLogy(); 
	cdiffy->SetGridx(); 
	cdiffy->SetGridy(); 
	//cdiffy->SetLogy(); 
	cdiffy->SaveAs( (outpath+"/diffy_nmu_perseg_dphi.png").Data() ); 
	cdiffy->SaveAs( (outpath+"/diffy_nmu_perseg_dphi.root").Data() ); 
	cdiffy->SetGridx(0); 
	cdiffy->SetGridy(0); 
      } 

      // Diff. phi: 1, 25, 26, 27, 28, 29 
      mgdiffphi->Add(tgaes[ 1], "lep"); label.Form("#Delta#phi < %.2f", params[ 1][4]); legdiffphi->AddEntry(tgaes[ 1], label.Data(), "LEP"); 
      mgdiffphi->Add(tgaes[25], "lep"); label.Form("#Delta#phi < %.2f", params[25][4]); legdiffphi->AddEntry(tgaes[25], label.Data(), "LEP"); 
      mgdiffphi->Add(tgaes[26], "lep"); label.Form("#Delta#phi < %.2f", params[26][4]); legdiffphi->AddEntry(tgaes[26], label.Data(), "LEP"); 
      mgdiffphi->Add(tgaes[27], "lep"); label.Form("#Delta#phi < %.2f", params[27][4]); legdiffphi->AddEntry(tgaes[27], label.Data(), "LEP"); 
      mgdiffphi->Add(tgaes[28], "lep"); label.Form("#Delta#phi < %.2f", params[28][4]); legdiffphi->AddEntry(tgaes[28], label.Data(), "LEP"); 
      mgdiffphi->Add(tgaes[29], "lep"); label.Form("#Delta#phi < %.2f", params[29][4]); legdiffphi->AddEntry(tgaes[29], label.Data(), "LEP"); 
      cdiffphi->cd(); 
      //cdiffphi->SetLogy(); 
      cdiffphi->SetGridx(); 
      cdiffphi->SetGridy(); 
      mgdiffphi->SetTitle(xytitle.Data()); 
      mgdiffphi->Draw("a"); 
      if(var[w].EqualTo("pt")) { 
	if(k==0) 
	  mgdiffphi->GetXaxis()->SetLimits(0.0,  35.0); 
	else { 
	  cdiffphi->SetLogx(); 
	  mgdiffphi->GetXaxis()->SetLimits(0.3, 200.0); 
	} 
      } 
      mgdiffphi->GetXaxis()->SetTitleOffset(0.84); 
      mgdiffphi->GetYaxis()->SetTitleOffset(0.60); 
      legdiffphi->Draw(); 
      //cdiffphi->SetLogy(); 
      cdiffphi->SaveAs( (outpath+"/diffphi_"  +finname+".png").Data() ); 
      cdiffphi->SaveAs( (outpath+"/diffphi_"  +finname+".root").Data() ); 
      cdiffphi->SetLogx(0); 
      cdiffphi->SetGridx(0); 
      cdiffphi->SetGridy(0); 

      if(k>0) { 
	mgdiffphi2->Add(tgaes2[ 1], "lep"); 
	mgdiffphi2->Add(tgaes2[25], "lep"); 
	mgdiffphi2->Add(tgaes2[26], "lep"); 
	mgdiffphi2->Add(tgaes2[27], "lep"); 
	mgdiffphi2->Add(tgaes2[28], "lep"); 
	mgdiffphi2->Add(tgaes2[29], "lep"); 
	cdiffphi->cd(); 
	cdiffphi->Clear(); 
	//cdiffphi->SetLogy(); 
	cdiffphi->SetGridx(); 
	cdiffphi->SetGridy(); 
	mgdiffphi2->SetTitle(xytitle2.Data()); 
	mgdiffphi2->Draw("a"); 
	if(var[w].EqualTo("pt")) { 
	  mgdiffphi2->GetYaxis()->SetRangeUser(10.0, 400.0); 
	  cdiffphi->SetLogy(); 
	  mgdiffphi2->GetYaxis()->UnZoom(); 
	  cdiffphi->Modified(); 
	  cdiffphi->Update(); 
	  cdiffphi->SetLogx(); 
	  mgdiffphi2->GetXaxis()->SetLimits(0.3, 200.0); 
	} 
	mgdiffphi2->GetXaxis()->SetTitleOffset(0.84); 
	mgdiffphi2->GetYaxis()->SetTitleOffset(0.60); 
	legdiffphi->Draw(); 
	//cdiffphi->SetLogy(); 
	cdiffphi->SaveAs( (outpath+"/diffphi_"  +finname+"_perevt.png").Data() ); 
	cdiffphi->SaveAs( (outpath+"/diffphi_"  +finname+"_perevt.root").Data() ); 
	cdiffphi->SetLogx(0); 
	cdiffphi->SetLogy(0); 
	cdiffphi->SetGridx(0); 
	cdiffphi->SetGridy(0); 
      } 

      if(k==0 && w==0) { 
	mgdiffphi3->Add(tgaes3[ 1], "lep"); 
	mgdiffphi3->Add(tgaes3[25], "lep"); 
	mgdiffphi3->Add(tgaes3[26], "lep"); 
	mgdiffphi3->Add(tgaes3[27], "lep"); 
	mgdiffphi3->Add(tgaes3[28], "lep"); 
	mgdiffphi3->Add(tgaes3[29], "lep"); 
	cdiffphi->cd(); 
	cdiffphi->Clear(); 
	mgdiffphi3->SetTitle(xytitle3.Data()); 
	mgdiffphi3->Draw("a"); 
	mgdiffphi3->GetXaxis()->SetTitleOffset(0.84); 
	mgdiffphi3->GetYaxis()->SetTitleOffset(0.60); 
	legdiffphi->Draw(); 
	//cdiffphi->SetLogy(); 
	cdiffphi->SaveAs( (outpath+"/diffphi_nmu_perseg_aeta.png").Data() ); 
	cdiffphi->SaveAs( (outpath+"/diffphi_nmu_perseg_aeta.root").Data() ); 
	cdiffphi->SetGridx(0); 
	cdiffphi->SetGridy(0); 

	mgdiffphi4->Add(tgaes4[ 1], "lep"); 
	mgdiffphi4->Add(tgaes4[25], "lep"); 
	mgdiffphi4->Add(tgaes4[26], "lep"); 
	mgdiffphi4->Add(tgaes4[27], "lep"); 
	mgdiffphi4->Add(tgaes4[28], "lep"); 
	mgdiffphi4->Add(tgaes4[29], "lep"); 
	cdiffphi->cd(); 
	cdiffphi->Clear(); 
	mgdiffphi4->SetTitle(xytitle4.Data()); 
	mgdiffphi4->Draw("a"); 
	mgdiffphi4->GetXaxis()->SetTitleOffset(0.84); 
	mgdiffphi4->GetYaxis()->SetTitleOffset(0.60); 
	legdiffphi->Draw(); 
	//cdiffphi->SetLogy(); 
	cdiffphi->SetGridx(); 
	cdiffphi->SetGridy(); 
	//cdiffphi->SetLogy(); 
	cdiffphi->SaveAs( (outpath+"/diffphi_nmu_perseg_dphi.png").Data() ); 
	cdiffphi->SaveAs( (outpath+"/diffphi_nmu_perseg_dphi.root").Data() ); 
	cdiffphi->SetGridx(0); 
	cdiffphi->SetGridy(0); 
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

  TLegend *leg = new TLegend(legxmin, legymin, legxmax, legymax); 
  leg->SetFillColor(kWhite);
  leg->SetFillStyle(1001);
  leg->SetLineStyle(0);
  leg->SetLineWidth(0);
  leg->SetBorderSize(0);
  //leg->SetTextFont(42);


  // PullX 
  TGraphAsymmErrors *tgpullx = new TGraphAsymmErrors(6); 
  size_t idx = 0; 
  tgpullx->SetPoint(idx, xxx[0], yyy[0]); tgpullx->SetPointError(idx++, exl[0], exh[0], eyl[0], eyh[0]); 
  tgpullx->SetPoint(idx, xxx[1], yyy[1]); tgpullx->SetPointError(idx++, exl[1], exh[1], eyl[1], eyh[1]); 
  tgpullx->SetPoint(idx, xxx[2], yyy[2]); tgpullx->SetPointError(idx++, exl[2], exh[2], eyl[2], eyh[2]); 
  tgpullx->SetPoint(idx, xxx[3], yyy[3]); tgpullx->SetPointError(idx++, exl[3], exh[3], eyl[3], eyh[3]); 
  tgpullx->SetPoint(idx, xxx[4], yyy[4]); tgpullx->SetPointError(idx++, exl[4], exh[4], eyl[4], eyh[4]); 
  tgpullx->SetPoint(idx, xxx[5], yyy[5]); tgpullx->SetPointError(idx++, exl[5], exh[5], eyl[5], eyh[5]); 
  tgpullx->SetLineColor(kBlack); 
  tgpullx->SetLineWidth(2); 
  tgpullx->SetMarkerStyle(20); 
  tgpullx->SetMarkerColor(kBlack); 
  tgpullx->SetMarkerSize(1.6); 
  tgpullx->Draw("PL"); 
  leg->AddEntry(tgpullx, "Pull x: #Deltax/#sigma_{x}", "PL"); 

  // DiffX 
  TGraphAsymmErrors *tgdiffx = new TGraphAsymmErrors(6); 
  idx = 0; 
  tgdiffx->SetPoint(idx, xxx[ 6], yyy[ 6]); tgdiffx->SetPointError(idx++, exl[ 6], exh[ 6], eyl[ 6], eyh[ 6]); 
  tgdiffx->SetPoint(idx, xxx[ 1], yyy[ 1]); tgdiffx->SetPointError(idx++, exl[ 1], exh[ 1], eyl[ 1], eyh[ 1]); 
  tgdiffx->SetPoint(idx, xxx[ 7], yyy[ 7]); tgdiffx->SetPointError(idx++, exl[ 7], exh[ 7], eyl[ 7], eyh[ 7]); 
  tgdiffx->SetPoint(idx, xxx[ 8], yyy[ 8]); tgdiffx->SetPointError(idx++, exl[ 8], exh[ 8], eyl[ 8], eyh[ 8]); 
  tgdiffx->SetPoint(idx, xxx[ 9], yyy[ 9]); tgdiffx->SetPointError(idx++, exl[ 9], exh[ 9], eyl[ 9], eyh[ 9]); 
  tgdiffx->SetPoint(idx, xxx[10], yyy[10]); tgdiffx->SetPointError(idx++, exl[10], exh[10], eyl[10], eyh[10]); 
  tgdiffx->SetLineColor(kBlue); 
  tgdiffx->SetLineWidth(2); 
  tgdiffx->SetMarkerStyle(21); 
  tgdiffx->SetMarkerColor(kBlue); 
  tgdiffx->SetMarkerSize(1.6); 
  tgdiffx->Draw("PL"); 
  leg->AddEntry(tgdiffx, "Dist. x: #Deltax, no #Deltax/#sigma_{x}", "PL"); 

  // PullY 
  TGraphAsymmErrors *tgpully = new TGraphAsymmErrors(8); 
  idx = 0; 
  tgpully->SetPoint(idx, xxx[ 1], yyy[ 1]); tgpully->SetPointError(idx++, exl[ 1], exh[ 1], eyl[ 1], eyh[ 1]); 
  tgpully->SetPoint(idx, xxx[11], yyy[11]); tgpully->SetPointError(idx++, exl[11], exh[11], eyl[11], eyh[11]); 
  tgpully->SetPoint(idx, xxx[12], yyy[12]); tgpully->SetPointError(idx++, exl[12], exh[12], eyl[12], eyh[12]); 
  tgpully->SetPoint(idx, xxx[13], yyy[13]); tgpully->SetPointError(idx++, exl[13], exh[13], eyl[13], eyh[13]); 
  tgpully->SetPoint(idx, xxx[14], yyy[14]); tgpully->SetPointError(idx++, exl[14], exh[14], eyl[14], eyh[14]); 
  tgpully->SetPoint(idx, xxx[15], yyy[15]); tgpully->SetPointError(idx++, exl[15], exh[15], eyl[15], eyh[15]); 
  tgpully->SetPoint(idx, xxx[16], yyy[16]); tgpully->SetPointError(idx++, exl[16], exh[16], eyl[16], eyh[16]); 
  tgpully->SetPoint(idx, xxx[17], yyy[17]); tgpully->SetPointError(idx++, exl[17], exh[17], eyl[17], eyh[17]); 
  tgpully->SetLineColor(kRed); 
  tgpully->SetLineWidth(2); 
  tgpully->SetMarkerStyle(22); 
  tgpully->SetMarkerColor(kRed); 
  tgpully->SetMarkerSize(1.6); 
  tgpully->Draw("PL"); 
  leg->AddEntry(tgpully, "Pull y: #Deltay/#sigma_{y}", "PL"); 

  // DiffY 
  TGraphAsymmErrors *tgdiffy = new TGraphAsymmErrors(8); 
  idx = 0; 
  tgdiffy->SetPoint(idx, xxx[ 1], yyy[ 1]); tgdiffy->SetPointError(idx++, exl[ 1], exh[ 1], eyl[ 1], eyh[ 1]); 
  tgdiffy->SetPoint(idx, xxx[18], yyy[18]); tgdiffy->SetPointError(idx++, exl[18], exh[18], eyl[18], eyh[18]); 
  tgdiffy->SetPoint(idx, xxx[19], yyy[19]); tgdiffy->SetPointError(idx++, exl[19], exh[19], eyl[19], eyh[19]); 
  tgdiffy->SetPoint(idx, xxx[20], yyy[20]); tgdiffy->SetPointError(idx++, exl[20], exh[20], eyl[20], eyh[20]); 
  tgdiffy->SetPoint(idx, xxx[21], yyy[21]); tgdiffy->SetPointError(idx++, exl[21], exh[21], eyl[21], eyh[21]); 
  tgdiffy->SetPoint(idx, xxx[22], yyy[22]); tgdiffy->SetPointError(idx++, exl[22], exh[22], eyl[22], eyh[22]); 
  tgdiffy->SetPoint(idx, xxx[23], yyy[23]); tgdiffy->SetPointError(idx++, exl[23], exh[23], eyl[23], eyh[23]); 
  tgdiffy->SetPoint(idx, xxx[24], yyy[24]); tgdiffy->SetPointError(idx++, exl[24], exh[24], eyl[24], eyh[24]); 
  tgdiffy->SetLineColor(kOrange); 
  tgdiffy->SetLineWidth(2); 
  tgdiffy->SetMarkerStyle(23); 
  tgdiffy->SetMarkerColor(kOrange); 
  tgdiffy->SetMarkerSize(1.6); 
  tgdiffy->Draw("PL"); 
  leg->AddEntry(tgdiffy, "Dist. y: #Deltay", "PL"); 

  // DiffPhi 
  TGraphAsymmErrors *tgdiffphi = new TGraphAsymmErrors(6); 
  idx = 0; 
  tgdiffphi->SetPoint(idx, xxx[ 1], yyy[ 1]); tgdiffphi->SetPointError(idx++, exl[ 1], exh[ 1], eyl[ 1], eyh[ 1]); 
  tgdiffphi->SetPoint(idx, xxx[25], yyy[25]); tgdiffphi->SetPointError(idx++, exl[25], exh[25], eyl[25], eyh[25]); 
  tgdiffphi->SetPoint(idx, xxx[26], yyy[26]); tgdiffphi->SetPointError(idx++, exl[26], exh[26], eyl[26], eyh[26]); 
  tgdiffphi->SetPoint(idx, xxx[27], yyy[27]); tgdiffphi->SetPointError(idx++, exl[27], exh[27], eyl[27], eyh[27]); 
  tgdiffphi->SetPoint(idx, xxx[28], yyy[28]); tgdiffphi->SetPointError(idx++, exl[28], exh[28], eyl[28], eyh[28]); 
  tgdiffphi->SetPoint(idx, xxx[29], yyy[29]); tgdiffphi->SetPointError(idx++, exl[29], exh[29], eyl[29], eyh[29]); 
  tgdiffphi->SetLineColor(kMagenta); 
  tgdiffphi->SetLineWidth(2); 
  tgdiffphi->SetMarkerStyle(24); 
  tgdiffphi->SetMarkerColor(kMagenta); 
  tgdiffphi->SetMarkerSize(1.6); 
  tgdiffphi->Draw("PL"); 
  leg->AddEntry(tgdiffphi, "Dist. #phi: #Delta#phi", "PL"); 

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
  cc->cd(); 
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

  leg = new TLegend(legxmin, legymin, legxmax, legymax); 
  leg->SetFillColor(kWhite);
  leg->SetFillStyle(1001);
  leg->SetLineStyle(0);
  leg->SetLineWidth(0);
  leg->SetBorderSize(0);
  //leg->SetTextFont(42);


  // PullX 
  tgpullx = new TGraphAsymmErrors(6); 
  idx = 0; 
  tgpullx->SetPoint(idx, xxx[0], zzz[0]); tgpullx->SetPointError(idx++, exl[0], exh[0], ezl[0], ezh[0]); 
  tgpullx->SetPoint(idx, xxx[1], zzz[1]); tgpullx->SetPointError(idx++, exl[1], exh[1], ezl[1], ezh[1]); 
  tgpullx->SetPoint(idx, xxx[2], zzz[2]); tgpullx->SetPointError(idx++, exl[2], exh[2], ezl[2], ezh[2]); 
  tgpullx->SetPoint(idx, xxx[3], zzz[3]); tgpullx->SetPointError(idx++, exl[3], exh[3], ezl[3], ezh[3]); 
  tgpullx->SetPoint(idx, xxx[4], zzz[4]); tgpullx->SetPointError(idx++, exl[4], exh[4], ezl[4], ezh[4]); 
  tgpullx->SetPoint(idx, xxx[5], zzz[5]); tgpullx->SetPointError(idx++, exl[5], exh[5], ezl[5], ezh[5]); 
  tgpullx->SetLineColor(kBlack); 
  tgpullx->SetLineWidth(2); 
  tgpullx->SetMarkerStyle(20); 
  tgpullx->SetMarkerColor(kBlack); 
  tgpullx->SetMarkerSize(1.6); 
  tgpullx->Draw("PL"); 
  leg->AddEntry(tgpullx, "Pull x: #Deltax/#sigma_{x}", "PL"); 

  // DiffX 
  tgdiffx = new TGraphAsymmErrors(6); 
  idx = 0; 
  tgdiffx->SetPoint(idx, xxx[ 6], zzz[ 6]); tgdiffx->SetPointError(idx++, exl[ 6], exh[ 6], ezl[ 6], ezh[ 6]); 
  tgdiffx->SetPoint(idx, xxx[ 1], zzz[ 1]); tgdiffx->SetPointError(idx++, exl[ 1], exh[ 1], ezl[ 1], ezh[ 1]); 
  tgdiffx->SetPoint(idx, xxx[ 7], zzz[ 7]); tgdiffx->SetPointError(idx++, exl[ 7], exh[ 7], ezl[ 7], ezh[ 7]); 
  tgdiffx->SetPoint(idx, xxx[ 8], zzz[ 8]); tgdiffx->SetPointError(idx++, exl[ 8], exh[ 8], ezl[ 8], ezh[ 8]); 
  tgdiffx->SetPoint(idx, xxx[ 9], zzz[ 9]); tgdiffx->SetPointError(idx++, exl[ 9], exh[ 9], ezl[ 9], ezh[ 9]); 
  tgdiffx->SetPoint(idx, xxx[10], zzz[10]); tgdiffx->SetPointError(idx++, exl[10], exh[10], ezl[10], ezh[10]); 
  tgdiffx->SetLineColor(kBlue); 
  tgdiffx->SetLineWidth(2); 
  tgdiffx->SetMarkerStyle(21); 
  tgdiffx->SetMarkerColor(kBlue); 
  tgdiffx->SetMarkerSize(1.6); 
  tgdiffx->Draw("PL"); 
  leg->AddEntry(tgdiffx, "Dist. x: #Deltax", "PL"); 

  // PullY 
  tgpully = new TGraphAsymmErrors(8); 
  idx = 0; 
  tgpully->SetPoint(idx, xxx[ 1], zzz[ 1]); tgpully->SetPointError(idx++, exl[ 1], exh[ 1], ezl[ 1], ezh[ 1]); 
  tgpully->SetPoint(idx, xxx[11], zzz[11]); tgpully->SetPointError(idx++, exl[11], exh[11], ezl[11], ezh[11]); 
  tgpully->SetPoint(idx, xxx[12], zzz[12]); tgpully->SetPointError(idx++, exl[12], exh[12], ezl[12], ezh[12]); 
  tgpully->SetPoint(idx, xxx[13], zzz[13]); tgpully->SetPointError(idx++, exl[13], exh[13], ezl[13], ezh[13]); 
  tgpully->SetPoint(idx, xxx[14], zzz[14]); tgpully->SetPointError(idx++, exl[14], exh[14], ezl[14], ezh[14]); 
  tgpully->SetPoint(idx, xxx[15], zzz[15]); tgpully->SetPointError(idx++, exl[15], exh[15], ezl[15], ezh[15]); 
  tgpully->SetPoint(idx, xxx[16], zzz[16]); tgpully->SetPointError(idx++, exl[16], exh[16], ezl[16], ezh[16]); 
  tgpully->SetPoint(idx, xxx[17], zzz[17]); tgpully->SetPointError(idx++, exl[17], exh[17], ezl[17], ezh[17]); 
  tgpully->SetLineColor(kRed); 
  tgpully->SetLineWidth(2); 
  tgpully->SetMarkerStyle(22); 
  tgpully->SetMarkerColor(kRed); 
  tgpully->SetMarkerSize(1.6); 
  tgpully->Draw("PL"); 
  leg->AddEntry(tgpully, "Dist. y: #Deltay/#sigma_{y}", "PL"); 

  // DiffY 
  tgdiffy = new TGraphAsymmErrors(8); 
  idx = 0; 
  tgdiffy->SetPoint(idx, xxx[ 1], zzz[ 1]); tgdiffy->SetPointError(idx++, exl[ 1], exh[ 1], ezl[ 1], ezh[ 1]); 
  tgdiffy->SetPoint(idx, xxx[18], zzz[18]); tgdiffy->SetPointError(idx++, exl[18], exh[18], ezl[18], ezh[18]); 
  tgdiffy->SetPoint(idx, xxx[19], zzz[19]); tgdiffy->SetPointError(idx++, exl[19], exh[19], ezl[19], ezh[19]); 
  tgdiffy->SetPoint(idx, xxx[20], zzz[20]); tgdiffy->SetPointError(idx++, exl[20], exh[20], ezl[20], ezh[20]); 
  tgdiffy->SetPoint(idx, xxx[21], zzz[21]); tgdiffy->SetPointError(idx++, exl[21], exh[21], ezl[21], ezh[21]); 
  tgdiffy->SetPoint(idx, xxx[22], zzz[22]); tgdiffy->SetPointError(idx++, exl[22], exh[22], ezl[22], ezh[22]); 
  tgdiffy->SetPoint(idx, xxx[23], zzz[23]); tgdiffy->SetPointError(idx++, exl[23], exh[23], ezl[23], ezh[23]); 
  tgdiffy->SetPoint(idx, xxx[24], zzz[24]); tgdiffy->SetPointError(idx++, exl[24], exh[24], ezl[24], ezh[24]); 
  tgdiffy->SetLineColor(kOrange); 
  tgdiffy->SetLineWidth(2); 
  tgdiffy->SetMarkerStyle(23); 
  tgdiffy->SetMarkerColor(kOrange); 
  tgdiffy->SetMarkerSize(1.6); 
  tgdiffy->Draw("PL"); 
  leg->AddEntry(tgdiffy, "Dist. y: #Deltay", "PL"); 

  // DiffPhi 
  tgdiffphi = new TGraphAsymmErrors(6); 
  idx = 0; 
  tgdiffphi->SetPoint(idx, xxx[ 1], zzz[ 1]); tgdiffphi->SetPointError(idx++, exl[ 1], exh[ 1], ezl[ 1], ezh[ 1]); 
  tgdiffphi->SetPoint(idx, xxx[25], zzz[25]); tgdiffphi->SetPointError(idx++, exl[25], exh[25], ezl[25], ezh[25]); 
  tgdiffphi->SetPoint(idx, xxx[26], zzz[26]); tgdiffphi->SetPointError(idx++, exl[26], exh[26], ezl[26], ezh[26]); 
  tgdiffphi->SetPoint(idx, xxx[27], zzz[27]); tgdiffphi->SetPointError(idx++, exl[27], exh[27], ezl[27], ezh[27]); 
  tgdiffphi->SetPoint(idx, xxx[28], zzz[28]); tgdiffphi->SetPointError(idx++, exl[28], exh[28], ezl[28], ezh[28]); 
  tgdiffphi->SetPoint(idx, xxx[29], zzz[29]); tgdiffphi->SetPointError(idx++, exl[29], exh[29], ezl[29], ezh[29]); 
  tgdiffphi->SetLineColor(kMagenta); 
  tgdiffphi->SetLineWidth(2); 
  tgdiffphi->SetMarkerStyle(24); 
  tgdiffphi->SetMarkerColor(kMagenta); 
  tgdiffphi->SetMarkerSize(1.6); 
  tgdiffphi->Draw("PL"); 
  leg->AddEntry(tgdiffphi, "Dist. #phi: #Delta#phi", "PL"); 

  leg->Draw("same"); 
  //cc->SetLogy(); 
  cc->SetGridx(); 
  cc->SetGridy(); 
  cc->SaveAs( (outpath+"/roc_eff_bkg.png").Data() ); 
  cc->SaveAs( (outpath+"/roc_eff_bkg.root").Data() ); 
  cc->SetGridx(0); 
  cc->SetGridy(0); 

  return; 
} 

