
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

//
// Global variables
//

TString allcuts="";                               // Cuts (data, MC)
TString allweights="";                            // Weight (MC)
std::map<TString, TString> allweightspersample;   // Weight (one per MC sample)
std::vector<TString> cutSet;                      // Set of all cuts (data, MC)
std::vector<TString> cutCascade;                  // Chain of all cuts (data, MC)
std::vector<TString> variables;                   // Variables to plot
std::vector<TString> xtitles;                     // Variables to plot
std::vector<TString> ytitles;                     // Variables to plot
std::vector<TString> binnings;                    // Binning and range limits
unsigned int nVars=0;
vector<double> puWeights(36, 1.);                 // Vector of weights for PU


//
// Implementation of some necessary variables
//

void initializeGlobalVariables() {
  allcuts="";
  allweights="";
  cutSet.clear();
  cutCascade.clear();
  variables.clear();
  xtitles.clear();
  ytitles.clear();
  binnings.clear();
  nVars=0;
  return;
}

void addVariable(TString newvar, TString newxtitle, TString newytitle, TString newbin="") {
  variables.push_back(newvar);
  xtitles.push_back(newxtitle);
  if(newytitle.Length()>0) 
    ytitles.push_back(newytitle);
  else
    ytitles.push_back("Events/");
  if(newbin.Length()>0) newbin="("+newbin+")";
  binnings.push_back(newbin);
  nVars=variables.size();
  return;
}


void addCut(TString newcutname, TString newcut) {

  if(allcuts.Length()==0) 
    allcuts+=("("+newcut+")");
  else 
    allcuts+=(" && ("+newcut+")");

  cutSet.push_back(newcutname);
  cutCascade.push_back(allcuts);

  return;
}


void addWeight(TString newweight) {

  if(allweights.Length()==0) 
    allweights+=("("+newweight+")");
  else 
    allweights+=("*("+newweight+")");

  return;
}

void addWeight(TString newweight, TString &oldweight) {

  if(oldweight.Length()==0) 
    oldweight+=("("+newweight+")");
  else 
    oldweight+=("*("+newweight+")");

  return;
}

void addWeight(TString newweight, const char *asample) {

  addWeight(newweight, allweightspersample[asample]); 

  return;
}

void addWeight(TString newweight, TString *samples, unsigned int nSmpsToUse) {

  for(unsigned int j=0; j<nSmpsToUse; ++j) {
    addWeight(newweight, allweightspersample[samples[j].Data()]); 
  }

  return;
}


//
// Other useful variables (could go to another file...)
//

// Colins-Soper angle (a la CMS)
double csCosThetaAbs(float pt1, float eta1, float phi1, float charge1,
		     float pt2, float eta2, float phi2) {

  TLorentzVector mu(0., 0., 0., 0.);
  TLorentzVector mubar(0., 0., 0., 0.);

  if(charge1<0) {
    mu.SetPtEtaPhiM(pt1, eta1, phi1, 0.);
    mubar.SetPtEtaPhiM(pt2, eta2, phi2, 0.);
  }
  else {
    mubar.SetPtEtaPhiM(pt1, eta1, phi1, 0.);
    mu.SetPtEtaPhiM(pt2, eta2, phi2, 0.);
  }

  TLorentzVector Q(mu+mubar);

  double muplus  = 1.0/sqrt(2.0) * (mu.E() + mu.Z());
  double muminus = 1.0/sqrt(2.0) * (mu.E() - mu.Z());

  double mubarplus  = 1.0/sqrt(2.0) * (mubar.E() + mubar.Z());
  double mubarminus = 1.0/sqrt(2.0) * (mubar.E() - mubar.Z());

  double costheta = 2.0 / Q.Mag() / sqrt(pow(Q.Mag(), 2) + pow(Q.Pt(), 2)) * 
    (muplus * mubarminus - muminus * mubarplus);

  return fabs(costheta);

}

double dileptMetDeltaPhi(float dileptX, float dileptY, float metX, float metY) {

  double dileptPhi, metPhi;

  // Find dilepton phi
  if(dileptX==0) {
    if(dileptY>0) {
      dileptPhi=1.570796327; // pi/2, 90°
    }
    else {
      dileptPhi=-1.570796327; // -pi/2, -90°
    }
  }
  else {
    dileptPhi = atan(dileptY/dileptX);
    if(dileptX<0) {
      if(dileptY>0) {
	dileptPhi+=3.141592654;
      }
      else {
	dileptPhi-=3.141592654;
      }
    }
  }

  // Find MET phi
  if(metX==0) {
    if(metY>0) {
      metPhi=1.570796327; // pi/2, 90°
    }
    else {
      metPhi=-1.570796327; // -pi/2, -90°
    }
  }
  else {
    metPhi = atan(metY/metX);
    if(metX<0) {
      if(metY>0) {
	metPhi+=3.141592654;
      }
      else {
	metPhi-=3.141592654;
      }
    }
  }

  double deltaPhi=fabs(dileptPhi-metPhi);
  if(deltaPhi>3.141592654)
    deltaPhi=6.283185307-deltaPhi;

  return deltaPhi;
}

double myDeltaPhi(float phi1, float phi2) {

  double deltaPhi=fabs(phi1-phi2);
  if(deltaPhi>3.141592654)
    deltaPhi=6.283185307-deltaPhi;

  return deltaPhi;

}

double getOverallNorm(float finalState) {    // inv.mass, lept.type
  if( fabs(finalState-1)<0.1 ) {       // mm
    return 0.894879;
  }
  else if( fabs(finalState-2)<0.1 ) {  // ee
    return 0.971806;
  }
  else if( fabs(finalState-3)<0.1 ) {  // em
    return 0.930538;
  }
  else {
    return 1.;
  }
}


double getExtraRecoVertexNorm(float nRecVert) {
  if     ( fabs(nRecVert -  1) < 0.1 ) return 0.0530869;
  else if( fabs(nRecVert -  2) < 0.1 ) return 0.318133;
  else if( fabs(nRecVert -  3) < 0.1 ) return 0.641236;
  else if( fabs(nRecVert -  4) < 0.1 ) return 0.960034;
  else if( fabs(nRecVert -  5) < 0.1 ) return 1.18911;
  else if( fabs(nRecVert -  6) < 0.1 ) return 1.31578;
  else if( fabs(nRecVert -  7) < 0.1 ) return 1.369;
  else if( fabs(nRecVert -  8) < 0.1 ) return 1.38435;
  else if( fabs(nRecVert -  9) < 0.1 ) return 1.39915;
  else if( fabs(nRecVert - 10) < 0.1 ) return 1.41679;
  else if( fabs(nRecVert - 11) < 0.1 ) return 1.46922;
  else if( fabs(nRecVert - 12) < 0.1 ) return 1.51979;
  else if( fabs(nRecVert - 13) < 0.1 ) return 1.59561;
  else if( fabs(nRecVert - 14) < 0.1 ) return 1.66335;
  else if( fabs(nRecVert - 15) < 0.1 ) return 1.74668;
  else if( fabs(nRecVert - 16) < 0.1 ) return 1.82898;
  else if( fabs(nRecVert - 17) < 0.1 ) return 1.89655;
  else if( fabs(nRecVert - 18) < 0.1 ) return 2.00866;
  else if( fabs(nRecVert - 19) < 0.1 ) return 2.10335;
  else if( fabs(nRecVert - 20) < 0.1 ) return 2.28661;
  else if( fabs(nRecVert - 21) < 0.1 ) return 2.56388;
  else if( fabs(nRecVert - 22) < 0.1 ) return 3.06602;
  else if( fabs(nRecVert - 23) < 0.1 ) return 3.38279;
  else if( fabs(nRecVert - 24) < 0.1 ) return 3.78195;
  else if( fabs(nRecVert - 25) < 0.1 ) return 4.63139;
  else if( fabs(nRecVert - 26) < 0.1 ) return 4.85099;
  else if( fabs(nRecVert - 27) < 0.1 ) return 4.91525;
  else if( fabs(nRecVert - 28) < 0.1 ) return 6.58376;
  else if( fabs(nRecVert - 29) < 0.1 ) return 6.70113;
  else if( fabs(nRecVert - 30) < 0.1 ) return 6.63118;
  else if( fabs(nRecVert - 31) < 0.1 ) return 6.24989;
  else if( fabs(nRecVert - 32) < 0.1 ) return 7.98107;
  else if( fabs(nRecVert - 33) < 0.1 ) return 6.62338;
  else if( fabs(nRecVert - 34) < 0.1 ) return 2.06608;
  else if( fabs(nRecVert - 35) < 0.1 ) return 11.3634;
  else                                 return 1.;
}

double getExtraGoodRecoVertexNorm(float nRecVert) {
  if     ( fabs(nRecVert -  1) < 0.1 ) return 0.0715549;
  else if( fabs(nRecVert -  2) < 0.1 ) return 0.383385;
  else if( fabs(nRecVert -  3) < 0.1 ) return 0.725732;
  else if( fabs(nRecVert -  4) < 0.1 ) return 1.03871;
  else if( fabs(nRecVert -  5) < 0.1 ) return 1.2404;
  else if( fabs(nRecVert -  6) < 0.1 ) return 1.34387;
  else if( fabs(nRecVert -  7) < 0.1 ) return 1.38838;
  else if( fabs(nRecVert -  8) < 0.1 ) return 1.41323;
  else if( fabs(nRecVert -  9) < 0.1 ) return 1.4463;
  else if( fabs(nRecVert - 10) < 0.1 ) return 1.47697;
  else if( fabs(nRecVert - 11) < 0.1 ) return 1.55544;
  else if( fabs(nRecVert - 12) < 0.1 ) return 1.62036;
  else if( fabs(nRecVert - 13) < 0.1 ) return 1.67083;
  else if( fabs(nRecVert - 14) < 0.1 ) return 1.73942;
  else if( fabs(nRecVert - 15) < 0.1 ) return 1.79357;
  else if( fabs(nRecVert - 16) < 0.1 ) return 1.80336;
  else if( fabs(nRecVert - 17) < 0.1 ) return 1.84698;
  else if( fabs(nRecVert - 18) < 0.1 ) return 1.92475;
  else if( fabs(nRecVert - 19) < 0.1 ) return 2.04169;
  else if( fabs(nRecVert - 20) < 0.1 ) return 2.2201;
  else if( fabs(nRecVert - 21) < 0.1 ) return 2.59883;
  else if( fabs(nRecVert - 22) < 0.1 ) return 2.99641;
  else if( fabs(nRecVert - 23) < 0.1 ) return 3.40613;
  else if( fabs(nRecVert - 24) < 0.1 ) return 3.95436;
  else if( fabs(nRecVert - 25) < 0.1 ) return 2.80663;
  else if( fabs(nRecVert - 26) < 0.1 ) return 4.69622;
  else if( fabs(nRecVert - 27) < 0.1 ) return 4.60286;
  else if( fabs(nRecVert - 28) < 0.1 ) return 5.91928;
  else if( fabs(nRecVert - 29) < 0.1 ) return 2.34438;
  else if( fabs(nRecVert - 30) < 0.1 ) return 10.112;
  else if( fabs(nRecVert - 31) < 0.1 ) return 4.94807;
  else if( fabs(nRecVert - 32) < 0.1 ) return 0.757291;
  else if( fabs(nRecVert - 33) < 0.1 ) return 4.54375;
  else if( fabs(nRecVert - 34) < 0.1 ) return 2.27187;
  else if( fabs(nRecVert - 35) < 0.1 ) return 0.;
  else                                 return 1.;
}

double getUnclProjLongNorm(float unclProjLong) { // ee, +-, Z mass, jet number <= 1
  if     ( fabs(unclProjLong-(-45))<=5. ) return 1.65135;
  else if( fabs(unclProjLong-(-35))<=5. ) return 1.35709;
  else if( fabs(unclProjLong-(-25))<=5. ) return 1.34884;
  else if( fabs(unclProjLong-(-15))<=5. ) return 1.17875;
  else if( fabs(unclProjLong-(-5))<=5.  ) return 1.01814;
  else if( fabs(unclProjLong-(5))<=5.   ) return 0.953373;
  else if( fabs(unclProjLong-(15))<=5.  ) return 0.997537;
  else if( fabs(unclProjLong-(25))<=5.  ) return 1.02733;
  else if( fabs(unclProjLong-(35))<=5.  ) return 1.02469;
  else if( fabs(unclProjLong-(45))<=5.  ) return 0.967728;
  else if( fabs(unclProjLong-(55))<=5.  ) return 0.937805;
  else if( fabs(unclProjLong-(65))<=5.  ) return 0.845883;
  else if( fabs(unclProjLong-(75))<=5.  ) return 0.936875;
  else if( fabs(unclProjLong-(85))<=5.  ) return 0.924292;
  else if( fabs(unclProjLong-(95))<=5.  ) return 0.846026;
  else if( fabs(unclProjLong-(105))<=5. ) return 0.869044;
  else if( fabs(unclProjLong-(115))<=5. ) return 0.835796;
  else if( fabs(unclProjLong-(125))<=5. ) return 1.08366;
  else if( fabs(unclProjLong-(135))<=5. ) return 1.08037;
  else if( fabs(unclProjLong-(145))<=5. ) return 0.855114;
  else                                    return 1.;
}

double getUnclProjPerpNorm(float unclProjPerp) { // ee, +-, Z mass, jet number <= 1, after getUnclProjLongNorm re-weighting
  if     ( fabs(unclProjPerp-(-55))<=5. ) return 1.28321;
  else if( fabs(unclProjPerp-(-45))<=5. ) return 1.22922;
  else if( fabs(unclProjPerp-(-35))<=5. ) return 1.35677;
  else if( fabs(unclProjPerp-(-25))<=5. ) return 1.26899;
  else if( fabs(unclProjPerp-(-15))<=5. ) return 1.16051;
  else if( fabs(unclProjPerp-( -5))<=5. ) return 1.01475;
  else if( fabs(unclProjPerp-(  5))<=5. ) return 0.956595;
  else if( fabs(unclProjPerp-( 15))<=5. ) return 1.00224;
  else if( fabs(unclProjPerp-( 25))<=5. ) return 1.03097;
  else if( fabs(unclProjPerp-( 35))<=5. ) return 1.0205;
  else if( fabs(unclProjPerp-( 45))<=5. ) return 1.01042;
  else if( fabs(unclProjPerp-( 55))<=5. ) return 0.971337;
  else if( fabs(unclProjPerp-( 65))<=5. ) return 0.994124;
  else if( fabs(unclProjPerp-( 75))<=5. ) return 0.953588;
  else if( fabs(unclProjPerp-( 85))<=5. ) return 0.9279;
  else if( fabs(unclProjPerp-( 95))<=5. ) return 0.930336;
  else if( fabs(unclProjPerp-(105))<=5. ) return 0.842078;
  else if( fabs(unclProjPerp-(115))<=5. ) return 0.874163;
  else if( fabs(unclProjPerp-(125))<=5. ) return 0.860305;
  else if( fabs(unclProjPerp-(135))<=5. ) return 0.92118;
  else if( fabs(unclProjPerp-(170))<=30.) return 0.741493;
  else if( fabs(unclProjPerp-(230))<=30.) return 0.710686;
  else if( fabs(unclProjPerp-(290))<=30.) return 0.747752;
  //else if( fabs(unclProjPerp-(350))<=30.) return 2.13501;
  else                                    return 1.;
}

double getRecoilOppJetsLongNorm(float recoilOppJetsLong) { // mm, +-, Z mass, jet number <= 1
  if     ( recoilOppJetsLong<(-157.5-2.5) ) return 1.;
  else if( recoilOppJetsLong<(-157.5+2.5) ) return 0.902441;
  else if( recoilOppJetsLong<(-152.5+2.5) ) return 1.91396;
  else if( recoilOppJetsLong<(-147.5+2.5) ) return 1.2542;
  else if( recoilOppJetsLong<(-142.5+2.5) ) return 0.850504;
  else if( recoilOppJetsLong<(-137.5+2.5) ) return 1.07945;
  else if( recoilOppJetsLong<(-132.5+2.5) ) return 1.19874;
  else if( recoilOppJetsLong<(-127.5+2.5) ) return 1.27225;
  else if( recoilOppJetsLong<(-122.5+2.5) ) return 0.946466;
  else if( recoilOppJetsLong<(-117.5+2.5) ) return 0.977706;
  else if( recoilOppJetsLong<(-112.5+2.5) ) return 0.898443;
  else if( recoilOppJetsLong<(-107.5+2.5) ) return 1.04547;
  else if( recoilOppJetsLong<(-102.5+2.5) ) return 0.858365;
  else if( recoilOppJetsLong<( -97.5+2.5) ) return 0.784349;
  else if( recoilOppJetsLong<( -92.5+2.5) ) return 1.03814;
  else if( recoilOppJetsLong<( -87.5+2.5) ) return 0.899264;
  else if( recoilOppJetsLong<( -82.5+2.5) ) return 0.990597;
  else if( recoilOppJetsLong<( -77.5+2.5) ) return 0.874477;
  else if( recoilOppJetsLong<( -72.5+2.5) ) return 0.988083;
  else if( recoilOppJetsLong<( -67.5+2.5) ) return 1.00418;
  else if( recoilOppJetsLong<( -62.5+2.5) ) return 0.972143;
  else if( recoilOppJetsLong<( -57.5+2.5) ) return 0.950922;
  else if( recoilOppJetsLong<( -52.5+2.5) ) return 0.930473;
  else if( recoilOppJetsLong<( -47.5+2.5) ) return 0.969212;
  else if( recoilOppJetsLong<( -42.5+2.5) ) return 1.0085;
  else if( recoilOppJetsLong<( -37.5+2.5) ) return 1.00657;
  else if( recoilOppJetsLong<( -32.5+2.5) ) return 1.02349;
  else if( recoilOppJetsLong<( -27.5+2.5) ) return 1.01408;
  else if( recoilOppJetsLong<( -22.5+2.5) ) return 1.02057;
  else if( recoilOppJetsLong<( -17.5+2.5) ) return 1.00415;
  else if( recoilOppJetsLong<( -12.5+2.5) ) return 0.989669;
  else if( recoilOppJetsLong<(  -7.5+2.5) ) return 0.96876;
  else if( recoilOppJetsLong<(  -2.5+2.5) ) return 0.949391;
  else if( recoilOppJetsLong<(   2.5+2.5) ) return 1.05383;
  else                                      return 1.; 
}

double getRecoilOppJetsPerpNorm(float recoilOppJetsPerp) { // mm, +-, Z mass, jet number <= 1, after getRecoilOppJetsLongNorm re-weighting
  if     ( recoilOppJetsPerp<(-155.-5.) ) return 1.;
  else if( recoilOppJetsPerp<(-155.+5.) ) return 1.02585;
  else if( recoilOppJetsPerp<(-145.+5.) ) return 0.927008;
  else if( recoilOppJetsPerp<(-135.+5.) ) return 0.892117;
  else if( recoilOppJetsPerp<(-125.+5.) ) return 0.917091;
  else if( recoilOppJetsPerp<(-115.+5.) ) return 1.00573;
  else if( recoilOppJetsPerp<(-105.+5.) ) return 0.936748;
  else if( recoilOppJetsPerp<( -95.+5.) ) return 0.994819;
  else if( recoilOppJetsPerp<( -85.+5.) ) return 1.00274;
  else if( recoilOppJetsPerp<( -75.+5.) ) return 0.956078;
  else if( recoilOppJetsPerp<( -65.+5.) ) return 0.954484;
  else if( recoilOppJetsPerp<( -55.+5.) ) return 0.96165;
  else if( recoilOppJetsPerp<( -45.+5.) ) return 1.02195;
  else if( recoilOppJetsPerp<( -35.+5.) ) return 1.02367;
  else if( recoilOppJetsPerp<( -25.+5.) ) return 1.0279;
  else if( recoilOppJetsPerp<( -15.+5.) ) return 0.995428;
  else if( recoilOppJetsPerp<(  -5.+5.) ) return 0.956045;
  else if( recoilOppJetsPerp<(   5.+5.) ) return 1.05303;
  else                                  return 1.;
}


double getD0RedMet(float flav, 
		   float dileptLong, float dileptUncLong, float metPlusDileptLong, float sumJetsLong, 
		   float dileptPerp, float dileptUncPerp, float metPlusDileptPerp, float sumJetsPerp,
		   int pickAFlav=0 ) {

  if( fabs(flav-3.)<0.1) { 
    if( pickAFlav!=1 && pickAFlav!=2 ) {
      cout << " *** ERROR *** " << endl;
      cout << "  You need to pick a flavor in getD0RedMet(...)! " << endl;
      throw std::exception();
      return -1.;
    }
    else {
      flav = pickAFlav;
    }
  }


  float wBisMu = 1.0;
  float wDilMu = 1.0;
  float wRecMu = 2.0;
  float wUncMu = 2.5;

  float wBisEl = 1.5;
  float wDilEl = 1.0;
  float wRecEl = 2.25;
  float wUncEl = 0.0;

  float wBis = 1.0;
  float wDil = 1.0;
  float wRec = 1.0;
  float wUnc = 1.0;

  if( fabs(flav-1.)<0.1 ) {        // mm
    wBis = wBisMu;
    wDil = wDilMu;
    wRec = wRecMu;
    wUnc = wUncMu;
  }
  else if ( fabs(flav-2.)<0.1 ) {  // ee
    wBis = wBisEl;
    wDil = wDilEl;
    wRec = wRecEl;
    wUnc = wUncEl;
  }
  else {}

  metPlusDileptLong *= -1.;
  metPlusDileptPerp *= -1.;

  float recoilLong = ( sumJetsLong<metPlusDileptLong ? sumJetsLong : metPlusDileptLong ); 
  recoilLong = ( recoilLong<0. ? recoilLong : 0. ); 
  float recoilPerp = ( sumJetsPerp<metPlusDileptPerp ? sumJetsPerp : metPlusDileptPerp ); 
  recoilPerp = ( recoilPerp<0. ? recoilPerp : 0. ); 

  float redMetLong = wDil*dileptLong + wRec*recoilLong + wUnc*dileptUncLong;
  redMetLong = ( redMetLong>=0. ? redMetLong : 0. );
  float redMetPerp = wDil*dileptPerp + wRec*recoilPerp + wUnc*dileptUncPerp;
  redMetPerp = ( redMetPerp>=0. ? redMetPerp : 0. );

  return sqrt( redMetLong*redMetLong + wBis*redMetPerp*redMetPerp ); 
}

double getCmsRedMet(float flav, 
		    float dileptLong, float dileptUncLong, float metPlusDileptLong, float sumJetsLong, 
		    float dileptPerp, float dileptUncPerp, float metPlusDileptPerp, float sumJetsPerp, 
		    int pickAFlav=0 ) {

  if( fabs(flav-3.)<0.1) { 
    if( pickAFlav!=1 && pickAFlav!=2 ) {
      cout << " *** ERROR *** " << endl;
      cout << "  You need to pick a flavor in getCmsRedMet(...)! " << endl;
      throw std::exception();
      return -1.;
    }
    else {
      flav = pickAFlav;
    }
  }

  
  float wBisMu = 1.00;
  float wDilMu = 1.00;
  float wRecMu = 1.00;
  float wUncMu = 1.00;

  float wBisEl = 0.75;
  float wDilEl = 1.00;
  float wRecEl = 1.00;
  float wUncEl = 0.00;

  float wBis = 1.0;
  float wDil = 1.0;
  float wRec = 1.0;
  float wUnc = 1.0;

  if( fabs(flav-1.)<0.1 ) {        // mm
    wBis = wBisMu;
    wDil = wDilMu;
    wRec = wRecMu;
    wUnc = wUncMu;
  }
  else if ( fabs(flav-2.)<0.1 ) {  // ee
    wBis = wBisEl;
    wDil = wDilEl;
    wRec = wRecEl;
    wUnc = wUncEl;
  }
  else {}

  float redMetLong_withJets = wDil*dileptLong + wRec*sumJetsLong       + wUnc*dileptUncLong;
  float redMetLong_withMet  = wDil*dileptLong - wRec*metPlusDileptLong + wUnc*dileptUncLong;

  float redMetPerp_withJets = wDil*dileptPerp + wRec*sumJetsPerp       + wUnc*dileptUncPerp;
  float redMetPerp_withMet  = wDil*dileptPerp - wRec*metPlusDileptPerp + wUnc*dileptUncPerp;

  float redMetLong_min = ( redMetLong_withJets<redMetLong_withMet ? redMetLong_withJets : redMetLong_withMet );
  float redMetPerp_min = ( redMetPerp_withJets<redMetPerp_withMet ? redMetPerp_withJets : redMetPerp_withMet );

  return sqrt( redMetLong_min*redMetLong_min + wBis*redMetPerp_min*redMetPerp_min ); 
}



void makeWeightDistribution(TString mcScenario, TString dataScenario, vector<double> & result) {

  vector<double> mcArr(result.size(), 1.);
  vector<double> dtArr(result.size(), 1.);

  if(mcScenario.CompareTo("mc_observed_summer11_s4_35")==0) {
    double mcArrTmp[36] = {1.45346E-01, 6.42802E-02, 6.95255E-02, 6.96747E-02, 6.92955E-02, 6.84997E-02, 6.69528E-02, 6.45515E-02, 6.09865E-02, 5.63323E-02, 5.07322E-02, 4.44681E-02, 3.79205E-02, 3.15131E-02, 2.54220E-02, 2.00184E-02, 1.53776E-02, 1.15387E-02, 8.47608E-03, 6.08715E-03, 4.28255E-03, 2.97185E-03, 2.01918E-03, 1.34490E-03, 8.81587E-04, 5.69954E-04, 3.61493E-04, 2.28692E-04, 1.40791E-04, 8.44606E-05, 5.10204E-05, 3.07802E-05, 1.81401E-05, 1.00201E-05, 5.80004E-06, 0};
    for(unsigned int i=0; i<mcArr.size(); ++i) {
      mcArr[i] = mcArrTmp[i];
    }
  }

  if(dataScenario.CompareTo("data_observed_prompt_all2011_35")==0) {
    double dtArrTmp[36] = {0.00285942, 0.0125603, 0.0299631, 0.051313, 0.0709713, 0.0847864, 0.0914627, 0.0919255, 0.0879994, 0.0814127, 0.0733995, 0.0647191, 0.0558327, 0.0470663, 0.0386988, 0.0309811, 0.0241175, 0.018241, 0.0133997, 0.00956071, 0.00662814, 0.00446735, 0.00292946, 0.00187057, 0.00116414, 0.000706805, 0.000419059, 0.000242856, 0.0001377, 7.64582e-05, 4.16101e-05, 2.22135e-05, 1.16416e-05, 5.9937e-06, 5.95541e-06, 1.70121e-12};
    for(unsigned int i=0; i<dtArr.size(); ++i) {
      dtArr[i] = dtArrTmp[i];
    }
  }

  if(dataScenario.CompareTo("data_observed_may10_160404_163869_35")==0) {
    double dtArrTmp[36] = {0.00654107, 0.0282192, 0.0658972, 0.108633, 0.140701, 0.151866, 0.141831, 0.117629, 0.0882906, 0.0608404, 0.0389243, 0.0233308, 0.0132003, 0.00709547, 0.00364394, 0.00179698, 0.00085481, 0.000393846, 0.000176402, 7.70529e-05, 3.29119e-05, 1.37768e-05, 5.66111e-06, 2.28622e-06, 9.08062e-07, 3.54848e-07, 1.36434e-07, 5.16043e-08, 1.91956e-08, 7.01957e-09, 2.52253e-09, 8.90441e-10, 3.0865e-10, 1.05021e-10, 5.19305e-11, 2.77965e-11};
    for(unsigned int i=0; i<dtArr.size(); ++i) {
      dtArr[i] = dtArrTmp[i];
    }
  }

  if(dataScenario.CompareTo("data_observed_prompt_165088_167913_35")==0) {
    double dtArrTmp[36] = {0.00719213, 0.0313564, 0.0723671, 0.116859, 0.147734, 0.155388, 0.141249, 0.113879, 0.0829699, 0.0554065, 0.0342899, 0.0198428, 0.0108156, 0.00558697, 0.0027494, 0.00129469, 0.000585616, 0.000255268, 0.000107529, 4.38764e-05, 1.73774e-05, 6.69158e-06, 2.50888e-06, 9.16976e-07, 3.27043e-07, 1.13916e-07, 3.87805e-08, 1.29107e-08, 4.20543e-09, 1.34087e-09, 4.1863e-10, 1.28019e-10, 3.83498e-11, 1.12636e-11, 4.49471e-12, 2.14544e-12};
    for(unsigned int i=0; i<dtArr.size(); ++i) {
      dtArr[i] = dtArrTmp[i];
    }
  }

  if(dataScenario.CompareTo("data_observed_aug05_170249_172619_35")==0) {
    double dtArrTmp[36] = {0.00563937, 0.0229973, 0.0517398, 0.0847196, 0.112657, 0.128767, 0.130759, 0.120404, 0.101874, 0.0799242, 0.0585292, 0.0402185, 0.0260466, 0.0159597, 0.00928456, 0.00514434, 0.00272259, 0.00137995, 0.000671453, 0.000314329, 0.000141851, 6.18205e-05, 2.6061e-05, 1.06425e-05, 4.21557e-06, 1.62162e-06, 6.06438e-07, 2.2069e-07, 7.82201e-08, 2.7023e-08, 9.10621e-09, 2.99511e-09, 9.62094e-10, 3.01987e-10, 1.31916e-10, 0};
    for(unsigned int i=0; i<dtArr.size(); ++i) {
      dtArr[i] = dtArrTmp[i];
    }
  }

  if(dataScenario.CompareTo("data_observed_prompt_172620_173692_35")==0) {
    double dtArrTmp[36] = {0.00412608, 0.0179609, 0.042596, 0.0723529, 0.0992189, 0.117116, 0.12351, 0.118964, 0.106031, 0.0881499, 0.0687174, 0.0504288, 0.0349536, 0.0229512, 0.0143164, 0.00850614, 0.00482603, 0.00262082, 0.00136533, 0.000683752, 0.000329812, 0.000153508, 6.90612e-05, 3.00792e-05, 1.2702e-05, 5.20793e-06, 2.07592e-06, 8.05468e-07, 3.04564e-07, 1.12351e-07, 4.04756e-08, 1.42541e-08, 4.91145e-09, 1.65718e-09, 8.07636e-10, 0};
    for(unsigned int i=0; i<dtArr.size(); ++i) {
      dtArr[i] = dtArrTmp[i];
    }
  }

  if(dataScenario.CompareTo("data_observed_prompt_175832_177515_35")==0) {
    double dtArrTmp[36] = {0.000166117, 0.00113211, 0.00420346, 0.0109629, 0.0224119, 0.0382129, 0.0565315, 0.0745463, 0.0893209, 0.0986259, 0.101425, 0.0979401, 0.0893745, 0.0774687, 0.0640481, 0.0506811, 0.0384942, 0.0281332, 0.0198261, 0.0134978, 0.00889231, 0.00567749, 0.00351801, 0.00211841, 0.00124119, 0.00070844, 0.000394371, 0.000214353, 0.000113879, 5.9197e-05, 3.01389e-05, 1.50432e-05, 7.36755e-06, 3.54354e-06, 3.0957e-06, 0};
    for(unsigned int i=0; i<dtArr.size(); ++i) {
      dtArr[i] = dtArrTmp[i];
    }
  }

  if(dataScenario.CompareTo("data_observed_prompt_177718_178078_35")==0) {
    double dtArrTmp[36] = {0.000335355, 0.00221611, 0.00774264, 0.018714, 0.0351157, 0.0546301, 0.0735665, 0.0884342, 0.0971149, 0.0991616, 0.0954114, 0.0873614, 0.0766446, 0.0647211, 0.0527504, 0.0415661, 0.0316972, 0.0234085, 0.0167518, 0.0116242, 0.007827, 0.00511818, 0.00325322, 0.00201189, 0.00121177, 0.000711534, 0.000407723, 0.000228219, 0.000124901, 6.68955e-05, 3.50921e-05, 1.80444e-05, 9.10131e-06, 4.50583e-06, 4.14068e-06, 0};
    for(unsigned int i=0; i<dtArr.size(); ++i) {
      dtArr[i] = dtArrTmp[i];
    }
  }

  if(dataScenario.CompareTo("data_observed_prompt_178098_180252_35")==0) {
    double dtArrTmp[36] = {0.000169652, 0.00112017, 0.00399151, 0.00996989, 0.0195683, 0.0322346, 0.0465274, 0.0606406, 0.0729092, 0.0820905, 0.0874425, 0.0887038, 0.0860508, 0.0800375, 0.0715024, 0.0614359, 0.0508299, 0.0405435, 0.0312141, 0.0232244, 0.0167204, 0.0116631, 0.00789212, 0.0051872, 0.00331564, 0.00206356, 0.00125195, 0.000741254, 0.000428768, 0.000242548, 0.000134312, 7.28737e-05, 3.87733e-05, 2.0246e-05, 2.06084e-05, 0};
    for(unsigned int i=0; i<dtArr.size(); ++i) {
      dtArr[i] = dtArrTmp[i];
    }
  }

  for(unsigned int k=0; k<result.size(); ++k) {
    if(mcArr[k]>0)
      result[k] = dtArr[k]/mcArr[k];
    else 
      result[k] = 0.;

    //cout << dtArr[k] << "/" << mcArr[k] << " = " << result[k] << endl;
  }

  return;

}


double getPuWeights(unsigned int nInt) {
  if(nInt>=puWeights.size()) return 1.;
  return puWeights[nInt];
}

double getWzInefficiency(float evttype, float flav, float run, float eta) {

  double runRanges[5] = {160329, 165071, 170053, 172620, 175832};
  double etaBins[9] = {-2.1, -1.6, -1.1, -0.6, 0., 0.6, 1.1, 1.6, 2.1};

  double acceptance_m_mm = 0.609138; 
  double acceptance_e_mm = 0.869258; 
  double acceptance_m_ee = 0.744716; 
  double acceptance_e_ee = 0.810419; 

  unsigned int irun(0), ieta(0);

  for(; irun<5; ++irun) {
    if( run<runRanges[irun] ) break;
  }

  for(; ieta<9; ++ieta) {
    if( eta<etaBins[ieta] ) break;
  }

  if(irun>5) {
    std::cout << " *** ERROR ***" << std::endl;
    std::cout << "    irun>5 !" << std::endl;
    return 1;
  }

  if(ieta>9) {
    std::cout << " *** ERROR ***" << std::endl;
    std::cout << "    ieta>9 !" << std::endl;
    return 1;
  }

  double corrections_mu[6][10] = 
    {
      {0.997584, 0.99708,  0.993657, 0.992452, 0.99172,  0.992198, 0.992336, 0.993421, 0.997389, 0.997304},
      {0.998929, 0.998697, 0.997065, 0.995065, 0.994154, 0.995496, 0.994854, 0.997653, 0.997718, 0.997766},
      {0.998713, 0.997658, 0.996679, 0.996187, 0.994774, 0.995086, 0.995279, 0.996717, 0.997765, 0.997849},
      {0.997817, 0.997813, 0.996363, 0.995094, 0.994658, 0.995369, 0.995144, 0.996646, 0.997981, 0.997497},
      {0.997709, 0.997468, 0.997484, 0.994442, 0.994166, 0.995559, 0.995076, 0.997492, 0.998003, 0.997578},
      {0.998261, 0.998958, 0.998451, 0.998201, 0.997763, 0.997716, 0.997871, 0.998543, 0.998899, 0.997428}
    };

  double corrections_e[6][10] = 
    {
      {0.995, 0.995, 0.995, 0.995, 0.995, 0.995, 0.995, 0.995, 0.995, 0.995},
      {0.995, 0.995, 0.995, 0.995, 0.995, 0.995, 0.995, 0.995, 0.995, 0.995},
      {0.995, 0.995, 0.995, 0.995, 0.995, 0.995, 0.995, 0.995, 0.995, 0.995},
      {0.995, 0.995, 0.995, 0.995, 0.995, 0.995, 0.995, 0.995, 0.995, 0.995},
      {0.995, 0.995, 0.995, 0.995, 0.995, 0.995, 0.995, 0.995, 0.995, 0.995},
      {0.995, 0.995, 0.995, 0.995, 0.995, 0.995, 0.995, 0.995, 0.995, 0.995}
    };

  double eff = 1.;

  if( fabs( fabs(flav) - 13 ) < 0.1 ) {
    eff = corrections_mu[irun][ieta];
    if( fabs(evttype-1)<0.1)       // mm
      eff *= acceptance_m_mm;
    else if( fabs(evttype-2)<0.1)  // ee
      eff *= acceptance_m_ee;
    else {}
  }
  else if( fabs( fabs(flav) - 11 ) < 0.1 ) {
    eff = corrections_e[irun][ieta];
    if( fabs(evttype-1)<0.1)       // mm
      eff *= acceptance_e_mm;
    else if( fabs(evttype-2)<0.1)  // ee
      eff *= acceptance_e_ee;
    else {}
  }
  else {}

  return ((1-eff)/eff);

}

