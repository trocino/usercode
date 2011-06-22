
/*
 *  See header file for a description of this class.
 *
 *  $Date: 2008/03/13 17:42:59 $
 *  $Revision: 1.3 $
 *  \author G. Cerminara - CERN
 *          D. Trocino - Northeastern University
 */

#include "LumiNorm.h"

#include <iostream>

#include "TH1F.h"
#include "TFile.h"
#include "TSystem.h"

using namespace std;



LumiNorm::LumiNorm(const TString& inputDir, const TString& lumiFileName,
		   bool rescaleLumi, bool rescaleXsect) 
  : nData(0),
    nMC(0),
    theLuminosity(0.),
    theInputDir(inputDir),
    theLumiFileName(lumiFileName),
    theLumiHistoName("PreEventCounter"),
    theNormHistoName("cut2_MassWindow/DileptKin_cut2_MassWindow_hMass"),
    additionalScale(1.),
    normalizeToZ(true), 
    rescaleLuminosity(rescaleLumi),
    rescaleCrossSect(rescaleXsect)
{
  // Initialize the converion map
  conversionMap["data"]="AllData";
  conversionMap["zz_llnn"]="ZZtoAnything_only2l2n_Spring11";
  conversionMap["zz_x"]="ZZtoAnything_allBut2l2n_Spring11";
  conversionMap["wz"]="WZtoAnything_Spring11";
  conversionMap["ww"]="WWtoAnything_Spring11";
  conversionMap["wjets"]="WJetsToLNu_Spring11";
  conversionMap["tt"]="TTJets_madgraph_Spring11";
  conversionMap["t"]="TToBLNu_Spring11";
  conversionMap["dy_ee"]="DYToEE_M-20_Spring11";
  conversionMap["dy_mm"]="DYToMuMu_M-20_Spring11";
  conversionMap["dy_tt"]="DYToTauTau_M-20_Spring11";
  conversionMap["qcd"]="QCD_Pt-30to1000_Spring11";
}


LumiNorm::~LumiNorm() {}


// Add data to the count of events to be used for data/MC normalization
// the histo used to get the count is zmass_restricted_1
// the default file for this histo is: histos_zzanalysis_data_FINALSTATE.root
/*
void LumiNorm::addData() {
  TString fileName = theInputDir+"/histos_zzanalysis_data_"+theFinalState+".root";
  addData(fileName);
}
*/
void LumiNorm::addData() {
  TString thisFileName=theLumiFileName;
  thisFileName.ReplaceAll("whatsample", conversionMap["data"]);
  TString fileName=theInputDir+"/"+thisFileName;
  addData(fileName);
}


// Add data to the count of events to be used for data/MC normalization
// the histo used to get the count is zmass_restricted_1
// This method allows to specify the name of the file which contains this histo
void LumiNorm::addData(const TString& filename) {
  FileStat_t buffer;
  if(gSystem->GetPathInfo(filename.Data(), buffer) != 0) {
    cout << "[LumiNorm] *** Error: Input file " <<  filename << " doesn't exist!" << endl;
    throw std::exception();
  }
  TFile file(filename.Data());
  TH1F *hZPeak=(TH1F*) file.Get(theNormHistoName.Data());
  if(hZPeak==0) {
    cout << "[LumiNorm] *** Error: histo " << theNormHistoName << " not valid!" << endl;
    throw std::exception();
  }
  nData+=hZPeak->Integral();
  cout << "# ev. in data is " << nData << endl;
}


// Add MC sample to the count of events to be used for data/MC normalization
// the histo used to get the count is zmass_restricted_1
// the default file for this histo is: theInputDir+"/histos_zzanalysis_mc_"+sampleName+"_"+theFinalState+".root"
void LumiNorm::addMC(const TString& sampleName) {
  TString thisFileName=theLumiFileName;
  thisFileName.ReplaceAll("whatsample", conversionMap[sampleName]);
  TString fileName=theInputDir+"/"+thisFileName;
  addMC(sampleName, fileName);
}


void LumiNorm::addMC(const TString& sampleName, const TString& filename) {
  FileStat_t buffer;
  if(gSystem->GetPathInfo(filename.Data(), buffer) != 0) {
    cout << "[LumiNorm] *** Error: Input file " <<  filename << " doesn't exist!" << endl;
    throw std::exception();
  }
  TFile file(filename.Data());
  TH1F *hZPeak=(TH1F*) file.Get(theNormHistoName.Data());

  // Get n. events, cross section, branching ratio, compute scale factor
  TH1F *hCount=(TH1F*) file.Get(theLumiHistoName.Data());
  initNEvMap[sampleName]=hCount->GetBinContent(1);
  xSectMap[sampleName]=hCount->GetBinContent(6);
  brRatioMap[sampleName]=hCount->GetBinContent(7);
  if(theLuminosity==0.) theLuminosity=hCount->GetBinContent(5);
  if(theLuminosity==0.) {
    cout << "[LumiNorm] *** Error: luminosity is 0!" << endl;
    throw std::exception();
  }

  double scale=additionalScale;
  if(rescaleLuminosity) 
    scale*=theLuminosity;
  if(rescaleCrossSect) 
    scale*=(xSectMap[sampleName]*brRatioMap[sampleName]/initNEvMap[sampleName]);

  double unscIntg=hZPeak->Integral();
  if(scale!=1.) hZPeak->Scale(scale);
  cout << "# ev. in sample " << sampleName << " is: " 
       << unscIntg << " (before rescaling) and " 
       << hZPeak->Integral() << " (after rescaling)" << endl;
  cout << "  # of entries: " << hZPeak->GetEntries() << endl;
  nMC+=hZPeak->Integral();
}


void  LumiNorm::addMC(const vector<TString>& sampleNames) {
  for(vector<TString>::const_iterator name=sampleNames.begin(); name!=sampleNames.end();
      ++name) {
    addMC(*name);
  }
}



double LumiNorm::getNormalizationFactor() const {
  double ret=1;
  if(normalizeToZ==true && nData!=0 && nMC!=0) {
    ret=nData/nMC;
  }
  //cout << "Normalization to Z peak factor: " << ret << endl;
  return ret;
}


double LumiNorm::getLuminosity() const {
  return theLuminosity;
}


float LumiNorm::getInitialNEv(const TString& sampleName) const {
  std::map<TString, float>::const_iterator it=initNEvMap.find(sampleName);
  return (it==initNEvMap.end() ? 0 : it->second);
}


float  LumiNorm::getCrossSection(const TString& sampleName) const {
  std::map<TString, float>::const_iterator it=xSectMap.find(sampleName);
  return (it==xSectMap.end() ? 0 : it->second);
}


float  LumiNorm::getBranchingRatio(const TString& sampleName) const {
  std::map<TString, float>::const_iterator it=brRatioMap.find(sampleName);
  return (it==brRatioMap.end() ? 0 : it->second);
}


float LumiNorm::getScaleFactor(const TString& sampleName) const {
  if(sampleName=="data")
    return 1.;
  if(initNEvMap.find(sampleName)==initNEvMap.end()) {
    cout << "[LumiNorm] *** Warning: sample " << sampleName.Data() 
	 << " not added yet!" << endl;
    return 1.;
  }
  double scale=additionalScale*getNormalizationFactor();
  if(rescaleLuminosity) 
    scale*=theLuminosity;
  if(rescaleCrossSect) 
    scale*=(xSectMap.find(sampleName)->second*brRatioMap.find(sampleName)->second/initNEvMap.find(sampleName)->second);
  return scale;
}


void LumiNorm::setAdditionalScale(double scale) {
  additionalScale=scale;
}


void LumiNorm::setNormalizHistoName(const TString& histoname) {
  theNormHistoName=histoname;
}


void LumiNorm::normalizeToZPeak(bool doNormalize) {
  normalizeToZ=doNormalize;
}
