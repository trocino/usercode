
/*
 *  See header file for a description of this class.
 *
 *  $Date: 2011/06/06 16:21:32 $
 *  $Revision: 1.1 $
 *  \author G. Cerminara - CERN
            D. Trocino - Northeastern University
 */

#include "LikelihoodDisc.h"
#include "TTree.h"
#include "TString.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "TMinuit.h"

#include <iostream>
#include "TCanvas.h"

using namespace std;


// Build the discriminat from signal and backgorund trees
LikelihoodDisc::LikelihoodDisc(TTree *sigTree, TTree *bckTree) : theSigTree(sigTree),
								 theBckTree(bckTree),
								 massMin(-9999),
								 massMax(9999),
								 doSeparateSample(true),
								 verbosityLevel(-1),
								 theChargeCut(0),
								 doMassCut(false) {}

LikelihoodDisc::~LikelihoodDisc() {}

// Add a variable to the discriminat
void LikelihoodDisc::addVariable(const TString& lName, const TString& varName, const TString& title,
				 const int nBins, double min, double max,
				 const TString& fitFuncLR) {
  TF1 *LRatioFit=new TF1(lName+"Fit", fitFuncLR.Data(), min, max);
  addVariable(lName, varName, title, nBins, min, max, LRatioFit);
}

void LikelihoodDisc::addVariable(const TString& lName, const TString& varName, const TString& title,
				 const int nBins, double min, double max, const double* bins, 
				 const TString& fitFuncLR) {
  TF1 *LRatioFit=new TF1(lName+"Fit", fitFuncLR.Data(), min, max);
  addVariable(lName, varName, title, nBins, min, max, bins, LRatioFit);
}

void LikelihoodDisc::addVariable(const TString& lName, const TString& varName, const TString& title,
				 const int nBins, double min, double max, TF1* fitFunc) {
  double *bins=new double[nBins+1];
  const double binSize=(max-min)/nBins;
  for(int i=0; i!=nBins+1; ++i) {
    bins[i]=min + i*binSize;
  }
  addVariable(lName, varName, title, nBins, min, max, bins,fitFunc);
  
}

void LikelihoodDisc::addVariable(const TString& lName, const TString& varName, const TString& title,
				 const int nBins, double min, double max,
				 const double* bins, TF1* fitFunc) {

  cout << "----------------- Variable: " << varName << endl;
  theVariableMap[lName]=varName;
  // Set the variable to be read from the trees
  float event, weight, mass, variable;
  theSigTree->SetBranchAddress("event", &event);
  theSigTree->SetBranchAddress("weight", &weight);
  theSigTree->SetBranchAddress(varName.Data(), &variable);
  if(doMassCut && varName!="diLeptInvMass") {
    theSigTree->SetBranchAddress("diLeptInvMass", &mass);
  }

  // Read the charge of the leptons only if needed
  float lept1Charge, lept2Charge;
  if(theChargeCut!=0) {
    theSigTree->SetBranchAddress("leadCharge", &lept1Charge);
    theSigTree->SetBranchAddress("subleadCharge", &lept2Charge);
  }

  // Build the relative histos for signal and bck
  TH1D *hSig=new TH1D(("h"+varName+"Sig").Data(), title.Data(), nBins, bins);
  hSig->Sumw2();
  TH1D *hBck=new TH1D(("h"+varName+"Bck").Data(), title.Data(), nBins, bins);
  hBck->Sumw2();

  // Loop over signal events
  for(int j=0; j!=theSigTree->GetEntries(); j++) {
    theSigTree->GetEntry(j);

    // Use odd MC events for the training and even events for the testing
    if(doSeparateSample && ((int)event+1)%2) continue;
    if(varName=="diLeptInvMass") mass=variable;

    if(doMassCut && (mass<massMin || mass>massMax)) continue;
    if(theChargeCut!=0 && lept1Charge*lept2Charge!=theChargeCut) continue;

    hSig->Fill(variable, weight);
  }

  // Normalize the histo to 1
  normalize(hSig);

  theBckTree->SetBranchAddress("weight", &weight);
  theBckTree->SetBranchAddress(varName.Data(), &variable);
  if(doMassCut && varName!="diLeptInvMass") {
    theBckTree->SetBranchAddress("diLeptInvMass", &mass);
  }

  // Read the charge of the leptons only if needed
  if(theChargeCut!=0) {
    theBckTree->SetBranchAddress("leadCharge", &lept1Charge);
    theBckTree->SetBranchAddress("subleadCharge", &lept2Charge);
  }


  // Loop over background events
  for(int j=0; j!=theBckTree->GetEntries(); ++j) {
    theBckTree->GetEntry(j);
    if(varName=="diLeptInvMass") mass=variable;

    if(doMassCut && (mass<massMin || mass>massMax)) continue;
    if(theChargeCut!=0 && lept1Charge*lept2Charge!=theChargeCut) continue;

    hBck->Fill(variable, weight);
  }

  // Normalize the histo to 1
  normalize(hBck);

  // Store the pdf histos
  thePdfMap[lName]=make_pair(hSig, hBck);
  // Store the likelihood ratio for this variable
  theLRMap[lName]=buildLikelihoodRatio(hSig, hBck, fitFunc->GetName(), fitFunc, min, max);

  theSigTree->ResetBranchAddresses();
  theBckTree->ResetBranchAddresses();
}

void LikelihoodDisc::normalize(TH1D *histo) {
  histo->Scale(1./histo->Integral());
}

pair<TH1D *, TF1 *> LikelihoodDisc::buildLikelihoodRatio(const TH1D *hSig, const TH1D *hBck,
							 TString fitName, TF1 *LRatioFit,
							 double funcMin, double funcMax) {
  TString ratioName=TString(hSig->GetName())+"_ratio";
  TH1D *LRatio=(TH1D*)hSig->Clone(ratioName.Data());
  LRatio->Divide(hBck);
  LRatio->Fit(fitName.Data(), "");
  return make_pair(LRatio,LRatioFit);
}

double LikelihoodDisc::evaluate(const vector<TString> &likelihoodComposition,
				const vector<float> &variables) const {
  if(likelihoodComposition.size()!=variables.size()) {
    cerr << "[LikelihoodDisc] *** Error: size of likelihood composition "
	 << "and variables don't match!" << endl;
    cerr << "      Likelihoods: ";
    for(vector<TString>::const_iterator llname=likelihoodComposition.begin();
	llname!=likelihoodComposition.end(); ++llname) {
      cerr << *llname << " ";
    }
    cerr << endl;
    return -1;
  }
  double ret=1;
  for(unsigned int i=0; i!=likelihoodComposition.size(); ++i) {
    map<TString, std::pair<TH1D *, TF1 *> >::const_iterator lRatio=theLRMap.find(likelihoodComposition[i]);
    if(lRatio!=theLRMap.end()) {
      ret*=((*lRatio).second.second)->Eval(variables[i]);
    } else {
      cerr << "[LikelihoodDisc] *** Error: Likelihood " << likelihoodComposition[i]
	   << " not trained!" << endl;
      return -1;
    }
  }
  return transformLikelihoodOutput(ret);
}

double LikelihoodDisc::transformLikelihoodOutput(double r) const {
  const double fEpsilon=1e-8;
  double ret=r/(r+1);
  if(ret<=0.0) ret=fEpsilon;
  else if(ret >= 1.0) ret=1.0-fEpsilon;

  return ret;
}

pair<TH1D *, TH1D *>  LikelihoodDisc::getPdf(const TString& lName) const {
  map<TString, std::pair<TH1D *, TH1D *> >::const_iterator pdfPair=thePdfMap.find(lName);
  if(pdfPair!=thePdfMap.end()) {
    return (*pdfPair).second;
  } else {
    cerr << "[LikelihoodDisc] *** Error: No pfd for likelihood " << lName << endl;
    return make_pair((TH1D*)0,(TH1D*)0);
  }
}
  
pair<TH1D *, TF1 *>  LikelihoodDisc::getLR(const TString& lName) const {
  map<TString, std::pair<TH1D *, TF1 *> >::const_iterator lRatio=theLRMap.find(lName);
  if(lRatio!=theLRMap.end()) {
    return (*lRatio).second;
  } else {
    cerr << "[LikelihoodDisc] *** Error: No LR for likelihood " << lName << endl;
    return make_pair((TH1D*)0,(TF1*)0);
  }
}
  
void LikelihoodDisc::setMassCut(double min, double max) {
  massMin=min;
  massMax=max;
}

TH1F *LikelihoodDisc::evaluate(const vector<TString>& likelihoodComposition,
				TTree *tree, 
				const TString& sampleName, 
				const int nBins) const {
  double *bins=new double[nBins+1];
  double max=1.;
  double min=0.;
  const double binSize=(max-min)/nBins;
  for(int i=0; i!=nBins+1; ++i) {
    bins[i]=min + i*binSize;
  }
  return evaluate(likelihoodComposition, tree, sampleName, nBins, bins);
}
  
TH1F *LikelihoodDisc::evaluate(const std::vector<TString>& likelihoodComposition,
			       TTree *tree, 
			       const TString& sampleName,
			       const int nBins, 
			       const double *bins) const {
  // Check that the variable binning is valid, otherwise run with constant binning
  if(bins==0) return evaluate(likelihoodComposition, tree, sampleName, nBins);

  // This is where we read the tree branches
  vector<float> theReadValues(likelihoodComposition.size(), -1);
  vector<float> dumpAddress(likelihoodComposition.size(), -1);

  /*
  // Switch for the kinematic fit machinery
  bool evalKinemFit=false;
  int kinemFitPlace=-1;  
  */

  TString hName="hTest_";
  float mass;
  int isMassInVar=-99;

  // Set the branch addresses
  for(vector<TString>::const_iterator likelih=likelihoodComposition.begin();
      likelih!=likelihoodComposition.end(); ++likelih) {
    TString varname=theVariableMap.find(*likelih)->second;
    if(varname=="diLeptInvMass") {
      isMassInVar=likelih-likelihoodComposition.begin();
    }

    tree->SetBranchAddress(varname.Data(), &theReadValues[likelih - likelihoodComposition.begin()]);
    /*
    if(varname!="Chi2") {
      tree->SetBranchAddress(varname.Data(), &theReadValues[likelih - likelihoodComposition.begin()]);
    } else {
      evalKinemFit=true;
      kinemFitPlace=likelih - likelihoodComposition.begin();
    }
    */

    if(likelih!=likelihoodComposition.begin()) {
      hName+="x";
    }
    hName+=*likelih;
  }
  float weight=1;
  tree->SetBranchAddress("weight", &weight);
  float event;
  tree->SetBranchAddress("event", &event);
  if(doMassCut && isMassInVar==-99) {
    tree->SetBranchAddress("diLeptInvMass", &mass);
  }
  
  // Read the charge of the leptons only if needed
  float lept1Charge, lept2Charge;
  if(theChargeCut!=0) {
    tree->SetBranchAddress("leadCharge", &lept1Charge);
    tree->SetBranchAddress("subleadCharge", &lept2Charge);
  }

  /*
  // Read the branches needed for the kinematic fit
  float Lept1Px, Lept1Py, Lept1Pz, Lept1E, dPt1,
    Lept2Px, Lept2Py, Lept2Pz, Lept2E, dPt2;
  if(evalKinemFit) {
    tree->SetBranchAddress("Lept1Px", &Lept1Px);
    tree->SetBranchAddress("Lept1Py", &Lept1Py);
    tree->SetBranchAddress("Lept1Pz", &Lept1Pz);
    tree->SetBranchAddress("Lept1E", &Lept1E);
    //tree->SetBranchAddress("leadCharge", &leadCharge);
    tree->SetBranchAddress("dPt1", &dPt1);
    tree->SetBranchAddress("Lept2Px", &Lept2Px);
    tree->SetBranchAddress("Lept2Py", &Lept2Py);
    tree->SetBranchAddress("Lept2Pz", &Lept2Pz);
    tree->SetBranchAddress("Lept2E", &Lept2E);
    //tree->SetBranchAddress("subleadCharge", &subleadCharge);
    tree->SetBranchAddress("dPt2", &dPt2);
  }
  */

  hName+=sampleName;
  TH1D hCount("h", "h", 1, 0.5, 1.5);
  TH1F *ret=new TH1F(hName.Data(), "Test statistics", nBins, bins);
  ret->Sumw2();
  for(unsigned int i=0; i!=tree->GetEntries(); ++i) {
    tree->GetEntry(i);

    if(isMassInVar!=-99) {
      mass=theReadValues[isMassInVar];
    }
    if(doMassCut && (mass<massMin || mass>massMax)) {
      continue;
    }

    if(theChargeCut!=0 && lept1Charge*lept2Charge!=theChargeCut) continue;

    /*
    // Compute the chi2 from the kineamtic fit
    if(evalKinemFit) {
      TLorentzVector lept1(Lept1Px, Lept1Py, Lept1Pz, Lept1E);
      TLorentzVector lept2(Lept2Px, Lept2Py, Lept2Pz, Lept2E);
      
      double chi2=computeChi2KinemFitV2(lept1, dPt1, lept2, dPt2);
      theReadValues[kinemFitPlace]=chi2;
    }
    */

    hCount.Fill(1, weight);
    // Use odd MC events for the training and even events for the testing
    if(doSeparateSample && sampleName!="data" && ((int)event)%2) continue;
    ret->Fill(evaluate(likelihoodComposition, theReadValues), weight);
  }
  if(doSeparateSample && sampleName!="data") {
    // Rescale the integral to match the expected number of events
    double integral_init=hCount.Integral();
    if(integral_init!=0) {
      double integral_new=ret->Integral();
      if(integral_new!=0) {
	ret->Scale(integral_init/integral_new);
      } else {
	cout << "[LikelihoodDisc] *** Warning: the evaluation histo "
	     << "has 0 entries but the training one doesn't!"
	     << endl; 
      }
    }
  }

  tree->ResetBranchAddresses();
  return ret;
}

void LikelihoodDisc::separateSamples(bool doSeparate) {
  doSeparateSample=doSeparate;
}


// Evaluate the Likelihood on a sample: an histogram (Test Statistic) is produced
// a cut can be applied on the likelihood.
// NOTE: for samples named "data", the weight branch is not read from the tree
TH1F *LikelihoodDisc::evaluateData(const vector<TString>& likelihoodComposition,
				   TTree *tree, 
				   const TString& sampleName, 
				   double llCut,
				   const int nBins) const {

  // This is where we read the tree branches
  vector<float> theReadValues(likelihoodComposition.size(), -1);
  vector<float> dumpAddress(likelihoodComposition.size(), -1);

  TString hName="hTest_";
  float mass;

  int isMassInVar=-99;
  // Set the branch addresses
  for(vector<TString>::const_iterator likelih=likelihoodComposition.begin();
      likelih!=likelihoodComposition.end(); ++likelih) {
    TString varname=theVariableMap.find(*likelih)->second;
    if(varname=="diLeptInvMass") {
      isMassInVar=likelih-likelihoodComposition.begin();
    }

    tree->SetBranchAddress(varname.Data(), &theReadValues[likelih-likelihoodComposition.begin()]);

    if(likelih!=likelihoodComposition.begin()) {
      hName+="x";
    }
    hName+=(*likelih);
  }

  // Read the charge of the leptons only if needed
  float lept1Charge, lept2Charge;
  if(theChargeCut!=0) {
    tree->SetBranchAddress("leadCharge", &lept1Charge);
    tree->SetBranchAddress("subleadCharge", &lept2Charge);
  }

  float weight=1;
  if(sampleName!="data")
    tree->SetBranchAddress("weight", &weight);
  float event;
  tree->SetBranchAddress("event",&event);
  if(doMassCut && isMassInVar==-99) {
    tree->SetBranchAddress("diLeptInvMass", &mass);
  }

  hName+=sampleName;
  TH1D hCount("h", "h", 1, 0.5, 1.5);
  TH1F *ret=new TH1F(hName.Data(), "Test statistics", nBins, 0, 1);
  ret->Sumw2();
  for(unsigned int i=0; i!=tree->GetEntries(); ++i) {
    tree->GetEntry(i);
    if(isMassInVar!=-99) {
      mass=theReadValues[isMassInVar];
    }

    if(doMassCut && (mass<massMin || mass>massMax)) {
      continue;
    }
    if(theChargeCut!=0 && lept1Charge*lept2Charge!=theChargeCut) continue;
    if(evaluate(likelihoodComposition, theReadValues)<llCut) continue;

    hCount.Fill(1,weight);
    ret->Fill(evaluate(likelihoodComposition, theReadValues), weight);
  }
  if(doSeparateSample && sampleName!="data") {
    // Rescale the integral to match the expected number of events
    double integral_init=hCount.Integral();
    if(integral_init!=0) {
      double integral_new=ret->Integral();
      if(integral_new!=0) {
	ret->Scale(integral_init/integral_new);
      } else {
	cout << "[LikelihoodDisc] *** Warning: the evaluation histo has 0 entries "
	     << "but the training one doesn't!" << endl; 
      }
    }
  }

  tree->ResetBranchAddresses();
  return ret;
}

void LikelihoodDisc::setVerbosityLevel(int level) {
  verbosityLevel=level;
}

// Switch for the cut on the lepton charge
//  0 -> no cut
// -1 -> opposite charge
//  1 -> same charge
void LikelihoodDisc::setChargeCut(int cut) {
  theChargeCut=cut;
}


void LikelihoodDisc::useMassCut(bool doit) {
  doMassCut=doit;
}


//
// Methods for chi2 and kinematic fit
// Will not be used for the time being
//
/*
void LikelihoodDisc::addChi2Variable(const TString& lName, const TString& title,
				     const int nBins, double min, double max, const double* bins, 
				     const TString& fitFuncLR, bool useProbability) {
  useChi2Probability=useProbability;
  TF1 * fitFunc=new TF1(lName+"Fit", fitFuncLR.Data(),min,max);
  cout << "----------------- Variable: Chi2 from kinematic fit" << endl;
  // Map the likelihood name to this variable
  TString varName="Chi2";
  theVariableMap[lName]=varName;
  
  // set the variable to be read from the trees
  float Lept1Px, Lept1Py, Lept1Pz, Lept1E, leadCharge,
    Lept2Px, Lept2Py, Lept2Pz, Lept2E, subleadCharge,
    dPt1,dPt2,
    weight, mass;
  float nSMT1, nSMT2;
  float event;
  // ---------------------------------------------------------------------
  // --- Signal tree
  // set the branches
  theSigTree->SetBranchAddress("Lept1Px", &Lept1Px);
  theSigTree->SetBranchAddress("Lept1Py", &Lept1Py);
  theSigTree->SetBranchAddress("Lept1Pz", &Lept1Pz);
  theSigTree->SetBranchAddress("Lept1E", &Lept1E);
  theSigTree->SetBranchAddress("leadCharge", &leadCharge);
  theSigTree->SetBranchAddress("dPt1", &dPt1);
  theSigTree->SetBranchAddress("Lept2Px", &Lept2Px);
  theSigTree->SetBranchAddress("Lept2Py", &Lept2Py);
  theSigTree->SetBranchAddress("Lept2Pz", &Lept2Pz);
  theSigTree->SetBranchAddress("Lept2E", &Lept2E);
  theSigTree->SetBranchAddress("subleadCharge", &subleadCharge);
  theSigTree->SetBranchAddress("dPt2", &dPt2);
  theSigTree->SetBranchAddress("weight", &weight);
  theSigTree->SetBranchAddress("diLeptInvMass", &mass);
  theSigTree->SetBranchAddress("event", &event);
  // Read the # of SMT hits only if needed
  if(nSMTHitsSwitch!=0) {
    theSigTree->SetBranchAddress("nSMT1", &nSMT1);
    theSigTree->SetBranchAddress("nSMT2", &nSMT2);
  }

  float alatCorrReducedT;
  // Read the alat value only if needed
  if(doAlatCut) {
    theSigTree->SetBranchAddress("alatCorrReduceT", &alatCorrReducedT);
  }

  float lept1Charge, lept2Charge;
  // Read the charge of the leptons only if needed
  if(theChargeCut!=0) {
    theSigTree->SetBranchAddress("leadCharge", &lept1Charge);
    theSigTree->SetBranchAddress("subleadCharge", &lept2Charge);
  }



  // Build the relative histos for signal and bck
  TH1D *hSig=new TH1D(("h"+varName+"Sig").Data(),title.Data(),nBins,bins);
  hSig->Sumw2();
  TH1D *hBck=new TH1D(("h"+varName+"Bck").Data(),title.Data(),nBins,bins);
  hBck->Sumw2();

  // loop over signal events
  for(int j=0; j!=theSigTree->GetEntries();j++) {
    theSigTree->GetEntry(j);
    // Use odd MC events for the training and even events for the testing
    if(doSeparateSample && ((int)event+1)%2) continue;
    if(mass < massMin || mass > massMax) continue;
    
    // Require at least 1 of the two tracks to have SMT hits
    if(nSMTHitsSwitch==1 && (nSMT1==0 && nSMT2==0)) continue;

    // Require both tracks to have SMT hits
    if(nSMTHitsSwitch==2 && (nSMT1==0 || nSMT2==0)) continue;
    
    // Require both tracks to have no SMT hits
    if(nSMTHitsSwitch==3 && (nSMT1!=0 || nSMT2!=0)) continue;

    if(doAlatCut && alatCorrReducedT < alatMin) continue;

    if(theChargeCut!=0 && lept1Charge*lept2Charge!=theChargeCut) continue;

    TLorentzVector lept1(Lept1Px, Lept1Py, Lept1Pz, Lept1E);
    TLorentzVector lept2(Lept2Px, Lept2Py, Lept2Pz, Lept2E);

    double variable=computeChi2KinemFitV2(lept1, dPt1, lept2, dPt2);
    hSig->Fill(variable, weight);
  }
  // normalize the histo to 1
  normalize(hSig);
  
  // ---------------------------------------------------------------------
  // --- Background tree
  theBckTree->SetBranchAddress("Lept1Px", &Lept1Px);
  theBckTree->SetBranchAddress("Lept1Py", &Lept1Py);
  theBckTree->SetBranchAddress("Lept1Pz", &Lept1Pz);
  theBckTree->SetBranchAddress("Lept1E", &Lept1E);
  theBckTree->SetBranchAddress("leadCharge", &leadCharge);
  theBckTree->SetBranchAddress("dPt1", &dPt1);
  theBckTree->SetBranchAddress("Lept2Px", &Lept2Px);
  theBckTree->SetBranchAddress("Lept2Py", &Lept2Py);
  theBckTree->SetBranchAddress("Lept2Pz", &Lept2Pz);
  theBckTree->SetBranchAddress("Lept2E", &Lept2E);
  theBckTree->SetBranchAddress("subleadCharge", &subleadCharge);
  theBckTree->SetBranchAddress("dPt2", &dPt2);
  theBckTree->SetBranchAddress("weight", &weight);
  theBckTree->SetBranchAddress("diLeptInvMass", &mass);
  // Read the # of SMT hits only if needed
  if(nSMTHitsSwitch!=0) {
    theBckTree->SetBranchAddress("nSMT1", &nSMT1);
    theBckTree->SetBranchAddress("nSMT2", &nSMT2);
  }


  // Read the alat value only if needed
  if(doAlatCut) {
    theBckTree->SetBranchAddress("alatCorrReduceT", &alatCorrReducedT);
  }

    // Read the charge of the leptons only if needed
  if(theChargeCut!=0) {
    theBckTree->SetBranchAddress("leadCharge", &lept1Charge);
    theBckTree->SetBranchAddress("subleadCharge", &lept2Charge);
  }




  // loop over bck events
  for(int j=0; j!=theBckTree->GetEntries();j++) {
    theBckTree->GetEntry(j);
    if(mass < massMin || mass > massMax) continue;

    // Require at least 1 of the two tracks to have SMT hits
    if(nSMTHitsSwitch==1 && (nSMT1==0 && nSMT2==0)) continue;

    // Require both tracks to have SMT hits
    if(nSMTHitsSwitch==2 && (nSMT1==0 || nSMT2==0)) continue;
    
    // Require both tracks to have no SMT hits
    if(nSMTHitsSwitch==3 && (nSMT1!=0 || nSMT2!=0)) continue;

    if(doAlatCut && alatCorrReducedT < alatMin) continue;

    if(theChargeCut!=0 && lept1Charge*lept2Charge!=theChargeCut) continue;

    TLorentzVector lept1(Lept1Px, Lept1Py, Lept1Pz, Lept1E);
    TLorentzVector lept2(Lept2Px, Lept2Py, Lept2Pz, Lept2E);

//     cout << "Lept1: " << lept1.Px() << " " << lept1.Py() << " " << lept1.Pz() << " " << lept1.E() << endl;
//     cout << "Lept2: " << lept2.Px() << " " << lept2.Py() << " " << lept2.Pz() << " " << lept2.E() << endl;

    double variable=computeChi2KinemFitV2(lept1, dPt1, lept2, dPt2);
    hBck->Fill(variable, weight);
  }

  // normalize the histo to 1
  normalize(hBck);


  // Store the pdf histos
  thePdfMap[lName]=make_pair(hSig, hBck);
  // Store the likelihood ratio for this variable
  theLRMap[lName]=buildLikelihoodRatio(hSig, hBck, fitFunc->GetName(),fitFunc,min,max);

  theSigTree->ResetBranchAddresses();
  theBckTree->ResetBranchAddresses();
  

}
*/

/*
double LikelihoodDisc::computeChi2KinemFit(const TLorentzVector& lept1, float dOneOverPt1,
					   const TLorentzVector& lept2, float dOneOverPt2) const {
  cout << " --- computeChi2KinemFit --------------------------------------------------" << endl;
  // Compute the error on 1/p starting from the error on 1/pt
  double sigma_1OverP1=dOneOverPt1*(1/sqrt(1+(1/TMath::Power(TMath::Tan(lept1.Theta()),2))));
  double sigma_1OverP2=dOneOverPt2*(1/sqrt(1+(1/TMath::Power(TMath::Tan(lept2.Theta()),2))));


//   cout << " Delta 1/P (1): " << dOneOverPt1 << endl;
//   cout << " Delta 1/P (2): " << dOneOverPt1 << endl;


  // Create the MINUIT object
  TMinuit minuit(6);

  // Set the MINUIT verbosity level
  minuit.SetPrintLevel(verbosityLevel);

  // Pass to MINUIT the pointer to the function to be minimized
  minuit.SetFCN(chi2function);

  double fval=0; //The result of the function that you want to minimize
  int nPar=6; //Number of variable parameters
  double *grad=0; //Gradient, not used
  


  double parameters[6];//Parameter array
  parameters[0]=1./lept1.P();
  parameters[1]=lept1.Angle(lept2.Vect());
  parameters[2]=1./lept1.P();
  parameters[3]=1./lept2.P();
  parameters[4]=sigma_1OverP1;
  parameters[5]=sigma_1OverP2;

//   double parameters[6];//Parameter array
//   parameters[0]=1./lept1.Pt();
//   parameters[1]=lept1.Angle(lept2.Vect());
//   parameters[2]=1./lept1.Pt();
//   parameters[3]=1./lept2.Pt();
//   parameters[4]=dOneOverPt1;
//   parameters[5]=dOneOverPt2;



  int iflag=1; //Switch to execute different operations inside the function (useless)

  //Evaluate the function
  minuit.Eval(nPar, grad, fval, parameters, iflag);
  cout << " chi2 init.=" << fval << endl;


  // Assign a value to a parameter
  // DefineParameter(Int_t parNo, const char *name, Double_t initVal, Double_t initErr, Double_t lowerLimit, Double_t upperLimit ) 
  minuit.DefineParameter(0, "1/P_lead (fit)", parameters[0], 1., 0, 0);
  minuit.DefineParameter(1, "angle",  parameters[1], 0, 0, 0);
  minuit.DefineParameter(2, "1/P_lead (obs)",  parameters[2], 0, 0, 0);
  minuit.DefineParameter(3, "1/P_trail (obs)",  parameters[3], 0, 0, 0);
  minuit.DefineParameter(4, "sigma(1/P_lead) (obs)",  parameters[4], 0, 0, 0);
  minuit.DefineParameter(5, "sigma(1/P_trail) (obs)",  parameters[5], 0, 0, 0);

  // Release a parameter (It becomes a variable parameter)
  minuit.Release(0);

  // Fix other parameters (No effect (IGNORED) on const parameters)
  minuit.FixParameter(1);
  minuit.FixParameter(2);
  minuit.FixParameter(3);
  minuit.FixParameter(4);
  minuit.FixParameter(5);


  if(verbosityLevel > 0) {
    //  returns the total number of parameters that have been defined.
    //  (fixed and free)
    cout << "    Number of Par: " <<  minuit.GetNumPars() << endl;
    //Returns the number of currently fixed parameters
    cout << "    Number of Fixed Par: " <<  minuit.GetNumFixedPars() << endl;
    // returns the number of currently free parameters
    cout << "    Number of Free Par: " << minuit.GetNumFreePars() << endl;
  }

  if(verbosityLevel > -1)
    cout << " MIGRAD ============================================================================" << endl;
  // Call MIGRAD routine
  int migradErrorCode=minuit.Migrad();
  cout << " MIGRAD error code: " << migradErrorCode << endl;
  
  cout << "===== Minimization results:" << endl;
  double fitParam=-1;
  double fitParamErr=-1;
  // Get parameter values and errors
  minuit.GetParameter(0, fitParam, fitParamErr);

  parameters[0]=fitParam;
  minuit.Eval(nPar, grad, fval, parameters, iflag);


  cout << "== Par0: " << fitParam << " +/- " << fitParamErr << endl;
  cout << "parameters[0]=" <<  parameters[0] << endl;
  cout << "== Chi2 final: " << fval << endl;

  double ret=fval;
  if(useChi2Probability) ret=TMath::Prob(fval,1);
  return ret;
}
*/

/*
void LikelihoodDisc::addChi2Variable(const TString& lName, const TString& title,
				     const int nBins, double min, double max, 
				     const TString& fitFuncLR, bool useProbability) {
  double *bins=new double[nBins+1];
  const double binSize=(max - min)/nBins;
  for(int i=0; i!=nBins+1; ++i) {
    bins[i]=min + i*binSize;
  }
  addChi2Variable(lName, title, nBins, min, max, bins,fitFuncLR,useProbability);
}
*/

/*
// The function to be minimized
void chi2function(int &npar, double *grad, double &fval, double *par, int flag) {
  double oneOverPt1_fit=par[0];
  double angle=par[1];
  double oneOverPt1_obs=par[2];
  double oneOverPt2_obs=par[3];
  double sigma_oneOverPt1=par[4];
  double sigma_oneOverPt2=par[5];

//   double massZ=91.1876;
  double massZ=89.0;
  double oneOverPt2_fit=2.*(1-TMath::Cos(angle))/(oneOverPt1_fit*massZ*massZ);


  double term1=(oneOverPt1_obs - oneOverPt1_fit)/sigma_oneOverPt1;
  double term2=(oneOverPt2_obs - oneOverPt2_fit)/sigma_oneOverPt2;

  fval=term1*term1 + term2*term2;
}
*/

/*
double LikelihoodDisc::computeChi2KinemFitV2(const TLorentzVector& lept1, float dOneOverPt1,
					   const TLorentzVector& lept2, float dOneOverPt2) const {
  if(verbosityLevel > -1)
    cout << " --- computeChi2KinemFit Version 2--------------------------------------------------" << endl;
  // Compute the error on 1/p starting from the error on 1/pt

  // Create the MINUIT object
  TMinuit minuit(11);

  // Set the MINUIT verbosity level
  minuit.SetPrintLevel(verbosityLevel);

  // Pass to MINUIT the pointer to the function to be minimized
  minuit.SetFCN(chi2functionV2);

  double fval=0; //The result of the function that you want to minimize
  int nPar=11; //Number of variable parameters
  double *grad=0; //Gradient, not used
  
  double parameters[11];//Parameter array
  parameters[0]=lept1.Px();
  parameters[1]=lept1.Py();
  parameters[2]=lept1.Pz();
  parameters[3]=lept1.E();
  parameters[4]=dOneOverPt1;
  parameters[5]=lept2.Px();
  parameters[6]=lept2.Py();
  parameters[7]=lept2.Pz();
  parameters[8]=lept2.E();
  parameters[9]=dOneOverPt2;
  parameters[10]=lept1.Pt(); // this is the parameter we want to fit




  int iflag=1; //Switch to execute different operations inside the function (useless)

  //Evaluate the function
  minuit.Eval(nPar, grad, fval, parameters, iflag);
  if(verbosityLevel > -1)
    cout << " chi2 init.=" << fval << endl;


  // Assign a value to a parameter
  // DefineParameter(Int_t parNo, const char *name, Double_t initVal, Double_t initErr, Double_t lowerLimit, Double_t upperLimit ) 
  minuit.DefineParameter(0, "px_1 (obs)", parameters[0], 0, 0, 0);
  minuit.DefineParameter(1, "py_1 (obs)", parameters[1], 0, 0, 0);
  minuit.DefineParameter(2, "pz_1 (obs)", parameters[2], 0, 0, 0);
  minuit.DefineParameter(3, "E_1 (obs)", parameters[3], 0, 0, 0);
  minuit.DefineParameter(4, "sigma 1/pt_1", parameters[4], 0, 0, 0);
  minuit.DefineParameter(5, "px_2 (obs)", parameters[5], 0, 0, 0);
  minuit.DefineParameter(6, "py_2 (obs)", parameters[6], 0, 0, 0);
  minuit.DefineParameter(7, "pz_2 (obs)", parameters[7], 0, 0, 0);
  minuit.DefineParameter(8, "E_2 (obs)", parameters[8], 0, 0, 0);
  minuit.DefineParameter(9, "sigma 1/pt_2", parameters[9], 0, 0, 0);
  minuit.DefineParameter(10, "pt_1 (fit)", parameters[10], 1., 0, 0);


//   // Release a parameter (It becomes a variable parameter)
//   minuit.Release(0);

//   // Fix other parameters (No effect (IGNORED) on const parameters)
//   minuit.FixParameter(1);
//   minuit.FixParameter(2);
//   minuit.FixParameter(3);
//   minuit.FixParameter(4);
//   minuit.FixParameter(5);


  if(verbosityLevel > 0) {
    //  returns the total number of parameters that have been defined.
    //  (fixed and free)
    cout << "    Number of Par: " <<  minuit.GetNumPars() << endl;
    //Returns the number of currently fixed parameters
    cout << "    Number of Fixed Par: " <<  minuit.GetNumFixedPars() << endl;
    // returns the number of currently free parameters
    cout << "    Number of Free Par: " << minuit.GetNumFreePars() << endl;
  }
  
  if(verbosityLevel > -1)
    cout << "=== MIGRAD =========================================================================" << endl;
  // Call MIGRAD routine
  int migradErrorCode=minuit.Migrad();
  if(migradErrorCode!=0) {
    cout << "[LikelihoodDisc::computeChi2KinemFitV2] Error, MIGRAD error code: " << migradErrorCode << endl;
    return -9999;    
  }


  double fitParam=-1;
  double fitParamErr=-1;
  // Get parameter values and errors
  minuit.GetParameter(10, fitParam, fitParamErr);

  parameters[10]=fitParam;
  minuit.Eval(nPar, grad, fval, parameters, iflag);

  if(verbosityLevel > -1) {
    cout << "=Minimization results:" << endl;
    cout << "   1/Pt (lead): " << fitParam << " +/- " << fitParamErr << endl;
    //   cout << "parameters[10]=" <<  parameters[10] << endl;
    cout << "   Chi2 final: " << fval << endl;
  }
  double ret=fval;
  if(useChi2Probability) ret=TMath::Prob(fval,1);
  return ret;
}
*/

/*
// The function to be minimized
void chi2functionV2(int &npar, double *grad, double &fval, double *par, int flag) {
  double p1x_obs=par[0];
  double p1y_obs=par[1];
  double p1z_obs=par[2];
  double E1_obs=par[3];
  double dOneOverPt1=par[4];
  double p2x_obs=par[5];
  double p2y_obs=par[6];
  double p2z_obs=par[7];
  double E2_obs=par[8];
  double dOneOverPt2=par[9];
  double pt1_fit=par[10]; // the free parameter


  // 4-momenta of the observed leptons
  TLorentzVector lept1_obs(p1x_obs, p1y_obs, p1z_obs, E1_obs);
  TLorentzVector lept2_obs(p2x_obs, p2y_obs, p2z_obs, E2_obs);

  // The fitted 4-momenta: 
  // lept1: assume the same angle and fit only pt
  TLorentzVector lept1_fit=(pt1_fit/lept1_obs.Pt())*lept1_obs;
  // lept2: assume the same angle and p from mass constraint
  const double massZ=89.0;
  double angle=lept1_obs.Angle(lept2_obs.Vect());
  double p2_fit=massZ*massZ/(2.*(1-TMath::Cos(angle))*lept1_fit.P());
  TLorentzVector lept2_fit=(p2_fit/lept2_obs.P())*lept2_obs;

  

//   double massZ=91.1876;



  double term1=(1./lept1_obs.Pt() - 1./pt1_fit)/dOneOverPt1;
  double term2=(1./lept2_obs.Pt() - 1./lept2_fit.Pt())/dOneOverPt2;

  fval=term1*term1 + term2*term2;
}
*/

/*
TH1F *LikelihoodDisc::plotChi2Prob(TTree *tree, const TString& sampleName, const int nBins) const {
  // set the prefix of the histo name
  TString hName="hChi2Prob_";


  // Variables needed for the selection cuts
  float mass;
  tree->SetBranchAddress("diLeptInvMass", &mass);
  float weight=1;
  tree->SetBranchAddress("weight", &weight);
  float event;
  tree->SetBranchAddress("event",&event);
  float nSMT1, nSMT2;
  // Read the # of SMT hits only if needed
  if(nSMTHitsSwitch!=0) {
    tree->SetBranchAddress("nSMT1", &nSMT1);
    tree->SetBranchAddress("nSMT2", &nSMT2);
  }

  float alatCorrReducedT;
  // Read the # of SMT hits only if needed
  if(doAlatCut) {
    tree->SetBranchAddress("alatCorrReduceT", &alatCorrReducedT);
  }

  float lept1Charge, lept2Charge;
  // Read the charge of the leptons only if needed
  if(theChargeCut!=0) {
    tree->SetBranchAddress("leadCharge", &lept1Charge);
    tree->SetBranchAddress("subleadCharge", &lept2Charge);
  }



  // Read the branches needed for the kinematic fit
  float Lept1Px, Lept1Py, Lept1Pz, Lept1E,dPt1,
    Lept2Px, Lept2Py, Lept2Pz, Lept2E,dPt2;
  tree->SetBranchAddress("Lept1Px", &Lept1Px);
  tree->SetBranchAddress("Lept1Py", &Lept1Py);
  tree->SetBranchAddress("Lept1Pz", &Lept1Pz);
  tree->SetBranchAddress("Lept1E", &Lept1E);
  //     tree->SetBranchAddress("leadCharge", &leadCharge);
  tree->SetBranchAddress("dPt1", &dPt1);
  tree->SetBranchAddress("Lept2Px", &Lept2Px);
  tree->SetBranchAddress("Lept2Py", &Lept2Py);
  tree->SetBranchAddress("Lept2Pz", &Lept2Pz);
  tree->SetBranchAddress("Lept2E", &Lept2E);
  //     tree->SetBranchAddress("subleadCharge", &subleadCharge);
  tree->SetBranchAddress("dPt2", &dPt2);


  hName += sampleName;

  TH1D hCount("h","h",1,0.5,1.5);
  TH1F *ret=new TH1F(hName.Data(),"Chi2 probability of the kineamtic fit", nBins, 0,1);
  ret->Sumw2();
  for(unsigned int i=0; i!=tree->GetEntries(); ++i) {
    tree->GetEntry(i);

    // mass cut
    if(mass < massMin || mass > massMax) {
      continue;
    }     
    // Require at least 1 of the two tracks to have SMT hits
    if(nSMTHitsSwitch==1 && (nSMT1==0 && nSMT2==0)) continue;
    // Require both tracks to have SMT hits
    if(nSMTHitsSwitch==2 && (nSMT1==0 || nSMT2==0)) continue;
    // Require both tracks to have no SMT hits
    if(nSMTHitsSwitch==3 && (nSMT1!=0 || nSMT2!=0)) continue;
    if(doAlatCut && alatCorrReducedT < alatMin) continue;
    if(theChargeCut!=0 && lept1Charge*lept2Charge!=theChargeCut) continue;

    // Compute the chi2 from the kineamtic fit
    TLorentzVector lept1(Lept1Px, Lept1Py, Lept1Pz, Lept1E);
    TLorentzVector lept2(Lept2Px, Lept2Py, Lept2Pz, Lept2E);
    
    double chi2=computeChi2KinemFitV2(lept1, dPt1, lept2, dPt2);
   
    hCount.Fill(1,weight);
    // Use odd MC events for the training and even events for the testing
    if(doSeparateSample && sampleName!="data" && ((int)event)%2 ) continue;
    ret->Fill(chi2, weight);
  }

  if(doSeparateSample && sampleName!="data") {
    // Rescale the integral to match the expected # of events
    double integral_init=hCount.Integral();
    if(integral_init!=0) {
      double integral_new=ret->Integral();
      if(integral_new!=0) {
	ret->Scale(integral_init/integral_new);
      } else {
	cout << "[LikelihoodDisc] Warning: the evaluation histo has 0 entries but the training one doesn't!"
	     << endl; 
      }
    }
  }

  tree->ResetBranchAddresses();
  return ret;
}
*/
