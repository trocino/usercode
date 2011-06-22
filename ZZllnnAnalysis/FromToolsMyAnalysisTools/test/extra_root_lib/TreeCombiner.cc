
/*
 *  See header file for a description of this class.
 *
 *  $Date: 2008/03/13 17:44:26 $
 *  $Revision: 1.3 $
 *  \author G. Cerminara - CERN
 *          D. Trocino - Northeastern University
 */

#include "TreeCombiner.h"
#include "LumiNorm.h"


#include "TTree.h"
#include "TNtuple.h"
#include "TString.h"
#include "TFile.h"

#include <iostream>

using namespace std;

TreeCombiner::TreeCombiner(const LumiNorm* lumiNorm, 
			   const TString& inputDir, 
			   const TString& fileName) 
  : theLumiNorm(lumiNorm),
    theFileName(fileName),
    theInputDir(inputDir),
    massMinCut(0.),
    massMaxCut(99999),
    fullSelectionReq(true) {}

TNtuple *TreeCombiner::getDataTree() const {
  TString thisFileName=theFileName;
  thisFileName.ReplaceAll("whatsample", theLumiNorm->conversionMap.find("data")->second);
  TString fileName=theInputDir+"/"+thisFileName;
  TFile *file=new TFile(fileName.Data());
  if(file==0) {
    cout << "[TreeCombiner] *** Warning: File: " << fileName << " not found" << endl;
    return 0;
  }

  TNtuple *origNtuple=(TNtuple *)file->Get("ntuple");

  return origNtuple;
}



TNtuple *TreeCombiner::getCombinedTree(const vector<TString>& sampleNames) {
  // Create the new tree
  /*
  TString ntupleVars;
  ntupleVars="run"; 
  ntupleVars+=":lumi"; 
  ntupleVars+=":event"; 
  ntupleVars+=":weight"; 
  ntupleVars+=":leadCharge";
  ntupleVars+=":subleadCharge";
  ntupleVars+=":leadPt";
  ntupleVars+=":subleadPt";
  ntupleVars+=":diLeptInvMass";
  ntupleVars+=":dileptLeadDeltaPhi";
  ntupleVars+=":leptMinusCmCosTheta";
  ntupleVars+=":Analyzed";
  ntupleVars+=":Dilepton";
  ntupleVars+=":MassWindow"; 
  ntupleVars+=":RedMETcut";
  ntupleVars+=":JetVeto";
  ntupleVars+=":LeptonVeto";
  ntupleVars+=":AllCutsPassed";
  */

  // Create the TNtuple to store the variables to be used by the NN
  TString name("NewDiLeptonTree");
  TNtuple *ret(0);
  // TNtuple *ret=new TNtuple(name.Data(), name.Data(), ntupleVars.Data());
  // ret->SetDirectory(0);
  bool initializeDestTree=true;
  int numberOfBranches=0;

  // Loop over all the samples and add them to the new tree with the correct normalization
  for(vector<TString>::const_iterator sampleName=sampleNames.begin();
      sampleName!=sampleNames.end(); ++sampleName) {
    // Retrieve the original tree from the file
    TString thisFileName=theFileName;
    thisFileName.ReplaceAll("whatsample", theLumiNorm->conversionMap.find(*sampleName)->second);
    TString fileName=theInputDir+"/"+thisFileName;
    TFile *file=new TFile(fileName.Data());
    if(file==0) {
      cout << "[TreeCombiner] *** Warning: file " 
	   << fileName.Data() << " not found!" << endl;
      continue;
    }
    cout << "File: " << file->GetName() << endl;
    TNtuple *origNtuple=(TNtuple*) file->Get("ntuple");
    if(origNtuple==0) {
      cout << "[TreeCombiner] *** Error: invalid ntuple in file " 
	   << fileName.Data() << "!" << endl;
      throw std::exception();
    }

    // If first time, initialize destination tree
    if(initializeDestTree) {
      TString ntupleVars;
      numberOfBranches=origNtuple->GetNbranches();
      for(int j=0; j<numberOfBranches; ++j) {
	if(j!=0) ntupleVars+=":";
	ntupleVars+=origNtuple->GetListOfBranches()->At(j)->GetName();
	treeVars[origNtuple->GetListOfBranches()->At(j)->GetName()]=-99999.;
	//treeVars[origNtuple->GetListOfBranches()->At(j)->GetName()];
	varNames.push_back(origNtuple->GetListOfBranches()->At(j)->GetName());
      }
      ret=new TNtuple(name.Data(), name.Data(), ntupleVars.Data());
      if(ret==0) {
	cout << "[TreeCombiner] *** Error: initialization of combined " 
	     << "tree failed!" << endl;
	throw std::exception();	
      }
      if(ret->GetNbranches()!=numberOfBranches) {
	cout << "[TreeCombiner] *** Error: wrong initialization of " 
	     << "combined tree: original tree contains " << numberOfBranches 
	     << " branches, combined tree " << ret->GetNbranches() 
	     << "!" << endl;
	throw std::exception();	
      }
      ret->SetDirectory(0);
      initializeDestTree=false;
    }

    // Set the branches
    for(map<TString, float>::iterator it=treeVars.begin();
	it!=treeVars.end(); ++it) {
      origNtuple->SetBranchAddress(it->first.Data(), &(it->second));
    }
    /*
    origNtuple->SetBranchAddress("run", &run);
    origNtuple->SetBranchAddress("lumi", &lumi);
    origNtuple->SetBranchAddress("event", &event);
    origNtuple->SetBranchAddress("weight", &weight);
    origNtuple->SetBranchAddress("leadCharge", &leadCharge);
    origNtuple->SetBranchAddress("subleadCharge", &subleadCharge);
    origNtuple->SetBranchAddress("leadPt", &leadPt);
    origNtuple->SetBranchAddress("subleadPt", &subleadPt);
    origNtuple->SetBranchAddress("diLeptInvMass", &diLeptInvMass);
    origNtuple->SetBranchAddress("dileptLeadDeltaPhi", &dileptLeadDeltaPhi);
    origNtuple->SetBranchAddress("leptMinusCmCosTheta", &leptMinusCmCosTheta);
    */

    // Get the scale factor for normalization
    double scale=theLumiNorm->getScaleFactor(*sampleName);
    
    // Loop over the entries and add the event to the new ntuple
    for(int entry=0; entry!=origNtuple->GetEntries(); entry++) {
      origNtuple->GetEntry(entry);

      // Selection 
      map<TString, float>::const_iterator selIt, endIt=treeVars.end();

      // Require full selection
      if(fullSelectionReq) {
	selIt=treeVars.find("AllCutsPassed");
	if(selIt!=endIt) {
	  if(selIt->second<0.5) continue;
	}
	else {
	  cout << "[TreeCombiner] *** Error: AllCutsPassed not found in the tree " 
	       << "for sample " << (*sampleName).Data() << endl;
	  throw std::exception();
	}
      }
      else {
	// Cut on the dilepton invariant mass
	selIt=treeVars.find("diLeptInvMass");
	if(selIt!=treeVars.end()) {
	  if(selIt->second<massMinCut || selIt->second>massMaxCut) continue;
	}
	else {
	  cout << "[TreeCombiner] *** Error: diLeptInvMass not found in "
	       << "the tree for sample " << (*sampleName).Data() << endl;
	  throw std::exception();
	}
      }

      const int dimBr(numberOfBranches);
      float values[dimBr];
      for(int h=0; h<numberOfBranches; ++h) {
	values[h]=treeVars[varNames[h]];  // to preserve the order
	if(varNames[h].CompareTo("weight")==0) 
	  values[h]*=scale;
      }
      /*
      float values[11] = {0.};
      int i=0;
      values[i++]=run; 
      values[i++]=lumi; 
      values[i++]=event; 
      values[i++]=weight*scale; 
      values[i++]=leadCharge;
      values[i++]=subleadCharge;
      values[i++]=leadPt;
      values[i++]=subleadPt;
      values[i++]=diLeptInvMass;
      values[i++]=dileptLeadDeltaPhi;
      values[i++]=leptMinusCmCosTheta;
      */
      ret->Fill(values);
    }
    origNtuple->ResetBranchAddresses();
  }
  return ret;
}

TreeCombiner::~TreeCombiner() {}

void TreeCombiner::setMassCut(double min, double max) {
  massMinCut = min;
  massMaxCut = max;
}

void TreeCombiner::requireFullSection(bool doit) {
  fullSelectionReq=doit;
}

