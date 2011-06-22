
/*
 *  See header file for a description of this class.
 *
 *  $Date: 2007/08/10 19:15:15 $
 *  $Revision: 1.1 $
 *  \author G. Cerminara - NEU Boston & INFN Torino
 */

#include "TreeSum.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TString.h"

#include <iostream>

using namespace std;

TreeSum::TreeSum(TTree *sigTree, TTree *bckTree) : theSigTree(sigTree),
						   theBckTree(bckTree),
						   sigScale(1),
						   bckScale(1) {

  theSigTree->SetBranchAddress("run", &run);
  theSigTree->SetBranchAddress("event", &event);
  theSigTree->SetBranchAddress("weight", &weight);
  theSigTree->SetBranchAddress("type", &type);
  theSigTree->SetBranchAddress("DiLeptMinv", &DiLeptMinv);
  theSigTree->SetBranchAddress("LeadLeptPt", &LeadLeptPt);
  theSigTree->SetBranchAddress("NoLeadLeptPt", &NoLeadLeptPt);
  theSigTree->SetBranchAddress("DeltaPhiLeadLeptZ", &DeltaPhiLeadLeptZ); // delta phi leading mu - Z boson
  theSigTree->SetBranchAddress("CosThetaStarLeptMinus", &CosThetaStarLeptMinus);  // cos(theta*)  mu -
  theSigTree->SetBranchAddress("DiLeptDeltaEta",&DiLeptDeltaEta);// delta eta di-leptons
  theSigTree->SetBranchAddress("TransMass4Body", &TransMass4Body); // Transverse mass of the 4-body system neglecting di-lepton mass
  theSigTree->SetBranchAddress("MinMass4Body", &MinMass4Body); // Min invariant mass of the 4-body system



  theBckTree->SetBranchAddress("run", &run);
  theBckTree->SetBranchAddress("event", &event);
  theBckTree->SetBranchAddress("weight", &weight);
  theBckTree->SetBranchAddress("type", &type);
  theBckTree->SetBranchAddress("DiLeptMinv", &DiLeptMinv);
  theBckTree->SetBranchAddress("LeadLeptPt", &LeadLeptPt);
  theBckTree->SetBranchAddress("NoLeadLeptPt", &NoLeadLeptPt);
  theBckTree->SetBranchAddress("DeltaPhiLeadLeptZ", &DeltaPhiLeadLeptZ); // delta phi leading mu - Z boson
  theBckTree->SetBranchAddress("CosThetaStarLeptMinus", &CosThetaStarLeptMinus);  // cos(theta*)  mu -
  theBckTree->SetBranchAddress("DiLeptDeltaEta",&DiLeptDeltaEta);// delta eta di-leptons
  theBckTree->SetBranchAddress("TransMass4Body", &TransMass4Body); // Transverse mass of the 4-body system neglecting di-lepton mass
  theBckTree->SetBranchAddress("MinMass4Body", &MinMass4Body); // Min invariant mass of the 4-body system


}


TreeSum::~TreeSum(){}




TNtuple *TreeSum::getSumTree() {
  TString ntupleVars;
  ntupleVars = "run"; // Run #
  ntupleVars += ":event"; // event #
  ntupleVars += ":weight"; // event weight
  ntupleVars += ":type"; // event type: 1 for signal 0 for background
  ntupleVars += ":DiLeptMinv"; // di-muon invariant mass
  ntupleVars += ":LeadLeptPt"; // pT leading mu
  ntupleVars += ":NoLeadLeptPt"; // pT non-leading mu
  ntupleVars += ":DeltaPhiLeadLeptZ"; // delta phi leading mu - Z boson
  ntupleVars += ":CosThetaStarLeptMinus";  // cos(theta*)  mu -
  ntupleVars += ":TransMass4Body"; // Transverse mass of the 4-body system neglecting di-lepton mass
  ntupleVars += ":MinMass4Body";// Min invariant mass of the 4-body system
  
  // Create the TNtuple to store the variables to be used by the NN
  TString name(theSigTree->GetName());
  name+=theBckTree->GetName();
  TNtuple *ret =  new TNtuple(name.Data(), name.Data(),ntupleVars.Data());

  // Signal
  for(int entry = 0; entry != theSigTree->GetEntries(); entry++) {
    theSigTree->GetEntry(entry);
//     // Cut on the di-lepton invariant mass....
//     if(DiLeptMinv < massMinCut || DiLeptMinv > massMaxCut) continue;

    float values[11] = {0.};    
    int i = 0;
    values[i++] = run; // Run #
    values[i++] = event; // event #
    values[i++] = weight*sigScale; // event weight
    values[i++] = type; // event type: 1 for signal 0 for background
    values[i++] = DiLeptMinv; // di-muon invariant mass
    values[i++] = LeadLeptPt; // pT leading mu
    values[i++] = NoLeadLeptPt; // pT non-leading mu
    values[i++] = DeltaPhiLeadLeptZ; // delta phi leading mu - Z boson
    values[i++] = CosThetaStarLeptMinus;  // cos(theta*)  mu -
    values[i++] = TransMass4Body;
    values[i++] = MinMass4Body;

//     cout << "[2a] run: " << run << " event: " << event << " weight: " << weight
// 	 << " minv: " << DiLeptMinv << endl;
    ret->Fill(values);
  }

  // Background
  for(int entry = 0; entry != theBckTree->GetEntries(); entry++) {
    theBckTree->GetEntry(entry);
//     // Cut on the di-lepton invariant mass....
//     if(DiLeptMinv < massMinCut || DiLeptMinv > massMaxCut) continue;

    float values[11] = {0.};    
    int i = 0;
    values[i++] = run; // Run #
    values[i++] = event; // event #
    values[i++] = weight*bckScale; // event weight
    values[i++] = type; // event type: 1 for signal 0 for background
    values[i++] = DiLeptMinv; // di-muon invariant mass
    values[i++] = LeadLeptPt; // pT leading mu
    values[i++] = NoLeadLeptPt; // pT non-leading mu
    values[i++] = DeltaPhiLeadLeptZ; // delta phi leading mu - Z boson
    values[i++] = CosThetaStarLeptMinus;  // cos(theta*)  mu -
    values[i++] = TransMass4Body;
    values[i++] = MinMass4Body;


//     cout << "[2b] run: " << run << " event: " << event << " weight: " << weight
// 	 << " minv: " << DiLeptMinv << endl;
    ret->Fill(values);
    
  }

//   file->cd();
//   allMcTree->Write();

  return ret;
  
}


void TreeSum::setScaleFactors(double scale1, double scale2) {
  sigScale = scale1;
  bckScale = scale2;
}
