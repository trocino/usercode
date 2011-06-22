
/*
 *  See header file for a description of this class.
 *
 *  $Date: 2007/12/04 00:11:38 $
 *  $Revision: 1.2 $
 *  \author G. Cerminara - NEU Boston & INFN Torino
 */

#include "EnsembleProducer.h"

#include "RndmGenerator.h"


#include "TTree.h"

#include "TFile.h"
#include "TNtuple.h"

#include <iostream>
#include <math.h>

using namespace std;


EnsembleProducer::EnsembleProducer(const TString& name, TTree *tree, const double expectedEvents,
				   double massMin, double massMax) : theTree(tree),
								     massMinCut(massMin),
								     massMaxCut(massMax) {
  int entries = theTree->GetEntries();
  cout << "The Tree has " << entries << " entries" << endl;
  theTree->SetBranchAddress("run", &run);
  theTree->SetBranchAddress("event", &event);
  theTree->SetBranchAddress("weight", &weight);
  theTree->SetBranchAddress("type", &type);
  theTree->SetBranchAddress("DiLeptMinv", &DiLeptMinv);
  theTree->SetBranchAddress("LeadLeptPt", &LeadLeptPt);
  theTree->SetBranchAddress("NoLeadLeptPt", &NoLeadLeptPt);
  theTree->SetBranchAddress("DeltaPhiLeadLeptZ", &DeltaPhiLeadLeptZ); // delta phi leading mu - Z boson
  theTree->SetBranchAddress("CosThetaStarLeptMinus", &CosThetaStarLeptMinus);  // cos(theta*)  mu -
  theTree->SetBranchAddress("DiLeptDeltaEta",&DiLeptDeltaEta);// delta eta di-leptons
  theTree->SetBranchAddress("TransMass4Body", &TransMass4Body); // Transverse mass of the 4-body system neglecting di-lepton mass
  theTree->SetBranchAddress("MinMass4Body", &MinMass4Body); // Min invariant mass of the 4-body system

  double yeld = getYeld();
  cout << "The Tree has " << entriesInMassWind << " entries in the selected mass range" << endl;
  cout << "yeld: " << yeld << endl;
  
//   nExpectedAverage = rounding(expectedEvents);
  nExpectedAverage = expectedEvents;

  nEns = entriesInMassWind*yeld/nExpectedAverage;
  cout << "# of ensembles: " << nEns << endl;
  
  // build the random generator
  rndm = new RndmGenerator(10); // FIXME: handling of the seed
  TString fileName = name+"EnsembleProducer.root";
  // The file which will store the trees
  theFile = new TFile(fileName.Data(),"RECREATE");
  
  initialize();
  produce();
  write();
  theTree->ResetBranchAddresses();
}

EnsembleProducer::~EnsembleProducer(){}


double EnsembleProducer::getYeld() {
  // loop over all entries and get maximum and mean event weight
  double sumWeight = 0;
  maxWeight = -99999;
  int counter = 0;
  for(int i = 0; i != theTree->GetEntries(); i++) {
    theTree->GetEntry(i);
    if(massOutWindow()) continue;
    counter++;
    sumWeight += weight;
    if(weight > maxWeight) {
      maxWeight = weight;
    }
  }
  entriesInMassWind = counter;
  double meanWeight = sumWeight/counter;
  cout << "[EnsembleProducer::getYeld] Mean weight: " << meanWeight
       << " max weight: " <<  maxWeight << endl;
  double ret = meanWeight/maxWeight;
  
  return ret;
  
}




int EnsembleProducer::rounding(double expectedEvents) const {
  return rint(expectedEvents);
}


TNtuple *EnsembleProducer::getEnsemble(int i) {
  if(theEnsambles.find(i) == theEnsambles.end()) {
    // build the tree
      // List the variables to be stored in the ntuple

    string ntupleVars;
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
    theFile->cd();
    TString name = "nsmb";
    name+= theTree->GetName();
    name+=i;
    theEnsambles[i] =  new TNtuple(name.Data(), name.Data(),ntupleVars.c_str());
  }
  return theEnsambles[i];
}


void EnsembleProducer::produce() {
  // loop over all the events
  for(int i = 0; i != theTree->GetEntries(); i++) {
    theTree->GetEntry(i);
    if(massOutWindow()) continue;
    float values[11] = {0.};    
    int i = 0;
    values[i++] = run; // Run #
    values[i++] = event; // event #
    values[i++] = weight; // event weight
    values[i++] = type; // event type: 1 for signal 0 for background
    values[i++] = DiLeptMinv; // di-muon invariant mass
    values[i++] = LeadLeptPt; // pT leading mu
    values[i++] = NoLeadLeptPt; // pT non-leading mu
    values[i++] = DeltaPhiLeadLeptZ; // delta phi leading mu - Z boson
    values[i++] = CosThetaStarLeptMinus;  // cos(theta*)  mu -
    values[i++] = TransMass4Body;
    values[i++] = MinMass4Body;
  
    // Discard or reject the event on the basis of its weight
    double random = rndm->extractFlat(0,maxWeight);
    if(weight < random) continue;
    
    // Extract the Ensemble index number
    int ensIndex = rndm->extractInteg(nEns);
//     cout << "[1] run: " << run << " event: " << event << " weight: " << weight
// 	 << " minv: " << DiLeptMinv << endl;
    //fill the corresponding ensemble
    getEnsemble(ensIndex)->Fill(values);
  }
}


int EnsembleProducer::getNEnsembles() const {
  return nEns;
}

TNtuple* EnsembleProducer::getEnsembleN(int i) {
  return theEnsambles[i];
}

void EnsembleProducer::setSeed(unsigned int seed) {
  rndm->setSeed(seed);
}


void EnsembleProducer::initialize() {
  // List the variables to be stored in the ntuple
  string ntupleVars;
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

    
  for(int i=0; i !=nEns; i++) {
    // Create the TNtuple to store the variables to be used by the NN
    theFile->cd();
    TString name = "nsmb";
    name+= theTree->GetName();
    name+=i;
    theEnsambles[i] =  new TNtuple(name.Data(), name.Data(),ntupleVars.c_str());
  }
}

void EnsembleProducer::write() {
  theFile->cd();
  for(map<int, TNtuple *>::iterator ens = theEnsambles.begin();
      ens != theEnsambles.end(); ens++) {
    (*ens).second->Write();
  }
}


int EnsembleProducer::getExpectedAverage() const {
  return nExpectedAverage;
}


bool EnsembleProducer::massOutWindow() {
  if(DiLeptMinv < massMinCut || DiLeptMinv > massMaxCut) return true;
  else return false;
}
