// -*- C++ -*-
//
// Package:    HLTPathFiring
// Class:      HLTPathFiring
// 
/**\class HLTPathFiring HLTPathFiring.cc DQM/HLTPathFiring/src/HLTPathFiring.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/

// system include files
#include <memory>
#include <iomanip>

// User include files
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

// HLT collection objects
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"

// Root files for inclusion
#include "TFile.h"
#include "TMath.h"

//
// class declaration
//

using namespace edm;
using namespace trigger;
using std::vector;
using std::string;


class HLTPathFiring : public edm::EDAnalyzer {
public:
  explicit HLTPathFiring(const edm::ParameterSet&);
  ~HLTPathFiring();

private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  bool checkTriggerName(std::string const);

  // ----------member data ---------------------------

  bool debugPrint;
  bool outputPrint;

  std::string hltProcessName;
  std::string lumiProcessName;
  std::string lumiLabel;

  HLTConfigProvider hltConfig_;

  // Do we need to print every run? Or only at the end of the job?
  bool printAtEndRun;

  vector<std::string> listOfPaths;
  
  bool isFirstRun;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
HLTPathFiring::HLTPathFiring(const edm::ParameterSet& iConfig) {

  if(debugPrint) std::cout << "Inside Constructor" << std::endl;

  hltProcessName = iConfig.getUntrackedParameter<std::string>("hltProcessName", "TEST");
  lumiProcessName = iConfig.getUntrackedParameter<std::string>("lumiProcessName", hltProcessName);
  lumiLabel = iConfig.getUntrackedParameter<std::string>("lumiLabel", "hltLumiScalers");

  debugPrint = iConfig.getUntrackedParameter<bool>("verbose", false);
  outputPrint = iConfig.getUntrackedParameter<bool>("printOutput", false);

  printAtEndRun = iConfig.getUntrackedParameter<bool>("printAtEndRun", false);

  listOfPaths = iConfig.getUntrackedParameter<std::vector<std::string> >("ListOfPaths", std::vector<std::string>());

  if(debugPrint) std::cout << "Now Printing verbose messages"  << std::endl;
  if(outputPrint) std::cout << "Now Printing printOutput messages" << std::endl;

}


HLTPathFiring::~HLTPathFiring() { 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called for each event  ------------
void HLTPathFiring::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using std::string;

  if(debugPrint || outputPrint) std::cout << "============================  Inside analyze "
					  << iEvent.id().run() << "  "
					  << iEvent.id().luminosityBlock() << "  "
					  << iEvent.id().event() << "  "
					  << " ============================" << std::endl;

  edm::Handle<edm::TriggerResults> triggerResults;
  iEvent.getByLabel(InputTag("TriggerResults","", hltProcessName), triggerResults);
   
  if(!triggerResults.isValid()) {
    if(debugPrint) std::cout << "Trigger results not valid" << std::endl;
    return;
  } 
  else {
    if(debugPrint) std::cout << "Found trigger results" << std::endl;
  }
   
  // Sanity check
  assert(triggerResults->size()==hltConfig_.size());

  for(unsigned int iHltPath=0; iHltPath<hltConfig_.size(); ++iHltPath) {

    if(listOfPaths.size()==0) {
      std::cout << " ####### Trigger path: " << hltConfig_.triggerName(iHltPath).c_str() << ": " << triggerResults->accept(iHltPath) << std::endl; 
    }
    else {
      if( !checkTriggerName( hltConfig_.triggerName(iHltPath) ) ) continue; 
      std::cout << " ####### Trigger path: " << hltConfig_.triggerName(iHltPath).c_str() << ": " << triggerResults->accept(iHltPath) << std::endl; 
    }

  }
}


// ------------ method called once each job just before starting event loop  ------------
void HLTPathFiring::beginJob() {

  if(debugPrint) std::cout << "Inside begin job" << std::endl; 

  isFirstRun = true;
}

// ------------ method called once each job just after ending the event loop  ------------
void HLTPathFiring::endJob() {

  if(debugPrint) std::cout << "Inside end job" << std::endl; 
  
}

// ------------ method called when starting to processes a run  ------------
void HLTPathFiring::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {

  if (debugPrint) std::cout << "Inside beginRun" << std::endl;

  bool changed = true;
  if(hltConfig_.init(iRun, iSetup, hltProcessName, changed)) {
    if(debugPrint)
      if(debugPrint) std::cout << "HLT config with process name " 
			       << hltProcessName << " successfully extracted" << std::endl;
  } 
  else {
    if(debugPrint)
      if(debugPrint) std::cout << "Warning, didn't find process " << hltProcessName << std::endl;

    // if you can't find the HLT process and init fails
    // then crash
    assert(false);
  }
} // end of beginRun


// ------------ method called when ending the processing of a run  ------------
void HLTPathFiring::endRun(edm::Run const&, edm::EventSetup const&) {
  using namespace std;

  if(debugPrint) std::cout << "Inside end run, print = " << printAtEndRun << std::endl;
  if(!printAtEndRun) return; 
}


// ------------ method called when starting to processes a luminosity block  ------------
void HLTPathFiring::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

// ------------ method called when ending the processing of a luminosity block  ------------
void HLTPathFiring::endLuminosityBlock(edm::LuminosityBlock const& iLumiBlock, edm::EventSetup const&) {
}

bool HLTPathFiring::checkTriggerName(std::string const thisTrigName) {
  if(listOfPaths.size()==0) return true; 
  for(unsigned int i=0; i<listOfPaths.size(); ++i) {
    if( thisTrigName.find( listOfPaths[i] ) != std::string::npos ) return true;
  }
  return false;
}


// define this as a plug-in
DEFINE_FWK_MODULE(HLTPathFiring);
