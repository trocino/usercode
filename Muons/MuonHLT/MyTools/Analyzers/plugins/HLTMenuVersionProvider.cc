/** \class HLTMenuVersionProvider
 *  Class to extract the HLT menu version from provenance
 *
 *  $Date: 2012/11/05 15:56:24 $
 *  $Revision: 1.0 $
 *  \authors D. Trocino - INFN Torino <daniele.trocino@cern.ch>
 */

// System include files
#include <memory>

// Framework include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// HLT objects
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

// ROOT include files

//
// Class declaration
//

using namespace edm;
using namespace trigger;
using std::vector;
using std::string;


class HLTMenuVersionProvider : public edm::EDAnalyzer {
public:
  explicit HLTMenuVersionProvider(const edm::ParameterSet&);
  ~HLTMenuVersionProvider();

private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  // ----------member data ---------------------------

  std::string hltProcessName;
  bool debugPrint;

  HLTConfigProvider hltConfig_;
};

//
// Constants, enums and typedefs
//

//
// Static data member definitions
//

//
// Constructors and destructor
//
HLTMenuVersionProvider::HLTMenuVersionProvider(const edm::ParameterSet& iConfig) {
  debugPrint = iConfig.getUntrackedParameter<bool>("verbose", false);
  hltProcessName = iConfig.getUntrackedParameter<std::string>("hltProcessName", "TEST");
  if(debugPrint) std::cout << "HLT process name: " << hltProcessName.c_str() << std::endl;
}


HLTMenuVersionProvider::~HLTMenuVersionProvider() { 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


//
// Member functions
//

// ------------ Method called for each event  ------------
void HLTMenuVersionProvider::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {}

// ------------ Method called once each job just before starting event loop  ------------
void HLTMenuVersionProvider::beginJob() {}

// ------------ Method called once each job just after ending the event loop  ------------
void HLTMenuVersionProvider::endJob() {}

// ------------ Method called when starting to processes a run  ------------
void HLTMenuVersionProvider::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {

  bool changed = true;
  if(hltConfig_.init(iRun, iSetup, hltProcessName, changed)) {
    if(debugPrint) 
      std::cout << "Run " << iRun.run() << ": HLT config with process name " 
		<< hltProcessName << " successfully extracted" << std::endl;
  } 
  else {  // If you can't find the HLT process, then crash
    std::cout << "Warning, didn't find process " << hltProcessName << std::endl;
    assert(false);
  }

  if(changed) {
    std::cout << "Run " << iRun.run() << ": HLT menu changed: " 
	      << hltConfig_.tableName() << std::endl;
  }
} // end of beginRun


// ------------ Method called when ending the processing of a run  ------------
void HLTMenuVersionProvider::endRun(edm::Run const&, edm::EventSetup const&) {}


// ------------ Method called when starting to processes a luminosity block  ------------
void HLTMenuVersionProvider::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

// ------------ Method called when ending the processing of a luminosity block  ------------
void HLTMenuVersionProvider::endLuminosityBlock(edm::LuminosityBlock const& iLumiBlock, edm::EventSetup const&) {}

// Define this as a plug-in
DEFINE_FWK_MODULE(HLTMenuVersionProvider);
