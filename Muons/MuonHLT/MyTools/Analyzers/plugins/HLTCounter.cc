// -*- C++ -*-
//
// Package:    HLTCounter
// Class:      HLTCounter
// 
/**\class HLTCounter HLTCounter.cc DQM/HLTCounter/src/HLTCounter.cc

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
#include "DataFormats/Scalers/interface/LumiScalers.h"

// Root files for inclusion
#include "TFile.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TMath.h"
#include "TStyle.h"

//
// class declaration
//

using namespace edm;
using namespace trigger;
using std::vector;
using std::string;


class HLTCounter : public edm::EDAnalyzer {
public:
  explicit HLTCounter(const edm::ParameterSet&);
  ~HLTCounter();

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

  std::string plotDirectoryName;

  HLTConfigProvider hltConfig_;

  // count the number of times each path fired
  // store it in this vector
  unsigned         countsAllPathByLS_OR;  // OR (with no prescales)
  vector<unsigned> countsPerPathByLS_EX;  // exclusive (with no prescales)

  vector<unsigned> countsPerPath;

  vector<unsigned> countsPerPathByLS;
  vector<unsigned> countsPerPathByLS_L1;
  vector<unsigned> countsPerPathByLS_L2;
  vector<unsigned> countsPerPathByLS_L3;
  vector<unsigned> countsPerPathByLS_L3iso;

  //vector<unsigned> LSPerPath;

  bool useInstLumi;  
  double averageInstLumi;
  int NumInstLumis;
  int NumInstLumisPrevLS;
  std::vector<double> InstLumiArray;

  int TotNumLumiSections;
  double NormIntLumi;
  int nLumisInAverage;

  // Do we need to print every run? Or only at the end of the job?
  bool printAtEndRun;

  // Do we want to use original HLT results?
  // If we do, what is the input tag for them?

  bool filterByOriginalHLTResults;
  std::string originalHLTName;
  std::string requiredOriginalPath;
  int originalPathIndex;
  HLTConfigProvider originalConfig_;

  vector<std::string> listOfPaths;

  bool doRatesByPtThreshold;

  void addLumiToAverage(double inLumi, int nLumi);
  void clearLumiAverage();

  vector<double> averagePrescale;  
  vector<double> averagePrescale_L1;  
  vector<double> averagePrescale_L2;  
  vector<double> averagePrescale_L3;  
  vector<double> averagePrescale_L3iso;  
  void addPrescaleToAverage(double evtPs, unsigned int pos, vector<double>& avgPs, vector<unsigned>& cntPerLS);
  
  bool isFirstRun;

  // Output file
  std::string fileName;
  TFile *fout;

  // Histograms
  TH1D *hRateVsPt;

  TH1D *hRatePerPath;
  TH1D *hRatePerPath_L1;
  TH1D *hRatePerPath_L2;
  TH1D *hRatePerPath_L3;
  TH1D *hRatePerPath_L3iso;

  TProfile *             hRate_OR;  // OR (with no prescales)
  std::vector<TProfile*> hRates_EX;  // exclusive (with no prescales)

  std::vector<TProfile*> hRates;
  std::vector<TProfile*> hRates_L1;
  std::vector<TProfile*> hRates_L2;
  std::vector<TProfile*> hRates_L3;
  std::vector<TProfile*> hRates_L3iso;

  unsigned int nRateBins;
  double firstRateBin;
  double lastRateBin;
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
HLTCounter::HLTCounter(const edm::ParameterSet& iConfig) {

  if(debugPrint) std::cout << "Inside Constructor" << std::endl;

  hltProcessName = iConfig.getUntrackedParameter<std::string>("hltProcessName", "TEST");
  lumiProcessName = iConfig.getUntrackedParameter<std::string>("lumiProcessName", hltProcessName);
  lumiLabel = iConfig.getUntrackedParameter<std::string>("lumiLabel", "hltLumiScalers");

  nLumisInAverage = iConfig.getUntrackedParameter<int>("nLumisInAverage", 1);
  NormIntLumi = iConfig.getUntrackedParameter<double>("NormIntLumi", 10000.);

  debugPrint = iConfig.getUntrackedParameter<bool>("verbose", false);
  outputPrint = iConfig.getUntrackedParameter<bool>("printOutput", false);
  useInstLumi = iConfig.getUntrackedParameter<bool>("useInstLumi", false);

  printAtEndRun = iConfig.getUntrackedParameter<bool>("printAtEndRun", false);

  filterByOriginalHLTResults = iConfig.getUntrackedParameter<bool>("filterByOriginalHLTResults", false);
  originalHLTName = iConfig.getUntrackedParameter<std::string>("originalHLTName", "HLT");
  requiredOriginalPath = iConfig.getUntrackedParameter<std::string>("requiredOriginalPath", "HLT_Mu20_v");

  listOfPaths = iConfig.getUntrackedParameter<std::vector<std::string> >("ListOfPaths", std::vector<std::string>());

  doRatesByPtThreshold = iConfig.getUntrackedParameter<bool>("RatesByPtThreshold", false);
  
  if( doRatesByPtThreshold==true && filterByOriginalHLTResults==false ) {
    std::cout << "ERROR: With doRatesByPtThreshold==true, filterByOriginalHLTResults must be true!" << std::endl;
    assert(false);
  }

  fileName = iConfig.getUntrackedParameter<std::string>("FileName", "output.root");

  nRateBins    = iConfig.getUntrackedParameter<unsigned int>("NRateBins", 100);
  firstRateBin = iConfig.getUntrackedParameter<double>("FirstRateBin", 2000.);
  lastRateBin  = iConfig.getUntrackedParameter<double>("LastRateBin", 7000.);

  if(debugPrint) std::cout << "Now Printing verbose messages"  << std::endl;
  if(outputPrint) std::cout << "Now Printing printOutput messages" << std::endl;

  if(debugPrint) std::cout << "Got plot dirname = " << plotDirectoryName << std::endl;
  if(debugPrint) std::cout << "Got plot hltProcessName = " << hltProcessName << std::endl;

  averageInstLumi = 0.;
  NumInstLumis = 0;
  NumInstLumisPrevLS = 0;
  TotNumLumiSections = 0;

  fout = new TFile(fileName.c_str(), "RECREATE");

  if(doRatesByPtThreshold) {
    //hRateVsPt = new TH1D("RateVsPt", "Rate vs p_{T} threshold", 501, -0.5, 500.5); 
    hRateVsPt = new TH1D("RateVsPt", "Rate vs p_{T} threshold", 100, 0.5, 100.5); 
    hRateVsPt->Sumw2();
  }
}


HLTCounter::~HLTCounter() { 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called for each event  ------------
void HLTCounter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using std::string;

  if(debugPrint || outputPrint) std::cout << "============================  Inside analyze "
					  << iEvent.id().run() << "  "
					  << iEvent.id().luminosityBlock() << "  "
					  << iEvent.id().event() << "  "
					  << " ============================" << std::endl;
  unsigned currentLS = iEvent.id().luminosityBlock();

  /////////////////////////////////////////////
  //
  //  Access All collections needed 
  //  for doing some checks
  //
  //////////////////////////////////////////////


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

  edm::Handle<edm::TriggerResults> originalTrigResults;
  if(filterByOriginalHLTResults) {
    iEvent.getByLabel(InputTag("TriggerResults","", originalHLTName), originalTrigResults);

    if(!originalTrigResults.isValid()) {
      std::cout << "Trigger results not valid for process " << originalHLTName << std::endl;
      return;
    } 
    else {   
      if(debugPrint) std::cout << "Found trigger results" << std::endl;
    }

    // Sanity check
    assert(originalTrigResults->size()==originalConfig_.size());

    if(debugPrint) std::cout << "Checking to see if you fired the original path" << std::endl;
     
    // If didn't fire the path, don't bother processing the event    
    if( !(originalTrigResults->accept(originalPathIndex)) ) {
      if(debugPrint) std::cout << "Didn't fire the path, skipping the event" << std::endl;
      return;
    }

    if(debugPrint) std::cout << "Fired original path, continue processing" << std::endl;     
  } // end if filterByOriginalHLTResults
   

  /////////////////////////////////////////////////
  //
  //   Get lumi and use it in an average
  //
  ////////////////////////////////////////////////
     
  edm::Handle<LumiScalersCollection> lumiScalers;
  double currentLumi = -1;

  if(useInstLumi) {
    //bool lumiHandleOK = iEvent.getByLabel( InputTag("hltScalersRawToDigi", "", lumiProcessName), lumiScalers );
    bool lumiHandleOK = iEvent.getByLabel( InputTag(lumiLabel, "", lumiProcessName), lumiScalers );

    if( !lumiHandleOK ) {
      if(debugPrint) std::cout << "Problem unpacking lumi, skipping event" << std::endl;
      return;
    } 
    else {
      if(debugPrint) std::cout << "Got the lumi" << std::endl;
    } // end if lumi handle OK
     
    if( lumiScalers->size() ) {
      LumiScalersCollection::const_iterator it3 = lumiScalers->begin();
      unsigned int lumisection = it3->sectionNumber();
      currentLumi = it3->instantLumi();
      if(currentLumi<=0) {
	if(debugPrint) std::cout << "Problem with lumi (l = " << currentLumi << "), skipping the event" << std::endl;
	return;
      }
      if(debugPrint) std::cout << ".:.:.: Current Lumi :.:.:. " << currentLumi << std::endl;
      if(lumisection) {
	addLumiToAverage(it3->instantLumi(), NumInstLumis);
      } // end if lumi section
    } // end if there are lumiScalers
  } // end if use lumi

  //
  // Prescales
  // 
  // This will store the prescale values of original triggers
  std::pair<int,int>  psValueOrigCombo;
  if(filterByOriginalHLTResults) {
    if( originalTrigResults->accept(originalPathIndex) ) {
      if(outputPrint) std::cout << "Reference path " << originalPathIndex << "  " << originalConfig_.triggerName(originalPathIndex)
				<< " fired." << std::endl;

      // Get the ps for this path
      psValueOrigCombo = originalConfig_.prescaleValues(iEvent, iSetup, originalConfig_.triggerName(originalPathIndex) );
      double psWeight = psValueOrigCombo.first * psValueOrigCombo.second;

      if(outputPrint) std::cout << "Reference path prescale L1 = " << psValueOrigCombo.first << "; HLT = " << psValueOrigCombo.second
				<< "; total = " << psWeight << std::endl;
    } 
  }
  else { // filterByOriginalHLTResults == false
    psValueOrigCombo.first = 1;
    psValueOrigCombo.second = 1;
  }

  // This will store the prescale values
  std::pair<int,int>  psValueCombo;

  unsigned int hltPathCount(0);
  unsigned int hltPathsMask(0);
  for(unsigned int iHltPath=0; iHltPath<hltConfig_.size(); ++iHltPath) {

    if( !checkTriggerName( hltConfig_.triggerName(iHltPath) ) ) continue; 

    // Get the ps for this path
    psValueCombo = hltConfig_.prescaleValues(iEvent, iSetup, hltConfig_.triggerName(iHltPath) );
    double psWeight = psValueCombo.first * psValueCombo.second;

    if(outputPrint) std::cout << "Prescale L1 = " << psValueCombo.first << "; HLT = " << psValueCombo.second
			      << "; total = " << psWeight << std::endl;

    // Total event weight:  (L1 ps * HLT ps)_ref * (HLT ps)_path / inst.lumi
    // N.B. path is seeded by same seed as reference path, so no need to multiply by path L1 ps
    double l1Weight    = psValueOrigCombo.first * psValueOrigCombo.second;
    double eventWeight = psValueOrigCombo.first * psValueOrigCombo.second * psValueCombo.second;

    // 
    // Full path
    // 

    if( triggerResults->accept(iHltPath) ) {
      if(outputPrint) std::cout << "Fired Path " << iHltPath << "  " << hltConfig_.triggerName(iHltPath)
				<< std::endl;

      hltPathsMask += std::pow(2, hltPathCount); 

      countsPerPath[hltPathCount]++;
      countsPerPathByLS[hltPathCount]++;

      addPrescaleToAverage(eventWeight, hltPathCount, averagePrescale, countsPerPathByLS);

      hRatePerPath->Fill(hltPathCount, eventWeight/currentLumi); 

      std::string trigName = hltConfig_.triggerName(iHltPath); 

      if( doRatesByPtThreshold ) {
	double trigThr = 0.;

	if( trigName.find("HLTriggerFirstPath") != std::string::npos ||
	    trigName.find("HLTriggerFinalPath") != std::string::npos   ) continue;

	if( originalConfig_.triggerName(originalPathIndex).find("SingleMuOpen") != std::string::npos ) {
	  if( trigName.find("Open") == std::string::npos && 
	      trigName.find("_2") == std::string::npos   &&
	      trigName.compare("HLT_L2SingleMu10") != 0  &&
	      trigName.compare("HLT_L2SingleMu15") != 0    ) continue;

	  trigName = trigName.substr( trigName.find("SingleMu")+8 ); 
	  if( trigName.find("_") != std::string::npos )
	    trigName = trigName.substr( 0, trigName.find("_") );
	  if( trigName.compare("Open")==0 )
	    trigThr = 0.;
	  else 
	    trigThr = atof( trigName.c_str() ); 
	}

	else if( originalConfig_.triggerName(originalPathIndex).find("SingleMu12") != std::string::npos ) {
	  if( trigName.find("Open") != std::string::npos ||
	      trigName.find("_2") != std::string::npos   ||
	      trigName.compare("HLT_L2SingleMu10") == 0  ||
	      trigName.compare("HLT_L2SingleMu15") == 0    ) continue;

	  trigName = trigName.substr( trigName.find("SingleMu")+8 ); 
	  if( trigName.find("_") != std::string::npos )
	    trigName = trigName.substr( 0, trigName.find("_") );
	  trigThr = atof( trigName.c_str() ); 
	}

	else if( originalConfig_.triggerName(originalPathIndex).find("HLT_IsoMu") != std::string::npos ) {
	  if( trigName.find("HLT_Mu") != std::string::npos ) continue;

	  trigName = trigName.substr( trigName.find("HLT_IsoMu")+9 ); 
	  if( trigName.find("_") != std::string::npos )
	    trigName = trigName.substr( 0, trigName.find("_") );
	  trigThr = atof( trigName.c_str() ); 
	}

	else if( originalConfig_.triggerName(originalPathIndex).find("HLT_Mu") != std::string::npos ) {
	  if( trigName.find("HLT_IsoMu") != std::string::npos ) continue;

	  trigName = trigName.substr( trigName.find("HLT_Mu")+6 ); 
	  if( trigName.find("_") != std::string::npos )
	    trigName = trigName.substr( 0, trigName.find("_") );
	  trigThr = atof( trigName.c_str() ); 
	  //std::cout << trigName.c_str() << ": " << trigThr << std::endl;
	}

	else {}

	//std::cout << "  Wgh: " << eventWeight << ",  trigger name: " << trigName.c_str() << ",  trigger Thr: " << trigThr << std::endl;
	hRateVsPt->Fill(trigThr, eventWeight); 
      }
    } 

    else {
      if(outputPrint) std::cout << "NOT fired Path " << iHltPath << "  " << hltConfig_.triggerName(iHltPath)
				<< std::endl;
    }

    // 
    // Now L1, L2, L3, L3iso separately
    // 

    if( triggerResults->accept(iHltPath) ) { // if HLT path fired, all intermediate steps fired as well... 
      // L1
      countsPerPathByLS_L1[hltPathCount]++;
      addPrescaleToAverage(l1Weight, hltPathCount, averagePrescale_L1, countsPerPathByLS_L1);
      hRatePerPath_L1->Fill(hltPathCount, l1Weight/currentLumi); 
      // L2
      countsPerPathByLS_L2[hltPathCount]++;
      addPrescaleToAverage(eventWeight, hltPathCount, averagePrescale_L2, countsPerPathByLS_L2);
      hRatePerPath_L2->Fill(hltPathCount, eventWeight/currentLumi); 
      // L3
      countsPerPathByLS_L3[hltPathCount]++;
      addPrescaleToAverage(eventWeight, hltPathCount, averagePrescale_L3, countsPerPathByLS_L3);
      hRatePerPath_L3->Fill(hltPathCount, eventWeight/currentLumi); 
      // L3 isolation
      countsPerPathByLS_L3iso[hltPathCount]++;
      addPrescaleToAverage(eventWeight, hltPathCount, averagePrescale_L3iso, countsPerPathByLS_L3iso);
      hRatePerPath_L3iso->Fill(hltPathCount, eventWeight/currentLumi); 
    }

    else {  // if HLT paths did not fire, check filter by filter
      // Get trigger summary 
      edm::Handle<trigger::TriggerEvent> triggerEvent;
      iEvent.getByLabel(edm::InputTag("hltTriggerSummaryAOD", "", hltProcessName), triggerEvent);

      if(!triggerEvent.isValid()) { 
	std::cout << "TriggerEvent not valid" << std::endl;
	return;
      }

      // Modules in trigger path
      const std::vector<std::string>& moduleLabels(hltConfig_.moduleLabels(iHltPath));
      if(debugPrint) std::cout << "moduleLabels.size() = " << moduleLabels.size() << std::endl;

      const unsigned int m(hltConfig_.size(iHltPath));
      if(debugPrint) std::cout << "hltConfig_.size(iHltPath) = " << m << std::endl;
      const unsigned int lastModuleIndex(triggerResults->index(iHltPath));
      if(debugPrint) std::cout << "lastModuleIndex = " << lastModuleIndex << std::endl;
      assert(lastModuleIndex<m);

      // Results from TriggerEvent product - Attention: must look only for
      // Modules actually run in this path for this event!
      for (unsigned int j=0; j<=lastModuleIndex; ++j) {
	const std::string& moduleLabel(moduleLabels[j]);
	const std::string  moduleType(hltConfig_.moduleType(moduleLabel));
	// Check whether the module is packed up in TriggerEvent product
	const unsigned int filterIndex(triggerEvent->filterIndex(edm::InputTag(moduleLabel, "", hltProcessName)));
	if(filterIndex<triggerEvent->sizeFilters()) {
	  const unsigned int nI = triggerEvent->filterIds(filterIndex).size();
	  const unsigned int nK = triggerEvent->filterKeys(filterIndex).size();
	  assert(nI==nK);
	  if( nI<1 ) continue; 

	  if( moduleType.compare("HLTLevel1GTSeed")==0 ) { // L1
	    countsPerPathByLS_L1[hltPathCount]++;
	    addPrescaleToAverage(l1Weight, hltPathCount, averagePrescale_L1, countsPerPathByLS_L1);
	    hRatePerPath_L1->Fill(hltPathCount, l1Weight/currentLumi); 
	  }
	  else if( moduleType.compare("HLTMuonL2PreFilter")==0 ) { // L2
	    countsPerPathByLS_L2[hltPathCount]++;
	    addPrescaleToAverage(eventWeight, hltPathCount, averagePrescale_L2, countsPerPathByLS_L2);
	    hRatePerPath_L2->Fill(hltPathCount, eventWeight/currentLumi); 
	  }
	  else if( moduleType.compare("HLTMuonL3PreFilter")==0 ) { // L3
	    countsPerPathByLS_L3[hltPathCount]++;
	    addPrescaleToAverage(eventWeight, hltPathCount, averagePrescale_L3, countsPerPathByLS_L3);
	    hRatePerPath_L3->Fill(hltPathCount, eventWeight/currentLumi); 
	  }
	  else if( moduleType.compare("HLTMuonIsoFilter")==0 ) { // L3 isolation
	    countsPerPathByLS_L3iso[hltPathCount]++;
	    addPrescaleToAverage(eventWeight, hltPathCount, averagePrescale_L3iso, countsPerPathByLS_L3iso);
	    hRatePerPath_L3iso->Fill(hltPathCount, eventWeight/currentLumi); 
	  }
	  else {}
	}
      }
    } // end else (HLT path was not fired)

    hltPathCount++;
  }

  // Total rate and exclusive rates
  if(hltPathsMask>0) {
    countsAllPathByLS_OR++; 

    hltPathCount = 0; 
    for(unsigned int iHltPath=0; iHltPath<hltConfig_.size(); ++iHltPath) {
      if( !checkTriggerName( hltConfig_.triggerName(iHltPath) ) ) continue; 
      if(hltPathsMask==std::pow(2, hltPathCount)) countsPerPathByLS_EX[hltPathCount]++;
      hltPathCount++;
    } 
  }

}


// ------------ method called once each job just before starting event loop  ------------
void HLTCounter::beginJob() {

  if(debugPrint) std::cout << "Inside begin job" << std::endl; 

  isFirstRun = true;
}

// ------------ method called once each job just after ending the event loop  ------------
void HLTCounter::endJob() {

  if(debugPrint) std::cout << "Inside end job" << std::endl; 
  
  using namespace std;
  
  if(debugPrint) std::cout << "End of JOB printing summary table " << std::endl;

  cout << "################### JOB TRIGGER SUMMARY ####################### "
       << endl << endl;

  //std::cout << "TotNumLumiSections: " << TotNumLumiSections << std::endl;
  fout->cd();

  if(doRatesByPtThreshold) {
    hRateVsPt->Scale( 1.0 * NormIntLumi / (TotNumLumiSections*23.31) );
    hRateVsPt->Write();
  }

  hRatePerPath->Scale( 1.0 * NormIntLumi / (TotNumLumiSections*23.31) );
  hRatePerPath_L1->Scale( 1.0 * NormIntLumi / (TotNumLumiSections*23.31) );
  hRatePerPath_L2->Scale( 1.0 * NormIntLumi / (TotNumLumiSections*23.31) );
  hRatePerPath_L3->Scale( 1.0 * NormIntLumi / (TotNumLumiSections*23.31) );
  hRatePerPath_L3iso->Scale( 1.0 * NormIntLumi / (TotNumLumiSections*23.31) );

  hRatePerPath->Write();
  hRatePerPath_L1->Write();
  hRatePerPath_L2->Write();
  hRatePerPath_L3->Write();
  hRatePerPath_L3iso->Write();

  hRate_OR->Write();

  unsigned int hltCnt(0);
  for(unsigned iHltPath = 0; iHltPath<hltConfig_.size(); ++iHltPath) {
    if( !checkTriggerName( hltConfig_.triggerName(iHltPath) ) ) continue; 

    string tempName = hltConfig_.removeVersion(hltConfig_.triggerName(iHltPath));

    // cut the string down to size
    // then pad it out with . if it was smaller
    if(tempName.size() > 52) {
      tempName.resize(52);
      tempName.resize(55, '.');
    }

    cout << left << setw(58) << tempName
	 << right << countsPerPath[hltCnt] << endl;

    hRates_EX[hltCnt]->Write();
    hRates[hltCnt]->Write();
    hRates_L1[hltCnt]->Write();
    hRates_L2[hltCnt]->Write();
    hRates_L3[hltCnt]->Write();
    hRates_L3iso[hltCnt]->Write();
    hltCnt++;
  }

  fout->Close();
}

// ------------ method called when starting to processes a run  ------------
void HLTCounter::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {

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

  ////////////// Look up the reference path, and save the index of that path
  if(filterByOriginalHLTResults) {
    if(debugPrint) std::cout << "Extracting original config" << std::endl;
    bool changedOrig = true;
    if(originalConfig_.init(iRun, iSetup, originalHLTName, changedOrig)) {
      if(debugPrint)
        std::cout << "Original HLT config with process name " 
                  << originalHLTName << " successfully extracted" << std::endl;
    } 
    else { 
      if(debugPrint) std::cout << "Warning, didn't find process " << originalHLTName << std::endl;

      // if you can't find the HLT process and init fails
      // then crash
      assert(false);
    }

    originalPathIndex = -1;
    for(unsigned iHltPath=0; iHltPath<originalConfig_.size(); ++iHltPath) {
      string tempName = originalConfig_.triggerName(iHltPath);
      if(tempName.find(requiredOriginalPath) != std::string::npos) {
        originalPathIndex = int(iHltPath);
        std::cout << "Found original path index. Path = "
                  << requiredOriginalPath
                  << " index = " << originalPathIndex
                  << std::endl;
      }
    } // end for each path
  } // end if filter by original 

  //
  // Initialize all counts to be zero
  //
  if(isFirstRun) {
    if(listOfPaths.size()==0) {
      hRatePerPath = new TH1D("RatePerPath", "Rate per path", hltConfig_.size(), -0.5, hltConfig_.size()-0.5); 
      hRatePerPath_L1 = new TH1D("RatePerPath_L1", "L1 rate per path", hltConfig_.size(), -0.5, hltConfig_.size()-0.5); 
      hRatePerPath_L2 = new TH1D("RatePerPath_L2", "L2 rate per path", hltConfig_.size(), -0.5, hltConfig_.size()-0.5); 
      hRatePerPath_L3 = new TH1D("RatePerPath_L3", "L3 rate per path", hltConfig_.size(), -0.5, hltConfig_.size()-0.5); 
      hRatePerPath_L3iso = new TH1D("RatePerPath_L3iso", "L3 isolation rate per path", hltConfig_.size(), -0.5, hltConfig_.size()-0.5); 
    }
    else {
      hRatePerPath = new TH1D("RatePerPath", "Rate per path", listOfPaths.size(), -0.5, listOfPaths.size()-0.5); 
      hRatePerPath_L1 = new TH1D("RatePerPath_L1", "L1 rate per path", listOfPaths.size(), -0.5, listOfPaths.size()-0.5); 
      hRatePerPath_L2 = new TH1D("RatePerPath_L2", "L2 rate per path", listOfPaths.size(), -0.5, listOfPaths.size()-0.5); 
      hRatePerPath_L3 = new TH1D("RatePerPath_L3", "L3 rate per path", listOfPaths.size(), -0.5, listOfPaths.size()-0.5); 
      hRatePerPath_L3iso = new TH1D("RatePerPath_L3iso", "L3 isolation rate per path", listOfPaths.size(), -0.5, listOfPaths.size()-0.5); 
    }
    hRatePerPath->Sumw2();
    hRatePerPath_L1->Sumw2();
    hRatePerPath_L2->Sumw2();
    hRatePerPath_L3->Sumw2();
    hRatePerPath_L3iso->Sumw2();

    hRate_OR = new TProfile( "TotalRate", "Total rate", nRateBins, firstRateBin, lastRateBin); 

    std::cout << "===== Initializing counters to zero ====================" << std::endl;
    countsAllPathByLS_OR = 0;
    unsigned int hltCnt(0);
    for(unsigned iHltPath=0; iHltPath<hltConfig_.size(); ++iHltPath) {
      if( !checkTriggerName( hltConfig_.triggerName(iHltPath) ) ) continue; 

      countsPerPath.push_back(0);
      countsPerPathByLS_EX.push_back(0);
      countsPerPathByLS.push_back(0);
      countsPerPathByLS_L1.push_back(0);
      countsPerPathByLS_L2.push_back(0);
      countsPerPathByLS_L3.push_back(0);
      countsPerPathByLS_L3iso.push_back(0);

      averagePrescale.push_back(0);
      averagePrescale_L1.push_back(0);
      averagePrescale_L2.push_back(0);
      averagePrescale_L3.push_back(0);
      averagePrescale_L3iso.push_back(0);

      hRatePerPath->GetXaxis()->SetBinLabel(hltCnt+1, hltConfig_.triggerName(iHltPath).c_str()); 
      hRatePerPath_L1->GetXaxis()->SetBinLabel(hltCnt+1, hltConfig_.triggerName(iHltPath).c_str()); 
      hRatePerPath_L2->GetXaxis()->SetBinLabel(hltCnt+1, hltConfig_.triggerName(iHltPath).c_str()); 
      hRatePerPath_L3->GetXaxis()->SetBinLabel(hltCnt+1, hltConfig_.triggerName(iHltPath).c_str()); 
      hRatePerPath_L3iso->GetXaxis()->SetBinLabel(hltCnt+1, hltConfig_.triggerName(iHltPath).c_str()); 

      std::string plotNameBase = "rate_"+hltConfig_.triggerName(iHltPath);
      hRates_EX.push_back(    new TProfile( ("EX_"   +plotNameBase).c_str(), ("EX_"   +plotNameBase).c_str(), nRateBins, firstRateBin, lastRateBin) ); 
      hRates.push_back(       new TProfile( (         plotNameBase).c_str(), (         plotNameBase.c_str()), nRateBins, firstRateBin, lastRateBin) ); 
      hRates_L1.push_back(    new TProfile( ("L1_"   +plotNameBase).c_str(), ("L1_"   +plotNameBase).c_str(), nRateBins, firstRateBin, lastRateBin) );
      hRates_L2.push_back(    new TProfile( ("L2_"   +plotNameBase).c_str(), ("L2_"   +plotNameBase).c_str(), nRateBins, firstRateBin, lastRateBin) );
      hRates_L3.push_back(    new TProfile( ("L3_"   +plotNameBase).c_str(), ("L3_"   +plotNameBase).c_str(), nRateBins, firstRateBin, lastRateBin) );
      hRates_L3iso.push_back( new TProfile( ("L3iso_"+plotNameBase).c_str(), ("L3iso_"+plotNameBase).c_str(), nRateBins, firstRateBin, lastRateBin) );
      hltCnt++;
    }
    isFirstRun = false;
  }

  ///////////////////////////////// BOOK PLOTS ///////////////////////////////
  
} // end of beginRun


// ------------ method called when ending the processing of a run  ------------
void HLTCounter::endRun(edm::Run const&, edm::EventSetup const&) {
  using namespace std;

  if(debugPrint) std::cout << "Inside end run, print = " << printAtEndRun << std::endl;
  if(!printAtEndRun) return; 
  if(debugPrint) std::cout << "End of run printing summary table " << std::endl;
  cout << "################### TRIGGER SUMMARY ####################### "
       << endl << endl;

  unsigned int hltCnt(0);
  for(unsigned iHltPath=0; iHltPath<hltConfig_.size(); ++iHltPath) {
    if( !checkTriggerName( hltConfig_.triggerName(iHltPath) ) ) continue; 
    string tempName = hltConfig_.removeVersion(hltConfig_.triggerName(iHltPath));

    // cut the string down to size
    // then pad it out with . if it was smaller
    if(tempName.size()>52) {
      tempName.resize(52);
      tempName.resize(55, '.');
    }
    cout << left << setw(58) << tempName
	 << right << countsPerPath[hltCnt] << endl;
    hltCnt++;
  }
}


// ------------ method called when starting to processes a luminosity block  ------------
void HLTCounter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
  TotNumLumiSections++;
}

// ------------ method called when ending the processing of a luminosity block  ------------
void HLTCounter::endLuminosityBlock(edm::LuminosityBlock const& iLumiBlock, edm::EventSetup const&) {

  int lumi = int(iLumiBlock.id().luminosityBlock());

  if(debugPrint) std::cout << "In LS " << lumi << ", NumInstLumis = " << NumInstLumis-NumInstLumisPrevLS << std::endl;
  for( int n=NumInstLumisPrevLS; n<NumInstLumis; ++n) {
    if(debugPrint) std::cout << "InstLumiArray[" << n << "] = " << InstLumiArray[n] << std::endl;
  }

  if( TotNumLumiSections%nLumisInAverage==0 ) {
    hRate_OR->Fill( averageInstLumi, countsAllPathByLS_OR / (nLumisInAverage*23.31) ); 

    unsigned int hltCnt(0);
    for(unsigned iHltPath=0; iHltPath<hltConfig_.size(); ++iHltPath) {
      if( !checkTriggerName( hltConfig_.triggerName(iHltPath) ) ) continue; 
      hRates_EX[hltCnt]   ->Fill( averageInstLumi, countsPerPathByLS_EX[hltCnt]                                    / (nLumisInAverage*23.31) ); 
      hRates[hltCnt]      ->Fill( averageInstLumi, countsPerPathByLS[hltCnt]       * averagePrescale[hltCnt]       / (nLumisInAverage*23.31) ); 
      hRates_L1[hltCnt]   ->Fill( averageInstLumi, countsPerPathByLS_L1[hltCnt]    * averagePrescale_L1[hltCnt]    / (nLumisInAverage*23.31) ); 
      hRates_L2[hltCnt]   ->Fill( averageInstLumi, countsPerPathByLS_L2[hltCnt]    * averagePrescale_L2[hltCnt]    / (nLumisInAverage*23.31) ); 
      hRates_L3[hltCnt]   ->Fill( averageInstLumi, countsPerPathByLS_L3[hltCnt]    * averagePrescale_L3[hltCnt]    / (nLumisInAverage*23.31) ); 
      hRates_L3iso[hltCnt]->Fill( averageInstLumi, countsPerPathByLS_L3iso[hltCnt] * averagePrescale_L3iso[hltCnt] / (nLumisInAverage*23.31) ); 
      hltCnt++;
    }

    clearLumiAverage();
  }
}

void HLTCounter::addLumiToAverage(double lumi, int nLumi) {

  averageInstLumi = (averageInstLumi*nLumi + lumi) / (nLumi+1);
  InstLumiArray.push_back(lumi);
  NumInstLumis = nLumi + 1;

  return;
}

void HLTCounter::clearLumiAverage() {

  averageInstLumi = 0;
  NumInstLumis = 0;
  InstLumiArray.clear();
  countsAllPathByLS_OR = 0;
  for(unsigned int h=0; h<countsPerPathByLS.size(); ++h) {
    countsPerPathByLS_EX[h] = 0;
    countsPerPathByLS[h] = 0;
    countsPerPathByLS_L1[h] = 0;
    countsPerPathByLS_L2[h] = 0;
    countsPerPathByLS_L3[h] = 0;
    countsPerPathByLS_L3iso[h] = 0;

    averagePrescale[h] = 0.;
    averagePrescale_L1[h] = 0.;
    averagePrescale_L2[h] = 0.;
    averagePrescale_L3[h] = 0.;
    averagePrescale_L3iso[h] = 0.;
  }

  return;  
}

void HLTCounter::addPrescaleToAverage(double evtPs, unsigned int pos, vector<double>& avgPs, vector<unsigned>& cntPerLS) {
  avgPs[pos] = ( avgPs[pos] * (cntPerLS[pos]-1) + evtPs ) / cntPerLS[pos]; 
  return; 
}

bool HLTCounter::checkTriggerName(std::string const thisTrigName) {
  if(listOfPaths.size()==0) return true; 
  for(unsigned int i=0; i<listOfPaths.size(); ++i) {
    if( thisTrigName.find( listOfPaths[i] ) != std::string::npos ) return true;
  }
  return false;
}


// define this as a plug-in
DEFINE_FWK_MODULE(HLTCounter);
