#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTEventAnalyzerAOD.h"

#include <map>
#include <string>
#include <memory>
#include <iomanip>

#include "TH1F.h"
#include "TH2F.h"


/**
   Class declaration

   \class MuonTriggerEfficiencyAnalyzer MuonTriggerEfficiencyAnalyzer.h "PhysicsTools/UtilAlgos/interface/MuonTriggerEfficiencyAnalyzer.h"
   \brief Example class that can be used both within FWLite and within the full framework

   This is an example for keeping classes that can be used both within FWLite and within the full 
   framework. The class is derived from the BasicAnalyzer base class, which is an interface for 
   the two wrapper classes EDAnalyzerWrapper and FWLiteAnalyzerWrapper. The latter provides basic 
   configuration file reading and event looping equivalent to the FWLiteHistograms executable of 
   this package. You can see the FWLiteAnalyzerWrapper class at work in the FWLiteWithBasicAnalyzer
   executable of this package.
*/

class MuonTriggerEfficiencyAnalyzer : public edm::EDAnalyzer {

 public:
  /// default constructor
  MuonTriggerEfficiencyAnalyzer(const edm::ParameterSet& cfg);
  /// default destructor
  virtual ~MuonTriggerEfficiencyAnalyzer() {};
  /// everything that needs to be done before the event loop

 private:
  virtual void beginJob();
  /// everything that needs to be done after the event loop
  virtual void endJob();
  /// everything that needs to be done before each run
  virtual void beginRun(const edm::Run & run, const edm::EventSetup & eventSetup);
  /// everything that needs to be done after each run
  virtual void endRun(const edm::Run & run, const edm::EventSetup & eventSetup);
  /// everything that needs to be done during the event loop
  virtual void analyze(const edm::Event& event, const edm::EventSetup& eventSetup);

  /// input tag for mouns
  edm::InputTag vertexes_;
  /// input tag for mouns
  edm::InputTag muons_;
  /// file service
  edm::Service<TFileService> outfile_;
  /// histograms
  std::map<std::string, TH1*> hists_;

  // Trigger process
  std::string triggerProcess_;
  // Trigger names
  std::string tagTriggerName_;
  std::string triggerName_;
  // Trigger indexes
  int tagTriggerIndex_;
  int triggerIndex_;
  // HLTConfig
  HLTConfigProvider hltConfig_;

  /// Max number of offline muons allowed in the event
  unsigned int nMaxMuons_;
};

/// default constructor
MuonTriggerEfficiencyAnalyzer::MuonTriggerEfficiencyAnalyzer(const edm::ParameterSet& cfg): 
  vertexes_(cfg.getParameter<edm::InputTag>("vertexes")), 
  muons_(cfg.getParameter<edm::InputTag>("muons")), 
  triggerProcess_(cfg.getParameter<std::string>("triggerProcess")), 
  tagTriggerName_(cfg.getParameter<std::string>("tagTriggerName")), 
  triggerName_(cfg.getParameter<std::string>("triggerName")), 
  nMaxMuons_(cfg.getUntrackedParameter<unsigned int>("maxNumberMuons", 2))
{}

void MuonTriggerEfficiencyAnalyzer::beginJob() {

  double eta_bins[] = {-2.4, -2.1, -1.6, -1.2, -0.9, -0.3, -0.2, 0.2, 0.3, 0.9, 1.2, 1.6, 2.1, 2.4}; 
  int eta_bin_n = sizeof(eta_bins)/sizeof(double); 

  hists_["muonPt_tag"]   = outfile_->make<TH1F>("muonPt_tag"  , "pt"  ,  100,  0., 300.);
  hists_["muonEta_tag"]  = outfile_->make<TH1F>("muonEta_tag" , "eta" ,  eta_bin_n-1, eta_bins);
  hists_["muonPhi_tag"]  = outfile_->make<TH1F>("muonPhi_tag" , "phi" ,  100, -5.,   5.); 
  hists_["muonNvtx_tag"] = outfile_->make<TH1F>("nvtx_tag",     "nvtx", 30, 0.5, 30.5);

  hists_["muonPt_tag"]->Sumw2();
  hists_["muonEta_tag"]->Sumw2();
  hists_["muonPhi_tag"]->Sumw2();
  hists_["muonNvtx_tag"]->Sumw2();

  hists_["muonPt_probe_den"]   = outfile_->make<TH1F>("muonPt_probe_den"  , "pt"  ,  100,  0., 300.);
  hists_["muonEta_probe_den"]  = outfile_->make<TH1F>("muonEta_probe_den" , "eta" ,  eta_bin_n-1, eta_bins);
  hists_["muonPhi_probe_den"]  = outfile_->make<TH1F>("muonPhi_probe_den" , "phi" ,  100, -5.,   5.); 
  hists_["muonNvtx_probe_den"] = outfile_->make<TH1F>("muonNvtx_probe_den", "nvtx", 30, 0.5, 30.5);

  hists_["muonPt_probe_den"]->Sumw2();
  hists_["muonEta_probe_den"]->Sumw2();
  hists_["muonPhi_probe_den"]->Sumw2();
  hists_["muonNvtx_probe_den"]->Sumw2();

  hists_["muonPt_probe_num"]   = outfile_->make<TH1F>("muonPt_probe_num"  , "pt"  ,  100,  0., 300.);
  hists_["muonEta_probe_num"]  = outfile_->make<TH1F>("muonEta_probe_num" , "eta" ,  eta_bin_n-1, eta_bins);
  hists_["muonPhi_probe_num"]  = outfile_->make<TH1F>("muonPhi_probe_num" , "phi" ,  100, -5.,   5.); 
  hists_["muonNvtx_probe_num"] = outfile_->make<TH1F>("muonNvtx_probe_num", "nvtx", 30, 0.5, 30.5);

  hists_["muonPt_probe_num"]->Sumw2();
  hists_["muonEta_probe_num"]->Sumw2();
  hists_["muonPhi_probe_num"]->Sumw2();
  hists_["muonNvtx_probe_num"]->Sumw2();

  hists_["muonPt_probe"]   = outfile_->make<TH1F>("muonPt_probe"  , "pt"  ,  100,  0., 300.);
  hists_["muonEta_probe"]  = outfile_->make<TH1F>("muonEta_probe" , "eta" ,  eta_bin_n-1, eta_bins);
  hists_["muonPhi_probe"]  = outfile_->make<TH1F>("muonPhi_probe" , "phi" ,  100, -5.,   5.); 
  hists_["muonNvtx_probe"] = outfile_->make<TH1F>("muonNvtx_probe", "nvtx", 30, 0.5, 30.5);

  hists_["muonPt_probe"]->Sumw2();
  hists_["muonEta_probe"]->Sumw2();
  hists_["muonPhi_probe"]->Sumw2();
  hists_["muonNvtx_probe"]->Sumw2();

  hists_["mumuMass_all"] = outfile_->make<TH1F>("mumuMass_all", "mass",   90, 30., 120.);
  hists_["mumuMass_den"] = outfile_->make<TH1F>("mumuMass_den", "mass",   90, 30., 120.);
  hists_["mumuMass_num"] = outfile_->make<TH1F>("mumuMass_num", "mass",   90, 30., 120.);
  hists_["mumuMass"]     = outfile_->make<TH1F>("mumuMass",     "mass",   90, 30., 120.);

  hists_["mumuMass_all"]->Sumw2();
  hists_["mumuMass_den"]->Sumw2();
  hists_["mumuMass_num"]->Sumw2();
  hists_["mumuMass"]->Sumw2();

  hists_["muonPt12_den"]  = outfile_->make<TH2F>("muonPt12_den"  , "p_{T,1} vs p_{T,2}",  80, 20., 100., 80, 20., 100.);
  hists_["muonEta12_den"] = outfile_->make<TH2F>("muonEta12_den" , "eta" , eta_bin_n-1, eta_bins, eta_bin_n-1, eta_bins);
  hists_["muonPhi12_den"] = outfile_->make<TH2F>("muonPhi12_den" , "phi" , 16, -3.2, 3.2, 16, -3.2, 3.2); 

  hists_["muonPt12_den"]->Sumw2();
  hists_["muonEta12_den"]->Sumw2();
  hists_["muonPhi12_den"]->Sumw2();

  hists_["muonPt12_num"]  = outfile_->make<TH2F>("muonPt12_num"  , "p_{T,1} vs p_{T,2}",  80, 20., 100., 80, 20., 100.);
  hists_["muonEta12_num"] = outfile_->make<TH2F>("muonEta12_num" , "eta" , eta_bin_n-1, eta_bins, eta_bin_n-1, eta_bins);
  hists_["muonPhi12_num"] = outfile_->make<TH2F>("muonPhi12_num" , "phi" , 16, -3.2, 3.2, 16, -3.2, 3.2); 

  hists_["muonPt12_num"]->Sumw2();
  hists_["muonEta12_num"]->Sumw2();
  hists_["muonPhi12_num"]->Sumw2();

  hists_["muonPt12"]  = outfile_->make<TH2F>("muonPt12"  , "p_{T,1} vs p_{T,2}",  80, 20., 100., 80, 20., 100.);
  hists_["muonEta12"] = outfile_->make<TH2F>("muonEta12" , "eta" , eta_bin_n-1, eta_bins, eta_bin_n-1, eta_bins);
  hists_["muonPhi12"] = outfile_->make<TH2F>("muonPhi12" , "phi" , 16, -3.2, 3.2, 16, -3.2, 3.2); 

  hists_["muonPt12"]->Sumw2();
  hists_["muonEta12"]->Sumw2();
  hists_["muonPhi12"]->Sumw2();

  hists_["deltaR_trobj_tag"]   = outfile_->make<TH1F>("deltaR_trobj_tag" ,   "#DeltaR(trig,#mu)" , 150, 0., 0.3); 
  hists_["deltaR_trobj_probe"] = outfile_->make<TH1F>("deltaR_trobj_probe" , "#DeltaR(trig,#mu)" , 150, 0., 0.3); 
}

void MuonTriggerEfficiencyAnalyzer::endJob() {
  hists_["muonPt_probe"]->Divide(   hists_["muonPt_probe_num"],   hists_["muonPt_probe_den"],   1.0, 1.0, "B" );
  hists_["muonEta_probe"]->Divide(  hists_["muonEta_probe_num"],  hists_["muonEta_probe_den"],  1.0, 1.0, "B" );
  hists_["muonPhi_probe"]->Divide(  hists_["muonPhi_probe_num"],  hists_["muonPhi_probe_den"],  1.0, 1.0, "B" );
  hists_["muonNvtx_probe"]->Divide( hists_["muonNvtx_probe_num"], hists_["muonNvtx_probe_den"], 1.0, 1.0, "B" ); 

  hists_["mumuMass"]->Divide(  hists_["mumuMass_num"],  hists_["mumuMass_den"],  1.0, 1.0, "B" ); 

  hists_["muonPt12"]->Divide(  hists_["muonPt12_num"],  hists_["muonPt12_den"],  1.0, 1.0, "B" ); 
  hists_["muonEta12"]->Divide( hists_["muonEta12_num"], hists_["muonEta12_den"], 1.0, 1.0, "B" ); 
  hists_["muonPhi12"]->Divide( hists_["muonPhi12_num"], hists_["muonPhi12_den"], 1.0, 1.0, "B" ); 
}

void MuonTriggerEfficiencyAnalyzer::beginRun(const edm::Run & run, const edm::EventSetup & eventSetup) {
  bool changed = true;
  if( hltConfig_.init(run, eventSetup, triggerProcess_, changed) ) {
  }
  else {
    std::cout << "Warning, didn't find process " << triggerProcess_.c_str() << std::endl;
    // Now crash
    assert(false);
  }

  triggerIndex_ = -1; 
  tagTriggerIndex_ = -1; 

  for(unsigned iHltPath=0; iHltPath<hltConfig_.size(); ++iHltPath) {
    std::string tempName = hltConfig_.triggerName(iHltPath);
    if(tempName.find(triggerName_) != std::string::npos) {
      triggerIndex_ = int(iHltPath);
    }
    if(tempName.find(tagTriggerName_) != std::string::npos) {
      tagTriggerIndex_ = int(iHltPath);
    }

    if( triggerIndex_>-1 && tagTriggerIndex_>-1 ) break; 
  } // end for each path

  if( triggerIndex_ == -1 ) {
    std::cout << "Warning, didn't find trigger " <<  triggerName_.c_str() << std::endl;
    // Now crash
    assert(false);    
  }
  if( tagTriggerIndex_ == -1 ) {
    std::cout << "Warning, didn't find tag trigger " <<  tagTriggerName_.c_str() << std::endl;
    // Now crash
    assert(false);    
  }
}

void MuonTriggerEfficiencyAnalyzer::endRun(const edm::Run & run, const edm::EventSetup & eventSetup) {}
 
void MuonTriggerEfficiencyAnalyzer::analyze(const edm::Event &event, const edm::EventSetup &eventSetup) {
  // define what muon you are using; this is necessary as FWLite is not 
  // capable of reading edm::Views
  using reco::Muon;

  edm::Handle<reco::VertexCollection> pvHandle; 
  event.getByLabel(vertexes_, pvHandle);
  const reco::VertexCollection & vertices = *pvHandle.product();
  unsigned int nGoodVtx = 0; 
  for(reco::VertexCollection::const_iterator it=vertices.begin(); it!=vertices.end(); ++it) {
    if( it->ndof()>4                     && 
	(fabs(it->z())<=24.)             && 
	(fabs(it->position().rho())<=2.)   ) 
      nGoodVtx++;
  }
  if( nGoodVtx==0 ) return;
  const reco::Vertex & pv = vertices[0];

  // Handle to the muon collection
  edm::Handle<std::vector<Muon> > muons;
  event.getByLabel(muons_, muons);

  if( nMaxMuons_>0 && muons->size()>nMaxMuons_ ) return; 

  // Get trigger results
  edm::Handle<edm::TriggerResults> triggerResults;
  event.getByLabel(edm::InputTag("TriggerResults", "", triggerProcess_), triggerResults);

  if(!triggerResults.isValid()) {
    std::cout << "Trigger results not valid" << std::endl;
    return;
  } 

  if( !triggerResults->accept(tagTriggerIndex_) ) return; // there are no tags

  // Get trigger summary 
  edm::Handle<trigger::TriggerEvent> triggerEvent;
  event.getByLabel(edm::InputTag("hltTriggerSummaryAOD", "", triggerProcess_), triggerEvent);

  if(!triggerEvent.isValid()) { 
    std::cout << "TriggerEvent not valid" << std::endl;
    return;
  }

  // Sanity check
  assert(triggerResults->size()==hltConfig_.size());

  // Get trigger objects from trigger summary
  const trigger::TriggerObjectCollection & toc = triggerEvent->getObjects();

  // Modules in tag trigger path
  const std::vector<std::string>& tagModuleLabels(hltConfig_.moduleLabels(tagTriggerIndex_));
  assert( tagModuleLabels.size()==hltConfig_.size(tagTriggerIndex_) );
  const unsigned int tagModuleIndex( hltConfig_.size(tagTriggerIndex_)-2 ); // index of last filter (excluding HLTEndBool)
  const unsigned int tagFilterIndex( triggerEvent->filterIndex( edm::InputTag( tagModuleLabels[tagModuleIndex], "", triggerProcess_) ) );
  assert( tagFilterIndex < triggerEvent->sizeFilters() );
  const trigger::Vids & tagVids( triggerEvent->filterIds(tagFilterIndex) );
  const trigger::Keys & tagKeys( triggerEvent->filterKeys(tagFilterIndex) );
  assert( tagVids.size()==tagKeys.size() );
  const unsigned int nTagTrig(tagVids.size());

  // Loop muon collection and fill histograms
  for(std::vector<Muon>::const_iterator mu1=muons->begin(); mu1!=muons->end(); ++mu1) {
    if( muon::isTightMuon( (*mu1), pv ) && (*mu1).pt()>20 ) {

      // Is this a suitable tag? 
      double maxTagDeltaR = 0.1; 
      double finTagDeltaR = 10.0; 
      bool isTagTrigMatch = false; 
      for(unsigned int i=0; i!=nTagTrig; ++i) {
	const trigger::TriggerObject & tagTo = toc[tagKeys[i]];
	double tmpTagDeltaR = deltaR( (*mu1), tagTo ); 
	if( tmpTagDeltaR<maxTagDeltaR ) {
	  isTagTrigMatch = true;
	  if( tmpTagDeltaR<finTagDeltaR ) {
	    finTagDeltaR = tmpTagDeltaR;
	  }
	  //break;
	}
      }

      hists_["deltaR_trobj_tag"]->Fill(finTagDeltaR); 

      if(isTagTrigMatch==false) continue; 

      // Go on and look for a probe
      for(std::vector<Muon>::const_iterator mu2=muons->begin(); mu2!=muons->end(); ++mu2) {
	if( mu2==mu1 ) continue; 
	if( muon::isTightMuon( (*mu2), pv ) && (*mu2).pt()>20 ) {
	  if( mu1->charge()*mu2->charge()<0 ) { // check only muon pairs of unequal charge 

	    double mumuMass = (mu1->p4()+mu2->p4()).mass();
	    hists_["mumuMass_all"]->Fill( mumuMass );

	    if( mumuMass>86. && mumuMass<96. ) {// check only muon pairs compatible with Z decay
	      // Ok, tag-probe pair found
	      // Go on and fill denominator plots
	      hists_["muonPt_tag" ]->Fill( mu1->pt () );
	      hists_["muonEta_tag"]->Fill( mu1->eta() );
	      hists_["muonPhi_tag"]->Fill( mu1->phi() );
	      hists_["muonNvtx_tag"]->Fill( nGoodVtx );

	      hists_["muonPt_probe_den" ]->Fill( mu2->pt () );
	      hists_["muonEta_probe_den"]->Fill( mu2->eta() );
	      hists_["muonPhi_probe_den"]->Fill( mu2->phi() );
	      hists_["muonNvtx_probe_den"]->Fill( nGoodVtx );

	      hists_["mumuMass_den"]->Fill( mumuMass );

	      hists_["muonPt12_den"]->Fill( mu1->pt(), mu2->pt() );
	      hists_["muonEta12_den"]->Fill( mu1->eta(), mu2->eta() );
	      hists_["muonPhi12_den"]->Fill( mu1->phi(), mu2->phi() );
  
	      // Has the probe fired the trigger?
	      if( !triggerResults->accept(triggerIndex_) ) continue; // no muon has fired the target trigger

	      // Modules in probe trigger path
	      const std::vector<std::string>& moduleLabels(hltConfig_.moduleLabels(triggerIndex_));
	      assert( moduleLabels.size()==hltConfig_.size(triggerIndex_) );
	      const unsigned int moduleIndex( hltConfig_.size(triggerIndex_)-2 ); // index of last filter (excluding HLTEndBool)
	      const unsigned int filterIndex( triggerEvent->filterIndex( edm::InputTag( moduleLabels[moduleIndex], "", triggerProcess_) ) );
	      assert( filterIndex < triggerEvent->sizeFilters() );
	      const trigger::Vids & vids( triggerEvent->filterIds(filterIndex) );
	      const trigger::Keys & keys( triggerEvent->filterKeys(filterIndex) );
	      assert( vids.size()==keys.size() );
	      const unsigned int nProbeTrig(vids.size());
 
	      double maxProbeDeltaR = 0.1; 
	      double finProbeDeltaR = 10.0; 
	      bool isProbeTrigMatch = false; 
	      for(unsigned int i=0; i!=nProbeTrig; ++i) {
		const trigger::TriggerObject & to = toc[keys[i]];
		double tmpProbeDeltaR = deltaR( (*mu2), to ); 
		if( tmpProbeDeltaR<maxProbeDeltaR ) {
		  isProbeTrigMatch = true;
		  if( tmpProbeDeltaR<finProbeDeltaR ) {
		    finProbeDeltaR = tmpProbeDeltaR;
		  }
		  //break;
		}
	      }

	      hists_["deltaR_trobj_probe"]->Fill(finProbeDeltaR); 

	      if(isProbeTrigMatch==false) continue; 

	      // Ok, probe passed
	      // Go on and fill numerator plots
	      hists_["muonPt_probe_num" ]->Fill( mu2->pt () );
	      hists_["muonEta_probe_num"]->Fill( mu2->eta() );
	      hists_["muonPhi_probe_num"]->Fill( mu2->phi() );
	      hists_["muonNvtx_probe_num"]->Fill( nGoodVtx );

	      hists_["mumuMass_num"]->Fill( mumuMass );

	      hists_["muonPt12_num"]->Fill( mu1->pt(), mu2->pt() );
	      hists_["muonEta12_num"]->Fill( mu1->eta(), mu2->eta() );
	      hists_["muonPhi12_num"]->Fill( mu1->phi(), mu2->phi() );
	    } 
	  }
	}
      }
    }
  }
}

// define this as a plug-in
DEFINE_FWK_MODULE(MuonTriggerEfficiencyAnalyzer);
