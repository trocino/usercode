#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/PatCandidates/interface/EventHypothesis.h"
#include "DataFormats/PatCandidates/interface/EventHypothesisLooper.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CMGTools/HtoZZ2l2nu/interface/ObjectFilters.h"
#include "CMGTools/HtoZZ2l2nu/interface/setStyle.h"

#include "TH1D.h"
#include "TString.h"
 

class DileptonPlusMETEventProducer : public edm::EDProducer {
public:
  DileptonPlusMETEventProducer(const edm::ParameterSet &iConfig) ;
  virtual void produce( edm::Event &iEvent, const edm::EventSetup &iSetup) ;
private:
  std::map<std::string, edm::ParameterSet > objConfig;
  std::map<TString,TH1D *> controlHistos_;
};


using namespace std;


//
DileptonPlusMETEventProducer::DileptonPlusMETEventProducer(const edm::ParameterSet &iConfig) {
  produces<std::vector<pat::EventHypothesis> >("selectedEvent");
  produces<reco::VertexCollection>("selectedVertices");
  produces<reco::MuonCollection>("selectedTracks");
  produces<std::vector<int> >("selectionInfo");
  std::string objs[]={"Generator", "Vertices", "Electrons", "Muons", "Dileptons", "Tracks", "Jets", "MET" };
  for(size_t iobj=0; iobj<sizeof(objs)/sizeof(string); iobj++)
    objConfig[ objs[iobj] ] = iConfig.getParameter<edm::ParameterSet>( objs[iobj] );

  edm::Service<TFileService> fs;
  TString cats[]={"electron","muon"};
  size_t ncats=sizeof(cats)/sizeof(TString);
  for(size_t icat=0; icat<ncats; icat++) {
    TFileDirectory newDir=fs->mkdir(cats[icat].Data());
    controlHistos_[cats[icat]+"_rho"] = (TH1D *) formatPlot( newDir.make<TH1F>(cats[icat]+"_rho", "; #rho; Events", 100, 0.,10.), 1,1,1,20,0,false,true,1,1,1 );
    controlHistos_[cats[icat]+"_ecaliso"] = (TH1D *) formatPlot( newDir.make<TH1F>(cats[icat]+"_ecaliso", ";ECAL isolation; Events", 100, 0.,10.), 1,1,1,20,0,false,true,1,1,1 );
    controlHistos_[cats[icat]+"_hcaliso"] = (TH1D *) formatPlot( newDir.make<TH1F>(cats[icat]+"_hcaliso", ";HCAL isolation; Events", 100, 0.,10.), 1,1,1,20,0,false,true,1,1,1 );
    controlHistos_[cats[icat]+"_caloiso"] = (TH1D *) formatPlot( newDir.make<TH1F>(cats[icat]+"_caloiso", ";Calorimeter isolation; Events", 100, 0.,10.), 1,1,1,20,0,false,true,1,1,1 );
    controlHistos_[cats[icat]+"_trackiso"] = (TH1D *) formatPlot( newDir.make<TH1F>(cats[icat]+"_trackiso", ";Tracker Isolation; Events", 100, 0.,10.), 1,1,1,20,0,false,true,1,1,1 );
    controlHistos_[cats[icat]+"_reliso"] = (TH1D *) formatPlot( newDir.make<TH1F>(cats[icat]+"_reliso", "; Isolation; Events", 100, 0.,10.), 1,1,1,20,0,false,true,1,1,1 );
  }
}

//
void DileptonPlusMETEventProducer::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  using namespace std;
  using namespace edm;
  using namespace pat::eventhypothesis;
  using reco::Candidate; 
  using reco::CandidatePtr;


  pat::EventHypothesis hyp;
  int selStep(0), selPath(0);
  
  // Pre-select vertices
  Handle<reco::VertexCollection> hVtx;
  iEvent.getByLabel(objConfig["Vertices"].getParameter<edm::InputTag>("source"), hVtx);  
  std::vector<reco::VertexRef> selVertices=vertex::filter(hVtx, objConfig["Vertices"]);
  const reco::Vertex *theSelVertex=0;
  if(selVertices.size()>0) selStep=1;
  else return;

  // Average energy density
  edm::Handle< double > rho;
  iEvent.getByLabel(edm::InputTag("kt6PFJets:rho"), rho);

  // Select muons (id + very loose uncorrected isolation)
  Handle<View<Candidate> > hMu; 
  iEvent.getByLabel(objConfig["Muons"].getParameter<edm::InputTag>("source"), hMu);
  CandidateWithVertexCollection selMuons=muon::filter(hMu, selVertices, objConfig["Muons"], iSetup);

  // Select electrons (id + conversion veto + very loose uncorrected isolation)
  Handle<View<Candidate> > hEle; 
  iEvent.getByLabel(objConfig["Electrons"].getParameter<edm::InputTag>("source"), hEle);
  CandidateWithVertexCollection selElectrons=electron::filter(hEle, hMu, selVertices, objConfig["Electrons"], iSetup);
  
  // Build inclusive collection
  CandidateWithVertexCollection selLeptons=selMuons;
  selLeptons.insert(selLeptons.end(), selElectrons.begin(), selElectrons.end());
  if(selLeptons.size()>0) selStep=2;

  // Build the dilepton (all tightly isolated leptons will be returned)
  CandidateWithVertexCollection isolLeptons;

  DileptonWithVertex dileptonWithVertex = dilepton::filter(selLeptons,
							   selVertices,
							   objConfig["Dileptons"],
							   iSetup,
							   *rho,
							   isolLeptons,
							   &controlHistos_);
  selPath=dilepton::classify(dileptonWithVertex.second);
  if(selPath>0) {
    selStep=3;

    std::vector<CandidatePtr> &dilepton=dileptonWithVertex.second;
    hyp.add(dilepton[0],"leg1");
    hyp.add(dilepton[1],"leg2");

    // Dilepton vertex
    theSelVertex=dileptonWithVertex.first.get();
    if(theSelVertex==0) return;

    //add the remaining isolated leptons now
    bool onlyMuFromPV=( objConfig["Muons"].existsAs<bool>("OnlyFromPV") ? 
			objConfig["Muons"].getParameter<bool>("OnlyFromPV") : false ); 
    bool onlyEFromPV=( objConfig["Electrons"].existsAs<bool>("OnlyFromPV") ? 
		       objConfig["Electrons"].getParameter<bool>("OnlyFromPV") : false ); 
    for(CandidateWithVertexCollection::iterator lIt=isolLeptons.begin(); 
	lIt!=isolLeptons.end(); ++lIt) {
      reco::VertexRef aVtx=(*lIt).first;
      CandidatePtr iLept=(*lIt).second;
      // Skip selected dilepton
      if(iLept.get()==dilepton[0].get() || iLept.get()==dilepton[1].get()) continue;

      string flav;
      bool skipIfNotSameVtx;
      if( fabs(dilepton::getLeptonId(iLept))==13 ) {
	flav="muon";
	skipIfNotSameVtx=onlyMuFromPV;
      }
      else {
	flav="electron";
	skipIfNotSameVtx=onlyEFromPV;
      }

      // Skip if not from PV
      if(skipIfNotSameVtx && aVtx.get()!=theSelVertex) continue;

      // Save lepton in event hypothesis
      hyp.add(iLept, flav.c_str());
    }
  }

  // Select all (tight-)isolated tracks (from generalTracks)
  Handle<std::vector<reco::Track> > hTrks; 
  iEvent.getByLabel(objConfig["Tracks"].getParameter<edm::InputTag>("source"), hTrks);
  TrackWithVertexCollection selTracks=track::filter(hTrks, 
						    isolLeptons, 
						    selVertices, 
						    objConfig["Tracks"], 
						    iSetup); 
  reco::MuonCollection selVtxTracks; 
  bool onlyTrackFromPV=( objConfig["Tracks"].existsAs<bool>("OnlyFromPV") ? 
			 objConfig["Tracks"].getParameter<bool>("OnlyFromPV") : false );
  for(TrackWithVertexCollection::iterator itTk=selTracks.begin(); itTk!=selTracks.end(); ++itTk) {
    if(onlyTrackFromPV && (*itTk).first.get()!=theSelVertex) continue;
    selVtxTracks.push_back( (*itTk).second );
  }

  // Add the jets
  Handle<View<Candidate> > hJet; 
  iEvent.getByLabel(objConfig["Jets"].getParameter<edm::InputTag>("source"), hJet);
  CandidateWithVertexCollection selJets=jet::filter(hJet, 
						    isolLeptons, 
						    selVertices, 
						    objConfig["Jets"]); 
  bool onlyJetFromPV=( objConfig["Jets"].existsAs<bool>("OnlyFromPV") ? 
		       objConfig["Jets"].getParameter<bool>("OnlyFromPV") : false );
  for(CandidateWithVertexCollection::iterator jIt=selJets.begin(); jIt!=selJets.end(); ++jIt) {
    if(onlyJetFromPV && (*jIt).first.get()!=theSelVertex) continue;
    hyp.add((*jIt).second, "jet");
  }

  //add the met
  Handle<View<Candidate> > hMET; 
  iEvent.getByLabel(objConfig["MET"].getParameter<edm::InputTag>("source"), hMET);
  CandidatePtr met=hMET->ptrAt(0);
  hyp.add(met, "met");

  //if event is MC filter out the genparticle collection also
  if(!iEvent.isRealData()) {
    Handle<View<Candidate> > hGen;
    iEvent.getByLabel(objConfig["Generator"].getParameter<edm::InputTag>("source"), hGen);
    std::map<std::string,std::vector<CandidatePtr> > genEvent=gen::filter(hGen, objConfig["Generator"]);
    for(std::map<std::string,std::vector<CandidatePtr> >::iterator it=genEvent.begin();
	it!=genEvent.end(); ++it) {
      for(std::vector<CandidatePtr>::iterator itt=it->second.begin(); 
	  itt!=it->second.end(); ++itt)
	hyp.add( *itt, it->first );
    }
  }
      
  // work done, save results
  auto_ptr<std::vector<pat::EventHypothesis> > hyps(new std::vector<pat::EventHypothesis>() );
  hyps->push_back(hyp);
  iEvent.put(hyps,"selectedEvent");

  auto_ptr<std::vector<int> > selectionInfo(new std::vector<int>());
  selectionInfo->push_back( selPath );
  selectionInfo->push_back( selStep );
  iEvent.put(selectionInfo,"selectionInfo");

  auto_ptr<reco::VertexCollection> selVertex(new reco::VertexCollection() );
  if(theSelVertex) selVertex->push_back( *theSelVertex );   
  iEvent.put(selVertex,"selectedVertices");

  auto_ptr<reco::MuonCollection> selAllTracks(new reco::MuonCollection() );
  for(reco::MuonCollection::iterator tIt=selVtxTracks.begin(); tIt!=selVtxTracks.end(); ++tIt) {
    selAllTracks->push_back( *tIt ); 
  }
  iEvent.put(selAllTracks,"selectedTracks");

}

//  
DEFINE_FWK_MODULE(DileptonPlusMETEventProducer);
