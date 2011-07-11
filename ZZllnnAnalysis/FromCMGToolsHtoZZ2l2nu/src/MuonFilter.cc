#include "CMGTools/HtoZZ2l2nu/interface/ObjectFilters.h"

using namespace std;

namespace muon{
  
  //
  CandidateWithVertexCollection filter(edm::Handle<edm::View<reco::Candidate> > &hMu, 
				       std::vector<reco::VertexRef> &selVtx, 
				       const edm::ParameterSet &iConfig,
				       const edm::EventSetup &iSetup)
  {
    CandidateWithVertexCollection selMuons;
    
    try {
    
      //config parameters
      double minPt = iConfig.getParameter<double>("minPt");
      double maxEta = iConfig.getParameter<double>("maxEta");
      std::string id = iConfig.getParameter<string>("id");
      double maxRelIso = iConfig.getParameter<double>("maxRelIso");
      double maxTrackChi2 = iConfig.getParameter<double>("maxTrackChi2");
      int minValidPixelHits = iConfig.getParameter<int>("minValidPixelHits");
      int minValidTrackerHits = iConfig.getParameter<int>("minValidTrackerHits");
      int minValidMuonHits = iConfig.getParameter<int>("minValidMuonHits");
      int minMatches = iConfig.getParameter<int>("minMatches");
      double maxDxy = iConfig.getParameter<double>("maxDxy");
      double maxDz  = iConfig.getParameter<double>("maxDz");

      //iterate over the muons
      for(size_t iMuon=0; iMuon< hMu.product()->size(); ++iMuon)      
	{
	  reco::CandidatePtr muonPtr = hMu->ptrAt(iMuon);
	  const pat::Muon *muon=dynamic_cast<const pat::Muon *>( muonPtr.get() );

	  //muon type
	  bool isTracker = muon->isTrackerMuon();
	  bool isGlobal = muon->isGlobalMuon();
	  if(!isTracker || !isGlobal) continue;

	  //kinematics
	  double mPt = muon->pt();
	  double mEta = muon->eta();
	  if( mPt<minPt || fabs(mEta)>maxEta) continue; 

	  //track selection
	  reco::TrackRef glbTrk=muon->globalTrack();
	  double chi2 = glbTrk->normalizedChi2();	  
	  int nValidPixelHits = glbTrk->hitPattern().numberOfValidPixelHits();
	  int nValidTrackerHits = glbTrk->hitPattern().numberOfValidTrackerHits();
	  int nValidMuonHits = glbTrk->hitPattern().numberOfValidMuonHits();
	  int nMatches = muon->numberOfMatches();
	  if( chi2>maxTrackChi2                     || 
	      nValidPixelHits<minValidPixelHits     || 
	      nValidTrackerHits<minValidTrackerHits || 
	      nValidMuonHits<minValidMuonHits       || 
	      nMatches<minMatches                    ) 
	    continue;

	  //id 
	  if(!id.empty()) {
	    bool hasId = (muon->muonID(id)>0) ;
	    if(!hasId) continue;
	  }
	  
	  //isolation
	  double norm = std::max((double)20.0, (double)mPt);
	  double ecalIso = muon->ecalIso();
	  double hcalIso = muon->hcalIso();
	  double trkIso = muon->trackIso();
	  double totalIso = (ecalIso+hcalIso+trkIso);
	  double relIso = totalIso/norm;
	  if(relIso>maxRelIso) continue;

	  //vertex selection and IP cuts
	  reco::TrackRef innTrk=muon->innerTrack();
	  reco::VertexRef vtx=vertex::getClosestVertexTo<reco::Track>(innTrk.get(), selVtx, iSetup, true);
	  if( fabs(innTrk->dxy(vtx->position()))>maxDxy ||
	      fabs(innTrk->dz(vtx->position()))>maxDz    ) continue;

	  //muon is selected
	  selMuons.push_back( CandidateWithVertex(vtx, muonPtr) );
	}
      
    } catch(std::exception &e) {
      std::cout << "[muon::filter] failed with " << e.what() << endl;
    } 
    
    return selMuons;
  }
}
