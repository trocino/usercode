#include "CMGTools/HtoZZ2l2nu/interface/ObjectFilters.h"

using namespace std;

namespace track{
  
  //
  TrackWithVertexCollection filter(edm::Handle<std::vector<reco::Track> > &hTrks, 
				   CandidateWithVertexCollection &isoLepts, 
				   std::vector<reco::VertexRef> &selVtx, 
				   const edm::ParameterSet &iConfig, 
				   const edm::EventSetup &iSetup)
  {
    TrackWithVertexCollection selTracks;

    try {
    
      // Config parameters
      double minPt = iConfig.getParameter<double>("minPt");
      double maxEta = iConfig.getParameter<double>("maxEta");
      double minDRFromLept = iConfig.getParameter<double>("minDeltaRtoLeptons");
      double maxRelIso = iConfig.getParameter<double>("maxRelIso");
      double maxTrackChi2 = iConfig.getParameter<double>("maxTrackChi2");
      int minValidPixelHits = iConfig.getParameter<int>("minValidPixelHits");
      int minValidTrackerHits = iConfig.getParameter<int>("minValidTrackerHits");
      double maxDxy = iConfig.getParameter<double>("maxDxy");
      double maxDz  = iConfig.getParameter<double>("maxDz");

      // Iterate over tracks
      for(size_t iTrk=0; iTrk<hTrks->size(); ++iTrk) {
	reco::TrackRef tRef(hTrks, iTrk);

	// Kinematics
	double tPt = tRef->pt();
	double tEta = tRef->eta();
	if( tPt<minPt || fabs(tEta)>maxEta) 
	  continue; 

	// Track selection
	double chi2 = tRef->normalizedChi2();
	int nValidPixelHits = tRef->hitPattern().numberOfValidPixelHits();
	int nValidTrackerHits = tRef->hitPattern().numberOfValidTrackerHits();
	if(chi2>maxTrackChi2                     || 
	   nValidPixelHits<minValidPixelHits     || 
	   nValidTrackerHits<minValidTrackerHits  ) 
	  continue;

	// Vertex
	reco::VertexRef vtx=vertex::getClosestVertexTo<reco::Track>(tRef.get(), selVtx, iSetup, true);
	double trkDxy=1.0e7;
	double trkDz=1.0e7;
	if(vtx.get()) {
	  trkDxy=fabs(tRef->dxy(vtx->position()));
	  trkDz=fabs(tRef->dz(vtx->position()));
	}
	if(trkDxy>maxDxy || trkDz>maxDz) continue;

	// Overlaps with selected leptons
	double minDR(1000.);
	for(CandidateWithVertexCollection::iterator lIt=isoLepts.begin(); 
	    lIt!=isoLepts.end(); ++lIt) {
	  double dR=deltaR( *tRef, *(lIt->second.get()) );
	  if(dR>minDR) continue; 
	  minDR=dR;
	}
	if(minDR<minDRFromLept) continue;

	// Create fake muon for track
	double fakeEnergy=sqrt(tRef->p()*tRef->p() + 0.011163691);  // using mass of a muon...
	math::XYZTLorentzVector fakeP4(tRef->px(), tRef->py(), tRef->pz(), fakeEnergy);
	//reco::Muon fakeMu(tRef->charge(), fakeP4, vtx->position());
	//fakeMu.setInnerTrack(tRef);

	// Isolation
	pair<int, double> trackIso=computeTrackIsolation(tRef.get(), hTrks, vtx, true, 0.5, maxDz);
	double relIso=trackIso.second/tRef->pt();
	if(relIso>maxRelIso) continue;

	reco::MuonIsolation isoR03, isoR05;
	isoR03.nTracks=trackIso.first;
	isoR03.sumPt=trackIso.second;
	//fakeMu.setIsolation(isoR03, isoR05);

	// Fill the collection
	selTracks.push_back( TrackWithVertex( vtx, reco::Muon(tRef->charge(), fakeP4, vtx->position()) ) );
	reco::Muon & refMu=selTracks.back().second; 
	refMu.setInnerTrack(tRef);
	refMu.setIsolation(isoR03, isoR05);

      } // end for(size_t iTrk=0; iTrk<hTrks->size(); ++iTrk) 

    } catch(std::exception &e) {
      std::cout << "[track::filter] failed with : " << e.what() << endl;
    }

    return selTracks;
  }


  pair<int, double> computeTrackIsolation(const reco::Track *trk, 
					  edm::Handle<std::vector<reco::Track> > &hTracks,
					  reco::VertexRef &vx, bool vxConstr, 
					  double minPtForIso,
					  double maxDzCut)
  {
    double iso(0.);
    int ntrks(0);
    //    for(size_t jTrk=0; jTrk<hTracks.product()->size(); ++jTrk) {
    for(reco::TrackCollection::const_iterator jTrk=hTracks->begin(); jTrk!=hTracks->end(); ++jTrk) {
      const reco::Track *aTrk=&*jTrk;
      if( aTrk->pt()<minPtForIso ) continue;
      if( (&*aTrk)==(&*trk) ) continue;
      if( deltaR( *trk, *aTrk )>0.3 ) continue;
      if(vx.get()!=0 && vxConstr) {
	double aDz=fabs(aTrk->dz(vx->position()));
	if(aDz>maxDzCut) continue;
      }
      iso+=aTrk->pt();
      ntrks++;
    }
    //iso/=trk->pt();
    return pair<int, double>(ntrks, iso);
  }

}
