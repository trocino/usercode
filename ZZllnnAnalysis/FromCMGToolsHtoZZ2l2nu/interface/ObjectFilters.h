#ifndef objectfilters_h
#define objectfilters_h

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "TVector3.h"
#include "TH1D.h"

#include "Math/LorentzVector.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;
typedef std::vector<LorentzVector> LorentzVectorCollection;
typedef std::pair<reco::VertexRef, reco::CandidatePtr> CandidateWithVertex;
typedef std::vector<CandidateWithVertex> CandidateWithVertexCollection;
typedef std::pair<reco::VertexRef, reco::Muon> TrackWithVertex;
typedef std::vector<TrackWithVertex> TrackWithVertexCollection;
typedef std::pair<reco::VertexRef, std::vector<reco::CandidatePtr> > DileptonWithVertex;

namespace gen
{
  std::map<std::string, std::vector<reco::CandidatePtr> > filter(edm::Handle<edm::View<reco::Candidate> > &hGen, const edm::ParameterSet &iConfig);
  const reco::Candidate *getFinalStateFor(const reco::Candidate *p);
}

namespace vertex
{
  std::vector<reco::VertexRef> filter(edm::Handle<reco::VertexCollection> &hVtx, const edm::ParameterSet &iConfig);
  float getVertexMomentumFlux(const reco::Vertex *vtx, float minWeight=0.5);

  template<class T>
  std::pair<bool,Measurement1D> getImpactParameter(const T &trk, reco::Vertex *vtx, const edm::EventSetup &iSetup, bool is3d=true)
    {
      edm::ESHandle<TransientTrackBuilder> trackBuilder;
      iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", trackBuilder);
      reco::TransientTrack tt = trackBuilder->build(trk);
      if(is3d) return IPTools::absoluteImpactParameter3D(tt, *vtx);
      else     return IPTools::absoluteTransverseImpactParameter(tt, *vtx);
    }
  
  template<class T>
  reco::VertexRef getClosestVertexTo(const T *trk, std::vector<reco::VertexRef> &selVertices, const edm::EventSetup &iSetup, bool is3d=true)
    {
      reco::VertexRef bestvtx;
      double bestDz(1.0e7);
      for(std::vector<reco::VertexRef>::iterator vIt = selVertices.begin(); vIt != selVertices.end(); ++vIt)
	{
	  double dz = fabs( trk->dz( vIt->get()->position() ) );
	  if( dz>bestDz ) continue;
	  bestvtx=*vIt;
	  bestDz=dz;
	}
      return bestvtx;
    }
}

namespace muon
{
  CandidateWithVertexCollection filter(edm::Handle<edm::View<reco::Candidate> > &hMu, 
				       std::vector<reco::VertexRef> &selVtx, 
				       const edm::ParameterSet &iConfig,
				       const edm::EventSetup &iSetup);
}

namespace electron
{
  CandidateWithVertexCollection filter(edm::Handle<edm::View<reco::Candidate> > &hEle, 
				       edm::Handle<edm::View<reco::Candidate> > &hMu, 
				       std::vector<reco::VertexRef> &selVtx, 
				       const edm::ParameterSet &iConfig,
				       const edm::EventSetup &iSetup);
}

namespace dilepton
{
  enum DileptonClassification {UNKNOWN=0,MUMU=1,EE=2,EMU=3};
  DileptonWithVertex filter(CandidateWithVertexCollection &selLeptons, 
			    std::vector<reco::VertexRef> &selVertices, 
			    const edm::ParameterSet &iConfig,
			    const edm::EventSetup &iSetup,
			    double rho,
			    CandidateWithVertexCollection &isolLeptons,
			    std::map<TString, TH1D *> *controlHistos_=0);

  int classify(std::vector<reco::CandidatePtr> &selDilepton);
  int getLeptonId(reco::CandidatePtr &lepton);
  double getPtErrorFor(reco::CandidatePtr &lepton);
  enum IsolType { ECAL_ISO=0, HCAL_ISO, TRACKER_ISO};
  std::vector<double> getLeptonIso(reco::CandidatePtr &lepton);
  const reco::GenParticle *getLeptonGenMatch(reco::CandidatePtr &lepton);
}

namespace track
{
  TrackWithVertexCollection filter(edm::Handle<std::vector<reco::Track> > &hTrks, 
				   CandidateWithVertexCollection &isoLepts, 
				   std::vector<reco::VertexRef> &selVtx, 
				   const edm::ParameterSet &iConfig, 
				   const edm::EventSetup &iSetup);
  std::pair<int, double> computeTrackIsolation(const reco::Track *trk, 
					       edm::Handle<std::vector<reco::Track> > &hTracks,
					       reco::VertexRef &vx, bool vxConstr=true, 
					       double minPtForIso=0.5,
					       double maxDzCut=0.1);
}

namespace jet
{
  CandidateWithVertexCollection filter(edm::Handle<edm::View<reco::Candidate> > &hJet, 
				       CandidateWithVertexCollection &selLeptons, 
				       std::vector<reco::VertexRef> &selVtx, 
				       const edm::ParameterSet &iConfig);
  double fAssoc(const pat::Jet *jet, const reco::Vertex *vtx);
}

namespace met
{
  enum MetTypes { RAWMET=0, TYPEIMET, CORRECTED_TYPEIMET };
  LorentzVectorCollection filter(const reco::PFMET &pfMET,  std::vector<const pat::Jet *> &assocJets,  std::vector<const pat::Jet *> &puJets);
}



#endif
