#include "CMGTools/HtoZZ2l2nu/interface/ObjectFilters.h"

using namespace std;

namespace jet{

  //
  CandidateWithVertexCollection filter(edm::Handle<edm::View<reco::Candidate> > &hJet, 
				       CandidateWithVertexCollection &selLeptons, 
				       std::vector<reco::VertexRef> &selVtx, 
				       const edm::ParameterSet &iConfig)
  {
    CandidateWithVertexCollection selJets;

    using namespace edm;
    try {

      // Config parameters
      double minPt = iConfig.getParameter<double>("minPt");
      double maxEta = iConfig.getParameter<double>("maxEta");
      double minDeltaRtoLepton = iConfig.getParameter<double>("minDeltaRtoLeptons");
      PFJetIDSelectionFunctor jetIdSelector( iConfig.getParameter<edm::ParameterSet>("jetId") );
      pat::strbitset hasId = jetIdSelector.getBitTemplate();

      // Iterate over the jets
      for(size_t iJet=0; iJet< hJet.product()->size(); ++iJet) {
	reco::CandidatePtr jetPtr=hJet->ptrAt(iJet);
	const pat::Jet *jet=dynamic_cast<const pat::Jet *>( jetPtr.get() );

	// Basic kinematics
	double pt=jet->pt();
	double eta=jet->eta();
	if(pt<minPt || fabs(eta)>maxEta) continue;

	// Jet id
	hasId.set(false);
	if( !jetIdSelector( *jet, hasId ) ) continue;

	// Vertex selection
	reco::VertexRef sVtx;
	double assoVar(0.);
	for(std::vector<reco::VertexRef>::iterator vIt=selVtx.begin(); 
	    vIt!=selVtx.end(); ++vIt) {
	  double tmpAss=fAssoc(jet, vIt->get());
	  if(tmpAss>assoVar) {
	    assoVar=tmpAss;
	    sVtx=(*vIt);
	  }
	}
	if(sVtx.get()==0) continue; // FIXME!!!

	// Check overlaps with selected leptons
	double minDR(1000);
	for(CandidateWithVertexCollection::iterator lIt=selLeptons.begin(); 
	    lIt!=selLeptons.end(); ++lIt) {
	  double dR=deltaR( *jet, *(lIt->second.get()) );
	  if(dR>minDR) continue; 
	  minDR=dR;
	}
	if(minDR<minDeltaRtoLepton) continue;

	// Jet is selected
	selJets.push_back( CandidateWithVertex(sVtx, jetPtr) );

      }
    } catch(std::exception &e) {
      cout << "[jet::filter] failed with " << e.what() << endl;
    }
    
    return selJets;
  }


  //
  double fAssoc(const pat::Jet *jet, const reco::Vertex *vtx) {
    double fassoc(-1.);
    if(jet==0 || vtx==0) return fassoc;

    // Iterate over the tracks associated to a jet
    double sumpttracks(0), assocsumpttracks(0);
    const reco::TrackRefVector &jtracks=jet->associatedTracks();
    for(reco::TrackRefVector::const_iterator jtIt=jtracks.begin();
	jtIt!=jtracks.end(); ++jtIt) {
      if( jtIt->isNull() ) continue;
      const reco::Track *jtrack=jtIt->get();
      sumpttracks+=jtrack->pt();

      // Find track match
      for(reco::Vertex::trackRef_iterator vtIt=vtx->tracks_begin(); 
	  vtIt!=vtx->tracks_end(); ++vtIt) {
	if( vtIt->isNull() ) continue;
	const reco::Track *vtrack=vtIt->get();
	if(vtrack!=jtrack) continue;
	assocsumpttracks+=jtrack->pt();
	break;
      }
    }
    if(sumpttracks>0) fassoc=assocsumpttracks/sumpttracks;
    return fassoc;
  }

}
