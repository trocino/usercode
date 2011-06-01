// -*- C++ -*-
//
// Package:    ZZllnnFilter
// Class:      ZZllnnFilter
//
/**\class ZZllnnFilter ZZllnnFilter.cc AnalysisCode/LeptonJetFilter/src/ZZllnnFilter.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Jeff Temple
//         Created:  Mon Aug  3 13:02:30 CEST 2009
// $Id: ZZllnnFilter.cc,v 1.6 2010/09/11 04:27:25 ferencek Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

//#include "DataFormats/PatCandidates/interface/Muon.h"
//#include "DataFormats/PatCandidates/interface/Electron.h"

//#include "FWCore/ServiceRegistry/interface/Service.h"
//#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "TH1.h"

using namespace edm;
using namespace std;
using namespace reco;
//
// class declaration
//

class ZZllnnFilter : public edm::EDFilter {
   public:
      explicit ZZllnnFilter(const edm::ParameterSet&);
      ~ZZllnnFilter();

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------
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
ZZllnnFilter::ZZllnnFilter(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed

  // specify labels
}


ZZllnnFilter::~ZZllnnFilter()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
ZZllnnFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  bool isZZllnn=false;
  unsigned int nGenZ(0), nZtoMM(0), nZtoEE(0), nZtoNN(0);
  unsigned int nZtoMp(0), nZtoMm(0), nZtoEp(0), nZtoEm(0), nZtoN(0), nZtoAN(0);

  //   Handle<HepMCProduct> hepMcProd;
  //   iEvent.getByLabel("generator", hepMcProd);
  Handle<GenParticleCollection> genParticles;
  iEvent.getByLabel("prunedGen", genParticles);
  unsigned int gensize=genParticles->size();

  //cout << "taken hepMcProd" << endl;
  //   for( SimTrackContainer::const_iterator simTrack=simMuTracks.begin(); simTrack!=simMuTracks.end(); ++simTrack ) {
  //  for( NewMuonTrackResolutionAnalyzer::SimTracksMap::const_iterator simTrackPair=simMuTracksMap.begin(); simTrackPair!=simMuTracksMap.end(); ++simTrackPair ) {

  for(unsigned int i=0; i<gensize; ++i) {

    //cout << "taken gpIt" << endl;
    //    HepMC::GenParticle* gp = hepMcProd->GetEvent()->barcode_to_particle((*simTrackPair).second->genpartIndex());
    const reco::GenParticle& gp=(*genParticles)[i];
    //cout << "taken gp:  gp=" << gp << " - status=" << gp->status() << " - id=" << gp->pdg_id()  << endl;

    if( gp.status()==3 && gp.pdgId()==13 ) nZtoMm++;
    if( gp.status()==3 && gp.pdgId()==-13 ) nZtoMp++;
    if( gp.status()==3 && gp.pdgId()==11 ) nZtoEm++;
    if( gp.status()==3 && gp.pdgId()==-11 ) nZtoEp++;
    if( gp.status()==3 && (gp.pdgId()==12 || gp.pdgId()==14 || gp.pdgId()==16) ) nZtoN++;
    if( gp.status()==3 && (gp.pdgId()==-12 || gp.pdgId()==-14 || gp.pdgId()==-16) ) nZtoAN++;

    if( ( (nZtoMm>=1 && nZtoMp>=1) || (nZtoEm>=1 && nZtoEp>=1) ) && nZtoN>=1 && nZtoAN>=1 ) {
      isZZllnn=true;
      break;
    }
//     else 
//       break;

    /*
    if( gp!=0 && gp->status()==2 && gp->pdgId()==23 ) {
      nGenZ++;
      cout << nGenZ << endl;
      for(HepMC::GenVertex::particle_iterator daughter=gp->end_vertex()->particles_begin(HepMC::children);
	  daughter!=gp->end_vertex()->particles_end(HepMC::children); ++daughter) {
	cout << "taken daughter" << endl;
	if( abs((*daughter)->pdg_id())==13 ) {
	  nZtoMM++;
	  break;
	}
	else if( abs((*daughter)->pdg_id())==11 ) {
	  nZtoEE++;
	  break;
	}
	else if( abs((*daughter)->pdg_id())==12 || abs((*daughter)->pdg_id())==14 || abs((*daughter)->pdg_id())==16 ) {
	  nZtoNN++;
	  break;
	}
	else {
	  break;
	}
      }
    }
    if(nGenZ==2) {
      if( (nZtoMM>=1 || nZtoEE>=1) && nZtoNN>=1 ) {
	isZZllnn=true;
	break;
      }
      else 
	break;
    }
    */
  } // end for

  return isZZllnn;
}

// ------------ method called once each job just before starting event loop  ------------
void
ZZllnnFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
ZZllnnFilter::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(ZZllnnFilter);
