//
#include <sstream>
#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
//#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "DataFormats/GEMRecHit/interface/ME0SegmentCollection.h"
#include "DataFormats/MuonDetId/interface/ME0DetId.h"

//#include "SimMuon/MCTruth/plugins/MuonTrackProducer.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "Geometry/GEMGeometry/interface/ME0Geometry.h"
#include "Geometry/GEMGeometry/interface/ME0EtaPartition.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"

#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixStateInfo.h"

#include "DataFormats/Math/interface/deltaPhi.h"

//
// Class declaration
//

using namespace edm;
using std::vector;
using std::map;
using std::string;

class TPMuonTrackProducer : public edm::stream::EDProducer<> {
public:
  explicit TPMuonTrackProducer(const edm::ParameterSet&);
  virtual ~TPMuonTrackProducer();

  typedef std::map<unsigned int, std::vector<unsigned int> > MuonSegmentMap; 

private:
  virtual void produce(edm::Event&, const edm::EventSetup&);
  
  edm::Handle<TrackingParticleCollection> trackParticles;
  edm::EDGetTokenT<TrackingParticleCollection> tpcToken; 
};


// 
// Class implementation 
// 

TPMuonTrackProducer::TPMuonTrackProducer(const edm::ParameterSet& parset) :
  tpcToken(consumes<TrackingParticleCollection>(parset.getParameter<edm::InputTag>("tpCollectionTag")))
{
  produces<TrackingParticleCollection>();
}

TPMuonTrackProducer::~TPMuonTrackProducer() {}

void TPMuonTrackProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) { 

  std::unique_ptr<TrackingParticleCollection> selectedTrackingParticles(new TrackingParticleCollection);

  TrackingParticleRefProd rTrackParts = iEvent.getRefBeforePut<TrackingParticleCollection>();

  // TrackingParticles 
  iEvent.getByToken(tpcToken, trackParticles); 
  if(!trackParticles.isValid()) {
    std::cout << "TrackingParticle collection not valid!" << std::endl;
    return;
  } 
  else {
    //if(debugPrint) std::cout << "Found TrackingParticle collection!" << std::endl;
  }
  for(size_t h=0; h<trackParticles->size(); ++h) { 
    if((*trackParticles)[h].status()!=1) continue; 
    if(abs((*trackParticles)[h].pdgId())!=13) continue; 
    //selectedTrackingParticles->push_back(TrackingParticleRef(trackParticles, h)); 
    selectedTrackingParticles->push_back((*trackParticles)[h]); 
  } 
  
  iEvent.put(std::move(selectedTrackingParticles));
}

// define this as a plug-in
DEFINE_FWK_MODULE(TPMuonTrackProducer);
