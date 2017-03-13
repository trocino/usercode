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
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

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
#include "DataFormats/MuonReco/interface/ME0Muon.h" 
#include "DataFormats/MuonReco/interface/ME0MuonCollection.h" 
#include "DataFormats/GEMRecHit/interface/ME0SegmentCollection.h"
#include "DataFormats/MuonDetId/interface/ME0DetId.h"

//#include "SimMuon/MCTruth/plugins/MuonTrackProducer.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "Geometry/GEMGeometry/interface/ME0Geometry.h"
#include "Geometry/GEMGeometry/interface/ME0EtaPartition.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"

#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
//#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
//#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixStateInfo.h"

#include "DataFormats/Math/interface/deltaPhi.h"

//
// Class declaration
//

using namespace edm;
using std::vector;
using std::map;
using std::string;

class ME0MuonTrackProducer : public edm::stream::EDProducer<> {
public:
  explicit ME0MuonTrackProducer(const edm::ParameterSet&);
  virtual ~ME0MuonTrackProducer();

  //typedef std::map<unsigned int, std::vector<unsigned int> > MuonSegmentMap; 

private:
  virtual void produce(edm::Event&, const edm::EventSetup&);
  
  bool isGoodME0Muon(edm::ESHandle<ME0Geometry>, const reco::ME0Muon&, double, double, double, double, double); 
  //bool isGoodME0Muon(edm::ESHandle<Propagator>, edm::ESHandle<MagneticField>, edm::ESHandle<ME0Geometry>, const reco::ME0Muon&, double, double, double, double, double); 

  edm::Handle<ME0MuonCollection> muonCollectionH;
  edm::EDGetTokenT<ME0MuonCollection> muonsToken;

  edm::Handle<ME0SegmentCollection> me0Segments; 
  edm::EDGetTokenT<ME0SegmentCollection> me0SegmentsToken;

  std::string trackType;

  // Track-ME0Segment selection cuts 
  double maxPullX; 
  double maxDiffX; 
  double maxPullY; 
  double maxDiffY; 
  double maxDiffPhiDir; 

  double minP; 
  double minPt; 

  // Temporary plots 
  edm::Service<TFileService> outfile_; 
  std::map<std::string, TH1*> hists_; 
};


// 
// Class implementation 
// 

ME0MuonTrackProducer::ME0MuonTrackProducer(const edm::ParameterSet& parset) :
  muonsToken(consumes<ME0MuonCollection>(parset.getParameter<edm::InputTag>("muonsTag"))),
  me0SegmentsToken(consumes<ME0SegmentCollection>(parset.getParameter<edm::InputTag>("me0segmentsTag"))), 
  trackType(parset.getParameter<std::string>("trackType")) 
{
  produces<reco::TrackCollection>();
  produces<reco::TrackExtraCollection>();
  produces<TrackingRecHitCollection>();

  maxPullX      = parset.getUntrackedParameter<double>("maxPullX"     ,  3.0  ); 
  maxDiffX      = parset.getUntrackedParameter<double>("maxDiffX"     ,  4.0  ); 
  maxPullY      = parset.getUntrackedParameter<double>("maxPullY"     , 20.0  ); 
  maxDiffY      = parset.getUntrackedParameter<double>("maxDiffY"     , 20.0  ); 
  maxDiffPhiDir = parset.getUntrackedParameter<double>("maxDiffPhiDir",  3.14 ); 
  minP          = parset.getUntrackedParameter<double>("minP"         , -1.0  ); 
  minPt         = parset.getUntrackedParameter<double>("minPt"        , -1.0  ); 

  // Temporary plots 
  hists_["segm_pos_x"] = outfile_->make<TH1F>("segm_pos_x", "segm_pos_x",  60, -30., 30.); 
  hists_["trck_pos_x"] = outfile_->make<TH1F>("trck_pos_x", "trck_pos_x",  60, -30., 30.); 
  hists_["diff_pos_x"] = outfile_->make<TH1F>("diff_pos_x", "diff_pos_x",  60, -30., 30.); 

  hists_["segm_pos_y"] = outfile_->make<TH1F>("segm_pos_y", "segm_pos_y", 120, -60., 60.); 
  hists_["trck_pos_y"] = outfile_->make<TH1F>("trck_pos_y", "trck_pos_y", 120, -60., 60.); 
  hists_["diff_pos_y"] = outfile_->make<TH1F>("diff_pos_y", "diff_pos_y", 120, -60., 60.); 

  hists_["segm_locdir_phi"] = outfile_->make<TH1F>("segm_locdir_phi", "segm_locdir_phi", 640, -3.2, 3.2); 
  hists_["trck_locdir_phi"] = outfile_->make<TH1F>("trck_locdir_phi", "trck_locdir_phi", 640, -3.2, 3.2); 
  hists_["diff_locdir_phi"] = outfile_->make<TH1F>("diff_locdir_phi", "diff_locdir_phi", 640, -3.2, 3.2); 

  hists_["segm_glbdir_phi"] = outfile_->make<TH1F>("segm_glbdir_phi", "segm_glbdir_phi", 640, -3.2, 3.2); 
  hists_["trck_glbdir_phi"] = outfile_->make<TH1F>("trck_glbdir_phi", "trck_glbdir_phi", 640, -3.2, 3.2); 
  hists_["diff_glbdir_phi"] = outfile_->make<TH1F>("diff_glbdir_phi", "diff_glbdir_phi", 640, -3.2, 3.2); 

}

ME0MuonTrackProducer::~ME0MuonTrackProducer() {}

void ME0MuonTrackProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) { 

  std::unique_ptr<reco::TrackCollection> selectedTracks(new reco::TrackCollection);
  std::unique_ptr<reco::TrackExtraCollection> selectedTrackExtras(new reco::TrackExtraCollection());
  std::unique_ptr<TrackingRecHitCollection> selectedTrackHits(new TrackingRecHitCollection());

  reco::TrackRefProd rTracks = iEvent.getRefBeforePut<reco::TrackCollection>();
  reco::TrackExtraRefProd rTrackExtras = iEvent.getRefBeforePut<reco::TrackExtraCollection>();
  TrackingRecHitRefProd rHits = iEvent.getRefBeforePut<TrackingRecHitCollection>();

  // ME0 muons 
  iEvent.getByToken(muonsToken, muonCollectionH);
  if(!muonCollectionH.isValid()) {
    std::cout << "Muon collection not valid!" << std::endl;
    return;
  } 
  else {
    //if(debugPrint) std::cout << "Found Muon collection!" << std::endl;
  }

  // ME0 segments 
  iEvent.getByToken(me0SegmentsToken, me0Segments); 
  if(!me0Segments.isValid()) {
    std::cout << "ME0Segment collection not valid!" << std::endl;
    return;
  } 
  else {
    //if(debugPrint) std::cout << "Found ME0Segment collection!" << std::endl;
  }

  edm::ESHandle<ME0Geometry> hGeom; 
  iSetup.get<MuonGeometryRecord>().get(hGeom); 


  // // SteppingHelix propagator 
  // edm::ESHandle<Propagator> propag;
  // iSetup.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorAlong", propag);

  // // Magnetic field 
  // edm::ESHandle<MagneticField> bField; 
  // iSetup.get<IdealMagneticFieldRecord>().get(bField); 

  edm::Ref<reco::TrackExtraCollection>::key_type idx = 0;
  edm::Ref<reco::TrackExtraCollection>::key_type hidx = 0;

  edm::LogVerbatim("ME0MuonTrackProducer") <<"\nThere are "<< muonCollectionH->size() <<" reco::Muons.";
  unsigned int muon_index = 0;
  for(ME0MuonCollection::const_iterator muon=muonCollectionH->begin();
      muon!=muonCollectionH->end(); ++muon, muon_index++) {
    edm::LogVerbatim("ME0MuonTrackProducer") <<"\n******* muon index : "<<muon_index;

    // Select ME0 muons 
    bool isGoodResult = (muon->innerTrack()->p()>minP || muon->innerTrack()->pt()>minPt); 

    if(isGoodResult) 
      //isGoodResult = isGoodME0Muon(propag, bField, hGeom, *muon, maxPullX, maxDiffX, maxPullY, maxDiffY, maxDiffPhiDir); 
      isGoodResult = isGoodME0Muon(hGeom, *muon, maxPullX, maxDiffX, maxPullY, maxDiffY, maxDiffPhiDir); 

    if(isGoodResult) {
      // new copy of Track
      reco::TrackRef trackref = muon->innerTrack();
      //const reco::Track *trk = &(*trackref); 
      const reco::Track *trk = trackref.get(); 
      reco::Track *newTrk = new reco::Track(*trk); 
      newTrk->setExtra( reco::TrackExtraRef(rTrackExtras, idx++) ); 
      PropagationDirection seedDir = trk->seedDirection(); 
      reco::TrackExtra *newExtra = new reco::TrackExtra( trk->outerPosition(), trk->outerMomentum(), 
							 trk->outerOk(), trk->innerPosition(), 
							 trk->innerMomentum(), trk->innerOk(),
							 trk->outerStateCovariance(), trk->outerDetId(),
							 trk->innerStateCovariance(), trk->innerDetId(), 
							 seedDir ); 

      // New copy of the silicon hits 
      // Add hit refs to Extra and hits to hit collection
      unsigned int nHitsToAdd = 0;
      for(trackingRecHit_iterator iHit=trk->recHitsBegin(); iHit!=trk->recHitsEnd(); ++iHit) {
        TrackingRecHit* hit = (*iHit)->clone(); 
        selectedTrackHits->push_back(hit); 
        ++nHitsToAdd; 
      }
      newExtra->setHits(rHits, hidx, nHitsToAdd);
      hidx += nHitsToAdd;
      
      selectedTracks->push_back(*newTrk);
      selectedTrackExtras->push_back(*newExtra);

    } // if (isGoodResult)
  }  // loop on reco::MuonCollection
  
  iEvent.put(std::move(selectedTracks));
  iEvent.put(std::move(selectedTrackExtras));
  iEvent.put(std::move(selectedTrackHits));
}



bool ME0MuonTrackProducer::isGoodME0Muon(edm::ESHandle<ME0Geometry> me0Geom, const reco::ME0Muon& me0muon, double MaxPullX, double MaxDiffX, double MaxPullY, double MaxDiffY, double MaxDiffPhiDir) 
{ 
  const ME0Segment &thisSegment = me0muon.me0segment(); 
  auto thisChamber = me0Geom->chamber(thisSegment.me0DetId()); 

  // Local positions (at ME0Chamber surface) for x and y cuts 
  const LocalPoint &segmLocPosition = thisSegment.localPosition(); 
  const LocalPoint &muonLocPosition = me0muon.localTrackPosAtSurface(); 

  // Global momentum (at ME0Chamber surface) for phi cut 
  const GlobalVector &muonGlbMomentum = me0muon.globalTrackMomAtSurface(); 

  // Local errors for x and u cuts 
  Double_t sigmax = sqrt( (me0muon.localTrackCov())[3][3] + thisSegment.localPositionError().xx() ); 
  Double_t sigmay = sqrt( (me0muon.localTrackCov())[4][4] + thisSegment.localPositionError().yy() ); 

  // Temporary plots 
  hists_["segm_pos_x"]->Fill(segmLocPosition.x()); 
  hists_["trck_pos_x"]->Fill(muonLocPosition.x()); 
  hists_["diff_pos_x"]->Fill(segmLocPosition.x() - muonLocPosition.x()); 

  hists_["segm_pos_y"]->Fill(segmLocPosition.y()); 
  hists_["trck_pos_y"]->Fill(muonLocPosition.y()); 
  hists_["diff_pos_y"]->Fill(segmLocPosition.y() - muonLocPosition.y()); 

  hists_["segm_locdir_phi"]->Fill(thisSegment.localDirection().phi()); 
  hists_["trck_locdir_phi"]->Fill(me0muon.localTrackMomAtSurface().phi()); 
  hists_["diff_locdir_phi"]->Fill(me0muon.localTrackMomAtSurface().phi() - thisSegment.localDirection().phi()); 

  hists_["segm_glbdir_phi"]->Fill(thisChamber->toGlobal(thisSegment.localDirection()).barePhi()); 
  hists_["trck_glbdir_phi"]->Fill(muonGlbMomentum.barePhi()); 
  hists_["diff_glbdir_phi"]->Fill(thisChamber->toGlobal(thisSegment.localDirection()).barePhi() - muonGlbMomentum.barePhi()); 

  bool xmatch(false), ymatch(false), phimatch(false); 

  if(( (std::abs(segmLocPosition.x() - muonLocPosition.x())/sigmax) < MaxPullX ) || 
     (  std::abs(segmLocPosition.x() - muonLocPosition.x())         < MaxDiffX )) 
    xmatch = true;

  if(MaxPullY<0. && MaxDiffY<0.) ymatch = true; 
  else if(( (std::abs(segmLocPosition.y() - muonLocPosition.y())/sigmay) < MaxPullY ) || 
	  (  std::abs(segmLocPosition.y() - muonLocPosition.y())         < MaxDiffY )) 
    ymatch = true;

  if(MaxDiffPhiDir<0.) phimatch = true; 
  else if(std::abs( reco::deltaPhi(muonGlbMomentum.barePhi(), thisChamber->toGlobal(thisSegment.localDirection()).barePhi()) ) < MaxDiffPhiDir) 
    phimatch = true; 

  return (xmatch && ymatch && phimatch); 
} 


/*
bool ME0MuonTrackProducer::isGoodME0Muon(edm::ESHandle<Propagator> prop, edm::ESHandle<MagneticField> bfield, edm::ESHandle<ME0Geometry> me0Geom, const reco::ME0Muon& me0muon, double MaxPullX, double MaxDiffX, double MaxPullY, double MaxDiffY, double MaxDiffPhiDir )
{
  ME0Segment thisSegment = me0muon.me0segment();
  LocalPoint thisPosition(thisSegment.localPosition());
  LocalVector thisDirection(thisSegment.localDirection().x(),thisSegment.localDirection().y(),thisSegment.localDirection().z());  //FIXME

  ME0DetId id = thisSegment.me0DetId(); 
  //auto roll = me0Geom->etaPartition(id); // Geometry/GEMGeometry/interface/ME0EtaPartition.h 
  auto roll = me0Geom->chamber(id); // Geometry/GEMGeometry/interface/ME0EtaPartition.h 
  reco::TrackRef thisTrk = me0muon.innerTrack(); 
  float zSign = thisTrk->outerZ()/fabs(thisTrk->outerZ());
  //if( zSign*me0muon.globalTrackPosAtSurface().z() < 0 ) return false;
  if( zSign*roll->toGlobal(thisPosition).z() < 0 ) return false;

  GlobalPoint r3reco(thisTrk->outerX(), thisTrk->outerY(), thisTrk->outerZ()); 
  GlobalVector p3reco(thisTrk->outerPx(), thisTrk->outerPy(), thisTrk->outerPz()); 
  GlobalTrajectoryParameters gtp(r3reco, p3reco, me0muon.charge(), bfield.product()); 
  CurvilinearTrajectoryError cov(thisTrk->outerStateCovariance()); 
  FreeTrajectoryState fts(gtp, cov); 
  TrajectoryStateOnSurface me0state = prop->propagate(fts, roll->surface()); 
  if(!me0state.isValid()) { 
    std::cout << "Propagation to ME0 eta partition failed!" << std::endl; 
    return false; 
  } 

  LocalPoint r3FinalReco = me0state.localPosition(); 
  LocalVector p3FinalReco = me0state.localMomentum(); 

  // AlgebraicMatrix thisCov(4,4,0);   
  // for(int i=1; i<=4; i++){
  //   for(int j=1; j<=4; j++){
  //     thisCov(i,j) = thisSegment.parametersError()(i,j); 
  //   }
  // }

  /////////////////////////////////////////////////////////////////////////////////////////


  // LocalTrajectoryParameters ltp(r3FinalReco,p3FinalReco,me0muon.trackCharge());
  // JacobianCartesianToLocal jctl(roll->surface(),ltp);
  // AlgebraicMatrix56 jacobGlbToLoc = jctl.jacobian(); 
  // AlgebraicMatrix55 Ctmp =  (jacobGlbToLoc * me0muon.globalTrackCov()) * ROOT::Math::Transpose(jacobGlbToLoc); 
  // AlgebraicMatrix55 Ctmp =  me0muon.localTrackCov(); 

  // AlgebraicSymMatrix55 C;  // I couldn't find any other way, so I resort to the brute force
  // for(int i=0; i<5; ++i) {
  //   for(int j=0; j<5; ++j) {
  //     C[i][j] = Ctmp[i][j]; 
  //   }
  // }  

  Double_t sigmax = sqrt( me0state.localError().positionError().xx() + thisSegment.localPositionError().xx() );      
  Double_t sigmay = sqrt( me0state.localError().positionError().yy() + thisSegment.localPositionError().yy() );

  bool X_MatchFound(false), Y_MatchFound(false), Dir_MatchFound(false);
  
  // Temporary plots 
  hists_["segm_pos_x"]->Fill(thisPosition.x()); 
  hists_["trck_pos_x"]->Fill(r3FinalReco.x()); 
  hists_["diff_pos_x"]->Fill(thisPosition.x()-r3FinalReco.x()); 

  hists_["segm_pos_y"]->Fill(thisPosition.y()); 
  hists_["trck_pos_y"]->Fill(r3FinalReco.y()); 
  hists_["diff_pos_y"]->Fill(thisPosition.y()-r3FinalReco.y()); 

  hists_["segm_locdir_phi"]->Fill(thisSegment.localDirection().phi()); 
  hists_["trck_locdir_phi"]->Fill(p3FinalReco.phi()); 
  hists_["diff_locdir_phi"]->Fill(p3FinalReco.phi()-thisSegment.localDirection().phi()); 

  hists_["segm_glbdir_phi"]->Fill(roll->toGlobal(thisSegment.localDirection()).phi()); 
  hists_["trck_glbdir_phi"]->Fill(roll->toGlobal(p3FinalReco).phi()); 
  hists_["diff_glbdir_phi"]->Fill(roll->toGlobal(p3FinalReco).phi()-roll->toGlobal(thisSegment.localDirection()).phi()); 


  if( ( (fabs(thisPosition.x()-r3FinalReco.x())/sigmax ) < MaxPullX ) || (fabs(thisPosition.x()-r3FinalReco.x()) < MaxDiffX ) ) X_MatchFound = true;

  if( MaxPullY < 0. && MaxDiffY < 0. ) Y_MatchFound = true;
  else if( ( (fabs(thisPosition.y()-r3FinalReco.y())/sigmay ) < MaxPullY ) || (fabs(thisPosition.y()-r3FinalReco.y()) < MaxDiffY ) ) Y_MatchFound = true;
  
  if( MaxDiffPhiDir < 0. ) Dir_MatchFound = true;
  else if( fabs( roll->toGlobal(p3FinalReco).phi() - roll->toGlobal(thisSegment.localDirection()).phi()) < MaxDiffPhiDir) Dir_MatchFound = true; // fix with deltaPhi

  return (X_MatchFound && Y_MatchFound && Dir_MatchFound);
}
*/ 

// define this as a plug-in
DEFINE_FWK_MODULE(ME0MuonTrackProducer);
