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
#include "DataFormats/GEMRecHit/interface/ME0SegmentCollection.h"
#include "DataFormats/MuonDetId/interface/ME0DetId.h"

//#include "SimMuon/MCTruth/plugins/MuonTrackProducer.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "Geometry/GEMGeometry/interface/ME0Geometry.h"
#include "Geometry/GEMGeometry/interface/ME0EtaPartition.h"
#include "Geometry/GEMGeometry/interface/ME0Chamber.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"

#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixStateInfo.h"

#include "DataFormats/Math/interface/deltaPhi.h"

// ROOT includes 
#include "TFile.h"
#include "TH1.h"

//
// Class declaration
//

using namespace edm;
using std::vector;
using std::map;
using std::string;

class RecoMuonTrackProducer : public edm::stream::EDProducer<> {
public:
  explicit RecoMuonTrackProducer(const edm::ParameterSet&);
  virtual ~RecoMuonTrackProducer();

  typedef std::map<unsigned int, std::vector<unsigned int> > MuonSegmentMap; 

private:
  virtual void produce(edm::Event&, const edm::EventSetup&);
  
  bool isGoodME0Muon(edm::ESHandle<Propagator>, edm::ESHandle<MagneticField>, edm::ESHandle<ME0Geometry>, const reco::Muon&, const ME0SegmentRef&, double, double, double, double, double); 

  edm::Handle<reco::MuonCollection> muonCollectionH;
  edm::EDGetTokenT<reco::MuonCollection> muonsToken;

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

RecoMuonTrackProducer::RecoMuonTrackProducer(const edm::ParameterSet& parset) :
  muonsToken(consumes<reco::MuonCollection>(parset.getParameter<edm::InputTag>("muonsTag"))),
  me0SegmentsToken( consumes<ME0SegmentCollection>(parset.getParameter<edm::InputTag>("me0segmentsTag")) ), 
  trackType(parset.getParameter<std::string>("trackType")) 
{
  produces<reco::TrackCollection>();
  produces<reco::TrackExtraCollection>();
  produces<TrackingRecHitCollection>();
  produces<MuonSegmentMap>(); 

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

RecoMuonTrackProducer::~RecoMuonTrackProducer() {}

void RecoMuonTrackProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) { 

  std::unique_ptr<reco::TrackCollection> selectedTracks(new reco::TrackCollection);
  std::unique_ptr<reco::TrackExtraCollection> selectedTrackExtras(new reco::TrackExtraCollection());
  std::unique_ptr<TrackingRecHitCollection> selectedTrackHits(new TrackingRecHitCollection());
  std::unique_ptr<MuonSegmentMap> muonSegmMap(new MuonSegmentMap());

  reco::TrackRefProd rTracks = iEvent.getRefBeforePut<reco::TrackCollection>();
  reco::TrackExtraRefProd rTrackExtras = iEvent.getRefBeforePut<reco::TrackExtraCollection>();
  TrackingRecHitRefProd rHits = iEvent.getRefBeforePut<TrackingRecHitCollection>();
  edm::RefProd<MuonSegmentMap> rMuSegmMap = iEvent.getRefBeforePut<MuonSegmentMap>();

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

  edm::ESHandle<TrackerTopology> httopo;
  iSetup.get<TrackerTopologyRcd>().get(httopo);
  const TrackerTopology& ttopo = *httopo;

  edm::ESHandle<ME0Geometry> hGeom; 
  iSetup.get<MuonGeometryRecord>().get(hGeom); 


  // SteppingHelix propagator 
  edm::ESHandle<Propagator> propag;
  iSetup.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorAlong", propag);

  // Magnetic field 
  edm::ESHandle<MagneticField> bField; 
  iSetup.get<IdealMagneticFieldRecord>().get(bField); 

  edm::Ref<reco::TrackExtraCollection>::key_type idx = 0;
  edm::Ref<reco::TrackExtraCollection>::key_type hidx = 0;

  edm::LogVerbatim("RecoMuonTrackProducer") <<"\nThere are "<< muonCollectionH->size() <<" reco::Muons.";
  unsigned int muon_index = 0;
  for(reco::MuonCollection::const_iterator muon=muonCollectionH->begin();
      muon!=muonCollectionH->end(); ++muon, muon_index++) {
    edm::LogVerbatim("RecoMuonTrackProducer") <<"\n******* muon index : "<<muon_index;

    // Select ME0 muons 
    bool isGoodResult = false; 
    if(muon->isME0Muon()) { 
      std::vector<unsigned int> segmlist; 
      //std::cout << "  ==> Muon " << muon_index << std::endl; 

      for(std::vector<reco::MuonChamberMatch>::const_iterator chamberMatch=muon->matches().begin(); 
	  chamberMatch!=muon->matches().end(); ++chamberMatch) {
	if(chamberMatch->detector()==MuonSubdetId::ME0) { 
	  //std::cout << "          Chamber match: "; 

	  for(std::vector<reco::MuonSegmentMatch>::const_iterator segmentMatch=chamberMatch->me0Matches.begin();
	      segmentMatch!=chamberMatch->me0Matches.end(); ++segmentMatch) { 
	    ME0SegmentRef me0segm = segmentMatch->me0SegmentRef; 
	    bool thissegmisgood = isGoodME0Muon(propag, bField, hGeom, *muon, me0segm, maxPullX, maxDiffX, maxPullY, maxDiffY, maxDiffPhiDir); 
	    //std::cout << " -- " << me0segm.get() << "  " << thissegmisgood; 
	    isGoodResult |= thissegmisgood; 
	    if(thissegmisgood) {  
	      unsigned int segm_index = 0;
	      for(ME0SegmentCollection::const_iterator segmit=me0Segments->begin(); segmit!=me0Segments->end(); ++segmit, segm_index++) { 
		if( me0segm.get()==(&(*segmit)) ) { 
		  segmlist.push_back(segm_index); 
		  break; 
		} 
	      }
	    }  
	  } // end for(...MuonSegmentMatch...) 
	  //std::cout << std::endl; 
	} // end if(chamberMatch->detector()==MuonSubdetId::ME0) 
      } // end for(...MuonChamberMatch...) 

      if(segmlist.size()>0) 
	(*muonSegmMap)[muon_index] = segmlist; 
    } // end if(muon->isME0Muon()) 

    if(isGoodResult) {
      // new copy of Track
      reco::TrackRef trackref;
      if (trackType == "innerTrack") {
        if (muon->innerTrack().isNonnull()) trackref = muon->innerTrack();
        else continue;
      } 
      else if (trackType == "outerTrack") {
        if (muon->outerTrack().isNonnull()) trackref = muon->outerTrack();
        else continue;
      } 
      else if (trackType == "globalTrack") {
        if (muon->globalTrack().isNonnull()) trackref = muon->globalTrack();
        else continue;
      }
      else if (trackType == "innerTrackPlusSegments") {
	if (muon->innerTrack().isNonnull()) trackref = muon->innerTrack();
	else continue;
      }

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
      if(trackType=="innerTrackPlusSegments") { 
	int wheel, station, sector;
	int endcap, /*station, */ ring, chamber;
	
	edm::LogVerbatim("RecoMuonTrackProducer") << "Number of chambers: " << muon->matches().size() 
						 << ", arbitrated: " << muon->numberOfMatches(reco::Muon::SegmentAndTrackArbitration); 
	unsigned int index_chamber = 0; 
	
	for(std::vector<reco::MuonChamberMatch>::const_iterator chamberMatch=muon->matches().begin(); 
	    chamberMatch!=muon->matches().end(); ++chamberMatch, index_chamber++) {
	  std::stringstream chamberStr;
	  chamberStr <<"\nchamber index: "<<index_chamber; 

	  int subdet = chamberMatch->detector(); 
	  DetId did = chamberMatch->id; 

	  if(subdet==MuonSubdetId::DT) { 
	    DTChamberId dtdetid = DTChamberId(did); 
	    wheel = dtdetid.wheel(); 
	    station = dtdetid.station(); 
	    sector = dtdetid.sector(); 
	    chamberStr << ", DT chamber Wh:"<<wheel<<",St:"<<station<<",Se:"<<sector; 
	  } 
	  else if(subdet == MuonSubdetId::CSC) { 
	    CSCDetId cscdetid = CSCDetId(did); 
	    endcap = cscdetid.endcap(); 
	    station = cscdetid.station(); 
	    ring = cscdetid.ring(); 
	    chamber = cscdetid.chamber(); 
	    chamberStr << ", CSC chamber End:"<<endcap<<",St:"<<station<<",Ri:"<<ring<<",Ch:"<<chamber; 
	  } 

	  chamberStr << ", Number of segments: " << chamberMatch->segmentMatches.size(); 
	  edm::LogVerbatim("RecoMuonTrackProducer") << chamberStr.str(); 

	  unsigned int index_segment = 0;
	  
	  for(std::vector<reco::MuonSegmentMatch>::const_iterator segmentMatch=chamberMatch->segmentMatches.begin();
	      segmentMatch!=chamberMatch->segmentMatches.end(); ++segmentMatch, index_segment++) {

	    float segmentX = segmentMatch->x; 
	    float segmentY = segmentMatch->y; 
	    float segmentdXdZ = segmentMatch->dXdZ; 
	    float segmentdYdZ = segmentMatch->dYdZ; 
	    float segmentXerr = segmentMatch->xErr; 
	    float segmentYerr = segmentMatch->yErr; 
	    float segmentdXdZerr = segmentMatch->dXdZErr; 
	    float segmentdYdZerr = segmentMatch->dYdZErr; 

	    CSCSegmentRef segmentCSC = segmentMatch->cscSegmentRef; 
	    DTRecSegment4DRef segmentDT = segmentMatch->dtSegmentRef; 

	    bool segment_arbitrated_Ok = (segmentMatch->isMask(reco::MuonSegmentMatch::BestInChamberByDR) && 
					  segmentMatch->isMask(reco::MuonSegmentMatch::BelongsToTrackByDR)); 

	    std::string ARBITRATED(" ***Arbitrated Off*** "); 
	    if(segment_arbitrated_Ok) ARBITRATED = " ***ARBITRATED OK*** "; 

	    if(subdet==MuonSubdetId::DT) {	      
	      edm::LogVerbatim("RecoMuonTrackProducer")
		<<"\n\t segment index: "<<index_segment << ARBITRATED
		<<"\n\t  Local Position (X,Y)=("<<segmentX<<","<<segmentY<<") +/- ("<<segmentXerr<<","<<segmentYerr<<"), " 
		<<"\n\t  Local Direction (dXdZ,dYdZ)=("<<segmentdXdZ<<","<<segmentdYdZ<<") +/- ("<<segmentdXdZerr<<","<<segmentdYdZerr<<")"; 
	      
	      if(!segment_arbitrated_Ok) continue;
	      
	      if(segmentDT.get()!=0) {
		const DTRecSegment4D* segment = segmentDT.get();
		
		edm::LogVerbatim("RecoMuonTrackProducer")<<"\t ===> MATCHING with DT segment with index = "<<segmentDT.key();
		
		if(segment->hasPhi()) {
		  const DTChamberRecSegment2D* phiSeg = segment->phiSegment();
		  std::vector<const TrackingRecHit*> phiHits = phiSeg->recHits();
                  unsigned int nHitsAdded = 0;
		  for(std::vector<const TrackingRecHit*>::const_iterator ihit=phiHits.begin();
		      ihit!=phiHits.end(); ++ihit) {
		    TrackingRecHit* seghit = (*ihit)->clone();
		    newTrk->appendHitPattern(*seghit, ttopo);
		    //		    edm::LogVerbatim("RecoMuonTrackProducer")<<"hit pattern for position "<<index_hit<<" set to:";
		    //		    newTrk->hitPattern().printHitPattern(index_hit, std::cout);
		    selectedTrackHits->push_back(seghit);
                    ++nHitsAdded;
		  }
                  newExtra->setHits(rHits, hidx, nHitsAdded);
                  hidx += nHitsAdded;
		}
		
		if(segment->hasZed()) {
		  const DTSLRecSegment2D* zSeg = (*segment).zSegment();
		  std::vector<const TrackingRecHit*> zedHits = zSeg->recHits();
                  unsigned int nHitsAdded = 0;
		  for(std::vector<const TrackingRecHit*>::const_iterator ihit=zedHits.begin();
		      ihit!=zedHits.end(); ++ihit) {
		    TrackingRecHit* seghit = (*ihit)->clone();
		    newTrk->appendHitPattern(*seghit, ttopo);
		    //		    edm::LogVerbatim("RecoMuonTrackProducer")<<"hit pattern for position "<<index_hit<<" set to:";
		    //		    newTrk->hitPattern().printHitPattern(index_hit, std::cout);
		    selectedTrackHits->push_back(seghit);
                    ++nHitsAdded;
		  }
                  newExtra->setHits(rHits, hidx, nHitsAdded);
                  hidx += nHitsAdded;
		}
	      } else edm::LogWarning("RecoMuonTrackProducer")<<"\n***WARNING: UNMATCHED DT segment ! \n";
	    } // if (subdet == MuonSubdetId::DT)

	    else if(subdet == MuonSubdetId::CSC) {
	      edm::LogVerbatim("RecoMuonTrackProducer")
		<<"\n\t segment index: "<<index_segment << ARBITRATED
		<<"\n\t  Local Position (X,Y)=("<<segmentX<<","<<segmentY<<") +/- ("<<segmentXerr<<","<<segmentYerr<<"), " 
		<<"\n\t  Local Direction (dXdZ,dYdZ)=("<<segmentdXdZ<<","<<segmentdYdZ<<") +/- ("<<segmentdXdZerr<<","<<segmentdYdZerr<<")"; 
	      
	      if(!segment_arbitrated_Ok) continue;
	      
	      if(segmentCSC.get() != 0) {
		const CSCSegment* segment = segmentCSC.get();
		
		edm::LogVerbatim("RecoMuonTrackProducer")<<"\t ===> MATCHING with CSC segment with index = "<<segmentCSC.key();
		
		std::vector<const TrackingRecHit*> hits = segment->recHits();
                unsigned int nHitsAdded = 0;
		for(std::vector<const TrackingRecHit*>::const_iterator ihit=hits.begin();
		    ihit!=hits.end(); ++ihit) {
		  TrackingRecHit* seghit = (*ihit)->clone();
		  newTrk->appendHitPattern(*seghit, ttopo);
		  //		    edm::LogVerbatim("RecoMuonTrackProducer")<<"hit pattern for position "<<index_hit<<" set to:";
		  //		    newTrk->hitPattern().printHitPattern(index_hit, std::cout);
		  selectedTrackHits->push_back(seghit);
                  ++nHitsAdded;
		}
                newExtra->setHits(rHits, hidx, nHitsAdded);
                hidx += nHitsAdded;
	      } else edm::LogWarning("RecoMuonTrackProducer")<<"\n***WARNING: UNMATCHED CSC segment ! \n";
	    }  //  else if (subdet == MuonSubdetId::CSC)

	  } // loop on vector<MuonSegmentMatch>	  
	} // loop on vector<MuonChamberMatch>	
      } // if(trackType == "innerTrackPlusSegments")
      
      //      edm::LogVerbatim("RecoMuonTrackProducer")<<"\n printing final hit_pattern";
      //      newTrk->hitPattern().print();
      
      selectedTracks->push_back(*newTrk);
      selectedTrackExtras->push_back(*newExtra);

    } // if (isGoodResult)
  }  // loop on reco::MuonCollection
  
  iEvent.put(std::move(selectedTracks));
  iEvent.put(std::move(selectedTrackExtras));
  iEvent.put(std::move(selectedTrackHits));
  iEvent.put(std::move(muonSegmMap)); 
}


bool RecoMuonTrackProducer::isGoodME0Muon(edm::ESHandle<Propagator> prop, edm::ESHandle<MagneticField> bfield, edm::ESHandle<ME0Geometry> me0Geom, const reco::Muon& me0muon, const ME0SegmentRef& thisSegment, double MaxPullX, double MaxDiffX, double MaxPullY, double MaxDiffY, double MaxDiffPhiDir )
{
  LocalPoint thisPosition(thisSegment->localPosition());
  LocalVector thisDirection(thisSegment->localDirection().x(),thisSegment->localDirection().y(),thisSegment->localDirection().z());  //FIXME

  int tmpcnt(0); 
  // std::cout << "  *** MINCHIA " << tmpcnt++ << std::endl; // 0 
  ME0DetId id = thisSegment->me0DetId(); 
  //auto roll = me0Geom->etaPartition(id); // Geometry/GEMGeometry/interface/ME0EtaPartition.h 
  // std::cout << "  *** MINCHIA " << tmpcnt++ << std::endl; // 1 
  auto roll = me0Geom->chamber(id); // Geometry/GEMGeometry/interface/ME0EtaPartition.h 
  if(roll==0) { 
    std::cout << " ~*~*~*~*~*~ No ME0 chamber found! 'sti cazzi!" << std::endl; 
    return false; 
  } 
  // std::cout << "  *** MINCHIA " << tmpcnt++ << std::endl; // 2 
  float zSign  = me0muon.pz()/fabs(me0muon.pz());
  // std::cout << "  *** MINCHIA " << tmpcnt++ << std::endl; // 3 
  if ( zSign * roll->toGlobal(thisPosition).z() < 0 ) return false;

  // std::cout << "  *** MINCHIA " << tmpcnt++ << std::endl; // 4 
  GlobalPoint r3reco(me0muon.innerTrack()->outerX(), me0muon.innerTrack()->outerY(), me0muon.innerTrack()->outerZ()); 
  GlobalVector p3reco(me0muon.innerTrack()->outerPx(), me0muon.innerTrack()->outerPy(), me0muon.innerTrack()->outerPz()); 
  GlobalTrajectoryParameters gtp(r3reco, p3reco, me0muon.charge(), bfield.product()); 
  CurvilinearTrajectoryError cov(me0muon.innerTrack()->outerStateCovariance()); 
  FreeTrajectoryState fts(gtp, cov); 
  TrajectoryStateOnSurface me0state = prop->propagate(fts, roll->surface()); 
  if(!me0state.isValid()) { 
    std::cout << "Propagation to ME0 eta partition failed!" << std::endl; 
    return false; 
  } 
  // std::cout << "  *** MINCHIA " << tmpcnt++ << std::endl; // 5 

  LocalPoint r3FinalReco = me0state.localPosition(); 
  LocalVector p3FinalReco = me0state.localMomentum(); 

  // AlgebraicMatrix thisCov(4,4,0);   
  // for(int i=1; i<=4; i++){
  //   for(int j=1; j<=4; j++){
  //     thisCov(i,j) = thisSegment->parametersError()(i,j); 
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

  Double_t sigmax = sqrt( me0state.localError().positionError().xx() + thisSegment->localPositionError().xx() );      
  Double_t sigmay = sqrt( me0state.localError().positionError().yy() + thisSegment->localPositionError().yy() );

  bool X_MatchFound(false), Y_MatchFound(false), Dir_MatchFound(false);
  
  // Temporary plots 
  hists_["segm_pos_x"]->Fill(thisPosition.x()); 
  hists_["trck_pos_x"]->Fill(r3FinalReco.x()); 
  hists_["diff_pos_x"]->Fill(thisPosition.x()-r3FinalReco.x()); 

  hists_["segm_pos_y"]->Fill(thisPosition.y()); 
  hists_["trck_pos_y"]->Fill(r3FinalReco.y()); 
  hists_["diff_pos_y"]->Fill(thisPosition.y()-r3FinalReco.y()); 

  hists_["segm_locdir_phi"]->Fill(thisSegment->localDirection().phi()); 
  hists_["trck_locdir_phi"]->Fill(p3FinalReco.phi()); 
  hists_["diff_locdir_phi"]->Fill(p3FinalReco.phi()-thisSegment->localDirection().phi()); 

  hists_["segm_glbdir_phi"]->Fill(roll->toGlobal(thisSegment->localDirection()).phi()); 
  hists_["trck_glbdir_phi"]->Fill(roll->toGlobal(p3FinalReco).phi()); 
  hists_["diff_glbdir_phi"]->Fill(roll->toGlobal(p3FinalReco).phi()-roll->toGlobal(thisSegment->localDirection()).phi()); 


  if( ( (fabs(thisPosition.x()-r3FinalReco.x())/sigmax ) < MaxPullX ) || (fabs(thisPosition.x()-r3FinalReco.x()) < MaxDiffX ) ) X_MatchFound = true;

  if( MaxPullY < 0. && MaxDiffY < 0. ) Y_MatchFound = true;
  else if( ( (fabs(thisPosition.y()-r3FinalReco.y())/sigmay ) < MaxPullY ) || (fabs(thisPosition.y()-r3FinalReco.y()) < MaxDiffY ) ) Y_MatchFound = true;
  
  if( MaxDiffPhiDir < 0. ) Dir_MatchFound = true;
  //else if( fabs( roll->toGlobal(p3FinalReco).phi() - roll->toGlobal(thisSegment->localDirection()).phi()) < MaxDiffPhiDir) Dir_MatchFound = true; // fix with deltaPhi
  else if( fabs( fabs( roll->toGlobal(p3FinalReco).phi() - roll->toGlobal(thisSegment->localDirection()).phi()) - 3.14159265359 ) < MaxDiffPhiDir) Dir_MatchFound = true; // fix with deltaPhi

  // if(X_MatchFound && Y_MatchFound) { 
  //   std::cout << " ~~~~~~ DeltaPhi:  " << roll->toGlobal(p3FinalReco).phi() << " - " << roll->toGlobal(thisSegment->localDirection()).phi() << " = " << roll->toGlobal(p3FinalReco).phi() - roll->toGlobal(thisSegment->localDirection()).phi() << std::endl; 
  // } 

  return (X_MatchFound && Y_MatchFound && Dir_MatchFound);
}

// define this as a plug-in
DEFINE_FWK_MODULE(RecoMuonTrackProducer);
