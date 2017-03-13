// -*- C++ -*-
//
// Package:    NewME0MuonGeneralAnalyzer
// Class:      NewME0MuonGeneralAnalyzer
// 
/**\class NewME0MuonGeneralAnalyzer 

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/

// System include files
#include <memory>
#include <iomanip>
#include <map>
#include <string>

// Framework and utilities include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// Data formats 
//#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/Provenance/interface/ProductID.h"
//#include "DataFormats/MuonReco/interface/Muon.h"
//#include "DataFormats/MuonReco/interface/MuonFwd.h"
//#include "DataFormats/MuonReco/interface/MuonSelectors.h"
//#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h" 
//#include "DataFormats/DetId/interface/DetId.h"
//#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
//#include "DataFormats/RPCRecHit/interface/RPCRecHit.h"
//#include "DataFormats/MuonDetId/interface/RPCDetId.h"
#include "DataFormats/MuonReco/interface/ME0Muon.h" 
#include "DataFormats/MuonReco/interface/ME0MuonCollection.h" 
#include "DataFormats/GEMRecHit/interface/ME0SegmentCollection.h"
// #include "DataFormats/TrajectoryState/interface/LocalTrajectoryParameters.h"
#include "DataFormats/MuonDetId/interface/ME0DetId.h"
#include "DataFormats/Math/interface/deltaPhi.h" 
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h" 
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h" 
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h" 

// Tracking tools 
// #include "TrackingTools/AnalyticalJacobians/interface/JacobianCartesianToLocal.h"
// #include "TrackingTools/AnalyticalJacobians/interface/JacobianLocalToCartesian.h"

// Muon geometry (some are probably not needed) 
#include "Geometry/GEMGeometry/interface/ME0Geometry.h"
#include "Geometry/GEMGeometry/interface/ME0EtaPartition.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"

// Muon ID 
//#include "RecoMuon/MuonIdentification/interface/ME0MuonSelector.h" 

// Associator by hits or chi2 
#include "SimTracker/TrackAssociatorProducers/plugins/TrackAssociatorByChi2Impl.h"
#include "SimTracker/TrackAssociatorProducers/plugins/TrackAssociatorByHitsImpl.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"

// ROOT includes 
#include "TFile.h"
#include "TH1.h"
//#include "TH2.h"
#include "TGraphAsymmErrors.h"
//#include "TCanvas.h"
#include "TString.h"

//
// class declaration
//

using namespace edm;
using std::vector;
using std::string;


class NewME0MuonGeneralAnalyzer : public edm::EDAnalyzer {
public:
  explicit NewME0MuonGeneralAnalyzer(const edm::ParameterSet&);
  ~NewME0MuonGeneralAnalyzer();

  //typedef std::map<unsigned int, std::vector<unsigned int> > MuonSegmentMap; 

private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  typedef math::XYZTLorentzVectorD FourMomentum; 
  typedef std::map<unsigned int, std::vector<const PSimHit*> > SimHitsMap;
  //typedef std::map<unsigned int, const SimTrack*> SimTracksMap;
  typedef std::map<unsigned int, std::vector<std::pair<unsigned int, unsigned int> > > SimRecoAssoMap;
  typedef std::map<unsigned int, unsigned int> SimRecoMap;

  void associateRecoToSimByHits(/*SimHitsMap&, SimHitsMap&,*/ SimHitsMap&, SimHitsMap&, SimHitsMap&, ME0MuonCollection&,
				edm::Handle<reco::GenParticleCollection>&, edm::Handle<SimTrackContainer>&); 
  void getSimHitsAndTracks(const edm::Event & event, edm::Handle<SimTrackContainer>&, /*SimHitsMap&, SimHitsMap&,*/ SimHitsMap&, SimHitsMap&, SimHitsMap& /*, SimTracksMap&*/); 
  void sortByNumberOfMatches(SimRecoAssoMap&); 

  // 
  // WARNING: all the following maps NEED TO BE CLEANED at the end of each event!!! 
  // 
  //SimRecoAssoMap simToRecoCompleteAssoMap;  // 1 SIM -> many RECO (1 RECO -> many SIM)
  //SimRecoAssoMap simToRecoCompleteTrackerAssoMap;  
  //SimRecoAssoMap simToRecoCompleteME0AssoMap;  

  //SimRecoAssoMap simToRecoAssoMap; // 1 SIM -> many RECO (1 RECO -> 1 SIM)
  //SimRecoAssoMap simToRecoTrackerAssoMap; 
  //SimRecoAssoMap simToRecoME0AssoMap; 

  SimRecoMap simToRecoMap; // 1 SIM  -> 1 RECO
  SimRecoMap recoToSimMap; // 1 RECO -> 1 SIM

  // This one is NOT a sim-to-reco map formally, just a map<int, int> 
  // It maps SimTrack::trackId() <-> index in SimTrackContainer 
  // Useful later, for easy access to SimTrackContainer 
  SimRecoMap trkIdIndexMap; 

  // This one is NOT a sim-to-reco map formally, just a map<int, int> 
  // It maps SimTrack::genpartIndex() <-> index in GenParticleCollection 
  // Useful later, for easy access to GenParticleCollection 
  SimRecoMap genpartIndexMap; 


  bool isGoodME0Muon(edm::ESHandle<ME0Geometry>, const reco::ME0Muon&, double, double, double, double, double); 

  // ----------member data ---------------------------

  // Output file and histos 
  std::string outputFileName; 
  edm::Service<TFileService> outfile_; 
  std::map<std::string, TH1*> hists_; 
  std::map<std::string, TGraphAsymmErrors*> graphs_; 

  // Debug options 
  bool debugPrint; 
  bool outputPrint; 

  // Associators 
  //std::vector<std::string> associatorNames; 
  //std::vector<const TrackAssociatorBase*> associators; 

  // Tokens 
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> > puInfoTag_;
  edm::EDGetTokenT<ME0SegmentCollection> me0SegmentsToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> gpToken_; 
  edm::EDGetTokenT<std::vector<int> > gpidxToken_; 
  //edm::EDGetTokenT<TrackingParticleCollection> tpcToken_; 
  edm::EDGetTokenT<ME0MuonCollection> muonToken_; 
  //edm::EDGetTokenT<reco::MuonCollection> muonToken_; 
  //edm::EDGetTokenT<edm::View<reco::Track> > tracksToken_; 
  //edm::EDGetTokenT<MuonSegmentMap> muSegmMapToken_; 
  //edm::EDGetTokenT<ME0DigiPreRecoCollection> me0digiToken_;
  edm::EDGetTokenT<edm::PSimHitContainer> pxeHiToken; 
  edm::EDGetTokenT<edm::PSimHitContainer> pxeLoToken; 
  edm::EDGetTokenT<edm::PSimHitContainer> me0Token; 
  edm::EDGetTokenT<SimTrackContainer> simTrackToken; 

  // Track-ME0Segment selection cuts 
  double minP; 
  double minPt; 
  double maxPullX; 
  double maxDiffX; 
  double maxPullY; 
  double maxDiffY; 
  double maxDiffPhiDir; 
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
NewME0MuonGeneralAnalyzer::NewME0MuonGeneralAnalyzer(const edm::ParameterSet& iConfig) :
  puInfoTag_( consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("puInfoTag")) ), 
  me0SegmentsToken_( consumes<ME0SegmentCollection>(iConfig.getParameter<edm::InputTag>("me0segmentsTag")) ), 
  //tpcToken_( consumes<TrackingParticleCollection>(iConfig.getParameter<edm::InputTag>("tpCollectionTag")) ), 
  gpToken_( consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("gpCollectionTag")) ), 
  gpidxToken_( consumes<std::vector<int> >(iConfig.getParameter<edm::InputTag>("gpCollectionTag")) ), 
  muonToken_( consumes<ME0MuonCollection>(iConfig.getParameter<edm::InputTag>("muonTag")) ), 
  //muonToken_( consumes<reco::MuonCollection>(iConfig.getParameter<edm::InputTag>("muonTag")) ), 
  //tracksToken_( consumes<edm::View<reco::Track> >(iConfig.getParameter<edm::InputTag>("tracksTag")) ), 
  //muSegmMapToken_( consumes<MuonSegmentMap>(iConfig.getParameter<edm::InputTag>("muSegmMapTag")) ), 
  //me0digiToken_( consumes<ME0DigiPreRecoCollection>(iConfig.getParameter<edm::InputTag>("digisTag")) ) 
  pxeHiToken( consumes<edm::PSimHitContainer>(edm::InputTag("g4SimHits", "TrackerHitsPixelEndcapHighTof")) ), 
  pxeLoToken( consumes<edm::PSimHitContainer>(edm::InputTag("g4SimHits", "TrackerHitsPixelEndcapLowTof")) ), 
  me0Token( consumes<edm::PSimHitContainer>(edm::InputTag("g4SimHits", "MuonME0Hits")) ), 
  simTrackToken( consumes<SimTrackContainer>(edm::InputTag("g4SimHits")) ) 
{
  if(debugPrint) std::cout << "Inside Constructor" << std::endl;

  debugPrint       = iConfig.getUntrackedParameter<bool>("verbose"          , false); 
  outputPrint      = iConfig.getUntrackedParameter<bool>("printOutput"      , false); 

  //associatorNames = iConfig.getParameter<std::vector<std::string> >("associatorNames"); 

  // for (auto const& thisassociator : associatorNames) {
  //   //std::cout << thisassociator << std::endl;
  //   consumes<reco::RecoToSimCollection>(edm::InputTag(thisassociator));
  //   consumes<reco::SimToRecoCollection>(edm::InputTag(thisassociator));
  // } 

  minP          = iConfig.getUntrackedParameter<double>("minP"         , -1.0  ); 
  minPt         = iConfig.getUntrackedParameter<double>("minPt"        , -1.0  ); 

  maxPullX      = iConfig.getUntrackedParameter<double>("maxPullX"     ,  3.0  ); 
  maxDiffX      = iConfig.getUntrackedParameter<double>("maxDiffX"     ,  4.0  ); 
  maxPullY      = iConfig.getUntrackedParameter<double>("maxPullY"     , 20.0  ); 
  maxDiffY      = iConfig.getUntrackedParameter<double>("maxDiffY"     , 20.0  ); 
  maxDiffPhiDir = iConfig.getUntrackedParameter<double>("maxDiffPhiDir",  3.14 ); 

  // -- Binnings -- 
  //double pt_bins[] = {0., 5., 10., 15., 20., 25., 30., 35., 40., 45., 50., 60., 70., 80., 90., 100., 150., 200.}; 
  double pt_bins[] = {0., 1., 2., 3., 4., 5., 7.5, 10., 15., 20., 35., 50., 100., 200.}; 
  int n_pt_bins = sizeof(pt_bins)/sizeof(double); 

  // //double eta_bins[] = {-2.4, -2.1, -1.6, -1.2, -0.9, -0.3, -0.2, 0.2, 0.3, 0.9, 1.2, 1.6, 2.1, 2.4}; 
  // double eta_bins[] = {-3.4, -3.2, -3.0, -2.8, -2.6, -2.4, -2.2, -2.0, -1.8, -1.6, -1.4, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4}; 
  // int n_eta_bins = sizeof(eta_bins)/sizeof(double); 

  double aeta_bins[] = {1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4}; 
  int n_aeta_bins = sizeof(aeta_bins)/sizeof(double); 

  static const double pig=3.1416; 

  // -- Association plots -- 
  hists_["hPxeHi2DHitResiduals"       ] = outfile_->make<TH1F>("hPxeHi2DHitResiduals"       , "hPxeHi2DHitResiduals"       , 200, 0., 8.);
  hists_["hPxeHiXHitResiduals"        ] = outfile_->make<TH1F>("hPxeHiXHitResiduals"        , "hPxeHiXHitResiduals"        , 200, -8., 8.);
  hists_["hPxeHiYHitResiduals"        ] = outfile_->make<TH1F>("hPxeHiYHitResiduals"        , "hPxeHiYHitResiduals"        , 200, -8., 8.);
  // hists_["hPxeHi2DValidHitResiduals"  ] = outfile_->make<TH1F>("hPxeHi2DValidHitResiduals"  , "hPxeHi2DValidHitResiduals"  , 200, 0., 20.);
  // hists_["hPxeHiXValidHitResiduals"   ] = outfile_->make<TH1F>("hPxeHiXValidHitResiduals"   , "hPxeHiXValidHitResiduals"   , 200, -10., 10.);
  // hists_["hPxeHiYValidHitResiduals"   ] = outfile_->make<TH1F>("hPxeHiYValidHitResiduals"   , "hPxeHiYValidHitResiduals"   , 200, -10., 10.);
  // hists_["hPxeHi2DInvalidHitResiduals"] = outfile_->make<TH1F>("hPxeHi2DInvalidHitResiduals", "hPxeHi2DInvalidHitResiduals", 200, 0., 20.);
  // hists_["hPxeHiXInvalidHitResiduals" ] = outfile_->make<TH1F>("hPxeHiXInvalidHitResiduals" , "hPxeHiXInvalidHitResiduals" , 200, -10., 10.);
  // hists_["hPxeHiYInvalidHitResiduals" ] = outfile_->make<TH1F>("hPxeHiYInvalidHitResiduals" , "hPxeHiYInvalidHitResiduals" , 200, -10., 10.);

  hists_["hPxeLo2DHitResiduals"       ] = outfile_->make<TH1F>("hPxeLo2DHitResiduals"       , "hPxeLo2DHitResiduals"       , 200, 0., 8.);
  hists_["hPxeLoXHitResiduals"        ] = outfile_->make<TH1F>("hPxeLoXHitResiduals"        , "hPxeLoXHitResiduals"        , 200, -8., 8.);
  hists_["hPxeLoYHitResiduals"        ] = outfile_->make<TH1F>("hPxeLoYHitResiduals"        , "hPxeLoYHitResiduals"        , 200, -8., 8.);
  // hists_["hPxeLo2DValidHitResiduals"  ] = outfile_->make<TH1F>("hPxeLo2DValidHitResiduals"  , "hPxeLo2DValidHitResiduals"  , 200, 0., 20.);
  // hists_["hPxeLoXValidHitResiduals"   ] = outfile_->make<TH1F>("hPxeLoXValidHitResiduals"   , "hPxeLoXValidHitResiduals"   , 200, -10., 10.);
  // hists_["hPxeLoYValidHitResiduals"   ] = outfile_->make<TH1F>("hPxeLoYValidHitResiduals"   , "hPxeLoYValidHitResiduals"   , 200, -10., 10.);
  // hists_["hPxeLo2DInvalidHitResiduals"] = outfile_->make<TH1F>("hPxeLo2DInvalidHitResiduals", "hPxeLo2DInvalidHitResiduals", 200, 0., 20.);
  // hists_["hPxeLoXInvalidHitResiduals" ] = outfile_->make<TH1F>("hPxeLoXInvalidHitResiduals" , "hPxeLoXInvalidHitResiduals" , 200, -10., 10.);
  // hists_["hPxeLoYInvalidHitResiduals" ] = outfile_->make<TH1F>("hPxeLoYInvalidHitResiduals" , "hPxeLoYInvalidHitResiduals" , 200, -10., 10.);

  hists_["hMe0Hi2DHitResiduals"       ] = outfile_->make<TH1F>("hMe0Hi2DHitResiduals"       , "hMe0Hi2DHitResiduals"       , 200, 0., 20.);
  hists_["hMe0HiXHitResiduals"        ] = outfile_->make<TH1F>("hMe0HiXHitResiduals"        , "hMe0HiXHitResiduals"        , 200, -20., 20.);
  hists_["hMe0HiYHitResiduals"        ] = outfile_->make<TH1F>("hMe0HiYHitResiduals"        , "hMe0HiYHitResiduals"        , 200, -20., 20.);
  // hists_["hMe0Hi2DValidHitResiduals"  ] = outfile_->make<TH1F>("hMe0Hi2DValidHitResiduals"  , "hMe0Hi2DValidHitResiduals"  , 200, 0., 200.);
  // hists_["hMe0HiXValidHitResiduals"   ] = outfile_->make<TH1F>("hMe0HiXValidHitResiduals"   , "hMe0HiXValidHitResiduals"   , 200, -100., 100.);
  // hists_["hMe0HiYValidHitResiduals"   ] = outfile_->make<TH1F>("hMe0HiYValidHitResiduals"   , "hMe0HiYValidHitResiduals"   , 200, -100., 100.);
  // hists_["hMe0Hi2DInvalidHitResiduals"] = outfile_->make<TH1F>("hMe0Hi2DInvalidHitResiduals", "hMe0Hi2DInvalidHitResiduals", 200, 0., 200.);
  // hists_["hMe0HiXInvalidHitResiduals" ] = outfile_->make<TH1F>("hMe0HiXInvalidHitResiduals" , "hMe0HiXInvalidHitResiduals" , 200, -100., 100.);
  // hists_["hMe0HiYInvalidHitResiduals" ] = outfile_->make<TH1F>("hMe0HiYInvalidHitResiduals" , "hMe0HiYInvalidHitResiduals" , 200, -100., 100.);


  hists_["numberMatchesPerTrack_tracker_all"] = outfile_->make<TH1F>("numberMatchesPerTrack_tracker_all", "numberMatchesPerTrack_tracker_all",  31, -0.5, 30.5); 
  hists_["fractnMatchesPerTrack_tracker_all"] = outfile_->make<TH1F>("fractnMatchesPerTrack_tracker_all", "fractnMatchesPerTrack_tracker_all", 100,  0.0,  1.0); 

  hists_["numberMatchesPerTrack_me0segm_all"] = outfile_->make<TH1F>("numberMatchesPerTrack_me0segm_all", "numberMatchesPerTrack_me0segm_all",  31, -0.5, 30.5); 
  hists_["fractnMatchesPerTrack_me0segm_all"] = outfile_->make<TH1F>("fractnMatchesPerTrack_me0segm_all", "fractnMatchesPerTrack_me0segm_all", 100,  0.0,  1.0); 

  hists_["numberMatchesPerTrack_wholemu_all"] = outfile_->make<TH1F>("numberMatchesPerTrack_wholemu_all", "numberMatchesPerTrack_wholemu_all",  31, -0.5, 30.5); 
  hists_["fractnMatchesPerTrack_wholemu_all"] = outfile_->make<TH1F>("fractnMatchesPerTrack_wholemu_all", "fractnMatchesPerTrack_wholemu_all", 100,  0.0,  1.0); 

  hists_["numberMatchesPerTrack_tracker_best"] = outfile_->make<TH1F>("numberMatchesPerTrack_tracker_best", "numberMatchesPerTrack_tracker_best",  31, -0.5, 30.5); 
  hists_["fractnMatchesPerTrack_tracker_best"] = outfile_->make<TH1F>("fractnMatchesPerTrack_tracker_best", "fractnMatchesPerTrack_tracker_best", 100,  0.0,  1.0); 

  hists_["numberMatchesPerTrack_me0segm_best"] = outfile_->make<TH1F>("numberMatchesPerTrack_me0segm_best", "numberMatchesPerTrack_me0segm_best",  31, -0.5, 30.5); 
  hists_["fractnMatchesPerTrack_me0segm_best"] = outfile_->make<TH1F>("fractnMatchesPerTrack_me0segm_best", "fractnMatchesPerTrack_me0segm_best", 100,  0.0,  1.0); 

  hists_["numberMatchesPerTrack_wholemu_best"] = outfile_->make<TH1F>("numberMatchesPerTrack_wholemu_best", "numberMatchesPerTrack_wholemu_best",  31, -0.5, 30.5); 
  hists_["fractnMatchesPerTrack_wholemu_best"] = outfile_->make<TH1F>("fractnMatchesPerTrack_wholemu_best", "fractnMatchesPerTrack_wholemu_best", 100,  0.0,  1.0); 

  // Signal only 
  hists_["numberMatchesPerTrack_sign_tracker_all"] = outfile_->make<TH1F>("numberMatchesPerTrack_sign_tracker_all", "numberMatchesPerTrack_sign_tracker_all",  31, -0.5, 30.5); 
  hists_["fractnMatchesPerTrack_sign_tracker_all"] = outfile_->make<TH1F>("fractnMatchesPerTrack_sign_tracker_all", "fractnMatchesPerTrack_sign_tracker_all", 100,  0.0,  1.0); 

  hists_["numberMatchesPerTrack_sign_me0segm_all"] = outfile_->make<TH1F>("numberMatchesPerTrack_sign_me0segm_all", "numberMatchesPerTrack_sign_me0segm_all",  31, -0.5, 30.5); 
  hists_["fractnMatchesPerTrack_sign_me0segm_all"] = outfile_->make<TH1F>("fractnMatchesPerTrack_sign_me0segm_all", "fractnMatchesPerTrack_sign_me0segm_all", 100,  0.0,  1.0); 

  hists_["numberMatchesPerTrack_sign_wholemu_all"] = outfile_->make<TH1F>("numberMatchesPerTrack_sign_wholemu_all", "numberMatchesPerTrack_sign_wholemu_all",  31, -0.5, 30.5); 
  hists_["fractnMatchesPerTrack_sign_wholemu_all"] = outfile_->make<TH1F>("fractnMatchesPerTrack_sign_wholemu_all", "fractnMatchesPerTrack_sign_wholemu_all", 100,  0.0,  1.0); 

  hists_["numberMatchesPerTrack_sign_tracker_best"] = outfile_->make<TH1F>("numberMatchesPerTrack_sign_tracker_best", "numberMatchesPerTrack_sign_tracker_best",  31, -0.5, 30.5); 
  hists_["fractnMatchesPerTrack_sign_tracker_best"] = outfile_->make<TH1F>("fractnMatchesPerTrack_sign_tracker_best", "fractnMatchesPerTrack_sign_tracker_best", 100,  0.0,  1.0); 

  hists_["numberMatchesPerTrack_sign_me0segm_best"] = outfile_->make<TH1F>("numberMatchesPerTrack_sign_me0segm_best", "numberMatchesPerTrack_sign_me0segm_best",  31, -0.5, 30.5); 
  hists_["fractnMatchesPerTrack_sign_me0segm_best"] = outfile_->make<TH1F>("fractnMatchesPerTrack_sign_me0segm_best", "fractnMatchesPerTrack_sign_me0segm_best", 100,  0.0,  1.0); 

  hists_["numberMatchesPerTrack_sign_wholemu_best"] = outfile_->make<TH1F>("numberMatchesPerTrack_sign_wholemu_best", "numberMatchesPerTrack_sign_wholemu_best",  31, -0.5, 30.5); 
  hists_["fractnMatchesPerTrack_sign_wholemu_best"] = outfile_->make<TH1F>("fractnMatchesPerTrack_sign_wholemu_best", "fractnMatchesPerTrack_sign_wholemu_best", 100,  0.0,  1.0); 


  // -- Efficiencies -- 
  hists_["eff_den_pt"      ] = outfile_->make<TH1F>("eff_den_pt"    , "eff_den_pt"    , n_pt_bins-1  , pt_bins  ); 
  // hists_["eff_den_eta"     ] = outfile_->make<TH1F>("eff_den_eta"   , "eff_den_eta"   , n_eta_bins-1 , eta_bins ); 
  hists_["eff_den_aeta"    ] = outfile_->make<TH1F>("eff_den_aeta"  , "eff_den_aeta"  , n_aeta_bins-1, aeta_bins); 
  hists_["eff_den_phi"     ] = outfile_->make<TH1F>("eff_den_phi"   , "eff_den_phi"   , 10, -pig   , pig     ); 
  // hists_["eff_den_truepu"  ] = outfile_->make<TH1F>("eff_den_truepu", "eff_den_truepu", 300, -0.5    , 299.5    ); 
  // hists_["eff_den_itpu"    ] = outfile_->make<TH1F>("eff_den_itpu"  , "eff_den_itpu"  , 300, -0.5    , 299.5    ); 
  hists_["eff_den_allpu"   ] = outfile_->make<TH1F>("eff_den_allpu" , "eff_den_allpu" , 300, -0.5    , 299.5    ); 

  hists_["eff_num_pt"      ] = outfile_->make<TH1F>("eff_num_pt"    , "eff_num_pt"    , n_pt_bins-1  , pt_bins  ); 
  // hists_["eff_num_eta"     ] = outfile_->make<TH1F>("eff_num_eta"   , "eff_num_eta"   , n_eta_bins-1 , eta_bins ); 
  hists_["eff_num_aeta"    ] = outfile_->make<TH1F>("eff_num_aeta"  , "eff_num_aeta"  , n_aeta_bins-1, aeta_bins); 
  hists_["eff_num_phi"     ] = outfile_->make<TH1F>("eff_num_phi"   , "eff_num_phi"   , 10, -pig   , pig     ); 
  // hists_["eff_num_truepu"  ] = outfile_->make<TH1F>("eff_num_truepu", "eff_num_truepu", 300, -0.5    , 299.5    ); 
  // hists_["eff_num_itpu"    ] = outfile_->make<TH1F>("eff_num_itpu"  , "eff_num_itpu"  , 300, -0.5    , 299.5    ); 
  hists_["eff_num_allpu"   ] = outfile_->make<TH1F>("eff_num_allpu" , "eff_num_allpu" , 300, -0.5    , 299.5    ); 

 //  hists_["eff_pt"          ] = outfile_->make<TH1F>("eff_pt"        , "eff_pt"        , n_pt_bins-1  , pt_bins  ); 
 // //  hists_["eff_eta"         ] = outfile_->make<TH1F>("eff_eta"       , "eff_eta"       , n_eta_bins-1 , eta_bins ); 
 //  hists_["eff_aeta"        ] = outfile_->make<TH1F>("eff_aeta"      , "eff_aeta"      , n_aeta_bins-1, aeta_bins); 
 //  hists_["eff_phi"         ] = outfile_->make<TH1F>("eff_phi"       , "eff_phi"       , 10, -pig   , pig     ); 
 //  // hists_["eff_truepu"      ] = outfile_->make<TH1F>("eff_truepu"    , "eff_truepu"    , 300, -0.5    , 299.5    ); 
 //  // hists_["eff_itpu"        ] = outfile_->make<TH1F>("eff_itpu"      , "eff_itpu"      , 300, -0.5    , 299.5    ); 
 //  hists_["eff_allpu"       ] = outfile_->make<TH1F>("eff_allpu"     , "eff_allpu"     , 300, -0.5    , 299.5    ); 

  graphs_["eff_pt_err"     ] = outfile_->make<TGraphAsymmErrors>(); 
  // graphs_["eff_eta_err"    ] = outfile_->make<TGraphAsymmErrors>(); 
  graphs_["eff_aeta_err"   ] = outfile_->make<TGraphAsymmErrors>(); 
  graphs_["eff_phi_err"    ] = outfile_->make<TGraphAsymmErrors>(); 
  // graphs_["eff_truepu_err" ] = outfile_->make<TGraphAsymmErrors>(); 
  // graphs_["eff_itpu_err"   ] = outfile_->make<TGraphAsymmErrors>(); 
  graphs_["eff_allpu_err"  ] = outfile_->make<TGraphAsymmErrors>(); 

  graphs_["eff_pt_err"     ]->SetNameTitle("eff_pt_err"    , "eff_pt_err"    ); 
  // graphs_["eff_eta_err"    ]->SetNameTitle("eff_eta_err"   , "eff_eta_err"   ); 
  graphs_["eff_aeta_err"   ]->SetNameTitle("eff_aeta_err"  , "eff_aeta_err"  ); 
  graphs_["eff_phi_err"    ]->SetNameTitle("eff_phi_err"   , "eff_phi_err"   ); 
  // graphs_["eff_truepu_err" ]->SetNameTitle("eff_truepu_err", "eff_truepu_err"); 
  // graphs_["eff_itpu_err"   ]->SetNameTitle("eff_itpu_err"  , "eff_itpu_err"  ); 
  graphs_["eff_allpu_err"  ]->SetNameTitle("eff_allpu_err" , "eff_allpu_err" ); 

  // -- Fakes -- 
  // hists_["fake_pertrack_den_pt"     ] = outfile_->make<TH1F>("fake_pertrack_den_pt"    , "fake_pertrack_den_pt"    , n_pt_bins-1  , pt_bins  );
  // // hists_["fake_pertrack_den_eta"    ] = outfile_->make<TH1F>("fake_pertrack_den_eta"   , "fake_pertrack_den_eta"   , n_eta_bins-1 , eta_bins );
  // hists_["fake_pertrack_den_aeta"   ] = outfile_->make<TH1F>("fake_pertrack_den_aeta"  , "fake_pertrack_den_aeta"  , n_aeta_bins-1, aeta_bins);
  // hists_["fake_pertrack_den_phi"    ] = outfile_->make<TH1F>("fake_pertrack_den_phi"   , "fake_pertrack_den_phi"   , 10, -pig   , pig     );
  // // hists_["fake_pertrack_den_truepu" ] = outfile_->make<TH1F>("fake_pertrack_den_truepu", "fake_pertrack_den_truepu", 300, -0.5    , 299.5    ); 
  // // hists_["fake_pertrack_den_itpu"   ] = outfile_->make<TH1F>("fake_pertrack_den_itpu"  , "fake_pertrack_den_itpu"  , 300, -0.5    , 299.5    ); 
  // hists_["fake_pertrack_den_allpu"  ] = outfile_->make<TH1F>("fake_pertrack_den_allpu" , "fake_pertrack_den_allpu" , 300, -0.5    , 299.5    ); 

  hists_["fake_den_pt"     ] = outfile_->make<TH1F>("fake_den_pt"    , "fake_den_pt"    , n_pt_bins-1  , pt_bins  );
  // hists_["fake_den_eta"    ] = outfile_->make<TH1F>("fake_den_eta"   , "fake_den_eta"   , n_eta_bins-1 , eta_bins );
  hists_["fake_den_aeta"   ] = outfile_->make<TH1F>("fake_den_aeta"  , "fake_den_aeta"  , n_aeta_bins-1, aeta_bins);
  hists_["fake_den_phi"    ] = outfile_->make<TH1F>("fake_den_phi"   , "fake_den_phi"   , 10, -pig   , pig     );
  // hists_["fake_den_truepu" ] = outfile_->make<TH1F>("fake_den_truepu", "fake_den_truepu", 300, -0.5    , 299.5    ); 
  // hists_["fake_den_itpu"   ] = outfile_->make<TH1F>("fake_den_itpu"  , "fake_den_itpu"  , 300, -0.5    , 299.5    ); 
  hists_["fake_den_allpu"  ] = outfile_->make<TH1F>("fake_den_allpu" , "fake_den_allpu" , 300, -0.5    , 299.5    ); 

  // hists_["fake_pertrack_num_pt"     ] = outfile_->make<TH1F>("fake_pertrack_num_pt"    , "fake_pertrack_num_pt"    , n_pt_bins-1  , pt_bins  );
  // // hists_["fake_pertrack_num_eta"    ] = outfile_->make<TH1F>("fake_pertrack_num_eta"   , "fake_pertrack_num_eta"   , n_eta_bins-1 , eta_bins );
  // hists_["fake_pertrack_num_aeta"   ] = outfile_->make<TH1F>("fake_pertrack_num_aeta"  , "fake_pertrack_num_aeta"  , n_aeta_bins-1, aeta_bins);
  // hists_["fake_pertrack_num_phi"    ] = outfile_->make<TH1F>("fake_pertrack_num_phi"   , "fake_pertrack_num_phi"   , 10, -pig   , pig     );
  // // hists_["fake_pertrack_num_truepu" ] = outfile_->make<TH1F>("fake_pertrack_num_truepu", "fake_pertrack_num_truepu", 300, -0.5    , 299.5    ); 
  // // hists_["fake_pertrack_num_itpu"   ] = outfile_->make<TH1F>("fake_pertrack_num_itpu"  , "fake_pertrack_num_itpu"  , 300, -0.5    , 299.5    ); 
  // hists_["fake_pertrack_num_allpu"  ] = outfile_->make<TH1F>("fake_pertrack_num_allpu" , "fake_pertrack_num_allpu" , 300, -0.5    , 299.5    ); 

  hists_["fake_num_pt"     ] = outfile_->make<TH1F>("fake_num_pt"    , "fake_num_pt"    , n_pt_bins-1  , pt_bins  );
  // hists_["fake_num_eta"    ] = outfile_->make<TH1F>("fake_num_eta"   , "fake_num_eta"   , n_eta_bins-1 , eta_bins );
  hists_["fake_num_aeta"   ] = outfile_->make<TH1F>("fake_num_aeta"  , "fake_num_aeta"  , n_aeta_bins-1, aeta_bins);
  hists_["fake_num_phi"    ] = outfile_->make<TH1F>("fake_num_phi"   , "fake_num_phi"   , 10, -pig   , pig     );
  // hists_["fake_num_truepu" ] = outfile_->make<TH1F>("fake_num_truepu", "fake_num_truepu", 300, -0.5    , 299.5    ); 
  // hists_["fake_num_itpu"   ] = outfile_->make<TH1F>("fake_num_itpu"  , "fake_num_itpu"  , 300, -0.5    , 299.5    ); 
  hists_["fake_num_allpu"  ] = outfile_->make<TH1F>("fake_num_allpu" , "fake_num_allpu" , 300, -0.5    , 299.5    ); 

  // hists_["fake_pt"         ] = outfile_->make<TH1F>("fake_pt"        , "fake_pt"        , n_pt_bins-1  , pt_bins  );
  // // hists_["fake_eta"        ] = outfile_->make<TH1F>("fake_eta"       , "fake_eta"       , n_eta_bins-1 , eta_bins );
  // hists_["fake_aeta"       ] = outfile_->make<TH1F>("fake_aeta"      , "fake_aeta"      , n_aeta_bins-1, aeta_bins);
  // hists_["fake_phi"        ] = outfile_->make<TH1F>("fake_phi"       , "fake_phi"       , 10, -pig   , pig     );
  // // hists_["fake_truepu"     ] = outfile_->make<TH1F>("fake_truepu"    , "fake_truepu"    , 300, -0.5    , 299.5    ); 
  // // hists_["fake_itpu"       ] = outfile_->make<TH1F>("fake_itpu"      , "fake_itpu"      , 300, -0.5    , 299.5    ); 
  // hists_["fake_allpu"      ] = outfile_->make<TH1F>("fake_allpu"     , "fake_allpu"     , 300, -0.5    , 299.5    ); 

  // graphs_["fake_pertrack_pt_err"    ] = outfile_->make<TGraphAsymmErrors>(); 
  // // graphs_["fake_pertrack_eta_err"   ] = outfile_->make<TGraphAsymmErrors>(); 
  // graphs_["fake_pertrack_aeta_err"  ] = outfile_->make<TGraphAsymmErrors>(); 
  // graphs_["fake_pertrack_phi_err"   ] = outfile_->make<TGraphAsymmErrors>(); 
  // // graphs_["fake_pertrack_truepu_err"] = outfile_->make<TGraphAsymmErrors>(); 
  // // graphs_["fake_pertrack_itpu_err"  ] = outfile_->make<TGraphAsymmErrors>(); 
  // graphs_["fake_pertrack_allpu_err" ] = outfile_->make<TGraphAsymmErrors>(); 

  // graphs_["fake_pertrack_pt_err"    ]->SetNameTitle("fake_pertrack_pt_err"    , "fake_pertrack_pt_err"    ); 
  // // graphs_["fake_pertrack_eta_err"   ]->SetNameTitle("fake_pertrack_eta_err"   , "fake_pertrack_eta_err"   ); 
  // graphs_["fake_pertrack_aeta_err"  ]->SetNameTitle("fake_pertrack_aeta_err"  , "fake_pertrack_aeta_err"  ); 
  // graphs_["fake_pertrack_phi_err"   ]->SetNameTitle("fake_pertrack_phi_err"   , "fake_pertrack_phi_err"   ); 
  // // graphs_["fake_pertrack_truepu_err"]->SetNameTitle("fake_pertrack_truepu_err", "fake_pertrack_truepu_err"); 
  // // graphs_["fake_pertrack_itpu_err"  ]->SetNameTitle("fake_pertrack_itpu_err"  , "fake_pertrack_itpu_err"  ); 
  // graphs_["fake_pertrack_allpu_err" ]->SetNameTitle("fake_pertrack_allpu_err" , "fake_pertrack_allpu_err" ); 

  graphs_["fake_pt_err"    ] = outfile_->make<TGraphAsymmErrors>(); 
  // graphs_["fake_eta_err"   ] = outfile_->make<TGraphAsymmErrors>(); 
  graphs_["fake_aeta_err"  ] = outfile_->make<TGraphAsymmErrors>(); 
  graphs_["fake_phi_err"   ] = outfile_->make<TGraphAsymmErrors>(); 
  // graphs_["fake_truepu_err"] = outfile_->make<TGraphAsymmErrors>(); 
  // graphs_["fake_itpu_err"  ] = outfile_->make<TGraphAsymmErrors>(); 
  graphs_["fake_allpu_err" ] = outfile_->make<TGraphAsymmErrors>(); 

  graphs_["fake_pt_err"    ]->SetNameTitle("fake_pt_err"    , "fake_pt_err"    ); 
  // graphs_["fake_eta_err"   ]->SetNameTitle("fake_eta_err"   , "fake_eta_err"   ); 
  graphs_["fake_aeta_err"  ]->SetNameTitle("fake_aeta_err"  , "fake_aeta_err"  ); 
  graphs_["fake_phi_err"   ]->SetNameTitle("fake_phi_err"   , "fake_phi_err"   ); 
  // graphs_["fake_truepu_err"]->SetNameTitle("fake_truepu_err", "fake_truepu_err"); 
  // graphs_["fake_itpu_err"  ]->SetNameTitle("fake_itpu_err"  , "fake_itpu_err"  ); 
  graphs_["fake_allpu_err" ]->SetNameTitle("fake_allpu_err" , "fake_allpu_err" ); 

  // -- Background (includes real muons, e.g. from decays in flight) -- 
  // Denominator not needed, same as for fakes 
  // hists_["bkg_den_pt"     ] = outfile_->make<TH1F>("bkg_den_pt"    , "bkg_den_pt"    , n_pt_bins-1  , pt_bins  );
  // // hists_["bkg_den_eta"    ] = outfile_->make<TH1F>("bkg_den_eta"   , "bkg_den_eta"   , n_eta_bins-1 , eta_bins );
  // hists_["bkg_den_aeta"   ] = outfile_->make<TH1F>("bkg_den_aeta"  , "bkg_den_aeta"  , n_aeta_bins-1, aeta_bins);
  // hists_["bkg_den_phi"    ] = outfile_->make<TH1F>("bkg_den_phi"   , "bkg_den_phi"   , 10, -pig   , pig     );
  // // hists_["bkg_den_truepu" ] = outfile_->make<TH1F>("bkg_den_truepu", "bkg_den_truepu", 300, -0.5    , 299.5    ); 
  // // hists_["bkg_den_itpu"   ] = outfile_->make<TH1F>("bkg_den_itpu"  , "bkg_den_itpu"  , 300, -0.5    , 299.5    ); 
  // hists_["bkg_den_allpu"  ] = outfile_->make<TH1F>("bkg_den_allpu" , "bkg_den_allpu" , 300, -0.5    , 299.5    ); 

  // hists_["bkg_pertrack_num_pt"     ] = outfile_->make<TH1F>("bkg_pertrack_num_pt"    , "bkg_pertrack_num_pt"    , n_pt_bins-1  , pt_bins  );
  // // hists_["bkg_pertrack_num_eta"    ] = outfile_->make<TH1F>("bkg_pertrack_num_eta"   , "bkg_pertrack_num_eta"   , n_eta_bins-1 , eta_bins );
  // hists_["bkg_pertrack_num_aeta"   ] = outfile_->make<TH1F>("bkg_pertrack_num_aeta"  , "bkg_pertrack_num_aeta"  , n_aeta_bins-1, aeta_bins);
  // hists_["bkg_pertrack_num_phi"    ] = outfile_->make<TH1F>("bkg_pertrack_num_phi"   , "bkg_pertrack_num_phi"   , 10, -pig   , pig     );
  // // hists_["bkg_pertrack_num_truepu" ] = outfile_->make<TH1F>("bkg_pertrack_num_truepu", "bkg_pertrack_num_truepu", 300, -0.5    , 299.5    ); 
  // // hists_["bkg_pertrack_num_itpu"   ] = outfile_->make<TH1F>("bkg_pertrack_num_itpu"  , "bkg_pertrack_num_itpu"  , 300, -0.5    , 299.5    ); 
  // hists_["bkg_pertrack_num_allpu"  ] = outfile_->make<TH1F>("bkg_pertrack_num_allpu" , "bkg_pertrack_num_allpu" , 300, -0.5    , 299.5    ); 

  hists_["bkg_num_pt"     ] = outfile_->make<TH1F>("bkg_num_pt"    , "bkg_num_pt"    , n_pt_bins-1  , pt_bins  );
  // hists_["bkg_num_eta"    ] = outfile_->make<TH1F>("bkg_num_eta"   , "bkg_num_eta"   , n_eta_bins-1 , eta_bins );
  hists_["bkg_num_aeta"   ] = outfile_->make<TH1F>("bkg_num_aeta"  , "bkg_num_aeta"  , n_aeta_bins-1, aeta_bins);
  hists_["bkg_num_phi"    ] = outfile_->make<TH1F>("bkg_num_phi"   , "bkg_num_phi"   , 10, -pig   , pig     );
  // hists_["bkg_num_truepu" ] = outfile_->make<TH1F>("bkg_num_truepu", "bkg_num_truepu", 300, -0.5    , 299.5    ); 
  // hists_["bkg_num_itpu"   ] = outfile_->make<TH1F>("bkg_num_itpu"  , "bkg_num_itpu"  , 300, -0.5    , 299.5    ); 
  hists_["bkg_num_allpu"  ] = outfile_->make<TH1F>("bkg_num_allpu" , "bkg_num_allpu" , 300, -0.5    , 299.5    ); 

  // hists_["bkg_pt"         ] = outfile_->make<TH1F>("bkg_pt"        , "bkg_pt"        , n_pt_bins-1  , pt_bins  );
  // // hists_["bkg_eta"        ] = outfile_->make<TH1F>("bkg_eta"       , "bkg_eta"       , n_eta_bins-1 , eta_bins );
  // hists_["bkg_aeta"       ] = outfile_->make<TH1F>("bkg_aeta"      , "bkg_aeta"      , n_aeta_bins-1, aeta_bins);
  // hists_["bkg_phi"        ] = outfile_->make<TH1F>("bkg_phi"       , "bkg_phi"       , 10, -pig   , pig     );
  // // hists_["bkg_truepu"     ] = outfile_->make<TH1F>("bkg_truepu"    , "bkg_truepu"    , 300, -0.5    , 299.5    ); 
  // // hists_["bkg_itpu"       ] = outfile_->make<TH1F>("bkg_itpu"      , "bkg_itpu"      , 300, -0.5    , 299.5    ); 
  // hists_["bkg_allpu"      ] = outfile_->make<TH1F>("bkg_allpu"     , "bkg_allpu"     , 300, -0.5    , 299.5    ); 

  // graphs_["bkg_pertrack_pt_err"    ] = outfile_->make<TGraphAsymmErrors>(); 
  // // graphs_["bkg_pertrack_eta_err"   ] = outfile_->make<TGraphAsymmErrors>(); 
  // graphs_["bkg_pertrack_aeta_err"  ] = outfile_->make<TGraphAsymmErrors>(); 
  // graphs_["bkg_pertrack_phi_err"   ] = outfile_->make<TGraphAsymmErrors>(); 
  // // graphs_["bkg_pertrack_truepu_err"] = outfile_->make<TGraphAsymmErrors>(); 
  // // graphs_["bkg_pertrack_itpu_err"  ] = outfile_->make<TGraphAsymmErrors>(); 
  // graphs_["bkg_pertrack_allpu_err" ] = outfile_->make<TGraphAsymmErrors>(); 

  // graphs_["bkg_pertrack_pt_err"    ]->SetNameTitle("bkg_pertrack_pt_err"    , "bkg_pertrack_pt_err"    ); 
  // // graphs_["bkg_pertrack_eta_err"   ]->SetNameTitle("bkg_pertrack_eta_err"   , "bkg_pertrack_eta_err"   ); 
  // graphs_["bkg_pertrack_aeta_err"  ]->SetNameTitle("bkg_pertrack_aeta_err"  , "bkg_pertrack_aeta_err"  ); 
  // graphs_["bkg_pertrack_phi_err"   ]->SetNameTitle("bkg_pertrack_phi_err"   , "bkg_pertrack_phi_err"   ); 
  // // graphs_["bkg_pertrack_truepu_err"]->SetNameTitle("bkg_pertrack_truepu_err", "bkg_pertrack_truepu_err"); 
  // // graphs_["bkg_pertrack_itpu_err"  ]->SetNameTitle("bkg_pertrack_itpu_err"  , "bkg_pertrack_itpu_err"  ); 
  // graphs_["bkg_pertrack_allpu_err" ]->SetNameTitle("bkg_pertrack_allpu_err" , "bkg_pertrack_allpu_err" ); 

  graphs_["bkg_pt_err"    ] = outfile_->make<TGraphAsymmErrors>(); 
  // graphs_["bkg_eta_err"   ] = outfile_->make<TGraphAsymmErrors>(); 
  graphs_["bkg_aeta_err"  ] = outfile_->make<TGraphAsymmErrors>(); 
  graphs_["bkg_phi_err"   ] = outfile_->make<TGraphAsymmErrors>(); 
  // graphs_["bkg_truepu_err"] = outfile_->make<TGraphAsymmErrors>(); 
  // graphs_["bkg_itpu_err"  ] = outfile_->make<TGraphAsymmErrors>(); 
  graphs_["bkg_allpu_err" ] = outfile_->make<TGraphAsymmErrors>(); 

  graphs_["bkg_pt_err"    ]->SetNameTitle("bkg_pt_err"    , "bkg_pt_err"    ); 
  // graphs_["bkg_eta_err"   ]->SetNameTitle("bkg_eta_err"   , "bkg_eta_err"   ); 
  graphs_["bkg_aeta_err"  ]->SetNameTitle("bkg_aeta_err"  , "bkg_aeta_err"  ); 
  graphs_["bkg_phi_err"   ]->SetNameTitle("bkg_phi_err"   , "bkg_phi_err"   ); 
  // graphs_["bkg_truepu_err"]->SetNameTitle("bkg_truepu_err", "bkg_truepu_err"); 
  // graphs_["bkg_itpu_err"  ]->SetNameTitle("bkg_itpu_err"  , "bkg_itpu_err"  ); 
  graphs_["bkg_allpu_err" ]->SetNameTitle("bkg_allpu_err" , "bkg_allpu_err" ); 


  /// --- Rates per event or per segment 
  // Fake rates per event 
  // hists_["fake_pertrack_pt_perevt"    ] = outfile_->make<TH1F>();
  // //hists_["fake_pertrack_eta_perevt"   ] = outfile_->make<TH1F>();
  // hists_["fake_pertrack_aeta_perevt"  ] = outfile_->make<TH1F>();
  // hists_["fake_pertrack_phi_perevt"   ] = outfile_->make<TH1F>();
  // //hists_["fake_pertrack_truepu_perevt"] = outfile_->make<TH1F>(); 
  // //hists_["fake_pertrack_itpu_perevt"  ] = outfile_->make<TH1F>(); 
  // hists_["fake_pertrack_allpu_perevt" ] = outfile_->make<TH1F>(); 

  hists_["fake_pt_perevt"    ] = outfile_->make<TH1F>("fake_pt_perevt"    , "fake_pt_perevt"    , n_pt_bins-1  , pt_bins     );
  //hists_["fake_eta_perevt"   ] = outfile_->make<TH1F>("fake_eta_perevt"   , "fake_eta_perevt"   , n_eta_bins-1 , eta_bins    );
  hists_["fake_aeta_perevt"  ] = outfile_->make<TH1F>("fake_aeta_perevt"  , "fake_aeta_perevt"  , n_aeta_bins-1, aeta_bins   );
  hists_["fake_phi_perevt"   ] = outfile_->make<TH1F>("fake_phi_perevt"   , "fake_phi_perevt"   , 10           , -pig, pig   );
  //hists_["fake_truepu_perevt"] = outfile_->make<TH1F>("fake_truepu_perevt", "fake_truepu_perevt", 300          , -0.5, 299.5 ); 
  //hists_["fake_itpu_perevt"  ] = outfile_->make<TH1F>("fake_itpu_perevt"  , "fake_itpu_perevt"  , 300          , -0.5, 299.5 ); 
  hists_["fake_allpu_perevt" ] = outfile_->make<TH1F>("fake_allpu_perevt" , "fake_allpu_perevt" , 300          , -0.5, 299.5 ); 

  // Fake rates per segment 
  // hists_["fake_pertrack_pt_perseg"    ] = outfile_->make<TH1F>();
  // //hists_["fake_pertrack_eta_perseg"   ] = outfile_->make<TH1F>();
  // hists_["fake_pertrack_aeta_perseg"  ] = outfile_->make<TH1F>();
  // hists_["fake_pertrack_phi_perseg"   ] = outfile_->make<TH1F>();
  // //hists_["fake_pertrack_truepu_perseg"] = outfile_->make<TH1F>(); 
  // //hists_["fake_pertrack_itpu_perseg"  ] = outfile_->make<TH1F>(); 
  // hists_["fake_pertrack_allpu_perseg" ] = outfile_->make<TH1F>(); 

  hists_["fake_pt_perseg"    ] = outfile_->make<TH1F>("fake_pt_perseg"    , "fake_pt_perseg"    , n_pt_bins-1  , pt_bins     );
  //hists_["fake_eta_perseg"   ] = outfile_->make<TH1F>("fake_eta_perseg"   , "fake_eta_perseg"   , n_eta_bins-1 , eta_bins    );
  hists_["fake_aeta_perseg"  ] = outfile_->make<TH1F>("fake_aeta_perseg"  , "fake_aeta_perseg"  , n_aeta_bins-1, aeta_bins   );
  hists_["fake_phi_perseg"   ] = outfile_->make<TH1F>("fake_phi_perseg"   , "fake_phi_perseg"   , 10           , -pig, pig   );
  //hists_["fake_truepu_perseg"] = outfile_->make<TH1F>("fake_truepu_perseg", "fake_truepu_perseg", 300          , -0.5, 299.5 ); 
  //hists_["fake_itpu_perseg"  ] = outfile_->make<TH1F>("fake_itpu_perseg"  , "fake_itpu_perseg"  , 300          , -0.5, 299.5 ); 
  hists_["fake_allpu_perseg" ] = outfile_->make<TH1F>("fake_allpu_perseg" , "fake_allpu_perseg" , 300          , -0.5, 299.5 ); 

  // Background rates per event 
  // hists_["bkg_pertrack_pt_perevt"     ] = outfile_->make<TH1F>();
  // //hists_["bkg_pertrack_eta_perevt"    ] = outfile_->make<TH1F>();
  // hists_["bkg_pertrack_aeta_perevt"   ] = outfile_->make<TH1F>();
  // hists_["bkg_pertrack_phi_perevt"    ] = outfile_->make<TH1F>();
  // //hists_["bkg_pertrack_truepu_perevt" ] = outfile_->make<TH1F>(); 
  // //hists_["bkg_pertrack_itpu_perevt"   ] = outfile_->make<TH1F>(); 
  // hists_["bkg_pertrack_allpu_perevt"  ] = outfile_->make<TH1F>(); 

  hists_["bkg_pt_perevt"    ] = outfile_->make<TH1F>("bkg_pt_perevt"    , "bkg_pt_perevt"    , n_pt_bins-1  , pt_bins     );
  //hists_["bkg_eta_perevt"   ] = outfile_->make<TH1F>("bkg_eta_perevt"   , "bkg_eta_perevt"   , n_eta_bins-1 , eta_bins    );
  hists_["bkg_aeta_perevt"  ] = outfile_->make<TH1F>("bkg_aeta_perevt"  , "bkg_aeta_perevt"  , n_aeta_bins-1, aeta_bins   );
  hists_["bkg_phi_perevt"   ] = outfile_->make<TH1F>("bkg_phi_perevt"   , "bkg_phi_perevt"   , 10           , -pig, pig   );
  //hists_["bkg_truepu_perevt"] = outfile_->make<TH1F>("bkg_truepu_perevt", "bkg_truepu_perevt", 300          , -0.5, 299.5 ); 
  //hists_["bkg_itpu_perevt"  ] = outfile_->make<TH1F>("bkg_itpu_perevt"  , "bkg_itpu_perevt"  , 300          , -0.5, 299.5 ); 
  hists_["bkg_allpu_perevt" ] = outfile_->make<TH1F>("bkg_allpu_perevt" , "bkg_allpu_perevt" , 300          , -0.5, 299.5 ); 


  // Background rates per segment 
  // hists_["bkg_pertrack_pt_perseg"     ] = outfile_->make<TH1F>();
  // //hists_["bkg_pertrack_eta_perseg"    ] = outfile_->make<TH1F>();
  // hists_["bkg_pertrack_aeta_perseg"   ] = outfile_->make<TH1F>();
  // hists_["bkg_pertrack_phi_perseg"    ] = outfile_->make<TH1F>();
  // //hists_["bkg_pertrack_truepu_perseg" ] = outfile_->make<TH1F>(); 
  // //hists_["bkg_pertrack_itpu_perseg"   ] = outfile_->make<TH1F>(); 
  // hists_["bkg_pertrack_allpu_perseg"  ] = outfile_->make<TH1F>(); 

  hists_["bkg_pt_perseg"    ] = outfile_->make<TH1F>("bkg_pt_perseg"    , "bkg_pt_perseg"    , n_pt_bins-1  , pt_bins     );
  //hists_["bkg_eta_perseg"   ] = outfile_->make<TH1F>("bkg_eta_perseg"   , "bkg_eta_perseg"   , n_eta_bins-1 , eta_bins    );
  hists_["bkg_aeta_perseg"  ] = outfile_->make<TH1F>("bkg_aeta_perseg"  , "bkg_aeta_perseg"  , n_aeta_bins-1, aeta_bins   );
  hists_["bkg_phi_perseg"   ] = outfile_->make<TH1F>("bkg_phi_perseg"   , "bkg_phi_perseg"   , 10           , -pig, pig   );
  //hists_["bkg_truepu_perseg"] = outfile_->make<TH1F>("bkg_truepu_perseg", "bkg_truepu_perseg", 300          , -0.5, 299.5 ); 
  //hists_["bkg_itpu_perseg"  ] = outfile_->make<TH1F>("bkg_itpu_perseg"  , "bkg_itpu_perseg"  , 300          , -0.5, 299.5 ); 
  hists_["bkg_allpu_perseg" ] = outfile_->make<TH1F>("bkg_allpu_perseg" , "bkg_allpu_perseg" , 300          , -0.5, 299.5 ); 

  // hists_["fake_pertrack_pt_perevt"    ]->SetNameTitle("fake_pertrack_pt_perevt"    ,"fake_pertrack_pt_perevt"    ); 
  // //hists_["fake_pertrack_eta_perevt"   ]->SetNameTitle("fake_pertrack_eta_perevt"   ,"fake_pertrack_eta_perevt"   ); 
  // hists_["fake_pertrack_aeta_perevt"  ]->SetNameTitle("fake_pertrack_aeta_perevt"  ,"fake_pertrack_aeta_perevt"  ); 
  // hists_["fake_pertrack_phi_perevt"   ]->SetNameTitle("fake_pertrack_phi_perevt"   ,"fake_pertrack_phi_perevt"   ); 
  // //hists_["fake_pertrack_truepu_perevt"]->SetNameTitle("fake_pertrack_truepu_perevt","fake_pertrack_truepu_perevt"); 
  // //hists_["fake_pertrack_itpu_perevt"  ]->SetNameTitle("fake_pertrack_itpu_perevt"  ,"fake_pertrack_itpu_perevt"  ); 
  // hists_["fake_pertrack_allpu_perevt" ]->SetNameTitle("fake_pertrack_allpu_perevt" ,"fake_pertrack_allpu_perevt" ); 

  // hists_["fake_pt_perevt"    ]->SetNameTitle("fake_pt_perevt"    ,"fake_pt_perevt"    ); 
  // //hists_["fake_eta_perevt"   ]->SetNameTitle("fake_eta_perevt"   ,"fake_eta_perevt"   ); 
  // hists_["fake_aeta_perevt"  ]->SetNameTitle("fake_aeta_perevt"  ,"fake_aeta_perevt"  ); 
  // hists_["fake_phi_perevt"   ]->SetNameTitle("fake_phi_perevt"   ,"fake_phi_perevt"   ); 
  // //hists_["fake_truepu_perevt"]->SetNameTitle("fake_truepu_perevt","fake_truepu_perevt"); 
  // //hists_["fake_itpu_perevt"  ]->SetNameTitle("fake_itpu_perevt"  ,"fake_itpu_perevt"  ); 
  // hists_["fake_allpu_perevt" ]->SetNameTitle("fake_allpu_perevt" ,"fake_allpu_perevt" ); 

  // hists_["fake_pertrack_pt_perseg"    ]->SetNameTitle("fake_pertrack_pt_perseg"    ,"fake_pertrack_pt_perseg"    ); 
  // //hists_["fake_pertrack_eta_perseg"   ]->SetNameTitle("fake_pertrack_eta_perseg"   ,"fake_pertrack_eta_perseg"   ); 
  // hists_["fake_pertrack_aeta_perseg"  ]->SetNameTitle("fake_pertrack_aeta_perseg"  ,"fake_pertrack_aeta_perseg"  ); 
  // hists_["fake_pertrack_phi_perseg"   ]->SetNameTitle("fake_pertrack_phi_perseg"   ,"fake_pertrack_phi_perseg"   ); 
  // //hists_["fake_pertrack_truepu_perseg"]->SetNameTitle("fake_pertrack_truepu_perseg","fake_pertrack_truepu_perseg"); 
  // //hists_["fake_pertrack_itpu_perseg"  ]->SetNameTitle("fake_pertrack_itpu_perseg"  ,"fake_pertrack_itpu_perseg"  ); 
  // hists_["fake_pertrack_allpu_perseg" ]->SetNameTitle("fake_pertrack_allpu_perseg" ,"fake_pertrack_allpu_perseg" ); 

  // hists_["fake_pt_perseg"    ]->SetNameTitle("fake_pt_perseg"    ,"fake_pt_perseg"    ); 
  // //hists_["fake_eta_perseg"   ]->SetNameTitle("fake_eta_perseg"   ,"fake_eta_perseg"   ); 
  // hists_["fake_aeta_perseg"  ]->SetNameTitle("fake_aeta_perseg"  ,"fake_aeta_perseg"  ); 
  // hists_["fake_phi_perseg"   ]->SetNameTitle("fake_phi_perseg"   ,"fake_phi_perseg"   ); 
  // //hists_["fake_truepu_perseg"]->SetNameTitle("fake_truepu_perseg","fake_truepu_perseg"); 
  // //hists_["fake_itpu_perseg"  ]->SetNameTitle("fake_itpu_perseg"  ,"fake_itpu_perseg"  ); 
  // hists_["fake_allpu_perseg" ]->SetNameTitle("fake_allpu_perseg" ,"fake_allpu_perseg" ); 

  // hists_["bkg_pertrack_pt_perevt"     ]->SetNameTitle("bkg_pertrack_pt_perevt"     ,"bkg_pertrack_pt_perevt"     ); 
  // //hists_["bkg_pertrack_eta_perevt"    ]->SetNameTitle("bkg_pertrack_eta_perevt"    ,"bkg_pertrack_eta_perevt"    ); 
  // hists_["bkg_pertrack_aeta_perevt"   ]->SetNameTitle("bkg_pertrack_aeta_perevt"   ,"bkg_pertrack_aeta_perevt"   ); 
  // hists_["bkg_pertrack_phi_perevt"    ]->SetNameTitle("bkg_pertrack_phi_perevt"    ,"bkg_pertrack_phi_perevt"    ); 
  // //hists_["bkg_pertrack_truepu_perevt" ]->SetNameTitle("bkg_pertrack_truepu_perevt" ,"bkg_pertrack_truepu_perevt" ); 
  // //hists_["bkg_pertrack_itpu_perevt"   ]->SetNameTitle("bkg_pertrack_itpu_perevt"   ,"bkg_pertrack_itpu_perevt"   ); 
  // hists_["bkg_pertrack_allpu_perevt"  ]->SetNameTitle("bkg_pertrack_allpu_perevt"  ,"bkg_pertrack_allpu_perevt"  ); 

  // hists_["bkg_pt_perevt"     ]->SetNameTitle("bkg_pt_perevt"     ,"bkg_pt_perevt"     ); 
  // //hists_["bkg_eta_perevt"    ]->SetNameTitle("bkg_eta_perevt"    ,"bkg_eta_perevt"    ); 
  // hists_["bkg_aeta_perevt"   ]->SetNameTitle("bkg_aeta_perevt"   ,"bkg_aeta_perevt"   ); 
  // hists_["bkg_phi_perevt"    ]->SetNameTitle("bkg_phi_perevt"    ,"bkg_phi_perevt"    ); 
  // //hists_["bkg_truepu_perevt" ]->SetNameTitle("bkg_truepu_perevt" ,"bkg_truepu_perevt" ); 
  // //hists_["bkg_itpu_perevt"   ]->SetNameTitle("bkg_itpu_perevt"   ,"bkg_itpu_perevt"   ); 
  // hists_["bkg_allpu_perevt"  ]->SetNameTitle("bkg_allpu_perevt"  ,"bkg_allpu_perevt"  ); 

  // hists_["bkg_pertrack_pt_perseg"     ]->SetNameTitle("bkg_pertrack_pt_perseg"     ,"bkg_pertrack_pt_perseg"     ); 
  // //hists_["bkg_pertrack_eta_perseg"    ]->SetNameTitle("bkg_pertrack_eta_perseg"    ,"bkg_pertrack_eta_perseg"    ); 
  // hists_["bkg_pertrack_aeta_perseg"   ]->SetNameTitle("bkg_pertrack_aeta_perseg"   ,"bkg_pertrack_aeta_perseg"   ); 
  // hists_["bkg_pertrack_phi_perseg"    ]->SetNameTitle("bkg_pertrack_phi_perseg"    ,"bkg_pertrack_phi_perseg"    ); 
  // //hists_["bkg_pertrack_truepu_perseg" ]->SetNameTitle("bkg_pertrack_truepu_perseg" ,"bkg_pertrack_truepu_perseg" ); 
  // //hists_["bkg_pertrack_itpu_perseg"   ]->SetNameTitle("bkg_pertrack_itpu_perseg"   ,"bkg_pertrack_itpu_perseg"   ); 
  // hists_["bkg_pertrack_allpu_perseg"  ]->SetNameTitle("bkg_pertrack_allpu_perseg"  ,"bkg_pertrack_allpu_perseg"  ); 

  // hists_["bkg_pt_perseg"     ]->SetNameTitle("bkg_pt_perseg"     ,"bkg_pt_perseg"     ); 
  // //hists_["bkg_eta_perseg"    ]->SetNameTitle("bkg_eta_perseg"    ,"bkg_eta_perseg"    ); 
  // hists_["bkg_aeta_perseg"   ]->SetNameTitle("bkg_aeta_perseg"   ,"bkg_aeta_perseg"   ); 
  // hists_["bkg_phi_perseg"    ]->SetNameTitle("bkg_phi_perseg"    ,"bkg_phi_perseg"    ); 
  // //hists_["bkg_truepu_perseg" ]->SetNameTitle("bkg_truepu_perseg" ,"bkg_truepu_perseg" ); 
  // //hists_["bkg_itpu_perseg"   ]->SetNameTitle("bkg_itpu_perseg"   ,"bkg_itpu_perseg"   ); 
  // hists_["bkg_allpu_perseg"  ]->SetNameTitle("bkg_allpu_perseg"  ,"bkg_allpu_perseg"  ); 


  // Number of muons per segment 
  hists_["nmu_perseg_den_aeta"] = outfile_->make<TH1F>("nmu_perseg_den_aeta", "nmu_perseg_den_aeta", n_aeta_bins-1, aeta_bins);
  hists_["nmu_perseg_den_dphi"] = outfile_->make<TH1F>("nmu_perseg_den_dphi", "nmu_perseg_den_dphi", 40, -2.0, 2.0);

  hists_["nmu_perseg_num_aeta"] = outfile_->make<TH1F>("nmu_perseg_num_aeta", "nmu_perseg_num_aeta", n_aeta_bins-1, aeta_bins);
  hists_["nmu_perseg_num_dphi"] = outfile_->make<TH1F>("nmu_perseg_num_dphi", "nmu_perseg_num_dphi", 40, -2.0, 2.0);

  hists_["nmu_perseg_aeta"    ] = outfile_->make<TH1F>("nmu_perseg_aeta"    , "nmu_perseg_aeta"    , n_aeta_bins-1, aeta_bins);
  hists_["nmu_perseg_dphi"    ] = outfile_->make<TH1F>("nmu_perseg_dphi"    , "nmu_perseg_dphi"    , 40, -2.0, 2.0);



  // Event and segment counter 
  hists_["evt_seg_counter"] = outfile_->make<TH1F>("evt_seg_counter", "evt_seg_counter", 2, 0.5, 2.5); 
  hists_["evt_seg_counter"]->GetXaxis()->SetBinLabel(1, "N. events"); 
  hists_["evt_seg_counter"]->GetXaxis()->SetBinLabel(2, "N. segments"); 

  // Sumw2 
  hists_["eff_den_pt"      ]->Sumw2(); 
  // hists_["eff_den_eta"     ]->Sumw2(); 
  hists_["eff_den_aeta"    ]->Sumw2(); 
  hists_["eff_den_phi"     ]->Sumw2(); 
  //hists_["eff_den_truepu"  ]->Sumw2(); 
  // hists_["eff_den_itpu"    ]->Sumw2(); 
  hists_["eff_den_allpu"   ]->Sumw2(); 

  hists_["eff_num_pt"      ]->Sumw2(); 
  // hists_["eff_num_eta"     ]->Sumw2(); 
  hists_["eff_num_aeta"    ]->Sumw2(); 
  hists_["eff_num_phi"     ]->Sumw2(); 
  //hists_["eff_num_truepu"  ]->Sumw2(); 
  // hists_["eff_num_itpu"    ]->Sumw2(); 
  hists_["eff_num_allpu"   ]->Sumw2(); 

  // hists_["eff_pt"          ]->Sumw2(); 
  // // hists_["eff_eta"         ]->Sumw2(); 
  // hists_["eff_aeta"        ]->Sumw2(); 
  // hists_["eff_phi"         ]->Sumw2(); 
  // //hists_["eff_truepu"      ]->Sumw2(); 
  // // hists_["eff_itpu"        ]->Sumw2(); 
  // hists_["eff_allpu"       ]->Sumw2(); 


  // hists_["fake_pertrack_den_pt"     ]->Sumw2(); 
  // // hists_["fake_pertrack_den_eta"    ]->Sumw2(); 
  // hists_["fake_pertrack_den_aeta"   ]->Sumw2(); 
  // hists_["fake_pertrack_den_phi"    ]->Sumw2(); 
  // //hists_["fake_pertrack_den_truepu" ]->Sumw2(); 
  // // hists_["fake_pertrack_den_itpu"   ]->Sumw2(); 
  // hists_["fake_pertrack_den_allpu"  ]->Sumw2(); 

  hists_["fake_den_pt"     ]->Sumw2(); 
  // hists_["fake_den_eta"    ]->Sumw2(); 
  hists_["fake_den_aeta"   ]->Sumw2(); 
  hists_["fake_den_phi"    ]->Sumw2(); 
  //hists_["fake_den_truepu" ]->Sumw2(); 
  // hists_["fake_den_itpu"   ]->Sumw2(); 
  hists_["fake_den_allpu"  ]->Sumw2(); 

  // hists_["fake_pertrack_num_pt"     ]->Sumw2(); 
  // // hists_["fake_pertrack_num_eta"    ]->Sumw2(); 
  // hists_["fake_pertrack_num_aeta"   ]->Sumw2(); 
  // hists_["fake_pertrack_num_phi"    ]->Sumw2(); 
  // //hists_["fake_pertrack_num_truepu" ]->Sumw2(); 
  // // hists_["fake_pertrack_num_itpu"   ]->Sumw2(); 
  // hists_["fake_pertrack_num_allpu"  ]->Sumw2(); 

  // hists_["fake_pt"         ]->Sumw2(); 
  hists_["fake_num_pt"     ]->Sumw2(); 
  // hists_["fake_num_eta"    ]->Sumw2(); 
  hists_["fake_num_aeta"   ]->Sumw2(); 
  hists_["fake_num_phi"    ]->Sumw2(); 
  //hists_["fake_num_truepu" ]->Sumw2(); 
  // hists_["fake_num_itpu"   ]->Sumw2(); 
  hists_["fake_num_allpu"  ]->Sumw2(); 

  // hists_["fake_pt"         ]->Sumw2(); 
  // // hists_["fake_eta"        ]->Sumw2(); 
  // hists_["fake_aeta"       ]->Sumw2(); 
  // hists_["fake_phi"        ]->Sumw2(); 
  // //hists_["fake_truepu"     ]->Sumw2(); 
  // // hists_["fake_itpu"       ]->Sumw2(); 
  // hists_["fake_allpu"      ]->Sumw2(); 


  // hists_["bkg_den_pt"     ]->Sumw2(); 
  // // hists_["bkg_den_eta"    ]->Sumw2(); 
  // hists_["bkg_den_aeta"   ]->Sumw2(); 
  // hists_["bkg_den_phi"    ]->Sumw2(); 
  // //hists_["bkg_den_truepu" ]->Sumw2(); 
  // // hists_["bkg_den_itpu"   ]->Sumw2(); 
  // hists_["bkg_den_allpu"  ]->Sumw2(); 

  // hists_["bkg_pertrack_num_pt"     ]->Sumw2(); 
  // // hists_["bkg_pertrack_num_eta"    ]->Sumw2(); 
  // hists_["bkg_pertrack_num_aeta"   ]->Sumw2(); 
  // hists_["bkg_pertrack_num_phi"    ]->Sumw2(); 
  // //hists_["bkg_pertrack_num_truepu" ]->Sumw2(); 
  // // hists_["bkg_pertrack_num_itpu"   ]->Sumw2(); 
  // hists_["bkg_pertrack_num_allpu"  ]->Sumw2(); 

  hists_["bkg_num_pt"     ]->Sumw2(); 
  // hists_["bkg_num_eta"    ]->Sumw2(); 
  hists_["bkg_num_aeta"   ]->Sumw2(); 
  hists_["bkg_num_phi"    ]->Sumw2(); 
  //hists_["bkg_num_truepu" ]->Sumw2(); 
  // hists_["bkg_num_itpu"   ]->Sumw2(); 
  hists_["bkg_num_allpu"  ]->Sumw2(); 

  // hists_["bkg_pt"         ]->Sumw2(); 
  // // hists_["bkg_eta"        ]->Sumw2(); 
  // hists_["bkg_aeta"       ]->Sumw2(); 
  // hists_["bkg_phi"        ]->Sumw2(); 
  // //hists_["bkg_truepu"     ]->Sumw2(); 
  // // hists_["bkg_itpu"       ]->Sumw2(); 
  // hists_["bkg_allpu"      ]->Sumw2(); 

  hists_["nmu_perseg_den_aeta"]->Sumw2(); 
  hists_["nmu_perseg_den_dphi"]->Sumw2(); 

  hists_["nmu_perseg_num_aeta"]->Sumw2(); 
  hists_["nmu_perseg_num_dphi"]->Sumw2(); 

  hists_["nmu_perseg_aeta"    ]->Sumw2(); 
  hists_["nmu_perseg_dphi"    ]->Sumw2(); 
}


NewME0MuonGeneralAnalyzer::~NewME0MuonGeneralAnalyzer() { 
  if(debugPrint) std::cout << "Inside Destructor" << std::endl;
}


//
// member functions
//

// ------------ method called for each event  ------------
void NewME0MuonGeneralAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using std::string;

  if(debugPrint || outputPrint) std::cout << "============================  Inside analyze "
					  << iEvent.id().run() << "  "
					  << iEvent.id().luminosityBlock() << "  "
					  << iEvent.id().event() << "  "
					  << " ============================" << std::endl;

  hists_["evt_seg_counter"]->Fill(1.); 

  /*
  edm::Handle<reco::GenParticleCollection> genParts; 
  iEvent.getByLabel(InputTag("genParticles"), genParts); 
  if(!genParts.isValid()) {
    if(debugPrint) std::cout << "GenParticle collection not valid!" << std::endl;
    return;
  } 
  else {
    if(debugPrint) std::cout << "Found GenParticle collection!" << std::endl;
  }
  */


  edm::Handle<std::vector<PileupSummaryInfo> > puInfoH;
  //event.getByLabel("addPileupInfo", puInfoH);
  //event.getByLabel("slimmedAddPileupInfo", puInfoH);
  iEvent.getByToken(puInfoTag_,puInfoH);
  int npuOOT(0),npuIT(0),npuOOTm1(0);
  //float truePU(0);
  if(puInfoH.isValid()) {
    for(std::vector<PileupSummaryInfo>::const_iterator it=puInfoH->begin(); it!=puInfoH->end(); ++it) {
      if(it->getBunchCrossing()==0) {
	npuIT += it->getPU_NumInteractions();
	//truePU = it->getTrueNumInteractions();
      }
      else {
	//npuOOT += it->getPU_NumInteractions();
	npuOOT = 0; // TMP: figure out how to use this!!! 
      } 

      if(it->getBunchCrossing()<0)
	npuOOTm1 += it->getPU_NumInteractions();
    }
  }

  // // Tracking particles (preselected) 
  // edm::Handle<TrackingParticleCollection> trackParts; // trackParticles; 
  // iEvent.getByToken(tpcToken_, trackParts); 
  // if(!trackParts.isValid()) {
  //   //if(outputPrint || debugPrint) std::cout << "TrackingParticle collection not valid!" << std::endl;
  //   throw cms::Exception("TrackingParticle collection not valid!"); 
  // } 
  // else {
  //   if(debugPrint) std::cout << "Found TrackingParticle collection!" << std::endl;
  // }

  // GenParticles 
  edm::Handle<reco::GenParticleCollection> genParts; 
  iEvent.getByToken(gpToken_, genParts); 
  if(!genParts.isValid()) {
    throw cms::Exception("GenParticleCollection not valid!"); 
  } 
  else {
    if(debugPrint) std::cout << "Found GenParticleCollection!" << std::endl;
  }

  edm::Handle<ME0MuonCollection> muons; 
  iEvent.getByToken(muonToken_, muons); 
  if(!muons.isValid()) {
    //if(outputPrint || debugPrint) std::cout << "ME0Muon collection not valid!" << std::endl;
    throw cms::Exception("ME0Muon collection not valid!"); 
  } 
  else {
    if(debugPrint) std::cout << "Found ME0Muon collection!" << std::endl;
  }

  // Muons (all) 
  // edm::Handle<reco::MuonCollection> muons; 
  // iEvent.getByToken(muonToken_, muons); 
  // if(!muons.isValid()) {
  //   if(outputPrint || debugPrint) std::cout << "Muon collection not valid!" << std::endl;
  //   return;
  // } 
  // else {
  //   if(debugPrint) std::cout << "Found ME0Muon collection!" << std::endl;
  // }

  // // Tracks of ME0 muons (preselected) 
  // edm::Handle<edm::View<reco::Track> > trackColl; // trackCollection; 
  // iEvent.getByToken(tracksToken_, trackColl); 
  // if(!trackColl.isValid()) {
  //   //if(outputPrint || debugPrint) std::cout << "Track collection not valid!" << std::endl;
  //   throw cms::Exception("Track collection not valid!"); 
  // } 
  // else {
  //   if(debugPrint) std::cout << "Found Track collection!" << std::endl;
  // }

  edm::ESHandle<ME0Geometry> hGeom; 
  iSetup.get<MuonGeometryRecord>().get(hGeom); 

  // unsigned int nGoodVtx = 0; 
  // for(reco::VertexCollection::const_iterator vertex=vertexes->begin(); vertex!=vertexes->end(); ++vertex) {
  //   if( vertex->ndof()>4                     && 
  // 	(fabs(vertex->z())<=24.)             && 
  // 	(fabs(vertex->position().rho())<=2.)   ) 
  //     nGoodVtx++;
  // }
  // if( nGoodVtx==0 ) return;
  // const reco::Vertex & pv = (*vertexes)[0];

  edm::Handle<ME0SegmentCollection> me0Segments; 
  iEvent.getByToken(me0SegmentsToken_, me0Segments); 
  if(!me0Segments.isValid()) {
    //if(outputPrint || debugPrint) std::cout << "ME0Segment collection not valid!" << std::endl;
    throw cms::Exception("ME0Segment collection not valid!"); 
  } 
  else {
    if(debugPrint) std::cout << "Found ME0Segment collection!" << std::endl;
  }

  // edm::Handle<ME0DigiPreRecoCollection> me0digis; 
  // iEvent.getByToken(me0digiToken_, me0digis);

  // std::map<const edm::RefToBase<reco::Track>, const edm::RefToBase<ME0Muon> > tpartTrackMap; 
  std::map<int, std::pair<float, float> > segmetaphi; 
  int segmcnt(0); 
  for(ME0SegmentCollection::const_iterator segm=me0Segments->begin(); segm!=me0Segments->end(); ++segm, ++segmcnt) { 
    hists_["evt_seg_counter"]->Fill(2.); 

    ME0DetId id = (*segm).me0DetId();
    //auto roll = hGeom->etaPartition(id); 
    auto roll = hGeom->chamber(id); 
    float etapos = (roll->toGlobal((*segm).localPosition())).eta(); 
    float phipos = (roll->toGlobal((*segm).localPosition())).phi(); 
    float phidir = (roll->toGlobal((*segm).localDirection())).phi(); 

    hists_["nmu_perseg_den_aeta"]->Fill(etapos); 
    hists_["nmu_perseg_den_dphi"]->Fill(deltaPhi(phidir, phipos)); 

    segmetaphi[segmcnt] = std::pair<float, float>(etapos, deltaPhi(phidir, phipos)); 
  } 


  //bool isselected = false; 
  //unsigned int muonidx = 0; 
  ME0MuonCollection selMuons; 
  for(ME0MuonCollection::const_iterator muon=muons->begin(); muon!=muons->end(); ++muon) { 
    //isselected = isGoodME0Muon(hGeom, *muon, maxPullX, maxDiffX, maxPullX, maxDiffX, maxDiffPhiDir); 
    if( isGoodME0Muon(hGeom, *muon, maxPullX, maxDiffX, maxPullX, maxDiffX, maxDiffPhiDir) ) { 
      selMuons.push_back(*muon); 
      hists_["nmu_perseg_num_aeta"]->Fill( segmetaphi[muon->me0segid()].first  ); 
      hists_["nmu_perseg_num_dphi"]->Fill( segmetaphi[muon->me0segid()].second );       
      //++muonidx; 
    } // end if( isGoodME0Muon(hGeom, *muon, maxPullX, maxDiffX, maxPullX, maxDiffX, maxDiffPhiDir) )
  } // end for(ME0MuonCollection::const_iterator muon=muons->begin(); muon!=muons->end(); ++muon, muonidx++) 


  // reco::RecoToSimCollection recSimColl; 
  // reco::SimToRecoCollection simRecColl; 
  //TrackingParticleRefVector trackParts; // a.k.a. edm::RefVector<TrackingParticleCollection> 
  //edm::RefToBaseVector<reco::Track> trackColl; 
  //std::vector<bool> passIds; 

  // for(size_t h=0; h<trackParticles->size(); ++h) { 
  //   if((*trackParticles)[h].status()!=1) continue; 
  //   if(selectMuonOnlyTP && abs((*trackParticles)[h].pdgId())!=13) continue; 
  //   trackParts.push_back(TrackingParticleRef(trackParticles, h)); 
  // } 

  // for(size_t j=0; j<trackColl->size(); ++j) { 
  //   //edm::ProductID tkid = trackColl->refAt(j).id(); 
  //   //const reco::Track *tkptr = trackColl->refAt(j).get(); 
  //   const edm::RefToBase<reco::Track> tkptr = trackColl->refAt(j); 

  //   bool isselected = false; 
  //   unsigned int muonidx = 0; 
  //   for(ME0MuonCollection::const_iterator muon=muons->begin(); muon!=muons->end(); ++muon, muonidx++) { 

  //     //const reco::Track *muptr = muon->innerTrack().get(); 
  //     const edm::RefToBase<reco::Track> muptr = muon->innerTrack(); 

  //     if(false && debugPrint) {
  // 	std::cout << "*** Track vs Muon: " << std::endl; 
  // 	std::cout << "    Muon " << muonidx << ": (" << muptr->px() << ", " << muptr->py() << ", " << muptr->pz() << ") " 
  // 		  << " - Track " << j       << ": (" << tkptr->px() << ", " << tkptr->py() << ", " << tkptr->pz() << ") " << std::endl; 
  // 	//std::cout << " Muon: " << muptr << " - Track " << j << ": " << tkptr << std::endl; 
  // 	//std::cout << " Muon: " << muon->innerTrack().id() << " - Track " << j << ": " << trackColl->refAt(j).id() << std::endl; 
  //     } 
  //     //if(muon->innerTrack() == trackColl->refAt(j)) { // this DOESN'T work! Different types of Ref! 
  //     isselected = 
  // 	(muptr->charge()==tkptr->charge())   && 
  // 	(fabs(muptr->px()-tkptr->px())<1E-4) && 
  // 	(fabs(muptr->py()-tkptr->py())<1E-4) && 
  // 	(fabs(muptr->pz()-tkptr->pz())<1E-4) && 
  // 	(fabs(muptr->vx()-tkptr->vx())<1E-4) && 
  // 	(fabs(muptr->vy()-tkptr->vy())<1E-4) && 
  // 	(fabs(muptr->vz()-tkptr->vz())<1E-4)  ; 
  //     if(isselected) { 
  // 	// tpartTrackMap[tkptr] = muons->refAt(muonidx); 
  // 	break; // with this break, muonidx does not increment 
  //     } 

  //     // if(muptr==tkptr) { 
  //     // 	isselected = true; 
  //     // 	break; // with this break, muonidx does not increment 
  //     // } 
  //   }

  //   if(!isselected) {  
  //     //if(outputPrint || debugPrint) std::cout << "ME0 pre-selected track " << j << " can't be matched to any reco::ME0Muon!" << std::endl;
  //     char buff[99]; 
  //     sprintf(buff, "ME0 pre-selected track %d can't be matched to any reco::ME0Muon!", int(j)); 
  //     throw cms::Exception(buff); 
  //   } 
  //   else { 
  //     assert(muonidx<(*muons).size()); 
  //     // muonidx is correct, thanks to the 'break' above 
  //     hists_["nmu_perseg_num_aeta"]->Fill( segmetaphi[(*muons)[muonidx].me0segid()].first  ); 
  //     hists_["nmu_perseg_num_dphi"]->Fill( segmetaphi[(*muons)[muonidx].me0segid()].second ); 
  //   } 
  // } // end for(size_t j=0; j<trackColl->size(); ++j) 

  // Invert order of loops -- slower... 
  // for(size_t j=0; j<trackCollection->size(); ++j) { 
  //   //edm::ProductID tkid = trackCollection->refAt(j).id(); 
  //   const reco::Track *tkptr = trackCollection->refAt(j).get(); 
  //   for(ME0MuonCollection::const_iterator muon=muons->begin(); muon!=muons->end(); ++muon) { 
  //     if(debugPrint) {
  // 	std::cout << "*** Track vs ME0Muon: " << std::endl; 
  // 	std::cout << j << " : " << tkptr << " - " << muon->innerTrack().get() << std::endl; 
  // 	//std::cout << j << " : " << tkid  << " - " << muon->innerTrack().id()  << std::endl; 
  //     } 
  //     //if(muon->innerTrack() == trackCollection->refAt(j)) { // this DOESN'T work! Different types of Ref! 
  //     if(muon->innerTrack().get()==tkptr) { 
  // 	trackColl.push_back(trackCollection->refAt(j)); 
  // 	// Muon ID 
  // 	passIds.push_back( isGoodME0Muon(hGeom, *muon, maxPullX, maxDiffX, maxPullX, maxDiffX, maxDiffPhiDir) ); 
  // 	break; 
  //     } 
  //   } 
  // } 

  // if((*trackParts).size()==0 || (*trackColl).size()==0) { 
  //   recSimColl.post_insert(); 
  //   simRecColl.post_insert(); 
  // } 
  // else { 
  //   edm::Handle<reco::SimToRecoCollection> simtorecoCollectionH;
  //   iEvent.getByLabel(associatorNames[0], simtorecoCollectionH);
  //   simRecColl = *(simtorecoCollectionH.product());

  //   edm::Handle<reco::RecoToSimCollection> recotosimCollectionH;
  //   iEvent.getByLabel(associatorNames[0], recotosimCollectionH);
  //   recSimColl = *(recotosimCollectionH.product());

  //   // recSimColl = associators[0]->associateRecoToSim(trackColl, 
  //   // 						    trackParts, 
  //   // 						    &iEvent, 
  //   // 						    &iSetup); 

  //   // simRecColl = associators[0]->associateSimToReco(trackColl, 
  //   // 						    trackParts, 
  //   // 						    &iEvent, 
  //   // 						    &iSetup); 
  // } 


  // SimTracks 
  edm::Handle<SimTrackContainer> simTracksH; 

  // SimHits and SimTracks 
  //SimHitsMap pxbhimap;  
  //SimHitsMap pxblomap; 
  SimHitsMap pxehimap; 
  SimHitsMap pxelomap; 
  SimHitsMap me0map; 
  getSimHitsAndTracks(iEvent, simTracksH, /*pxbhimap, pxblomap*/ pxehimap, pxelomap, me0map); 

  // SIM-RECO association 
  associateRecoToSimByHits(/*pxbhimap, pxblomap,*/ pxehimap, pxelomap, me0map, selMuons, genParts, simTracksH); 

  // Efficiencies: Sim-Reco 
  //for(SimRecoMap::const_iterator simRecIt=simToRecoMap.begin(); simRecIt!=simToRecoMap.end(); ++simRecIt) { 
  for(SimTrackContainer::const_iterator stit=simTracksH->begin(); stit!=simTracksH->end(); ++stit) { 

    // ONLY signal muons! 
    // 
    // 1. has a GenParticle 
    int gpindex = stit->genpartIndex(); 
    if(gpindex<0) continue; 

    // 2. is a muon 
    if(abs(stit->type())!=13) continue; 

    assert(genpartIndexMap.count(gpindex)>0); 
    assert(genpartIndexMap[gpindex]<genParts->size()); 
    const reco::GenParticle& gpr = genParts->at(genpartIndexMap[gpindex]); 

    // 3. is a prompt, stable particle  
    if(gpr.isPromptFinalState()==false) continue; 

    // 4. is roughly in the ME0 |eta| region 
    const FourMomentum thisTp = stit->momentum(); // just recycling the name from TrackingParticle... 
    if(fabs(thisTp.eta())<1.4 || fabs(thisTp.eta())>3.4) continue; 

    //SimTrack& thisSt = (*simTracksH)[simRecIt->first]; 
    //SimTrackRef thisSt = simTracksH.refAt(simRecIt->first); 

    // Fill plots 
    if(thisTp.P()>minP || thisTp.Pt()>minPt) { 
      // hists_["eff_den_eta" ]->Fill( thisTp.eta() );
      hists_["eff_den_aeta"]->Fill( fabs(thisTp.eta()) );
    } 
    if(fabs(thisTp.eta())>2.0 && fabs(thisTp.eta())<2.8) { // only within ME0 |eta| fiducial region 
      hists_["eff_den_pt" ]->Fill( thisTp.Pt()  ); 
      if(thisTp.P()>minP || thisTp.Pt()>minPt) { 
	hists_["eff_den_phi"]->Fill( thisTp.Phi() ); 
	//hists_["eff_den_truepu"]->Fill( truePU ); 
	// hists_["eff_den_itpu"]->Fill( npuIT ); 
	hists_["eff_den_allpu"]->Fill( npuIT+npuOOT ); 
      } 
    } 
		       
    bool matchFound(simToRecoMap.count(stit->trackId())>0); 

    // Match found! 
    if(matchFound) { 
      // Fill plots 
      if(thisTp.P()>minP || thisTp.Pt()>minPt) { 
	// hists_["eff_num_eta" ]->Fill( thisTp.eta() ); 
	hists_["eff_num_aeta"]->Fill( fabs(thisTp.eta()) ); 
      } 
      if(fabs(thisTp.eta())>2.0 && fabs(thisTp.eta())<2.8) { 
	hists_["eff_num_pt" ]->Fill( thisTp.Pt()  ); 
	if(thisTp.P()>minP || thisTp.Pt()>minPt) { 
	  hists_["eff_num_phi"]->Fill( thisTp.Phi() ); 
	  //hists_["eff_num_truepu"]->Fill( truePU ); 
	  // hists_["eff_num_itpu"]->Fill( npuIT ); 
	  hists_["eff_num_allpu"]->Fill( npuIT+npuOOT ); 
	} 
      } 

      // // Segment composition 
      // SimTrack &stk = thisTp.g4Tracks()[0]; 
      // const edm::RefToBase<ME0Muon> me0muref = tpartTrackMap[assoTrack]; 
      // auto me0segm = me0muref->me0segment(); 
      // auto me0rhs = me0segm.specificRecHits(); 
      // for (auto rh=me0rhs.begin(); rh!=me0rhs.end(); ++rh) {
      // 	ME0DetId me0id = rh->me0Id(); 
      // 	int region = me0id.region(); 
      // 	int chamber = me0id.chamber(); 
      // 	int roll = me0id.roll(); 
      // 	int layer = me0id.layer(); 
      // 	me0digis; 
      // } 
    } // end if(matchFound) 
  } // end for(SimTrackContainer::const_iterator stit=simTracksH->begin(); stit!=simTracksH->end(); ++stit)

  // Fakes: Reco-Sim 
  for(unsigned int h=0; h<selMuons.size(); ++h) { 
    reco::TrackRef thisTrack = selMuons[h].innerTrack(); 

    if(debugPrint) printf("--- %u) Muon track -- (%.3f, %.3f, %.3f)\n", h, thisTrack->p(), thisTrack->eta(), thisTrack->phi()); 

    if(fabs(thisTrack->eta())>2.0 && fabs(thisTrack->eta())<2.8) { 
      hists_["fake_den_pt"    ]->Fill( thisTrack->pt()        ); 
      hists_["fake_den_aeta"  ]->Fill( fabs(thisTrack->eta()) ); 
      hists_["fake_den_phi"   ]->Fill( thisTrack->phi()       ); 
      //hists_["fake_den_truepu"]->Fill( truePU                 ); 
      // hists_["fake_den_itpu"  ]->Fill( npuIT                  ); 
      hists_["fake_den_allpu" ]->Fill( npuIT+npuOOT           ); 
    }

    bool matchFound(recoToSimMap.count(h)>0); 
    bool matchOnlySignalFound(false); 
    if(matchFound) { 
      assert( trkIdIndexMap[recoToSimMap[h]]<simTracksH->size() ); 
      const SimTrack& stref = simTracksH->at( trkIdIndexMap[recoToSimMap[h]] ); 
      int gpindex = stref.genpartIndex(); 
      if(gpindex>=0) { 
	assert(genpartIndexMap.count(gpindex)>0); 
	assert(genpartIndexMap[gpindex]<genParts->size()); 
	const reco::GenParticle& gpr = genParts->at(genpartIndexMap[gpindex]); 
	matchOnlySignalFound = gpr.isPromptFinalState(); 
      } 
    } 

    // NO match! 
    if(!matchFound) { 
      // // Fill plots SINGLE TRACK (i.e. NOT counting ghosts from single track matched to multiple segments)  
      // // hists_["fake_pertrack_num_eta"   ]->Fill( thisTrack->eta()       ); 
      // // hists_["fake_pertrack_eta_perevt"]->Fill( thisTrack->eta()       ); 
      // // hists_["fake_pertrack_eta_perseg"]->Fill( thisTrack->eta()       ); 

      // if(fabs(thisTrack->eta())>2.0 && fabs(thisTrack->eta())<2.8) { 
      // 	hists_["fake_pertrack_num_pt"      ]->Fill( thisTrack->pt()        ); 
      // 	hists_["fake_pertrack_pt_perevt"   ]->Fill( thisTrack->pt()        ); 
      // 	hists_["fake_pertrack_pt_perseg"   ]->Fill( thisTrack->pt()        ); 

      // 	hists_["fake_pertrack_num_aeta"    ]->Fill( fabs(thisTrack->eta()) ); 
      // 	hists_["fake_pertrack_aeta_perevt" ]->Fill( fabs(thisTrack->eta()) ); 
      // 	hists_["fake_pertrack_aeta_perseg" ]->Fill( fabs(thisTrack->eta()) ); 

      // 	hists_["fake_pertrack_num_phi"     ]->Fill( thisTrack->phi()       ); 
      // 	hists_["fake_pertrack_phi_perevt"  ]->Fill( thisTrack->phi()       ); 
      // 	hists_["fake_pertrack_phi_perseg"  ]->Fill( thisTrack->phi()       ); 

      // 	// hists_["fake_pertrack_num_truepu"   ]->Fill( truePU                 ); 
      // 	// hists_["fake_pertrack_truepu_perevt"]->Fill( truePU                 ); 
      // 	// hists_["fake_pertrack_truepu_perseg"]->Fill( truePU                 ); 

      // 	// hists_["fake_pertrack_num_itpu"     ]->Fill( npuIT                  ); 
      // 	// hists_["fake_pertrack_itpu_perevt"  ]->Fill( npuIT                  ); 
      // 	// hists_["fake_pertrack_itpu_perseg"  ]->Fill( npuIT                  ); 

      // 	hists_["fake_pertrack_num_allpu"    ]->Fill( npuIT+npuOOT           );
      // 	hists_["fake_pertrack_allpu_perevt" ]->Fill( npuIT+npuOOT           );
      // 	hists_["fake_pertrack_allpu_perseg" ]->Fill( npuIT+npuOOT           );
      // }

      // Fill plots ALL TRACKS (i.e. counting ghosts from single track matched to multiple segments)  
      // hists_["fake_num_eta"   ]->Fill( thisTrack->eta()       ); 
      // hists_["fake_eta_perevt"]->Fill( thisTrack->eta()       ); 
      // hists_["fake_eta_perseg"]->Fill( thisTrack->eta()       ); 

      if(fabs(thisTrack->eta())>2.0 && fabs(thisTrack->eta())<2.8) { 
	hists_["fake_num_pt"      ]->Fill( thisTrack->pt()        ); 
	hists_["fake_pt_perevt"   ]->Fill( thisTrack->pt()        ); 
	hists_["fake_pt_perseg"   ]->Fill( thisTrack->pt()        ); 

	hists_["fake_num_aeta"    ]->Fill( fabs(thisTrack->eta()) ); 
	hists_["fake_aeta_perevt" ]->Fill( fabs(thisTrack->eta()) ); 
	hists_["fake_aeta_perseg" ]->Fill( fabs(thisTrack->eta()) ); 

	hists_["fake_num_phi"     ]->Fill( thisTrack->phi()       ); 
	hists_["fake_phi_perevt"  ]->Fill( thisTrack->phi()       ); 
	hists_["fake_phi_perseg"  ]->Fill( thisTrack->phi()       ); 

	// hists_["fake_num_truepu"   ]->Fill( truePU                 ); 
	// hists_["fake_truepu_perevt"]->Fill( truePU                 ); 
	// hists_["fake_truepu_perseg"]->Fill( truePU                 ); 

	// hists_["fake_num_itpu"     ]->Fill( npuIT                  ); 
	// hists_["fake_itpu_perevt"  ]->Fill( npuIT                  ); 
	// hists_["fake_itpu_perseg"  ]->Fill( npuIT                  ); 

	hists_["fake_num_allpu"    ]->Fill( npuIT+npuOOT           );
	hists_["fake_allpu_perevt" ]->Fill( npuIT+npuOOT           );
	hists_["fake_allpu_perseg" ]->Fill( npuIT+npuOOT           );
      }
    } // end if(!matchFound) 

    // NO match with a SIGNAL particle (e.g. Z->mumu)! 
    // (But still possible matches with a background particle, e.g. pi->munu) 
    if(!matchOnlySignalFound) { 
      // // Fill plots 
      // // hists_["bkg_pertrack_num_eta"      ]->Fill( thisTrack->eta()       ); 
      // // hists_["bkg_pertrack_eta_perevt"   ]->Fill( thisTrack->eta()       ); 
      // // hists_["bkg_pertrack_eta_perseg"   ]->Fill( thisTrack->eta()       ); 

      // if(fabs(thisTrack->eta())>2.0 && fabs(thisTrack->eta())<2.8) { 
      // 	hists_["bkg_pertrack_num_pt"       ]->Fill( thisTrack->pt()        ); 
      // 	hists_["bkg_pertrack_pt_perevt"    ]->Fill( thisTrack->pt()        ); 
      // 	hists_["bkg_pertrack_pt_perseg"    ]->Fill( thisTrack->pt()        ); 

      // 	hists_["bkg_pertrack_num_aeta"     ]->Fill( fabs(thisTrack->eta()) ); 
      // 	hists_["bkg_pertrack_aeta_perevt"  ]->Fill( fabs(thisTrack->eta()) ); 
      // 	hists_["bkg_pertrack_aeta_perseg"  ]->Fill( fabs(thisTrack->eta()) ); 

      // 	hists_["bkg_pertrack_num_phi"      ]->Fill( thisTrack->phi()       ); 
      // 	hists_["bkg_pertrack_phi_perevt"   ]->Fill( thisTrack->phi()       ); 
      // 	hists_["bkg_pertrack_phi_perseg"   ]->Fill( thisTrack->phi()       ); 

      // 	// hists_["bkg_pertrack_num_truepu"   ]->Fill( truePU                 ); 
      // 	// hists_["bkg_pertrack_truepu_perevt"]->Fill( truePU                 ); 
      // 	// hists_["bkg_pertrack_truepu_perseg"]->Fill( truePU                 ); 

      // 	// hists_["bkg_pertrack_num_itpu"     ]->Fill( npuIT                  ); 
      // 	// hists_["bkg_pertrack_itpu_perevt"  ]->Fill( npuIT                  ); 
      // 	// hists_["bkg_pertrack_itpu_perseg"  ]->Fill( npuIT                  ); 

      // 	hists_["bkg_pertrack_num_allpu"    ]->Fill( npuIT+npuOOT           );
      // 	hists_["bkg_pertrack_allpu_perevt" ]->Fill( npuIT+npuOOT           );
      // 	hists_["bkg_pertrack_allpu_perseg" ]->Fill( npuIT+npuOOT           );
      // }

      // Fill plots 
      // hists_["bkg_num_eta"      ]->Fill( thisTrack->eta()       ); 
      // hists_["bkg_eta_perevt"   ]->Fill( thisTrack->eta()       ); 
      // hists_["bkg_eta_perseg"   ]->Fill( thisTrack->eta()       ); 

      if(fabs(thisTrack->eta())>2.0 && fabs(thisTrack->eta())<2.8) { 
	hists_["bkg_num_pt"       ]->Fill( thisTrack->pt()        ); 
	hists_["bkg_pt_perevt"    ]->Fill( thisTrack->pt()        ); 
	hists_["bkg_pt_perseg"    ]->Fill( thisTrack->pt()        ); 

	hists_["bkg_num_aeta"     ]->Fill( fabs(thisTrack->eta()) ); 
	hists_["bkg_aeta_perevt"  ]->Fill( fabs(thisTrack->eta()) ); 
	hists_["bkg_aeta_perseg"  ]->Fill( fabs(thisTrack->eta()) ); 

	hists_["bkg_num_phi"      ]->Fill( thisTrack->phi()       ); 
	hists_["bkg_phi_perevt"   ]->Fill( thisTrack->phi()       ); 
	hists_["bkg_phi_perseg"   ]->Fill( thisTrack->phi()       ); 

	// hists_["bkg_num_truepu"   ]->Fill( truePU                 ); 
	// hists_["bkg_truepu_perevt"]->Fill( truePU                 ); 
	// hists_["bkg_truepu_perseg"]->Fill( truePU                 ); 

	// hists_["bkg_num_itpu"     ]->Fill( npuIT                  ); 
	// hists_["bkg_itpu_perevt"  ]->Fill( npuIT                  ); 
	// hists_["bkg_itpu_perseg"  ]->Fill( npuIT                  ); 

	hists_["bkg_num_allpu"    ]->Fill( npuIT+npuOOT           );
	hists_["bkg_allpu_perevt" ]->Fill( npuIT+npuOOT           );
	hists_["bkg_allpu_perseg" ]->Fill( npuIT+npuOOT           );
      } 
    } // end if(!matchOnlySignalFound) 
  } 

  // Clear maps 
  simToRecoMap.clear(); 
  recoToSimMap.clear(); 
  trkIdIndexMap.clear(); 
  genpartIndexMap.clear(); 

  return;
}


// ------------ method called once each job just before starting event loop  ------------
void NewME0MuonGeneralAnalyzer::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void NewME0MuonGeneralAnalyzer::endJob() { 
  // hists_["eff_pt"      ]->Divide( hists_["eff_num_pt"    ], hists_["eff_den_pt"    ], 1.0, 1.0, "B" ); 
  // // hists_["eff_eta"     ]->Divide( hists_["eff_num_eta"   ], hists_["eff_den_eta"   ], 1.0, 1.0, "B" ); 
  // hists_["eff_aeta"    ]->Divide( hists_["eff_num_aeta"  ], hists_["eff_den_aeta"  ], 1.0, 1.0, "B" ); 
  // hists_["eff_phi"     ]->Divide( hists_["eff_num_phi"   ], hists_["eff_den_phi"   ], 1.0, 1.0, "B" ); 
  // // hists_["eff_truepu"  ]->Divide( hists_["eff_num_truepu"], hists_["eff_den_truepu"], 1.0, 1.0, "B" ); 
  // // hists_["eff_itpu"    ]->Divide( hists_["eff_num_itpu"]  , hists_["eff_den_itpu"]  , 1.0, 1.0, "B" ); 
  // hists_["eff_allpu"   ]->Divide( hists_["eff_num_allpu"] , hists_["eff_den_allpu"] , 1.0, 1.0, "B" ); 

  try{ graphs_["eff_pt_err"    ]->Divide( hists_["eff_num_pt"    ], hists_["eff_den_pt"    ], "cl=0.683 b(1,1) mode" ); } catch(cms::Exception& ex) {} 
  // try{ graphs_["eff_eta_err"   ]->Divide( hists_["eff_num_eta"   ], hists_["eff_den_eta"   ], "cl=0.683 b(1,1) mode" ); } catch(cms::Exception& ex) {} 
  try{ graphs_["eff_aeta_err"  ]->Divide( hists_["eff_num_aeta"  ], hists_["eff_den_aeta"  ], "cl=0.683 b(1,1) mode" ); } catch(cms::Exception& ex) {} 
  try{ graphs_["eff_phi_err"   ]->Divide( hists_["eff_num_phi"   ], hists_["eff_den_phi"   ], "cl=0.683 b(1,1) mode" ); } catch(cms::Exception& ex) {} 
  // try{ graphs_["eff_truepu_err"]->Divide( hists_["eff_num_truepu"], hists_["eff_den_truepu"], "cl=0.683 b(1,1) mode" ); } catch(cms::Exception& ex) {} 
  // try{ graphs_["eff_itpu_err"  ]->Divide( hists_["eff_num_itpu"  ], hists_["eff_den_itpu"  ], "cl=0.683 b(1,1) mode" ); } catch(cms::Exception& ex) {} 
  try{ graphs_["eff_allpu_err" ]->Divide( hists_["eff_num_allpu" ], hists_["eff_den_allpu" ], "cl=0.683 b(1,1) mode" ); } catch(cms::Exception& ex) {} 

  // hists_["fake_pt"     ]->Divide( hists_["fake_num_pt"    ], hists_["fake_den_pt"    ], 1.0, 1.0, "B" ); 
  // // hists_["fake_eta"    ]->Divide( hists_["fake_num_eta"   ], hists_["fake_den_eta"   ], 1.0, 1.0, "B" ); 
  // hists_["fake_aeta"   ]->Divide( hists_["fake_num_aeta"  ], hists_["fake_den_aeta"  ], 1.0, 1.0, "B" ); 
  // hists_["fake_phi"    ]->Divide( hists_["fake_num_phi"   ], hists_["fake_den_phi"   ], 1.0, 1.0, "B" ); 
  // // hists_["fake_truepu" ]->Divide( hists_["fake_num_truepu"], hists_["fake_den_truepu"], 1.0, 1.0, "B" ); 
  // // hists_["fake_itpu"   ]->Divide( hists_["fake_num_itpu"  ], hists_["fake_den_itpu"  ], 1.0, 1.0, "B" ); 
  // hists_["fake_allpu"  ]->Divide( hists_["fake_num_allpu" ], hists_["fake_den_allpu" ], 1.0, 1.0, "B" ); 

  // try{ graphs_["fake_pertrack_pt_err"    ]->Divide( hists_["fake_pertrack_num_pt"    ], hists_["fake_pertrack_den_pt"    ], "cl=0.683 b(1,1) mode" ); } catch(cms::Exception& ex) {} 
  // // try{ graphs_["fake_pertrack_eta_err"   ]->Divide( hists_["fake_pertrack_num_eta"   ], hists_["fake_pertrack_den_eta"   ], "cl=0.683 b(1,1) mode" ); } catch(cms::Exception& ex) {} 
  // try{ graphs_["fake_pertrack_aeta_err"  ]->Divide( hists_["fake_pertrack_num_aeta"  ], hists_["fake_pertrack_den_aeta"  ], "cl=0.683 b(1,1) mode" ); } catch(cms::Exception& ex) {} 
  // try{ graphs_["fake_pertrack_phi_err"   ]->Divide( hists_["fake_pertrack_num_phi"   ], hists_["fake_pertrack_den_phi"   ], "cl=0.683 b(1,1) mode" ); } catch(cms::Exception& ex) {} 
  // // try{ graphs_["fake_pertrack_truepu_err"]->Divide( hists_["fake_pertrack_num_truepu"], hists_["fake_pertrack_den_truepu"], "cl=0.683 b(1,1) mode" ); } catch(cms::Exception& ex) {} 
  // // try{ graphs_["fake_pertrack_itpu_err"  ]->Divide( hists_["fake_pertrack_num_itpu"  ], hists_["fake_pertrack_den_itpu¯"  ], "cl=0.683 b(1,1) mode"); } catch(cms::Exception& ex) {} 
  // try{ graphs_["fake_pertrack_allpu_err" ]->Divide( hists_["fake_pertrack_num_allpu" ], hists_["fake_pertrack_den_allpu" ], "cl=0.683 b(1,1) mode" ); } catch(cms::Exception& ex) {} 

  try{ graphs_["fake_pt_err"    ]->Divide( hists_["fake_num_pt"    ], hists_["fake_den_pt"    ], "cl=0.683 b(1,1) mode" ); } catch(cms::Exception& ex) {} 
  // try{ graphs_["fake_eta_err"   ]->Divide( hists_["fake_num_eta"   ], hists_["fake_den_eta"   ], "cl=0.683 b(1,1) mode" ); } catch(cms::Exception& ex) {} 
  try{ graphs_["fake_aeta_err"  ]->Divide( hists_["fake_num_aeta"  ], hists_["fake_den_aeta"  ], "cl=0.683 b(1,1) mode" ); } catch(cms::Exception& ex) {} 
  try{ graphs_["fake_phi_err"   ]->Divide( hists_["fake_num_phi"   ], hists_["fake_den_phi"   ], "cl=0.683 b(1,1) mode" ); } catch(cms::Exception& ex) {} 
  // try{ graphs_["fake_truepu_err"]->Divide( hists_["fake_num_truepu"], hists_["fake_den_truepu"], "cl=0.683 b(1,1) mode" ); } catch(cms::Exception& ex) {} 
  // try{ graphs_["fake_itpu_err"  ]->Divide( hists_["fake_num_itpu"  ], hists_["fake_den_itpu¯"  ], "cl=0.683 b(1,1) mode"); } catch(cms::Exception& ex) {} 
  try{ graphs_["fake_allpu_err" ]->Divide( hists_["fake_num_allpu" ], hists_["fake_den_allpu" ], "cl=0.683 b(1,1) mode" ); } catch(cms::Exception& ex) {} 

  // try{ graphs_["bkg_pertrack_pt_err"    ]->Divide( hists_["bkg_pertrack_num_pt"    ], hists_["fake_den_pt"    ], "cl=0.683 b(1,1) mode" ); } catch(cms::Exception& ex) {} 
  // // try{ graphs_["bkg_pertrack_eta_err"   ]->Divide( hists_["bkg_pertrack_num_eta"   ], hists_["fake_den_eta"   ], "cl=0.683 b(1,1) mode" ); } catch(cms::Exception& ex) {} 
  // try{ graphs_["bkg_pertrack_aeta_err"  ]->Divide( hists_["bkg_pertrack_num_aeta"  ], hists_["fake_den_aeta"  ], "cl=0.683 b(1,1) mode" ); } catch(cms::Exception& ex) {} 
  // try{ graphs_["bkg_pertrack_phi_err"   ]->Divide( hists_["bkg_pertrack_num_phi"   ], hists_["fake_den_phi"   ], "cl=0.683 b(1,1) mode" ); } catch(cms::Exception& ex) {} 
  // // try{ graphs_["bkg_pertrack_truepu_err"]->Divide( hists_["bkg_pertrack_num_truepu"], hists_["fake_den_truepu"], "cl=0.683 b(1,1) mode" ); } catch(cms::Exception& ex) {} 
  // // try{ graphs_["bkg_pertrack_itpu_err"  ]->Divide( hists_["bkg_pertrack_num_itpu"  ], hists_["fake_den_itpu¯"  ], "cl=0.683 b(1,1) mode"); } catch(cms::Exception& ex) {} 
  // try{ graphs_["bkg_pertrack_allpu_err" ]->Divide( hists_["bkg_pertrack_num_allpu" ], hists_["fake_den_allpu" ], "cl=0.683 b(1,1) mode" ); } catch(cms::Exception& ex) {} 

  try{ graphs_["bkg_pt_err"    ]->Divide( hists_["bkg_num_pt"    ], hists_["fake_den_pt"    ], "cl=0.683 b(1,1) mode" ); } catch(cms::Exception& ex) {} 
  // try{ graphs_["bkg_eta_err"   ]->Divide( hists_["bkg_num_eta"   ], hists_["fake_den_eta"   ], "cl=0.683 b(1,1) mode" ); } catch(cms::Exception& ex) {} 
  try{ graphs_["bkg_aeta_err"  ]->Divide( hists_["bkg_num_aeta"  ], hists_["fake_den_aeta"  ], "cl=0.683 b(1,1) mode" ); } catch(cms::Exception& ex) {} 
  try{ graphs_["bkg_phi_err"   ]->Divide( hists_["bkg_num_phi"   ], hists_["fake_den_phi"   ], "cl=0.683 b(1,1) mode" ); } catch(cms::Exception& ex) {} 
  // try{ graphs_["bkg_truepu_err"]->Divide( hists_["bkg_num_truepu"], hists_["fake_den_truepu"], "cl=0.683 b(1,1) mode" ); } catch(cms::Exception& ex) {} 
  // try{ graphs_["bkg_itpu_err"  ]->Divide( hists_["bkg_num_itpu"  ], hists_["fake_den_itpu¯"  ], "cl=0.683 b(1,1) mode"); } catch(cms::Exception& ex) {} 
  try{ graphs_["bkg_allpu_err" ]->Divide( hists_["bkg_num_allpu" ], hists_["fake_den_allpu" ], "cl=0.683 b(1,1) mode" ); } catch(cms::Exception& ex) {} 


  // hists_["fake_pertrack_pt_perevt"    ] ->Scale( 1./hists_["evt_seg_counter"]->GetBinContent(1) ); 
  // //hists_["fake_pertrack_eta_perevt"   ] ->Scale( 1./hists_["evt_seg_counter"]->GetBinContent(1) ); 
  // hists_["fake_pertrack_aeta_perevt"  ] ->Scale( 1./hists_["evt_seg_counter"]->GetBinContent(1) ); 
  // hists_["fake_pertrack_phi_perevt"   ] ->Scale( 1./hists_["evt_seg_counter"]->GetBinContent(1) ); 
  // //hists_["fake_pertrack_truepu_perevt"] ->Scale( 1./hists_["evt_seg_counter"]->GetBinContent(1) ); 
  // //hists_["fake_pertrack_itpu_perevt"  ] ->Scale( 1./hists_["evt_seg_counter"]->GetBinContent(1) ); 
  // hists_["fake_pertrack_allpu_perevt" ] ->Scale( 1./hists_["evt_seg_counter"]->GetBinContent(1) ); 
  
  hists_["fake_pt_perevt"    ] ->Scale( 1./hists_["evt_seg_counter"]->GetBinContent(1) ); 
  //hists_["fake_eta_perevt"   ] ->Scale( 1./hists_["evt_seg_counter"]->GetBinContent(1) ); 
  hists_["fake_aeta_perevt"  ] ->Scale( 1./hists_["evt_seg_counter"]->GetBinContent(1) ); 
  hists_["fake_phi_perevt"   ] ->Scale( 1./hists_["evt_seg_counter"]->GetBinContent(1) ); 
  //hists_["fake_truepu_perevt"] ->Scale( 1./hists_["evt_seg_counter"]->GetBinContent(1) ); 
  //hists_["fake_itpu_perevt"  ] ->Scale( 1./hists_["evt_seg_counter"]->GetBinContent(1) ); 
  hists_["fake_allpu_perevt" ] ->Scale( 1./hists_["evt_seg_counter"]->GetBinContent(1) ); 
  
  // hists_["fake_pertrack_pt_perseg"    ] ->Scale( 1./hists_["evt_seg_counter"]->GetBinContent(2) ); 
  // //hists_["fake_pertrack_eta_perseg"   ] ->Scale( 1./hists_["evt_seg_counter"]->GetBinContent(2) ); 
  // hists_["fake_pertrack_aeta_perseg"  ] ->Scale( 1./hists_["evt_seg_counter"]->GetBinContent(2) ); 
  // hists_["fake_pertrack_phi_perseg"   ] ->Scale( 1./hists_["evt_seg_counter"]->GetBinContent(2) ); 
  // //hists_["fake_pertrack_truepu_perseg"] ->Scale( 1./hists_["evt_seg_counter"]->GetBinContent(2) ); 
  // //hists_["fake_pertrack_itpu_perseg"  ] ->Scale( 1./hists_["evt_seg_counter"]->GetBinContent(2) ); 
  // hists_["fake_pertrack_allpu_perseg" ] ->Scale( 1./hists_["evt_seg_counter"]->GetBinContent(2) ); 
  
  hists_["fake_pt_perseg"    ] ->Scale( 1./hists_["evt_seg_counter"]->GetBinContent(2) ); 
  //hists_["fake_eta_perseg"   ] ->Scale( 1./hists_["evt_seg_counter"]->GetBinContent(2) ); 
  hists_["fake_aeta_perseg"  ] ->Scale( 1./hists_["evt_seg_counter"]->GetBinContent(2) ); 
  hists_["fake_phi_perseg"   ] ->Scale( 1./hists_["evt_seg_counter"]->GetBinContent(2) ); 
  //hists_["fake_truepu_perseg"] ->Scale( 1./hists_["evt_seg_counter"]->GetBinContent(2) ); 
  //hists_["fake_itpu_perseg"  ] ->Scale( 1./hists_["evt_seg_counter"]->GetBinContent(2) ); 
  hists_["fake_allpu_perseg" ] ->Scale( 1./hists_["evt_seg_counter"]->GetBinContent(2) ); 
  
  // hists_["bkg_pertrack_pt_perevt"     ] ->Scale( 1./hists_["evt_seg_counter"]->GetBinContent(1) ); 
  // //hists_["bkg_pertrack_eta_perevt"    ] ->Scale( 1./hists_["evt_seg_counter"]->GetBinContent(1) ); 
  // hists_["bkg_pertrack_aeta_perevt"   ] ->Scale( 1./hists_["evt_seg_counter"]->GetBinContent(1) ); 
  // hists_["bkg_pertrack_phi_perevt"    ] ->Scale( 1./hists_["evt_seg_counter"]->GetBinContent(1) ); 
  // //hists_["bkg_pertrack_truepu_perevt" ] ->Scale( 1./hists_["evt_seg_counter"]->GetBinContent(1) ); 
  // //hists_["bkg_pertrack_itpu_perevt"   ] ->Scale( 1./hists_["evt_seg_counter"]->GetBinContent(1) ); 
  // hists_["bkg_pertrack_allpu_perevt"  ] ->Scale( 1./hists_["evt_seg_counter"]->GetBinContent(1) ); 
  
  hists_["bkg_pt_perevt"     ] ->Scale( 1./hists_["evt_seg_counter"]->GetBinContent(1) ); 
  //hists_["bkg_eta_perevt"    ] ->Scale( 1./hists_["evt_seg_counter"]->GetBinContent(1) ); 
  hists_["bkg_aeta_perevt"   ] ->Scale( 1./hists_["evt_seg_counter"]->GetBinContent(1) ); 
  hists_["bkg_phi_perevt"    ] ->Scale( 1./hists_["evt_seg_counter"]->GetBinContent(1) ); 
  //hists_["bkg_truepu_perevt" ] ->Scale( 1./hists_["evt_seg_counter"]->GetBinContent(1) ); 
  //hists_["bkg_itpu_perevt"   ] ->Scale( 1./hists_["evt_seg_counter"]->GetBinContent(1) ); 
  hists_["bkg_allpu_perevt"  ] ->Scale( 1./hists_["evt_seg_counter"]->GetBinContent(1) ); 
  
  // hists_["bkg_pertrack_pt_perseg"     ] ->Scale( 1./hists_["evt_seg_counter"]->GetBinContent(2) ); 
  // //hists_["bkg_pertrack_eta_perseg"    ] ->Scale( 1./hists_["evt_seg_counter"]->GetBinContent(2) ); 
  // hists_["bkg_pertrack_aeta_perseg"   ] ->Scale( 1./hists_["evt_seg_counter"]->GetBinContent(2) ); 
  // hists_["bkg_pertrack_phi_perseg"    ] ->Scale( 1./hists_["evt_seg_counter"]->GetBinContent(2) ); 
  // //hists_["bkg_pertrack_truepu_perseg" ] ->Scale( 1./hists_["evt_seg_counter"]->GetBinContent(2) ); 
  // //hists_["bkg_pertrack_itpu_perseg"   ] ->Scale( 1./hists_["evt_seg_counter"]->GetBinContent(2) ); 
  // hists_["bkg_pertrack_allpu_perseg"  ] ->Scale( 1./hists_["evt_seg_counter"]->GetBinContent(2) ); 

  hists_["bkg_pt_perseg"     ] ->Scale( 1./hists_["evt_seg_counter"]->GetBinContent(2) ); 
  //hists_["bkg_eta_perseg"    ] ->Scale( 1./hists_["evt_seg_counter"]->GetBinContent(2) ); 
  hists_["bkg_aeta_perseg"   ] ->Scale( 1./hists_["evt_seg_counter"]->GetBinContent(2) ); 
  hists_["bkg_phi_perseg"    ] ->Scale( 1./hists_["evt_seg_counter"]->GetBinContent(2) ); 
  //hists_["bkg_truepu_perseg" ] ->Scale( 1./hists_["evt_seg_counter"]->GetBinContent(2) ); 
  //hists_["bkg_itpu_perseg"   ] ->Scale( 1./hists_["evt_seg_counter"]->GetBinContent(2) ); 
  hists_["bkg_allpu_perseg"  ] ->Scale( 1./hists_["evt_seg_counter"]->GetBinContent(2) ); 

  try{ hists_["nmu_perseg_aeta"]->Divide(hists_["nmu_perseg_num_aeta"], hists_["nmu_perseg_den_aeta"]); } catch(cms::Exception& ex) {} 
  try{ hists_["nmu_perseg_dphi"]->Divide(hists_["nmu_perseg_num_dphi"], hists_["nmu_perseg_den_dphi"]); } catch(cms::Exception& ex) {} 
} 

// ------------ method called when starting to processes a run  ------------
void NewME0MuonGeneralAnalyzer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {
  // edm::ESHandle<TrackAssociatorBase> theAssociator;
  // for(size_t w=0; w<associatorNames.size(); ++w) { 
  //   iSetup.get<TrackAssociatorRecord>().get(associatorNames[w], theAssociator); 
  //   if(theAssociator.isValid()) associators.push_back(theAssociator.product()); 
  // } 
  // if(associators.size()==0) { 
  //   TString assonames; 
  //   for(size_t w=0; w<associatorNames.size(); ++w) assonames += (associatorNames[w]+" "); 
  //   throw cms::Exception("No Associators found")
  //     << "Cannot find any of the following associators: " << assonames.Data(); 
  // } 
} // end of beginRun


// ------------ method called when ending the processing of a run  ------------
void NewME0MuonGeneralAnalyzer::endRun(edm::Run const&, edm::EventSetup const&) {
}


// ------------ method called when starting to processes a luminosity block  ------------
void NewME0MuonGeneralAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

// ------------ method called when ending the processing of a luminosity block  ------------
void NewME0MuonGeneralAnalyzer::endLuminosityBlock(edm::LuminosityBlock const& iLumiBlock, edm::EventSetup const&) {
}


bool NewME0MuonGeneralAnalyzer::isGoodME0Muon(edm::ESHandle<ME0Geometry> me0Geom, const reco::ME0Muon& me0muon, double MaxPullX, double MaxDiffX, double MaxPullY, double MaxDiffY, double MaxDiffPhiDir )
{
  ME0Segment thisSegment = me0muon.me0segment();
  
  ME0DetId id = thisSegment.me0DetId();

  //auto roll = me0Geom->etaPartition(id); 
  auto roll = me0Geom->chamber(id); 
  
  float zSign  = me0muon.globalTrackMomAtSurface().z()/fabs(me0muon.globalTrackMomAtSurface().z());
  if ( zSign * roll->toGlobal(thisSegment.localPosition()).z() < 0 ) return false;
	
  //GlobalPoint r3FinalReco_glob(r3FinalReco_globv.x(),r3FinalReco_globv.y(),r3FinalReco_globv.z());

  LocalPoint r3FinalReco = roll->toLocal(me0muon.globalTrackPosAtSurface());
  //LocalVector p3FinalReco=roll->toLocal(me0muon.globalTrackMomAtSurface());

  LocalPoint thisPosition(thisSegment.localPosition());
  LocalVector thisDirection(thisSegment.localDirection().x(),thisSegment.localDirection().y(),thisSegment.localDirection().z());  //FIXME

  //The same goes for the error
  AlgebraicMatrix thisCov(4,4,0);   
  for (int i = 1; i <=4; i++){
    for (int j = 1; j <=4; j++){
      thisCov(i,j) = thisSegment.parametersError()(i,j);
    }
  }

  /////////////////////////////////////////////////////////////////////////////////////////


  // LocalTrajectoryParameters ltp(r3FinalReco,p3FinalReco,me0muon.trackCharge());
  // JacobianCartesianToLocal jctl(roll->surface(),ltp);
  // AlgebraicMatrix56 jacobGlbToLoc = jctl.jacobian(); 
  // AlgebraicMatrix55 Ctmp =  (jacobGlbToLoc * me0muon.globalTrackCov()) * ROOT::Math::Transpose(jacobGlbToLoc); 
  AlgebraicMatrix55 Ctmp =  me0muon.localTrackCov(); 

  AlgebraicSymMatrix55 C;  // I couldn't find any other way, so I resort to the brute force
  for(int i=0; i<5; ++i) {
    for(int j=0; j<5; ++j) {
      C[i][j] = Ctmp[i][j]; 
      
    }
  }  

  Double_t sigmax = sqrt(C[3][3]+thisSegment.localPositionError().xx() );      
  Double_t sigmay = sqrt(C[4][4]+thisSegment.localPositionError().yy() );

  bool X_MatchFound = false, Y_MatchFound = false, Dir_MatchFound = false;
  
  
  if ( ( (fabs(thisPosition.x()-r3FinalReco.x())/sigmax ) < MaxPullX ) || (fabs(thisPosition.x()-r3FinalReco.x()) < MaxDiffX ) ) X_MatchFound = true;

  if ( MaxPullY < 0. && MaxDiffY < 0. ) Y_MatchFound = true;
  else if ( ( (fabs(thisPosition.y()-r3FinalReco.y())/sigmay ) < MaxPullY ) || (fabs(thisPosition.y()-r3FinalReco.y()) < MaxDiffY ) ) Y_MatchFound = true;
  
  if ( MaxDiffPhiDir < 0. ) Dir_MatchFound = true;
  else if ( fabs(me0muon.globalTrackMomAtSurface().phi()-roll->toGlobal(thisSegment.localDirection()).phi()) < MaxDiffPhiDir) Dir_MatchFound = true;

  return (X_MatchFound && Y_MatchFound && Dir_MatchFound);

}



void NewME0MuonGeneralAnalyzer::getSimHitsAndTracks( const Event& event, 
						     edm::Handle<SimTrackContainer>& simTracks, 
						     //SimHitsMap& pxbHiMap, 
						     //SimHitsMap& pxbLoMap, 
						     SimHitsMap& pxeHiMap, 
						     SimHitsMap& pxeLoMap, 
						     SimHitsMap& me0Map 
						     /*, SimTracksMap& simTracksMap*/ ) {
  
  // Get the SimHit collection from the event
  // //edm::EDGetTokenT<edm::PSimHitContainer> pxbHiToken = consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("pxbHiTag")); 
  // edm::EDGetTokenT<edm::PSimHitContainer> pxbHiToken = consumes<edm::PSimHitContainer>(edm::InputTag("g4SimHits", "TrackerHitsPixelBarrelHighTof")); 
  // edm::Handle<edm::PSimHitContainer> pxbHiH;
  // event.getByToken(pxbHiToken, pxbHiH);
  
  // edm::EDGetTokenT<edm::PSimHitContainer> pxbLoToken = consumes<edm::PSimHitContainer>(edm::InputTag("g4SimHits", "TrackerHitsPixelBarrelLowTof")); 
  // edm::Handle<edm::PSimHitContainer> pxbLoH;
  // event.getByToken(pxbLoToken, pxbLoH);

  //edm::EDGetTokenT<edm::PSimHitContainer> pxeHiToken = consumes<edm::PSimHitContainer>(edm::InputTag("g4SimHits", "TrackerHitsPixelEndcapHighTof")); 
  edm::Handle<edm::PSimHitContainer> pxeHiH;
  event.getByToken(pxeHiToken, pxeHiH);
  
  //edm::EDGetTokenT<edm::PSimHitContainer> pxeLoToken = consumes<edm::PSimHitContainer>(edm::InputTag("g4SimHits", "TrackerHitsPixelEndcapLowTof")); 
  edm::Handle<edm::PSimHitContainer> pxeLoH;
  event.getByToken(pxeLoToken, pxeLoH);

  //edm::EDGetTokenT<edm::PSimHitContainer> me0Token = consumes<edm::PSimHitContainer>(edm::InputTag("g4SimHits", "MuonME0Hits")); 
  edm::Handle<edm::PSimHitContainer> me0H;
  event.getByToken(me0Token, me0H);

  //edm::EDGetTokenT<SimTrackContainer> simTrackToken = consumes<SimTrackContainer>(edm::InputTag("g4SimHits")); 
  //edm::Handle<SimTrackContainer> simTracks;
  event.getByToken(simTrackToken, simTracks);


  // for(edm::PSimHitContainer::const_iterator simhit=pxbHiH->begin(); simhit!=pxbHiH->end(); ++simhit) {
  //   if(abs(simhit->particleType())!=13) continue;
  //   pxbHiMap[simhit->trackId()].push_back(&*simhit);
  // }

  // for(edm::PSimHitContainer::const_iterator simhit=pxbLoH->begin(); simhit!=pxbLoH->end(); ++simhit) {
  //   if(abs(simhit->particleType())!=13) continue;
  //   pxbLoMap[simhit->trackId()].push_back(&*simhit);
  // }

  for(edm::PSimHitContainer::const_iterator simhit=pxeHiH->begin(); simhit!=pxeHiH->end(); ++simhit) {
    if(abs(simhit->particleType())!=13) continue;
    pxeHiMap[simhit->trackId()].push_back(&*simhit);
  }

  for(edm::PSimHitContainer::const_iterator simhit=pxeLoH->begin(); simhit!=pxeLoH->end(); ++simhit) {
    if(abs(simhit->particleType())!=13) continue;
    pxeLoMap[simhit->trackId()].push_back(&*simhit);
  }

  for(edm::PSimHitContainer::const_iterator simhit=me0H->begin(); simhit!=me0H->end(); ++simhit) {
    if(abs(simhit->particleType())!=13) continue;
    me0Map[simhit->trackId()].push_back(&*simhit);
  }


  // GenParticle indexes 
  edm::Handle<std::vector<int> > genPartIdxs; 
  event.getByToken(gpidxToken_, genPartIdxs); 
  if(!genPartIdxs.isValid()) {
    throw cms::Exception("GenParticle indexes collection not valid!"); 
  } 
  else {
    if(debugPrint) std::cout << "Found GenParticle indexes collection!" << std::endl;
  }

  unsigned int simtrackIdx = 0; 
  for(SimTrackContainer::const_iterator simtrack=simTracks->begin(); simtrack!=simTracks->end(); ++simtrack, ++simtrackIdx) {
    if(abs(simtrack->type())!=13) continue;
    unsigned int trId = simtrack->trackId();
    if( pxeHiMap.count(trId)>0 || pxeLoMap.count(trId)>0 || me0Map.count(trId)>0 ) {
      trkIdIndexMap[trId] = simtrackIdx; 
      //simTracksMap[trId] = &*simtrack; 
    }

    int genpidx = simtrack->genpartIndex(); 
    if(genpidx<0) continue; 
    for(unsigned int gpcnt=0; gpcnt<genPartIdxs->size(); ++gpcnt) { 
      if(genPartIdxs->at(gpcnt)==genpidx) { 
	genpartIndexMap[genpidx] = gpcnt; 
	break; 
      } 
    } 
  }

  return;
}



void NewME0MuonGeneralAnalyzer::associateRecoToSimByHits( //SimHitsMap& pxbHiMap, SimHitsMap& pxbLoMap, 
							 SimHitsMap& pxeHiMap, SimHitsMap& pxeLoMap, 
							 SimHitsMap& me0Map,   ME0MuonCollection& recomus,
							 edm::Handle<reco::GenParticleCollection>& gpcoll, 
							 edm::Handle<SimTrackContainer>& stcoll) {

  SimRecoAssoMap simToRecoCompleteAssoMap;  // 1 SIM -> many RECO (1 RECO -> many SIM)
  SimRecoAssoMap simToRecoAssoMap;          // 1 SIM -> many RECO (1 RECO -> 1 SIM)

  unsigned int recoCnt = 0; 
  for(ME0MuonCollection::const_iterator itmu=recomus.begin(); itmu!=recomus.end(); ++itmu, ++recoCnt) {

    const reco::TrackRef recoTrack = itmu->innerTrack(); 
    const ME0Segment& muonSegm = itmu->me0segment(); 

    // bool matchFound = false;
    // unsigned int matchTrack = 0;
    std::map<unsigned int, unsigned int> assoMap;
    std::map<unsigned int, unsigned int> trkAssoMap;
    std::map<unsigned int, unsigned int> me0AssoMap;

    for(trackingRecHit_iterator recHit=recoTrack->recHitsBegin(); recHit!=recoTrack->recHitsEnd(); ++recHit) {
      if((*recHit)->isValid()==false) continue; 

      DetId detId((*recHit)->geographicalId());

      // --------------- Pixel Endcap --------------- 
      if(detId.det()==DetId::Tracker && detId.subdetId()==PixelSubdetector::PixelEndcap) { 
	PXFDetId pxeId((*recHit)->geographicalId()); 

	// ============= High TOF ============= 
	for(SimHitsMap::const_iterator mapIt=pxeHiMap.begin(); mapIt!=pxeHiMap.end(); ++mapIt) {
	  for(std::vector<const PSimHit*>::const_iterator shIt=(*mapIt).second.begin(); shIt!=(*mapIt).second.end(); ++shIt) {
	    PXFDetId simPxeDetId((*shIt)->detUnitId());
	    if( pxeId.side()   == simPxeDetId.side()   && 
		pxeId.disk()   == simPxeDetId.disk()   && 
		pxeId.blade()  == simPxeDetId.blade()  && 
		pxeId.panel()  == simPxeDetId.panel()  && 
		pxeId.module() == simPxeDetId.module() 
		) { 
	      double distX = (*recHit)->localPosition().x() - (*shIt)->localPosition().x();
	      double distY = (*recHit)->localPosition().y() - (*shIt)->localPosition().y();
	      double dist = sqrt(distX*distX + distY*distY);
	      hists_["hPxeHi2DHitResiduals"]->Fill(dist);
	      hists_["hPxeHiXHitResiduals" ]->Fill(distX);
	      hists_["hPxeHiYHitResiduals" ]->Fill(distY);
	      // if((*recHit)->isValid()) {
	      // 	hists_["hPxeHi2DValidHitResiduals"]->Fill(dist);
	      // 	hists_["hPxeHiXValidHitResiduals" ]->Fill(distX);
	      // 	hists_["hPxeHiYValidHitResiduals" ]->Fill(distY);
	      // }
	      // else {
	      // 	hists_["hPxeHi2DInvalidHitResiduals"]->Fill(dist);
	      // 	hists_["hPxeHiXInvalidHitResiduals" ]->Fill(distX);
	      // 	hists_["hPxeHiYInvalidHitResiduals" ]->Fill(distY);
	      // }
	      if(dist<3.5) {
		trkAssoMap[(*mapIt).first]++;
		assoMap[(*mapIt).first]++;
	      }
	    }
	  } // end for(std::vector<const PSimHit*>::const_iterator shIt=(*mapIt).second.begin(); shIt!=(*mapIt).second.end(); ++shIt)
	} // end for(SimHitsMap::const_iterator mapIt=pxeHiMap.begin(); mapIt!=pxeHiMap.end(); ++mapIt)

	// ============= Low TOF ============= 
	for(SimHitsMap::const_iterator mapIt=pxeLoMap.begin(); mapIt!=pxeLoMap.end(); ++mapIt) {
	  for(std::vector<const PSimHit*>::const_iterator shIt=(*mapIt).second.begin(); shIt!=(*mapIt).second.end(); ++shIt) {
	    PXFDetId simPxeDetId((*shIt)->detUnitId());
	    if( pxeId.side()   == simPxeDetId.side()   && 
		pxeId.disk()   == simPxeDetId.disk()   && 
		pxeId.blade()  == simPxeDetId.blade()  && 
		pxeId.panel()  == simPxeDetId.panel()  && 
		pxeId.module() == simPxeDetId.module() 
		) { 
	      double distX = (*recHit)->localPosition().x() - (*shIt)->localPosition().x();
	      double distY = (*recHit)->localPosition().y() - (*shIt)->localPosition().y();
	      double dist = sqrt(distX*distX + distY*distY);
	      hists_["hPxeLo2DHitResiduals"]->Fill(dist);
	      hists_["hPxeLoXHitResiduals" ]->Fill(distX);
	      hists_["hPxeLoYHitResiduals" ]->Fill(distY);
	      // if((*recHit)->isValid()) {
	      // 	hists_["hPxeLo2DValidHitResiduals"]->Fill(dist);
	      // 	hists_["hPxeLoXValidHitResiduals" ]->Fill(distX);
	      // 	hists_["hPxeLoYValidHitResiduals" ]->Fill(distY);
	      // }
	      // else {
	      // 	hists_["hPxeLo2DInvalidHitResiduals"]->Fill(dist);
	      // 	hists_["hPxeLoXInvalidHitResiduals" ]->Fill(distX);
	      // 	hists_["hPxeLoYInvalidHitResiduals" ]->Fill(distY);
	      // }
	      if(dist<3.5) {
		trkAssoMap[(*mapIt).first]++;
		assoMap[(*mapIt).first]++;
	      }
	    }
	  } // end for(std::vector<const PSimHit*>::const_iterator shIt=(*mapIt).second.begin(); shIt!=(*mapIt).second.end(); ++shIt)
	} // end for(SimHitsMap::const_iterator mapIt=pxeHiMap.begin(); mapIt!=pxeHiMap.end(); ++mapIt)

      } // end if(detId.subdetId()==PixelSubdetector::PixelEndcap)
    } // end for(trackingRecHit_iterator recHit=recoTrack->recHitsBegin(); recHit!=recoTrack->recHitsEnd(); ++recHit)


    auto me0rhs = muonSegm.specificRecHits(); 
    for(auto recHit=me0rhs.begin(); recHit!=me0rhs.end(); ++recHit) {
      if(recHit->isValid()==false) continue; 

      // DetId detId((*recHit)->me0Id()); 

      // // --------------- ME0 --------------- 
      // if(detId.det()==DetId::Muon && detId.subdetId()==MuonSubdetId::ME0) {
      ME0DetId me0Id(recHit->me0Id()); 

      // ============= ME0 ============= 
      for(SimHitsMap::const_iterator mapIt=me0Map.begin(); mapIt!=me0Map.end(); ++mapIt) {
	for(std::vector<const PSimHit*>::const_iterator shIt=(*mapIt).second.begin(); shIt!=(*mapIt).second.end(); ++shIt) {
	  ME0DetId simMe0DetId((*shIt)->detUnitId()); 
	  if( me0Id.region()  == simMe0DetId.region()  && 
	      me0Id.layer()   == simMe0DetId.layer()   && 
	      me0Id.chamber() == simMe0DetId.chamber() && 
	      me0Id.roll()    == simMe0DetId.roll() 
	      ) { 
	    double distX = recHit->localPosition().x() - (*shIt)->localPosition().x();
	    double distY = recHit->localPosition().y() - (*shIt)->localPosition().y();
	    double dist = sqrt(distX*distX + distY*distY);
	    hists_["hMe0Hi2DHitResiduals"]->Fill(dist);
	    hists_["hMe0HiXHitResiduals" ]->Fill(distX);
	    hists_["hMe0HiYHitResiduals" ]->Fill(distY);
	    // if(recHit->isValid()) {
	    //   hists_["hMe0Hi2DValidHitResiduals"]->Fill(dist);
	    //   hists_["hMe0HiXValidHitResiduals" ]->Fill(distX);
	    //   hists_["hMe0HiYValidHitResiduals" ]->Fill(distY);
	    // }
	    // else {
	    //   hists_["hMe0Hi2DInvalidHitResiduals"]->Fill(dist);
	    //   hists_["hMe0HiXInvalidHitResiduals" ]->Fill(distX);
	    //   hists_["hMe0HiYInvalidHitResiduals" ]->Fill(distY);
	    // }
	    if(dist<1.0) {
	      me0AssoMap[(*mapIt).first]++;
	      assoMap[(*mapIt).first]++;
	    }
	  }
	} // end for(std::vector<const PSimHit*>::const_iterator shIt=(*mapIt).second.begin(); shIt!=(*mapIt).second.end(); ++shIt)
      } // end for(SimHitsMap::const_iterator mapIt=me0HiMap.begin(); mapIt!=me0HiMap.end(); ++mapIt)
	// } // end if(detId.det()==DetId::Muon && detId.subdetId()==MuonSubdetId::ME0)
    } // end for(auto recHit=me0rhs.begin(); recHit!=me0rhs.end(); ++rh) 

    // Various plots 
    for(std::map<unsigned int, unsigned int>::const_iterator itAssoMap=trkAssoMap.begin(); itAssoMap!=trkAssoMap.end(); ++itAssoMap) { 
      hists_["numberMatchesPerTrack_tracker_all"]->Fill((*itAssoMap).second); 
      hists_["fractnMatchesPerTrack_tracker_all"]->Fill(float((*itAssoMap).second)/recoTrack->numberOfValidHits()); 
      if(trkIdIndexMap.count((*itAssoMap).first)>0) {
	const SimTrack& ast = stcoll->at( trkIdIndexMap[(*itAssoMap).first] ); 
	if(ast.genpartIndex()>=0) { 
	  const reco::GenParticle& agp = gpcoll->at( genpartIndexMap[ast.genpartIndex()] ); 
	  if(agp.isPromptFinalState()) { 
	    hists_["numberMatchesPerTrack_sign_tracker_all"]->Fill((*itAssoMap).second); 
	    hists_["fractnMatchesPerTrack_sign_tracker_all"]->Fill(float((*itAssoMap).second)/recoTrack->numberOfValidHits()); 
	  } 
	} 
      }
    } 

    for(std::map<unsigned int, unsigned int>::const_iterator itAssoMap=me0AssoMap.begin(); itAssoMap!=me0AssoMap.end(); ++itAssoMap) { 
      hists_["numberMatchesPerTrack_me0segm_all"]->Fill((*itAssoMap).second); 
      hists_["fractnMatchesPerTrack_me0segm_all"]->Fill(float((*itAssoMap).second)/muonSegm.nRecHits()); 
      if(trkIdIndexMap.count((*itAssoMap).first)>0) {
	const SimTrack& ast = stcoll->at( trkIdIndexMap[(*itAssoMap).first] ); 
	if(ast.genpartIndex()>=0) { 
	  const reco::GenParticle& agp = gpcoll->at( genpartIndexMap[ast.genpartIndex()] ); 
	  if(agp.isPromptFinalState()) { 
	    hists_["numberMatchesPerTrack_sign_me0segm_all"]->Fill((*itAssoMap).second); 
	    hists_["fractnMatchesPerTrack_sign_me0segm_all"]->Fill(float((*itAssoMap).second)/muonSegm.nRecHits()); 
	  } 
	} 
      }
    } 

    for(std::map<unsigned int, unsigned int>::const_iterator itAssoMap=assoMap.begin(); itAssoMap!=assoMap.end(); ++itAssoMap) { 
      hists_["numberMatchesPerTrack_wholemu_all"]->Fill((*itAssoMap).second); 
      hists_["fractnMatchesPerTrack_wholemu_all"]->Fill(float((*itAssoMap).second)/(recoTrack->numberOfValidHits()+muonSegm.nRecHits())); 
      if(trkIdIndexMap.count((*itAssoMap).first)>0) {
	const SimTrack& ast = stcoll->at( trkIdIndexMap[(*itAssoMap).first] ); 
	if(ast.genpartIndex()>=0) { 
	  const reco::GenParticle& agp = gpcoll->at( genpartIndexMap[ast.genpartIndex()] ); 
	  if(agp.isPromptFinalState()) { 
	    hists_["numberMatchesPerTrack_sign_wholemu_all"]->Fill((*itAssoMap).second); 
	    hists_["fractnMatchesPerTrack_sign_wholemu_all"]->Fill(float((*itAssoMap).second)/(recoTrack->numberOfValidHits()+muonSegm.nRecHits())); 
	  } 
	} 
      }
    } 

    if(assoMap.size()!=0) { 
      unsigned int bestSim = 0;
      unsigned int bestRec = 0;
      unsigned int bestMatch = 0;
      //for(map<unsigned int, unsigned int>::const_iterator itAssoMap=trkAssoMap.begin(); itAssoMap!=trkAssoMap.end(); ++itAssoMap) {
      //for(map<unsigned int, unsigned int>::const_iterator itAssoMap=me0AssoMap.begin(); itAssoMap!=me0AssoMap.end(); ++itAssoMap) {
      for(std::map<unsigned int, unsigned int>::const_iterator itAssoMap=assoMap.begin(); itAssoMap!=assoMap.end(); ++itAssoMap) {
	simToRecoCompleteAssoMap[(*itAssoMap).first].push_back(std::make_pair(recoCnt, (*itAssoMap).second)); 
	if((*itAssoMap).second>bestMatch) {
	  bestSim = (*itAssoMap).first;
	  bestRec = recoCnt;
	  bestMatch = (*itAssoMap).second;
	}
      } // end for(map<unsigned int, unsigned int>::const_iterator itAssoMap=assoMap.begin(); itAssoMap!=assoMap.end(); ++itAssoMap) 
      simToRecoAssoMap[bestSim].push_back( std::make_pair( bestRec, bestMatch) ); 

      hists_["numberMatchesPerTrack_tracker_best"]->Fill(trkAssoMap[bestSim]); 
      hists_["fractnMatchesPerTrack_tracker_best"]->Fill(float(trkAssoMap[bestSim])/recoTrack->numberOfValidHits()); 

      hists_["numberMatchesPerTrack_me0segm_best"]->Fill(me0AssoMap[bestSim]); 
      hists_["fractnMatchesPerTrack_me0segm_best"]->Fill(float(me0AssoMap[bestSim])/muonSegm.nRecHits()); 

      hists_["numberMatchesPerTrack_wholemu_best"]->Fill(bestMatch); 
      hists_["fractnMatchesPerTrack_wholemu_best"]->Fill(float(bestMatch)/(recoTrack->numberOfValidHits()+muonSegm.nRecHits())); 

      if(trkIdIndexMap.count(bestSim)>0) {
	const SimTrack& ast = stcoll->at( trkIdIndexMap[bestSim] ); 
	if(ast.genpartIndex()>=0) { 
	  const reco::GenParticle& agp = gpcoll->at( genpartIndexMap[ast.genpartIndex()] ); 
	  if(agp.isPromptFinalState()) { 
	    hists_["numberMatchesPerTrack_sign_tracker_best"]->Fill(trkAssoMap[bestSim]); 
	    hists_["fractnMatchesPerTrack_sign_tracker_best"]->Fill(float(trkAssoMap[bestSim])/recoTrack->numberOfValidHits()); 

	    hists_["numberMatchesPerTrack_sign_me0segm_best"]->Fill(me0AssoMap[bestSim]); 
	    hists_["fractnMatchesPerTrack_sign_me0segm_best"]->Fill(float(me0AssoMap[bestSim])/muonSegm.nRecHits()); 

	    hists_["numberMatchesPerTrack_sign_wholemu_best"]->Fill(bestMatch); 
	    hists_["fractnMatchesPerTrack_sign_wholemu_best"]->Fill(float(bestMatch)/(recoTrack->numberOfValidHits()+muonSegm.nRecHits())); 
	  } 
	} 
      }
    } // end if(assoMap.size()!=0) 

  } // end for(ME0MuonCollection::const_iterator itmu=recomus.begin(); itmu!=recomus.end(); ++itmu, ++recoCnt)

  if(simToRecoAssoMap.empty()) 
    return;

  sortByNumberOfMatches(simToRecoCompleteAssoMap); 
  sortByNumberOfMatches(simToRecoAssoMap); 

  for(SimRecoAssoMap::const_iterator simRecIt=simToRecoAssoMap.begin(); simRecIt!=simToRecoAssoMap.end(); ++simRecIt) { 
    //simToRecoMap[ trkIdIndexMap[(*simRecIt).first] ] = (simRecIt->second)[0].first; 
    simToRecoMap[(*simRecIt).first] = (simRecIt->second)[0].first; 
    recoToSimMap[(simRecIt->second)[0].first] = (*simRecIt).first; 
  } 

  return;
}



void NewME0MuonGeneralAnalyzer::sortByNumberOfMatches(SimRecoAssoMap& mapToSort) {

  SimRecoAssoMap tmpMap;
  //std::cout << " [sortByNumberOfMatches]: start: " << mapToSort.size() << std::endl;

  for(SimRecoAssoMap::const_iterator simRecIt=mapToSort.begin(); simRecIt!=mapToSort.end(); ++simRecIt) {
    //std::cout << " [sortByNumberOfMatches]:        before " << (*simRecIt).first << " - " << (*simRecIt).second.size() << std::endl;
    std::vector<std::pair<unsigned int, unsigned int> > tmpVec = (*simRecIt).second;

    //for(unsigned int ii=0; ii!=tmpVec.size(); ++ii) {
    //  std::cout << " [sortByNumberOfMatches]:                   recoCnt " << tmpVec[ii].first << "  matches " << tmpVec[ii].second << std::endl;
    //}

    for(unsigned int i1=0; i1<tmpVec.size()-1; ++i1) {
      for(unsigned int i2=i1+1; i2<tmpVec.size(); ++i2) {
	if(tmpVec[i1].second<tmpVec[i2].second) {
	  std::pair<unsigned int, unsigned int> auxPair = tmpVec[i1];
	  tmpVec[i1] = tmpVec[i2];
	  tmpVec[i2] = auxPair;
	}
      }
    }
    tmpMap[(*simRecIt).first] = tmpVec;
    //std::cout << " [sortByNumberOfMatches]:        after  " << tmpVec.size() << " - " << tmpMap[(*simRecIt).first].size() << std::endl;
  }
  mapToSort.clear();
  mapToSort = tmpMap;
  //std::cout << " [sortByNumberOfMatches]:   end: " << mapToSort.size() << std::endl;
  return;
} 

// define this as a plug-in
DEFINE_FWK_MODULE(NewME0MuonGeneralAnalyzer);
