// -*- C++ -*-
//
// Package:    RecoMuonGeneralAnalyzer
// Class:      RecoMuonGeneralAnalyzer
// 
/**\class RecoMuonGeneralAnalyzer 

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
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/Provenance/interface/ProductID.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
//#include "DataFormats/MuonReco/interface/MuonSelectors.h"
//#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h" 
//#include "DataFormats/DetId/interface/DetId.h"
//#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
//#include "DataFormats/RPCRecHit/interface/RPCRecHit.h"
//#include "DataFormats/MuonDetId/interface/RPCDetId.h"
//#include "DataFormats/MuonReco/interface/ME0Muon.h" 
//#include "DataFormats/MuonReco/interface/ME0MuonCollection.h" 
#include "DataFormats/GEMRecHit/interface/ME0SegmentCollection.h"
//#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
//#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/TrajectoryState/interface/LocalTrajectoryParameters.h"
#include "DataFormats/MuonDetId/interface/ME0DetId.h"
#include "DataFormats/Math/interface/deltaPhi.h" 

// Tracking tools 
#include "TrackingTools/AnalyticalJacobians/interface/JacobianCartesianToLocal.h"
#include "TrackingTools/AnalyticalJacobians/interface/JacobianLocalToCartesian.h"

// Muon geometry (some are probably not needed) 
#include "Geometry/GEMGeometry/interface/ME0Geometry.h"
#include "Geometry/GEMGeometry/interface/ME0EtaPartition.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"

// Muon ID 
//#include "RecoMuon/MuonIdentification/interface/ME0MuonSelector.h" 

// Associator by hits or chi2 
//#include "SimTracker/TrackAssociation/interface/TrackAssociatorByChi2.h"
//#include "SimTracker/TrackAssociation/interface/TrackAssociatorByHits.h"
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


class RecoMuonGeneralAnalyzer : public edm::EDAnalyzer {
public:
  explicit RecoMuonGeneralAnalyzer(const edm::ParameterSet&);
  ~RecoMuonGeneralAnalyzer();

  typedef std::map<unsigned int, std::vector<unsigned int> > MuonSegmentMap; 

private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  //bool isGoodME0Muon(edm::ESHandle<ME0Geometry>, const reco::ME0Muon&, double, double, double, double, double); 

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
  std::vector<std::string> associatorNames; 
  //std::vector<const TrackAssociatorBase*> associators; 

  // Tokens 
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> > puInfoTag_;
  edm::EDGetTokenT<ME0SegmentCollection> me0SegmentsToken_;
  edm::EDGetTokenT<TrackingParticleCollection> tpcToken_; 
  //edm::EDGetTokenT<ME0MuonCollection> me0MuonToken_; 
  edm::EDGetTokenT<reco::MuonCollection> muonToken_; 
  edm::EDGetTokenT<edm::View<reco::Track> > tracksToken_; 
  edm::EDGetTokenT<MuonSegmentMap> muSegmMapToken_; 

  // Track-ME0Segment selection cuts 
  double minP; 
  double minPt; 
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
RecoMuonGeneralAnalyzer::RecoMuonGeneralAnalyzer(const edm::ParameterSet& iConfig) :
  puInfoTag_( consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("puInfoTag")) ), 
  me0SegmentsToken_( consumes<ME0SegmentCollection>(iConfig.getParameter<edm::InputTag>("me0segmentsTag")) ), 
  tpcToken_( consumes<TrackingParticleCollection>(iConfig.getParameter<edm::InputTag>("tpCollectionTag")) ), 
  //me0MuonToken_( consumes<ME0MuonCollection>(iConfig.getParameter<edm::InputTag>("me0muonTag")) ), 
  muonToken_( consumes<reco::MuonCollection>(iConfig.getParameter<edm::InputTag>("muonTag")) ), 
  tracksToken_( consumes<edm::View<reco::Track> >(iConfig.getParameter<edm::InputTag>("tracksTag")) ), 
  muSegmMapToken_( consumes<MuonSegmentMap>(iConfig.getParameter<edm::InputTag>("muSegmMapTag")) ) 
{
  if(debugPrint) std::cout << "Inside Constructor" << std::endl;

  debugPrint       = iConfig.getUntrackedParameter<bool>("verbose"          , false); 
  outputPrint      = iConfig.getUntrackedParameter<bool>("printOutput"      , false); 

  associatorNames = iConfig.getParameter<std::vector<std::string> >("associatorNames"); 

  for (auto const& thisassociator : associatorNames) {
    //std::cout << thisassociator << std::endl;
    consumes<reco::RecoToSimCollection>(edm::InputTag(thisassociator));
    consumes<reco::SimToRecoCollection>(edm::InputTag(thisassociator));
  } 

  minP          = iConfig.getUntrackedParameter<double>("minP"         , -1.0  ); 
  minPt         = iConfig.getUntrackedParameter<double>("minPt"        , -1.0  ); 


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


RecoMuonGeneralAnalyzer::~RecoMuonGeneralAnalyzer() { 
  if(debugPrint) std::cout << "Inside Destructor" << std::endl;
}


//
// member functions
//

// ------------ method called for each event  ------------
void RecoMuonGeneralAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
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


  // Tracking particles (preselected) 
  edm::Handle<TrackingParticleCollection> trackParts; // trackParticles; 
  iEvent.getByToken(tpcToken_, trackParts); 
  if(!trackParts.isValid()) {
    if(outputPrint || debugPrint) std::cout << "TrackingParticle collection not valid!" << std::endl;
    return;
  } 
  else {
    if(debugPrint) std::cout << "Found TrackingParticle collection!" << std::endl;
  }

  // edm::Handle<ME0MuonCollection> muons; 
  // iEvent.getByToken(me0MuonToken_, muons); 
  // if(!muons.isValid()) {
  //   if(outputPrint || debugPrint) std::cout << "ME0Muon collection not valid!" << std::endl;
  //   return;
  // } 
  // else {
  //   if(debugPrint) std::cout << "Found ME0Muon collection!" << std::endl;
  // }

  // Muons (all) 
  edm::Handle<reco::MuonCollection> muons; 
  iEvent.getByToken(muonToken_, muons); 
  if(!muons.isValid()) {
    if(outputPrint || debugPrint) std::cout << "Muon collection not valid!" << std::endl;
    return;
  } 
  else {
    if(debugPrint) std::cout << "Found ME0Muon collection!" << std::endl;
  }

  // Tracks of ME0 muons (preselected) 
  edm::Handle<edm::View<reco::Track> > trackColl; // trackCollection; 
  iEvent.getByToken(tracksToken_, trackColl); 
  if(!trackColl.isValid()) {
    if(outputPrint || debugPrint) std::cout << "Track collection not valid!" << std::endl;
    return;
  } 
  else {
    if(debugPrint) std::cout << "Found Track collection!" << std::endl;
  }

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
    if(outputPrint || debugPrint) std::cout << "ME0Segment collection not valid!" << std::endl;
    return;
  } 
  else {
    if(debugPrint) std::cout << "Found ME0Segment collection!" << std::endl;
  }


  edm::Handle<MuonSegmentMap> muonSegmMapH; 
  iEvent.getByToken(muSegmMapToken_, muonSegmMapH); 
  if(!muonSegmMapH.isValid()) {
    if(outputPrint || debugPrint) std::cout << "MuonSegmentMapH not valid!" << std::endl;
    return;
  } 
  else {
    if(debugPrint) std::cout << "Found MuonSegmentMapH!" << std::endl;
  }
  MuonSegmentMap muonSegmMap = (*muonSegmMapH); 

  assert( muonSegmMap.size()==(*trackColl).size() ); 

  std::map<unsigned int, unsigned int> trackToMuonMap; 

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


  reco::RecoToSimCollection recSimColl; 
  reco::SimToRecoCollection simRecColl; 
  //TrackingParticleRefVector trackParts; // a.k.a. edm::RefVector<TrackingParticleCollection> 
  //edm::RefToBaseVector<reco::Track> trackColl; 
  //std::vector<bool> passIds; 

  // for(size_t h=0; h<trackParticles->size(); ++h) { 
  //   if((*trackParticles)[h].status()!=1) continue; 
  //   if(selectMuonOnlyTP && abs((*trackParticles)[h].pdgId())!=13) continue; 
  //   trackParts.push_back(TrackingParticleRef(trackParticles, h)); 
  // } 

  for(size_t j=0; j<trackColl->size(); ++j) { 
    //edm::ProductID tkid = trackColl->refAt(j).id(); 
    const reco::Track *tkptr = trackColl->refAt(j).get(); 

    bool isselected = false; 
    unsigned int muonidx = 0; 
    for(reco::MuonCollection::const_iterator muon=muons->begin(); muon!=muons->end(); ++muon, muonidx++) { 
      // Muon ID 
      if(!muon->isME0Muon()) continue; 

      const reco::Track *muptr = muon->innerTrack().get(); 

      if(false && debugPrint) {
  	std::cout << "*** Track vs Muon: " << std::endl; 
	std::cout << "    Muon " << muonidx << ": (" << muptr->px() << ", " << muptr->py() << ", " << muptr->pz() << ") " 
		  << " - Track " << j       << ": (" << tkptr->px() << ", " << tkptr->py() << ", " << tkptr->pz() << ") " << std::endl; 
  	//std::cout << " Muon: " << muptr << " - Track " << j << ": " << tkptr << std::endl; 
  	//std::cout << " Muon: " << muon->innerTrack().id() << " - Track " << j << ": " << trackColl->refAt(j).id() << std::endl; 
      } 
      //if(muon->innerTrack() == trackColl->refAt(j)) { // this DOESN'T work! Different types of Ref! 
      isselected = 
	(muptr->charge()==tkptr->charge())   && 
	(fabs(muptr->px()-tkptr->px())<1E-4) && 
	(fabs(muptr->py()-tkptr->py())<1E-4) && 
	(fabs(muptr->pz()-tkptr->pz())<1E-4) && 
	(fabs(muptr->vx()-tkptr->vx())<1E-4) && 
	(fabs(muptr->vy()-tkptr->vy())<1E-4) && 
	(fabs(muptr->vz()-tkptr->vz())<1E-4)  ; 
      if(isselected) break; // with this break, muonidx does not increment 

      // if(muptr==tkptr) { 
      // 	isselected = true; 
      // 	trackToMuonMap[j] = muonidx; 
      // 	break; 
      // } 
    } 

    if(!isselected) {  
	if(outputPrint || debugPrint) std::cout << "ME0 pre-selected track " << j << " can't be matched to any reco::Muon!" << std::endl;
	return;      
    } 
    else { 
      if(muonSegmMap.count(muonidx)==0) { 
	if(outputPrint || debugPrint) std::cout << "Muon " << muonidx << " is an ME0 muon, yet it can't be found in the MuonSegmentMap!" << std::endl;
	return;
      } 

      std::vector<unsigned int> thismusegms = muonSegmMap[muonidx]; 
      for(unsigned int ii=0; ii<thismusegms.size(); ++ii) { 
	hists_["nmu_perseg_num_aeta"]->Fill( segmetaphi[thismusegms[ii]].first  ); 
	hists_["nmu_perseg_num_dphi"]->Fill( segmetaphi[thismusegms[ii]].second ); 
      } 
    } 
  } 

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

  if((*trackParts).size()==0 || (*trackColl).size()==0) { 
    recSimColl.post_insert(); 
    simRecColl.post_insert(); 
  } 
  else { 
    edm::Handle<reco::SimToRecoCollection> simtorecoCollectionH;
    iEvent.getByLabel(associatorNames[0], simtorecoCollectionH);
    simRecColl = *(simtorecoCollectionH.product());

    edm::Handle<reco::RecoToSimCollection> recotosimCollectionH;
    iEvent.getByLabel(associatorNames[0], recotosimCollectionH);
    recSimColl = *(recotosimCollectionH.product());

    // recSimColl = associators[0]->associateRecoToSim(trackColl, 
    // 						    trackParts, 
    // 						    &iEvent, 
    // 						    &iSetup); 

    // simRecColl = associators[0]->associateSimToReco(trackColl, 
    // 						    trackParts, 
    // 						    &iEvent, 
    // 						    &iSetup); 
  } 

  // Efficiencies: Sim-Reco 
  for(size_t h=0; h<(*trackParts).size(); ++h) { 
    //TrackingParticleRef thisTp = (*trackParts)[h]; // a.k.a. edm::Ref<TrackingParticleCollection> 
    TrackingParticleRef thisTp(trackParts, h); 

    // Select only muons of status 1 (unless already requested) 
    //if(abs(thisTp->pdgId())!=13 || thisTp->status()!=1) continue; 

    // Large |eta| region for efficiency vs eta 
    if(fabs(thisTp->eta())<1.4 || fabs(thisTp->eta())>3.4) continue; 

    // int motherId(0), grandmaId(0), greatgrandmaId(0); 
    // if( thisTp->genParticles().size()>0 && (*thisTp->genParticle_begin())->numberOfMothers()>0 ) { 
    //   motherId = abs( (*thisTp->genParticle_begin())->mother()->pdgId() ); 
    //   if( (*thisTp->genParticle_begin())->mother()->numberOfMothers()>0 ) { 
    // 	grandmaId = abs( (*thisTp->genParticle_begin())->mother()->mother()->pdgId() ); 
    // 	if( (*thisTp->genParticle_begin())->mother()->mother()->numberOfMothers()>0 ) { 
    // 	  greatgrandmaId = abs( (*thisTp->genParticle_begin())->mother()->mother()->mother()->pdgId() ); 
    // 	} 
    //   } 
    // } 
    // // Is it a signal muon? 
    // bool isSignal = ( grandmaId==theMotherId || greatgrandmaId==theMotherId || 
    // 		      (motherId==13 && grandmaId==13 && greatgrandmaId==13)   ); 
    // if(!isSignal) continue; 

    // ONLY signal muons! 
    if(thisTp->genParticles().size()==0) 
      continue; 
    if((*thisTp->genParticle_begin())->isPromptFinalState()==false) 
      continue; 

    // Fill plots 
    if(thisTp->p()>minP || thisTp->pt()>minPt) { 
      // hists_["eff_den_eta" ]->Fill( thisTp->eta() );
      hists_["eff_den_aeta"]->Fill( fabs(thisTp->eta()) );
    } 
    if(fabs(thisTp->eta())>2.0 && fabs(thisTp->eta())<2.8) { // only within ME0 |eta| fiducial region 
      hists_["eff_den_pt" ]->Fill( thisTp->pt()  ); 
      if(thisTp->p()>minP || thisTp->pt()>minPt) { 
	hists_["eff_den_phi"]->Fill( thisTp->phi() ); 
	//hists_["eff_den_truepu"]->Fill( truePU ); 
	// hists_["eff_den_itpu"]->Fill( npuIT ); 
	hists_["eff_den_allpu"]->Fill( npuIT+npuOOT ); 
      } 
    } 
		       
    bool matchFound(false); 
    if(simRecColl.find(thisTp)!=simRecColl.end()) { 
      if(simRecColl[thisTp].size()>0) { 
	edm::RefToBase<reco::Track> assoTrack = simRecColl[thisTp].begin()->first; 

	// // Muon ID 
	// size_t itrk = 0; 
	// for(; itrk<trackColl.size(); ++itrk) {
	//   if(trackColl[itrk] == assoTrack) break; 
	//   //if(trackColl[itrk].get() == assoTrack.get()) break; // this works too 
	// } 
	// if(itrk>=trackColl.size()) { 
	//   std::cout << " *********************************************************\n" 
	//             << " *** ERROR: Couldn't find the associated ME0Muon track ***\n" 
	// 	    << " ***        inside the ME0Muon track collection!       ***\n" 
	// 	    << " ***        This should never happen!!                 ***\n" 
	//             << " *********************************************************" << std::endl; 
	//   continue; 
	// } 	
	// if(!passIds[itrk]) continue; 

	if(recSimColl.find(assoTrack)!=recSimColl.end()) { 
	  if(recSimColl[assoTrack].size()>0) { 
	    TrackingParticleRef assoTp = recSimColl[assoTrack].begin()->first; 
	    //matchFound = (assoTp.id() == thisTp.id()); 
	    matchFound = (assoTp == thisTp); 
	  } 
	} 
      } 
    } 

    // Match found! 
    if(matchFound) { 
      // Fill plots 
      if(thisTp->p()>minP || thisTp->pt()>minPt) { 
	// hists_["eff_num_eta" ]->Fill( thisTp->eta() ); 
	hists_["eff_num_aeta"]->Fill( fabs(thisTp->eta()) ); 
      } 
      if(fabs(thisTp->eta())>2.0 && fabs(thisTp->eta())<2.8) { 
	hists_["eff_num_pt" ]->Fill( thisTp->pt()  ); 
	if(thisTp->p()>minP || thisTp->pt()>minPt) { 
	  hists_["eff_num_phi"]->Fill( thisTp->phi() ); 
	  //hists_["eff_num_truepu"]->Fill( truePU ); 
	  // hists_["eff_num_itpu"]->Fill( npuIT ); 
	  hists_["eff_num_allpu"]->Fill( npuIT+npuOOT ); 
	} 
      } 
    } 
  } 

  // Fakes: Reco-Sim 
  for(unsigned int h=0; h<(*trackColl).size(); ++h) { 
    //edm::RefToBase<reco::Track> thisTrack = (*trackColl)[h]; 
    edm::RefToBase<reco::Track> thisTrack = trackColl->refAt(h); 
    unsigned int numberOfMatches = (muonSegmMap[ trackToMuonMap[h] ]).size(); // number of ME0 segments matched to a single track 

    //char buff[99]; 
    //printf("Muon track %u, %x -- (%.3f, %.3f, %.3f, %.3f)\n", h, (unsigned int)(thisTrack.get()), thisTrack->p(), thisTrack->px(), thisTrack->py(), thisTrack->pz()); 
    //printf("--- %u) Muon track -- (%.3f, %.3f, %.3f, %.3f)\n", h, thisTrack->p(), thisTrack->px(), thisTrack->py(), thisTrack->pz()); 
    if(debugPrint) printf("--- %u) Muon track -- (%.3f, %.3f, %.3f) -- matches: %d\n", h, thisTrack->p(), thisTrack->eta(), thisTrack->phi(), numberOfMatches); 
    //if(debug) std::cout << buff; 

    //if(fabs(thisTrack->eta())<2.0 || fabs(thisTrack->eta())>2.8) continue; 

    // // Muon ID 
    // if(!passIds[h]) continue; 

    // Fill plots 
    // // hists_["fake_pertrack_den_eta"   ]->Fill( thisTrack->eta()       ); 
    // if(fabs(thisTrack->eta())>2.0 && fabs(thisTrack->eta())<2.8) { 
    //   hists_["fake_pertrack_den_pt"    ]->Fill( thisTrack->pt()        ); 
    //   hists_["fake_pertrack_den_aeta"  ]->Fill( fabs(thisTrack->eta()) ); 
    //   hists_["fake_pertrack_den_phi"   ]->Fill( thisTrack->phi()       ); 
    //   //hists_["fake_pertrack_den_truepu"]->Fill( truePU                 ); 
    //   // hists_["fake_pertrack_den_itpu"  ]->Fill( npuIT                  ); 
    //   hists_["fake_pertrack_den_allpu" ]->Fill( npuIT+npuOOT           ); 
    // }

    for(unsigned int k=0; k<numberOfMatches; ++k) { 
      if(fabs(thisTrack->eta())>2.0 && fabs(thisTrack->eta())<2.8) { 
	hists_["fake_den_pt"    ]->Fill( thisTrack->pt()        ); 
	hists_["fake_den_aeta"  ]->Fill( fabs(thisTrack->eta()) ); 
	hists_["fake_den_phi"   ]->Fill( thisTrack->phi()       ); 
	//hists_["fake_den_truepu"]->Fill( truePU                 ); 
	// hists_["fake_den_itpu"  ]->Fill( npuIT                  ); 
	hists_["fake_den_allpu" ]->Fill( npuIT+npuOOT           ); 
      }
    } 

    bool matchFound(false); 
    bool matchOnlySignalFound(false); 
    if(recSimColl.find(thisTrack)!=recSimColl.end()) { 
      if(recSimColl[thisTrack].size()>0) { 
	TrackingParticleRef assoTp = recSimColl[thisTrack].begin()->first; 
	//printf("       Gen track -- (%.3f, %.3f, %.3f, %.3f)\n", assoTp->p(), assoTp->px(), assoTp->py(), assoTp->pz()); 

	// This is not enough to identify a prompt muon! 
	// Need to check number of mothers and pdgId of mothers! 
	//  (Cross-check with Cesare) 
	//if(abs(assoTp->pdgId())!=13 || assoTp->status()!=1) continue; 
	//if(fabs(assoTp->eta())<2.0 || fabs(assoTp->eta())>2.8) continue; 

	if(simRecColl.find(assoTp)!=simRecColl.end()) { 
	  if(simRecColl[assoTp].size()>0) { 
	    edm::RefToBase<reco::Track> assoTrack = simRecColl[assoTp].begin()->first; 
	    //printf("       Asso track -- (%.3f, %.3f, %.3f, %.3f)\n", assoTrack->p(), assoTrack->px(), assoTrack->py(), assoTrack->pz()); 
	    if(debugPrint) printf("       Asso track -- (%.3f, %.3f, %.3f)\n", assoTrack->p(), assoTrack->eta(), assoTrack->phi()); 
	    //matchFound = (assoTrack.id() == thisTrack.id()); // this DOESN'T work! All tracks have the same ID... 
	    matchFound = (assoTrack == thisTrack); 
	  } 
	} 
	// int motherId(0), grandmaId(0), greatgrandmaId(0); 
	// if( assoTp->genParticles().size()>0 && (*assoTp->genParticle_begin())->numberOfMothers()>0 ) { 
	//   motherId = (*assoTp->genParticle_begin())->mother()->pdgId(); 
	//   if( (*assoTp->genParticle_begin())->mother()->numberOfMothers()>0 ) { 
	//     grandmaId = (*assoTp->genParticle_begin())->mother()->mother()->pdgId(); 
	//     if( (*assoTp->genParticle_begin())->mother()->mother()->numberOfMothers()>0 ) { 
	//       greatgrandmaId = (*assoTp->genParticle_begin())->mother()->mother()->mother()->pdgId(); 
	//     } 
	//   } 
	// } 
	// if(debugPrint) { 
	//   //printf("        Simu track -- (%.3f, %.3f, %.3f, %.3f) <-- %d <-- %d\n", assoTp->p(), assoTp->px(), assoTp->py(), assoTp->pz(), motherId, grandmaId); 
	//   printf("       Simu track -- (%.3f, %.3f, %.3f) <-- %d <-- %d <-- %d\n", assoTp->p(), assoTp->eta(), assoTp->phi(), motherId, grandmaId, greatgrandmaId); 
	// } // end if(debugPrint) 

	// // Is it a signal muon? 
	// matchOnlySignalFound = ( grandmaId==theMotherId || greatgrandmaId==theMotherId || 
	// 			 (motherId==13 && grandmaId==13 && greatgrandmaId==13)   ); 

	matchOnlySignalFound = ( assoTp->genParticles().size()>0 && (*assoTp->genParticle_begin())->isPromptFinalState() ); 
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
      for(unsigned int k=0; k<numberOfMatches; ++k) { 
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
      } // for(...numberOfMatches...)
    } 

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

      for(unsigned int k=0; k<numberOfMatches; ++k) { 
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
      } // for(...numberOfMatches...)
    } 
  } 

  return;
}


// ------------ method called once each job just before starting event loop  ------------
void RecoMuonGeneralAnalyzer::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void RecoMuonGeneralAnalyzer::endJob() { 
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
void RecoMuonGeneralAnalyzer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {
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
void RecoMuonGeneralAnalyzer::endRun(edm::Run const&, edm::EventSetup const&) {
}


// ------------ method called when starting to processes a luminosity block  ------------
void RecoMuonGeneralAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

// ------------ method called when ending the processing of a luminosity block  ------------
void RecoMuonGeneralAnalyzer::endLuminosityBlock(edm::LuminosityBlock const& iLumiBlock, edm::EventSetup const&) {
}


/*
bool RecoMuonGeneralAnalyzer::isGoodME0Muon(edm::ESHandle<ME0Geometry> me0Geom, const reco::ME0Muon& me0muon, double MaxPullX, double MaxDiffX, double MaxPullY, double MaxDiffY, double MaxDiffPhiDir )
{
  ME0Segment thisSegment = me0muon.me0segment();
  
  ME0DetId id = thisSegment.me0DetId();

  auto roll = me0Geom->etaPartition(id); 
  
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
*/ 


// define this as a plug-in
DEFINE_FWK_MODULE(RecoMuonGeneralAnalyzer);
