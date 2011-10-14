/*
 *  See header file for a description of this class.
 *
 *  $Date: 2011/07/20 11:33:20 $
 *  $Revision: 1.5 $
 *  \author G. Cerminara - CERN
 *          D. Trocino   - Northeastern University
 */
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "Tools/MyAnalysisTools/src/ZZllvvAnalyzer.h"
#include "Tools/MyAnalysisTools/interface/Histograms.h"

#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 
#include "DataFormats/PatCandidates/interface/EventHypothesis.h"
#include "DataFormats/PatCandidates/interface/EventHypothesisLooper.h"
#include "DataFormats/Common/interface/MergeableCounter.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include <iostream>
#include <algorithm>
#include "TFile.h"
#include "TString.h"
#include "TNtuple.h"

#include "CMGTools/HtoZZ2l2nu/interface/Utils.h"
#include "CMGTools/HtoZZ2l2nu/interface/ObjectFilters.h"

using namespace std;
using namespace edm;

int nGeneratedEvents;
int nEventsBaseFilterEE;
int nEventsBaseFilterMM;
int nEventsBaseFilterEM;
int nEventsSkimEE;
int nEventsSkimMM;
int nEventsSkimEM;

vector<double> puweight;
float theWeight=1.;

// Analysis cuts
TString allSteps[]={"Analyzed", "Dilepton", "MassWindow", 
		    "JetVeto", "JetMETDeltaPhi", "AntiBTagging", 
		    "LeptonVeto", "RedMETcut"};
const unsigned int nSteps=sizeof(allSteps)/sizeof(TString);
map<string, int> maps;
map<string, int> passedCuts;
//int gShift=3;

TH1F *hPreEventCounter;
TH1F *hEventCounter;
TH1F *hNVertexAll;
TH1F *hNVertexGood;
vector<HistoLept*>       hLeptLead;
vector<HistoLept*>       hLeptSubLead;
vector<HistoLept*>       hLeptThird;
// vector<HistoTrack*>      hSelectedIsoTracks; 
// vector<HistoTrack*>      hSelectedIsoTrackThird; 
vector<HistoDilept*>     hDileptKin;
vector<HistoKin*>        hSelJetKin;
vector<HistoKin*>        hLeadJetKin;
vector<HistoKin*>        hMETKin;
vector<HistoRedMET*>     hRedMetStd;

vector<HistoObjectPair*> hDileptMET;
vector<HistoObjectPair*> hLeadLeptMET;
vector<HistoObjectPair*> hSubleadLeptMET;
vector<HistoObjectPair*> hDileptLeadSelJet;
vector<HistoObjectPair*> hMETLeadSelJet;
vector<HistoObjectPair*> hMETLeadLeptThird;
// vector<HistoObjectPair*> hMETLeadIsoTrackThird;

ReducedMETComputer *redMETComputer_std;

TString varsForLL[]={"crossSection", "branchingRatio", "intLumi", 
		     "run", "lumi", "event", "weight", 
		     "Flavor", 
		     "leadCharge", "subleadCharge", 
		     "leadPt", "subleadPt",
		     "leadEta","subleadEta",
		     "leadPhi","subleadPhi",
		     "dileptPt", 
		     "dileptEta", 
		     "dileptPhi", 
		     "dileptInvMass", 
		     "dileptLeadDeltaPhi",
		     "leptMinusCmCosTheta",
		     "patMet",
		     "patMetPhi",
		     "dileptMetTransMass", 
		     "dileptMetTransMassZ", 
		     "dileptMetTransMassZZ", 
		     "leadMetTransMass", 
		     "leadMetTransMassZ", 
		     "leadMetTransMassZZ", 
		     "subleadMetTransMass", 
		     "subleadMetTransMassZ", 
		     "subleadMetTransMassZZ", 
		     "dilepPROJLong",
		     "dilepPROJPerp",
		     "uncertPROJLong",
		     "uncertPROJPerp",
		     "sumjetPROJLong",
		     "sumjetPROJPerp",
		     "sumAllJetPROJLong",
		     "sumAllJetPROJPerp",
		     "sumOppositeJetPROJLong",
		     "sumOppositeJetPROJPerp",
		     "METPROJLong",
		     "METPROJPerp",
		     "unclPROJLong",
		     "unclPROJPerp",
		     "recoilPROJLong",
		     "recoilPROJPerp",
		     "recoilAllJetPROJLong",
		     "recoilAllJetPROJPerp",
		     "recoilOppositeJetPROJLong",
		     "recoilOppositeJetPROJPerp",
		     "D0redMet",
		     "D0redMet_long",
		     "D0redMet_perp",
		     "CMSredMet",
		     "CMSredMet_long",
		     "CMSredMet_perp",
		     "numberOfLeptons",
		     "thirdLept_flavor",
		     "thirdLept_isGlobalMuon",
		     "thirdLept_isGlobalMuonPT",
		     "thirdLept_isTrackerMuon",
		     "thirdLept_isTrackerMuonLSL",
		     "thirdLept_isTrackerMuonLST",
		     "thirdLept_isTrackerMuonLSAL",
		     "thirdLept_isTrackerMuonLSAT",
		     "thirdLept_isElectronVBTF85",
		     "thirdLept_isElectronVBTF95",
		     "thirdLept_pt",
		     "thirdLept_eta",
		     "thirdLept_phi",
		     "thirdLept_trkChi2",
		     "thirdLept_glbChi2",
		     "thirdLept_corrRelIso",
		     "thirdLept_pixHits",
		     "thirdLept_trkHits",
		     "thirdLept_trkLostInnerHits",
		     "thirdLept_trkLostOuterHits",
		     "thirdLept_muHits",
		     "thirdLept_muMatches",
		     "thirdLept_dxy",
		     "thirdLept_dz",
		     "thirdLept_METleptTransvMass",
		     "fourthLept_flavor",
		     "fourthLept_isGlobalMuon",
		     "fourthLept_isGlobalMuonPT",
		     "fourthLept_isTrackerMuon",
		     "fourthLept_isTrackerMuonLSL",
		     "fourthLept_isTrackerMuonLST",
		     "fourthLept_isTrackerMuonLSAL",
		     "fourthLept_isTrackerMuonLSAT",
		     "fourthLept_isElectronVBTF85",
		     "fourthLept_isElectronVBTF95",
		     "fourthLept_pt",
		     "fourthLept_eta",
		     "fourthLept_phi",
		     "fourthLept_trkChi2",
		     "fourthLept_glbChi2",
		     "fourthLept_corrRelIso",
		     "fourthLept_pixHits",
		     "fourthLept_trkHits",
		     "fourthLept_trkLostInnerHits",
		     "fourthLept_trkLostOuterHits",
		     "fourthLept_muHits",
		     "fourthLept_muMatches",
		     "fourthLept_dxy",
		     "fourthLept_dz",
		     "fourthLept_METleptTransvMass",
		     "CosmicMuon",
		     "Jet_number",
		     "Jet_numberFromPV",
		     "Jet_numberFromPU",
		     "jet_pt_0",
		     "jet_eta_0",
		     "jet_phi_0",
		     "jet_mass_0",
		     "jet_isFromPV_0",
		     "jet_TCHE_0",
		     "jet_SSVHE_0",
		     "jet_JBP_0",
		     "jet_pt_1",
		     "jet_eta_1",
		     "jet_phi_1",
		     "jet_mass_1",
		     "jet_isFromPV_1",
		     "jet_TCHE_1",
		     "jet_SSVHE_1",
		     "jet_JBP_1",
		     "jet_pt_2",
		     "jet_eta_2",
		     "jet_phi_2",
		     "jet_mass_2",
		     "jet_isFromPV_2",
		     "jet_TCHE_2",
		     "jet_SSVHE_2",
		     "jet_JBP_2",
		     "jet_pt_3",
		     "jet_eta_3",
		     "jet_phi_3",
		     "jet_mass_3",
		     "jet_isFromPV_3",
		     "jet_TCHE_3",
		     "jet_SSVHE_3",
		     "jet_JBP_3",
		     "jet_pt_4",
		     "jet_eta_4",
		     "jet_phi_4",
		     "jet_mass_4",
		     "jet_isFromPV_4",
		     "jet_TCHE_4",
		     "jet_SSVHE_4",
		     "jet_JBP_4",
		     "jet_pt_5",
		     "jet_eta_5",
		     "jet_phi_5",
		     "jet_mass_5",
		     "jet_isFromPV_5",
		     "jet_TCHE_5",
		     "jet_SSVHE_5",
		     "jet_JBP_5",
		     "jet_pt_6",
		     "jet_eta_6",
		     "jet_phi_6",
		     "jet_mass_6",
		     "jet_isFromPV_6",
		     "jet_TCHE_6",
		     "jet_SSVHE_6",
		     "jet_JBP_6",
		     "jet_pt_7",
		     "jet_eta_7",
		     "jet_phi_7",
		     "jet_mass_7",
		     "jet_isFromPV_7",
		     "jet_TCHE_7",
		     "jet_SSVHE_7",
		     "jet_JBP_7",
		     "jet_pt_8",
		     "jet_eta_8",
		     "jet_phi_8",
		     "jet_mass_8",
		     "jet_isFromPV_8",
		     "jet_TCHE_8",
		     "jet_SSVHE_8",
		     "jet_JBP_8",
		     "jet_pt_9",
		     "jet_eta_9",
		     "jet_phi_9",
		     "jet_mass_9",
		     "jet_isFromPV_9",
		     "jet_TCHE_9",
		     "jet_SSVHE_9",
		     "jet_JBP_9",
		     "jet_min_deltaPhiJetMET",
		     "jet_min_deltaPhiJetMET_pos",
		     "jet_min_deltaPhiJetMET_isFromPV",
		     "jet_max_CSV",       // CombinedSecondaryVertexBJetTags
		     "jet_max_CSV_pos",
		     "jet_max_CSV_isFromPV",
		     "jet_max_CSVMVA",    // CombinedSecondaryVertexMVABJetTags
		     "jet_max_CSVMVA_pos",
		     "jet_max_CSVMVA_isFromPV",
		     "jet_max_JBP",	  // JetBProbabilityBJetTags
		     "jet_max_JBP_pos",
		     "jet_max_JBP_isFromPV",
		     "jet_max_JP",	  // JetProbabilityBJetTags
		     "jet_max_JP_pos",
		     "jet_max_JP_isFromPV",
		     "jet_max_SSVHE",     // SimpleSecondaryVertexHighEffBJetTags
		     "jet_max_SSVHE_pos",
		     "jet_max_SSVHE_isFromPV",
		     "jet_max_SSVHP",     // SimpleSecondaryVertexHighPurBJetTags
		     "jet_max_SSVHP_pos",
		     "jet_max_SSVHP_isFromPV",
		     "jet_max_SEBPT",     // SoftElectronByPtBJetTags                
		     "jet_max_SEBPT_pos",
		     "jet_max_SEBPT_isFromPV",
		     "jet_max_SEBIP3D",	  // SoftElectronByIP3dBJetTags
		     "jet_max_SEBIP3D_pos",
		     "jet_max_SEBIP3D_isFromPV",
		     "jet_max_SM",	  // SoftMuonBJetTags
		     "jet_max_SM_pos",
		     "jet_max_SM_isFromPV",
		     "jet_max_SMBPT",     // SoftMuonByPtBJetTags                
		     "jet_max_SMBPT_pos",
		     "jet_max_SMBPT_isFromPV",
		     "jet_max_SMBIP3D",	  // SoftMuonByIP3dBJetTags
		     "jet_max_SMBIP3D_pos",
		     "jet_max_SMBIP3D_isFromPV",
		     "jet_max_TCHE",	  // TrackCountingHighEffBJetTags
		     "jet_max_TCHE_pos",
		     "jet_max_TCHE_isFromPV",
		     "jet_max_TCHP",      // TrackCountingHighPurBJetTags        
		     "jet_max_TCHP_pos",
		     "jet_max_TCHP_isFromPV",
		     "totalNVertex",
		     "goodNVertex"};
const unsigned int nVarsOpt=sizeof(varsForLL)/sizeof(TString);

ZZllvvAnalyzer::ZZllvvAnalyzer(const ParameterSet& pSet) : totNEvents(0) 
{
  //theFile = new TFile(pSet.getUntrackedParameter<string>("fileName", "ZZllvvAnalyzer.root").c_str(),"RECREATE");
  theOutFileName = pSet.getUntrackedParameter<string>("fileName", "ZZllvvAnalyzer.root");
  source = pSet.getUntrackedParameter<InputTag>("source");
  isMC = pSet.getUntrackedParameter<bool>("isMC", false);
  theXsect = pSet.getUntrackedParameter<double>("xSection", 1.);
  theBrRatio = pSet.getUntrackedParameter<double>("branchingRatio", 1.);
  theLumi = pSet.getUntrackedParameter<double>("luminosity", 1.);
  debug = pSet.getUntrackedParameter<bool>("debug", false);
  nGeneratedEvents=0;
  nEventsBaseFilterEE=0;
  nEventsBaseFilterMM=0;
  nEventsBaseFilterEM=0;
  nEventsSkimEE=0;
  nEventsSkimMM=0;
  nEventsSkimEM=0;
  flavorCombo = pSet.getUntrackedParameter<int>("FlavorCombination", 0);
  chargeCombo = pSet.getUntrackedParameter<int>("ChargeCombination", -1);
  kRecoilLongWeight = pSet.getUntrackedParameter<double>("RecoilLongWeight", 2.);
  kRecoilPerpWeight = pSet.getUntrackedParameter<double>("RecoilPerpWeight", 2.);
  kSigmaPtLongWeight = pSet.getUntrackedParameter<double>("SigmaPtLongWeight", 2.5);
  kSigmaPtPerpWeight = pSet.getUntrackedParameter<double>("SigmaPtPerpWeight", 2.5);
  kPerpComponentWeight = pSet.getUntrackedParameter<double>("PerpComponentWeight", 1.);
  theRedMETMinCut = pSet.getUntrackedParameter<double>("RedMETMinCut", 50.);
  useAllJets = pSet.getUntrackedParameter<bool>("UseAllJets", true);
  redMETComputer_std   = new ReducedMETComputer(kRecoilLongWeight,
						kRecoilPerpWeight,
						kSigmaPtLongWeight,
						kSigmaPtPerpWeight,
						kPerpComponentWeight);
  vertexSelection =  pSet.getParameter<ParameterSet>("Vertices");
  if(false && isMC) {
    puweight=makePuDistr(pSet);
  }
}

ZZllvvAnalyzer::~ZZllvvAnalyzer(){
  cout << "destructor" << endl;
}



void ZZllvvAnalyzer::beginJob() {
  // Output file
  theFile = new TFile(theOutFileName.c_str(),"RECREATE");

  // Preliminary steps
  hPreEventCounter=new TH1F("PreEventCounter", "Unscaled number of events and more", 12, 0., 12.);
  hPreEventCounter->GetXaxis()->SetBinLabel( 1, "Generated");
  hPreEventCounter->GetXaxis()->SetBinLabel( 2, "BaseFilters-ee");
  hPreEventCounter->GetXaxis()->SetBinLabel( 3, "BaseFilters-mumu");
  hPreEventCounter->GetXaxis()->SetBinLabel( 4, "BaseFilters-emu");
  hPreEventCounter->GetXaxis()->SetBinLabel( 5, "Skim");
  hPreEventCounter->GetXaxis()->SetBinLabel( 6, "Skim-ee");
  hPreEventCounter->GetXaxis()->SetBinLabel( 7, "Skim-mumu");
  hPreEventCounter->GetXaxis()->SetBinLabel( 8, "Skim-emu");
  hPreEventCounter->GetXaxis()->SetBinLabel( 9, "Analyzed");
  hPreEventCounter->GetXaxis()->SetBinLabel(10, "Integrated-luminosity");
  hPreEventCounter->GetXaxis()->SetBinLabel(11, "Cross-section");
  hPreEventCounter->GetXaxis()->SetBinLabel(12, "Branching-ratio");

  // Book the histograms
  initializePlots();

  // Ntuple
  TString varsTmp;
  for(unsigned int z=0; z<nVarsOpt; ++z) {
    if(z) varsTmp+=":";
    varsTmp+=varsForLL[z];
  }
  for(unsigned int y=0; y<nSteps; ++y) {
    varsTmp+=(":"+allSteps[y]);
    passedCuts[allSteps[y].Data()]=1;
  }
  varsTmp+=":AllCutsPassed";
  passedCuts["AllCutsPassed"]=1;
  finalNtpl=new TNtuple("ntuple", "Final results", varsTmp.Data());
}


// Operations
void ZZllvvAnalyzer::beginRun(const Run& run, const EventSetup& eSetup) {
}
  
void ZZllvvAnalyzer::analyze(const Event& event, const EventSetup& eSetup) {

  // Count all the events (the EventCounter histo is filled few lines below...)
  totNEvents++;

  // Reset all flags
  resetFlags();

  if(isMC==event.isRealData()) {
    string decl(isMC ? "MC" : "data");
    string real(event.isRealData() ? "data" : "MC");
    cout << "[ZZllvvAnalyzer] *** Error: sample is declared as " << decl.c_str() 
	 << ", but is " << real.c_str() << endl;
    throw std::exception();
  }

  using reco::Candidate; 
  using reco::CandidatePtr;
  
  if(isMC) {  // MC
    if(false) {
      // Retrieve Pile-Up information
      Handle<vector<PileupSummaryInfo> > PupInfo;
      event.getByLabel("addPileupInfo", PupInfo);
      std::vector<PileupSummaryInfo>::const_iterator pvi;
      unsigned int nPuInt(0), nPuIntIT(0);
      for(pvi = PupInfo->begin(); pvi!=PupInfo->end(); ++pvi) {
	if(debug) {
	  cout << " Pileup Information: bunchXing: " << pvi->getBunchCrossing() << ", nvtx: " << pvi->getPU_NumInteractions() << endl;
	}
	nPuInt+=pvi->getPU_NumInteractions();
	if(pvi->getBunchCrossing()==0) nPuIntIT=pvi->getPU_NumInteractions();
      }

      if(debug) {
	cout << " Total number of PU interaction: " << nPuInt << endl;
	cout << " Total number of PU interaction in time: " << nPuIntIT << endl;
      }

      if(nPuIntIT>=puweight.size()) {
	cout << "[ZZllvvAnalyzer] *** Error: n. PU interactions in time (" << nPuIntIT 
	     << ") greater than the size of weights vector (" << puweight.size() 
	     << ")! " << endl;
	throw std::exception();
      }
      theWeight=puweight[nPuIntIT];
    }
    else {
      edm::Handle<float> puWeightHandle;
      event.getByLabel("puWeights", "puWeight", puWeightHandle);
      if(puWeightHandle.isValid()) theWeight=*(puWeightHandle.product());
    }
  } // end if(isMC)
  else {  // Real data
    theWeight=1.;   // already 1. at initialization, but you never know...
  }
  if(debug) cout << "event weight = "  << theWeight << endl;

  // Count all the events (EventCounter)
  hEventCounter->Fill(0., theWeight);

  // Pre-select vertices
  Handle<reco::VertexCollection> offlinePrimVertices;
  event.getByLabel("offlinePrimaryVertices", offlinePrimVertices);  
  unsigned int totNvtx=offlinePrimVertices->size(); //TOTAL NUMBER OF VERTICES
  if(debug) cout << "Total # of vertices: " << totNvtx << endl;;
  std::vector<reco::VertexRef> selVertices = vertex::filter(offlinePrimVertices,vertexSelection);
  unsigned int goodNvtx=selVertices.size(); //GOOD NUMBER OF VERTICES
  if(debug) cout << "# of good vertices: " << goodNvtx << endl;;
  float totNvtx_float = 1.0*totNvtx;
  float goodNvtx_float = 1.0*goodNvtx;
  

  // Average energy density
  edm::Handle< double > rhoH;
  event.getByLabel(InputTag("kt6PFJetsPFlow:rho"), rhoH);
  double rho=(*rhoH);
  if(debug) cout << "rho = "  << rho << endl;

  // Retrieve the Event Hypothesis
  Handle<vector<pat::EventHypothesis> > hyps;
  edm::InputTag selEvtTag(source.label()+":selectedEvent");
  event.getByLabel(selEvtTag, hyps);

  if(hyps->size()==0) return;
  const pat::EventHypothesis &h = (*hyps)[0];

  // Retrieve the Selection Path
  Handle<std::vector<int> > selectionInfo;
  edm::InputTag selInfoTag(source.label()+":selectionInfo");
  event.getByLabel(selInfoTag, selectionInfo);

  // Retrieve the selected Vertex
  Handle<reco::VertexCollection> selectedVertex;
  edm::InputTag selVtxTag(source.label()+":selectedVertices");
  event.getByLabel(selVtxTag, selectedVertex);
  const reco::Vertex *vtx = &(*selectedVertex)[0];

  // Dump selected event 
  //   SelectionInfo [0] -> path: {UNKNOWN=0,MUMU=1,EE=2,EMU=3};
  //   SelectionInfo [1] -> step: 
  if(debug) {
    cout << "Retrieve an event hypothesis selected at step=" << (*selectionInfo)[1]
	 << " for path " << (*selectionInfo)[0] 
	 << " with " << selectedVertex->size() << " vertices" << endl;
  }


  // Get the MET
  const pat::MET *met=h.getAs<pat::MET>("met");
  if(debug) cout << "\t MET:" << met->pt() << ";" << met->phi() << endl;


  // Check if there are Dileptons
  if((*selectionInfo)[0]==0) {
    if(debug) cout << "\t No dilepton has been selected" << endl;
    return;
  }

  // Skip unwanted flavor combinations
  switch(flavorCombo) {
  case 0:
    if( (*selectionInfo)[0]!=1 && (*selectionInfo)[0]!=2 ) return;
    break;
  case 1:
  case 2:
  case 3:
    if( (*selectionInfo)[0]!=flavorCombo ) return;
    break;
  case 4:
    if(debug) cout << "[ZZllvvAnalyzer] Any flavor combination will be considered." << endl;
    break;
  default:
    cout << "[ZZllvvAnalyzer] *** Error: wrong flavor combination selected!" << endl;
    throw std::exception(); 
  }

  // Get the two best lepton candidates... 
  CandidatePtr lep1=h["leg1"];
  CandidatePtr lep2=h["leg2"];

  if(lep1.get()==0 || lep2.get()==0) {
    cout << "[ZZllvvAnalyzer] *** Warning: no dilepton candidate!" << endl;
    return;
  }

  // Let lep1 be the leading lepton
  if(lep1->pt()<lep2->pt()) {
    CandidatePtr lepTmp=lep1;
    lep1=lep2;
    lep2=lepTmp;
  }

  // ... try casting to Muons...
  const pat::Muon *muonLead    = dynamic_cast<const pat::Muon *>(lep1.get());
  const pat::Muon *muonSubLead = dynamic_cast<const pat::Muon *>(lep2.get());

  // ... and to Electrons
  const pat::Electron *electronLead    = dynamic_cast<const pat::Electron *>(lep1.get());
  const pat::Electron *electronSubLead = dynamic_cast<const pat::Electron *>(lep2.get());

  if( (muonLead==0 && electronLead==0) || (muonSubLead==0 && electronSubLead==0) ) {
    cout << "[ZZllvvAnalyzer] *** Error: dilepton found, "
	 << "but something went wrong with the dynamc_cast!" << endl;
    cout << " Leading lepton:    muon=" << muonLead << ", electron=" << electronLead << endl;
    cout << " Subleading lepton: muon=" << muonSubLead << ", electron=" << electronSubLead << endl;
    throw std::exception();
  }

  // TEMPORARY!!!!
  double muonEffectiveArea=0.112;
  double electronEffectiveArea=0.24;
  //double muonEffectiveArea=0.;
  //double electronEffectiveArea=0.;
  double Aeff1=(muonLead ? muonEffectiveArea : electronEffectiveArea);
  double Aeff2=(muonSubLead ? muonEffectiveArea : electronEffectiveArea);
  double Aeff3=0., Aeff4=0.;

  double relIso1=lepton::getLeptonIso(lep1, 1., rho*Aeff1)[lepton::REL_ISO];
  double relIso2=lepton::getLeptonIso(lep2, 1., rho*Aeff2)[lepton::REL_ISO];
  
  if(debug) {
    cout << "\t dilepton leg1: type=" << (muonLead ? "MUON" : (electronLead ? "ELECTRON" : "UNKNOWN")) 
	 << "; charge=" << lep1->charge() << "; pt=" << lep1->pt() 
	 << "; eta=" << lep1->eta() << "; phi=" << lep1->phi() << "; relIsoCorr=" << relIso1 << endl;
    cout << "\t          leg2: type=" << (muonSubLead ? "MUON" : (electronSubLead ? "ELECTRON" : "UNKNOWN")) 
	 << "; charge=" << lep2->charge() << "; pt=" << lep2->pt() 
	 << "; eta=" << lep2->eta() << "; phi=" << lep2->phi() << "; relIsoCorr=" << relIso2 << endl;
  }

  //   // Skip unwanted flavor combinations
  //   if( flavorCombo==0 && (muonLead==0 || muonSubLead==0) && (electronLead==0 || electronSubLead==0) ) 
  //     return; 
  //   if( flavorCombo==1 && (muonLead==0 || muonSubLead==0) ) 
  //     return; 
  //   if( flavorCombo==2 && (electronLead==0 || electronSubLead==0) )
  //     return;
  //   if( flavorCombo==3 && (muonLead==0 || electronSubLead==0) && (electronLead==0 || muonSubLead==0) ) 
  //     return; 

  // =====================================================================
  //  ? ask for the dilepton to exist (right flavour, charge)
  // =====================================================================

  typedef reco::Candidate::LorentzVector LorentzVector;
  LorentzVector leadLeptMom = lep1->p4();
  LorentzVector subLeadLeptMom = lep2->p4();
  LorentzVector diLeptonMom = leadLeptMom + subLeadLeptMom;

  // FLAG COSMIC MUONS (and "cosmic" electrons)
  float lep1x = leadLeptMom.Px();
  float lep1y = leadLeptMom.Py();
  float lep1z = leadLeptMom.Pz();
  float lep2x = subLeadLeptMom.Px();
  float lep2y = subLeadLeptMom.Py();
  float lep2z = subLeadLeptMom.Pz();
  float cosang = (lep1x*lep2x+lep1y*lep2y+lep1z*lep2z)/(sqrt(lep1x*lep1x+lep1y*lep1y+lep1z*lep1z)*sqrt(lep2x*lep2x+lep2y*lep2y+lep2z*lep2z));
  float CosMu = 0.0;
  if(cosang <= -0.9998){
    CosMu = 1.0;
  }

  // =====================================================================
  //  ? apply mass window on the Z mass
  // =====================================================================


  // Get the leading "extra lepton" ("third lepton")
  CandidatePtr leadExtraLept;
  float leadExtraLeptId=0.;
  string leadExtraLeptType;
  float leadExtraLeptIsGM=-1.;
  float leadExtraLeptIsGMPT=-1.;
  float leadExtraLeptIsTM=-1.;
  float leadExtraLeptIsTMLSL=-1.;
  float leadExtraLeptIsTMLST=-1.;
  float leadExtraLeptIsTMLSAL=-1.;
  float leadExtraLeptIsTMLSAT=-1.;
  float leadExtraLeptIdVBTF85=-1.;
  float leadExtraLeptIdVBTF95=-1.;
  float leadExtraLeptPt=0.;
  float leadExtraLeptEta=-9999.;
  float leadExtraLeptPhi=-9999.;
  float leadExtraLeptTrkChi2=-1.;
  float leadExtraLeptGlbChi2=-1.;
  float leadExtraLeptRelIso=-1.;
  float leadExtraLeptNPixHits=-1.;
  float leadExtraLeptNTrkHits=-1.;
  float leadExtraLeptNTrkLostInnHits=-1.;
  float leadExtraLeptNTrkLostOutHits=-1.;
  float leadExtraLeptNMuHits=-1.;
  float leadExtraLeptNMuMatch=-1.;
  float leadExtraLeptDxy=-1.;
  float leadExtraLeptDz=-1.;
  float leadExtraLeptMetTransMass=-1.;

  // Get the subleading "extra lepton" ("fourth lepton")
  CandidatePtr subleadExtraLept;
  float subleadExtraLeptId=0.;
  string subleadExtraLeptType;
  float subleadExtraLeptIsGM=-1.;
  float subleadExtraLeptIsGMPT=-1.;
  float subleadExtraLeptIsTM=-1.;
  float subleadExtraLeptIsTMLSL=-1.;
  float subleadExtraLeptIsTMLST=-1.;
  float subleadExtraLeptIsTMLSAL=-1.;
  float subleadExtraLeptIsTMLSAT=-1.;
  float subleadExtraLeptIdVBTF85=-1.;
  float subleadExtraLeptIdVBTF95=-1.;
  float subleadExtraLeptPt=0.;
  float subleadExtraLeptEta=-9999.;
  float subleadExtraLeptPhi=-9999.;
  float subleadExtraLeptTrkChi2=-1.;
  float subleadExtraLeptGlbChi2=-1.;
  float subleadExtraLeptRelIso=-1.;
  float subleadExtraLeptNPixHits=-1.;
  float subleadExtraLeptNTrkHits=-1.;
  float subleadExtraLeptNTrkLostInnHits=-1.;
  float subleadExtraLeptNTrkLostOutHits=-1.;
  float subleadExtraLeptNMuHits=-1.;
  float subleadExtraLeptNMuMatch=-1.;
  float subleadExtraLeptDxy=-1.;
  float subleadExtraLeptDz=-1.;
  float subleadExtraLeptMetTransMass=-1.;


  // Get the collection of selected, isolated muons
  vector<const pat::Muon *> allSelMuons;
  for(pat::eventhypothesis::Looper<pat::Muon> muon=h.loopAs<pat::Muon>("muon"); muon; ++muon) {
    if(debug) 
      cout << "\t muon: " << muon->pt() << "; " << muon->eta() << "; " << muon->phi() << endl;
    const pat::Muon *tmpMu=dynamic_cast<const pat::Muon *>(muon.get());
    allSelMuons.push_back(tmpMu);
    if(muon->pt()>leadExtraLeptPt) {

      if(leadExtraLeptPt>0.) {
	subleadExtraLept=leadExtraLept;
	subleadExtraLeptPt=leadExtraLeptPt;
	subleadExtraLeptEta=leadExtraLeptEta;
	subleadExtraLeptPhi=leadExtraLeptPhi;
	subleadExtraLeptType=leadExtraLeptType;
	subleadExtraLeptId=leadExtraLeptId;
      }

      leadExtraLept=muon.ref();
      leadExtraLeptPt=muon->pt();
      leadExtraLeptEta=muon->eta();
      leadExtraLeptPhi=muon->phi();
      leadExtraLeptType="muon";
      leadExtraLeptId=13.*muon->charge();
    }
  }
  if(debug) 
    cout << "# of selected, isolated muons: " << allSelMuons.size() << endl;

  // Get the collection of selected, isolated electrons
  vector<const pat::Electron *> allSelElectrons;
  for(pat::eventhypothesis::Looper<pat::Electron> ele=h.loopAs<pat::Electron>("electron"); ele; ++ele) {
    if(debug) 
      cout << "\t electron: " << ele->pt() << "; " << ele->eta() << "; " << ele->phi() << endl;
    const pat::Electron *tmpEle=dynamic_cast<const pat::Electron *>(ele.get());
    allSelElectrons.push_back(tmpEle);
 
    if(ele->pt()>leadExtraLeptPt) {

      if(leadExtraLeptPt>0.) {
	subleadExtraLept=leadExtraLept;
	subleadExtraLeptPt=leadExtraLeptPt;
	subleadExtraLeptEta=leadExtraLeptEta;
	subleadExtraLeptPhi=leadExtraLeptPhi;
	subleadExtraLeptType=leadExtraLeptType;
	subleadExtraLeptId=leadExtraLeptId;
      }

      leadExtraLept=ele.ref();
      leadExtraLeptPt=ele->pt();
      leadExtraLeptEta=ele->eta();
      leadExtraLeptPhi=ele->phi();
      leadExtraLeptType="electron";
      leadExtraLeptId=11.*ele->charge();
    }
  }
  if(debug) {
    cout << "# of selected, isolated electrons: " << allSelElectrons.size() << endl;
  }

  const pat::Muon *mu3=dynamic_cast<const pat::Muon *>(leadExtraLept.get());
  const pat::Electron *el3=dynamic_cast<const pat::Electron *>(leadExtraLept.get());

  const pat::Muon *mu4=dynamic_cast<const pat::Muon *>(subleadExtraLept.get());
  const pat::Electron *el4=dynamic_cast<const pat::Electron *>(subleadExtraLept.get());

  if(leadExtraLept.get()!=0) {
    if( abs(abs(leadExtraLeptId)-13)<0.1 ) {    // if( abs(leadExtraLeptId)==13 ) {
      if( mu3==0 ) {
	cout << "[ZZllvvAnalyzer] *** Error: Leading extra lepton is a muon, "
	     << "but something went wrong with the dynamic_cast!" << endl;
	throw std::exception();
      }
      leadExtraLeptIsGM    =float(mu3->isGlobalMuon());
      leadExtraLeptIsGMPT  =float(mu3->muonID("GlobalMuonPromptTight"));
      leadExtraLeptIsTM    =float(mu3->isTrackerMuon());
      leadExtraLeptIsTMLSL =float(mu3->muonID("TMLastStationLoose"));
      leadExtraLeptIsTMLST =float(mu3->muonID("TMLastStationTight"));
      leadExtraLeptIsTMLSAL=float(mu3->muonID("TMLastStationAngLoose"));
      leadExtraLeptIsTMLSAT=float(mu3->muonID("TMLastStationAngTight"));
      leadExtraLeptTrkChi2=mu3->innerTrack()->normalizedChi2();
      if(mu3->isGlobalMuon()) leadExtraLeptGlbChi2=mu3->globalTrack()->normalizedChi2();
      leadExtraLeptNPixHits=mu3->innerTrack()->hitPattern().numberOfValidPixelHits();
      leadExtraLeptNTrkHits=mu3->innerTrack()->hitPattern().numberOfValidTrackerHits();
      leadExtraLeptNTrkLostInnHits=mu3->innerTrack()->trackerExpectedHitsInner().numberOfLostHits();
      leadExtraLeptNTrkLostOutHits=mu3->innerTrack()->trackerExpectedHitsOuter().numberOfLostHits();
      if(mu3->isGlobalMuon()) leadExtraLeptNMuHits=mu3->globalTrack()->hitPattern().numberOfValidMuonHits();
      leadExtraLeptNMuMatch=mu3->numberOfMatches();
      leadExtraLeptDxy=fabs( mu3->innerTrack()->dxy(vtx->position()) );
      leadExtraLeptDz=fabs( mu3->innerTrack()->dz(vtx->position()) );
      Aeff3=muonEffectiveArea;
    }
    else if( abs(abs(leadExtraLeptId)-11)<0.1 ) {    // else if( abs(leadExtraLeptId)==11 ) {
      if( el3==0 ) {
	cout << "[ZZllvvAnalyzer] *** Error: Leading extra lepton is an electron, "
	     << "but something went wrong with the dynamic_cast!" << endl;
	throw std::exception();
      }
      leadExtraLeptIdVBTF85=float( int(el3->electronID("eidVBTF85")) & 0x1 );
      leadExtraLeptIdVBTF95=float( int(el3->electronID("eidVBTF95")) & 0x1 );
      leadExtraLeptTrkChi2=el3->gsfTrack()->normalizedChi2();
      leadExtraLeptNPixHits=el3->gsfTrack()->hitPattern().numberOfValidPixelHits();
      leadExtraLeptNTrkHits=el3->gsfTrack()->hitPattern().numberOfValidTrackerHits();
      leadExtraLeptNTrkLostInnHits=el3->gsfTrack()->trackerExpectedHitsInner().numberOfLostHits();
      leadExtraLeptNTrkLostOutHits=el3->gsfTrack()->trackerExpectedHitsOuter().numberOfLostHits();
      leadExtraLeptDxy=fabs( el3->gsfTrack()->dxy(vtx->position()) );
      leadExtraLeptDz=fabs( el3->gsfTrack()->dz(vtx->position()) );
      Aeff3=electronEffectiveArea;
    }
    leadExtraLeptRelIso=lepton::getLeptonIso(leadExtraLept, 1., rho*Aeff3)[lepton::REL_ISO];
    double leadMetPtTot=sqrt((leadExtraLept->momentum()+met->momentum()).perp2());
    double leadMetEtTot=leadExtraLeptPt+met->pt();
    leadExtraLeptMetTransMass=sqrt(leadMetEtTot*leadMetEtTot - leadMetPtTot*leadMetPtTot);

    if(debug) {
      cout << " Leading extra lepton: " << leadExtraLeptType.c_str() 
	   << ", pt=" << leadExtraLeptPt << "; relIsoCorr=" << leadExtraLeptRelIso << endl;
    }
  }

  if(subleadExtraLept.get()!=0) {
    if( abs(abs(subleadExtraLeptId)-13)<0.1 ) {    // if( abs(subleadExtraLeptId)==13 ) {
      if( mu4==0 ) {
	cout << "[ZZllvvAnalyzer] *** Error: Subleading extra lepton is a muon, "
	     << "but something went wrong with the dynamic_cast!" << endl;
	throw std::exception();
      }
      subleadExtraLeptIsGM    =float(mu4->isGlobalMuon());
      subleadExtraLeptIsGMPT  =float(mu4->muonID("GlobalMuonPromptTight"));
      subleadExtraLeptIsTM    =float(mu4->isTrackerMuon());
      subleadExtraLeptIsTMLSL =float(mu4->muonID("TMLastStationLoose"));
      subleadExtraLeptIsTMLST =float(mu4->muonID("TMLastStationTight"));
      subleadExtraLeptIsTMLSAL=float(mu4->muonID("TMLastStationAngLoose"));
      subleadExtraLeptIsTMLSAT=float(mu4->muonID("TMLastStationAngTight"));
      subleadExtraLeptTrkChi2=mu4->innerTrack()->normalizedChi2();
      if(mu4->isGlobalMuon()) subleadExtraLeptGlbChi2=mu4->globalTrack()->normalizedChi2();
      subleadExtraLeptNPixHits=mu4->innerTrack()->hitPattern().numberOfValidPixelHits();
      subleadExtraLeptNTrkHits=mu4->innerTrack()->hitPattern().numberOfValidTrackerHits();
      subleadExtraLeptNTrkLostInnHits=mu4->innerTrack()->trackerExpectedHitsInner().numberOfLostHits();
      subleadExtraLeptNTrkLostOutHits=mu4->innerTrack()->trackerExpectedHitsOuter().numberOfLostHits();
      if(mu4->isGlobalMuon()) subleadExtraLeptNMuHits=mu4->globalTrack()->hitPattern().numberOfValidMuonHits();
      subleadExtraLeptNMuMatch=mu4->numberOfMatches();
      subleadExtraLeptDxy=fabs( mu4->innerTrack()->dxy(vtx->position()) );
      subleadExtraLeptDz=fabs( mu4->innerTrack()->dz(vtx->position()) );
      Aeff4=muonEffectiveArea;
    }
    else if( abs(abs(subleadExtraLeptId)-11)<0.1 ) {    // else if( abs(subleadExtraLeptId)==11 ) {
      if( el4==0 ) {
	cout << "[ZZllvvAnalyzer] *** Error: Subleading extra lepton is an electron, "
	     << "but something went wrong with the dynamic_cast!" << endl;
	throw std::exception();
      }
      subleadExtraLeptIdVBTF85=float( int(el4->electronID("eidVBTF85")) & 0x1 );
      subleadExtraLeptIdVBTF95=float( int(el4->electronID("eidVBTF95")) & 0x1 );
      subleadExtraLeptTrkChi2=el4->gsfTrack()->normalizedChi2();
      subleadExtraLeptNPixHits=el4->gsfTrack()->hitPattern().numberOfValidPixelHits();
      subleadExtraLeptNTrkHits=el4->gsfTrack()->hitPattern().numberOfValidTrackerHits();
      subleadExtraLeptNTrkLostInnHits=el4->gsfTrack()->trackerExpectedHitsInner().numberOfLostHits();
      subleadExtraLeptNTrkLostOutHits=el4->gsfTrack()->trackerExpectedHitsOuter().numberOfLostHits();
      subleadExtraLeptDxy=fabs( el4->gsfTrack()->dxy(vtx->position()) );
      subleadExtraLeptDz=fabs( el4->gsfTrack()->dz(vtx->position()) );
      Aeff4=electronEffectiveArea;
    }
    subleadExtraLeptRelIso=lepton::getLeptonIso(subleadExtraLept, 1., rho*Aeff4)[lepton::REL_ISO];
    double subleadMetPtTot=sqrt((leadExtraLept->momentum()+met->momentum()).perp2());
    double subleadMetEtTot=leadExtraLeptPt+met->pt();
    subleadExtraLeptMetTransMass=sqrt(subleadMetEtTot*subleadMetEtTot - subleadMetPtTot*subleadMetPtTot);

    if(debug) {
      cout << " Subleading extra lepton: " << subleadExtraLeptType.c_str() 
	   << ", pt=" << subleadExtraLeptPt << "; relIsoCorr=" << subleadExtraLeptRelIso << endl;
    }
  }

  int totNumberOfSelElectrons=allSelElectrons.size();
  int totNumberOfSelMuons=allSelMuons.size();
  int totNumberOfSelLeptons=totNumberOfSelElectrons+totNumberOfSelMuons;

  // =====================================================================
  //  ? veto on the third lepton
  // =====================================================================

  // Pt, eta, phi, mass of the highest-pt jets
  unsigned int numberOfJets(0), numberOfJetsFromPV(0), numberOfJetsFromPU(0);
  const unsigned int nMaxJets=10;
  vector<float> jet_pt_vect(nMaxJets, -1.);
  vector<float> jet_eta_vect(nMaxJets, 9999.);
  vector<float> jet_phi_vect(nMaxJets, 4.);
  vector<float> jet_mass_vect(nMaxJets, -1.);
  vector<float> jet_isFromPV_vect(nMaxJets, -1.);
  vector<float> jet_TCHE_vect(nMaxJets, -10.);
  vector<float> jet_SSVHE_vect(nMaxJets, -10.);
  vector<float> jet_JBP_vect(nMaxJets, -10.);

  // Min DeltaPhi between a jet and MET
  float jet_min_deltaPhiJetMET=4.;
  float jet_min_deltaPhiJetMET_pos=-1.;
  float jet_min_deltaPhiJetMET_isFromPV=-1.;

  // Max value of several b-tagging discriminators
  float jet_max_CSV=-10.;
  float jet_max_CSV_pos=-1.;
  float jet_max_CSV_isFromPV=-1.;
  float jet_max_CSVMVA=-10.;
  float jet_max_CSVMVA_pos=-1.;
  float jet_max_CSVMVA_isFromPV=-1.;
  float jet_max_JBP=-10.;
  float jet_max_JBP_pos=-1.;
  float jet_max_JBP_isFromPV=-1.;
  float jet_max_JP=-10.;
  float jet_max_JP_pos=-1.;
  float jet_max_JP_isFromPV=-1.;
  float jet_max_SSVHE=-10.;
  float jet_max_SSVHE_pos=-1.;
  float jet_max_SSVHE_isFromPV=-1.;
  float jet_max_SSVHP=-10.;
  float jet_max_SSVHP_pos=-1.;
  float jet_max_SSVHP_isFromPV=-1.;
  float jet_max_SEBPT=-10.;
  float jet_max_SEBPT_pos=-1.;
  float jet_max_SEBPT_isFromPV=-1.;
  float jet_max_SEBIP3D=-10.;
  float jet_max_SEBIP3D_pos=-1.;
  float jet_max_SEBIP3D_isFromPV=-1.;
  float jet_max_SM=-10.;
  float jet_max_SM_pos=-1.;
  float jet_max_SM_isFromPV=-1.;
  float jet_max_SMBPT=-10.;
  float jet_max_SMBPT_pos=-1.;
  float jet_max_SMBPT_isFromPV=-1.;
  float jet_max_SMBIP3D=-10.;
  float jet_max_SMBIP3D_pos=-1.;
  float jet_max_SMBIP3D_isFromPV=-1.;
  float jet_max_TCHE=-10.;
  float jet_max_TCHE_pos=-1.;
  float jet_max_TCHE_isFromPV=-1.;
  float jet_max_TCHP=-10.;
  float jet_max_TCHP_pos=-1.;
  float jet_max_TCHP_isFromPV=-1.;

  // Get the Jet momenta
  vector<const pat::Jet *> moreSelJets;
  vector<LorentzVector> moreSelJetMomenta;
  vector<bool> moreSelJetIsFromPV;

  for(pat::eventhypothesis::Looper<const pat::Jet> jet=h.loopAs<const pat::Jet>("jet"); jet; ++jet, ++numberOfJets, ++numberOfJetsFromPV) {
    const pat::Jet *tmpJet=dynamic_cast<const pat::Jet *>(jet.get());
    moreSelJets.push_back(tmpJet);
    moreSelJetMomenta.push_back(jet->p4());
    moreSelJetIsFromPV.push_back(true);
  }

  for(pat::eventhypothesis::Looper<const pat::Jet> jet=h.loopAs<const pat::Jet>("pujet"); jet; ++jet, ++numberOfJets, ++numberOfJetsFromPU) {
    const pat::Jet *tmpJet=dynamic_cast<const pat::Jet *>(jet.get());
    moreSelJets.push_back(tmpJet);
    moreSelJetMomenta.push_back(jet->p4());
    moreSelJetIsFromPV.push_back(false);
  }

  // Sort by pt
  for(unsigned int ii=0; ii<moreSelJetMomenta.size(); ++ii) {
    for(unsigned int jj=ii+1; jj<moreSelJetMomenta.size(); ++jj) {
      if(moreSelJetMomenta[jj].pt()>moreSelJetMomenta[ii].pt()) {
	const pat::Jet *tmpJet=moreSelJets[ii];
	moreSelJets[ii]=moreSelJets[jj];
	moreSelJets[jj]=tmpJet;

	LorentzVector tmpJetVec=moreSelJetMomenta[ii];
	moreSelJetMomenta[ii]=moreSelJetMomenta[jj];
	moreSelJetMomenta[jj]=tmpJetVec;

	bool tmpJetIsFromPV=moreSelJetIsFromPV[ii];
	moreSelJetIsFromPV[ii]=moreSelJetIsFromPV[jj];
	moreSelJetIsFromPV[jj]=tmpJetIsFromPV;
      } // end if
    } // end for(...jj...)

    if(debug) cout << "\t jet: " << moreSelJetMomenta[ii].pt() << ";" << moreSelJetMomenta[ii].eta() << ";" << moreSelJetMomenta[ii].phi() << std::endl;

    if(ii<nMaxJets) {
      jet_pt_vect[ii]=moreSelJets[ii]->pt();
      jet_eta_vect[ii]=moreSelJets[ii]->eta();
      jet_phi_vect[ii]=moreSelJets[ii]->phi();
      jet_mass_vect[ii]=moreSelJets[ii]->mass();
      jet_isFromPV_vect[ii]=(float)moreSelJetIsFromPV[ii];
      jet_TCHE_vect[ii]=moreSelJets[ii]->bDiscriminator("trackCountingHighEffBJetTags");
      jet_SSVHE_vect[ii]=moreSelJets[ii]->bDiscriminator("simpleSecondaryVertexHighEffBJetTags");
      jet_JBP_vect[ii]=moreSelJets[ii]->bDiscriminator("jetBProbabilityBJetTags");
    }

    if( deltaPhi(met->phi(), moreSelJetMomenta[ii].phi()) < jet_min_deltaPhiJetMET ) {
      jet_min_deltaPhiJetMET=deltaPhi(met->phi(), moreSelJetMomenta[ii].phi());
      jet_min_deltaPhiJetMET_pos=(float)ii;
      jet_min_deltaPhiJetMET_isFromPV=(float)moreSelJetIsFromPV[ii];
    }
    if( moreSelJets[ii]->bDiscriminator("combinedSecondaryVertexBJetTags")>jet_max_CSV ) {
      jet_max_CSV=moreSelJets[ii]->bDiscriminator("combinedSecondaryVertexBJetTags");
      jet_max_CSV_pos=(float)ii;
      jet_max_CSV_isFromPV=(float)moreSelJetIsFromPV[ii];
    }
    if( moreSelJets[ii]->bDiscriminator("combinedSecondaryVertexMVABJetTags")>jet_max_CSVMVA ) {
      jet_max_CSVMVA=moreSelJets[ii]->bDiscriminator("combinedSecondaryVertexMVABJetTags");
      jet_max_CSVMVA_pos=(float)ii;
      jet_max_CSVMVA_isFromPV=(float)moreSelJetIsFromPV[ii];
    }
    if( moreSelJets[ii]->bDiscriminator("jetBProbabilityBJetTags")>jet_max_JBP ) {
      jet_max_JBP=moreSelJets[ii]->bDiscriminator("jetBProbabilityBJetTags");
      jet_max_JBP_pos=(float)ii;
      jet_max_JBP_isFromPV=(float)moreSelJetIsFromPV[ii];
    }
    if( moreSelJets[ii]->bDiscriminator("jetProbabilityBJetTags")>jet_max_JP ) {
      jet_max_JP=moreSelJets[ii]->bDiscriminator("jetProbabilityBJetTags");
      jet_max_JP_pos=(float)ii;
      jet_max_JP_isFromPV=(float)moreSelJetIsFromPV[ii];
    }
    if( moreSelJets[ii]->bDiscriminator("simpleSecondaryVertexHighEffBJetTags")>jet_max_SSVHE ) {
      jet_max_SSVHE=moreSelJets[ii]->bDiscriminator("simpleSecondaryVertexHighEffBJetTags");
      jet_max_SSVHE_pos=(float)ii;
      jet_max_SSVHE_isFromPV=(float)moreSelJetIsFromPV[ii];
    }
    if( moreSelJets[ii]->bDiscriminator("simpleSecondaryVertexHighPurBJetTags")>jet_max_SSVHP ) {
      jet_max_SSVHP=moreSelJets[ii]->bDiscriminator("simpleSecondaryVertexHighPurBJetTags");
      jet_max_SSVHP_pos=(float)ii;
      jet_max_SSVHP_isFromPV=(float)moreSelJetIsFromPV[ii];
    }
    if( moreSelJets[ii]->bDiscriminator("softElectronByPtBJetTags")>jet_max_SEBPT ) {
      jet_max_SEBPT=moreSelJets[ii]->bDiscriminator("softElectronByPtBJetTags");
      jet_max_SEBPT_pos=(float)ii;
      jet_max_SEBPT_isFromPV=(float)moreSelJetIsFromPV[ii];
    }
    if( moreSelJets[ii]->bDiscriminator("softElectronByIP3dBJetTags")>jet_max_SEBIP3D ) {
      jet_max_SEBIP3D=moreSelJets[ii]->bDiscriminator("softElectronByIP3dBJetTags");
      jet_max_SEBIP3D_pos=(float)ii;
      jet_max_SEBIP3D_isFromPV=(float)moreSelJetIsFromPV[ii];
    }
    if( moreSelJets[ii]->bDiscriminator("softMuonBJetTags")>jet_max_SM ) {
      jet_max_SM=moreSelJets[ii]->bDiscriminator("softMuonBJetTags");
      jet_max_SM_pos=(float)ii;
      jet_max_SM_isFromPV=(float)moreSelJetIsFromPV[ii];
    }
    if( moreSelJets[ii]->bDiscriminator("softMuonByPtBJetTags")>jet_max_SMBPT ) {
      jet_max_SMBPT=moreSelJets[ii]->bDiscriminator("softMuonByPtBJetTags");
      jet_max_SMBPT_pos=(float)ii;
      jet_max_SMBPT_isFromPV=(float)moreSelJetIsFromPV[ii];
    }
    if( moreSelJets[ii]->bDiscriminator("softMuonByIP3dBJetTags")>jet_max_SMBIP3D ) {
      jet_max_SMBIP3D=moreSelJets[ii]->bDiscriminator("softMuonByIP3dBJetTags");
      jet_max_SMBIP3D_pos=(float)ii;
      jet_max_SMBIP3D_isFromPV=(float)moreSelJetIsFromPV[ii];
    }
    if( moreSelJets[ii]->bDiscriminator("trackCountingHighEffBJetTags")>jet_max_TCHE ) {
      jet_max_TCHE=moreSelJets[ii]->bDiscriminator("trackCountingHighEffBJetTags");
      jet_max_TCHE_pos=(float)ii;
      jet_max_TCHE_isFromPV=(float)moreSelJetIsFromPV[ii];
    }
    if( moreSelJets[ii]->bDiscriminator("trackCountingHighPurBJetTags")>jet_max_TCHP ) {
      jet_max_TCHP=moreSelJets[ii]->bDiscriminator("trackCountingHighPurBJetTags");
      jet_max_TCHP_pos=(float)ii;
      jet_max_TCHP_isFromPV=(float)moreSelJetIsFromPV[ii];
    }
  } // end for(...ii...)

  // // Get the "leading" jet
  // int leadJet=-1, leadJetCnt=0; //number of jets
  // float leadExtraJet_Pt=0.;
  // float leadExtraJet_Eta=-9999.0;
  // float leadExtraJet_Phi=6.5;
  // float leadExtraJet_Mass=0.;
  // float deltaphimetjet=3.5;
  // float deltaphimetjet_Pt=0.;
  // float deltaphimetjet_Eta=-9999.0;
  // float deltaphimetjet_Phi=6.5;
  // float deltaphimetjet_Mass=0.;
  // int deltaphimetjet_pos=-1;
  // float lead_metjet_equal=0.0;

  // // Get the Jet momenta
  // vector<CandidatePtr> allSelJets;
  // vector<LorentzVector> moreSelJetMomenta;
  // for(pat::eventhypothesis::Looper<pat::Jet> jet=h.loopAs<pat::Jet>("jet"); jet; ++jet, ++leadJetCnt) {
  //   double tmpJetPt=jet->pt();
  //   if(debug) cout << "\t jet: " << tmpJetPt << ";" << jet->eta() << ";" << jet->phi() << std::endl;
  //   moreSelJetMomenta.push_back(jet->p4());
  //   if(tmpJetPt>leadExtraJet_Pt) {
  //     leadExtraJet_Pt = tmpJetPt; //Pt maximize
  //     leadExtraJet_Eta = jet->eta();
  //     leadExtraJet_Phi = jet->phi();
  //     leadExtraJet_Mass = jet->mass();
  //     leadJet = leadJetCnt; //position of leading jet inside collection
  //   }
  //   double tmpdeltaphimetjet=deltaPhi(met->phi(), jet->phi());
  //   if(tmpdeltaphimetjet<deltaphimetjet) {
  //     deltaphimetjet = tmpdeltaphimetjet; //minimum deltaphi between a jet and met
  //     deltaphimetjet_Pt = jet->pt();
  //     deltaphimetjet_Eta = jet->eta();
  //     deltaphimetjet_Phi = jet->phi();
  //     deltaphimetjet_Mass = jet->mass();
  //     deltaphimetjet_pos = leadJetCnt; //position of closest jet to met inside collection
  //   }
  // }
  // if(leadJet==deltaphimetjet_pos) {
  //   lead_metjet_equal = 1.0;
  // }


  if(debug) {
    cout << "# of selected jets: " << moreSelJetMomenta.size() << endl;
    if(moreSelJetMomenta.size()!=0) {
      cout << " Leading jet: pt=" << jet_pt_vect[0] << endl;
      cout << " Closest jet to MET: n=" << jet_min_deltaPhiJetMET_pos << "; deltaPhi(jet, MET)=" << jet_min_deltaPhiJetMET << endl;
    }
  }


  // =====================================================================
  //  ? veto on the number of jets (< 2) and deltaPhi(jet, MET)
  // =====================================================================


  // Compute the redMET 
  //double lep1Err=( muonLead ? muonLead->track()->ptError() : 
  //		   electronLead->electronMomentumError()*sin(electronLead->theta()) );  // FIXME!!!
  //double lep2Err=( muonSubLead ? muonSubLead->track()->ptError() : 
  //		   electronSubLead->electronMomentumError()*sin(electronSubLead->theta()) );  // FIXME!!!
  double lep1Err=lepton::getPtErrorFor(lep1);
  double lep2Err=lepton::getPtErrorFor(lep2);
  redMETComputer_std->compute(lep1->p4(), lep1Err, 
			      lep2->p4(), lep2Err, 
			      moreSelJetMomenta,
			      met->p4(), 
			      useAllJets);


  // =====================================================================
  //  ? cut on the RedMET
  // =====================================================================



  // // //    __      _______ ________ / \_____      _______  __        ___ ___    ___   // // //
  // // //   |  |    |   ___ | _    _ |\ / ____\     \      \|  |      /   \\  |  |  /   // // //
  // // //   |  |    |  |__ \|/ |  | \| Y \____\|     |  D   |  |     /  O  \\  \/  /    // // //
  // // //   |  |    |   __|    |  |     \____ \      |   __/|  |    |  ___  \\    /     // // //
  // // //   |  |___/|  |___/|  |  |    |\____| |     |  |   |  |___/| /   \  \|  |      // // //
  // // //   |_______|_______| /____\   |______/     /____\  |_______|_\   /___\___\     // // //
  // // // 								                 // // //

  ///////////////////////////////////////////////////////////////////////////////////////////////
  //\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//
  //\\///////////////////////////////////////////////////////////////////////////////////////\\//
  //\\//                                                                                   //\\//
  //\\//       |\                                                                          //\\//
  //\\//       | |                                                                         //\\//
  //\\//       |/                                                                          //\\//
  //\\//       /                                                      /|      |\           //\\//
  //\\//      /|             |\                     0                //|      |\\          //\\//
  //\\// ----/-|-------------|\\---------/|--------|-------0--------|/-|------|-\|-------- //\\//
  //\\// ---/--|-------------|-\|-------//|------0-|------|---------|--|-----0---|-------- //\\//
  //\\// --|--/|--\---------0---|------|/-|-----|-/|------|---0-----|-0----------|-------- //\\//
  //\\// --|-|-|---|------------|------|-0------|//-------|\-|-----0------------O--------- //\\//
  //\\// ---\--|--/------------0-------|--------|/---------\\|---------------------------- //\\//
  //\\//      -+-                     0                     \|                             //\\//
  //\\//     _ |	                                                                   //\\//
  //\\//     \/ 	                                                                   //\\//
  //\\//                                                                                   //\\//
  //\\///////////////////////////////////////////////////////////////////////////////////////\\//
  //\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//
  ///////////////////////////////////////////////////////////////////////////////////////////////

  TString thiscut="Analyzed";
  fillPlots(thiscut.Data(), vtx, totNvtx, goodNvtx, 
	    totNumberOfSelElectrons, totNumberOfSelMuons, 
	    lep1, lep2, leadExtraLept, 
	    rho*Aeff1, rho*Aeff2, rho*Aeff3, 
	    //selIsoTracks, leadExtraTrkV, 
	    met, redMETComputer_std, 
	    moreSelJetMomenta/*, leadJet*/);

  ////////////
  //
  // =====================================================================
  //  Cut1: ask for the dilepton to exist (right flavour, charge)
  // =====================================================================
  //
  thiscut="Dilepton";
  if( chargeCombo!=0 && lep1->charge()*lep2->charge()!=chargeCombo &&  // wrong charge
      CosMu>0.5                                                        // cosmic muon
      ) {
    passedCuts[thiscut.Data()]=0;
    passedCuts["AllCutsPassed"]=0;
  }
  //
  if(passedCuts["AllCutsPassed"]>0.5) 
    fillPlots(thiscut.Data(), vtx, totNvtx, goodNvtx, 
	      totNumberOfSelElectrons, totNumberOfSelMuons, 
	      lep1, lep2, leadExtraLept, 
	      rho*Aeff1, rho*Aeff2, rho*Aeff3, 
	      //selIsoTracks, leadExtraTrkV, 
	      met, redMETComputer_std, 
	      moreSelJetMomenta/*, leadJet*/);
  //
  ////////////


  ////////////
  //
  // =====================================================================
  //  Cut2: apply mass window on the Z mass
  // =====================================================================
  //
  thiscut="MassWindow";
  if(fabs(diLeptonMom.mass()-91.1876)>15.) {
    passedCuts[thiscut.Data()]=0;
    passedCuts["AllCutsPassed"]=0;
  }
  //
  if(passedCuts["AllCutsPassed"]>0.5) 
    fillPlots(thiscut.Data(), vtx, totNvtx, goodNvtx, 
	      totNumberOfSelElectrons, totNumberOfSelMuons, 
	      lep1, lep2, leadExtraLept, 
	      rho*Aeff1, rho*Aeff2, rho*Aeff3, 
	      //selIsoTracks, leadExtraTrkV, 
	      met, redMETComputer_std, 
	      moreSelJetMomenta/*, leadJet*/);
  //
  ////////////


  ////////////
  //
  // =====================================================================
  //  Cut3: veto on the number of jets (< 2)
  // =====================================================================
  //
  thiscut="JetVeto";
  if(moreSelJetMomenta.size()>1) {
    passedCuts[thiscut.Data()]=0;
    passedCuts["AllCutsPassed"]=0;
  }
  //
  if(passedCuts["AllCutsPassed"]>0.5) 
    fillPlots(thiscut.Data(), vtx, totNvtx, goodNvtx, 
	      totNumberOfSelElectrons, totNumberOfSelMuons, 
	      lep1, lep2, leadExtraLept, 
	      rho*Aeff1, rho*Aeff2, rho*Aeff3, 
	      //selIsoTracks, leadExtraTrkV, 
	      met, redMETComputer_std, 
	      moreSelJetMomenta/*, leadJet*/);
  //
  ////////////


  ////////////
  //
  // =====================================================================
  //  Cut4: cut on DeltaPhi(jet, MET)
  // =====================================================================
  //
  thiscut="JetMETDeltaPhi";
  if(jet_min_deltaPhiJetMET<0.5) {
    passedCuts[thiscut.Data()]=0;
    passedCuts["AllCutsPassed"]=0;
  }
  //
  if(passedCuts["AllCutsPassed"]>0.5) 
    fillPlots(thiscut.Data(), vtx, totNvtx, goodNvtx, 
	      totNumberOfSelElectrons, totNumberOfSelMuons, 
	      lep1, lep2, leadExtraLept, 
	      rho*Aeff1, rho*Aeff2, rho*Aeff3, 
	      //selIsoTracks, leadExtraTrkV, 
	      met, redMETComputer_std, 
	      moreSelJetMomenta/*, leadJet*/);
  //
  ////////////



  ////////////
  //
  // =====================================================================
  //  Cut5: anti b-tagging
  // =====================================================================
  //
  thiscut="AntiBTagging";
  if(jet_max_TCHE>2.) {
    passedCuts[thiscut.Data()]=0;
    passedCuts["AllCutsPassed"]=0;
  }
  //
  if(passedCuts["AllCutsPassed"]>0.5) 
    fillPlots(thiscut.Data(), vtx, totNvtx, goodNvtx, 
	      totNumberOfSelElectrons, totNumberOfSelMuons, 
	      lep1, lep2, leadExtraLept, 
	      rho*Aeff1, rho*Aeff2, rho*Aeff3, 
	      //selIsoTracks, leadExtraTrkV, 
	      met, redMETComputer_std, 
	      moreSelJetMomenta/*, leadJet*/);
  //
  ////////////


  ////////////
  //
  // =====================================================================
  //  Cut6: veto on the third lepton
  // =====================================================================
  //
  thiscut="LeptonVeto";
  //  if(totNumberOfSelLeptons>0) {  // N.B. totNumberOfSelLeptons is the number of extra leptons, so must be = 0. 
  bool leptonVetoPass(true);
  if( leadExtraLept.get()!=0 ) {
    // if( abs(leadExtraLeptId)==13 ) {
    //   if( leadExtraLeptNPixHits>0.5     &&
    // 	  leadExtraLeptNTrkHits>10.5    && 
    // 	  ( leadExtraLeptIsGMPT>0.5     || 
    // 	    (leadExtraLeptIsTM>0.5      && 
    // 	     leadExtraLeptNMuMatch>1.5  && 
    // 	     leadExtraLeptTrkChi2<4.) ) && 
    // 	  ( (leadExtraLeptPt>10.        && 
    // 	     leadExtraLeptRelIso<0.3)   || 
    // 	    (leadExtraLeptPt>5.         && 
    // 	     leadExtraLeptRelIso<0.05) )    ) {
    // 	leptonVetoPass=false;
    //   }
    // }
    // else if( abs(leadExtraLeptId)==11 ) {
    //   if( (leadExtraLeptPt>10.        && 
    // 	   leadExtraLeptIdVBTF95>0.5  && 
    // 	   leadExtraLeptRelIso<0.3)   || 
    // 	  (leadExtraLeptPt>5.         && 
    // 	   leadExtraLeptIdVBTF85>0.5  && 
    // 	   leadExtraLeptRelIso<0.05) ) {
    // 	leptonVetoPass=false;
    //   }
    // }

    leptonVetoPass=false;

  }
  if(!leptonVetoPass) {
    passedCuts[thiscut.Data()]=0;
    passedCuts["AllCutsPassed"]=0;
  }
  //
  if(passedCuts["AllCutsPassed"]>0.5) 
    fillPlots(thiscut.Data(), vtx, totNvtx, goodNvtx, 
	      totNumberOfSelElectrons, totNumberOfSelMuons, 
	      lep1, lep2, leadExtraLept, 
	      rho*Aeff1, rho*Aeff2, rho*Aeff3, 
	      //selIsoTracks, leadExtraTrkV, 
	      met, redMETComputer_std, 
	      moreSelJetMomenta/*, leadJet*/);
  //
  ////////////


  ////////////
  //
  // =====================================================================
  //  Cut7: cut on the RedMET
  // =====================================================================
  //
  thiscut="RedMETcut";
  if(redMETComputer_std->reducedMET()<theRedMETMinCut) {
    passedCuts[thiscut.Data()]=0;
    passedCuts["AllCutsPassed"]=0;
  }
  //
  if(passedCuts["AllCutsPassed"]>0.5) 
    fillPlots(thiscut.Data(), vtx, totNvtx, goodNvtx, 
	      totNumberOfSelElectrons, totNumberOfSelMuons, 
	      lep1, lep2, leadExtraLept, 
	      rho*Aeff1, rho*Aeff2, rho*Aeff3, 
	      //selIsoTracks, leadExtraTrkV, 
	      met, redMETComputer_std, 
	      moreSelJetMomenta/*, leadJet*/);
  //
  ////////////


  // ********************* //
  // ********************* //
  // **                 ** //
  // **  Fill the tree  ** //
  // **                 ** //
  // ********************* //
  // ********************* //

  // Some extra variables
  //leadExtraLeptMetTransMass=(leadExtraLept.get()!=0 ? sqrt( 2*met->pt()*leadExtraLept->pt()*(1. - cos(deltaPhi(met->phi(), leadExtraLept->phi()))) ) : -1.);
  //leadExtraTrkMetTransMass=(leadExtraTrk.get()!=0 ? sqrt( 2*met->pt()*leadExtraTrk->pt()*(1. - cos(deltaPhi(met->phi(), leadExtraTrk->phi()))) ) : -1.);
  //leadExtraLeptMetTransMass=(leadExtraLept.get()!=0 ? hMETLeadLeptThird[0]->transverseMass : -1.);
  //leadExtraTrkMetTransMass=(leadExtraTrk.get()!=0 ? hMETLeadIsoTrackThird[0]->transverseMass : -1.);

  float tmpVars[nVarsOpt+nSteps+1];
  unsigned int tmpCnt=0;
  //
  // Variables
  tmpVars[tmpCnt++]=theXsect;
  tmpVars[tmpCnt++]=theBrRatio;
  tmpVars[tmpCnt++]=theLumi;
  tmpVars[tmpCnt++]=event.run();
  tmpVars[tmpCnt++]=event.luminosityBlock();
  tmpVars[tmpCnt++]=event.id().event();
  tmpVars[tmpCnt++]=theWeight;  // only PU weight
  tmpVars[tmpCnt++]=(*selectionInfo)[0]; // flavor combination
  tmpVars[tmpCnt++]=lep1->charge();
  tmpVars[tmpCnt++]=lep2->charge();
  tmpVars[tmpCnt++]=lep1->pt();
  tmpVars[tmpCnt++]=lep2->pt();
  tmpVars[tmpCnt++]=lep1->eta();
  tmpVars[tmpCnt++]=lep2->eta();
  tmpVars[tmpCnt++]=lep1->phi();
  tmpVars[tmpCnt++]=lep2->phi();
  tmpVars[tmpCnt++]=diLeptonMom.pt();
  tmpVars[tmpCnt++]=diLeptonMom.eta();
  tmpVars[tmpCnt++]=diLeptonMom.phi();
  tmpVars[tmpCnt++]=diLeptonMom.mass();
  tmpVars[tmpCnt++]=hDileptKin[0]->dileptLeadDeltaPhi;
  tmpVars[tmpCnt++]=hDileptKin[0]->leptMinusCmCosTheta;
  tmpVars[tmpCnt++]=met->pt();
  tmpVars[tmpCnt++]=met->phi();
  tmpVars[tmpCnt++]=hDileptMET[0]->transverseMass;
  tmpVars[tmpCnt++]=hDileptMET[0]->transverseMassZ;
  tmpVars[tmpCnt++]=hDileptMET[0]->transverseMassZZ;
  tmpVars[tmpCnt++]=hLeadLeptMET[0]->transverseMass;
  tmpVars[tmpCnt++]=hLeadLeptMET[0]->transverseMassZ;
  tmpVars[tmpCnt++]=hLeadLeptMET[0]->transverseMassZZ;
  tmpVars[tmpCnt++]=hSubleadLeptMET[0]->transverseMass;
  tmpVars[tmpCnt++]=hSubleadLeptMET[0]->transverseMassZ;
  tmpVars[tmpCnt++]=hSubleadLeptMET[0]->transverseMassZZ;
  tmpVars[tmpCnt++]=redMETComputer_std->dileptonProjComponents().first; 
  tmpVars[tmpCnt++]=redMETComputer_std->dileptonProjComponents().second;
  tmpVars[tmpCnt++]=redMETComputer_std->dileptonPtCorrComponents().first; 
  tmpVars[tmpCnt++]=redMETComputer_std->dileptonPtCorrComponents().second;
  tmpVars[tmpCnt++]=redMETComputer_std->sumJetProjComponents().first; 
  tmpVars[tmpCnt++]=redMETComputer_std->sumJetProjComponents().second;
  tmpVars[tmpCnt++]=redMETComputer_std->sumAllJetProjComponents().first; 
  tmpVars[tmpCnt++]=redMETComputer_std->sumAllJetProjComponents().second;
  tmpVars[tmpCnt++]=redMETComputer_std->sumOppositeJetProjComponents().first; 
  tmpVars[tmpCnt++]=redMETComputer_std->sumOppositeJetProjComponents().second;
  tmpVars[tmpCnt++]=redMETComputer_std->metProjComponents().first; 
  tmpVars[tmpCnt++]=redMETComputer_std->metProjComponents().second;
  tmpVars[tmpCnt++]=redMETComputer_std->unclProjComponents().first;
  tmpVars[tmpCnt++]=redMETComputer_std->unclProjComponents().second;
  tmpVars[tmpCnt++]=redMETComputer_std->recoilProjComponents().first; 
  tmpVars[tmpCnt++]=redMETComputer_std->recoilProjComponents().second;
  tmpVars[tmpCnt++]=redMETComputer_std->recoilAllJetProjComponents().first; 
  tmpVars[tmpCnt++]=redMETComputer_std->recoilAllJetProjComponents().second;
  tmpVars[tmpCnt++]=redMETComputer_std->recoilOppositeJetProjComponents().first; 
  tmpVars[tmpCnt++]=redMETComputer_std->recoilOppositeJetProjComponents().second;
  tmpVars[tmpCnt++]=redMETComputer_std->reducedMET();
  tmpVars[tmpCnt++]=redMETComputer_std->reducedMETComponents().first;
  tmpVars[tmpCnt++]=redMETComputer_std->reducedMETComponents().second;
  tmpVars[tmpCnt++]=redMETComputer_std->reducedMET(ReducedMETComputer::INDEPENDENTLYMINIMIZED);
  tmpVars[tmpCnt++]=redMETComputer_std->reducedMETComponents(ReducedMETComputer::INDEPENDENTLYMINIMIZED).first;
  tmpVars[tmpCnt++]=redMETComputer_std->reducedMETComponents(ReducedMETComputer::INDEPENDENTLYMINIMIZED).second; 
  tmpVars[tmpCnt++]=totNumberOfSelLeptons;
  tmpVars[tmpCnt++]=leadExtraLeptId;
  tmpVars[tmpCnt++]=leadExtraLeptIsGM;
  tmpVars[tmpCnt++]=leadExtraLeptIsGMPT;
  tmpVars[tmpCnt++]=leadExtraLeptIsTM;
  tmpVars[tmpCnt++]=leadExtraLeptIsTMLSL;
  tmpVars[tmpCnt++]=leadExtraLeptIsTMLST;
  tmpVars[tmpCnt++]=leadExtraLeptIsTMLSAL;
  tmpVars[tmpCnt++]=leadExtraLeptIsTMLSAT;
  tmpVars[tmpCnt++]=leadExtraLeptIdVBTF85;
  tmpVars[tmpCnt++]=leadExtraLeptIdVBTF95;
  tmpVars[tmpCnt++]=leadExtraLeptPt;
  tmpVars[tmpCnt++]=leadExtraLeptEta;
  tmpVars[tmpCnt++]=leadExtraLeptPhi;
  tmpVars[tmpCnt++]=leadExtraLeptTrkChi2;
  tmpVars[tmpCnt++]=leadExtraLeptGlbChi2;
  tmpVars[tmpCnt++]=leadExtraLeptRelIso;
  tmpVars[tmpCnt++]=leadExtraLeptNPixHits;
  tmpVars[tmpCnt++]=leadExtraLeptNTrkHits;
  tmpVars[tmpCnt++]=leadExtraLeptNTrkLostInnHits;
  tmpVars[tmpCnt++]=leadExtraLeptNTrkLostOutHits;
  tmpVars[tmpCnt++]=leadExtraLeptNMuHits;
  tmpVars[tmpCnt++]=leadExtraLeptNMuMatch;
  tmpVars[tmpCnt++]=leadExtraLeptDxy;
  tmpVars[tmpCnt++]=leadExtraLeptDz;
  tmpVars[tmpCnt++]=leadExtraLeptMetTransMass;
  tmpVars[tmpCnt++]=subleadExtraLeptId;
  tmpVars[tmpCnt++]=subleadExtraLeptIsGM;
  tmpVars[tmpCnt++]=subleadExtraLeptIsGMPT;
  tmpVars[tmpCnt++]=subleadExtraLeptIsTM;
  tmpVars[tmpCnt++]=subleadExtraLeptIsTMLSL;
  tmpVars[tmpCnt++]=subleadExtraLeptIsTMLST;
  tmpVars[tmpCnt++]=subleadExtraLeptIsTMLSAL;
  tmpVars[tmpCnt++]=subleadExtraLeptIsTMLSAT;
  tmpVars[tmpCnt++]=subleadExtraLeptIdVBTF85;
  tmpVars[tmpCnt++]=subleadExtraLeptIdVBTF95;
  tmpVars[tmpCnt++]=subleadExtraLeptPt;
  tmpVars[tmpCnt++]=subleadExtraLeptEta;
  tmpVars[tmpCnt++]=subleadExtraLeptPhi;
  tmpVars[tmpCnt++]=subleadExtraLeptTrkChi2;
  tmpVars[tmpCnt++]=subleadExtraLeptGlbChi2;
  tmpVars[tmpCnt++]=subleadExtraLeptRelIso;
  tmpVars[tmpCnt++]=subleadExtraLeptNPixHits;
  tmpVars[tmpCnt++]=subleadExtraLeptNTrkHits;
  tmpVars[tmpCnt++]=subleadExtraLeptNTrkLostInnHits;
  tmpVars[tmpCnt++]=subleadExtraLeptNTrkLostOutHits;
  tmpVars[tmpCnt++]=subleadExtraLeptNMuHits;
  tmpVars[tmpCnt++]=subleadExtraLeptNMuMatch;
  tmpVars[tmpCnt++]=subleadExtraLeptDxy;
  tmpVars[tmpCnt++]=subleadExtraLeptDz;
  tmpVars[tmpCnt++]=subleadExtraLeptMetTransMass;
  tmpVars[tmpCnt++]=CosMu;
  tmpVars[tmpCnt++]=numberOfJets;       // number of jets
  tmpVars[tmpCnt++]=numberOfJetsFromPV; // number of jets from PV
  tmpVars[tmpCnt++]=numberOfJetsFromPU; // number of jets from PU
  tmpVars[tmpCnt++]=jet_pt_vect[0];
  tmpVars[tmpCnt++]=jet_eta_vect[0];
  tmpVars[tmpCnt++]=jet_phi_vect[0];
  tmpVars[tmpCnt++]=jet_mass_vect[0];
  tmpVars[tmpCnt++]=jet_isFromPV_vect[0];
  tmpVars[tmpCnt++]=jet_TCHE_vect[0];
  tmpVars[tmpCnt++]=jet_SSVHE_vect[0];
  tmpVars[tmpCnt++]=jet_JBP_vect[0];
  tmpVars[tmpCnt++]=jet_pt_vect[1];
  tmpVars[tmpCnt++]=jet_eta_vect[1];
  tmpVars[tmpCnt++]=jet_phi_vect[1];
  tmpVars[tmpCnt++]=jet_mass_vect[1];
  tmpVars[tmpCnt++]=jet_isFromPV_vect[1];
  tmpVars[tmpCnt++]=jet_TCHE_vect[1];
  tmpVars[tmpCnt++]=jet_SSVHE_vect[1];
  tmpVars[tmpCnt++]=jet_JBP_vect[1];
  tmpVars[tmpCnt++]=jet_pt_vect[2];
  tmpVars[tmpCnt++]=jet_eta_vect[2];
  tmpVars[tmpCnt++]=jet_phi_vect[2];
  tmpVars[tmpCnt++]=jet_mass_vect[2];
  tmpVars[tmpCnt++]=jet_isFromPV_vect[2];
  tmpVars[tmpCnt++]=jet_TCHE_vect[2];
  tmpVars[tmpCnt++]=jet_SSVHE_vect[2];
  tmpVars[tmpCnt++]=jet_JBP_vect[2];
  tmpVars[tmpCnt++]=jet_pt_vect[3];
  tmpVars[tmpCnt++]=jet_eta_vect[3];
  tmpVars[tmpCnt++]=jet_phi_vect[3];
  tmpVars[tmpCnt++]=jet_mass_vect[3];
  tmpVars[tmpCnt++]=jet_isFromPV_vect[3];
  tmpVars[tmpCnt++]=jet_TCHE_vect[3];
  tmpVars[tmpCnt++]=jet_SSVHE_vect[3];
  tmpVars[tmpCnt++]=jet_JBP_vect[3];
  tmpVars[tmpCnt++]=jet_pt_vect[4];
  tmpVars[tmpCnt++]=jet_eta_vect[4];
  tmpVars[tmpCnt++]=jet_phi_vect[4];
  tmpVars[tmpCnt++]=jet_mass_vect[4];
  tmpVars[tmpCnt++]=jet_isFromPV_vect[4];
  tmpVars[tmpCnt++]=jet_TCHE_vect[4];
  tmpVars[tmpCnt++]=jet_SSVHE_vect[4];
  tmpVars[tmpCnt++]=jet_JBP_vect[4];
  tmpVars[tmpCnt++]=jet_pt_vect[5];
  tmpVars[tmpCnt++]=jet_eta_vect[5];
  tmpVars[tmpCnt++]=jet_phi_vect[5];
  tmpVars[tmpCnt++]=jet_mass_vect[5];
  tmpVars[tmpCnt++]=jet_isFromPV_vect[5];
  tmpVars[tmpCnt++]=jet_TCHE_vect[5];
  tmpVars[tmpCnt++]=jet_SSVHE_vect[5];
  tmpVars[tmpCnt++]=jet_JBP_vect[5];
  tmpVars[tmpCnt++]=jet_pt_vect[6];
  tmpVars[tmpCnt++]=jet_eta_vect[6];
  tmpVars[tmpCnt++]=jet_phi_vect[6];
  tmpVars[tmpCnt++]=jet_mass_vect[6];
  tmpVars[tmpCnt++]=jet_isFromPV_vect[6];
  tmpVars[tmpCnt++]=jet_TCHE_vect[6];
  tmpVars[tmpCnt++]=jet_SSVHE_vect[6];
  tmpVars[tmpCnt++]=jet_JBP_vect[6];
  tmpVars[tmpCnt++]=jet_pt_vect[7];
  tmpVars[tmpCnt++]=jet_eta_vect[7];
  tmpVars[tmpCnt++]=jet_phi_vect[7];
  tmpVars[tmpCnt++]=jet_mass_vect[7];
  tmpVars[tmpCnt++]=jet_isFromPV_vect[7];
  tmpVars[tmpCnt++]=jet_TCHE_vect[7];
  tmpVars[tmpCnt++]=jet_SSVHE_vect[7];
  tmpVars[tmpCnt++]=jet_JBP_vect[7];
  tmpVars[tmpCnt++]=jet_pt_vect[8];
  tmpVars[tmpCnt++]=jet_eta_vect[8];
  tmpVars[tmpCnt++]=jet_phi_vect[8];
  tmpVars[tmpCnt++]=jet_mass_vect[8];
  tmpVars[tmpCnt++]=jet_isFromPV_vect[8];
  tmpVars[tmpCnt++]=jet_TCHE_vect[8];
  tmpVars[tmpCnt++]=jet_SSVHE_vect[8];
  tmpVars[tmpCnt++]=jet_JBP_vect[8];
  tmpVars[tmpCnt++]=jet_pt_vect[9];
  tmpVars[tmpCnt++]=jet_eta_vect[9];
  tmpVars[tmpCnt++]=jet_phi_vect[9];
  tmpVars[tmpCnt++]=jet_mass_vect[9];
  tmpVars[tmpCnt++]=jet_isFromPV_vect[9];
  tmpVars[tmpCnt++]=jet_TCHE_vect[9];
  tmpVars[tmpCnt++]=jet_SSVHE_vect[9];
  tmpVars[tmpCnt++]=jet_JBP_vect[9];
  tmpVars[tmpCnt++]=jet_min_deltaPhiJetMET;
  tmpVars[tmpCnt++]=jet_min_deltaPhiJetMET_pos;
  tmpVars[tmpCnt++]=jet_min_deltaPhiJetMET_isFromPV;
  tmpVars[tmpCnt++]=jet_max_CSV;
  tmpVars[tmpCnt++]=jet_max_CSV_pos;
  tmpVars[tmpCnt++]=jet_max_CSV_isFromPV;
  tmpVars[tmpCnt++]=jet_max_CSVMVA;
  tmpVars[tmpCnt++]=jet_max_CSVMVA_pos;
  tmpVars[tmpCnt++]=jet_max_CSVMVA_isFromPV;
  tmpVars[tmpCnt++]=jet_max_JBP;
  tmpVars[tmpCnt++]=jet_max_JBP_pos;
  tmpVars[tmpCnt++]=jet_max_JBP_isFromPV;
  tmpVars[tmpCnt++]=jet_max_JP;
  tmpVars[tmpCnt++]=jet_max_JP_pos;
  tmpVars[tmpCnt++]=jet_max_JP_isFromPV;
  tmpVars[tmpCnt++]=jet_max_SSVHE;
  tmpVars[tmpCnt++]=jet_max_SSVHE_pos;
  tmpVars[tmpCnt++]=jet_max_SSVHE_isFromPV;
  tmpVars[tmpCnt++]=jet_max_SSVHP;
  tmpVars[tmpCnt++]=jet_max_SSVHP_pos;
  tmpVars[tmpCnt++]=jet_max_SSVHP_isFromPV;
  tmpVars[tmpCnt++]=jet_max_SEBPT;
  tmpVars[tmpCnt++]=jet_max_SEBPT_pos;
  tmpVars[tmpCnt++]=jet_max_SEBPT_isFromPV;
  tmpVars[tmpCnt++]=jet_max_SEBIP3D;
  tmpVars[tmpCnt++]=jet_max_SEBIP3D_pos;
  tmpVars[tmpCnt++]=jet_max_SEBIP3D_isFromPV;
  tmpVars[tmpCnt++]=jet_max_SM;
  tmpVars[tmpCnt++]=jet_max_SM_pos;
  tmpVars[tmpCnt++]=jet_max_SM_isFromPV;
  tmpVars[tmpCnt++]=jet_max_SMBPT;
  tmpVars[tmpCnt++]=jet_max_SMBPT_pos;
  tmpVars[tmpCnt++]=jet_max_SMBPT_isFromPV;
  tmpVars[tmpCnt++]=jet_max_SMBIP3D;
  tmpVars[tmpCnt++]=jet_max_SMBIP3D_pos;
  tmpVars[tmpCnt++]=jet_max_SMBIP3D_isFromPV;
  tmpVars[tmpCnt++]=jet_max_TCHE;
  tmpVars[tmpCnt++]=jet_max_TCHE_pos;
  tmpVars[tmpCnt++]=jet_max_TCHE_isFromPV;
  tmpVars[tmpCnt++]=jet_max_TCHP;
  tmpVars[tmpCnt++]=jet_max_TCHP_pos;
  tmpVars[tmpCnt++]=jet_max_TCHP_isFromPV;
  tmpVars[tmpCnt++]=totNvtx_float;
  tmpVars[tmpCnt++]=goodNvtx_float;
  
  
  if(tmpCnt!=nVarsOpt) {
    cout << "[ZZllvvAnalyzer] *** Error: wrong number of "
	 << "variables used to fill the tree!" << endl;
    throw std::exception();
  }
  //
  // Flags
  for(unsigned int jj=0; jj<nSteps; ++jj) {
    tmpVars[tmpCnt++]=float(passedCuts[allSteps[jj].Data()]);
  }
  if(tmpCnt!=nVarsOpt+nSteps) {
    cout << "[ZZllvvAnalyzer] *** Error: wrong number of "
	 << "flags used to fill the tree!" << endl;
    throw std::exception();
  }
  //
  // Global accept
  tmpVars[tmpCnt++]=float(passedCuts["AllCutsPassed"]);

  //
  // Fill the ntpule
  finalNtpl->Fill(tmpVars);

  //
  // ************************************************************
  // If event is MC debug gen level event
  // ************************************************************
  //

  if(!event.isRealData()) {
    if(debug) cout << "\t Generator level event " << flush;
    int igenpart(0);
    for (pat::eventhypothesis::Looper<reco::GenParticle> genpart = h.loopAs<reco::GenParticle>("genparticle"); genpart; ++genpart) {
      if(debug) cout << "\t" << genpart->pdgId() << " -> " << flush;  
	  
      int igenpartdau(0);
      char buf[20];
      sprintf(buf,"gendaughter_%d",igenpart);
      for(pat::eventhypothesis::Looper<reco::GenParticle> genpartdau = h.loopAs<reco::GenParticle>(buf); genpartdau; ++genpartdau) {
	if(debug) cout << genpartdau->pdgId() << " (" << flush;

	char buf[20];
	sprintf(buf,"gendaughter_%d_%d",igenpart,igenpartdau);
	for(pat::eventhypothesis::Looper<reco::GenParticle> genpartgdau = h.loopAs<reco::GenParticle>(buf); genpartgdau; ++genpartgdau) {
	  if(debug) cout << genpartgdau->pdgId() << " " << flush;
	}
	      
	cout << ") " << flush;
	igenpartdau++;
      }
      igenpart++;
    }
    if(debug) cout << endl;
  }
}


void ZZllvvAnalyzer::endJob() {
  cout << "Tot. # of generated events:               " << nGeneratedEvents       << endl;
  cout << "Tot. # of events after base filters (ee): " << nEventsBaseFilterEE    << endl;
  cout << "Tot. # of events after base filters (mm): " << nEventsBaseFilterMM    << endl;
  cout << "Tot. # of events after base filters (em): " << nEventsBaseFilterEM    << endl;
  cout << "Tot. # of events after skim:              " << nEventsSkimEE+nEventsSkimMM+nEventsSkimEM << endl;
  cout << "Tot. # of events after skim (ee):         " << nEventsSkimEE	         << endl;
  cout << "Tot. # of events after skim (mm):         " << nEventsSkimMM	         << endl;
  cout << "Tot. # of events after skim (em):         " << nEventsSkimEM          << endl;
  cout << "Tot. # of analyzed events:                " << totNEvents             << endl;

  hPreEventCounter->SetBinContent( 1, nGeneratedEvents);
  hPreEventCounter->SetBinContent( 2, nEventsBaseFilterEE);
  hPreEventCounter->SetBinContent( 3, nEventsBaseFilterMM);
  hPreEventCounter->SetBinContent( 4, nEventsBaseFilterEM);
  hPreEventCounter->SetBinContent( 5, nEventsSkimEE+nEventsSkimMM+nEventsSkimEM);
  hPreEventCounter->SetBinContent( 6, nEventsSkimEE);
  hPreEventCounter->SetBinContent( 7, nEventsSkimMM);
  hPreEventCounter->SetBinContent( 8, nEventsSkimEM);
  hPreEventCounter->SetBinContent( 9, totNEvents);
  hPreEventCounter->SetBinContent(10, theLumi);
  hPreEventCounter->SetBinContent(11, theXsect);
  hPreEventCounter->SetBinContent(12, theBrRatio);

  if(isMC) {
    double scaleFact=theLumi*theXsect*theBrRatio/nGeneratedEvents;
    cout << "Luminosity*cross-section*branching-ratio/#events [pb]: " << scaleFact << endl;
    scalePlots(scaleFact);
  }

  // Write the histograms
  theFile->cd();
  hPreEventCounter->Write();
  writePlots();
  finalNtpl->Write();
  theFile->Close();
}



void ZZllvvAnalyzer::beginLuminosityBlock(const edm::LuminosityBlock & iLumi, const edm::EventSetup & iSetup) {}

void ZZllvvAnalyzer::endLuminosityBlock(const edm::LuminosityBlock & iLumi, const edm::EventSetup & iSetup) {
  edm::Handle<edm::MergeableCounter> nGeneratedEvents_lumi;
  edm::Handle<edm::MergeableCounter> nEventsBaseFilterEE_lumi;
  edm::Handle<edm::MergeableCounter> nEventsBaseFilterMM_lumi;
  edm::Handle<edm::MergeableCounter> nEventsBaseFilterEM_lumi;
  edm::Handle<edm::MergeableCounter> nEventsSkimEE_lumi;
  edm::Handle<edm::MergeableCounter> nEventsSkimMM_lumi;
  edm::Handle<edm::MergeableCounter> nEventsSkimEM_lumi;

  iLumi.getByLabel("startCounter",          nGeneratedEvents_lumi);
  iLumi.getByLabel("eePreselectionCounter", nEventsBaseFilterEE_lumi);
  iLumi.getByLabel("mumuPreselectionCounter", nEventsBaseFilterMM_lumi);
  iLumi.getByLabel("emuPreselectionCounter", nEventsBaseFilterEM_lumi);
  iLumi.getByLabel("eeSelectionCounter",    nEventsSkimEE_lumi);
  iLumi.getByLabel("mumuSelectionCounter",    nEventsSkimMM_lumi);
  iLumi.getByLabel("emuSelectionCounter",    nEventsSkimEM_lumi);

  nGeneratedEvents      += nGeneratedEvents_lumi->value;
  nEventsBaseFilterEE   += nEventsBaseFilterEE_lumi->value;
  nEventsBaseFilterMM   += nEventsBaseFilterMM_lumi->value;
  nEventsBaseFilterEM   += nEventsBaseFilterEM_lumi->value;
  nEventsSkimEE         += nEventsSkimEE_lumi->value;
  nEventsSkimMM         += nEventsSkimMM_lumi->value;
  nEventsSkimEM         += nEventsSkimEM_lumi->value;
}


void ZZllvvAnalyzer::initializePlots() {

  if(debug) {
    cout << " Initializing plots..." << endl;
    cout << " - Number of steps = " << nSteps << endl;
  }
  hEventCounter = new TH1F("EventCounter", "Number of events", nSteps, -0.5, float(nSteps)-0.5);
  hNVertexAll   = new TH1F("hNVertexAll", "# of vertices", 100, -0.5, 99.5);
  hNVertexGood  = new TH1F("hNVertexGood", "# of good vertices", 100, -0.5, 99.5);
  for(unsigned int k=0; k<nSteps; ++k) {
    hEventCounter->GetXaxis()->SetBinLabel(k+1, allSteps[k]);
    maps[allSteps[k].Data()]=k;
    if(debug) cout << "   - step " << k << " -> " << allSteps[k].Data() << endl;
    TString cut("cut");
    cut+=k;
    cut+=("_"+allSteps[k]);
    //if(k!=0) theFile->mkdir(cut.Data());  // Do not do it for events without 2 leptons
    theFile->mkdir(cut.Data()); 
    hLeptLead.push_back( new HistoLept( ("LeptonLead_"+cut).Data() ) );
    hLeptSubLead.push_back( new HistoLept( ("LeptonSublead_"+cut).Data() ) );
    hLeptThird.push_back( new HistoLept( ("LeptonThird_"+cut).Data() ) );
    // hSelectedIsoTracks.push_back( new HistoTrack( ("SelectedIsoTracks_"+cut).Data() ) );
    // hSelectedIsoTrackThird.push_back( new HistoTrack( ("SelectedIsoTrackThird_"+cut).Data() ) );
    hDileptKin.push_back( new HistoDilept( ("DileptKin_"+cut).Data() ) );
    hSelJetKin.push_back( new HistoKin( ("SelJetKin_"+cut).Data() ) );
    hLeadJetKin.push_back( new HistoKin( ("LeadJetKin_"+cut).Data() ) );
    hMETKin.push_back( new HistoKin( ("METKin_"+cut).Data() ) );
    hRedMetStd.push_back( new HistoRedMET( ("RedMetStd_"+cut).Data() ) );

    hDileptMET.push_back( new HistoObjectPair( ("DileptMET_"+cut).Data() ) );
    hLeadLeptMET.push_back( new HistoObjectPair( ("LeadLeptMET_"+cut).Data() ) );
    hSubleadLeptMET.push_back( new HistoObjectPair( ("SubleadLeptMET_"+cut).Data() ) );
    hDileptLeadSelJet.push_back( new HistoObjectPair( ("DileptLeadSelJet_"+cut).Data() ) );
    hMETLeadSelJet.push_back( new HistoObjectPair( ("METLeadSelJet_"+cut).Data() ) );
    hMETLeadLeptThird.push_back( new HistoObjectPair( ("METLeadLeptThird_"+cut).Data() ) );
    // hMETLeadIsoTrackThird.push_back( new HistoObjectPair( ("METLeadIsoTrackThird_"+cut).Data() ) );

  }
  hEventCounter->Sumw2();
  hNVertexAll->Sumw2();
  hNVertexGood->Sumw2();
}


void ZZllvvAnalyzer::fillPlots(string cutString, const reco::Vertex *pv, unsigned int nVtx, unsigned int nGoodVtx, 
			       int nElectrons, int nMuons, 
			       reco::CandidatePtr leadLept, reco::CandidatePtr subleadLept, 
			       reco::CandidatePtr thirdLept, 
			       double puOffLead, double puOffSublead, double puOffThird, 
			       //reco::MuonRefVector &selectIsoTrks, reco::MuonRefVector &selectIsoTrkThird, 
			       const pat::MET *missEt, ReducedMETComputer *redmet_std, 
			       vector<LorentzVector> &selJetV/*, int nLeadJ*/) {

  int cutN=maps[cutString];
  hEventCounter->Fill(cutN, theWeight);
  hNVertexAll->Fill(nVtx, theWeight);
  hNVertexGood->Fill(nGoodVtx, theWeight);

  LorentzVector dileptMom=leadLept->p4()+subleadLept->p4();
  hLeptLead[cutN]->FillNLept(nElectrons, theWeight);
  hLeptLead[cutN]->Fill(leadLept, pv, nVtx, puOffLead, theWeight);
  hLeptSubLead[cutN]->FillNLept(nMuons, theWeight);
  hLeptSubLead[cutN]->Fill(subleadLept, pv, nVtx, puOffSublead, theWeight);
  if(thirdLept.get()!=0) {
    hLeptThird[cutN]->FillNLept(nElectrons+nMuons, theWeight);
    hLeptThird[cutN]->Fill(thirdLept, pv, nVtx, puOffThird, theWeight);
    hMETLeadLeptThird[cutN]->Fill(missEt->p4(), thirdLept->p4(), theWeight);
  }
  // hSelectedIsoTracks[cutN]->Fill(selectIsoTrks, pv, nVtx, theWeight);
  // if(selectIsoTrkThird.size()!=0) {
  //   hSelectedIsoTrackThird[cutN]->Fill(selectIsoTrkThird, pv, nVtx, theWeight);
  //   hMETLeadIsoTrackThird[cutN]->Fill(missEt->p4(), selectIsoTrkThird[0]->p4(), theWeight);
  // }
  hDileptKin[cutN]->Fill(leadLept, subleadLept, nVtx, theWeight);
  hSelJetKin[cutN]->FillNObj(selJetV.size(), theWeight);
  for(unsigned int nj=0; nj<selJetV.size(); ++nj) {
    hSelJetKin[cutN]->Fill(selJetV[nj].pt(), selJetV[nj].eta(), selJetV[nj].phi(), selJetV[nj].mass(), nVtx, theWeight);
  }
  if(selJetV.size()!=0) hLeadJetKin[cutN]->Fill(selJetV[0].pt(), selJetV[0].eta(), selJetV[0].phi(), selJetV[0].mass(), nVtx, theWeight);
  hMETKin[cutN]->Fill(missEt->pt(), missEt->eta(), missEt->phi(), missEt->mass(), nVtx, theWeight);  
  hRedMetStd[cutN]->Fill(redmet_std, missEt->pt(), nVtx, theWeight);

  hDileptMET[cutN]->Fill(missEt->p4(), dileptMom, theWeight);
  hLeadLeptMET[cutN]->Fill(missEt->p4(), leadLept->p4(), theWeight);
  hSubleadLeptMET[cutN]->Fill(missEt->p4(), subleadLept->p4(), theWeight);
  if(selJetV.size()!=0) {
    hDileptLeadSelJet[cutN]->Fill(dileptMom, selJetV[0], theWeight);
    hMETLeadSelJet[cutN]->Fill(missEt->p4(), selJetV[0], theWeight);
  }

  if(debug) {
    cout << " - Cut n. " << cutN << ": " << cutString.c_str() << endl;
    cout << " - Weight: " << theWeight << endl;
    cout << " - Lead Lept (pt, eta, phi) = (" 
	 << leadLept->pt() << ", " << leadLept->eta() << ", " << leadLept->phi() << ")" << endl;
    cout << " - Sublead Lept (pt, eta, phi) = (" 
	 << subleadLept->pt() << ", " << subleadLept->eta() << ", " << subleadLept->phi() << ")" << endl;
    if(thirdLept.get()!=0) {
      cout << " - Third Lept (pt, eta, phi) = (" << thirdLept->pt() << ", " 
	   << thirdLept->eta() << ", " << thirdLept->phi() << ")" << endl;
    }
    else cout << " - No Third Lept" << endl;
    // if(selectIsoTrkThird.size()!=0) {
    //   cout << " - Third Track (pt, eta, phi) = (" << selectIsoTrkThird[0]->pt() << ", " 
    // 	   << selectIsoTrkThird[0]->eta() << ", " << selectIsoTrkThird[0]->phi() << ")" 
    // 	   << endl;
    // }
    // else cout << " - No Third Track" << endl;
    cout << " - MET (pt, eta, phi, mass) = (" 
	 << missEt->pt() << ", " << missEt->eta() << ", " << missEt->phi() << ", " << missEt->mass() << ")" << endl;
    cout << " - RedMET [GeV] = " << missEt->pt() << endl;
    cout << " - Number of jets: " << selJetV.size() << endl;
    cout << "                   (pt, eta, phi, mass)" << endl;
    for(unsigned int nj=0; nj<selJetV.size(); ++nj) {
      cout << "       sel jet #" << nj << ": (" << selJetV[nj].pt() << ", " << selJetV[nj].eta() 
	   << ", " << selJetV[nj].phi() << ", " << selJetV[nj].mass() << ")";
      //if(int(nj)==nLeadJ) cout << " ---> LEADING JET" << endl; 
      if(nj==0) cout << " ---> LEADING JET" << endl; 
      else cout << endl;
    }
  }
}


void ZZllvvAnalyzer::scalePlots(double fact) {

  hEventCounter->Scale(fact);
  hNVertexAll->Scale(fact);
  hNVertexGood->Scale(fact);
  for(unsigned int k=1; k<nSteps; ++k) {
    hLeptLead[k]->Scale(fact);
    hLeptSubLead[k]->Scale(fact);
    hLeptThird[k]->Scale(fact);
    // hSelectedIsoTracks[k]->Scale(fact);
    // hSelectedIsoTrackThird[k]->Scale(fact);
    hDileptKin[k]->Scale(fact);
    hSelJetKin[k]->Scale(fact);
    hLeadJetKin[k]->Scale(fact);
    hMETKin[k]->Scale(fact);
    hRedMetStd[k]->Scale(fact);

    hDileptMET[k]->Scale(fact);
    hLeadLeptMET[k]->Scale(fact);
    hSubleadLeptMET[k]->Scale(fact);
    hDileptLeadSelJet[k]->Scale(fact);
    hMETLeadSelJet[k]->Scale(fact);
    hMETLeadLeptThird[k]->Scale(fact);
    // hMETLeadIsoTrackThird[k]->Scale(fact);
  }
}

void ZZllvvAnalyzer::writePlots() {

  theFile->cd();

  hEventCounter->Write();
  hNVertexAll->Write();
  hNVertexGood->Write();
  for(unsigned int k=0; k<nSteps; ++k) {
    TString cut("cut");
    cut+=k;
    cut+=("_"+allSteps[k]);
    theFile->cd(cut.Data());

    hLeptLead[k]->Write();
    hLeptSubLead[k]->Write();
    hLeptThird[k]->Write();
    // hSelectedIsoTracks[k]->Write();
    // hSelectedIsoTrackThird[k]->Write();
    hDileptKin[k]->Write();
    hSelJetKin[k]->Write();
    hLeadJetKin[k]->Write();
    hMETKin[k]->Write();
    hRedMetStd[k]->Write();

    hDileptMET[k]->Write();
    hLeadLeptMET[k]->Write();
    hSubleadLeptMET[k]->Write();
    hDileptLeadSelJet[k]->Write();
    hMETLeadSelJet[k]->Write();
    hMETLeadLeptThird[k]->Write();
    // hMETLeadIsoTrackThird[k]->Write();

    theFile->cd();
  }
}


vector<double> ZZllvvAnalyzer::makePuDistr(const ParameterSet& pSet) {
  vector<double> genPuDistr=pSet.getParameter<vector<double> >("GeneratedPU");
  string pufilename=pSet.getUntrackedParameter<string>("PuDistrFile");
  TFile *file=TFile::Open(pufilename.c_str());
  TH1D *h=0;
  int npuend=genPuDistr.size();
  if(file!=0 && file->IsOpen()) {
    h=(TH1D*)file->Get("pileup");
  }
  vector<double> result(npuend);
  double s=0.;
  for(int npu=0; npu<npuend; ++npu) {
    if(h) {
      if(npu<h->GetNbinsX()) {
	double npu_estimated=h->GetBinContent(h->GetXaxis()->FindBin(npu));
	result[npu]=npu_estimated/genPuDistr[npu];
	s+=npu_estimated;
      }
      else {
	result[npu]=0.;
      }
    }
    else {
      result[npu]=1.;
    }
  }
  // Normalize weights such that the total sum of weights over thw whole sample is 1.0, 
  // i.e., sum_i  result[i] * npu_probs[i] should be 1.0 (!)
  if(h) {
    for(int npu=0; npu<npuend; ++npu){
      result[npu]/=s;
    }
  }
  return result;
}

void ZZllvvAnalyzer::resetFlags() {
  for(map<string, int>::iterator it=passedCuts.begin(); 
      it!=passedCuts.end(); ++it) {
    it->second=1;
  }
}

