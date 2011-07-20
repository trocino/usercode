/*
 *  See header file for a description of this class.
 *
 *  $Date: 2011/07/11 19:28:22 $
 *  $Revision: 1.4 $
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

int nEventPreSkim;
int nEventBaseFilter;
int nEventSkim;

vector<double> puweight;
float theWeight=1.;

// Analysis cuts
TString allSteps[]={"Analyzed", "Dilepton", "MassWindow", 
		    "JetVeto", "LeptonVeto", /*"TrackVeto",*/ 
		    "RedMETcut"};
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
vector<HistoTrack*>      hSelectedIsoTracks; 
vector<HistoTrack*>      hSelectedIsoTrackThird; 
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
vector<HistoObjectPair*> hMETLeadIsoTrackThird;

ReducedMETComputer *redMETComputer_std;

TString varsForLL[]={"run", "lumi", "event", "weight", 
		     "leadCharge", "subleadCharge", 
		     "leadPt", "subleadPt", 
		     "dileptPt", 
		     "dileptInvMass", 
		     "dileptLeadDeltaPhi",
		     "leptMinusCmCosTheta",
		     "redMet",
		     "redMet_long",
		     "redMet_perp",
		     "redMet_dilept_long",
		     "redMet_dilept_perp",
		     "redMet_ptCorr_long",
		     "redMet_ptCorr_perp",
		     "redMet_sumJet_long",
		     "redMet_sumJet_perp",
		     "redMet_met_long",
		     "redMet_met_perp",
		     "redMet_recoil_long",
		     "redMet_recoil_perp",
		     "thirdLept_flavor",
		     "thirdLept_isGlobalMuon",
		     "thirdLept_isGlobalMuonPT",
		     "thirdLept_isTrackerMuon",
		     "thirdLept_isTrackerMuonLSL",
		     "thirdLept_isTrackerMuonLST",
		     "thirdLept_isTrackerMuonLSAL",
		     "thirdLept_isTrackerMuonLSAT",
		     "thirdLept_pt",
		     "thirdLept_eta",
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
		     "thirdTrack_pt",
		     "thirdTrack_eta",
		     "thirdTrack_chi2",
		     "thirdTrack_nTracksIso",
		     "thirdTrack_trkIso",
		     "thirdTrack_relTrkIso",
		     "thirdTrack_nPixHits",
		     "thirdTrack_nTrkHits",
		     "thirdTrack_nTrkLostInnHits",
		     "thirdTrack_nTrkLostOutHits",
		     "thirdTrack_dxy",
		     "thirdTrack_dz",
		     "thirdTrack_METtrackTransvMass"};
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
  nEventPreSkim = 0;
  nEventBaseFilter = 0;
  nEventSkim = 0;
  flavorCombo = pSet.getUntrackedParameter<int>("FlavorCombination", 0);
  chargeCombo = pSet.getUntrackedParameter<int>("ChargeCombination", -1);
  kRecoilLongWeight = pSet.getUntrackedParameter<double>("RecoilLongWeight", 2.);
  kRecoilPerpWeight = pSet.getUntrackedParameter<double>("RecoilPerpWeight", 2.);
  kSigmaPtLongWeight = pSet.getUntrackedParameter<double>("SigmaPtLongWeight", 2.5);
  kSigmaPtPerpWeight = pSet.getUntrackedParameter<double>("SigmaPtPerpWeight", 2.5);
  kPerpComponentWeight = pSet.getUntrackedParameter<double>("PerpComponentWeight", 1.);
  theRedMETMinCut = pSet.getUntrackedParameter<double>("RedMETMinCut", 50.);
  redMETComputer_std   = new ReducedMETComputer(kRecoilLongWeight,
						kRecoilPerpWeight,
						kSigmaPtLongWeight,
						kSigmaPtPerpWeight,
						kPerpComponentWeight);
  vertexSelection =  pSet.getParameter<ParameterSet>("Vertices");
  if(isMC) {
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
  hPreEventCounter=new TH1F("PreEventCounter", "Unscaled number of events and more", 7, 0., 7.);
  hPreEventCounter->GetXaxis()->SetBinLabel(1, "Before-skim");
  hPreEventCounter->GetXaxis()->SetBinLabel(2, "Pre-filters");
  hPreEventCounter->GetXaxis()->SetBinLabel(3, "After-skim");
  hPreEventCounter->GetXaxis()->SetBinLabel(4, "Analyzed");
  hPreEventCounter->GetXaxis()->SetBinLabel(5, "Integrated-luminosity");
  hPreEventCounter->GetXaxis()->SetBinLabel(6, "Cross-section");
  hPreEventCounter->GetXaxis()->SetBinLabel(7, "Branching-ratio");

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
    // Retrieve Pile-Up information
    Handle<vector<PileupSummaryInfo> > PupInfo;
    event.getByLabel("addPileupInfo", PupInfo);
    std::vector<PileupSummaryInfo>::const_iterator pvi;
    unsigned int nPuInt=0;
    for(pvi = PupInfo->begin(); pvi!=PupInfo->end(); ++pvi) {
      if(debug) {
	cout << " Pileup Information: bunchXing: " << pvi->getBunchCrossing() << ", nvtx: " << pvi->getPU_NumInteractions() << endl;
      }
      nPuInt+=pvi->getPU_NumInteractions();
    }

    if(debug) {
      cout << " Total number of PU interaction: " << nPuInt << endl;
    }

    if(nPuInt>=puweight.size()) {
      cout << "[ZZllvvAnalyzer] *** Error: n. PU interactions (" << nPuInt 
	   << ") greater than the size of weights vector (" << puweight.size() 
	   << ")! " << endl;
      throw std::exception();
    }

    theWeight=puweight[nPuInt];
  } // end if(isMC)
  else {  // Real data
    theWeight=1.;   // already 1. at initialization, but you never know...
  }

  // Count all the events (EventCounter)
  hEventCounter->Fill(0., theWeight);

  // Pre-select vertices
  Handle<reco::VertexCollection> offlinePrimVertices;
  event.getByLabel("offlinePrimaryVerticesDA", offlinePrimVertices);  
  unsigned int totNvtx=offlinePrimVertices->size();
  if(debug) cout << "Total # of vertices: " << totNvtx << endl;;
  std::vector<reco::VertexRef> selVertices = vertex::filter(offlinePrimVertices,vertexSelection);
  unsigned int goodNvtx=selVertices.size();
  if(debug) cout << "# of good vertices: " << goodNvtx << endl;;

  // Average energy density
  edm::Handle< double > rhoH;
  event.getByLabel(InputTag("kt6PFJets:rho"), rhoH);
  double rho=(*rhoH);

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
  CandidatePtr lep1 = h["leg1"];
  CandidatePtr lep2 = h["leg2"];

  // ... try casting to Muons...
  const pat::Muon *muonLead    = dynamic_cast<const pat::Muon *>(lep1.get());
  const pat::Muon *muonSubLead = dynamic_cast<const pat::Muon *>(lep2.get());

  // ... and to Electrons
  const pat::Electron *electronLead    = dynamic_cast<const pat::Electron *>(lep1.get());
  const pat::Electron *electronSubLead = dynamic_cast<const pat::Electron *>(lep2.get());

  // TEMPORARY!!!!
  double muonEffectiveArea=0.112;
  double electronEffectiveArea=0.24;
  double Aeff1=(muonLead ? muonEffectiveArea : electronEffectiveArea);
  double Aeff2=(muonSubLead ? muonEffectiveArea : electronEffectiveArea);

  vector<double> lepIso1=dilepton::getLeptonIso(lep1);
  double relIso1=(lepIso1[dilepton::TRACKER_ISO]
		  +max(lepIso1[dilepton::ECAL_ISO]
		       +lepIso1[dilepton::HCAL_ISO]
		       -Aeff1*rho, 0.0)) / lep1->pt();
  vector<double> lepIso2=dilepton::getLeptonIso(lep2);
  double relIso2=(lepIso2[dilepton::TRACKER_ISO]
		  +max(lepIso2[dilepton::ECAL_ISO]
		       +lepIso2[dilepton::HCAL_ISO]
		       -Aeff2*rho, 0.0)) / lep2->pt();

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
  float leadExtraLeptPt=0.;
  float leadExtraLeptEta=-9999.;
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

  // Get the collection of selected, isolated muons
  vector<const pat::Muon *> allSelMuons;
  for(pat::eventhypothesis::Looper<pat::Muon> muon=h.loopAs<pat::Muon>("muon"); muon; ++muon) {
    if(debug) 
      cout << "\t muon: " << muon->pt() << "; " << muon->eta() << "; " << muon->phi() << endl;
    const pat::Muon *tmpMu=dynamic_cast<const pat::Muon *>(muon.get());
    allSelMuons.push_back(tmpMu);
    if(muon->pt()>leadExtraLeptPt) {
      leadExtraLept=muon.ref();
      leadExtraLeptPt=muon->pt();
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
      leadExtraLept=ele.ref();
      leadExtraLeptPt=ele->pt();
      leadExtraLeptType="electron";
      leadExtraLeptId=11.*ele->charge();
    }
  }
  if(debug) {
    cout << "# of selected, isolated electrons: " << allSelElectrons.size() << endl;
  }

  const pat::Muon *mu3=dynamic_cast<const pat::Muon *>(leadExtraLept.get());
  const pat::Electron *el3=dynamic_cast<const pat::Electron *>(leadExtraLept.get());

  if(leadExtraLept.get()!=0) {
    leadExtraLeptEta=leadExtraLept->eta();
    double Aeff3;
    if( abs(leadExtraLeptId)==13 ) {
      if( mu3==0 ) {
	cout << "[ZZllvvAnalyzer] *** Error: Leading extra lepton is a muon, "
	     << "but dynamic_cast failed!" << endl;
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
    else if( abs(leadExtraLeptId)==11 ) {
      if( el3==0 ) {
	cout << "[ZZllvvAnalyzer] *** Error: Leading extra lepton is an electron, "
	     << "but dynamic_cast failed!" << endl;
	throw std::exception();
      }
      leadExtraLeptTrkChi2=el3->gsfTrack()->normalizedChi2();
      leadExtraLeptNPixHits=el3->gsfTrack()->hitPattern().numberOfValidPixelHits();
      leadExtraLeptNTrkHits=el3->gsfTrack()->hitPattern().numberOfValidTrackerHits();
      leadExtraLeptNTrkLostInnHits=el3->gsfTrack()->trackerExpectedHitsInner().numberOfLostHits();
      leadExtraLeptNTrkLostOutHits=el3->gsfTrack()->trackerExpectedHitsOuter().numberOfLostHits();
      leadExtraLeptDxy=fabs( el3->gsfTrack()->dxy(vtx->position()) );
      leadExtraLeptDz=fabs( el3->gsfTrack()->dz(vtx->position()) );
      Aeff3=electronEffectiveArea;
    }
    vector<double> lepIso3=dilepton::getLeptonIso(leadExtraLept);
    leadExtraLeptRelIso=(lepIso3[dilepton::TRACKER_ISO]
			 +max(lepIso3[dilepton::ECAL_ISO]
			      +lepIso3[dilepton::HCAL_ISO]
			      -Aeff3*rho, 0.0)) / leadExtraLeptPt;    
    if(debug) {
      cout << " Leading extra lepton: " << leadExtraLeptType.c_str() 
	   << ", pt=" << leadExtraLeptPt << "; relIsoCorr=" << leadExtraLeptRelIso << endl;
    }
  }

  int totNumberOfSelElectrons=allSelElectrons.size();
  int totNumberOfSelMuons=allSelMuons.size();
  int totNumberOfSelLeptons=totNumberOfSelElectrons+totNumberOfSelMuons;

  // =====================================================================
  //  ? veto on the third lepton
  // =====================================================================


  // Get the "leading" extra track
  reco::MuonRef leadExtraTrk;
  reco::MuonRefVector leadExtraTrkV; // to comply with Histograms.h format
  float leadExtraTrkPt=0.;
  float leadExtraTrkEta=-9999.;
  float leadExtraTrkTrkChi2=-1.;
  float leadExtraTrkIsoNTrks=-1.;
  float leadExtraTrkIsoSumPt=-1.;
  float leadExtraTrkRelIso=-1.;
  float leadExtraTrkNPixHits=-1.;
  float leadExtraTrkNTrkHits=-1.;
  float leadExtraTrkNTrkLostInnHits=-1.;
  float leadExtraTrkNTrkLostOutHits=-1.;
  float leadExtraTrkDxy=-1.;
  float leadExtraTrkDz=-1.;
  float leadExtraTrkMetTransMass=-1.;
  // Get the collection of selected, isolated tracks
  edm::Handle<reco::MuonCollection> hTracks;
  reco::MuonRefVector selIsoTracks;
  event.getByLabel(edm::InputTag(source.label()+":selectedTracks"), hTracks);
  for(reco::MuonCollection::size_type itTrk=0; itTrk<hTracks->size(); ++itTrk) {
    reco::MuonRef mRef(hTracks, itTrk);
    if(debug) 
      cout << "\t track: " << mRef->pt() << "; " << mRef->eta() << "; " << mRef->phi() << endl
	   << " Iso (nTracks, sumPt): " << mRef->isolationR03().nTracks << ", " << mRef->isolationR03().sumPt;

    selIsoTracks.push_back(mRef);
    if(mRef->pt()>leadExtraTrkPt) {
      leadExtraTrk=mRef;
      leadExtraTrkPt=mRef->pt();
    }

    // // Accept only if from primary vertex
    // if(mRef->innerTrack()->dz(vtx->position())<0.1) {
    //   if(debug) 
    //     cout << "; from primary vertex: accepted." << endl;
    //   selIsoTracks.push_back(mRef);
    //   if(mRef->pt()>leadExtraTrkPt) {
    //     leadExtraTrk=mRef;
    //     leadExtraTrkPt=mRef->pt();
    //   }
    // }
    // else 
    //   if(debug) cout << "; not from primary vertex: rejected." << endl;

  } // end for(reco::MuonCollection::size_type itTrk=0; itTrk<hTracks->size(); ++itTrk)

  if(leadExtraTrk.get()!=0) {
    leadExtraTrkV.push_back(leadExtraTrk);
    leadExtraTrkEta=leadExtraTrk->eta();
    leadExtraTrkTrkChi2=leadExtraTrk->innerTrack()->normalizedChi2();
    leadExtraTrkIsoNTrks=leadExtraTrk->isolationR03().nTracks;
    leadExtraTrkIsoSumPt=leadExtraTrk->isolationR03().sumPt;
    leadExtraTrkRelIso=leadExtraTrk->isolationR03().sumPt/leadExtraTrkPt;
    leadExtraTrkNPixHits=leadExtraTrk->innerTrack()->hitPattern().numberOfValidPixelHits();
    leadExtraTrkNTrkHits=leadExtraTrk->innerTrack()->hitPattern().numberOfValidTrackerHits();
    leadExtraTrkNTrkLostInnHits=leadExtraTrk->innerTrack()->trackerExpectedHitsInner().numberOfLostHits();
    leadExtraTrkNTrkLostOutHits=leadExtraTrk->innerTrack()->trackerExpectedHitsOuter().numberOfLostHits();
    leadExtraTrkDxy=fabs( leadExtraTrk->innerTrack()->dxy(vtx->position()) );
    leadExtraTrkDz=fabs( leadExtraTrk->innerTrack()->dz(vtx->position()) );
  }

  if(debug) {
    cout << "# of selected tracks: " << selIsoTracks.size() << endl;
    if(leadExtraTrk.get()!=0) 
      cout << " Leading extra track: pt=" << leadExtraTrkPt << endl;
  }

  int totNumberOfSelTracks=selIsoTracks.size();
  int totNumberOfSelLeptonsTracks=totNumberOfSelLeptons+totNumberOfSelTracks;


  // =====================================================================
  //  ? veto on the third track
  // =====================================================================


  // Get the "leading" jet
  int leadJet=-1, leadJetCnt=0;
  double leadExtraJetPt=0.;
  // Get the Jet momenta
  vector<CandidatePtr> allSelJets;
  vector<LorentzVector> moreSelJetMomenta;
  for(pat::eventhypothesis::Looper<pat::Jet> jet=h.loopAs<pat::Jet>("jet"); jet; ++jet, ++leadJetCnt) {
    double tmpJetPt=jet->pt();
    if(debug) 
      cout << "\t jet: " << tmpJetPt << ";" << jet->eta() << ";" << jet->phi() << std::endl;
    moreSelJetMomenta.push_back(jet->p4());
    if(tmpJetPt>leadExtraJetPt) {
      leadExtraJetPt=tmpJetPt;
      leadJet=leadJetCnt;
    }
  }
  if(debug) {
    cout << "# of selected jets: " << moreSelJetMomenta.size() << endl;
    if(leadJet!=-1) 
      cout << " Leading jet: n= " << leadJet << "; pt=" << leadExtraJetPt << endl;
  }


  // =====================================================================
  //  ? veto on the number of jets (< 2)
  // =====================================================================


  // Get the MET
  const pat::MET *met = h.getAs<pat::MET>("met");
  if(debug) cout << "\t MET:" << met->pt() << ";" << met->phi() << endl;

  // Compute the redMET 
  double lep1Err=( muonLead ? muonLead->track()->ptError() : 
		   electronLead->electronMomentumError()*sin(electronLead->theta()) );  // FIXME!!!
  double lep2Err=( muonSubLead ? muonSubLead->track()->ptError() : 
		   electronSubLead->electronMomentumError()*sin(electronSubLead->theta()) );  // FIXME!!!
  redMETComputer_std->compute(lep1->p4(), lep1Err, 
			      lep2->p4(), lep2Err, 
			      moreSelJetMomenta,
			      met->p4());


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

  TString thiscut="";

  ////////////
  //
  // =====================================================================
  //  Cut0: ask for the dilepton to exist (right flavour, charge)
  // =====================================================================
  //
  thiscut="Dilepton";
  if( chargeCombo!=0 && lep1->charge()*lep2->charge()!=chargeCombo ) {
    passedCuts[thiscut.Data()]=0;
    passedCuts["AllCutsPassed"]=0;
  }
  //
  if(passedCuts["AllCutsPassed"]>0.5) 
    fillPlots(thiscut.Data(), vtx, totNvtx, goodNvtx, 
	      totNumberOfSelElectrons, totNumberOfSelMuons, 
	      lep1, lep2, leadExtraLept, 
	      selIsoTracks, leadExtraTrkV, 
	      met, redMETComputer_std, 
	      moreSelJetMomenta, leadJet);
  //
  ////////////


  ////////////
  //
  // =====================================================================
  //  Cut1: apply mass window on the Z mass
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
	      selIsoTracks, leadExtraTrkV, 
	      met, redMETComputer_std, 
	      moreSelJetMomenta, leadJet);
  //
  ////////////


  ////////////
  //
  // =====================================================================
  //  Cut2: veto on the number of jets (< 2)
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
	      selIsoTracks, leadExtraTrkV, 
	      met, redMETComputer_std, 
	      moreSelJetMomenta, leadJet);
  //
  ////////////


  ////////////
  //
  // =====================================================================
  //  Cut3: veto on the third lepton
  // =====================================================================
  //
  thiscut="LeptonVeto";
  //  if(totNumberOfSelLeptons>0) {  // N.B. totNumberOfSelLeptons is the number of extra leptons, so must be = 0. 
  bool leptonVetoPass(true);
  if( leadExtraLept.get()!=0 ) {
    if( abs(leadExtraLeptId)==13 ) {
      if( leadExtraLeptNPixHits>0.5     &&
	  leadExtraLeptNTrkHits>10.5    && 
	  ( leadExtraLeptIsGMPT>0.5     || 
	    (leadExtraLeptIsTM>0.5      && 
	     leadExtraLeptNMuMatch>1.5  && 
	     leadExtraLeptTrkChi2<4.) ) && 
	  ( (leadExtraLeptPt>10.        && 
	     leadExtraLeptRelIso<0.3)   || 
	    (leadExtraLeptPt>5.         && 
	     leadExtraLeptRelIso<0.05) )    ) {
	leptonVetoPass=false;
      }
    }
    else if( abs(leadExtraLeptId)==11 ) {
      if( (leadExtraLeptPt>10.        && 
	   leadExtraLeptRelIso<0.3)   || 
	  (leadExtraLeptPt>5.         && 
	   leadExtraLeptRelIso<0.05) ) {
	leptonVetoPass=false;
      }
    }
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
	      selIsoTracks, leadExtraTrkV, 
	      met, redMETComputer_std, 
	      moreSelJetMomenta, leadJet);
  //
  ////////////


  /*
  ////////////
  //
  // =====================================================================
  //  Cut4: veto on the third isolated track
  // =====================================================================
  //
  thiscut="TrackVeto";
  if(totNumberOfSelTracks>0) {  
    passedCuts[thiscut.Data()]=0;
    passedCuts["AllCutsPassed"]=0;
  }
  //
  if(passedCuts["AllCutsPassed"]>0.5) 
    fillPlots(thiscut.Data(), vtx, totNvtx, goodNvtx, 
	      totNumberOfSelElectrons, totNumberOfSelMuons, 
	      lep1, lep2, leadExtraLept, 
	      selIsoTracks, leadExtraTrkV, 
	      met, redMETComputer_std, 
	      moreSelJetMomenta, leadJet);
  //
  ////////////
  */

  ////////////
  //
  // =====================================================================
  //  Cut5: cut on the RedMET
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
	      selIsoTracks, leadExtraTrkV, 
	      met, redMETComputer_std, 
	      moreSelJetMomenta, leadJet);
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
  leadExtraLeptMetTransMass=(leadExtraLept.get()!=0 ? sqrt( 2*met->pt()*leadExtraLept->pt()*(1. - cos(deltaPhi(met->phi(), leadExtraLept->phi()))) ) : -1.);
  leadExtraTrkMetTransMass=(leadExtraTrk.get()!=0 ? sqrt( 2*met->pt()*leadExtraTrk->pt()*(1. - cos(deltaPhi(met->phi(), leadExtraTrk->phi()))) ) : -1.);

  float tmpVars[nVarsOpt+nSteps+1];
  unsigned int tmpCnt=0;
  //
  // Variables
  tmpVars[tmpCnt++]=event.run();
  tmpVars[tmpCnt++]=event.luminosityBlock();
  tmpVars[tmpCnt++]=event.id().event();
  tmpVars[tmpCnt++]=theWeight;  // only PU weight
  tmpVars[tmpCnt++]=lep1->charge();
  tmpVars[tmpCnt++]=lep2->charge();
  tmpVars[tmpCnt++]=lep1->pt();
  tmpVars[tmpCnt++]=lep2->pt();
  tmpVars[tmpCnt++]=diLeptonMom.pt();
  tmpVars[tmpCnt++]=diLeptonMom.mass();
  tmpVars[tmpCnt++]=hDileptKin.back()->dileptLeadDeltaPhi;
  tmpVars[tmpCnt++]=hDileptKin.back()->leptMinusCmCosTheta;
  tmpVars[tmpCnt++]=redMETComputer_std->reducedMET();
  tmpVars[tmpCnt++]=redMETComputer_std->reducedMETComponents().first;
  tmpVars[tmpCnt++]=redMETComputer_std->reducedMETComponents().second;
  tmpVars[tmpCnt++]=redMETComputer_std->dileptonProjComponents().first; 
  tmpVars[tmpCnt++]=redMETComputer_std->dileptonProjComponents().second;
  tmpVars[tmpCnt++]=redMETComputer_std->dileptonPtCorrComponents().first; 
  tmpVars[tmpCnt++]=redMETComputer_std->dileptonPtCorrComponents().second;
  tmpVars[tmpCnt++]=redMETComputer_std->sumJetProjComponents().first; 
  tmpVars[tmpCnt++]=redMETComputer_std->sumJetProjComponents().second;
  tmpVars[tmpCnt++]=redMETComputer_std->metProjComponents().first; 
  tmpVars[tmpCnt++]=redMETComputer_std->metProjComponents().second;
  tmpVars[tmpCnt++]=redMETComputer_std->recoilProjComponents().first; 
  tmpVars[tmpCnt++]=redMETComputer_std->recoilProjComponents().second;
  tmpVars[tmpCnt++]=leadExtraLeptId;
  tmpVars[tmpCnt++]=leadExtraLeptIsGM;
  tmpVars[tmpCnt++]=leadExtraLeptIsGMPT;
  tmpVars[tmpCnt++]=leadExtraLeptIsTM;
  tmpVars[tmpCnt++]=leadExtraLeptIsTMLSL;
  tmpVars[tmpCnt++]=leadExtraLeptIsTMLST;
  tmpVars[tmpCnt++]=leadExtraLeptIsTMLSAL;
  tmpVars[tmpCnt++]=leadExtraLeptIsTMLSAT;
  tmpVars[tmpCnt++]=leadExtraLeptPt;
  tmpVars[tmpCnt++]=leadExtraLeptEta;
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
  tmpVars[tmpCnt++]=leadExtraTrkPt;
  tmpVars[tmpCnt++]=leadExtraTrkEta;
  tmpVars[tmpCnt++]=leadExtraTrkTrkChi2;
  tmpVars[tmpCnt++]=leadExtraTrkIsoNTrks;
  tmpVars[tmpCnt++]=leadExtraTrkIsoSumPt;
  tmpVars[tmpCnt++]=leadExtraTrkRelIso;
  tmpVars[tmpCnt++]=leadExtraTrkNPixHits;
  tmpVars[tmpCnt++]=leadExtraTrkNTrkHits;
  tmpVars[tmpCnt++]=leadExtraTrkNTrkLostInnHits;
  tmpVars[tmpCnt++]=leadExtraTrkNTrkLostOutHits;
  tmpVars[tmpCnt++]=leadExtraTrkDxy;
  tmpVars[tmpCnt++]=leadExtraTrkDz;
  tmpVars[tmpCnt++]=leadExtraTrkMetTransMass;
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
  cout << "Tot. # of events pre skim: "           << nEventPreSkim    << endl;
  cout << "Tot. # of events after base filters: " << nEventBaseFilter << endl;
  cout << "Tot. # of events after skim: "         << nEventSkim       << endl;
  cout << "Tot. # of events in analysis: "        << totNEvents       << endl;

  hPreEventCounter->SetBinContent(1, nEventPreSkim);
  hPreEventCounter->SetBinContent(2, nEventBaseFilter);
  hPreEventCounter->SetBinContent(3, nEventSkim);
  hPreEventCounter->SetBinContent(4, totNEvents);
  hPreEventCounter->SetBinContent(5, theLumi);
  hPreEventCounter->SetBinContent(6, theXsect);
  hPreEventCounter->SetBinContent(7, theBrRatio);

  if(isMC) {
    double scaleFact=theLumi*theXsect*theBrRatio/nEventPreSkim;
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
  edm::Handle<edm::MergeableCounter> nEventPreSkim_lumi;
  iLumi.getByLabel("startCounter", nEventPreSkim_lumi);
  nEventPreSkim += nEventPreSkim_lumi->value;

  edm::Handle<edm::MergeableCounter> nEventBaseFilter_lumi;
  iLumi.getByLabel("preFilterCounter", nEventBaseFilter_lumi);
  nEventBaseFilter += nEventBaseFilter_lumi->value;

  edm::Handle<edm::MergeableCounter> nEventSkim_lumi;
  iLumi.getByLabel("mumuCounter", nEventSkim_lumi);
  nEventSkim += nEventSkim_lumi->value;
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
    if(k!=0) theFile->mkdir(cut.Data());  // Do not do it for events without 2 leptons
    hLeptLead.push_back( new HistoLept( ("LeptonLead_"+cut).Data() ) );
    hLeptSubLead.push_back( new HistoLept( ("LeptonSublead_"+cut).Data() ) );
    hLeptThird.push_back( new HistoLept( ("LeptonThird_"+cut).Data() ) );
    hSelectedIsoTracks.push_back( new HistoTrack( ("SelectedIsoTracks_"+cut).Data() ) );
    hSelectedIsoTrackThird.push_back( new HistoTrack( ("SelectedIsoTrackThird_"+cut).Data() ) );
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
    hMETLeadIsoTrackThird.push_back( new HistoObjectPair( ("METLeadIsoTrackThird_"+cut).Data() ) );

  }
  hEventCounter->Sumw2();
  hNVertexAll->Sumw2();
  hNVertexGood->Sumw2();
}


void ZZllvvAnalyzer::fillPlots(string cutString, const reco::Vertex *pv, unsigned int nVtx, unsigned int nGoodVtx, 
			       int nElectrons, int nMuons, 
			       const reco::CandidatePtr leadLept, const reco::CandidatePtr subleadLept, 
			       reco::CandidatePtr thirdLept, 
			       reco::MuonRefVector &selectIsoTrks, reco::MuonRefVector &selectIsoTrkThird, 
			       const pat::MET *missEt, ReducedMETComputer *redmet_std, 
			       vector<LorentzVector> &selJetV, int nLeadJ) {

  int cutN=maps[cutString];
  hEventCounter->Fill(cutN, theWeight);
  hNVertexAll->Fill(nVtx, theWeight);
  hNVertexGood->Fill(nGoodVtx, theWeight);

  LorentzVector dileptMom=leadLept->p4()+subleadLept->p4();
  hLeptLead[cutN]->FillNLept(nElectrons, theWeight);
  hLeptLead[cutN]->Fill(leadLept, pv, nVtx, theWeight);
  hLeptSubLead[cutN]->FillNLept(nMuons, theWeight);
  hLeptSubLead[cutN]->Fill(subleadLept, pv, nVtx, theWeight);
  if(thirdLept.get()!=0) {
    hLeptThird[cutN]->FillNLept(nElectrons+nMuons, theWeight);
    hLeptThird[cutN]->Fill(thirdLept, pv, nVtx, theWeight);
    hMETLeadLeptThird[cutN]->Fill(missEt->p4(), thirdLept->p4(), theWeight);
  }
  hSelectedIsoTracks[cutN]->Fill(selectIsoTrks, pv, nVtx, theWeight);
  if(selectIsoTrkThird.size()!=0) {
    hSelectedIsoTrackThird[cutN]->Fill(selectIsoTrkThird, pv, nVtx, theWeight);
    hMETLeadIsoTrackThird[cutN]->Fill(missEt->p4(), selectIsoTrkThird[0]->p4(), theWeight);
  }
  hDileptKin[cutN]->Fill(leadLept, subleadLept, nVtx, theWeight);
  hSelJetKin[cutN]->FillNObj(selJetV.size(), theWeight);
  for(unsigned int nj=0; nj<selJetV.size(); ++nj) {
    hSelJetKin[cutN]->Fill(selJetV[nj].pt(), selJetV[nj].eta(), selJetV[nj].phi(), selJetV[nj].mass(), nVtx, theWeight);
  }
  if(nLeadJ!=-1) hLeadJetKin[cutN]->Fill(selJetV[nLeadJ].pt(), selJetV[nLeadJ].eta(), selJetV[nLeadJ].phi(), selJetV[nLeadJ].mass(), nVtx, theWeight);
  hMETKin[cutN]->Fill(missEt->pt(), missEt->eta(), missEt->phi(), missEt->mass(), nVtx, theWeight);  
  hRedMetStd[cutN]->Fill(redmet_std, missEt->pt(), nVtx, theWeight);

  hDileptMET[cutN]->Fill(dileptMom, missEt->p4(), theWeight);
  hLeadLeptMET[cutN]->Fill(leadLept->p4(), missEt->p4(), theWeight);
  hSubleadLeptMET[cutN]->Fill(subleadLept->p4(), missEt->p4(), theWeight);
  if(nLeadJ!=-1) {
    hDileptLeadSelJet[cutN]->Fill(dileptMom, selJetV[nLeadJ], theWeight);
    hMETLeadSelJet[cutN]->Fill(missEt->p4(), selJetV[nLeadJ], theWeight);
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
    if(selectIsoTrkThird.size()!=0) {
      cout << " - Third Track (pt, eta, phi) = (" << selectIsoTrkThird[0]->pt() << ", " 
	   << selectIsoTrkThird[0]->eta() << ", " << selectIsoTrkThird[0]->phi() << ")" 
	   << endl;
    }
    else cout << " - No Third Track" << endl;
    cout << " - MET (pt, eta, phi, mass) = (" 
	 << missEt->pt() << ", " << missEt->eta() << ", " << missEt->phi() << ", " << missEt->mass() << ")" << endl;
    cout << " - RedMET [GeV] = " << missEt->pt() << endl;
    cout << " - Number of jets: " << selJetV.size() << endl;
    cout << "                   (pt, eta, phi, mass)" << endl;
    for(unsigned int nj=0; nj<selJetV.size(); ++nj) {
      cout << "       sel jet #" << nj << ": (" << selJetV[nj].pt() << ", " << selJetV[nj].eta() 
	   << ", " << selJetV[nj].phi() << ", " << selJetV[nj].mass() << ")";
      if(int(nj)==nLeadJ) cout << " ---> LEADING JET" << endl; 
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
    hSelectedIsoTracks[k]->Scale(fact);
    hSelectedIsoTrackThird[k]->Scale(fact);
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
    hMETLeadIsoTrackThird[k]->Scale(fact);
  }
}

void ZZllvvAnalyzer::writePlots() {

  theFile->cd();

  hEventCounter->Write();
  hNVertexAll->Write();
  hNVertexGood->Write();
  for(unsigned int k=1; k<nSteps; ++k) {
    TString cut("cut");
    cut+=k;
    cut+=("_"+allSteps[k]);
    theFile->cd(cut.Data());

    hLeptLead[k]->Write();
    hLeptSubLead[k]->Write();
    hLeptThird[k]->Write();
    hSelectedIsoTracks[k]->Write();
    hSelectedIsoTrackThird[k]->Write();
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
    hMETLeadIsoTrackThird[k]->Write();

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

