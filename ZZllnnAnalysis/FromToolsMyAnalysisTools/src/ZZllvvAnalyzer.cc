/*
 *  See header file for a description of this class.
 *
 *  $Date: 2011/04/18 17:39:18 $
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
		    "JetVeto", "LeptonVeto", "RedMETcut"};
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
vector<HistoTrack*>      hSelectedTracks; 
vector<HistoTrack*>      hSelectedIsoTracks; 
vector<HistoDilept*>     hDileptKin;
vector<HistoKin*>        hJetKin;
vector<HistoKin*>        hSelJetKin;
vector<HistoObjectPair*> hDileptLeadJet;
vector<HistoObjectPair*> hDileptLeadSelJet;
vector<HistoKin*>        hMETKin;
vector<HistoObjectPair*> hDileptMET;
vector<HistoObjectPair*> hMETLeadJet;
vector<HistoObjectPair*> hMETLeadSelJet;
vector<HistoRedMET*>     hRedMetStd;

ReducedMETComputer *redMETComputer_std;

TString varsForLL[]={"run", "lumi", "event", "weight", 
		     "leadCharge", "subleadCharge", 
		     "leadPt", "subleadPt", 
		     "diLeptInvMass", 
		     "dileptLeadDeltaPhi",
		     "leptMinusCmCosTheta"};
const unsigned int nVarsOpt=sizeof(varsForLL)/sizeof(TString);

ZZllvvAnalyzer::ZZllvvAnalyzer(const ParameterSet& pSet) : totNEvents(0) 
{
  //theFile = new TFile(pSet.getUntrackedParameter<string>("fileName", "ZZllvvAnalyzer.root").c_str(),"RECREATE");
  theOutFileName = pSet.getUntrackedParameter<string>("fileName", "ZZllvvAnalyzer.root");
  source = pSet.getUntrackedParameter<InputTag>("source");
  zmmInput = pSet.getUntrackedParameter<InputTag>("zmmInput");
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
  if((*selectionInfo)[0] == 0) {
    if(debug)  cout << "\t No dilepton has been selected" << endl;
    return;
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

  if(debug) {
    cout << "\t dilepton leg1: type=" << (muonLead ? "MUON" : (electronLead ? "ELECTRON" : "UNKNOWN")) 
	 << "; charge=" << lep1->charge() << "; pt=" << lep1->pt() 
	 << "; eta=" << lep1->eta() << "; phi=" << lep1->phi() << endl;
    cout << "\t          leg2: type=" << (muonSubLead ? "MUON" : (electronSubLead ? "ELECTRON" : "UNKNOWN")) 
	 << "; charge=" << lep2->charge() << "; pt=" << lep2->pt() 
	 << "; eta=" << lep2->eta() << "; phi=" << lep2->phi() << endl;
  }


  // Skip unwanted flavor combinations
  if( flavorCombo==0 && (electronLead==0 || electronSubLead==0) && (muonLead==0 || muonSubLead==0) ) 
    return; 
  if( flavorCombo==1 && (electronLead==0 || electronSubLead==0) )
    return;
  if( flavorCombo==2 && (muonLead==0 || muonSubLead==0) ) 
    return; 
  if( flavorCombo==3 && (electronLead==0 || muonSubLead==0) && (muonLead==0 || electronSubLead==0) ) 
    return; 

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


  // Get the collection of selected, isolated muons
  vector<const pat::Muon *> allSelMuons;
  for(pat::eventhypothesis::Looper<pat::Muon> muon=h.loopAs<pat::Muon>("muon"); muon; ++muon) {
    if(debug) cout << "\t muon: " << muon->pt() << "; " << muon->eta() << "; " << muon->phi() << endl;
    const pat::Muon *tmpMu=dynamic_cast<const pat::Muon *>(muon.get());
    if(tmpMu) allSelMuons.push_back(tmpMu);
  }
  if(debug) cout << "# of selected, isolated muons: " << allSelMuons.size() << endl;

  // Get the collection of selected, isolated electrons
  vector<const pat::Electron *> allSelElectrons;
  for(pat::eventhypothesis::Looper<pat::Electron> ele=h.loopAs<pat::Electron>("electron"); ele; ++ele) {
    if(debug) cout << "\t electron: " << ele->pt() << "; " << ele->eta() << "; " << ele->phi() << endl;
    const pat::Electron *tmpEle=dynamic_cast<const pat::Electron *>(ele.get());
    if(tmpEle) allSelElectrons.push_back(tmpEle);
  }
  if(debug) cout << "# of selected, isolated electrons: " << allSelElectrons.size() << endl;

  // Get the collection of selected, isolated tracks
  edm::Handle<reco::MuonCollection> hTracks;
  reco::MuonRefVector allSelTracks, selIsoTracks;
  event.getByLabel(edm::InputTag(source.label()+":selectedTracks"), hTracks);
  for(reco::MuonCollection::size_type itTrk=0; itTrk<hTracks->size(); ++itTrk) {
    reco::MuonRef mRef(hTracks, itTrk);
    allSelTracks.push_back(mRef);
    if(debug) 
      cout << "\t selected track: " << mRef->pt() << "; " << mRef->eta() << "; " << mRef->phi() 
	   << "; Iso (nTracks, sumPt): " << mRef->isolationR03().nTracks << ", " << mRef->isolationR03().sumPt;
    if(mRef->innerTrack()->dz(vtx->position())<0.1) {
      if(debug) cout << "; from primary vertex. " << endl;
      if(mRef->isolationR03().sumPt/mRef->pt() < 0.05)
	selIsoTracks.push_back(mRef);
    }
    else 
      if(debug) cout << "; not from primary vertex. " << endl;
  }
  if(debug) cout << "# of selected tracks: " << allSelTracks.size() << endl;
  if(debug) cout << "# of selected, isolated tracks: " << selIsoTracks.size() << endl;

  int totNumberOfSelLeptons=allSelMuons.size()+allSelElectrons.size()+selIsoTracks.size();


  // =====================================================================
  //  ? veto on the third lepton
  // =====================================================================


  // Get the Jet momenta
  vector<CandidatePtr> allSelJets;
  vector<LorentzVector> selJetMomenta, moreSelJetMomenta;
  for (pat::eventhypothesis::Looper<pat::Jet> jet=h.loopAs<pat::Jet>("jet"); jet; ++jet) {
    if(debug) cout << "\t jet: " << jet->pt() << ";" << jet->eta() << ";" << jet->phi() << std::endl;
    selJetMomenta.push_back(jet->p4());
    if(jet->pt()>25.)
      moreSelJetMomenta.push_back(jet->p4());
  }
  if(debug) {
    cout << "# of selected jets: " << selJetMomenta.size() << endl;
    cout << "# of more selected jets: " << moreSelJetMomenta.size() << endl;
  }


  // =====================================================================
  //  ? veto on the number of jets (< 2)
  // =====================================================================


  // Get the MET
  const pat::MET *met = h.getAs<pat::MET>("met");
  if(debug) cout << "\t MET:" << met->pt() << ";" << met->phi() << endl;

  // Compute the redMET
  double lep1Err=( muonLead ? muonLead->track()->ptError() : 
		   electronLead->electronMomentumError()*sin(electronLead->theta()) );
  double lep2Err=( muonSubLead ? muonSubLead->track()->ptError() : 
		   electronSubLead->electronMomentumError()*sin(electronSubLead->theta()) );
  redMETComputer_std->compute(lep1->p4(), lep1Err, 
			      lep2->p4(), lep2Err, 
			      selJetMomenta,
			      met->p4());

  /*
  if(muonLead!=0 && muonSubLead!=0) {
    // Muon case
    redMETComputer_std->compute(muonLead->p4(), muonLead->track()->ptError(),
				muonSubLead->p4(), muonSubLead->track()->ptError(),
				selJetMomenta,
				met->p4());
  }
  else if(electronLead!=0 && electronSubLead!=0) {
    // Electron case
    // Nothing for the time being...
    return;
  }
  else {
    // Neither electrons nor Muons... You're done!
    // It should never happen at this point, but just in case...
    cout << " *** Neither an electron nor a muon valid pair was found! ***" << endl;
    return;
  }
  */

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
  //if((lep1->charge()*lep2->charge()>0) || ( (muonLead==0 || muonSubLead==0) && (electronLead==0 || electronSubLead==0) )) 
  if( chargeCombo!=0 && lep1->charge()*lep2->charge()!=chargeCombo ) {
    passedCuts[thiscut.Data()]=0;
    passedCuts["AllCutsPassed"]=0;
    //return;
  }
  //
  if(passedCuts["AllCutsPassed"]>0.5) 
    fillPlots(thiscut.Data(), vtx, totNvtx, goodNvtx, totNumberOfSelLeptons, lep1, lep2, allSelTracks, selIsoTracks, met, redMETComputer_std, selJetMomenta, moreSelJetMomenta);
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
    //return;
  }
  //
  if(passedCuts["AllCutsPassed"]>0.5) 
    fillPlots(thiscut.Data(), vtx, totNvtx, goodNvtx, totNumberOfSelLeptons, lep1, lep2, allSelTracks, selIsoTracks, met, redMETComputer_std, selJetMomenta, moreSelJetMomenta);
  //
  ////////////


  ////////////
  //
  // =====================================================================
  //  Cut2: veto on the number of jets (< 2)
  // =====================================================================
  //
  thiscut="JetVeto";
  if(selJetMomenta.size()>1) {
    passedCuts[thiscut.Data()]=0;
    passedCuts["AllCutsPassed"]=0;
    //return;
  }
  //
  if(passedCuts["AllCutsPassed"]>0.5) 
    fillPlots(thiscut.Data(), vtx, totNvtx, goodNvtx, totNumberOfSelLeptons, lep1, lep2, allSelTracks, selIsoTracks, met, redMETComputer_std, selJetMomenta, moreSelJetMomenta);
  //
  ////////////


  ////////////
  //
  // =====================================================================
  //  Cut3: veto on the third lepton
  // =====================================================================
  //
  thiscut="LeptonVeto";
  if(totNumberOfSelLeptons>0) {  // N.B. totNumberOfSelLeptons is the number of extra vertices, so must be = 0. 
    passedCuts[thiscut.Data()]=0;
    passedCuts["AllCutsPassed"]=0;
    //return;
  }
  //
  if(passedCuts["AllCutsPassed"]>0.5) 
    fillPlots(thiscut.Data(), vtx, totNvtx, goodNvtx, totNumberOfSelLeptons, lep1, lep2, allSelTracks, selIsoTracks, met, redMETComputer_std, selJetMomenta, moreSelJetMomenta);
  //
  ////////////


  ////////////
  //
  // =====================================================================
  //  Cut4: cut on the RedMET
  // =====================================================================
  //
  thiscut="RedMETcut";
  if(redMETComputer_std->reducedMET()<theRedMETMinCut) {
    passedCuts[thiscut.Data()]=0;
    passedCuts["AllCutsPassed"]=0;
    //return;
  }
  //
  if(passedCuts["AllCutsPassed"]>0.5) 
    fillPlots(thiscut.Data(), vtx, totNvtx, goodNvtx, totNumberOfSelLeptons, lep1, lep2, allSelTracks, selIsoTracks, met, redMETComputer_std, selJetMomenta, moreSelJetMomenta);
  //
  ////////////


  // ********************* //
  // ********************* //
  // **                 ** //
  // **  Fill the tree  ** //
  // **                 ** //
  // ********************* //
  // ********************* //

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
  tmpVars[tmpCnt++]=diLeptonMom.mass();
  tmpVars[tmpCnt++]=hDileptKin.back()->dileptLeadDeltaPhi;
  tmpVars[tmpCnt++]=hDileptKin.back()->leptMinusCmCosTheta;
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
    hSelectedTracks.push_back( new HistoTrack( ("SelectedTracks_"+cut).Data() ) );
    hSelectedIsoTracks.push_back( new HistoTrack( ("SelectedIsoTracks_"+cut).Data() ) );
    hDileptKin.push_back( new HistoDilept( ("DileptKin_"+cut).Data() ) );
    hJetKin.push_back( new HistoKin( ("JetKin_"+cut).Data() ) );
    hSelJetKin.push_back( new HistoKin( ("SelJetKin_"+cut).Data() ) );
    hDileptLeadJet.push_back( new HistoObjectPair( ("DileptLeadJet_"+cut).Data() ) );
    hDileptLeadSelJet.push_back( new HistoObjectPair( ("DileptLeadSelJet_"+cut).Data() ) );
    hMETKin.push_back( new HistoKin( ("METKin_"+cut).Data() ) );
    hDileptMET.push_back( new HistoObjectPair( ("DileptMET_"+cut).Data() ) );
    hMETLeadJet.push_back( new HistoObjectPair( ("METLeadJet_"+cut).Data() ) );
    hMETLeadSelJet.push_back( new HistoObjectPair( ("METLeadSelJet_"+cut).Data() ) );
    hRedMetStd.push_back( new HistoRedMET( ("RedMetStd_"+cut).Data() ) );
  }
  hEventCounter->Sumw2();
  hNVertexAll->Sumw2();
  hNVertexGood->Sumw2();
}


// void ZZllvvAnalyzer::fillPlots(string cutString, const reco::Vertex *pv, 
// 			       unsigned int nVtx, unsigned int nGoodVtx, int nLeptons, 
// 			       const pat::Muon *leadMu, const pat::Muon *subleadMu, 
// 			       reco::MuonRefVector &selectTrks, reco::MuonRefVector &selectIsoTrks, 
// 			       const pat::MET *missEt, ReducedMETComputer *redmet_std, 
// 			       vector<LorentzVector> &jetV, vector<LorentzVector> &selJetV) {
void ZZllvvAnalyzer::fillPlots(string cutString, const reco::Vertex *pv, 
			       unsigned int nVtx, unsigned int nGoodVtx, int nLeptons, 
			       const reco::CandidatePtr leadLept, const reco::CandidatePtr subleadLept, 
			       reco::MuonRefVector &selectTrks, reco::MuonRefVector &selectIsoTrks, 
			       const pat::MET *missEt, ReducedMETComputer *redmet_std, 
			       vector<LorentzVector> &jetV, vector<LorentzVector> &selJetV) {

  int cutN=maps[cutString];
  hEventCounter->Fill(cutN, theWeight);
  hNVertexAll->Fill(nVtx, theWeight);
  hNVertexGood->Fill(nGoodVtx, theWeight);
  LorentzVector dileptMom=leadLept->p4()+subleadLept->p4();
  hLeptLead[cutN]->FillNLept(nLeptons, theWeight);
  hLeptLead[cutN]->Fill(leadLept, pv, nVtx, theWeight);
  hLeptSubLead[cutN]->FillNLept(nLeptons, theWeight);
  hLeptSubLead[cutN]->Fill(subleadLept, pv, nVtx, theWeight);
  hSelectedTracks[cutN]->Fill(selectTrks, pv, nVtx, theWeight);
  hSelectedIsoTracks[cutN]->Fill(selectIsoTrks, pv, nVtx, theWeight);
  hDileptKin[cutN]->Fill(leadLept, subleadLept, nVtx, theWeight);
  hMETKin[cutN]->Fill(missEt->pt(), missEt->eta(), missEt->phi(), missEt->mass(), nVtx, theWeight);  
  hDileptMET[cutN]->Fill(dileptMom, missEt->p4(), theWeight);
  hRedMetStd[cutN]->Fill(redmet_std, missEt->pt(), nVtx, theWeight);
  hJetKin[cutN]->FillNObj(jetV.size(), theWeight);
  hSelJetKin[cutN]->FillNObj(selJetV.size(), theWeight);
  if(debug) {
    cout << " - Cut n. " << cutN << ": " << cutString.c_str() << endl;
    cout << " - Weight: " << theWeight << endl;
    cout << " - Lead Lept (pt, eta, phi) = (" 
	 << leadLept->pt() << ", " << leadLept->eta() << ", " << leadLept->phi() << ")" << endl;
    cout << " - Sublead Lept (pt, eta, phi) = (" 
	 << subleadLept->pt() << ", " << subleadLept->eta() << ", " << subleadLept->phi() << ")" << endl;
    cout << " - MET (pt, eta, phi, mass) = (" 
	 << missEt->pt() << ", " << missEt->eta() << ", " << missEt->phi() << ", " << missEt->mass() << ")" << endl;
    cout << " - RedMET [GeV] = " << missEt->pt() << endl;
    cout << " - Number of jets: " << jetV.size() << endl;
    cout << "                   (pt, eta, phi, mass)" << endl;
  }
  double leadJetPt=0;
  int leadJet=-1;
  for(unsigned int nj=0; nj<jetV.size(); ++nj) {
    hJetKin[cutN]->Fill(jetV[nj].pt(), jetV[nj].eta(), jetV[nj].phi(), jetV[nj].mass(), nVtx, theWeight);
    if(jetV[nj].pt()>leadJetPt) {
      leadJet=nj;
      leadJetPt=jetV[nj].pt();
    }
    if(debug) cout << "           jet #" << nj << ": (" << jetV[nj].pt() << ", " << jetV[nj].eta() 
		   << ", " << jetV[nj].phi() << ", " << jetV[nj].mass() << ")" << endl;
  }
  if(debug) {
    cout << " - Number of selected jets: " << selJetV.size() << endl;
    cout << "                   (pt, eta, phi, mass)" << endl;
  }
  for(unsigned int nj=0; nj<selJetV.size(); ++nj) {
    hSelJetKin[cutN]->Fill(selJetV[nj].pt(), selJetV[nj].eta(), selJetV[nj].phi(), selJetV[nj].mass(), nVtx, theWeight);
    if(debug) cout << "       sel jet #" << nj << ": (" << selJetV[nj].pt() << ", " << selJetV[nj].eta() 
  		   << ", " << selJetV[nj].phi() << ", " << selJetV[nj].mass() << ")" << endl;    
  }
  if(leadJet!=-1) {
    hMETLeadJet[cutN]->Fill(missEt->p4(), jetV[leadJet], theWeight);
    hDileptLeadJet[cutN]->Fill(dileptMom, jetV[leadJet], theWeight);
    if(selJetV.size()!=0) {
      hMETLeadSelJet[cutN]->Fill(missEt->p4(), jetV[leadJet], theWeight);
      hDileptLeadSelJet[cutN]->Fill(dileptMom, jetV[leadJet], theWeight);
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
    hSelectedTracks[k]->Scale(fact);
    hSelectedIsoTracks[k]->Scale(fact);
    hDileptKin[k]->Scale(fact);
    hJetKin[k]->Scale(fact);
    hSelJetKin[k]->Scale(fact);
    hDileptLeadJet[k]->Scale(fact);
    hDileptLeadSelJet[k]->Scale(fact);
    hMETKin[k]->Scale(fact);
    hDileptMET[k]->Scale(fact);
    hMETLeadJet[k]->Scale(fact);
    hMETLeadSelJet[k]->Scale(fact);
    hRedMetStd[k]->Scale(fact);
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
    hSelectedTracks[k]->Write();
    hSelectedIsoTracks[k]->Write();
    hDileptKin[k]->Write();
    hJetKin[k]->Write();
    hSelJetKin[k]->Write();
    hDileptLeadJet[k]->Write();
    hDileptLeadSelJet[k]->Write();
    hMETKin[k]->Write();
    hDileptMET[k]->Write();
    hMETLeadJet[k]->Write();
    hMETLeadSelJet[k]->Write();
    hRedMetStd[k]->Write();
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

