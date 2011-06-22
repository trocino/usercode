/*
 *  See header file for a description of this class.
 *
 *  $Date: 2011/04/18 17:39:18 $
 *  $Revision: 1.5 $
 *  \author G. Cerminara - CERN
 */
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "Tools/MyAnalysisTools/src/ZZllvvOptimizationAnalyzer.h"
#include "Tools/MyAnalysisTools/interface/Histograms.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/VertexReco/interface/Vertex.h"


#include "FWCore/Framework/interface/EventSetup.h"
// #include "AnalysisDataFormats/CMGTools/interface/Muon.h"
// #include "AnalysisDataFormats/CMGTools/interface/PFJet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/PatCandidates/interface/EventHypothesis.h"
#include "DataFormats/PatCandidates/interface/EventHypothesisLooper.h"
#include "DataFormats/Common/interface/MergeableCounter.h"
#include "TFile.h"

//#include "CMGTools/HtoZZ2l2nu/interface/Utils.h"
#include "CMGTools/HtoZZ2l2nu/interface/ReducedMETComputer.h"
#include "CMGTools/HtoZZ2l2nu/interface/ObjectFilters.h"

#include <iostream>
#include <stdio.h>

using namespace std;
using namespace edm;

int nEvtPreSkim;
int nEvtBaseFilter;
int nEvtSkim;
TH1F *hEvtCounter;

HistoLept *hMuonLead__cut0;    
HistoLept *hMuonSubLead__cut0;
HistoKin *hDiLeptKin__cut0;
HistoKin *hJetKin__cut0;
HistoKin *hMETKin__cut0;
vector<HistoRedMET*> hRedMetStd__cut0;

HistoLept *hMuonLead__cut1;    
HistoLept *hMuonSubLead__cut1;
HistoKin *hDiLeptKin__cut1;
HistoKin *hJetKin__cut1;
HistoKin *hMETKin__cut1;
vector<HistoRedMET*> hRedMetStd__cut1;

TH1F *hNVtxAll;

vector<ReducedMETComputer*> redMETComputer__std;
vector<string> kAllSteps;
unsigned int nCmb=0;

ZZllvvOptimizationAnalyzer::ZZllvvOptimizationAnalyzer(const ParameterSet& pSet) : totNEvents(0),
							   weight(1) {

  theFile = new TFile(pSet.getUntrackedParameter<string>("fileName", "ZZllvvOptimizationAnalyzer.root").c_str(),"RECREATE");
  source = pSet.getUntrackedParameter<InputTag>("source");
  zmmInput = pSet.getUntrackedParameter<InputTag>("zmmInput");
  debug = pSet.getUntrackedParameter<bool>("debug","False");
  nEvtPreSkim = 0;
  nEvtBaseFilter = 0;
  nEvtSkim = 0;
  vertexSelection =  pSet.getUntrackedParameter<ParameterSet>("Vertices");

  kRecoilSteps = pSet.getUntrackedParameter<vector<double> >("RecoilSteps", vector<double>());
  kSigmaPtSteps = pSet.getUntrackedParameter<vector<double> >("SigmaPtSteps", vector<double>());
  kPerpSteps = pSet.getUntrackedParameter<vector<double> >("PerpSteps", vector<double>());
  if(kRecoilSteps.size()==0) kRecoilSteps = vector<double>(1, 2.);
  if(kSigmaPtSteps.size()==0) kSigmaPtSteps = vector<double>(1, 2.5);
  if(kPerpSteps.size()==0) kPerpSteps = vector<double>(1, 1.);
  for(unsigned int ii=0; ii<kRecoilSteps.size(); ++ii) {
    for(unsigned int jj=0; jj<kSigmaPtSteps.size(); ++jj) {
      for(unsigned int kk=0; kk<kPerpSteps.size(); ++kk) {
	char tmpStr[50];
	sprintf(tmpStr, "%f_%f_%f", kRecoilSteps[ii], kSigmaPtSteps[jj], kPerpSteps[kk]);
	kAllSteps.push_back(string(tmpStr));
	// cout << nCmb << kAllSteps[nCmb].c_str() << endl;
	redMETComputer__std.push_back(new ReducedMETComputer( kRecoilSteps[ii],  kRecoilSteps[ii],
							     kSigmaPtSteps[jj], kSigmaPtSteps[jj],
							     kPerpSteps[kk]
							     ));
	++nCmb;
      }
    }
  }
}

ZZllvvOptimizationAnalyzer::~ZZllvvOptimizationAnalyzer(){
  cout << "destructor" << endl;
}

void ZZllvvOptimizationAnalyzer::beginJob() {
  
  
  theFile->cd();

  hEvtCounter = new TH1F("hEvtCounter", "Event counters", 10,0,10);
  hEvtCounter->GetXaxis()->SetBinLabel(1,"before skim");
  hEvtCounter->GetXaxis()->SetBinLabel(2,"after pre-filters");
  hEvtCounter->GetXaxis()->SetBinLabel(3,"after skim");
  hEvtCounter->GetXaxis()->SetBinLabel(4,"analyzed");
  hEvtCounter->GetXaxis()->SetBinLabel(5,"pre sel.");
  hEvtCounter->GetXaxis()->SetBinLabel(6,"sel dilepton");
  hEvtCounter->GetXaxis()->SetBinLabel(7,"mass window");

  hEvtCounter->Sumw2();

  // book the histograms
  hNVtxAll = new TH1F("hNVtxAll", "# of vertices", 100,0,100);
  hNVtxAll->Sumw2();

  hMuonLead__cut0    = new HistoLept("MuonLead__cut0");
  hMuonSubLead__cut0 = new HistoLept("MuonSubLead__cut0");
  hDiLeptKin__cut0   = new HistoKin("DiLeptKin__cut0");

  hMuonLead__cut1    = new HistoLept("MuonLead__cut1");
  hMuonSubLead__cut1 = new HistoLept("MuonSubLead__cut1");
  hDiLeptKin__cut1   = new HistoKin("DiLeptKin__cut1");

  hJetKin__cut0 = new HistoKin("JetKin__cut0");
  hMETKin__cut0 = new HistoKin("METKin__cut0");
  hJetKin__cut1 = new HistoKin("JetKin__cut1");
  hMETKin__cut1 = new HistoKin("METKin__cut1");
  //   hMETKin_J0__cut1 = new HistoKin("METKin_J0__cut1");
  //   hMETKin_J1__cut1 = new HistoKin("METKin_J1__cut1");

  for(unsigned int hh=0; hh<nCmb; ++hh) {
    // cout << " ciao: " << hh << " -- " 
    // << ("RedMetStd__"+kAllSteps[hh]+"___cut0").c_str() << " -- " 
    // << ("RedMetStd__"+kAllSteps[hh]+"___cut1").c_str();
    hRedMetStd__cut0.push_back(new HistoRedMET(("RedMetStd__"+kAllSteps[hh]+"___cut0").c_str()));
    hRedMetStd__cut1.push_back(new HistoRedMET(("RedMetStd__"+kAllSteps[hh]+"___cut1").c_str()));
    // cout << " ... Done!" << endl;
  }

}


// Operations
void ZZllvvOptimizationAnalyzer::beginRun(const Run& run, const EventSetup& eSetup) {
}
  
void ZZllvvOptimizationAnalyzer::analyze(const Event& event, const EventSetup& eSetup) {

  totNEvents++;
  // count all the events

  using reco::Candidate; 
  using reco::CandidatePtr;
  
  //pre-select vertices
  Handle<reco::VertexCollection> offlinePrimVertices;
  event.getByLabel("offlinePrimaryVertices", offlinePrimVertices);  
  std::vector<reco::VertexRef> selVertices = vertex::filter(offlinePrimVertices,vertexSelection);
  int nVert=selVertices.size();
  if(debug) cout << "# of vertices: " << nVert << endl;;



  //retrieve the event hypothesis
  Handle<vector<pat::EventHypothesis> > hyps;
  edm::InputTag selEvtTag(source.label()+":selectedEvent");
  event.getByLabel(selEvtTag, hyps);

  //retrieve the selection path
  Handle<std::vector<int> > selectionInfo;
  edm::InputTag selInfoTag(source.label()+":selectionInfo");
  event.getByLabel(selInfoTag, selectionInfo);

  //retrieve the selected vertex
  Handle<reco::VertexCollection> selectedVertex;
  edm::InputTag selVtxTag(source.label()+":selectedVertices");
  event.getByLabel(selVtxTag, selectedVertex);

  // Get the candidate collection
  Handle<View<reco::CompositeCandidate> > diMuonCands;
  event.getByLabel(string("zMMCand"), diMuonCands);  
  if(diMuonCands.isValid()) {
    if(debug) cout << "# of dimuons: " << diMuonCands->size() << endl;
  }

  // Get the candidate collection
  Handle<View<reco::CompositeCandidate> > diElectrCands;
  event.getByLabel("zEECand", diElectrCands);  
  
  if(debug) cout << "# of dielectrons: " << diElectrCands->size() << endl;
  

  // Added by Daniele
  if(diMuonCands->size()==0) return;

  if(hyps->size()==0) return;
  hEvtCounter->Fill(4);
  
  const pat::EventHypothesis &h = (*hyps)[0];

  //dump selected event 
  // SlectionInfo [0] -> path: {UNKNOWN=0,MUMU=1,EE=2,EMU=3};
  // SlectionInfo [1] -> step: 
  if(debug) {
    cout << "Retrieve an event hypothesis selected at step=" << (*selectionInfo)[1]
	 << " for path " << (*selectionInfo)[0] 
	 << " with " << selectedVertex->size() << " vertices" << std::endl;
  }

  if((*selectionInfo)[0] == 0) {
    if(debug)  cout << "\t no dilepton has been selected" << std::endl;
    return;
  }

  const pat::MET *met = h.getAs<pat::MET>("met");
  if(debug) cout << "\t met:" << met->pt() << ";" << met->phi() << endl;
  hMETKin__cut0->Fill(met->pt(), met->eta(), met->phi(), met->mass(), nVert, weight);

  vector<LorentzVector> jetMomenta;

  for (pat::eventhypothesis::Looper<pat::Jet> jet = h.loopAs<pat::Jet>("jet"); jet; ++jet) {
    if(debug) cout << "\t jet: " << jet->pt() << ";" << jet->eta() << ";" << jet->phi() << std::endl;
    jetMomenta.push_back(jet->p4());
    hJetKin__cut0->Fill(jet->pt(), jet->eta(), jet->phi(), jet->mass(), nVert, weight);
  }
  hJetKin__cut0->FillNObj(jetMomenta.size(), weight);
  if(debug) cout << "# of jets: " << jetMomenta.size() << endl;

  // =====================================================================
  // cut0: ask for the dilepton to exist
  
  hEvtCounter->Fill(5);

  CandidatePtr lep1 = h["leg1"];
  CandidatePtr lep2 = h["leg2"];
  if(debug) {
    cout << "\t dilepton leg1: " << lep1->pt() << ";" << lep1->eta() << ";" << lep1->phi() << endl
	 << "\t          leg2: " << lep2->pt() << ";" << lep2->eta() << ";" << lep2->phi() << endl;
  }
  const pat::Muon *muonLead    = dynamic_cast<const pat::Muon *>(lep1.get());
  const pat::Muon *muonSubLead = dynamic_cast<const pat::Muon *>(lep2.get());
  const reco::Vertex *vtx = &(*selectedVertex)[0];

  // Added by Daniele
  if(muonLead==0 || muonSubLead==0) return;

  typedef reco::Candidate::LorentzVector LorentzVector;
  LorentzVector leadLeptMom = lep1->p4();
  LorentzVector subLeadLeptMom = lep2->p4();
  LorentzVector diLeptonMom = leadLeptMom + subLeadLeptMom;


  hMuonLead__cut0->Fill(lep1, vtx, nVert, weight);
  hMuonSubLead__cut0->Fill(lep2, vtx, nVert, weight);
  hDiLeptKin__cut0->Fill(diLeptonMom.pt(), diLeptonMom.eta(), diLeptonMom.phi(), diLeptonMom.mass(), nVert, weight);

  

  for(unsigned int hh=0; hh<nCmb; ++hh) {
    redMETComputer__std[hh]->compute(muonLead->p4(), muonLead->track()->ptError(),
				muonSubLead->p4(), muonSubLead->track()->ptError(),
				jetMomenta,
				met->p4());
  
//     hRedMetStd__cut0[hh]->Fill(redMETComputer__std[hh]->reducedMET(),
// 			      redMETComputer__std[hh]->reducedMETComponents().first, redMETComputer__std[hh]->reducedMETComponents().second, 
// 			      redMETComputer__std[hh]->recoilProjComponents().first, redMETComputer__std[hh]->recoilProjComponents().second,
// 			      redMETComputer__std[hh]->metProjComponents().first, redMETComputer__std[hh]->metProjComponents().second,
// 			      redMETComputer__std[hh]->sumJetProjComponents().first, redMETComputer__std[hh]->sumJetProjComponents().second,
// 			      redMETComputer__std[hh]->dileptonProjComponents().first, redMETComputer__std[hh]->dileptonProjComponents().second,
// 			      redMETComputer__std[hh]->recoilType().first, redMETComputer__std[hh]->recoilType().second,
// 			      weight);
    hRedMetStd__cut0[hh]->Fill(redMETComputer__std[hh], met->pt(), offlinePrimVertices->size(), weight);
  }

  if(fabs(diLeptonMom.mass()-91.)>15.) return;
  
  // =====================================================================
  // cut1: apply mass window on the Z mass
  hEvtCounter->Fill(6);
  
  hMuonLead__cut1->Fill(lep1, vtx, nVert, weight);
  hMuonSubLead__cut1->Fill(lep2, vtx, nVert, weight);
  hDiLeptKin__cut1->Fill(diLeptonMom.pt(), diLeptonMom.eta(), diLeptonMom.phi(), diLeptonMom.mass(), nVert, weight);
  hMETKin__cut1->Fill(met->pt(), met->eta(), met->phi(), met->mass(), nVert, weight);
  for(unsigned int hh=0; hh<nCmb; ++hh) {
//     hRedMetStd__cut1[hh]->Fill(redMETComputer__std[hh]->reducedMET(),
// 			      redMETComputer__std[hh]->reducedMETComponents().first, redMETComputer__std[hh]->reducedMETComponents().second, 
// 			      redMETComputer__std[hh]->recoilProjComponents().first, redMETComputer__std[hh]->recoilProjComponents().second,
// 			      redMETComputer__std[hh]->metProjComponents().first, redMETComputer__std[hh]->metProjComponents().second,
// 			      redMETComputer__std[hh]->sumJetProjComponents().first, redMETComputer__std[hh]->sumJetProjComponents().second,
// 			      redMETComputer__std[hh]->dileptonProjComponents().first, redMETComputer__std[hh]->dileptonProjComponents().second,
// 			      redMETComputer__std[hh]->recoilType().first, redMETComputer__std[hh]->recoilType().second,
// 			      weight);
    hRedMetStd__cut1[hh]->Fill(redMETComputer__std[hh], met->pt(), offlinePrimVertices->size(), weight);
  }

  for (pat::eventhypothesis::Looper<pat::Jet> jet = h.loopAs<pat::Jet>("jet"); jet; ++jet) {
    if(debug) cout << "\t jet: " << jet->pt() << ";" << jet->eta() << ";" << jet->phi() << std::endl;
    hJetKin__cut1->Fill(jet->pt(), jet->eta(), jet->phi(), jet->mass(), nVert, weight);
  }
  hJetKin__cut1->FillNObj(jetMomenta.size(), weight);


  for (pat::eventhypothesis::Looper<pat::Electron> elec = h.loopAs<pat::Electron>("electron"); elec; ++elec) {
    if(debug) cout << "\t e: " << elec->pt() << ";" << elec->eta() << ";" << elec->phi() << std::endl;
  }

  for (pat::eventhypothesis::Looper<pat::Muon> muon = h.loopAs<pat::Muon>("muon"); muon; ++muon) {
    if(debug) cout << "\t mu: " << muon->pt() << ";" << muon->eta() << ";" << muon->phi() << std::endl;
  }


  //if event is MC debug gen level event
  if(!event.isRealData())
    {
      if(debug) cout << "\t Generator level event " << flush;
      int igenpart(0);
      for (pat::eventhypothesis::Looper<reco::GenParticle> genpart = h.loopAs<reco::GenParticle>("genparticle"); genpart; ++genpart) 
	{
	  if(debug) cout << "\t" << genpart->pdgId() << " -> " << flush;  

	  int igenpartdau(0);
	  char buf[20];
	  sprintf(buf,"gendaughter_%d",igenpart);
	  for(pat::eventhypothesis::Looper<reco::GenParticle> genpartdau = h.loopAs<reco::GenParticle>(buf); genpartdau; ++genpartdau)
	    {
	      if(debug) cout << genpartdau->pdgId() << " (" << flush;

	      char buf[20];
	      sprintf(buf,"gendaughter_%d_%d",igenpart,igenpartdau);
	      for(pat::eventhypothesis::Looper<reco::GenParticle> genpartgdau = h.loopAs<reco::GenParticle>(buf); genpartgdau; ++genpartgdau)
		if(debug) cout << genpartgdau->pdgId() << " " << flush;
	      
	      cout << ") " << flush;
	      igenpartdau++;
	    }
	  igenpart++;
	}
      if(debug) cout << endl;
    }
  
  



}


void ZZllvvOptimizationAnalyzer::endJob() {
  cout << "Tot. # of events pre skim: " << nEvtPreSkim << endl;
  cout << "Tot. # of events after base filters: " << nEvtBaseFilter << endl;
  cout << "Tot. # of events after skim: " << nEvtSkim << endl;
  cout << "Tot. # of events: " << totNEvents << endl;
  hEvtCounter->SetBinContent(1,nEvtPreSkim);
  hEvtCounter->SetBinContent(2,nEvtBaseFilter);
  hEvtCounter->SetBinContent(3,nEvtSkim);
  hEvtCounter->SetBinContent(4,totNEvents);

// //   // Write the histograms
  theFile->cd();
  hEvtCounter->Write();
  hNVtxAll->Write();
  hMuonLead__cut0->Write();    
  hMuonSubLead__cut0->Write();
  hDiLeptKin__cut0->Write();

  hMuonLead__cut1->Write();    
  hMuonSubLead__cut1->Write();
  hDiLeptKin__cut1->Write();
  hJetKin__cut0->Write();
  hMETKin__cut0->Write();
  hJetKin__cut1->Write();
  hMETKin__cut1->Write();
  for(unsigned int hh=0; hh<nCmb; ++hh) {
    TDirectory *subdir=theFile->mkdir(kAllSteps[hh].c_str());
    subdir->cd();
    hRedMetStd__cut0[hh]->Write();
    hRedMetStd__cut1[hh]->Write();
    theFile->cd();
  }

  theFile->Close();
}



void ZZllvvOptimizationAnalyzer::beginLuminosityBlock(const edm::LuminosityBlock & iLumi, const edm::EventSetup & iSetup) {

}

void ZZllvvOptimizationAnalyzer::endLuminosityBlock(const edm::LuminosityBlock & iLumi, const edm::EventSetup & iSetup) {
  edm::Handle<edm::MergeableCounter> nEvtPreSkim_lumi;
  iLumi.getByLabel("startCounter", nEvtPreSkim_lumi);
  nEvtPreSkim += nEvtPreSkim_lumi->value;

  edm::Handle<edm::MergeableCounter> nEvtBaseFilter_lumi;
  iLumi.getByLabel("preFilterCounter", nEvtBaseFilter_lumi);
  nEvtBaseFilter += nEvtBaseFilter_lumi->value;

  edm::Handle<edm::MergeableCounter> nEvtSkim_lumi;
  iLumi.getByLabel("mumuCounter", nEvtSkim_lumi);
  nEvtSkim += nEvtSkim_lumi->value;


}
