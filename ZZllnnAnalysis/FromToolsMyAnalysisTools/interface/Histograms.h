#ifndef Histograms_H
#define Histograms_H

/** \class Histograms
 *  No description available.
 *
 *  $Date: 2011/06/22 13:55:11 $
 *  $Revision: 1.2 $
 *  \author G. Cerminara - CERN
 */


#include <string>
#include <vector>
#include <iostream>

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
// #include "TProfile.h"
#include "TString.h"
#include "TMath.h"
#include "Math/GenVector/Boost.h"

#ifndef ROOTANALYSIS
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
//#include "PhysicsTools/CandUtils/interface/Booster.h"
#include "CMGTools/HtoZZ2l2nu/interface/Utils.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "CMGTools/HtoZZ2l2nu/interface/ReducedMETComputer.h"
#endif

using namespace std;


class HistoKin {
public:
  HistoKin(std::string name) : theName(name) {
    hPt=new TH1F(theName+"_hPt", "p_{T}", 100, 0. , 200.);
    hPtOverNVtx=new TH1F(theName+"_hPtOverNVtx", "p_{T}/n. vertices", 100, 0., 100.);
    hPtVsNVtx=new TH2F(theName+"_hPtVsNVtx", "p_{T} vs. n. vertices", 50, -0.5, 49.5, 100, 0., 200.);
    hPtOverNVtxVsNVtx=new TH2F(theName+"_hPtOverNVtxVsNVtx", "p_{T}/n. vertices vs. n. vertices", 50, -0.5, 49.5, 100, 0., 100.);
    hEta=new TH1F(theName+"_hEta","Eta",50,-5,5);
    hPhi=new TH1F(theName+"_hPhi","Phi (rad)",100, -3.15, 3.15);
    hMass=new TH1F(theName+"_hMass","Mass (GeV)", 200, 0., 500.);
    hNObj=new TH1F(theName+"_hNObj","# objects", 50, -0.5, 49.5);
    hPt->Sumw2();
    hEta->Sumw2();
    hPhi->Sumw2();
    hMass->Sumw2();
    hNObj->Sumw2();
  }


  HistoKin() : theName("") {
    hPt = 0;
    hPtOverNVtx = 0;
    hPtVsNVtx = 0;
    hPtOverNVtxVsNVtx = 0;
    hEta = 0;
    hPhi = 0;
    hMass = 0;
    hNObj = 0;
  }



  HistoKin(std::string name, TFile *file) : theName(name) {
    hPt = (TH1F *) file->Get(theName+"_hPt");
    hPtOverNVtx = (TH1F *) file->Get(theName+"_hPtOverNVtx");
    hPtVsNVtx = (TH2F *) file->Get(theName+"_hPtVsNVtx");
    hPtOverNVtxVsNVtx = (TH2F *) file->Get(theName+"_hPtOverNVtxVsNVtx");
    hEta = (TH1F *) file->Get(theName+"_hEta");
    hPhi = (TH1F *) file->Get(theName+"_hPhi");
    hMass = (TH1F *) file->Get(theName+"_hMass");
    hNObj = (TH1F *) file->Get(theName+"_hNObj");

  }


  HistoKin * Clone(std::string name) {
    HistoKin *ret = new HistoKin();
    ret->theName = name;

    if(hPt != 0) hPt->Clone((ret->theName+"_hPt").Data());
    if(hPtOverNVtx != 0) hPtOverNVtx->Clone((ret->theName+"_hPtOverNVtx").Data());
    if(hPtVsNVtx != 0)  hPtVsNVtx->Clone((ret->theName+"_hPtVsNVtx").Data());
    if(hPtOverNVtxVsNVtx != 0)  hPtOverNVtxVsNVtx->Clone((ret->theName+"_hPtOverNVtxVsNVtx").Data());
    if(hEta != 0) hEta->Clone((ret->theName+"_hEta").Data());
    if(hPhi != 0) hPhi->Clone((ret->theName+"_hPhi").Data());
    if(hMass != 0) hPhi->Clone((ret->theName+"_hMass").Data());
    if(hNObj != 0) hNObj->Clone((ret->theName+"_hNObj").Data());

    return ret;
  }



  void Add(const HistoKin* histSet) {
    if(hPt != 0) hPt->Add(histSet->hPt);
    if(hPtOverNVtx != 0) hPtOverNVtx->Add(hPtOverNVtx);
    if(hPtVsNVtx != 0) hPtVsNVtx->Add(hPtVsNVtx);
    if(hPtOverNVtxVsNVtx != 0) hPtOverNVtxVsNVtx->Add(hPtOverNVtxVsNVtx);
    if(hEta != 0) hEta->Add(histSet->hEta);
    if(hPhi != 0) hPhi->Add(histSet->hPhi);
    if(hMass != 0) hPhi->Add(histSet->hMass);
    if(hNObj != 0) hPhi->Add(histSet->hNObj);

  }



  void Scale(double scaleFact) {
    if(hPt != 0) hPt->Scale(scaleFact);
    if(hPtOverNVtx != 0) hPtOverNVtx->Scale(scaleFact);
    if(hPtVsNVtx != 0) hPtVsNVtx->Scale(scaleFact);
    if(hPtOverNVtxVsNVtx != 0) hPtOverNVtxVsNVtx->Scale(scaleFact);
    if(hEta != 0) hEta->Scale(scaleFact);
    if(hPhi != 0) hPhi->Scale(scaleFact);
    if(hMass != 0) hMass->Scale(scaleFact);
    if(hNObj != 0) hNObj->Scale(scaleFact);

  }


  void Write() {
    if(hPt != 0) hPt->Write();
    if(hPtOverNVtx != 0) hPtOverNVtx->Write();
    if(hPtVsNVtx != 0) hPtVsNVtx->Write();
    if(hPtOverNVtxVsNVtx != 0) hPtOverNVtxVsNVtx->Write();
    if(hEta != 0) hEta->Write();
    if(hPhi != 0) hPhi->Write();
    if(hMass != 0) hMass->Write();
    if(hNObj != 0) hNObj->Write();

  }
  
  void Fill(double pt, double eta, double phi, double mass, int nVtx, double weight) {
    hPt->Fill(pt, weight);
    if(nVtx!=0) hPtOverNVtx->Fill(pt/nVtx, weight);
    hPtVsNVtx->Fill(nVtx, pt, weight);
    if(nVtx!=0) hPtOverNVtxVsNVtx->Fill(nVtx, pt/nVtx, weight);
    hEta->Fill(eta, weight);
    hPhi->Fill(phi, weight);
    hMass->Fill(mass, weight);
  }

  void FillNObj(int nObj, double weight) {
    hNObj->Fill(nObj, weight);
  }


  /// Destructor
  virtual ~HistoKin() {}

  // Operations
  TString theName;

  TH1F *hPt;
  TH1F *hPtOverNVtx;
  TH2F *hPtVsNVtx;
  TH2F *hPtOverNVtxVsNVtx;
  TH1F *hEta;
  TH1F *hPhi;
  TH1F *hMass;
  TH1F *hNObj;
};




class HistoLept {
public:
  HistoLept(std::string name) : theName(name) {
    hKin = new HistoKin(name);
    hRelIso    = new TH1F(theName+"_hRelIso","Relative isolation",100,0,1);
    hDxy       = new TH1F(theName+"_hDxy","#Delta_{xy}", 100, 0., 0.1);
    hDz        = new TH1F(theName+"_hDz","#Delta_{z}", 100, 0., 0.5);
    hType      = new TH1F(theName+"_hType","Lepton type", 4, -1.5, 2.5);
    hNLept     = new TH1F(theName+"_hNLept","# leptons", 12, -0.5, 11.5);
    hNPixelHits = new TH1F(theName+"_hNPixelHits","# pixel hits", 10, -0.5, 9.5);
    hNTkHits    = new TH1F(theName+"_hNTkHits","# Tk hits", 40, -0.5, 39.5);
    hNMuonHits  = new TH1F(theName+"_hNMuonHits","# mu hits", 60, -0.5, 59.5);
    hNMatches  = new TH1F(theName+"_hNMatches","# matches", 15, -0.5, 14.5);


    hRelIso->Sumw2();
    hDxy->Sumw2();
    hDz->Sumw2();
    hType->Sumw2();
    hNLept->Sumw2();
    hNPixelHits->Sumw2();
    hNTkHits->Sumw2();
    hNMuonHits->Sumw2();
    hNMatches->Sumw2();

 }


  HistoLept() : theName("") {
    hKin = 0;
    hRelIso = 0;
    hDxy = 0;
    hDz = 0;
    hType = 0;
    hNLept = 0;
    hNPixelHits = 0;
    hNTkHits = 0;
    hNMuonHits = 0;
    hNMatches = 0;



  }


  HistoLept(std::string name, TFile *file) : theName(name) {
    hKin       = new HistoKin(name, file);
    hPt        = hKin->hPt;
    hEta       = hKin->hEta;
    hPhi       = hKin->hPhi;
    hMass      = hKin->hMass;
    hRelIso    = (TH1F *) file->Get(theName+"_hRelIso");
    hDxy       = (TH1F *) file->Get(theName+"_hDxy");
    hDz        = (TH1F *) file->Get(theName+"_hDz");
    hType      = (TH1F *) file->Get(theName+"_hType");
    hNLept     = (TH1F *) file->Get(theName+"_hNLept");
    hNPixelHits = (TH1F *) file->Get(theName+"_hNPixelHits");
    hNTkHits    = (TH1F *) file->Get(theName+"_hNTkHits");
    hNMuonHits  = (TH1F *) file->Get(theName+"_hNMuonHits");
    hNMatches   = (TH1F *) file->Get(theName+"_hNMatches");



 }


//   HistoLept * Clone(std::string name) {
//     HistoLept *ret = new HistoLept();
//     ret->theName = name;
//     hKin
//     if(hPt != 0) hPt->Clone((ret->theName+"_hPt").Data());
//     if(hEta != 0) hEta->Clone((ret->theName+"_hEta").Data());
//     if(hPhi != 0) hPhi->Clone((ret->theName+"_hPhi").Data());

//     return ret;
//   }



  void Add(const HistoLept* histSet) {
    if(hKin != 0) hKin->Add(histSet->hKin);
    if(hRelIso != 0) hRelIso->Add(histSet->hRelIso);
    if(hDxy != 0) hDxy->Add(histSet->hDxy);
    if(hDz != 0) hDz->Add(histSet->hDz);
    if(hType != 0) hType->Add(histSet->hType);
    if(hNLept != 0) hNLept->Add(histSet->hNLept);
    if(hNPixelHits != 0) hNPixelHits->Add(histSet->hNPixelHits);
    if(hNTkHits != 0) hNTkHits->Add(histSet->hNTkHits);
    if(hNMuonHits != 0) hNMuonHits->Add(histSet->hNMuonHits);
    if(hNMatches != 0) hNMatches->Add(histSet->hNMatches);


//     if(hPt != 0) hPt->Add(histSet->hPt);
//     if(hEta != 0) hEta->Add(histSet->hEta);
//     if(hPhi != 0) hPhi->Add(histSet->hPhi);
  }



  void Scale(double scaleFact) {
    if(hKin != 0) hKin->Scale(scaleFact);
    if(hRelIso != 0) hRelIso->Scale(scaleFact);
    if(hDxy != 0) hDxy->Scale(scaleFact);
    if(hDz != 0) hDz->Scale(scaleFact);
    if(hType != 0) hType->Scale(scaleFact);
    if(hNLept != 0) hNLept->Scale(scaleFact);
    if(hNPixelHits != 0) hNPixelHits->Scale(scaleFact);
    if(hNTkHits != 0) hNTkHits->Scale(scaleFact);
    if(hNMuonHits != 0) hNMuonHits->Scale(scaleFact);
    if(hNMatches != 0) hNMatches->Scale(scaleFact);


  }


  void Write() {
    if(hKin != 0) hKin->Write();
    if(hRelIso != 0) hRelIso->Write();
    if(hDxy != 0) hDxy->Write();
    if(hDz != 0) hDz->Write();
    if(hType != 0) hType->Write();
    if(hNLept != 0) hNLept->Write();
    if(hNPixelHits != 0) hNPixelHits->Write();
    if(hNTkHits != 0) hNTkHits->Write();
    if(hNMuonHits != 0) hNMuonHits->Write();
    if(hNMatches != 0) hNMatches->Write();


 }
  
  
#ifndef ROOTANALYSIS
  void Fill(const reco::CandidatePtr lept, const reco::Vertex *vtx, int nVtx, double weight) {

    const pat::Muon *muon=dynamic_cast<const pat::Muon*>(lept.get());
    const pat::Electron *ele=dynamic_cast<const pat::Electron*>(lept.get());
    if(muon==0 && ele==0) {
      cout << "[Histograms] *** Error: lepton is neither a muon nor an electron!" << endl;
      throw std::exception();
    }
    double relIso = (muon ? Utils::computeRelIsolation(muon) : Utils::computeRelIsolation(ele));
    int type = (muon ? Utils::muonType(muon) : Utils::muonType(ele));
    const reco::Track *track=(muon ? muon->track().get() : ele->gsfTrack().get());
    double dz=fabs(track->dz(vtx->position()));
    double dxy=fabs(track->dxy(vtx->position()));

    double nPxhits=track->hitPattern().numberOfValidPixelHits();
    double nTkhits=track->hitPattern().numberOfValidTrackerHits();
    double nMuonHits=(muon ? track->hitPattern().numberOfValidMuonHits() : 0);
    double nMatches=(muon ? muon->numberOfMatches() : 0);

    Fill(lept->pt(), lept->eta(), lept->phi(),
	 relIso,
	 dz, dxy,
	 type,
	 nPxhits, nTkhits, nMuonHits, nMatches,
	 nVtx, 
	 weight);

  }

#endif


  void Fill(double pt, double eta, double phi,
	    double relIso,
	    double dxy, double dz,
	    double type,
	    double nPxhits, double nTkhits, double nMuonHits, double nMatches,
	    int nVtx, 
	    double weight) {
    hKin->Fill(pt, eta, phi, 1, nVtx, weight);
    hRelIso->Fill(relIso, weight);
    hDxy->Fill(dxy, weight);
    hDz->Fill(dz, weight);
    hType->Fill(type, weight);
    hNPixelHits->Fill(nPxhits, weight);
    hNTkHits->Fill(nTkhits, weight);
    hNMuonHits->Fill(nMuonHits, weight);
    hNMatches->Fill(nMatches, weight);
  }

  void FillNLept(int num, int weight) {
    hNLept->Fill(num, weight);
  }

  /// Destructor
  virtual ~HistoLept() {}

  // Operations
  TString theName;

  HistoKin *hKin;
  TH1F *hPt;
  TH1F *hEta;
  TH1F *hPhi;
  TH1F *hMass;

  TH1F *hRelIso;
  TH1F *hDxy;
  TH1F *hDz;
  TH1F *hType;
  TH1F *hNLept;
  TH1F *hNPixelHits;
  TH1F *hNTkHits;
  TH1F *hNMuonHits;
  TH1F *hNMatches;

};


class HistoTrack {
public:
  HistoTrack(std::string name) : theName(name) {
    hKin         = new HistoKin(name);
    hIsoNTracks  = new TH1F(theName+"_hIsoNTracks","# tracks in iso cone", 50, -0.5, 49.5);
    hTrackIso    = new TH1F(theName+"_hTrackIso","Isolation", 100, 0., 20.);
    hTrackRelIso = new TH1F(theName+"_hTrackRelIso","Relative isolation", 100, 0., 0.1);
    hDxy         = new TH1F(theName+"_hDxy","#Delta_{xy}", 100, 0., 0.1);
    hDz          = new TH1F(theName+"_hDz","#Delta_{z}", 100, 0., 0.5);
    hNTrack      = new TH1F(theName+"_hNTrack","# tracks", 12, -0.5, 11.5);
    hNPixelHits  = new TH1F(theName+"_hNPixelHits","# pixel hits", 10, -0.5, 9.5);
    hNTkHits     = new TH1F(theName+"_hNTkHits","# Tk hits", 40, -0.5, 39.5);

    hIsoNTracks->Sumw2();
    hTrackIso->Sumw2();
    hTrackRelIso->Sumw2();
    hDxy->Sumw2();
    hDz->Sumw2();
    hNTrack->Sumw2();
    hNPixelHits->Sumw2();
    hNTkHits->Sumw2();
  }


  HistoTrack() : theName("") {
    hKin = 0;
    hIsoNTracks = 0;
    hTrackIso = 0;
    hTrackRelIso = 0;
    hDxy = 0;
    hDz = 0;
    hNTrack = 0;
    hNPixelHits = 0;
    hNTkHits = 0;
  }


  HistoTrack(std::string name, TFile *file) : theName(name) {
    hKin         = new HistoKin(name, file);
    hPt          = hKin->hPt;
    hEta         = hKin->hEta;
    hPhi         = hKin->hPhi;
    hMass        = hKin->hMass;
    hIsoNTracks  = (TH1F *) file->Get(theName+"_hIsoNTracks");
    hTrackIso    = (TH1F *) file->Get(theName+"_hTrackIso");
    hTrackRelIso = (TH1F *) file->Get(theName+"_hTrackRelIso");
    hDxy         = (TH1F *) file->Get(theName+"_hDxy");
    hDz          = (TH1F *) file->Get(theName+"_hDz");
    hNTrack      = (TH1F *) file->Get(theName+"_hNTrack");
    hNPixelHits  = (TH1F *) file->Get(theName+"_hNPixelHits");
    hNTkHits     = (TH1F *) file->Get(theName+"_hNTkHits");
  }


//   HistoTrack * Clone(std::string name) {
//     HistoTrack *ret = new HistoTrack();
//     ret->theName = name;
//     hKin
//     if(hPt != 0) hPt->Clone((ret->theName+"_hPt").Data());
//     if(hEta != 0) hEta->Clone((ret->theName+"_hEta").Data());
//     if(hPhi != 0) hPhi->Clone((ret->theName+"_hPhi").Data());

//     return ret;
//   }



  void Add(const HistoTrack* histSet) {
    if(hKin != 0)         hKin->Add(histSet->hKin);
    if(hIsoNTracks != 0)  hIsoNTracks->Add(histSet->hIsoNTracks);
    if(hTrackIso != 0)    hTrackIso->Add(histSet->hTrackIso);
    if(hTrackRelIso != 0) hTrackRelIso->Add(histSet->hTrackRelIso);
    if(hDxy != 0)         hDxy->Add(histSet->hDxy);
    if(hDz != 0)          hDz->Add(histSet->hDz);
    if(hNTrack != 0)      hNTrack->Add(histSet->hNTrack);
    if(hNPixelHits != 0)  hNPixelHits->Add(histSet->hNPixelHits);
    if(hNTkHits != 0)     hNTkHits->Add(histSet->hNTkHits);
  }



  void Scale(double scaleFact) {
    if(hKin != 0)         hKin->Scale(scaleFact);
    if(hIsoNTracks != 0)  hIsoNTracks->Scale(scaleFact);
    if(hTrackIso != 0)    hTrackIso->Scale(scaleFact);
    if(hTrackRelIso != 0) hTrackRelIso->Scale(scaleFact);
    if(hDxy != 0)         hDxy->Scale(scaleFact);
    if(hDz != 0)          hDz->Scale(scaleFact);
    if(hNTrack != 0)      hNTrack->Scale(scaleFact);
    if(hNPixelHits != 0)  hNPixelHits->Scale(scaleFact);
    if(hNTkHits != 0)     hNTkHits->Scale(scaleFact);
  }


  void Write() {
    if(hKin != 0)         hKin->Write();
    if(hIsoNTracks != 0)  hIsoNTracks->Write();
    if(hTrackIso != 0)    hTrackIso->Write();
    if(hTrackRelIso != 0) hTrackRelIso->Write();
    if(hDxy != 0)         hDxy->Write();
    if(hDz != 0)          hDz->Write();
    if(hNTrack != 0)      hNTrack->Write();
    if(hNPixelHits != 0)  hNPixelHits->Write();
    if(hNTkHits != 0)     hNTkHits->Write();
  }
  
  
#ifndef ROOTANALYSIS
  void Fill(reco::MuonRefVector &tracks, const reco::Vertex *vtx, int nVtx, double weight) {
    for(reco::muon_iterator it=tracks.begin(); it!=tracks.end(); ++it) {
      hKin->Fill((*it)->pt(), (*it)->eta(), (*it)->phi(), 1, nVtx, weight);
      hIsoNTracks->Fill( (*it)->isolationR03().nTracks, weight );
      hTrackIso->Fill( (*it)->isolationR03().sumPt, weight );
      hTrackRelIso->Fill( (*it)->isolationR03().sumPt/(*it)->pt(), weight );
      hDxy->Fill( fabs( (*it)->innerTrack()->dxy(vtx->position()) ), weight );
      hDz->Fill( fabs( (*it)->innerTrack()->dz(vtx->position()) ), weight );
      hNPixelHits->Fill( (*it)->innerTrack()->hitPattern().numberOfValidPixelHits(), weight );
      hNTkHits->Fill( (*it)->innerTrack()->hitPattern().numberOfValidTrackerHits(), weight );
    }
    FillNTrack(tracks.size(), weight);
  }
#endif

  void FillNTrack(int num, int weight) {
    hNTrack->Fill(num, weight);
  }

  /// Destructor
  virtual ~HistoTrack() {}

  // Operations
  TString theName;

  HistoKin *hKin;
  TH1F *hPt;
  TH1F *hEta;
  TH1F *hPhi;
  TH1F *hMass;

  TH1F *hIsoNTracks;
  TH1F *hTrackIso;
  TH1F *hTrackRelIso;
  TH1F *hDxy;
  TH1F *hDz;
  TH1F *hNTrack;
  TH1F *hNPixelHits;
  TH1F *hNTkHits;
};



class HistoObjectPair {
public:
  HistoObjectPair(std::string name) : theName(name) {
    hDeltaPhi    = new TH1F(theName+"_hDeltaPhi","#Delta#phi(1, 2)", 100, -3.15, 3.15);
    hDeltaEta    = new TH1F(theName+"_hDeltaEta","#Delta#eta(1, 2)", 100, -5., 5.);
    hDeltaTheta  = new TH1F(theName+"_hDeltaTheta","#Delta#theta(1, 2)", 100, -3.15, 3.15);
    hDeltaR      = new TH1F(theName+"_hDeltaR","#DeltaR(1, 1)", 100, 0., 6.);
    hAngle       = new TH1F(theName+"_hAngle","angle(1, 2)", 100, 0., 3.15);
    hCosAngle    = new TH1F(theName+"_hCosAngle","cosine of angle(1, 2)", 100, -1., 1.);
    hTransMass   = new TH1F(theName+"_hTransMass", "Transverse mass", 100, 0., 500.);
    hTransMassZZ = new TH1F(theName+"_hTransMassZZ", "ZZ transverse mass", 100, 0., 500.);

    hDeltaPhi->Sumw2();
    hDeltaEta->Sumw2();
    hDeltaTheta->Sumw2();
    hDeltaR->Sumw2();
    hAngle->Sumw2();
    hCosAngle->Sumw2();
    hTransMass->Sumw2();
    hTransMassZZ->Sumw2();
  }


  HistoObjectPair() : theName("") {
    hDeltaPhi    = 0;
    hDeltaEta    = 0;
    hDeltaTheta  = 0;
    hDeltaR      = 0;
    hAngle       = 0;
    hCosAngle    = 0;
    hTransMass   = 0;
    hTransMassZZ = 0;
  }


  HistoObjectPair(std::string name, TFile *file) : theName(name) {
    hDeltaPhi    = (TH1F *)file->Get(theName+"_hDeltaPhi");
    hDeltaEta    = (TH1F *)file->Get(theName+"_hDeltaEta");
    hDeltaTheta  = (TH1F *)file->Get(theName+"_hDeltaTheta");
    hDeltaR      = (TH1F *)file->Get(theName+"_hDeltaR");
    hAngle       = (TH1F *)file->Get(theName+"_hAngle");
    hCosAngle    = (TH1F *)file->Get(theName+"_hCosAngle");
    hTransMass   = (TH1F *)file->Get(theName+"_hTransMass");
    hTransMassZZ = (TH1F *)file->Get(theName+"_hTransMassZZ");
 }


  void Add(const HistoObjectPair* histSet) {
    if(hDeltaPhi    != 0) hDeltaPhi->Add(histSet->hDeltaPhi);
    if(hDeltaEta    != 0) hDeltaEta->Add(histSet->hDeltaEta);
    if(hDeltaTheta  != 0) hDeltaTheta->Add(histSet->hDeltaTheta);
    if(hDeltaR      != 0) hDeltaR->Add(histSet->hDeltaR);
    if(hAngle       != 0) hAngle->Add(histSet->hAngle);
    if(hCosAngle    != 0) hCosAngle->Add(histSet->hCosAngle);
    if(hTransMass   != 0) hTransMass->Add(histSet->hTransMass);
    if(hTransMassZZ != 0) hTransMassZZ->Add(histSet->hTransMassZZ);
  }

  void Scale(double scaleFact) {
    if(hDeltaPhi    != 0) hDeltaPhi->Scale(scaleFact);
    if(hDeltaEta    != 0) hDeltaEta->Scale(scaleFact);
    if(hDeltaTheta  != 0) hDeltaTheta->Scale(scaleFact);
    if(hDeltaR      != 0) hDeltaR->Scale(scaleFact);
    if(hAngle       != 0) hAngle->Scale(scaleFact);
    if(hCosAngle    != 0) hCosAngle->Scale(scaleFact);
    if(hTransMass   != 0) hTransMass->Scale(scaleFact);
    if(hTransMassZZ != 0) hTransMassZZ->Scale(scaleFact);
  }


  void Write() {
    if(hDeltaPhi    != 0) hDeltaPhi->Write();
    if(hDeltaEta    != 0) hDeltaEta->Write();
    if(hDeltaTheta  != 0) hDeltaTheta->Write();
    if(hDeltaR      != 0) hDeltaR->Write();
    if(hAngle       != 0) hAngle->Write();
    if(hCosAngle    != 0) hCosAngle->Write();
    if(hTransMass   != 0) hTransMass->Write();
    if(hTransMassZZ != 0) hTransMassZZ->Write();
  }
  
  
#ifndef ROOTANALYSIS
  void Fill(reco::Candidate::LorentzVector fm1, reco::Candidate::LorentzVector fm2, double weight) {
    reco::Candidate::Vector tm1(fm1.px(), fm1.py(), fm1.pz());
    reco::Candidate::Vector tm2(fm2.px(), fm2.py(), fm2.pz());

    hDeltaPhi->Fill( deltaPhi(tm1.phi(), tm2.phi()), weight );
    hDeltaEta->Fill( tm1.eta() - tm2.eta(), weight );
    hDeltaTheta->Fill( tm1.theta() - tm2.theta(), weight );
    hDeltaR->Fill( deltaR( tm1.eta(), tm1.phi(), tm2.eta(), tm2.phi() ), weight );
    hAngle->Fill( TMath::ACos( tm1.Unit().Dot(tm2.Unit()) ), weight );
    hCosAngle->Fill( tm1.Unit().Dot(tm2.Unit()), weight );
    double pt1=sqrt(tm1.perp2());
    double pt2=sqrt(tm2.perp2());
    double etTot=pt1+pt2;
    double ptTot=sqrt((tm1+tm2).perp2());
    double etZ1=sqrt( pt1*pt1 + 91.1876*91.1876 );
    double etZ2=sqrt( pt2*pt2 + 91.1876*91.1876 );
    double etZtot=etZ1+etZ2;
    hTransMass->Fill( sqrt(etTot*etTot - ptTot*ptTot), weight );
    hTransMassZZ->Fill( sqrt(etZtot*etZtot - ptTot*ptTot), weight );
  }
#endif

  /// Destructor
  virtual ~HistoObjectPair() {}

  // Operations
  TString theName;

  TH1F *hDeltaPhi;
  TH1F *hDeltaEta;
  TH1F *hDeltaTheta;
  TH1F *hDeltaR;
  TH1F *hAngle;
  TH1F *hCosAngle;
  TH1F *hTransMass;
  TH1F *hTransMassZZ; // Hypothesis: ZZ
};


class HistoDilept {
public:
  HistoDilept(std::string name) : theName(name) {

    hKin = new HistoKin(name);
    // Legend:
    //  D:  dilepton
    //  L:  leading
    //  S:  subleading
    //  P:  positive
    //  N:  negative
    //  Cm: center-of-mass
    hDeltaPhiLS  = new TH1F(theName+"_hDeltaPhiLS","#Delta#phi(lead, sublead)",      100, -3.15, 3.15);
    hDeltaPhiDL  = new TH1F(theName+"_hDeltaPhiDL","#Delta#phi(dilepton, lead)",     100, -3.15, 3.15);
    hDeltaPhiDS  = new TH1F(theName+"_hDeltaPhiDS","#Delta#phi(dilepton, sublead)",  100, -3.15, 3.15);
    hDeltaPhiDP  = new TH1F(theName+"_hDeltaPhiDP","#Delta#phi(dilepton, positive)", 100, -3.15, 3.15);
    hDeltaPhiDN  = new TH1F(theName+"_hDeltaPhiDN","#Delta#phi(dilepton, negative)", 100, -3.15, 3.15);

    hDeltaEtaLS  = new TH1F(theName+"_hDeltaEtaLS","#Delta#eta(lead, sublead)",      100, -5., 5.);
    hDeltaEtaDL  = new TH1F(theName+"_hDeltaEtaDL","#Delta#eta(dilepton, lead)",     100, -5., 5.);
    hDeltaEtaDS  = new TH1F(theName+"_hDeltaEtaDS","#Delta#eta(dilepton, sublead)",  100, -5., 5.);
    hDeltaEtaDP  = new TH1F(theName+"_hDeltaEtaDP","#Delta#eta(dilepton, positive)", 100, -5., 5.);
    hDeltaEtaDN  = new TH1F(theName+"_hDeltaEtaDN","#Delta#eta(dilepton, negative)", 100, -5., 5.);

    hDeltaRLS  = new TH1F(theName+"_hDeltaRLS","#DeltaR(lead, sublead)",      100, 0., 6.);
    hDeltaRDL  = new TH1F(theName+"_hDeltaRDL","#DeltaR(dilepton, lead)",     100, 0., 6.);
    hDeltaRDS  = new TH1F(theName+"_hDeltaRDS","#DeltaR(dilepton, sublead)",  100, 0., 6.);
    hDeltaRDP  = new TH1F(theName+"_hDeltaRDP","#DeltaR(dilepton, positive)", 100, 0., 6.);
    hDeltaRDN  = new TH1F(theName+"_hDeltaRDN","#DeltaR(dilepton, negative)", 100, 0., 6.);

    hCosAngleLS  = new TH1F(theName+"_hCosAngleLS","angle(lead, sublead)",      100, -1., 1.);
    hCosAngleDL  = new TH1F(theName+"_hCosAngleDL","angle(dilepton, lead)",     100, -1., 1.);
    hCosAngleDS  = new TH1F(theName+"_hCosAngleDS","angle(dilepton, sublead)",  100, -1., 1.);
    hCosAngleDP  = new TH1F(theName+"_hCosAngleDP","angle(dilepton, positive)", 100, -1., 1.);
    hCosAngleDN  = new TH1F(theName+"_hCosAngleDN","angle(dilepton, negative)", 100, -1., 1.);

    hCmCosThetaL  = new TH1F(theName+"_hCmCosThetaL","cos#theta* lead",     100, -1., 1.);
    hCmCosThetaS  = new TH1F(theName+"_hCmCosThetaS","cos#theta* sublead",  100, -1., 1.);
    hCmCosThetaP  = new TH1F(theName+"_hCmCosThetaP","cos#theta* positive", 100, -1., 1.);
    hCmCosThetaN  = new TH1F(theName+"_hCmCosThetaN","cos#theta* negative", 100, -1., 1.);

    hDeltaPhiLS->Sumw2();
    hDeltaPhiDL->Sumw2();
    hDeltaPhiDS->Sumw2();
    hDeltaPhiDP->Sumw2();
    hDeltaPhiDN->Sumw2();
    hDeltaEtaLS->Sumw2();
    hDeltaEtaDL->Sumw2();
    hDeltaEtaDS->Sumw2();
    hDeltaEtaDP->Sumw2();
    hDeltaEtaDN->Sumw2();
    hDeltaRLS->Sumw2();
    hDeltaRDL->Sumw2();
    hDeltaRDS->Sumw2();
    hDeltaRDP->Sumw2();
    hDeltaRDN->Sumw2();
    hCosAngleLS->Sumw2();
    hCosAngleDL->Sumw2();
    hCosAngleDS->Sumw2();
    hCosAngleDP->Sumw2();
    hCosAngleDN->Sumw2();
    hCmCosThetaL->Sumw2();
    hCmCosThetaS->Sumw2();
    hCmCosThetaP->Sumw2();
    hCmCosThetaN->Sumw2();

    // floats
    dileptLeadDeltaPhi=0.;
    leptMinusCmCosTheta=0.;
 }


  HistoDilept() : theName("") {
    hKin = 0;
    hDeltaPhiLS = 0;
    hDeltaPhiDL = 0;
    hDeltaPhiDS = 0;
    hDeltaPhiDP = 0;
    hDeltaPhiDN = 0;
    hDeltaEtaLS = 0;
    hDeltaEtaDL = 0;
    hDeltaEtaDS = 0;
    hDeltaEtaDP = 0;
    hDeltaEtaDN = 0;
    hDeltaRLS = 0;
    hDeltaRDL = 0;
    hDeltaRDS = 0;
    hDeltaRDP = 0;
    hDeltaRDN = 0;
    hCosAngleLS = 0;
    hCosAngleDL = 0;
    hCosAngleDS = 0;
    hCosAngleDP = 0;
    hCosAngleDN = 0;
    hCmCosThetaL = 0;
    hCmCosThetaS = 0;
    hCmCosThetaP = 0;
    hCmCosThetaN = 0;

    // floats
    dileptLeadDeltaPhi=0.;
    leptMinusCmCosTheta=0.;

  }


  HistoDilept(std::string name, TFile *file) : theName(name) {
    hKin       = new HistoKin(name, file);
    hPt          = hKin->hPt;
    hEta         = hKin->hEta;
    hPhi         = hKin->hPhi;
    hMass        = hKin->hMass;
    hDeltaPhiLS = (TH1F *)file->Get(theName+"_hDeltaPhiLS");
    hDeltaPhiDL = (TH1F *)file->Get(theName+"_hDeltaPhiDL");
    hDeltaPhiDS = (TH1F *)file->Get(theName+"_hDeltaPhiDS");
    hDeltaPhiDP = (TH1F *)file->Get(theName+"_hDeltaPhiDP");
    hDeltaPhiDN = (TH1F *)file->Get(theName+"_hDeltaPhiDN");
    hDeltaEtaLS = (TH1F *)file->Get(theName+"_hDeltaEtaLS");
    hDeltaEtaDL = (TH1F *)file->Get(theName+"_hDeltaEtaDL");
    hDeltaEtaDS = (TH1F *)file->Get(theName+"_hDeltaEtaDS");
    hDeltaEtaDP = (TH1F *)file->Get(theName+"_hDeltaEtaDP");
    hDeltaEtaDN = (TH1F *)file->Get(theName+"_hDeltaEtaDN");
    hDeltaRLS = (TH1F *)file->Get(theName+"_hDeltaRLS");
    hDeltaRDL = (TH1F *)file->Get(theName+"_hDeltaRDL");
    hDeltaRDS = (TH1F *)file->Get(theName+"_hDeltaRDS");
    hDeltaRDP = (TH1F *)file->Get(theName+"_hDeltaRDP");
    hDeltaRDN = (TH1F *)file->Get(theName+"_hDeltaRDN");
    hCosAngleLS = (TH1F *)file->Get(theName+"_hCosAngleLS");
    hCosAngleDL = (TH1F *)file->Get(theName+"_hCosAngleDL");
    hCosAngleDS = (TH1F *)file->Get(theName+"_hCosAngleDS");
    hCosAngleDP = (TH1F *)file->Get(theName+"_hCosAngleDP");
    hCosAngleDN = (TH1F *)file->Get(theName+"_hCosAngleDN");
    hCmCosThetaL = (TH1F *)file->Get(theName+"_hCmCosThetaL");
    hCmCosThetaS = (TH1F *)file->Get(theName+"_hCmCosThetaS");
    hCmCosThetaP = (TH1F *)file->Get(theName+"_hCmCosThetaP");
    hCmCosThetaN = (TH1F *)file->Get(theName+"_hCmCosThetaN");
 }


  void Add(const HistoDilept* histSet) {
    if(hKin != 0) hKin->Add(histSet->hKin);
    if(hDeltaPhiLS != 0) hDeltaPhiLS->Add(histSet->hDeltaPhiLS);
    if(hDeltaPhiDL != 0) hDeltaPhiDL->Add(histSet->hDeltaPhiDL);
    if(hDeltaPhiDS != 0) hDeltaPhiDS->Add(histSet->hDeltaPhiDS);
    if(hDeltaPhiDP != 0) hDeltaPhiDP->Add(histSet->hDeltaPhiDP);
    if(hDeltaPhiDN != 0) hDeltaPhiDN->Add(histSet->hDeltaPhiDN);
    if(hDeltaEtaLS != 0) hDeltaEtaLS->Add(histSet->hDeltaEtaLS);
    if(hDeltaEtaDL != 0) hDeltaEtaDL->Add(histSet->hDeltaEtaDL);
    if(hDeltaEtaDS != 0) hDeltaEtaDS->Add(histSet->hDeltaEtaDS);
    if(hDeltaEtaDP != 0) hDeltaEtaDP->Add(histSet->hDeltaEtaDP);
    if(hDeltaEtaDN != 0) hDeltaEtaDN->Add(histSet->hDeltaEtaDN);
    if(hDeltaRLS != 0) hDeltaRLS->Add(histSet->hDeltaRLS);
    if(hDeltaRDL != 0) hDeltaRDL->Add(histSet->hDeltaRDL);
    if(hDeltaRDS != 0) hDeltaRDS->Add(histSet->hDeltaRDS);
    if(hDeltaRDP != 0) hDeltaRDP->Add(histSet->hDeltaRDP);
    if(hDeltaRDN != 0) hDeltaRDN->Add(histSet->hDeltaRDN);
    if(hCosAngleLS != 0) hCosAngleLS->Add(histSet->hCosAngleLS);
    if(hCosAngleDL != 0) hCosAngleDL->Add(histSet->hCosAngleDL);
    if(hCosAngleDS != 0) hCosAngleDS->Add(histSet->hCosAngleDS);
    if(hCosAngleDP != 0) hCosAngleDP->Add(histSet->hCosAngleDP);
    if(hCosAngleDN != 0) hCosAngleDN->Add(histSet->hCosAngleDN);
    if(hCmCosThetaL != 0) hCmCosThetaL->Add(histSet->hCmCosThetaL);
    if(hCmCosThetaS != 0) hCmCosThetaS->Add(histSet->hCmCosThetaS);
    if(hCmCosThetaP != 0) hCmCosThetaP->Add(histSet->hCmCosThetaP);
    if(hCmCosThetaN != 0) hCmCosThetaN->Add(histSet->hCmCosThetaN);
  }



  void Scale(double scaleFact) {
    if(hKin != 0) hKin->Scale(scaleFact);
    if(hDeltaPhiLS != 0) hDeltaPhiLS->Scale(scaleFact);
    if(hDeltaPhiDL != 0) hDeltaPhiDL->Scale(scaleFact);
    if(hDeltaPhiDS != 0) hDeltaPhiDS->Scale(scaleFact);
    if(hDeltaPhiDP != 0) hDeltaPhiDP->Scale(scaleFact);
    if(hDeltaPhiDN != 0) hDeltaPhiDN->Scale(scaleFact);
    if(hDeltaEtaLS != 0) hDeltaEtaLS->Scale(scaleFact);
    if(hDeltaEtaDL != 0) hDeltaEtaDL->Scale(scaleFact);
    if(hDeltaEtaDS != 0) hDeltaEtaDS->Scale(scaleFact);
    if(hDeltaEtaDP != 0) hDeltaEtaDP->Scale(scaleFact);
    if(hDeltaEtaDN != 0) hDeltaEtaDN->Scale(scaleFact);
    if(hDeltaRLS != 0) hDeltaRLS->Scale(scaleFact);
    if(hDeltaRDL != 0) hDeltaRDL->Scale(scaleFact);
    if(hDeltaRDS != 0) hDeltaRDS->Scale(scaleFact);
    if(hDeltaRDP != 0) hDeltaRDP->Scale(scaleFact);
    if(hDeltaRDN != 0) hDeltaRDN->Scale(scaleFact);
    if(hCosAngleLS != 0) hCosAngleLS->Scale(scaleFact);
    if(hCosAngleDL != 0) hCosAngleDL->Scale(scaleFact);
    if(hCosAngleDS != 0) hCosAngleDS->Scale(scaleFact);
    if(hCosAngleDP != 0) hCosAngleDP->Scale(scaleFact);
    if(hCosAngleDN != 0) hCosAngleDN->Scale(scaleFact);
    if(hCmCosThetaL != 0) hCmCosThetaL->Scale(scaleFact);
    if(hCmCosThetaS != 0) hCmCosThetaS->Scale(scaleFact);
    if(hCmCosThetaP != 0) hCmCosThetaP->Scale(scaleFact);
    if(hCmCosThetaN != 0) hCmCosThetaN->Scale(scaleFact);
  }


  void Write() {
    if(hKin != 0) hKin->Write();
    if(hDeltaPhiLS != 0) hDeltaPhiLS->Write();
    if(hDeltaPhiDL != 0) hDeltaPhiDL->Write();
    if(hDeltaPhiDS != 0) hDeltaPhiDS->Write();
    if(hDeltaPhiDP != 0) hDeltaPhiDP->Write();
    if(hDeltaPhiDN != 0) hDeltaPhiDN->Write();
    if(hDeltaEtaLS != 0) hDeltaEtaLS->Write();
    if(hDeltaEtaDL != 0) hDeltaEtaDL->Write();
    if(hDeltaEtaDS != 0) hDeltaEtaDS->Write();
    if(hDeltaEtaDP != 0) hDeltaEtaDP->Write();
    if(hDeltaEtaDN != 0) hDeltaEtaDN->Write();
    if(hDeltaRLS != 0) hDeltaRLS->Write();
    if(hDeltaRDL != 0) hDeltaRDL->Write();
    if(hDeltaRDS != 0) hDeltaRDS->Write();
    if(hDeltaRDP != 0) hDeltaRDP->Write();
    if(hDeltaRDN != 0) hDeltaRDN->Write();
    if(hCosAngleLS != 0) hCosAngleLS->Write();
    if(hCosAngleDL != 0) hCosAngleDL->Write();
    if(hCosAngleDS != 0) hCosAngleDS->Write();
    if(hCosAngleDP != 0) hCosAngleDP->Write();
    if(hCosAngleDN != 0) hCosAngleDN->Write();
    if(hCmCosThetaL != 0) hCmCosThetaL->Write();
    if(hCmCosThetaS != 0) hCmCosThetaS->Write();
    if(hCmCosThetaP != 0) hCmCosThetaP->Write();
    if(hCmCosThetaN != 0) hCmCosThetaN->Write();
 }
  
  
#ifndef ROOTANALYSIS
  void Fill(const reco::CandidatePtr llep, const reco::CandidatePtr slep, int nVtx, double weight) {

    reco::Candidate::LorentzVector lfm=llep->p4();
    reco::Candidate::LorentzVector sfm=slep->p4();
    reco::Candidate::LorentzVector dfm=lfm+sfm;

    hKin->Fill(dfm.pt(), dfm.eta(), dfm.phi(), dfm.mass(), nVtx, weight );

    reco::Candidate::Vector ltm(lfm.px(), lfm.py(), lfm.pz());
    reco::Candidate::Vector stm(sfm.px(), sfm.py(), sfm.pz());
    reco::Candidate::Vector dtm=ltm+stm;

    reco::Candidate::Vector betaVect(dfm.BoostToCM());
    ROOT::Math::Boost booster(betaVect);
    reco::Candidate::LorentzVector lcfm=booster(lfm);
    reco::Candidate::LorentzVector scfm=booster(sfm);
    reco::Candidate::Vector lctm(lcfm.px(), lcfm.py(), lcfm.pz());
    reco::Candidate::Vector sctm(scfm.px(), scfm.py(), scfm.pz());
    // Booster booster(betaVect);
    // pat::Muon *lcmu=(pat::Muon *)lmu->clone();
    // pat::Muon *scmu=(pat::Muon *)smu->clone();
    // booster.set(*lcmu);
    // booster.set(*scmu);
    // reco::Candidate::Vector lctm(lcmu->px(), lcmu->py(), lcmu->pz());
    // reco::Candidate::Vector sctm(scmu->px(), scmu->py(), scmu->pz());

    reco::Candidate::Vector ptm, ntm, pctm, nctm;
    if(llep->charge()>0) {ptm=ltm; ntm=stm; pctm=lctm; nctm=sctm;}
    else {ptm=stm; ntm=ltm; pctm=sctm; nctm=lctm;}

    hDeltaPhiLS->Fill( deltaPhi(ltm.phi(), stm.phi()), weight );
    hDeltaPhiDL->Fill( deltaPhi(dtm.phi(), ltm.phi()), weight );
    hDeltaPhiDS->Fill( deltaPhi(dtm.phi(), stm.phi()), weight );
    hDeltaPhiDP->Fill( deltaPhi(dtm.phi(), ptm.phi()), weight );
    hDeltaPhiDN->Fill( deltaPhi(dtm.phi(), ntm.phi()), weight );
    hDeltaEtaLS->Fill( ltm.eta() - stm.eta(), weight );
    hDeltaEtaDL->Fill( dtm.eta() - ltm.eta(), weight );
    hDeltaEtaDS->Fill( dtm.eta() - stm.eta(), weight );
    hDeltaEtaDP->Fill( dtm.eta() - ptm.eta(), weight );
    hDeltaEtaDN->Fill( dtm.eta() - ntm.eta(), weight );
    hDeltaRLS->Fill( deltaR( ltm.eta(), ltm.phi(), stm.eta(), stm.phi() ), weight );
    hDeltaRDL->Fill( deltaR( dtm.eta(), dtm.phi(), ltm.eta(), ltm.phi() ), weight );
    hDeltaRDS->Fill( deltaR( dtm.eta(), dtm.phi(), stm.eta(), stm.phi() ), weight );
    hDeltaRDP->Fill( deltaR( dtm.eta(), dtm.phi(), ptm.eta(), ptm.phi() ), weight );
    hDeltaRDN->Fill( deltaR( dtm.eta(), dtm.phi(), ntm.eta(), ntm.phi() ), weight );

    hCosAngleLS->Fill( ltm.Unit().Dot(stm.Unit()), weight );
    hCosAngleDL->Fill( dtm.Unit().Dot(ltm.Unit()), weight );
    hCosAngleDS->Fill( dtm.Unit().Dot(stm.Unit()), weight );
    hCosAngleDP->Fill( dtm.Unit().Dot(ptm.Unit()), weight );
    hCosAngleDN->Fill( dtm.Unit().Dot(ntm.Unit()), weight );

    hCmCosThetaL->Fill(lctm.Unit().z(), weight );
    hCmCosThetaS->Fill(sctm.Unit().z(), weight );
    hCmCosThetaP->Fill(pctm.Unit().z(), weight );
    hCmCosThetaN->Fill(nctm.Unit().z(), weight );

    // Some (public) data member for easy access
    dileptLeadDeltaPhi=deltaPhi(dtm.phi(), ltm.phi());
    leptMinusCmCosTheta=nctm.Unit().z();
  }
#endif

  /// Destructor
  virtual ~HistoDilept() {}

  // Operations
  TString theName;

  HistoKin *hKin;
  TH1F *hPt;
  TH1F *hEta;
  TH1F *hPhi;
  TH1F *hMass;
  TH1F *hDeltaPhiLS;
  TH1F *hDeltaPhiDL;
  TH1F *hDeltaPhiDS;
  TH1F *hDeltaPhiDP;
  TH1F *hDeltaPhiDN;
  TH1F *hDeltaEtaLS;
  TH1F *hDeltaEtaDL;
  TH1F *hDeltaEtaDS;
  TH1F *hDeltaEtaDP;
  TH1F *hDeltaEtaDN;
  TH1F *hDeltaRLS;
  TH1F *hDeltaRDL;
  TH1F *hDeltaRDS;
  TH1F *hDeltaRDP;
  TH1F *hDeltaRDN;
  TH1F *hCosAngleLS;
  TH1F *hCosAngleDL;
  TH1F *hCosAngleDS;
  TH1F *hCosAngleDP;
  TH1F *hCosAngleDN;
  TH1F *hCmCosThetaL;
  TH1F *hCmCosThetaS;
  TH1F *hCmCosThetaP;
  TH1F *hCmCosThetaN;

  // Some (public) data member for easy access
  float dileptLeadDeltaPhi;
  float leptMinusCmCosTheta;

};


class HistoRedMET {
public:
  HistoRedMET(std::string name) : theName(name) {
    hRedMET         = new TH1F(theName+"_hRedMET",         "title", 100, 0., 100.);
    hRedMETCompLong = new TH1F(theName+"_hRedMETCompLong", "title", 100, 0., 100.);
    hRedMETCompPerp = new TH1F(theName+"_hRedMETCompPerp", "title", 100, 0., 100.);

    hRedMETOverNVtx         = new TH1F(theName+"_hRedMETOverNVtx",         "title", 100, 0., 50.);
    hRedMETCompLongOverNVtx = new TH1F(theName+"_hRedMETCompLongOverNVtx", "title", 100, 0., 50.);
    hRedMETCompPerpOverNVtx = new TH1F(theName+"_hRedMETCompPerpOverNVtx", "title", 100, 0., 50.);

    hRedMETVsNVtx         = new TH2F(theName+"_hRedMETVsNVtx",         "title", 50, -0.5, 49.5, 100, 0., 100.);
    hRedMETCompLongVsNVtx = new TH2F(theName+"_hRedMETCompLongVsNVtx", "title", 50, -0.5, 49.5, 100, 0., 100.);
    hRedMETCompPerpVsNVtx = new TH2F(theName+"_hRedMETCompPerpVsNVtx", "title", 50, -0.5, 49.5, 100, 0., 100.);

    hRedMETOverNVtxVsNVtx         = new TH2F(theName+"_hRedMETOverNVtxVsNVtx",         "title", 50, -0.5, 49.5, 100, 0., 50.);
    hRedMETCompLongOverNVtxVsNVtx = new TH2F(theName+"_hRedMETCompLongOverNVtxVsNVtx", "title", 50, -0.5, 49.5, 100, 0., 50.);
    hRedMETCompPerpOverNVtxVsNVtx = new TH2F(theName+"_hRedMETCompPerpOverNVtxVsNVtx", "title", 50, -0.5, 49.5, 100, 0., 50.);

    hRecoilCompLong = new TH1F(theName+ "_hRecoilCompLong","title",150,-150,150);
    hRecoilCompPerp = new TH1F(theName+ "_hRecoilCompPerp","title",150,-150,150);
    hMetCompLong = new TH1F(theName+ "_hMetCompLong","title",150,-150,150);
    hMetCompPerp = new TH1F(theName+ "_hMetCompPerp","title",150,-150,150);
    hSumJetCompLong = new TH1F(theName+ "_hSumJetCompLong","title",150,-150,150);
    hSumJetCompPerp = new TH1F(theName+ "_hSumJetCompPerp","title",150,-150,150);
    hDileptonCompLong = new TH1F(theName+ "_hDileptonCompLong","title",150,0,150);
    hDileptonCompPerp = new TH1F(theName+ "_hDileptonCompPerp","title",150,0,150);

    hDileptSigmaPtCompLong  = new TH1F(theName+ "_hDileptSigmaPtCompLong","title",150,0,150);
    hDileptSigmaPtCompPerp = new TH1F(theName+ "_hDileptSigmaPtCompPerp","title",150,0,150);

    hDileptSigmaPtCompLong  = new TH1F(theName+ "_hDileptSigmaPtCompLong","title",100,-50,0);
    hDileptSigmaPtCompPerp = new TH1F(theName+ "_hDileptSigmaPtCompPerp","title",100,-50,0.);

    hDileptSigmaPtCompLongVsMET  = new TH2F(theName+ "_hDileptSigmaPtCompLongVsMET","title",100,0.,200.,50,-50.,0.);
    hDileptSigmaPtCompPerpVsMET = new TH2F(theName+ "_hDileptSigmaPtCompPerpVsMET","title",100,0.,200.,50,-50,0.);
    hDileptonVsRecoilCompLong = new TH2F(theName+ "_hDileptonVsRecoilCompLong","title",150,-150.,0.,150,0.,150.);
    hDileptonVsRecoilCompPerp = new TH2F(theName+ "_hDileptonVsRecoilCompPerp","title",150,-150.,0.,150,0.,150.);

    hRecoilTypeLong = new TH1F(theName+ "_hRecoilTypeLong","Type",2,0,2);
    hRecoilTypeLong->GetXaxis()->SetBinLabel(1,"Jet-Like");
    hRecoilTypeLong->GetXaxis()->SetBinLabel(2,"MET-Like");
    hRecoilTypePerp = new TH1F(theName+ "_hRecoilTypePerp","Type",2,0,2);
    hRecoilTypePerp->GetXaxis()->SetBinLabel(1,"Jet-Like");
    hRecoilTypePerp->GetXaxis()->SetBinLabel(2,"MET-Like");


    hRedMET->Sumw2();
    hRedMETCompLong->Sumw2();
    hRedMETCompPerp->Sumw2();
    hRedMETOverNVtx->Sumw2();
    hRedMETCompLongOverNVtx->Sumw2();
    hRedMETCompPerpOverNVtx->Sumw2();
    hRedMETVsNVtx->Sumw2();
    hRedMETCompLongVsNVtx->Sumw2();
    hRedMETCompPerpVsNVtx->Sumw2();
    hRedMETOverNVtxVsNVtx->Sumw2();
    hRedMETCompLongOverNVtxVsNVtx->Sumw2();
    hRedMETCompPerpOverNVtxVsNVtx->Sumw2();
    hRecoilCompLong->Sumw2();
    hRecoilCompPerp->Sumw2();
    hMetCompLong->Sumw2();
    hMetCompPerp->Sumw2();
    hSumJetCompLong->Sumw2();
    hSumJetCompPerp->Sumw2();
    hDileptonCompLong->Sumw2();
    hDileptonCompPerp->Sumw2();
    hDileptSigmaPtCompLong->Sumw2();
    hDileptSigmaPtCompPerp->Sumw2();
    hDileptSigmaPtCompLongVsMET->Sumw2();
    hDileptSigmaPtCompPerpVsMET->Sumw2();
    hDileptonVsRecoilCompLong->Sumw2();
    hDileptonVsRecoilCompPerp->Sumw2();

    hRecoilTypeLong->Sumw2();
    hRecoilTypePerp->Sumw2();

  }


  HistoRedMET() : theName("") {
    hRedMET = 0;
    hRedMETCompLong = 0;
    hRedMETCompPerp = 0;
    hRedMETOverNVtx = 0;
    hRedMETCompLongOverNVtx = 0;
    hRedMETCompPerpOverNVtx = 0;
    hRedMETVsNVtx = 0;
    hRedMETCompLongVsNVtx = 0;
    hRedMETCompPerpVsNVtx = 0;
    hRedMETOverNVtxVsNVtx = 0;
    hRedMETCompLongOverNVtxVsNVtx = 0;
    hRedMETCompPerpOverNVtxVsNVtx = 0;
    hRecoilCompLong = 0;
    hRecoilCompPerp = 0;
    hMetCompLong = 0;
    hMetCompPerp = 0;
    hSumJetCompLong = 0;
    hSumJetCompPerp = 0;
    hDileptonCompLong = 0;
    hDileptonCompPerp = 0;
    hRecoilTypeLong = 0;
    hRecoilTypePerp = 0;

    hDileptSigmaPtCompLong = 0;
    hDileptSigmaPtCompPerp = 0;
    hDileptSigmaPtCompLongVsMET = 0;
    hDileptSigmaPtCompPerpVsMET = 0;
    hDileptonVsRecoilCompLong = 0;
    hDileptonVsRecoilCompPerp = 0;

  }



  HistoRedMET(std::string name, TFile *file) : theName(name) {
    hRedMET         = (TH1F *) file->Get(theName+"_hRedMET");
    hRedMETCompLong = (TH1F *) file->Get(theName+"_hRedMETCompLong");
    hRedMETCompPerp = (TH1F *) file->Get(theName+"_hRedMETCompPerp");
    hRedMETOverNVtx         = (TH1F *) file->Get(theName+"_hRedMETOverNVtx");	     
    hRedMETCompLongOverNVtx = (TH1F *) file->Get(theName+"_hRedMETCompLongOverNVtx");
    hRedMETCompPerpOverNVtx = (TH1F *) file->Get(theName+"_hRedMETCompPerpOverNVtx");
    hRedMETVsNVtx         = (TH2F *) file->Get(theName+"_hRedMETVsNVtx");
    hRedMETCompLongVsNVtx = (TH2F *) file->Get(theName+"_hRedMETCompLongVsNVtx");
    hRedMETCompPerpVsNVtx = (TH2F *) file->Get(theName+"_hRedMETCompPerpVsNVtx");
    hRedMETOverNVtxVsNVtx         = (TH2F *) file->Get(theName+"_hRedMETOverNVtxVsNVtx");	 
    hRedMETCompLongOverNVtxVsNVtx = (TH2F *) file->Get(theName+"_hRedMETCompLongOverNVtxVsNVtx");
    hRedMETCompPerpOverNVtxVsNVtx = (TH2F *) file->Get(theName+"_hRedMETCompPerpOverNVtxVsNVtx");
    hRecoilCompLong = (TH1F *) file->Get(theName+"_hRecoilCompLong");
    hRecoilCompPerp = (TH1F *) file->Get(theName+"_hRecoilCompPerp");
    hMetCompLong = (TH1F *) file->Get(theName+"_hMetCompLong");
    hMetCompPerp = (TH1F *) file->Get(theName+"_hMetCompPerp");
    hSumJetCompLong = (TH1F *) file->Get(theName+"_hSumJetCompLong");
    hSumJetCompPerp = (TH1F *) file->Get(theName+"_hSumJetCompPerp");
    hDileptonCompLong = (TH1F *) file->Get(theName+"_hDileptonCompLong");
    hDileptonCompPerp = (TH1F *) file->Get(theName+"_hDileptonCompPerp");
    hRecoilTypeLong = (TH1F *) file->Get(theName+"_hRecoilTypeLong");
    hRecoilTypePerp = (TH1F *) file->Get(theName+"_hRecoilTypePerp");
    hRecoilTypeLong->GetXaxis()->SetBinLabel(1,"Jet-Like");
    hRecoilTypeLong->GetXaxis()->SetBinLabel(2,"MET-Like");
    hRecoilTypePerp->GetXaxis()->SetBinLabel(1,"Jet-Like");
    hRecoilTypePerp->GetXaxis()->SetBinLabel(2,"MET-Like");
    hDileptSigmaPtCompLong = (TH1F *) file->Get(theName+"_hDileptSigmaPtCompLong");
    hDileptSigmaPtCompPerp = (TH1F *) file->Get(theName+"_hDileptSigmaPtCompPerp");
    hDileptSigmaPtCompLongVsMET = (TH2F *) file->Get(theName+"_hDileptSigmaPtCompLongVsMET");
    hDileptSigmaPtCompPerpVsMET = (TH2F *) file->Get(theName+"_hDileptSigmaPtCompPerpVsMET");
    hDileptonVsRecoilCompLong = (TH2F *) file->Get(theName+"_hDileptonVsRecoilCompLong");
    hDileptonVsRecoilCompPerp = (TH2F *) file->Get(theName+"_hDileptonVsRecoilCompPerp");

  }


  HistoRedMET * Clone(std::string name) {
    HistoRedMET *ret = new HistoRedMET();
    ret->theName = name;

    if(hRedMET!=0)         hRedMET->Clone((ret->theName+"_hRedMET").Data());
    if(hRedMETCompLong!=0) hRedMETCompLong->Clone((ret->theName+"_hRedMETCompLong").Data());
    if(hRedMETCompPerp!=0) hRedMETCompPerp->Clone((ret->theName+"_hRedMETCompPerp").Data());
    if(hRedMETOverNVtx!=0)         hRedMETOverNVtx->Clone((ret->theName+"_hRedMETOverNVtx").Data());	    
    if(hRedMETCompLongOverNVtx!=0) hRedMETCompLongOverNVtx->Clone((ret->theName+"_hRedMETCompLongOverNVtx").Data());
    if(hRedMETCompPerpOverNVtx!=0) hRedMETCompPerpOverNVtx->Clone((ret->theName+"_hRedMETCompPerpOverNVtx").Data());
    if(hRedMETVsNVtx!=0)         hRedMETVsNVtx->Clone((ret->theName+"_hRedMETVsNVtx").Data());	    
    if(hRedMETCompLongVsNVtx!=0) hRedMETCompLongVsNVtx->Clone((ret->theName+"_hRedMETCompLongVsNVtx").Data());
    if(hRedMETCompPerpVsNVtx!=0) hRedMETCompPerpVsNVtx->Clone((ret->theName+"_hRedMETCompPerpVsNVtx").Data());
    if(hRedMETOverNVtxVsNVtx!=0)         hRedMETOverNVtxVsNVtx->Clone((ret->theName+"_hRedMETOverNVtxVsNVtx").Data());	
    if(hRedMETCompLongOverNVtxVsNVtx!=0) hRedMETCompLongOverNVtxVsNVtx->Clone((ret->theName+"_hRedMETCompLongOverNVtxVsNVtx").Data());
    if(hRedMETCompPerpOverNVtxVsNVtx!=0) hRedMETCompPerpOverNVtxVsNVtx->Clone((ret->theName+"_hRedMETCompPerpOverNVtxVsNVtx").Data());
    if(hRecoilCompLong !=0) hRecoilCompLong->Clone((ret->theName+"_hRecoilCompLong").Data());
    if(hRecoilCompPerp !=0) hRecoilCompPerp->Clone((ret->theName+"_hRecoilCompPerp").Data());
    if(hMetCompLong !=0) hMetCompLong->Clone((ret->theName+"_hMetCompLong").Data());
    if(hMetCompPerp !=0) hMetCompPerp->Clone((ret->theName+"_hMetCompPerp").Data());
    if(hSumJetCompLong !=0) hSumJetCompLong->Clone((ret->theName+"_hSumJetCompLong").Data());
    if(hSumJetCompPerp !=0) hSumJetCompPerp->Clone((ret->theName+"_hSumJetCompPerp").Data());
    if(hDileptonCompLong !=0) hDileptonCompLong->Clone((ret->theName+"_hDileptonCompLong").Data());
    if(hDileptonCompPerp !=0) hDileptonCompPerp->Clone((ret->theName+"_hDileptonCompPerp").Data());
    if(hRecoilTypeLong !=0) hRecoilTypeLong->Clone((ret->theName+"_hRecoilTypeLong").Data());
    if(hRecoilTypePerp !=0) hRecoilTypePerp->Clone((ret->theName+"_hRecoilTypePerp").Data());

    if(hDileptSigmaPtCompLong != 0) hDileptSigmaPtCompLong->Clone((ret->theName+"_hDileptSigmaPtCompLong").Data());
    if(hDileptSigmaPtCompPerp != 0) hDileptSigmaPtCompPerp->Clone((ret->theName+"_hDileptSigmaPtCompPerp").Data());
    if(hDileptSigmaPtCompLongVsMET != 0) hDileptSigmaPtCompLongVsMET->Clone((ret->theName+"_hDileptSigmaPtCompLongVsMET").Data());
    if(hDileptSigmaPtCompPerpVsMET != 0) hDileptSigmaPtCompPerpVsMET->Clone((ret->theName+"_hDileptSigmaPtCompPerpVsMET").Data());
    if(hDileptonVsRecoilCompLong != 0) hDileptonVsRecoilCompLong->Clone((ret->theName+"_hDileptonVsRecoilCompLong").Data());
    if(hDileptonVsRecoilCompPerp != 0) hDileptonVsRecoilCompPerp->Clone((ret->theName+"_hDileptonVsRecoilCompPerp").Data());
    
    return ret;
  }



  void Add(const HistoRedMET* histSet) {
    if(hRedMET!=0)         hRedMET->Add(histSet->hRedMET);
    if(hRedMETCompLong!=0) hRedMETCompLong->Add(histSet->hRedMETCompLong);
    if(hRedMETCompPerp!=0) hRedMETCompPerp->Add(histSet->hRedMETCompPerp);
    if(hRedMETOverNVtx!=0)         hRedMETOverNVtx->Add(histSet->hRedMETOverNVtx);	    
    if(hRedMETCompLongOverNVtx!=0) hRedMETCompLongOverNVtx->Add(histSet->hRedMETCompLongOverNVtx);
    if(hRedMETCompPerpOverNVtx!=0) hRedMETCompPerpOverNVtx->Add(histSet->hRedMETCompPerpOverNVtx);
    if(hRedMETVsNVtx!=0)         hRedMETVsNVtx->Add(histSet->hRedMETVsNVtx);	    
    if(hRedMETCompLongVsNVtx!=0) hRedMETCompLongVsNVtx->Add(histSet->hRedMETCompLongVsNVtx);
    if(hRedMETCompPerpVsNVtx!=0) hRedMETCompPerpVsNVtx->Add(histSet->hRedMETCompPerpVsNVtx);
    if(hRedMETOverNVtxVsNVtx!=0)         hRedMETOverNVtxVsNVtx->Add(histSet->hRedMETOverNVtxVsNVtx);	
    if(hRedMETCompLongOverNVtxVsNVtx!=0) hRedMETCompLongOverNVtxVsNVtx->Add(histSet->hRedMETCompLongOverNVtxVsNVtx);
    if(hRedMETCompPerpOverNVtxVsNVtx!=0) hRedMETCompPerpOverNVtxVsNVtx->Add(histSet->hRedMETCompPerpOverNVtxVsNVtx);
    if(hRecoilCompLong != 0) hRecoilCompLong->Add(histSet->hRecoilCompLong);
    if(hRecoilCompPerp != 0) hRecoilCompPerp->Add(histSet->hRecoilCompPerp);
    if(hMetCompLong != 0) hMetCompLong->Add(histSet->hMetCompLong);
    if(hMetCompPerp != 0) hMetCompPerp->Add(histSet->hMetCompPerp);
    if(hSumJetCompLong != 0) hSumJetCompLong->Add(histSet->hSumJetCompLong);
    if(hSumJetCompPerp != 0) hSumJetCompPerp->Add(histSet->hSumJetCompPerp);
    if(hDileptonCompLong != 0) hDileptonCompLong->Add(histSet->hDileptonCompLong);
    if(hDileptonCompPerp != 0) hDileptonCompPerp->Add(histSet->hDileptonCompPerp);
    if(hRecoilTypeLong != 0) hRecoilTypeLong->Add(histSet->hRecoilTypeLong);
    if(hRecoilTypePerp != 0) hRecoilTypePerp->Add(histSet->hRecoilTypePerp);
    if(hDileptSigmaPtCompLong != 0) hDileptSigmaPtCompLong->Add(histSet->hDileptSigmaPtCompLong);
    if(hDileptSigmaPtCompPerp != 0) hDileptSigmaPtCompPerp->Add(histSet->hDileptSigmaPtCompPerp);
    if(hDileptSigmaPtCompLongVsMET != 0) hDileptSigmaPtCompLongVsMET->Add(histSet->hDileptSigmaPtCompLongVsMET);
    if(hDileptSigmaPtCompPerpVsMET != 0) hDileptSigmaPtCompPerpVsMET->Add(histSet->hDileptSigmaPtCompPerpVsMET);
    if(hDileptonVsRecoilCompLong != 0) hDileptonVsRecoilCompLong->Add(histSet->hDileptonVsRecoilCompLong);
    if(hDileptonVsRecoilCompPerp != 0) hDileptonVsRecoilCompPerp->Add(histSet->hDileptonVsRecoilCompPerp);

  }



  void Scale(double scaleFact) {
    if(hRedMET!=0)         hRedMET->Scale(scaleFact);
    if(hRedMETCompLong!=0) hRedMETCompLong->Scale(scaleFact);
    if(hRedMETCompPerp!=0) hRedMETCompPerp->Scale(scaleFact);
    if(hRedMETOverNVtx!=0)         hRedMETOverNVtx->Scale(scaleFact);
    if(hRedMETCompLongOverNVtx!=0) hRedMETCompLongOverNVtx->Scale(scaleFact);
    if(hRedMETCompPerpOverNVtx!=0) hRedMETCompPerpOverNVtx->Scale(scaleFact);
    if(hRedMETVsNVtx!=0)         hRedMETVsNVtx->Scale(scaleFact);
    if(hRedMETCompLongVsNVtx!=0) hRedMETCompLongVsNVtx->Scale(scaleFact);
    if(hRedMETCompPerpVsNVtx!=0) hRedMETCompPerpVsNVtx->Scale(scaleFact);
    if(hRedMETOverNVtxVsNVtx!=0)         hRedMETOverNVtxVsNVtx->Scale(scaleFact);
    if(hRedMETCompLongOverNVtxVsNVtx!=0) hRedMETCompLongOverNVtxVsNVtx->Scale(scaleFact);
    if(hRedMETCompPerpOverNVtxVsNVtx!=0) hRedMETCompPerpOverNVtxVsNVtx->Scale(scaleFact);
    if(hRecoilCompLong != 0) hRecoilCompLong->Scale(scaleFact);
    if(hRecoilCompPerp != 0) hRecoilCompPerp->Scale(scaleFact);
    if(hMetCompLong != 0) hMetCompLong->Scale(scaleFact);
    if(hMetCompPerp != 0) hMetCompPerp->Scale(scaleFact);
    if(hSumJetCompLong != 0) hSumJetCompLong->Scale(scaleFact);
    if(hSumJetCompPerp != 0) hSumJetCompPerp->Scale(scaleFact);
    if(hDileptonCompLong != 0) hDileptonCompLong->Scale(scaleFact);
    if(hDileptonCompPerp != 0) hDileptonCompPerp->Scale(scaleFact);
    if(hRecoilTypeLong != 0) hRecoilTypeLong->Scale(scaleFact);
    if(hRecoilTypePerp != 0) hRecoilTypePerp->Scale(scaleFact);
    if(hDileptSigmaPtCompLong != 0) hDileptSigmaPtCompLong->Scale(scaleFact);
    if(hDileptSigmaPtCompPerp != 0) hDileptSigmaPtCompPerp->Scale(scaleFact);
    if(hDileptSigmaPtCompLongVsMET != 0) hDileptSigmaPtCompLongVsMET->Scale(scaleFact);
    if(hDileptSigmaPtCompPerpVsMET != 0) hDileptSigmaPtCompPerpVsMET->Scale(scaleFact);
    if(hDileptonVsRecoilCompLong != 0) hDileptonVsRecoilCompLong->Scale(scaleFact);
    if(hDileptonVsRecoilCompPerp != 0) hDileptonVsRecoilCompPerp->Scale(scaleFact);

  }


  void Write() {
    if(hRedMET!=0)         hRedMET->Write();
    if(hRedMETCompLong!=0) hRedMETCompLong->Write();
    if(hRedMETCompPerp!=0) hRedMETCompPerp->Write();
    if(hRedMETOverNVtx!=0)         hRedMETOverNVtx->Write();
    if(hRedMETCompLongOverNVtx!=0) hRedMETCompLongOverNVtx->Write();
    if(hRedMETCompPerpOverNVtx!=0) hRedMETCompPerpOverNVtx->Write();
    if(hRedMETVsNVtx!=0)         hRedMETVsNVtx->Write();
    if(hRedMETCompLongVsNVtx!=0) hRedMETCompLongVsNVtx->Write();
    if(hRedMETCompPerpVsNVtx!=0) hRedMETCompPerpVsNVtx->Write();
    if(hRedMETOverNVtxVsNVtx!=0)         hRedMETOverNVtxVsNVtx->Write();
    if(hRedMETCompLongOverNVtxVsNVtx!=0) hRedMETCompLongOverNVtxVsNVtx->Write();
    if(hRedMETCompPerpOverNVtxVsNVtx!=0) hRedMETCompPerpOverNVtxVsNVtx->Write();
    if(hRecoilCompLong != 0) hRecoilCompLong->Write();
    if(hRecoilCompPerp != 0) hRecoilCompPerp->Write();
    if(hMetCompLong != 0) hMetCompLong->Write();
    if(hMetCompPerp != 0) hMetCompPerp->Write();
    if(hSumJetCompLong != 0) hSumJetCompLong->Write();
    if(hSumJetCompPerp != 0) hSumJetCompPerp->Write();
    if(hDileptonCompLong != 0) hDileptonCompLong->Write();
    if(hDileptonCompPerp != 0) hDileptonCompPerp->Write();
    if(hRecoilTypeLong != 0) hRecoilTypeLong->Write();
    if(hRecoilTypePerp != 0) hRecoilTypePerp->Write();
    if(hDileptSigmaPtCompLong != 0) hDileptSigmaPtCompLong->Write();
    if(hDileptSigmaPtCompPerp != 0) hDileptSigmaPtCompPerp->Write();
    if(hDileptSigmaPtCompLongVsMET != 0) hDileptSigmaPtCompLongVsMET->Write();
    if(hDileptSigmaPtCompPerpVsMET != 0) hDileptSigmaPtCompPerpVsMET->Write();
    if(hDileptonVsRecoilCompLong != 0) hDileptonVsRecoilCompLong->Write();
    if(hDileptonVsRecoilCompPerp != 0) hDileptonVsRecoilCompPerp->Write();

  }

#ifndef ROOTANALYSIS
  void Fill(const ReducedMETComputer* redMedComputer, double met, unsigned int nVtx, double weight) {
    Fill(redMedComputer->reducedMET(),
	 redMedComputer->reducedMETComponents().first, redMedComputer->reducedMETComponents().second, 
	 redMedComputer->recoilProjComponents().first, redMedComputer->recoilProjComponents().second,
	 redMedComputer->metProjComponents().first, redMedComputer->metProjComponents().second,
	 redMedComputer->sumJetProjComponents().first, redMedComputer->sumJetProjComponents().second,
	 redMedComputer->dileptonProjComponents().first, redMedComputer->dileptonProjComponents().second,
	 redMedComputer->dileptonPtCorrComponents().first, redMedComputer->dileptonPtCorrComponents().second,
	 redMedComputer->recoilType().first, redMedComputer->recoilType().second,
	 met, 
	 nVtx, 
	 weight);
  }
  
#endif

  void Fill(double redmet,
	    double redmet_long, double redmet_perp,
	    double recoil_long, double recoil_perp,
	    double met_long, double met_perp,
	    double sumjet_long, double sumjet_perp,
	    double dilepton_long, double dilepton_perp,
	    double dileptonSigma_long, double dileptonSigma_perp,
	    double type_long, double type_perp,
	    double met,
	    unsigned int nVtx,
	    double weight) {
    hRedMET->Fill(redmet, weight);
    hRedMETCompLong->Fill(redmet_long, weight);
    hRedMETCompPerp->Fill(redmet_perp, weight);

    if(nVtx!=0) hRedMETOverNVtx->Fill(redmet/nVtx, weight);
    if(nVtx!=0) hRedMETCompLongOverNVtx->Fill(redmet_long/nVtx, weight);
    if(nVtx!=0) hRedMETCompPerpOverNVtx->Fill(redmet_perp/nVtx, weight);
    hRedMETVsNVtx->Fill(nVtx, redmet, weight);
    hRedMETCompLongVsNVtx->Fill(nVtx, redmet_long, weight);
    hRedMETCompPerpVsNVtx->Fill(nVtx, redmet_perp, weight);
    if(nVtx!=0) hRedMETOverNVtxVsNVtx->Fill(nVtx, redmet/nVtx, weight);
    if(nVtx!=0) hRedMETCompLongOverNVtxVsNVtx->Fill(nVtx, redmet_long/nVtx, weight);
    if(nVtx!=0) hRedMETCompPerpOverNVtxVsNVtx->Fill(nVtx, redmet_perp/nVtx, weight);
    hRecoilCompLong->Fill(recoil_long, weight);
    hRecoilCompPerp->Fill(recoil_perp, weight);
    hMetCompLong->Fill(met_long, weight);
    hMetCompPerp->Fill(met_perp, weight);
    hSumJetCompLong->Fill(sumjet_long, weight);
    hSumJetCompPerp->Fill(sumjet_perp, weight);
    hDileptonCompLong->Fill(dilepton_long, weight);
    hDileptonCompPerp->Fill(dilepton_perp, weight);
    hRecoilTypeLong->Fill(type_long, weight);
    hRecoilTypePerp->Fill(type_perp, weight);
    hDileptSigmaPtCompLong->Fill(dileptonSigma_long, weight);
    hDileptSigmaPtCompPerp->Fill(dileptonSigma_perp, weight);
    hDileptSigmaPtCompLongVsMET->Fill(met, dileptonSigma_long, weight);
    hDileptSigmaPtCompPerpVsMET->Fill(met, dileptonSigma_perp, weight);
    hDileptonVsRecoilCompLong->Fill(recoil_long, dilepton_long, weight);
    hDileptonVsRecoilCompPerp->Fill(recoil_perp, dilepton_perp, weight);

  }

  /// Destructor
  virtual ~HistoRedMET() {}

  // Operations
  TString theName;

  TH1F *hRedMET;
  TH1F *hRedMETCompLong;
  TH1F *hRedMETCompPerp;
  TH1F *hRedMETOverNVtx;
  TH1F *hRedMETCompLongOverNVtx;
  TH1F *hRedMETCompPerpOverNVtx;
  TH2F *hRedMETVsNVtx;
  TH2F *hRedMETCompLongVsNVtx;
  TH2F *hRedMETCompPerpVsNVtx;
  TH2F *hRedMETOverNVtxVsNVtx;
  TH2F *hRedMETCompLongOverNVtxVsNVtx;
  TH2F *hRedMETCompPerpOverNVtxVsNVtx;
  TH1F *hRecoilCompLong;
  TH1F *hRecoilCompPerp;
  TH1F *hMetCompLong;
  TH1F *hMetCompPerp;
  TH1F *hSumJetCompLong;
  TH1F *hSumJetCompPerp;
  TH1F *hDileptonCompLong;
  TH1F *hDileptonCompPerp;
  TH1F *hDileptSigmaPtCompLong;
  TH1F *hDileptSigmaPtCompPerp;
  TH2F *hDileptSigmaPtCompLongVsMET;
  TH2F *hDileptSigmaPtCompPerpVsMET;
  TH2F *hDileptonVsRecoilCompLong;
  TH2F *hDileptonVsRecoilCompPerp;
  TH1F *hRecoilTypeLong;
  TH1F *hRecoilTypePerp;

};



#endif
