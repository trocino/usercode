
/*
 *  See header file for a description of this class.
 *
 *  $Date: 2011/05/17 14:02:25 $
 *  $Revision: 1.1 $
 *  \author G. Cerminara - CERN
 */

#include "CMGTools/HtoZZ2l2nu/interface/Utils.h"
//#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

Utils::Utils(){}

Utils::~Utils(){}

double Utils::computeRelIsolation(const pat::Muon *muon) {
  //isolation
  double norm = std::max((double)20.0,(double)muon->pt());
  double ecalIso = muon->ecalIso();
  double hcalIso = muon->hcalIso();
  double trkIso = muon->trackIso();
  double relIso = (ecalIso+hcalIso+trkIso)/norm;
      
  return relIso;
}

double Utils::computeRelIsolation(const pat::Electron *ele) {
  //isolation
  double norm = std::max((double)20.0,(double)ele->pt());
  double ecalIso = ele->ecalIso();
  double hcalIso = ele->hcalIso();
  double trkIso = ele->trackIso();
  double relIso = (ecalIso+hcalIso+trkIso)/norm;
      
  return relIso;
}

/*
double Utils::computeRelIsolation(const reco::CandidatePtr lept) {
  //isolation
  double norm = std::max((double)20.0,(double)lept->pt());
  double ecalIso = lept->ecalIso();
  double hcalIso = lept->hcalIso();
  double trkIso = lept->trackIso();
  double relIso = (ecalIso+hcalIso+trkIso)/norm;
      
  return relIso;
}
*/

int Utils::muonType(const pat::Muon *muon) {
  int type = -1;
  if(muon->isTrackerMuon() && !muon->isGlobalMuon()) {
    type = 0;
  } else if(!muon->isTrackerMuon() && muon->isGlobalMuon()) {
    type = 1;
  } else if(muon->isTrackerMuon() && muon->isGlobalMuon()) {
    type = 2;
  }
  return type;
}

int Utils::muonType(const pat::Electron *ele) {
  return -1;
}

/*
int  Utils::muonType(const reco::CandidatePtr lept) {
  int type = -1;
  const pat::Muon *muon=dynamic_cast<const pat::Muon*>(lept.get());
  if(muon!=0) {
    if(muon->isTrackerMuon() && !muon->isGlobalMuon()) {
      type = 0;
      //muonTrack = muon->globalTrack();
    } else if(!muon->isTrackerMuon() && muon->isGlobalMuon()) {
      type = 1;
    } else if(muon->isTrackerMuon() && muon->isGlobalMuon()) {
      type = 2;
    }
  }
  return type;
}
*/

