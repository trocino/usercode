#ifndef ZZllvvAnalyzer_H
#define ZZllvvAnalyzer_H

/** \class ZZllvvAnalyzer
 *  No description available.
 *
 *  $Date: 2011/07/11 19:28:23 $
 *  $Revision: 1.4 $
 *  \author G. Cerminara - CERN
 *          D. Trocino   - Northeastern University
 */
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"

#include "CMGTools/HtoZZ2l2nu/interface/ReducedMETComputer.h"

class TFile;
class TNtuple;
class TString;

namespace reco {
  class Vertex;
  class Muon;
}

namespace pat {
  //class Muon;
  class MET;
}

class ZZllvvAnalyzer   : public edm::EDAnalyzer {
public:
  /// Constructor
  ZZllvvAnalyzer(const edm::ParameterSet& pSet);

  /// Destructor
  virtual ~ZZllvvAnalyzer();

protected:

  // Operations
  virtual void beginJob();

  virtual void beginRun(const edm::Run& run, const edm::EventSetup& eSetup);

  virtual void beginLuminosityBlock(const edm::LuminosityBlock & iLumi, const edm::EventSetup & iSetup);
  virtual void endLuminosityBlock(const edm::LuminosityBlock & iLumi, const edm::EventSetup & iSetup);

  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();



private:

  int totNEvents;
  bool isMC;
  float theXsect;
  float theBrRatio;
  float theLumi;

  TFile *theFile;
  TNtuple *finalNtpl;
  std::string theOutFileName;
  edm::InputTag source;
  bool debug;
  edm::ParameterSet vertexSelection;

  int flavorCombo; // 0: EE/MM (default);  1=EE;  2=MM;  3=EM;  -1=ANY
  int chargeCombo; // -1: OPPOSITE (default);  0: ANY;  1: SAME

  // Parameters for RedMET
  double kRecoilLongWeight;
  double kRecoilPerpWeight;
  double kSigmaPtLongWeight;
  double kSigmaPtPerpWeight;
  double kPerpComponentWeight;
  double theRedMETMinCut;
  bool useAllJets;

  void initializePlots();
  void fillPlots(std::string, 
		 const reco::Vertex *, 
		 unsigned int, 
		 unsigned int, 
		 int, 
		 int, 
		 reco::CandidatePtr, 
		 reco::CandidatePtr, 
		 reco::CandidatePtr, 
		 double, double, double, 
		 //edm::RefVector<std::vector<reco::Muon> > &, 
		 //edm::RefVector<std::vector<reco::Muon> > &, 
		 const pat::MET *, 
		 ReducedMETComputer *, 
		 std::vector<LorentzVector> &//,
		 //int
		 );
  void scalePlots(double);
  void writePlots();
  std::vector<double> makePuDistr(const edm::ParameterSet&);
  void resetFlags();

};
#endif

