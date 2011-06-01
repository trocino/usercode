#ifndef ZZllvvAnalyzer_H
#define ZZllvvAnalyzer_H

/** \class ZZllvvAnalyzer
 *  No description available.
 *
 *  $Date: 2011/04/07 15:21:48 $
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

#include "CMGTools/HtoZZ2l2nu/interface/ReducedMETComputer.h"

class TFile;
class TString;

namespace reco {
  class Vertex;
  class Muon;
}

namespace pat {
  class Muon;
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
  edm::InputTag source;
  edm::InputTag zmmInput;
  bool debug;
  edm::ParameterSet vertexSelection;

  void initializePlots();
  void fillPlots(std::string, 
		 const reco::Vertex *, 
		 unsigned int, 
		 unsigned int, 
		 int, 
		 const pat::Muon *, 
		 const pat::Muon *, 
		 edm::RefVector<std::vector<reco::Muon> > &, 
		 edm::RefVector<std::vector<reco::Muon> > &, 
		 const pat::MET *, 
		 ReducedMETComputer *, 
		 std::vector<LorentzVector> &,
		 std::vector<LorentzVector> &);
  void scalePlots(double);
  void writePlots();
  std::vector<double> makePuDistr(const edm::ParameterSet&);

};
#endif

