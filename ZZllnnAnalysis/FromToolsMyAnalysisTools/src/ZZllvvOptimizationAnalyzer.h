#ifndef ZZllvvOptimizationAnalyzer_H
#define ZZllvvOptimizationAnalyzer_H

/** \class ZZllvvOptimizationAnalyzer
 *  No description available.
 *
 *  $Date: 2011/04/07 15:21:48 $
 *  $Revision: 1.4 $
 *  \author G. Cerminara - CERN
 */
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"


class TFile;
class HistoLept;

class ZZllvvOptimizationAnalyzer   : public edm::EDAnalyzer {
public:
  /// Constructor
  ZZllvvOptimizationAnalyzer(const edm::ParameterSet& pSet);

  /// Destructor
  virtual ~ZZllvvOptimizationAnalyzer();

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
  float weight;

  TFile *theFile;
  edm::InputTag source;
  edm::InputTag zmmInput;
  bool debug;
  edm::ParameterSet vertexSelection;
  std::vector<double> kRecoilSteps;
  std::vector<double> kSigmaPtSteps;
  std::vector<double> kPerpSteps;

};
#endif

