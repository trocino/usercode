#ifndef TreeCombiner_H
#define TreeCombiner_H

/** \class TreeCombiner
 *  Combines the various TTree for differnt samples into one tree weighting the
 *  events correctly. The event weight is computed using the LumiNormalization class.
 *  The trees DiLeptNtuple (and optionally KinemNtuple) are read from the files DiLeptDistr_*.root
 *
 *  $Date: 2008/03/13 17:44:26 $
 *  $Revision: 1.3 $
 *  \author G. Cerminara - NEU Boston & INFN Torino
 */

class TNtuple;
class LumiNorm;
class TFile;

#include "TString.h"
#include <vector>
#include <map>

class TreeCombiner {
public:
  /// Constructor
  // Input parameters are:
  //    - lumiNorm: pointer to the class handling the normalization
  //    - finalState: "diem" or "dimu" (used to choose the correct file)
  //    - inputDir: directory containing the root files
  TreeCombiner(const LumiNorm *lumiNorm, const TString &inputDir, const TString &fileName);

  /// Destructor
  virtual ~TreeCombiner();
  
  // Get the tree for the data sample
  TNtuple *getDataTree() const;

  // Get the tree for a combination of MC samples (specified with a vector of string)
  TNtuple *getCombinedTree(const std::vector<TString>& sampleNames);

  // Reject events using a mass window (default is no cut)
  void setMassCut(double min, double max);
  void requireFullSection(bool doit);

protected:

private:
  // Pointer to the normalization tool (not owned by this class)
  const LumiNorm *theLumiNorm;
  const TString theFileName;
  const TString theInputDir;
  double massMinCut;
  double massMaxCut;
  bool fullSelectionReq;

  // The tree variables 
  std::map<TString, float> treeVars;
  std::vector<TString> varNames;
  //   float run; 
  //   float lumi;
  //   float event; 
  //   float weight;
  //   float leadCharge;
  //   float subleadCharge;
  //   float leadPt;
  //   float subleadPt;
  //   float diLeptInvMass;
  //   float dileptLeadDeltaPhi;
  //   float leptMinusCmCosTheta;
  //   float passAnalyzed;
  //   float passDilepton;
  //   float passMassWindow; 
  //   float passRedMETcut;
  //   float passJetVeto;
  //   float passLeptonVeto;
  //   float passAllCuts;

};
#endif

