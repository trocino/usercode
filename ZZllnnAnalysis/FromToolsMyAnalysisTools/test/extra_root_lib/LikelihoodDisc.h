#ifndef LikelihoodDisc_H
#define LikelihoodDisc_H

/** \class LikelihoodDisc
 *  Implementation of a likelihood discriminant.
 *  The class accept new variables from a tree, builds the correspondig pdfs 
 *  and stores them. The pdf can be retrieved by name and combined by the 
 *  evaluate method. 
 *  The output of the likelihood is computed as pdf(signal)/[pdf(signal)+pdf(background)].
 *
 *  $Date: 2008/03/13 17:41:12 $
 *  $Revision: 1.3 $
 *  \author G. Cerminara - NEU Boston & INFN Torino
 */

class TTree;
class TString;
class TH1D;
class TH1F;
class TF1;
class TLorentzVector;

#include <map>
#include <vector>

class LikelihoodDisc {
public:
  /// Constructor
  // Accept as arguments the trees which store the variables for the training
  LikelihoodDisc(TTree *sigTree, TTree *bckTree);

  /// Destructor
  virtual ~LikelihoodDisc();

  // Operations
  // Add a variable to the likelihood. The variable is read from the tree (branch name = varName)
  // and a PDF histo is filled (# bins = nBins and axis limits = min and max) for signal and background.
  // The ratio of signal and backgorund PDF is computed (same binning of the PDF histo) and fitted with the 
  // function specified by fitFuncLR (see ROOT convetion for function names).
  // Both the PDF histos and the likelihood ratio are stored with the sub-likelihood name lName
  void addVariable(const TString& lName, const TString& varName, const TString& title,
		   const int nBins, double min, double max, const TString& fitFuncLR);

  // Add a variable to the likelihood. The variable is read from the tree (branch name = varName)
  // and a PDF histo is filled (# bins = nBins and axis limits = min and max) for signal and background.
  // The ratio of signal and backgorund PDF is computed (same binning of the PDF histo) and fitted with the 
  // function specified by fitFuncLR which must be instantiated by the user.
  // Both the PDF histos and the likelihood ratio are stored with the sub-likelihood name lName
  void addVariable(const TString& lName, const TString& varName, const TString& title,
		   const int nBins, double min, double max, TF1* fitFunc);

  // Add a variable to the likelihood. The variable is read from the tree (branch name = varName)
  // and a PDF histo is filled (# bins = nBins, axis limits = min and max, variable binning = bins)
  // for signal and background.
  // The ratio of signal and backgorund PDF is computed (same binning of the PDF histo) and fitted with the 
  // function specified by fitFuncLR which must be instantiated by the user.
  // Both the PDF histos and the likelihood ratio are stored with the sub-likelihood name lName
  void addVariable(const TString& lName, const TString& varName, const TString& title,
		   const int nBins, double min, double max, const double* bins, TF1* fitFuncLR);

  // Add a variable to the likelihood. The variable is read from the tree (branch name = varName)
  // and a PDF histo is filled (# bins = nBins, axis limits = min and max, variable binning = bins)
  // for signal and background.
  // The ratio of signal and backgorund PDF is computed (same binning of the PDF histo) and fitted with the 
  // function specified by fitFuncLR (see ROOT convetion for function names).
  // Both the PDF histos and the likelihood ratio are stored with the sub-likelihood name lName
  void addVariable(const TString& lName, const TString& varName, const TString& title,
		   const int nBins, double min, double max, const double* bins, 
		   const TString& fitFuncLR);

  /*
  // Add the chi2 of the kinematic fit as a variable of the likelihood.
  // The variable name is Chi2 in this case. The flag useProbability allows to use the
  // chi2 probability (chi2 for 1dof) instead of the chi2 itself 
  void addChi2Variable(const TString& lName, const TString& title,
		       const int nBins, double min, double max, 
		       const TString& fitFuncLR, bool useProbability = true);

  // Add the chi2 of the kinematic fit as a variable of the likelihood.
  // The variable name is Chi2 in this case. The flag useProbability allows to use the
  // chi2 probability (chi2 for 1dof) instead of the chi2 itself 
  void addChi2Variable(const TString& lName, const TString& title,
		       const int nBins, double min, double max, const double* bins, 
		       const TString& fitFuncLR, bool useProbability = true);

  // Given the 4-momenta of the 2 leptons and the errors on 1/pt minimize the chi2
  // for the kinematic fit built using the mass constraint
  double computeChi2KinemFit(const TLorentzVector& lept1, float dOneOverPt1,
			     const TLorentzVector& lept2, float dOneOverPt2) const;

  double computeChi2KinemFitV2(const TLorentzVector& lept1, float dOneOverPt1,
			       const TLorentzVector& lept2, float dOneOverPt2) const;
  */



  // Retrieve the signal and backgorund PDF for a particular sub-likelihood
  std::pair<TH1D *, TH1D *> getPdf(const TString& lName) const;
  
  // Retrieve the Likelihood ratio (histo and fitted function) for a particular sub-likelihood
  std::pair<TH1D *, TF1 *> getLR(const TString& lName) const;

  // Evaluate the Likelihood Ratio for a particular likelihood composition.
  // The sublikelihhod to be used are passed by a vector of names: likelihoodComposition
  // The value of each variable corresponding to the sub-likelihood used is passed in the vector: variables
  // Thi must follow the same order specified by likelihoodComposition
  double evaluate(const std::vector<TString>& likelihoodComposition,
		  const std::vector<float>& variables) const;

  // Loop on the given tree and fill an histo with the likelihood output
  // The sublikelihhod to be used are passed by a vector of names: likelihoodComposition
  TH1F * evaluate(const std::vector<TString>& likelihoodComposition,
		  TTree *tree, const TString& sampleName, const int nBins = 10) const;


  // Loop on the given tree and fill an histo with the likelihood output
  // The sublikelihhod to be used are passed by a vector of names: likelihoodComposition
  // Uses variable binning for the output histo
  TH1F * evaluate(const std::vector<TString>& likelihoodComposition,
		  TTree *tree, const TString& sampleName,
		  const int nBins, const double *bins) const;


  // Loop on the given DATA tree and fill an histo with the likelihood output
  // The sublikelihhod to be used are passed by a vector of names: likelihoodComposition
  TH1F * evaluateData(const std::vector<TString>& likelihoodComposition,
		      TTree *tree, const TString& sampleName, double llCut,
		      const int nBins = 10) const;
  
  // Set the mass cut to be used building the PDF and evaluating the likelihood output
  // default: no cut
  void setMassCut(double min, double max);
  
  // Specify to separate the samples for likelihood building and likelihood evaluation
  // default: doSeparate = true
  void separateSamples(bool doSeparate);

  /*
  // If set to true both leptons are required to have SMT hits:
  // SMTswitch = 0 -> no requirements (default)
  // SMTswitch = 1 -> at most 1 track without SMT
  // SMTswitch = 2 -> both tracks with SMT hits
  // SMTswitch = 3 -> both tracks without SMT hits
  void selectTracksWithSMT(unsigned int SMTswitch);

  // Switch on the cut on the alatCorrReduceT branch and set the value of the cut
  void setAlatCut(double minAlat);

  // Switch on the cut on the alatCorrReduceT branch and reset the value of the cut to -9999
  void noAlatCut();
  */
  
  // set the verbosity level (particularly important for MINUIT verbosity
  // Valid levels are:
  //    -1 quiet (also suppresse all warnings)
  //    0  normal
  //    1  verbose
  void setVerbosityLevel(int level);

  // Switch on the cut on the lepton charge
  // 0 -> no cut
  // -1 -> opposite charge
  // 1 -> same charge
  void setChargeCut(int cut);

  void useMassCut(bool doit);

  /*
  // Evaluate the chi2 probability of the kinematic fit (useful to plot the likelihood input)
  TH1F * plotChi2Prob(TTree *tree, const TString& sampleName, const int nBins = 10) const;
  */

protected:

private:

  void normalize(TH1D *histo);
  std::pair<TH1D *, TF1 *> buildLikelihoodRatio(const TH1D* hSig, const TH1D* hBck,
						TString fitName, TF1* LRatioFit,
						double funcMin = 0, double funcMax = 1);
  double transformLikelihoodOutput(double r) const;


  // The trees used for the training
  TTree *theSigTree;
  TTree *theBckTree;

  std::map<TString, std::pair<TH1D *, TH1D *> > thePdfMap;
  std::map<TString, std::pair<TH1D *, TF1 *> > theLRMap;
  std::map<TString, TString> theVariableMap;

  double massMin;
  double massMax;
  // float mass;
  
  bool doSeparateSample;

  /*  
  unsigned int nSMTHitsSwitch ;
  // Switch to cut on alat_corr_reduced_t (the relative branches are read only if necessary)
  bool doAlatCut;
  // Value of the cut on alat_corr_reduced_t (if used)
  double alatMin;
  // Flag to use the chi2 prob distribution instead of the chi2 distribution
  bool useChi2Probability;
  */

  // switch for the verbosity
  bool debug;
  // set the verbosity level
  int verbosityLevel;
  // the switch for the charge cut
  int theChargeCut;
  bool doMassCut;


};

/*
// The function to be minimized
void chi2function(int &npar, double *grad, double &fval, double *par, int flag);

// The function to be minimized
void chi2functionV2(int &npar, double *grad, double &fval, double *par, int flag);
*/

#endif

