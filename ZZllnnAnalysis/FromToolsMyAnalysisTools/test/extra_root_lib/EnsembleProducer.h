#ifndef EnsembleProducer_H
#define EnsembleProducer_H

/** \class EnsembleProducer
 *  This class produces ensemble trees starting from a bigger MC sample.
 *  The event weight is correctly taken into account in the ensemble generation,
 *  this can lead to possible inefficiencies.
 *
 *  $Date: 2007/08/10 19:15:15 $
 *  $Revision: 1.1 $
 *  \author G. Cerminara - NEU Boston & INFN Torino
 */



class TTree;
class TFile;
class TNtuple;
class RndmGenerator;

#include "TString.h"
#include <map>


class EnsembleProducer {
public:
  /// Constructor
  EnsembleProducer(const TString& name, TTree *tree, const double expectedEvents,
		   double massMin = 0, double massMax = 99999);

  /// Destructor
  virtual ~EnsembleProducer();

  // Operations
  int getNEnsembles() const;

  TNtuple* getEnsembleN(int i);

  void setSeed(unsigned int seed);

  int getExpectedAverage() const;

protected:

private:

  // Estimates the yeld of the sample on the basis of the weight
  double getYeld();
  // rounding the expected number of events to get the mean of the poisson distribution
  int rounding(double expectedEvents) const;
  // Get the ensemble i: if it is not in the map it will build it
  TNtuple *getEnsemble(int i);

  bool massOutWindow();

  // this does the job
  void produce();
  void initialize();
  void write();

  TTree *theTree;

  double massMinCut;
  double massMaxCut;

  // The tree variables 
  float run; 
  float event; 
  float type;
  float DiLeptMinv;
  float weight;
  float LeadLeptPt;
  float NoLeadLeptPt;
  float DeltaPhiLeadLeptZ;
  float CosThetaStarLeptMinus;
  float DiLeptDeltaEta;
  float TransMass4Body;
  float MinMass4Body;



  // the number of ensembles which will be produced
  int nEns;
  int entriesInMassWind;
  double maxWeight;

  int nExpectedAverage;

  TFile *theFile;
  std::map<int, TNtuple *> theEnsambles;
  
  RndmGenerator *rndm;

};
#endif

