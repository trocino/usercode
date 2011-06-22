#ifndef RndmGenerator_H
#define RndmGenerator_H


/** \class RndmGenerator
 *  Interface to the Root Mersenne Twistor random generator.
 *  
 *
 *  $Date: 2007/08/10 19:15:15 $
 *  $Revision: 1.1 $
 *  \author G. Cerminara - NEU Boston & INFN Torino
 */

#if !defined(__CINT__)||  defined(__MAKECINT__)
#include "TRandom3.h"
#endif


class RndmGenerator {
public:
  /// Constructor
  RndmGenerator(unsigned int seed = 1);

  /// Destructor
  virtual ~RndmGenerator();

  //Operations
  void setSeed(unsigned int newSeed);
  unsigned int seed() const;
  double extractFlat(double min = 0, double max = 1) const;
  double extractThetaDistrib(double n) const;
  int extractPoisson(double mean) const;
  
  int extractInteg(int max) const;


protected:

private:

  TRandom3* rRandom;
  unsigned int theSeed;

  double thetaDistrib(double theta, double n) const;

  ClassDef(RndmGenerator,1);
};
#endif

