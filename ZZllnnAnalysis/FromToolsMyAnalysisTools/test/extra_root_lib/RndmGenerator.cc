#include "RndmGenerator.h"

#if !defined(__CINT__)||  defined(__MAKECINT__)
#include "TMath.h"
#endif


ClassImp(RndmGenerator);

// TRandom3* RndmGenerator::rRandom = new TRandom3(1);


RndmGenerator::RndmGenerator(unsigned int seed) : theSeed(seed) {
  rRandom = new TRandom3(theSeed);
 }

RndmGenerator::~RndmGenerator() {}

void RndmGenerator::setSeed(unsigned int newSeed) {
  theSeed = newSeed;
  rRandom->SetSeed(newSeed);
}
unsigned int RndmGenerator::seed() const {
  return theSeed;
}

double RndmGenerator::extractFlat(double min, double max) const {
  //FIXME: Argomento Rndm = luxury level???
  return min+rRandom->Rndm()*(max-min); 
}



double RndmGenerator::thetaDistrib(double theta, double n) const {
  return TMath::Power(TMath::Cos(theta), n);
}

double RndmGenerator::extractThetaDistrib(double n) const {

  static const double PI = TMath::Pi();
  static const double epsilon = 0.1;
  double thetaMax = PI/2;
  double maxDistrib = 1+epsilon;
  if(n<0) {
    thetaMax = 0.6; //Cut-off from geometry
    maxDistrib = thetaDistrib(thetaMax, n);
  }
  
  for(;;) {
    double u1 = extractFlat(0, thetaMax);
    double u2 = extractFlat();
    if(thetaDistrib(u1,n) > maxDistrib*u2)
	return u1;
  }
}

int RndmGenerator::extractInteg(int max) const {
  return rRandom->Integer(max);
}

int RndmGenerator::extractPoisson(double mean) const {
  return rRandom->Poisson(mean);
}
