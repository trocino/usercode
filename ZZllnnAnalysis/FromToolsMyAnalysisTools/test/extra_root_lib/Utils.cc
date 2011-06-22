
/*
 *  See header file for a description of this class.
 *
 *  $Date: 2007/12/04 00:11:38 $
 *  $Revision: 1.2 $
 *  \author G. Cerminara - NEU Boston & INFN Torino
 */

#include "Utils.h"
#include "TMath.h"

#include <iostream>

using namespace std;

Utils::Utils(){}

Utils::~Utils(){}

double Utils::getSignifFromMinNLL(double minNLLSig, double minNLLBck) {
  return getSignifFromMaxL(TMath::Exp(-minNLLSig),TMath::Exp(-minNLLBck));
}


double Utils::getSignifFromMaxL(double maxLSig, double maxLBck) {

  double ratio = 2*TMath::Log(maxLSig/maxLBck);
  if(TMath::IsNaN(TMath::Sqrt(ratio))) {
//     cout << "#############################################" << endl;
//     cout << "maxLSig: " <<  maxLSig << " maxLBck: " << maxLBck
// 	 << " TMath::Log(maxLSig/maxLBck): " << TMath::Log(maxLSig/maxLBck)
// 	 << " ratio: " << ratio << " TMath::Sqrt(ratio): " << TMath::Sqrt(ratio) << endl;
    return 0.;
  }


  return TMath::Sqrt(ratio);

}
  


