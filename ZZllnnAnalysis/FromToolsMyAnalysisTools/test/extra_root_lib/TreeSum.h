#ifndef TreeSum_H
#define TreeSum_H

/** \class TreeSum
 *  No description available.
 *
 *  $Date: 2007/08/10 19:15:15 $
 *  $Revision: 1.1 $
 *  \author G. Cerminara - NEU Boston & INFN Torino
 */

class TTree;
class TNtuple;

class TreeSum {
public:
  /// Constructor
  TreeSum(TTree *sigTree, TTree *bckTree);
  
  

  /// Destructor
  virtual ~TreeSum();

  // Operations
  TNtuple *getSumTree();

  void setScaleFactors(double scale1, double scale2);

protected:

private:
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

  TTree *theSigTree;
  TTree *theBckTree;

  double sigScale;
  double bckScale;

};
#endif

