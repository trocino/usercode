#ifndef Utils_H
#define Utils_H

/** \class Utils
 *  Utility functions.
 *
 *  $Date: 2007/08/10 19:15:15 $
 *  $Revision: 1.1 $
 *  \author G. Cerminara - NEU Boston & INFN Torino
 */

class Utils {
public:
  /// Constructor
  Utils();

  /// Destructor
  virtual ~Utils();

  // Operations
  static double getSignifFromMinNLL(double minNLLSig, double minNLLBck);
  static double getSignifFromMaxL(double maxLSig, double maxLBck);
  

protected:

private:

};
#endif

