#ifndef LumiNorm_H
#define LumiNorm_H

/** \class LumiNorm
 *  Class to get the normalization of each sample.
 *  The relative normalization is get from MC xsections while the overall normalization
 *  can be derived from the number of events under the Z peak.
 *  The samples used in this luminosity normalization must be added explicitly throught the add method.
 *  Look to the XSection.txt and Luminosity.txt files for naming conventions.
 *
 *  $Date: 2008/03/13 17:41:58 $
 *  $Revision: 1.3 $
 *  \author G. Cerminara - CERN
 *          D. Trocino   - Northeastern University
 */

#include "TString.h"
#include <vector>
#include <map>

class XSecReader;

class LumiNorm {
public:
  /// Constructor
  // Create a tool for normalization. You need to specify the XSection file and the lumi file
  // and the elements to look for the right lumi in the map.
  // also the inputDir for the histo to be used for data/MC normalization in required (che be a null string)
  /*
  LumiNorm(const TString& inputDir, 
	   const TString& lumiFileName, 
	   const TString& finalState);
  */
  LumiNorm(const TString& inputDir, const TString& lumiFileName,
	   bool rescaleLumi=false, bool rescaleXsect=false);

  /// Destructor
  virtual ~LumiNorm();

  // Operations
  // Add data to the count of events to be used for data/MC normalization
  // the histo used to get the count is zmass_restricted_1
  // the default file for this histo is: histos_zzanalysis_data_FINALSTATE.root
  void addData();

  // Add data to the count of events to be used for data/MC normalization
  // the histo used to get the count is zmass_restricted_1
  // This method allows to specify the name of the file which contains this histo
  void addData(const TString& filename);

  // Add MC sample to the count of events to be used for data/MC normalization
  // the histo used to get the count is zmass_restricted_1
  // the default file for this histo is: theInputDir+"/histos_zzanalysis_mc_"+sampleName+"_"+theFinalState+".root"
  void addMC(const TString& sampleName);

  // Add MC sample to the count of events to be used for data/MC normalization
  // the histo used to get the count is zmass_restricted_1
  // This method allows to specify the name of the file which contains this histo
  void addMC(const TString& sampleName, const TString& filename);

  // Add MC samples to the count of events to be used for data/MC normalization
  // the histo used to get the count is zmass_restricted_1
  // the default file for this histo is: theInputDir+"/histos_zzanalysis_mc_"+sampleName+"_"+theFinalState+".root"
  void addMC(const std::vector<TString>& sampleNames);


  // Get the relative normalization of MC/data under the Z peak
  double getNormalizationFactor() const;

  double getLuminosity() const;

  // Get initial # of events (as in .root file)  // FIXME
  float getInitialNEv(const TString& sampleName) const;
  
  // Get cross section (as in .root file)  // FIXME
  float getCrossSection(const TString& sampleName) const;
  
  // Get branching ratio (as in .root file)  // FIXME
  float getBranchingRatio(const TString& sampleName) const;
  
  // Get the overall scale factor (relative normalization of MC/data under the Z peak + theo. xsec)
  float getScaleFactor(const TString& sampleName) const;

  // add an additional scale facto
  void setAdditionalScale(double scale);

  // set the name of the file histo to be used for the normalization of MC/data under the Z peak
  void setNormalizHistoName(const TString& histoname);

  // switch off the normalization of MC/data under the Z peak
  void normalizeToZPeak(bool doNormalize);

  // Maps
  // making them public for easier access
  std::map<TString, TString> conversionMap;
  std::map<TString, float> initNEvMap;
  std::map<TString, float> xSectMap;
  std::map<TString, float> brRatioMap;

protected:

private:
  double nData;
  double nMC;
  double theLuminosity;
  TString theInputDir;
  TString theLumiFileName;
  /*TString theFinalState;*/
  TString theLumiHistoName;
  TString theNormHistoName;
  double additionalScale;
  bool normalizeToZ;
  bool rescaleLuminosity;
  bool rescaleCrossSect;

};
#endif

