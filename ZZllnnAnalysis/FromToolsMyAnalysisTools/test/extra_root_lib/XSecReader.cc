
/*
 *  See header file for a description of this class.
 *
 *  $Date: 2008/01/11 15:26:10 $
 *  $Revision: 1.3 $
 *  \author G. Cerminara - NEU Boston & INFN Torino
 */

#include "XSecReader.h"
#include "TString.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;

XSecReader::XSecReader(const TString& xsecFileName, const TString& lumiFileName) {
  ifstream xsecFile(xsecFileName.Data());
  string line;
  // Read the xsection file and store numbers in a map (per sample)
  while (getline(xsecFile,line)) {
    if( line == "" || line[0] == '#' ) continue; // Skip comments and empty lines
    stringstream linestr;
    linestr << line;
    TString finalState, sampleName;
    int nEvents;
    double xsec, br;
    linestr >> finalState >> sampleName >> nEvents >> xsec >> br;
    double weight = xsec * br / nEvents;
    cout << "sample: " << sampleName << " #ev: " << nEvents
	 << " xsec: " << xsec << " br: " << br << " weight: " << weight << endl;
    if(finalState == "diem") {
      theWeightMapDiem[sampleName] = weight;
      theInitialNEvDiem[sampleName] = nEvents;
    } else if(finalState == "dimu") {
      theWeightMapDimu[sampleName] = weight;
      theInitialNEvDimu[sampleName] = nEvents;
    } else {
      cerr << "[XSecReader]*** Error: final state in the map is invalid: " << finalState << endl;
    }
  }
  xsecFile.close();

  ifstream lumiFile(lumiFileName.Data());

  // Read the lumi map and store numbers in a map (per sample)
  while (getline(lumiFile,line)) {
    if( line == "" || line[0] == '#' ) continue; // Skip comments and empty lines
    stringstream linestr;
    linestr << line;
    TString finalState;
    TString epoch;
    bool dqApplied;
    double lumi = -1;
    linestr >> finalState >> epoch >> dqApplied >> lumi;
    cout << "final state: " << finalState << " epoch: " << epoch << " DQ applied: " <<
	 dqApplied << " lumi: " << lumi << endl;
    if(dqApplied) {
      theLumiMapDQ[epoch][finalState] = lumi;
    } else {
      theLumiMapNoDQ[epoch][finalState] = lumi;
    }
  }
  lumiFile.close();

}


XSecReader::~XSecReader(){}


double XSecReader::getWeight(const TString& sampleName, const TString& epoch,
			     const TString& finalState, const bool dqApplied) const {
  const map<TString, double> *theWeightMap = 0;
  if(finalState == "diem") {
    theWeightMap = &theWeightMapDiem;
  } else if(finalState == "dimu") {
    theWeightMap = &theWeightMapDimu;
  } else {
    cerr << "[XSecReader]*** Error: final state is invalid: " << finalState << endl;
  }


  if(theWeightMap->find(sampleName) == theWeightMap->end() ||
     theLumiMapDQ.find(epoch) == theLumiMapDQ.end() ||
     theLumiMapDQ.find(epoch)->second.find(finalState) == theLumiMapDQ.find(epoch)->second.end()) {
//     if(theWeightMap->find(sampleName) == theWeightMap->end()) 
//       cout << "[XSecReader]***ERROR: invalid sample name: " << sampleName << endl;
    cout << "[XSecReader]***ERROR: invalid labels; sampleName: " << sampleName
	 << " epoch: " << epoch
	 << " finalState: " << finalState
	 << " dqApplied: " << dqApplied << endl;
    return -1;
  }
  double lumi = -1;
  if(dqApplied) {
    lumi = ((theLumiMapDQ.find(epoch)->second).find(finalState))->second;
  } else {
    lumi = ((theLumiMapNoDQ.find(epoch)->second).find(finalState))->second;
  }
//   cout << "[XSecReader] lumi: " << lumi
//        << " weight: " << theWeightMap->find(sampleName)->second << endl;

  return lumi*theWeightMap->find(sampleName)->second;
}


double XSecReader::getInitNEv(const TString& sampleName, const TString& finalState) const {
  const map<TString, double> *theNEvMap = 0;
  if(finalState == "diem") {
    theNEvMap = &theInitialNEvDiem;
  } else if(finalState == "dimu") {
    theNEvMap = &theInitialNEvDimu;
  } else {
    cerr << "[XSecReader]*** Error: final state is invalid: " << finalState << endl;
  }
  
  map<TString, double>::const_iterator elem = theNEvMap->find(sampleName);
  
  if(elem == theNEvMap->end()) {
    cout << "[XSecReader]***ERROR: invalid labels; sampleName: " << sampleName << endl;
    return -1;
  } else {
    return elem->second;
  }
}


double XSecReader::getLuminosity(const TString& epoch, const TString& finalState,
				 const bool dqApplied) const {
  double ret = -1;
  if(dqApplied) {
    ret = theLumiMapDQ.find(epoch)->second.find(finalState)->second;  
  } else {
    ret = theLumiMapNoDQ.find(epoch)->second.find(finalState)->second;
  }
  return ret;
}
