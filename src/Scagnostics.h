#include <iostream>
#include <fstream>
#include <string>
#include <cfloat>
#include "StringTokenizer.h"
#include "GraphMeasures.h"
using namespace std;

class Scagnostics {
  int MAX_ROWS;
  int numVars, numRows;
  int numScagnostics;
  bool *isScagnosticOutlier; // = null;

  double *dataMin, *dataMax;

  void getData(int argc, char * const argv[]);
  bool getFileData(string fname);
  void computeOnFileData(int argc, char * const argv[]);
  void computeScagnosticsOutliers();
  bool *computeMSTOutliers(double **pts);
  
  inline void initializeMinMax() {
    dataMin = new double[numVars];
    dataMax = new double[numVars];
    for (int i = 0; i < numVars; ++i) {
      dataMin[i] = DBL_MAX;
      dataMax[i] = -dataMin[i];
    }
  };

  inline void updateMinMax(double d, int i) {
    if (d < dataMin[i])
      dataMin[i] = d;
    if (d > dataMax[i])
      dataMax[i] = d;
  };

  inline void normalizeData(){
    for (int i = 0; i < numVars; ++i) {
      for (int j = 0; j < numRows; ++j) {
        data[i][j] = (data[i][j] - dataMin[i]) / (dataMax[i] - dataMin[i]);
      }
    }
  };

  public:
  double **data;
  double **scagnostics;
  string  *variableLabels;
  static const char* scagnosticsLabels[9];

  Scagnostics(int argc, char * const argv[]) {
    data = NULL;
    variableLabels = NULL;
    scagnostics = NULL;
    isScagnosticOutlier = NULL;

    MAX_ROWS = 100;
    numVars = 0;
    numRows = 0; 
    numScagnostics = Triangulation::numScagnostics;
    computeOnFileData(argc, argv);
  }
};


