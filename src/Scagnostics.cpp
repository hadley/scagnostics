#include <cstdio>
#include <fstream>
#include <string>
#include <cfloat>
#include <ctime>
#include <math.h>

#include "StringTokenizer.h"
#include "Scagnostics.h"
#include "Binner.h"
#include "GraphMeasures.h"

using namespace std;

const char *Scagnostics::scagnosticsLabels[9] = {"Outlying", "Skewed", "Clumpy", "Sparse", "Striated", "Convex", "Skinny", "Stringy", "Monotonic"};
// == Triangulation::scagnosticsLabels;


void Scagnostics::computeScagnosticsOutliers() {
  isScagnosticOutlier = computeMSTOutliers(scagnostics);
}

// for using qsort
double *X;    // needs to be set before using icompare
int icompare (const void * a, const void * b)
{
  if (X[*(int*)a] < X[*(int*)b]) return -1;
  if (X[*(int*)a] > X[*(int*)b]) return 1;
  
  return 0;
}


bool *Scagnostics::computeMSTOutliers(double **pts) {

  // Prim's algorithm on simple integer arrays

  int nVar = Triangulation::numScagnostics;                  // 9 Scagnostic dimensions
  int nPts = numVars * (numVars - 1) / 2;               // p*(p-1)/2 points representing scatterplots
  int *edges[2];
    edges[0] = new int[nPts - 1];
    edges[1] = new int[nPts - 1];
  int *list = new int[nPts];
  int *degrees = new int[nPts];
  double *cost = new double[nPts];
  double *lengths = new double[nPts - 1];

  list[0] = 0;
  cost[0] = DBL_MAX;
  int cheapest = 0;

  for (int i = 1; i < nPts; ++i) {
    for (int j = 0; j < nVar; ++j) {
      double d = pts[j][i] - pts[j][0];
      cost[i] += d * d;
    }
    if (cost[i] < cost[cheapest])
      cheapest = i;
  }
  for (int j = 1; j < nPts; ++j) {
    int end = list[cheapest];
    int jp = j - 1;
    edges[0][jp] = cheapest;
    edges[1][jp] = end;
    lengths[jp] = cost[cheapest];
    ++degrees[cheapest];
    ++degrees[end];
    cost[cheapest] = DBL_MAX;
    end = cheapest;

    for (int i = 1; i < nPts; ++i) {
      if (cost[i] != DBL_MAX) {
        double dist = 0.;
        for (int k = 0; k < nVar; ++k) {
          double d = (pts[k][i] - pts[k][end]);
          dist += d * d;
        }
        if (dist < cost[i]) {
          list[i] = end;
          cost[i] = dist;
        }
        if (cost[i] < cost[cheapest]) cheapest = i;
      }
    }
  }

  // find cutoff value for identifying extremely long edges

  bool *outlier = new bool[nPts];

//  int *index = Sort.indexedDoubleArraySort(lengths, 0, 0);
int *index = new int[nPts - 1];  
for (int k=0; k < nPts-1; ++k) index[k] = k;
X = lengths;
qsort(index, nPts-1, sizeof(int), icompare);

  int n50 = (nPts-1) / 2;
  int n25 = n50 / 2;
  int n75 = n50 + n50 / 2;
  double cutoff = lengths[index[n75]] + 1.5 * (lengths[index[n75]] - lengths[index[n25]]);
  for (int i = 0; i < nPts - 1; ++i) {
    if (lengths[i] > cutoff) {
      for (int k = 0; k < 2; ++k) {
        int node = edges[k][i];
        if (degrees[node] == 1) { // outliers are nodes of degree 1 with long adjacent edges
          outlier[node] = true;
        }
      }
    }
  }

  return outlier;
}

