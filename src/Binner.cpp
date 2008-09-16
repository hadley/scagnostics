#include "Binner.h"

BinnedData Binner::binHex(int n, double *x, double *y, const int nBins) {

  // scaling constants

  double con1 = .25;
  double con2 = 1. / 3.;
  double c1 = (double) (nBins - 1);
  double c2 = c1 / sqrt(3.);
  int jinc = nBins;
  int iinc = 2 * nBins;
  int nBin = (nBins + 20) * (nBins + 20);

  int *count = new int[nBin];
  double *xbin = new double[nBin];
  double *ybin = new double[nBin];

  // fill bins
  for(int i = 0; i < nBin; ++i) {
    count[i] = 0;
    xbin[i] = 0;
    ybin[i] = 0;
  }

  for (int i = 0; i < n; ++i) {
    double sx = c1 * x[i];
    double sy = c2 * y[i];
    int i1 = (int) (sy + .5);
    int j1 = (int) (sx + .5);
    double dy = sy - ((double) i1);
    double dx = sx - ((double) j1);
    double dist1 = dx * dx + 3. * dy * dy;
    int m = 0;
    if (dist1 < con1) {
      m = i1 * iinc + j1;
    } else if (dist1 > con2) {
      m = ((int) sy) * iinc + ((int) sx) + jinc;
    } else {
      int i2 = (int) sy;
      int j2 = (int) sx;
      dy = sy - ((double) i2) - .5;
      dx = sx - ((double) j2) - .5;
      double dist2 = dx * dx + 3. * dy * dy;
      if (dist1 <= dist2) {
        m = i1 * iinc + j1;
      } else {
        m = i2 * iinc + j2 + jinc;
      }
    }
    ++count[m];
    xbin[m] += (x[i] - xbin[m]) / count[m];
    ybin[m] += (y[i] - ybin[m]) / count[m];
  }


  nBin = deleteEmptyBins(nBin, count, xbin, ybin);
  if (nBin > 250) {
    BinnedData b = binHex(n, x, y, nBins/2);
    return b;
  }

  int *tcount = new int[nBin];
  double *xtbin = new double[nBin];
  double *ytbin = new double[nBin];

  for (int i =0; i < nBin; ++i) {    // copy arrays by value
    tcount[i] = count[i];
    xtbin[i] = xbin[i];
    ytbin[i] = ybin[i];
  }
  BinnedData bdata(nBin, xtbin, ytbin, tcount); 
  return bdata;
}

