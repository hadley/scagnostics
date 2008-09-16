/*
 * Binner
 *
 * Leland Wilkinson (SPSS, Inc.)
 *
 * Permission to use, copy, modify, and distribute this software for any
 * purpose without fee is hereby granted, provided that this entire notice
 * is included in all copies of any software which is or includes a copy
 * or modification of this software and in all copies of the supporting
 * documentation for such software.
 * THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED
 * WARRANTY.  IN PARTICULAR, THE AUTHORS MAKE NO
 * REPRESENTATION OR WARRANTY OF ANY KIND CONCERNING THE MERCHANTABILITY
 * OF THIS SOFTWARE OR ITS FITNESS FOR ANY PARTICULAR PURPOSE.
 */
#ifndef BINNER_H
#define BINNER_H

#include <cmath>
#include <iostream>

using namespace std;

class BinnedData
{

  static inline int *integerizeData (int n, double *x)
  {
    int *xd = new int[n];
    for (int i = 0; i < n; ++i)
    {
      xd[i] = (int) (1000 * x[i]);
    }
    return xd;
  }

public:
  int n;
  int *x, *y;
  int *counts;

  static const double RESOLUTION = 1000;
  static const int BINS = 40;

  inline BinnedData (int n, double *x, double *y, int *counts)
  {
    this->x = integerizeData (n, x);
    this->y = integerizeData (n, y);
    this->counts = counts;
    this->n = n;
  }

  inline int getNBins ()
  {
    return n;
  }

  inline int *getXData ()
  {
    return x;
  }

  inline int *getYData ()
  {
    return y;
  }

  inline int *getCounts ()
  {
    return counts;
  }

  inline void show ()
  {
    for (int i = 0; i < n; i++) {
      printf("(%d, %d, %d)\n", x[i], y[i], counts[i]);
    }
  };

};                              //end of class

class Binner
{
  static inline int deleteEmptyBins (int nBin, int *count, double *xbin,
                                     double *ybin)
  {
    int k = 0;
    for (int i = 0; i < nBin; ++i)
    {
      if (count[i] > 0) {
        count[k] = count[i];
        xbin[k] = xbin[i];
        ybin[k] = ybin[i];
        ++k;
      }
    }
    return k;
  };


public:
  Binner () {
  }

  static BinnedData binHex (int n, double *x, double *y, const int nBins);
};

#endif
