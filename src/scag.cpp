#include "Scagnostics.h"

#include <Rdefines.h>

extern "C" {
  void
  scagnostics(double* x, double* y, int* length, int* bins, double* results) {
    Binner b;

    BinnedData bdata = b.binHex(length[0], x, y, bins[0]);
    Triangulation dt;

    // bdata.show();

    double* r = dt.compute(bdata, false);
    results[9] = bdata.n;
    memcpy(results, r, sizeof(double) * 9);
    for (int i = 0; i < bdata.n; i++) {
      results[10 + 0 * bdata.n + i] = bdata.x[i];
      results[10 + 1 * bdata.n + i] = bdata.y[i];
      results[10 + 2 * bdata.n + i] = bdata.counts[i];
    }
  }

}
