#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "cdi.h"

#include "dmemory.h"

#include "simple_model_helper.h"

static double my_gamma_dist(double x);

#define missValue (-50.0)

int
main(int argc, const char **argv)
{
  // todo: handle optional arguments here to increase test coverage
  const char *fname = (argc > 1) ? argv[1] : "test.nc";

  int streamID = streamOpenWrite(fname, CDI_FILETYPE_NC);
  if (streamID < 0)
    {
      fprintf(stderr, "Open failed on %s: %s\n", fname, cdiStringError(streamID));
      return EXIT_FAILURE;
    }

  enum
  {
    sizey = 40,
    sizex = 2 * sizey,
  };

  int gridID = createLocalCurvilinearGrid(sizex, sizey);
  int zaxisID = zaxisCreate(ZAXIS_SURFACE, 1);

  int vlistID = vlistCreate();
  int varID = vlistDefVar(vlistID, gridID, zaxisID, TIME_VARYING);
  vlistDefVarMissval(vlistID, varID, missValue);
  {
    static const char creatorText[] = "CDI test_cdf_write";
    cdiDefAttTxt(vlistID, varID, "CDI Text Attribute test, created by", sizeof(creatorText) - 1, creatorText);
  }

  int taxisID = taxisCreate(TAXIS_ABSOLUTE);
  vlistDefTaxis(vlistID, taxisID);

  streamDefVlist(streamID, vlistID);

  (void) streamDefTimestep(streamID, 0);

  {
    double(*data)[sizex] = (double(*)[sizex]) Malloc(sizeof(**data) * sizex * sizey);
    for (size_t j = 0; j < sizey; ++j)
      for (size_t i = 0; i < sizex; ++i)
        {
          data[j][i] = my_gamma_dist((double) i / (double) (sizex - 1));
        }
    data[sizey / 3][sizex / 2] = missValue;
    streamWriteVar(streamID, 0, (const double *) data, 1);
    Free(data);
  }

  streamClose(streamID);

  return EXIT_SUCCESS;
}

static double
my_gamma_dist(double x)
{
  enum
  {
    k = 9,
  };
  const double theta = 0.5;
  x *= 20.0;
  double pdf_x = 1.0 / (tgamma(k) * pow(theta, k)) * pow(x, k - 1) * exp(-x / theta);
  return pdf_x;
}

/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
