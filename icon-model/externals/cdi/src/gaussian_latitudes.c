#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>

#ifndef M_SQRT2
#define M_SQRT2 1.41421356237309504880168872420969808
#endif

static void
cpledn(size_t kn, size_t kodd, double *pfn, double pdx, int kflag, double *pw, double *pdxn, double *pxmod)
{
  // 1.0 Newton iteration step

  double zdlx = pdx;
  double zdlldn = 0.0;

  size_t ik = 1;

  if (kflag == 0)
    {
      double zdlk = 0.5 * pfn[0];
      for (size_t jn = 2 - kodd; jn <= kn; jn += 2)
        {
          // normalised ordinary Legendre polynomial == \overbar{p_n}^0
          zdlk = zdlk + pfn[ik] * cos((double) (jn) *zdlx);
          // normalised derivative == d/d\theta(\overbar{p_n}^0)
          zdlldn = zdlldn - pfn[ik] * (double) (jn) *sin((double) (jn) *zdlx);
          ik++;
        }
      // Newton method
      double zdlmod = -(zdlk / zdlldn);
      double zdlxn = zdlx + zdlmod;
      *pdxn = zdlxn;
      *pxmod = zdlmod;
    }

  // 2.0 Compute weights

  if (kflag == 1)
    {
      for (size_t jn = 2 - kodd; jn <= kn; jn += 2)
        {
          // normalised derivative
          zdlldn = zdlldn - pfn[ik] * (double) (jn) *sin((double) (jn) *zdlx);
          ik++;
        }
      *pw = (double) (2 * kn + 1) / (zdlldn * zdlldn);
    }

  return;
}

static void
gawl(double *pfn, double *pl, double *pw, size_t kn)
{
  double pmod = 0.0;
  double zw = 0.0;
  double zdlxn = 0.0;

  // 1.0 Initizialization

  int iflag = 0;
  int itemax = 20;

  size_t iodd = (kn % 2);

  double zdlx = *pl;

  // 2.0 Newton iteration

  for (int jter = 1; jter <= itemax + 1; ++jter)
    {
      cpledn(kn, iodd, pfn, zdlx, iflag, &zw, &zdlxn, &pmod);
      zdlx = zdlxn;
      if (iflag == 1) break;
      if (fabs(pmod) <= DBL_EPSILON * 1000.0) iflag = 1;
    }

  *pl = zdlxn;
  *pw = zw;

  return;
}

static void
gauaw(size_t kn, double *restrict pl, double *restrict pw)
{
  /*
   * 1.0 Initialize Fourier coefficients for ordinary Legendre polynomials
   *
   * Belousov, Swarztrauber, and ECHAM use zfn(0,0) = sqrt(2)
   * IFS normalisation chosen to be 0.5*Integral(Pnm**2) = 1 (zfn(0,0) = 2.0)
   */
  double *zfn = (double *) malloc((kn + 1) * (kn + 1) * sizeof(double));
  double *zfnlat = (double *) malloc((kn / 2 + 1 + 1) * sizeof(double));

  zfn[0] = M_SQRT2;
  for (size_t jn = 1; jn <= kn; ++jn)
    {
      double zfnn = zfn[0];
      for (size_t jgl = 1; jgl <= jn; ++jgl)
        {
          zfnn *= sqrt(1.0 - 0.25 / ((double) (jgl * jgl)));
        }

      zfn[jn * (kn + 1) + jn] = zfnn;

      size_t iodd = jn % 2;
      for (size_t jgl = 2; jgl <= jn - iodd; jgl += 2)
        {
          zfn[jn * (kn + 1) + jn - jgl] = zfn[jn * (kn + 1) + jn - jgl + 2] * ((double) ((jgl - 1) * (2 * jn - jgl + 2)))
                                          / ((double) (jgl * (2 * jn - jgl + 1)));
        }
    }

  // 2.0 Gaussian latitudes and weights

  size_t iodd = kn % 2;
  size_t ik = iodd;
  for (size_t jgl = iodd; jgl <= kn; jgl += 2)
    {
      zfnlat[ik] = zfn[kn * (kn + 1) + jgl];
      ik++;
    }

  // 2.1 Find first approximation of the roots of the Legendre polynomial of degree kn

  size_t ins2 = kn / 2 + (kn % 2);

  for (size_t jgl = 1; jgl <= ins2; ++jgl)
    {
      double z = ((double) (4 * jgl - 1)) * M_PI / ((double) (4 * kn + 2));
      pl[jgl - 1] = z + 1.0 / (tan(z) * ((double) (8 * kn * kn)));
    }

  // 2.2 Computes roots and weights for transformed theta

  for (size_t jgl = ins2; jgl >= 1; --jgl)
    {
      size_t jglm1 = jgl - 1;
      gawl(zfnlat, &(pl[jglm1]), &(pw[jglm1]), kn);
    }

  // convert to physical latitude

  for (size_t jgl = 0; jgl < ins2; ++jgl) pl[jgl] = cos(pl[jgl]);

  for (size_t jgl = 1; jgl <= kn / 2; ++jgl)
    {
      size_t jglm1 = jgl - 1;
      size_t isym = kn - jgl;
      pl[isym] = -pl[jglm1];
      pw[isym] = pw[jglm1];
    }

  free(zfnlat);
  free(zfn);

  return;
}

void
gaussianLatitudes(size_t nlats, double *latitudes, double *weights)
{
  gauaw(nlats, latitudes, weights);
}

bool
isGaussianLatitudes(size_t nlats, const double *latitudes)
{
  bool is_gauss_lats = false;

  if (nlats > 2)  // check if gaussian
    {
      size_t i;
      double *yv = (double *) malloc(nlats * sizeof(double));
      double *yw = (double *) malloc(nlats * sizeof(double));
      gaussianLatitudes(nlats, yv, yw);
      free(yw);

      for (i = 0; i < nlats; i++) yv[i] = asin(yv[i]) / M_PI * 180.0;

      for (i = 0; i < nlats; i++)
        if (fabs(yv[i] - latitudes[i]) > ((yv[0] - yv[1]) / 500.0)) break;

      if (i == nlats) is_gauss_lats = true;

      // check S->N
      if (is_gauss_lats == false)
        {
          for (i = 0; i < nlats; i++)
            if (fabs(yv[i] - latitudes[nlats - i - 1]) > ((yv[0] - yv[1]) / 500.0)) break;

          if (i == nlats) is_gauss_lats = true;
        }

      free(yv);
    }

  return is_gauss_lats;
}
