#ifndef SIMPLE_MODEL_HELPER_H
#define SIMPLE_MODEL_HELPER_H

#if defined(HAVE_CONFIG_H)
#include "config.h"
#endif

#include <inttypes.h>
#include <math.h>
#include <time.h>

enum
{
  /* 129 is the first user-defined code of the GRIB format */
  GRIB_USERDEF = 129,
};

void var_scale(int datatype, double *mscale, double *mrscale);

#define DEG2RAD(phi) ((M_PI / 180.0) * (phi))
#define RAD2DEG(phi) ((180. / M_PI) * (phi))

static inline double
sign_flat(double v)
{
#ifdef _CRAYC
  if (v == 0 && copysign(1.0, v) < 0.0) return 0.0;
#else
  if (!(v < 0.0 || v > 0.0)) return 0.0;
#endif
  //    if (v == 0.0) return 0.0;
  return v;
}

/* data generator function */
static inline double
dg_wobble(double frac_x, double frac_y, double mscale, double mrscale)
{
  double r = sign_flat(round((cos(2.0 * M_PI * frac_x) * sin(2.0 * M_PI * frac_y)) * mscale)) * mrscale;
  return r;
}

time_t cditime2time_t(int date, int timeofday);
void time_t2cditime(time_t t, int *date, int *timeofday);

void composeFilename(char **buf, const char *fname_prefix, int tfID, const char *suffix);

int composeStream(char **buf, const char *fname_prefix, int tfID, const char *suffix, int filetype);

int createGlobalLatLonGrid(int nlon, int nlat);

int createLocalCurvilinearGrid(int sizex, int sizey);

struct cart_coord
{
  double lat, lon;
};

/*
 * coords_: points to array of size xyRange[0][0] * xyRange[0][1] * 2, where
 * the first half of xyRange[0][0] * xyRange[0][1] elements will be
 * initialized to longitude values and the second half to the latitude values
 * a: cartesian coordinates of rectangular part of globe in
 * anti-clockwise order
 * sizex, sizey: number of grid points in x and y directions respectively
 * xyRange: subset of the ranges [0,sizex) and [0,sizey) to compute,
 * given as start in xyRange[0] and count values in xyRange[1] for x
 * at index 0 and y at index 1
 */
void computeCurvilinearChunk(double *coords_, const struct cart_coord a[4], size_t sizex, size_t sizey, const size_t xyRange[2][2]);

#if defined USE_MPI && !defined HAVE_PPM_CORE
struct PPM_extent
{
  int32_t first, size;
};

struct PPM_extent PPM_uniform_partition(struct PPM_extent set_interval, int nparts, int part_idx);

int PPM_prime_factorization_32(uint32_t n, uint32_t **factors);

#endif

#if defined USE_MPI && defined HAVE_PPM_CORE
#include <ppm/ppm.h>
#include <core/ppm_random.h>
#ifdef _OPENMP
#define cdi_omp_pragma(prg) _Pragma(prg)
#else
#define cdi_omp_pragma(prg)
#endif
#define cdi_seed_repeatable_random(seed)                                                                                \
  do                                                                                                                    \
    {                                                                                                                   \
      int seed_ = (int) seed;                                                                                           \
      MPI_Comm comm_ = MPI_COMM_WORLD;                                                                                  \
      cdi_omp_pragma("omp parallel") if ((!PPM_initialized() || PPM_finalized())) PPM_initialize(&comm_, &seed_, NULL); \
      else PPM_ya_rand_init(MPI_COMM_WORLD, seed_);                                                                     \
    }                                                                                                                   \
  while (0)
#define cdi_repeatable_random() ((long) (PPM_irandp()))
#define cdi_repeatable_finalize() PPM_finalize()

#else
/* add random-wrapper to make call results repeatable, if PPM_random
 * can't be used */
void cdi_seed_repeatable_random(unsigned seed);
long cdi_repeatable_random(void);
#define cdi_repeatable_finalize() \
  do                              \
    {                             \
    }                             \
  while (0)
#endif

#endif
