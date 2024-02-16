#if defined(HAVE_CONFIG_H)
#include "config.h"
#endif

#include <assert.h>
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined USE_MPI && defined HAVE_PPM_CORE
#include <ppm/ppm_uniform_partition.h>
#include <core/ppm_combinatorics.h>
#endif

#include "cdi.h"
#include "dmemory.h"
#include "error.h"

#include "simple_model_helper.h"

void
var_scale(int datatype, double *mscale, double *mrscale)
{
  int mant_bits;
  switch (datatype)
    {
    case CDI_DATATYPE_PACK8: mant_bits = 7; break;
    case CDI_DATATYPE_PACK16: mant_bits = 15; break;
    case CDI_DATATYPE_PACK24: mant_bits = 23; break;
    case CDI_DATATYPE_FLT32: mant_bits = 24; break;
    case CDI_DATATYPE_FLT64: mant_bits = 53; break;
    case CDI_DATATYPE_INT8:
    case CDI_DATATYPE_INT16:
    case CDI_DATATYPE_INT32:
    default: fprintf(stderr, "Unexpected or unusable content format: %d\n", datatype); exit(EXIT_FAILURE);
    }
  *mscale = (double) (INT64_C(1) << mant_bits);
  *mrscale = 1.0 / *mscale;
}

/**
 * Compute UNIX epoch-based time_t from CDI's decimal encoding of date.
 */
time_t
cditime2time_t(int date, int timeofday)
{
  struct tm t_s;
  time_t t;
  t_s.tm_year = date / 10000;
  t_s.tm_mon = (date - t_s.tm_year * 10000) / 100;
  t_s.tm_mday = date % 100;
  t_s.tm_year -= 1900;
  t_s.tm_hour = timeofday / 10000;
  t_s.tm_min = (timeofday % 10000) / 100;
  t_s.tm_sec = timeofday % 100;
  t_s.tm_isdst = 0;
  t = mktime(&t_s);
  return t;
}

/**
 * Build decimal encoding of date from UNIX epoch-based time_t.
 */
void
time_t2cditime(time_t t, int *date, int *timeofday)
{
  struct tm *t_s;
  t_s = localtime(&t);
  *date = (t_s->tm_year + 1900) * 10000 + t_s->tm_mon * 100 + t_s->tm_mday;
  *timeofday = t_s->tm_hour * 10000 + t_s->tm_min * 100 + t_s->tm_sec;
}

void
composeFilename(char **buf, const char *fname_prefix, int tfID, const char *suffix)
{
  size_t fnSize = strlen(fname_prefix) + (tfID >= 0 ? 1 + sizeof(int) * 3 : 0) + 1 + strlen(suffix) + 1;
  char *filename = Realloc(*buf, fnSize);
  int plen = tfID >= 0 ? snprintf(filename, fnSize, "%s_%d.%s", fname_prefix, tfID, suffix)
                       : snprintf(filename, fnSize, "%s.%s", fname_prefix, suffix);
  if ((size_t) plen >= fnSize)
    {
      fprintf(stderr,
              "unexpected error: printing to string of size %zu"
              " results in %d chars to write\n",
              fnSize, plen);
      abort();
    }
  filename[fnSize - 1] = 0;
  *buf = filename;
}

int
composeStream(char **buf, const char *fname_prefix, int tfID, const char *suffix, int filetype)
{
  char *temp = buf ? *buf : NULL;
  composeFilename(&temp, fname_prefix, tfID, suffix);
  int streamID = streamOpenWrite(temp, filetype);
  if (streamID < 0)
    {
      fprintf(stderr, "Failed to open stream: %s\n", cdiStringError(streamID));
      abort();
    }
  if (!buf)
    free(temp);
  else
    *buf = temp;
  return streamID;
}

int
createGlobalLatLonGrid(int nlon, int nlat)
{
  int gridID = gridCreate(GRID_LONLAT, nlon * nlat);
  gridDefXsize(gridID, nlon);
  gridDefYsize(gridID, nlat);
  size_t maxAxisSize = (size_t) (nlon > nlat ? nlon : nlat);
  double *restrict coords = (double *) Malloc(maxAxisSize * sizeof(coords[0]));
  {
    double step = 360.0 / (double) nlon, ofs = step * 0.5;
    for (size_t i = 0; i < (size_t) nlon; ++i) coords[i] = (double) i * step + ofs;
    gridDefXvals(gridID, coords);
  }
  {
    double step = 180.0 / (double) nlat, ofs = step * 0.5;
    for (size_t i = 0; i < (size_t) nlat; ++i) coords[i] = (double) i * step - 90.0 + ofs;
    gridDefYvals(gridID, coords);
  }
  Free(coords);
  return gridID;
}

int
createLocalCurvilinearGrid(int sizex, int sizey)
{
  size_t gridsize = (size_t) sizex * (size_t) sizey;
  int gridID = gridCreate(GRID_CURVILINEAR, (int) gridsize);
  gridDefXsize(gridID, sizex);
  gridDefYsize(gridID, sizey);
  {
    /* anti-clockwise coordinates around Amazonia */
    static const struct cart_coord region[4]
#ifdef __cplusplus
        = { { DEG2RAD(-25.0), DEG2RAD(-85.0) },
            { DEG2RAD(-18.0), DEG2RAD(-44.0) },
            { DEG2RAD(-7.0), DEG2RAD(-50.0) },
            { DEG2RAD(10.0), DEG2RAD(-80.0) } };
#else
        = { { .lon = DEG2RAD(-85.0), .lat = DEG2RAD(-25.0) },
            { .lon = DEG2RAD(-44.0), .lat = DEG2RAD(-18.0) },
            { .lon = DEG2RAD(-50.0), .lat = DEG2RAD(7.0) },
            { .lon = DEG2RAD(-80.0), .lat = DEG2RAD(10.0) } };
#endif
    double(*gridCoords)[sizey][sizex] = (double(*)[sizey][sizex]) Malloc(sizeof(*gridCoords) * gridsize * 2);
    {
      const size_t xyRange[2][2] = { { 0, 0 }, { (size_t) sizex, (size_t) sizey } };
      computeCurvilinearChunk((double *) gridCoords, region, (size_t) sizex, (size_t) sizey, xyRange);
    }
    gridDefXvals(gridID, (double *) (gridCoords[1]));
    gridDefYvals(gridID, (double *) (gridCoords[0]));
    Free(gridCoords);
  }
  return gridID;
}

static inline double
cartDistance(struct cart_coord p1, struct cart_coord p2)
{
  double d_lat = sin((p1.lat - p2.lat) / 2.0), d_lon = sin((p1.lon - p2.lon) / 2.0),
         d = 2.0 * asin(sqrt(d_lat * d_lat + cos(p1.lat) * cos(p2.lat) * (d_lon * d_lon)));
  return d;
}

static inline struct cart_coord
intermediateCoord(struct cart_coord p1, struct cart_coord p2, double f)
{
  double d = cartDistance(p1, p2), sine_of_d = sin(d), A = sin((1 - f) * d) / sine_of_d, B = sin(f * d) / sine_of_d,
         x = A * cos(p1.lat) * cos(p1.lon) + B * cos(p2.lat) * cos(p2.lon),
         y = A * cos(p1.lat) * sin(p1.lon) + B * cos(p2.lat) * sin(p2.lon), z = A * sin(p1.lat) + B * sin(p2.lat);
  struct cart_coord ic = { .lat = atan2(z, sqrt(x * x + y * y)), .lon = atan2(y, x) };
  return ic;
}

void
computeCurvilinearChunk(double *coords_, const struct cart_coord a[4], size_t sizex, size_t sizey, const size_t xyRange[2][2])
{
  size_t startx = xyRange[0][0], starty = xyRange[0][1], chunkSizeX = xyRange[1][0], chunkSizeY = xyRange[1][1],
         endx = startx + chunkSizeX, endy = starty + chunkSizeY;
  assert(startx <= endx && endx <= sizex && starty <= endy && endy <= sizey);
#ifdef __cplusplus
  auto coords = (double(*)[chunkSizeY][chunkSizeX]) coords_;
#else
  double(*coords)[chunkSizeY][chunkSizeX] = (double(*)[sizey][sizex]) coords_;
#endif
  for (size_t j = starty; j < endy; ++j)
    {
      double g = (double) j / (double) (sizey - 1);
      /* compute start/end coordinates of great circle in x direction */
      struct cart_coord gc_left = intermediateCoord(a[0], a[3], g), gc_right = intermediateCoord(a[1], a[2], g);
      for (size_t i = startx; i < endx; ++i)
        {
          double f = (double) i / (double) (sizex - 1);
          struct cart_coord pij = intermediateCoord(gc_left, gc_right, f);
          coords[0][j - starty][i - startx] = RAD2DEG(pij.lat);
          coords[1][j - starty][i - startx] = RAD2DEG(pij.lon);
        }
    }
}

#if defined USE_MPI && !defined HAVE_PPM_CORE
static int32_t uniform_partition_start(struct PPM_extent set_interval, int nparts, int part_idx);

struct PPM_extent
PPM_uniform_partition(struct PPM_extent set_interval, int nparts, int part_idx)
{
  struct PPM_extent range;
  range.first = uniform_partition_start(set_interval, nparts, part_idx);
  range.size = uniform_partition_start(set_interval, nparts, part_idx + 1) - range.first;
  return range;
}

static int32_t
uniform_partition_start(struct PPM_extent set_interval, int nparts, int part_idx)
{
  int32_t part_offset = ((int64_t) set_interval.size * (int64_t) part_idx) / (int64_t) nparts;
  int32_t start = set_interval.first + part_offset;
  return start;
}

int
PPM_prime_factorization_32(uint32_t n, uint32_t **factors)
{
  if (n <= 1) return 0;
  uint32_t *restrict pfactors = Realloc(*factors, 32 * sizeof(pfactors[0]));
  size_t numFactors = 0;
  uint32_t unfactored = n;
  while (!(unfactored & 1))
    {
      pfactors[numFactors] = 2;
      ++numFactors;
      unfactored >>= 1;
    }
  uint32_t divisor = 3;
  while (unfactored > 1)
    {
      while (unfactored % divisor == 0)
        {
          unfactored /= divisor;
          pfactors[numFactors] = divisor;
          ++numFactors;
        }
      divisor += 1;
    }
  *factors = pfactors;
  return numFactors;
}

#endif

#if !defined USE_MPI || !defined HAVE_PPM_CORE
static char repeatable_rand_state[31 * sizeof(long)];

long
cdi_repeatable_random(void)
{
  char *caller_rand_state = setstate(repeatable_rand_state);
  long retval = random();
  setstate(caller_rand_state);
  return retval;
}

void
cdi_seed_repeatable_random(unsigned seed)
{
  char *caller_rand_state = initstate(seed, repeatable_rand_state, sizeof(repeatable_rand_state));
  setstate(caller_rand_state);
}
#endif

/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
