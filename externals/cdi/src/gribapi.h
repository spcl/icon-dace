#ifndef GRIBAPI_H
#define GRIBAPI_H

#ifdef HAVE_LIBGRIB_API
#include <grib_api.h>
#ifndef ERROR_H
#include "error.h"
#endif
#endif

#ifndef CDI_INT_H
#include "cdi_int.h"
#endif

// clang-format off

#define  GRIBAPI_MISSVAL  -9.E33

// GRIB2 Level Types
#define  GRIB2_LTYPE_SURFACE               1
#define  GRIB2_LTYPE_CLOUD_BASE            2
#define  GRIB2_LTYPE_CLOUD_TOP             3
#define  GRIB2_LTYPE_ISOTHERM0             4
#define  GRIB2_LTYPE_TROPOPAUSE            7
#define  GRIB2_LTYPE_TOA                   8
#define  GRIB2_LTYPE_SEA_BOTTOM            9
#define  GRIB2_LTYPE_ATMOSPHERE           10
#define  GRIB2_LTYPE_ISOBARIC            100
#define  GRIB2_LTYPE_MEANSEA             101
#define  GRIB2_LTYPE_ALTITUDE            102
#define  GRIB2_LTYPE_HEIGHT              103
#define  GRIB2_LTYPE_SIGMA               104
#define  GRIB2_LTYPE_HYBRID              105
#define  GRIB2_LTYPE_LANDDEPTH           106
#define  GRIB2_LTYPE_ISENTROPIC          107
#define  GRIB2_LTYPE_SNOW                114
#define  GRIB2_LTYPE_REFERENCE           150
#define  GRIB2_LTYPE_SEADEPTH            160  // Depth Below Sea Level
#define  GRIB2_LTYPE_LAKE_BOTTOM         162  // Lake or River Bottom
#define  GRIB2_LTYPE_SEDIMENT_BOTTOM     163  // Bottom Of Sediment Layer
#define  GRIB2_LTYPE_SEDIMENT_BOTTOM_TA  164  // Bottom Of Thermally Active Sediment Layer
#define  GRIB2_LTYPE_SEDIMENT_BOTTOM_TW  165  // Bottom Of Sediment Layer Penetrated By Thermal Wave
#define  GRIB2_LTYPE_MIX_LAYER           166  // Mixing Layer

// GRIB2 Data representation type (Grid Type)
#define  GRIB2_GTYPE_LATLON                0  // Latitude/longitude (or equidistant cylindrical, or Plate Carree)
#define  GRIB2_GTYPE_LATLON_ROT            1  // Rotated Latitude/longitude
#define  GRIB2_GTYPE_LATLON_STR            2  // Stretched Latitude/longitude
#define  GRIB2_GTYPE_LATLON_ROTSTR         3  // Stretched and Rotated Latitude/longitude
#define  GRIB2_GTYPE_STERE                20  // Polar stereographic projection
#define  GRIB2_GTYPE_LCC                  30  // Lambert conformal
#define  GRIB2_GTYPE_GAUSSIAN             40  // Gaussian latitude/longitude
#define  GRIB2_GTYPE_GAUSSIAN_ROT         41  // Rotated Gaussian latitude/longitude
#define  GRIB2_GTYPE_GAUSSIAN_STR         42  // Stretched Gaussian latitude/longitude
#define  GRIB2_GTYPE_GAUSSIAN_ROTSTR      43  // Stretched and rotated Gaussian latitude/longitude
#define  GRIB2_GTYPE_SPECTRAL             50  // Spherical harmonic coefficients
#define  GRIB2_GTYPE_GME                 100  // Triangular grid based on an icosahedron (GME)
#define  GRIB2_GTYPE_UNSTRUCTURED        101  // General Unstructured Grid

const char *gribapiLibraryVersionString(void);
void gribContainersNew(stream_t *streamptr);
void gribContainersDelete(stream_t *streamptr);

#ifdef HAVE_LIBGRIB_API

#ifdef ECCODES_VERSION
#if ECCODES_VERSION >= 23000
#define HAVE_GRIBAPI_FLOAT_INTERFACE 1
#endif
#endif

static inline int have_gribapi_float_interface(void)
{
#ifdef HAVE_GRIBAPI_FLOAT_INTERFACE
  return 1;
#else
  return 0;
#endif
}

static inline int my_grib_set_double(grib_handle* h, const char* key, double val)
{
  if (CDI_gribapi_debug)
    fprintf(stderr, "grib_set_double(\tgrib_handle* h, \"%s\", %f)\n", key, val);

  int retVal = grib_set_double(h, key, val);
  if (retVal != 0)
    fprintf(stderr, "!!! failed call to grib_set_double(\tgrib_handle* h, \"%s\", %f) !!!\n", key, val);
  return retVal;
}

static inline int my_grib_set_long(grib_handle* h, const char* key, long val)
{
  if (CDI_gribapi_debug)
    fprintf(stderr, "grib_set_long(  \tgrib_handle* h, \"%s\", %ld)\n", key, val);

  int retVal = grib_set_long(h, key, val);
  if (retVal != 0)
    fprintf(stderr, "!!! failed call to grib_set_long(  \tgrib_handle* h, \"%s\", %ld) !!!\n", key, val);
  return retVal;
}

static inline int my_grib_set_string(grib_handle* h, const char* key, const char* val, size_t* length)
{
  if (CDI_gribapi_debug)
    fprintf(stderr, "grib_set_string(\tgrib_handle* h, \"%s\", \"%s\")\n", key, val);

  int ret_val = grib_set_string(h, key, val, length);
  if (ret_val != 0)
    fprintf(stderr, "!!! grib_set_string(\tgrib_handle* h, \"%s\", \"%s\") !!!\n", key, val);
  return ret_val;
}
#endif

typedef struct {
  bool init;
  void *gribHandle;
}
gribContainer_t;

// clang-format on

#endif /* GRIBAPI_H */
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
