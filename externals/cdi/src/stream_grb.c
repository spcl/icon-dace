#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "dmemory.h"
#include "cdi.h"
#include "cdi_int.h"
#include "stream_cgribex.h"
#include "stream_grb.h"
#include "stream_gribapi.h"
#include "gribapi.h"
#include "file.h"
#include "cgribex.h" /* gribZip gribGetZip gribGinfo */

int cdiDebugExt = 0;  //  Debug level for the KNMI extensions
#ifdef HIRLAM_EXTENSIONS
// *** RELATED to GRIB only ***
int cdiGribUseTimeRangeIndicator = 0;  // normaly cdo looks in grib for attribute called "stepType"
                                       // but NWP models such as Harmonie 37h1.2, use "timeRangeIndicator"
                                       // where:  0: for instanteneous fields; 4: for accumulated fields
#endif                                 // HIRLAM_EXTENSIONS

double
zaxis_units_to_centimeter(int zaxisID)
{
  char units[CDI_MAX_NAME];
  int length = CDI_MAX_NAME;
  cdiInqKeyString(zaxisID, CDI_GLOBAL, CDI_KEY_UNITS, units, &length);

  double sf = 100.0;  // default: meter
  // clang-format off
  if (units[1] == 'm' && !units[2])
    {
      if      (units[0] == 'm') sf =   0.1;
      else if (units[0] == 'c') sf =   1.0;
      else if (units[0] == 'd') sf =  10.0;
    }
  // clang-format on

  return sf;
}

double
zaxis_units_to_meter(int zaxisID)
{
  char units[CDI_MAX_NAME];
  int length = CDI_MAX_NAME;
  cdiInqKeyString(zaxisID, CDI_GLOBAL, CDI_KEY_UNITS, units, &length);

  double sf = 1.0;  // default: meter
  // clang-format off
  if (units[1] == 'm' && !units[2])
    {
      if      (units[0] == 'm') sf /= 1000.0;
      else if (units[0] == 'c') sf /= 100.0;
      else if (units[0] == 'd') sf /= 10.0;
      else if (units[0] == 'k') sf *= 1000.0;
    }
  // clang-format on

  return sf;
}

bool
zaxis_units_is_Pa(int zaxisID)
{
  char units[CDI_MAX_NAME] = { 0 };
  int length = CDI_MAX_NAME;
  cdiInqKeyString(zaxisID, CDI_GLOBAL, CDI_KEY_UNITS, units, &length);

  return (units[0] == 'P' && units[1] == 'a');
}

void
ensureBufferSize(size_t requiredSize, size_t *curSize, void **buffer)
{
  if (*curSize < requiredSize)
    {
      *curSize = requiredSize;
      *buffer = Realloc(*buffer, *curSize);
    }
}

int
grbDecompress(size_t recsize, size_t *buffersize, void **gribbuffer)
{
  int comptype = CDI_COMPRESS_NONE;

  size_t unzipsize;
  if (gribGetZip(recsize, (unsigned char *) *gribbuffer, &unzipsize) > 0)
    {
      comptype = CDI_COMPRESS_SZIP;
      unzipsize += 100;  // need 0 to 1 bytes for rounding of bds
      ensureBufferSize(unzipsize, buffersize, gribbuffer);
    }

  return comptype;
}

// Regarding operation to change parameter identification:
// change if cdiGribChangeParameterID.active
struct cdiGribParamChange cdiGribChangeParameterID;

// Used only for CDO module Selmulti
void
streamGrbChangeParameterIdentification(int code, int ltype, int lev)
{
  // NOTE this is a "PROXY" function for gribapiChangeParameterIdentification();
  // This just sets the globals. There are probably better solutions to this.
  // The parameter change is done by function  gribapiChangeParameterIdentification() in stream_gribapi.c
  // Setting this control variable to true will cause calling fnc. gribapiChangeParameterIdentification later.
  // After grib attributes have been changed this variable goes to false.
  cdiGribChangeParameterID.active = true;
  cdiGribChangeParameterID.code = code;
  cdiGribChangeParameterID.ltype = ltype;
  cdiGribChangeParameterID.lev = lev;
}

struct cdiGribScanModeChange cdiGribDataScanningMode;

void
streamGrbDefDataScanningMode(int scanmode)
{
  cdiGribDataScanningMode.active = true;
  cdiGribDataScanningMode.value = scanmode;
}

enum
{
  MapKey = 0,
  MapValue = 1
};

static const int grib1ltypeMap[][2] = {
  // clang-format off
  {  GRIB1_LTYPE_SURFACE,             ZAXIS_SURFACE            },
  {  GRIB1_LTYPE_CLOUD_BASE,          ZAXIS_CLOUD_BASE         },
  {  GRIB1_LTYPE_CLOUD_TOP,           ZAXIS_CLOUD_TOP          },
  {  GRIB1_LTYPE_ISOTHERM0,           ZAXIS_ISOTHERM_ZERO      },
  {  GRIB1_LTYPE_TROPOPAUSE,          ZAXIS_TROPOPAUSE         },
  {  GRIB1_LTYPE_TOA,                 ZAXIS_TOA                },
  {  GRIB1_LTYPE_SEA_BOTTOM,          ZAXIS_SEA_BOTTOM         },
  {  GRIB1_LTYPE_ATMOSPHERE,          ZAXIS_ATMOSPHERE         },
  {  GRIB1_LTYPE_ISOBARIC,            ZAXIS_PRESSURE           },
  {  GRIB1_LTYPE_99,                  ZAXIS_PRESSURE           },
  {  GRIB1_LTYPE_ISOBARIC_PA,         ZAXIS_PRESSURE           },
  {  GRIB1_LTYPE_MEANSEA,             ZAXIS_MEANSEA            },
  {  GRIB1_LTYPE_ALTITUDE,            ZAXIS_ALTITUDE           },
  {  GRIB1_LTYPE_HEIGHT,              ZAXIS_HEIGHT             },
  {  GRIB1_LTYPE_SIGMA,               ZAXIS_SIGMA              },
  {  GRIB1_LTYPE_SIGMA_LAYER,         ZAXIS_SIGMA              },
  {  GRIB1_LTYPE_HYBRID,              ZAXIS_HYBRID             },
  {  GRIB1_LTYPE_HYBRID_LAYER,        ZAXIS_HYBRID             },
  {  GRIB1_LTYPE_LANDDEPTH,           ZAXIS_DEPTH_BELOW_LAND   },
  {  GRIB1_LTYPE_LANDDEPTH_LAYER,     ZAXIS_DEPTH_BELOW_LAND   },
  {  GRIB1_LTYPE_ISENTROPIC,          ZAXIS_ISENTROPIC         },
  {  GRIB1_LTYPE_SEADEPTH,            ZAXIS_DEPTH_BELOW_SEA    },
  {  GRIB1_LTYPE_LAKE_BOTTOM,         ZAXIS_LAKE_BOTTOM        },
  {  GRIB1_LTYPE_SEDIMENT_BOTTOM,     ZAXIS_SEDIMENT_BOTTOM    },
  {  GRIB1_LTYPE_SEDIMENT_BOTTOM_TA,  ZAXIS_SEDIMENT_BOTTOM_TA },
  {  GRIB1_LTYPE_SEDIMENT_BOTTOM_TW,  ZAXIS_SEDIMENT_BOTTOM_TW },
  {  GRIB1_LTYPE_MIX_LAYER,           ZAXIS_MIX_LAYER          },
  // clang-format on
};

static const size_t grib1ltypeMapSize = sizeof(grib1ltypeMap) / (2 * sizeof(int));

static const int grib2ltypeMap[][2] = {
  // clang-format off
  {  GRIB2_LTYPE_SURFACE,             ZAXIS_SURFACE            },
  {  GRIB2_LTYPE_CLOUD_BASE,          ZAXIS_CLOUD_BASE         },
  {  GRIB2_LTYPE_CLOUD_TOP,           ZAXIS_CLOUD_TOP          },
  {  GRIB2_LTYPE_ISOTHERM0,           ZAXIS_ISOTHERM_ZERO      },
  {  GRIB2_LTYPE_TROPOPAUSE,          ZAXIS_TROPOPAUSE         },
  {  GRIB2_LTYPE_TOA,                 ZAXIS_TOA                },
  {  GRIB2_LTYPE_SEA_BOTTOM,          ZAXIS_SEA_BOTTOM         },
  {  GRIB2_LTYPE_ATMOSPHERE,          ZAXIS_ATMOSPHERE         },
  {  GRIB2_LTYPE_ISOBARIC,            ZAXIS_PRESSURE           },
  {  GRIB2_LTYPE_MEANSEA,             ZAXIS_MEANSEA            },
  {  GRIB2_LTYPE_ALTITUDE,            ZAXIS_ALTITUDE           },
  {  GRIB2_LTYPE_HEIGHT,              ZAXIS_HEIGHT             },
  {  GRIB2_LTYPE_SIGMA,               ZAXIS_SIGMA              },
  {  GRIB2_LTYPE_HYBRID,              ZAXIS_HYBRID             },
  {  GRIB2_LTYPE_HYBRID,              ZAXIS_HYBRID_HALF        },
  {  GRIB2_LTYPE_LANDDEPTH,           ZAXIS_DEPTH_BELOW_LAND   },
  {  GRIB2_LTYPE_ISENTROPIC,          ZAXIS_ISENTROPIC         },
  {  GRIB2_LTYPE_SEADEPTH,            ZAXIS_DEPTH_BELOW_SEA    },
  {  GRIB2_LTYPE_LAKE_BOTTOM,         ZAXIS_LAKE_BOTTOM        },
  {  GRIB2_LTYPE_SEDIMENT_BOTTOM,     ZAXIS_SEDIMENT_BOTTOM    },
  {  GRIB2_LTYPE_SEDIMENT_BOTTOM_TA,  ZAXIS_SEDIMENT_BOTTOM_TA },
  {  GRIB2_LTYPE_SEDIMENT_BOTTOM_TW,  ZAXIS_SEDIMENT_BOTTOM_TW },
  {  GRIB2_LTYPE_MIX_LAYER,           ZAXIS_MIX_LAYER          },
  {  GRIB2_LTYPE_SNOW,                ZAXIS_SNOW               },
  {  GRIB2_LTYPE_REFERENCE,           ZAXIS_REFERENCE          },
  // clang-format on
};

static const size_t grib2ltypeMapSize = sizeof(grib2ltypeMap) / (2 * sizeof(int));

static int
getInt2IntMap(int searchKey, bool keyValue, size_t mapSize, const int gribltypeMap[][2], int defaultValue)
{
  for (size_t i = 0; i < mapSize; ++i)
    if (gribltypeMap[i][keyValue] == searchKey) return gribltypeMap[i][!keyValue];

  return defaultValue;
}

int
grib1ltypeToZaxisType(int grib_ltype)
{
  return getInt2IntMap(grib_ltype, MapKey, grib1ltypeMapSize, grib1ltypeMap, ZAXIS_GENERIC);
}

int
zaxisTypeToGrib1ltype(int zaxistype)
{
  return getInt2IntMap(zaxistype, MapValue, grib1ltypeMapSize, grib1ltypeMap, -1);
}

int
grib2ltypeToZaxisType(int grib_ltype)
{
  return getInt2IntMap(grib_ltype, MapKey, grib2ltypeMapSize, grib2ltypeMap, ZAXIS_GENERIC);
}

int
zaxisTypeToGrib2ltype(int zaxistype)
{
  return getInt2IntMap(zaxistype, MapValue, grib2ltypeMapSize, grib2ltypeMap, -1);
}

int
grbBitsPerValue(int datatype)
{
  int bitsPerValue = 16;

  if (datatype == CDI_DATATYPE_CPX32 || datatype == CDI_DATATYPE_CPX64) Error("CDI/GRIB library does not support complex numbers!");

  if (datatype != CDI_UNDEFID)
    {
      if (datatype > 0 && datatype <= 32)
        bitsPerValue = datatype;
      else if (datatype == CDI_DATATYPE_FLT64)
        bitsPerValue = 24;
      else
        bitsPerValue = 16;
    }

  return bitsPerValue;
}

/*
int grbInqRecord(stream_t * streamptr, int *varID, int *levelID)
{
  int status;

  status = cgribexInqRecord(streamptr, varID, levelID);

  return (status);
}
*/

void
grbDefRecord(stream_t *streamptr)
{
  UNUSED(streamptr);
}

static int
grbScanTimestep1(stream_t *streamptr)
{
  int status = CDI_EUFTYPE;

#ifdef HAVE_LIBCGRIBEX
  int filetype = streamptr->filetype;

  if (filetype == CDI_FILETYPE_GRB && !CDI_gribapi_grib1)
    status = cgribexScanTimestep1(streamptr);
  else
#endif
#ifdef HAVE_LIBGRIB_API
    status = gribapiScanTimestep1(streamptr);
#else
  Error("GRIB_API support unavailable!");
#endif

  return status;
}

static int
grbScanTimestep2(stream_t *streamptr)
{
  int status = CDI_EUFTYPE;

#ifdef HAVE_LIBCGRIBEX
  int filetype = streamptr->filetype;

  if (filetype == CDI_FILETYPE_GRB && !CDI_gribapi_grib1)
    status = cgribexScanTimestep2(streamptr);
  else
#endif
#ifdef HAVE_LIBGRIB_API
    status = gribapiScanTimestep2(streamptr);
#else
  Error("GRIB_API support unavailable!");
#endif

  return status;
}

static int
grbScanTimestep(stream_t *streamptr)
{
  int status = CDI_EUFTYPE;

#ifdef HAVE_LIBCGRIBEX
  int filetype = streamptr->filetype;

  if (filetype == CDI_FILETYPE_GRB && !CDI_gribapi_grib1)
    status = cgribexScanTimestep(streamptr);
  else
#endif
#ifdef HAVE_LIBGRIB_API
    status = gribapiScanTimestep(streamptr);
#else
  Error("GRIB_API support unavailable!");
#endif

  return status;
}

#ifdef HAVE_LIBGRIB
int
grbInqContents(stream_t *streamptr)
{
  streamptr->curTsID = 0;

  int status = grbScanTimestep1(streamptr);
  if (status == 0 && streamptr->ntsteps == -1) status = grbScanTimestep2(streamptr);

  int fileID = streamptr->fileID;
  fileSetPos(fileID, 0, SEEK_SET);

  return status;
}

int
fdbInqContents(stream_t *streamptr)
{
  streamptr->curTsID = 0;
#ifdef HAVE_LIBGRIB_API
  return fdbScanTimesteps(streamptr);
#else
  return -1;
#endif
}
#endif

int
fdbInqTimestep(stream_t *streamptr, int tsID)
{
  if (tsID == 0 && streamptr->rtsteps == 0) Error("Call to cdiInqContents missing!");

  if (CDI_Debug) Message("tsid = %d rtsteps = %d", tsID, streamptr->rtsteps);

  int nrecs;
  if (tsID >= streamptr->ntsteps && streamptr->ntsteps != CDI_UNDEFID)
    {
      nrecs = 0;
    }
  else
    {
      streamptr->curTsID = tsID;
      nrecs = streamptr->tsteps[tsID].nrecs;
    }

  return nrecs;
}

int
grbInqTimestep(stream_t *streamptr, int tsID)
{
  if (tsID == 0 && streamptr->rtsteps == 0) Error("Call to cdiInqContents missing!");

  if (CDI_Debug) Message("tsid = %d rtsteps = %d", tsID, streamptr->rtsteps);

  int ntsteps = CDI_UNDEFID;
  while ((tsID + 1) > streamptr->rtsteps && ntsteps == CDI_UNDEFID)
    {
      ntsteps = grbScanTimestep(streamptr);
      if (ntsteps == CDI_EUFSTRUCT)
        {
          streamptr->ntsteps = streamptr->rtsteps;
          break;
        }
    }

  int nrecs;
  if (tsID >= streamptr->ntsteps && streamptr->ntsteps != CDI_UNDEFID)
    {
      nrecs = 0;
    }
  else
    {
      streamptr->curTsID = tsID;
      nrecs = streamptr->tsteps[tsID].nrecs;
    }

  return nrecs;
}

// used in CDO!!!
void
streamInqGRIBinfo(int streamID, int *intnum, float *fltnum, off_t *bignum)
{
  stream_t *streamptr = stream_to_pointer(streamID);

  int filetype = streamptr->filetype;

  if (filetype == CDI_FILETYPE_GRB)
    {
      int tsID = streamptr->curTsID;
      int vrecID = streamptr->tsteps[tsID].curRecID;
      int recID = streamptr->tsteps[tsID].recIDs[vrecID];
      off_t recpos = streamptr->tsteps[tsID].records[recID].position;
      int zip = streamptr->tsteps[tsID].records[recID].zip;

      void *gribbuffer = streamptr->record->buffer;
      size_t gribbuffersize = streamptr->record->buffersize;

      if (zip > 0)
        Error("Compressed GRIB records unsupported!");
      else
        grib_info_for_grads(recpos, (long) gribbuffersize, (unsigned char *) gribbuffer, intnum, fltnum, bignum);
    }
}

int
grbGetGridtype(int *gridID, size_t gridsize, bool *gridIsRotated, bool *gridIsCurvilinear)
{
  int gridtype = gridInqType(*gridID);

  if (gridtype == GRID_GENERIC)
    {
      int xsize = (int) gridInqXsize(*gridID);
      int ysize = (int) gridInqYsize(*gridID);

      if ((ysize == 32 || ysize == 48 || ysize == 64 || ysize == 96 || ysize == 160 || ysize == 192 || ysize == 240 || ysize == 320
           || ysize == 384 || ysize == 480 || ysize == 768)
          && (xsize == 2 * ysize || xsize == 1))
        {
          gridtype = GRID_GAUSSIAN;
        }
      else if (gridsize == 1)
        {
          gridtype = GRID_LONLAT;
        }
      else if (gridInqXvals(*gridID, NULL) && gridInqYvals(*gridID, NULL))
        {
          gridtype = GRID_LONLAT;
        }
    }
  else if (gridtype == GRID_CURVILINEAR)
    {
      int projID = gridInqProj(*gridID);
      if (projID != CDI_UNDEFID && gridInqType(projID) == GRID_PROJECTION)
        {
          *gridID = projID;
          gridtype = GRID_PROJECTION;
        }
      else
        {
          static bool lwarning = true;
          if (lwarning && gridsize > 1)
            {
              lwarning = false;
              Warning("Curvilinear grid is unsupported in GRIB format! Created wrong Grid Description Section!");
            }
          *gridIsCurvilinear = true;
          gridtype = GRID_LONLAT;
        }
    }

  if (gridtype == GRID_PROJECTION)
    {
      int projtype = gridInqProjType(*gridID);
      if (projtype == CDI_PROJ_RLL)
        {
          gridtype = GRID_LONLAT;
          *gridIsRotated = true;
        }
      else if (projtype == CDI_PROJ_LCC)
        {
          gridtype = CDI_PROJ_LCC;
        }
      else if (projtype == CDI_PROJ_STERE)
        {
          gridtype = CDI_PROJ_STERE;
        }
      else if (projtype == CDI_PROJ_HEALPIX)
        {
          gridtype = CDI_PROJ_HEALPIX;
        }
    }

  return gridtype;
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
