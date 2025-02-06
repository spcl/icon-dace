#include <string.h>
#include <math.h>
#include <float.h>

#include "dmemory.h"

#include "cdi.h"
#include "cdi_cksum.h"
#include "cdi_int.h"
#include "cdi_uuid.h"
#include "resource_handle.h"
#include "resource_unpack.h"
#include "namespace.h"
#include "serialize.h"
#include "zaxis.h"

#define LevelUp 1
#define LevelDown 2

// clang-format off
static const struct
{
  unsigned char positive;   // 1: up;  2: down
  const char *name;
  const char *longname;
  const char *stdname;
  const char *units;
}
ZaxistypeEntry[] = {
  { /*  0 */ 0, "sfc",               "surface",                     "",               ""},
  { /*  1 */ 0, "lev",               "generic",                     "",               ""},
  { /*  2 */ 2, "lev",               "hybrid",                      "",               "level"},
  { /*  3 */ 2, "lev",               "hybrid_half",                 "",               "level"},
  { /*  4 */ 2, "plev",              "pressure",                    "air_pressure",   "Pa"},
  { /*  5 */ 1, "height",            "height",                      "height",         "m"},
  { /*  6 */ 2, "depth",             "depth_below_sea",             "depth",          "m"},
  { /*  7 */ 2, "depth",             "depth_below_land",            "",               "cm"},
  { /*  8 */ 0, "lev",               "isentropic",                  "",               "K"},
  { /*  9 */ 0, "lev",               "trajectory",                  "",               ""},
  { /* 10 */ 1, "alt",               "height above mean sea level", "altitude",       "m"},
  { /* 11 */ 0, "lev",               "sigma",                       "",               "level"},
  { /* 12 */ 0, "lev",               "meansea",                     "",               "level"},
  { /* 13 */ 0, "toa",               "top_of_atmosphere",           "",               ""},
  { /* 14 */ 0, "seabottom",         "sea_bottom",                  "",               ""},
  { /* 15 */ 0, "atmosphere",        "atmosphere",                  "",               ""},
  { /* 16 */ 0, "cloudbase",         "cloud_base",                  "",               ""},
  { /* 17 */ 0, "cloudtop",          "cloud_top",                   "",               ""},
  { /* 18 */ 0, "isotherm0",         "isotherm_zero",               "",               ""},
  { /* 19 */ 0, "snow",              "snow",                        "",               ""},
  { /* 20 */ 0, "lakebottom",        "lake_bottom",                 "",               ""},
  { /* 21 */ 0, "sedimentbottom",    "sediment_bottom",             "",               ""},
  { /* 22 */ 0, "sedimentbottomta",  "sediment_bottom_ta",          "",               ""},
  { /* 23 */ 0, "sedimentbottomtw",  "sediment_bottom_tw",          "",               ""},
  { /* 24 */ 0, "mixlayer",          "mix_layer",                   "",               ""},
  { /* 25 */ 0, "height",            "generalized_height",          "height",         ""},
  { /* 26 */ 0, "character",         "area_type",                   "",               ""},
  { /* 27 */ 0, "tropopause",        "tropopause",                  "",               ""},
};
// clang-format on

enum
{
  CDI_NumZaxistype = sizeof(ZaxistypeEntry) / sizeof(ZaxistypeEntry[0]),
};

static int zaxisCompareP(zaxis_t *z1, zaxis_t *z2);
static void zaxisDestroyP(void *zaxisptr);
static void zaxisPrintP(void *zaxisptr, FILE *fp);
static int zaxisGetPackSize(void *zaxisptr, void *context);
static void zaxisPack(void *zaxisptr, void *buffer, int size, int *pos, void *context);
static int zaxisTxCode(void *zaxisptr);

static const resOps zaxisOps
    = { (int (*)(void *, void *)) zaxisCompareP, zaxisDestroyP, zaxisPrintP, zaxisGetPackSize, zaxisPack, zaxisTxCode };

const resOps *
getZaxisOps(void)
{
  return &zaxisOps;
}

void
zaxisGetTypeDescription(int zaxisType, int *outPositive, const char **outName, const char **outLongName, const char **outStdName,
                        const char **outUnit)
{
  if (zaxisType < 0 || zaxisType >= CDI_NumZaxistype)
    {
      if (outPositive) *outPositive = 0;
      if (outName) *outName = NULL;
      if (outLongName) *outLongName = NULL;
      if (outStdName) *outStdName = NULL;
      if (outUnit) *outUnit = NULL;
    }
  else
    {
      if (outPositive) *outPositive = ZaxistypeEntry[zaxisType].positive;
      if (outName) *outName = ZaxistypeEntry[zaxisType].name;
      if (outLongName && zaxisType != ZAXIS_GENERIC) *outLongName = ZaxistypeEntry[zaxisType].longname;
      if (outStdName) *outStdName = ZaxistypeEntry[zaxisType].stdname;
      if (outUnit) *outUnit = ZaxistypeEntry[zaxisType].units;
    }
}

zaxis_t *
zaxis_to_pointer(int id)
{
  return (zaxis_t *) reshGetVal(id, &zaxisOps);
}

static void
zaxis_init(zaxis_t *zaxisptr)
{
  zaxisptr->self = CDI_UNDEFID;
  zaxisptr->vals = NULL;
#ifndef USE_MPI
  zaxisptr->cvals = NULL;
  zaxisptr->clength = 0;
#endif
  zaxisptr->ubounds = NULL;
  zaxisptr->lbounds = NULL;
  zaxisptr->weights = NULL;
  zaxisptr->type = CDI_UNDEFID;
  zaxisptr->positive = 0;
  zaxisptr->scalar = 0;
  zaxisptr->direction = 0;
  zaxisptr->size = 0;
  zaxisptr->vctsize = 0;
  zaxisptr->vct = NULL;

  cdiInitKeys(&zaxisptr->keys);
  zaxisptr->atts.nalloc = MAX_ATTRIBUTES;
  zaxisptr->atts.nelems = 0;

  cdiDefVarKeyInt(&zaxisptr->keys, CDI_KEY_DATATYPE, CDI_DATATYPE_FLT64);
}

static zaxis_t *
zaxisNewEntry(int id)
{
  zaxis_t *zaxisptr = (zaxis_t *) Malloc(sizeof(zaxis_t));
  zaxis_init(zaxisptr);

  if (id == CDI_UNDEFID)
    zaxisptr->self = reshPut(zaxisptr, &zaxisOps);
  else
    {
      zaxisptr->self = id;
      reshReplace(id, zaxisptr, &zaxisOps);
    }

  return zaxisptr;
}

static void
zaxisInit(void)
{
  static bool zaxisInitialized = false;
  if (zaxisInitialized) return;
  zaxisInitialized = true;
}

static void
zaxis_copy(zaxis_t *zaxisptr2, zaxis_t *zaxisptr1)
{
  int zaxisID2 = zaxisptr2->self;
  memcpy(zaxisptr2, zaxisptr1, sizeof(zaxis_t));
  zaxisptr2->self = zaxisID2;
  cdiInitKeys(&zaxisptr2->keys);
  cdiCopyVarKeys(&zaxisptr1->keys, &zaxisptr2->keys);
}

unsigned
cdiZaxisCount(void)
{
  return reshCountType(&zaxisOps);
}

static int
zaxisCreate_(int zaxistype, int size, int id)
{
  zaxis_t *zaxisptr = zaxisNewEntry(id);

  xassert(size >= 0);
  zaxisptr->type = zaxistype;
  zaxisptr->size = size;

  int zaxisID = zaxisptr->self;

  if (zaxistype >= 0 && zaxistype < CDI_NumZaxistype)
    {
      cdiDefKeyString(zaxisID, CDI_GLOBAL, CDI_KEY_NAME, ZaxistypeEntry[zaxistype].name);
      if (zaxistype != ZAXIS_GENERIC) zaxisDefLongname(zaxisID, ZaxistypeEntry[zaxistype].longname);
      cdiDefKeyString(zaxisID, CDI_GLOBAL, CDI_KEY_UNITS, ZaxistypeEntry[zaxistype].units);

      const char *stdname = ZaxistypeEntry[zaxistype].stdname;
      if (*stdname) cdiDefVarKeyBytes(&zaxisptr->keys, CDI_KEY_STDNAME, (const unsigned char *) stdname, (int) strlen(stdname) + 1);

      zaxisptr->positive = ZaxistypeEntry[zaxistype].positive;
    }
  else
    {
      Error("Internal problem! zaxistype=%d out of range (min=0/max=%d)!", zaxistype, CDI_NumZaxistype - 1);
    }

  return zaxisID;
}

/*
@Function  zaxisCreate
@Title     Create a vertical Z-axis

@Prototype int zaxisCreate(int zaxistype, int size)
@Parameter
    @Item  zaxistype  The type of the Z-axis, one of the set of predefined CDI Z-axis types.
                      The valid CDI Z-axis types are @func{ZAXIS_GENERIC}, @func{ZAXIS_SURFACE},
                      @func{ZAXIS_HYBRID}, @func{ZAXIS_SIGMA}, @func{ZAXIS_PRESSURE}, @func{ZAXIS_HEIGHT},
                      @func{ZAXIS_ISENTROPIC}, @func{ZAXIS_ALTITUDE}, @func{ZAXIS_MEANSEA}, @func{ZAXIS_TOA},
                      @func{ZAXIS_SEA_BOTTOM}, @func{ZAXIS_ATMOSPHERE}, @func{ZAXIS_CLOUD_BASE},
                      @func{ZAXIS_CLOUD_TOP}, @func{ZAXIS_ISOTHERM_ZERO}, @func{ZAXIS_SNOW},
                      @func{ZAXIS_LAKE_BOTTOM}, @func{ZAXIS_SEDIMENT_BOTTOM}, @func{ZAXIS_SEDIMENT_BOTTOM_TA},
                      @func{ZAXIS_SEDIMENT_BOTTOM_TW}, @func{ZAXIS_MIX_LAYER},
                      @func{ZAXIS_DEPTH_BELOW_SEA} and @func{ZAXIS_DEPTH_BELOW_LAND}.
    @Item  size       Number of levels.

@Description
The function @func{zaxisCreate} creates a vertical Z-axis.

@Result
@func{zaxisCreate} returns an identifier to the Z-axis.

@Example
Here is an example using @func{zaxisCreate} to create a pressure level Z-axis:

@Source
#include "cdi.h"
   ...
#define  nlev    5
   ...
double levs[nlev] = {101300, 92500, 85000, 50000, 20000};
int zaxisID;
   ...
zaxisID = zaxisCreate(ZAXIS_PRESSURE, nlev);
zaxisDefLevels(zaxisID, levs);
   ...
@EndSource
@EndFunction
*/
int
zaxisCreate(int zaxistype, int size)
{
  if (CDI_Debug) Message("zaxistype: %d size: %d ", zaxistype, size);

  xassert(size);
  zaxisInit();

  return zaxisCreate_(zaxistype, size, CDI_UNDEFID);
}

static void
zaxisDestroyKernel(zaxis_t *zaxisptr)
{
  xassert(zaxisptr);

  if (zaxisptr->vals) Free(zaxisptr->vals);
#ifndef USE_MPI
  if (zaxisptr->cvals)
    {
      for (int i = 0; i < zaxisptr->size; i++) Free(zaxisptr->cvals[i]);
      Free(zaxisptr->cvals);
    }
#endif
  if (zaxisptr->lbounds) Free(zaxisptr->lbounds);
  if (zaxisptr->ubounds) Free(zaxisptr->ubounds);
  if (zaxisptr->weights) Free(zaxisptr->weights);
  if (zaxisptr->vct) Free(zaxisptr->vct);

  int zaxisID = zaxisptr->self;
  cdiDeleteKeys(zaxisID, CDI_GLOBAL);
  cdiDeleteAtts(zaxisID, CDI_GLOBAL);

  Free(zaxisptr);
}

/*
@Function  zaxisDestroy
@Title     Destroy a vertical Z-axis

@Prototype void zaxisDestroy(int zaxisID)
@Parameter
    @Item  zaxisID  Z-axis ID, from a previous call to @fref{zaxisCreate}.

@EndFunction
*/
void
zaxisDestroy(int zaxisID)
{
  zaxis_t *zaxisptr = zaxis_to_pointer(zaxisID);
  zaxisDestroyKernel(zaxisptr);
  reshRemove(zaxisID, &zaxisOps);
}

static void
zaxisDestroyP(void *zaxisptr)
{
  zaxisDestroyKernel((zaxis_t *) zaxisptr);
}

const char *
zaxisNamePtr(int zaxistype)
{
  const char *name = (zaxistype >= 0 && zaxistype < CDI_NumZaxistype) ? ZaxistypeEntry[zaxistype].longname
                                                                      : ZaxistypeEntry[ZAXIS_GENERIC].longname;
  return name;
}

void
zaxisName(int zaxistype, char *zaxisname)
{
  strcpy(zaxisname, zaxisNamePtr(zaxistype));
}

// obsolete function
void
zaxisDefLtype(int zaxisID, int ltype)
{
  static bool printInfo = true;
  if (printInfo) printInfo = cdiObsoleteInfo(__func__, "cdiDefKeyInt");

  (void) cdiDefKeyInt(zaxisID, CDI_GLOBAL, CDI_KEY_TYPEOFFIRSTFIXEDSURFACE, ltype);
}

/*
@Function  zaxisDefName
@Title     Define the name of a Z-axis

@Prototype void zaxisDefName(int zaxisID, const char *name)
@Parameter
    @Item  zaxisID  Z-axis ID, from a previous call to @fref{zaxisCreate}.
    @Item  name     Name of the Z-axis.

@Description
The function @func{zaxisDefName} defines the name of a Z-axis.

@EndFunction
*/
void
zaxisDefName(int zaxisID, const char *name)
{
  (void) cdiDefKeyString(zaxisID, CDI_GLOBAL, CDI_KEY_NAME, name);
}

/*
@Function  zaxisDefLongname
@Title     Define the longname of a Z-axis

@Prototype void zaxisDefLongname(int zaxisID, const char *longname)
@Parameter
    @Item  zaxisID  Z-axis ID, from a previous call to @fref{zaxisCreate}.
    @Item  longname Longname of the Z-axis.

@Description
The function @func{zaxisDefLongname} defines the longname of a Z-axis.

@EndFunction
*/
void
zaxisDefLongname(int zaxisID, const char *longname)
{
  (void) cdiDefKeyString(zaxisID, CDI_GLOBAL, CDI_KEY_LONGNAME, longname);
}

/*
@Function  zaxisDefUnits
@Title     Define the units of a Z-axis

@Prototype void zaxisDefUnits(int zaxisID, const char *units)
@Parameter
    @Item  zaxisID  Z-axis ID, from a previous call to @fref{zaxisCreate}.
    @Item  units    Units of the Z-axis.

@Description
The function @func{zaxisDefUnits} defines the units of a Z-axis.

@EndFunction
*/
void
zaxisDefUnits(int zaxisID, const char *units)
{
  (void) cdiDefKeyString(zaxisID, CDI_GLOBAL, CDI_KEY_UNITS, units);
}

/*
@Function  zaxisInqName
@Title     Get the name of a Z-axis

@Prototype void zaxisInqName(int zaxisID, char *name)
@Parameter
    @Item  zaxisID  Z-axis ID, from a previous call to @fref{zaxisCreate} or @fref{vlistInqVarZaxis}.
    @Item  name     Name of the Z-axis. The caller must allocate space for the
                    returned string. The maximum possible length, in characters, of
                    the string is given by the predefined constant @func{CDI_MAX_NAME}.

@Description
The function @func{zaxisInqName} returns the name of a Z-axis.

@Result
@func{zaxisInqName} returns the name of the Z-axis to the parameter name.

@EndFunction
*/
void
zaxisInqName(int zaxisID, char *name)
{
  int length = CDI_MAX_NAME;
  (void) cdiInqKeyString(zaxisID, CDI_GLOBAL, CDI_KEY_NAME, name, &length);
}

const char *
zaxisInqNamePtr(int zaxisID)
{
  zaxis_t *zaxisptr = zaxis_to_pointer(zaxisID);
  return cdiInqVarKeyString(&zaxisptr->keys, CDI_KEY_NAME);
}

/*
@Function  zaxisInqLongname
@Title     Get the longname of a Z-axis

@Prototype void zaxisInqLongname(int zaxisID, char *longname)
@Parameter
    @Item  zaxisID  Z-axis ID, from a previous call to @fref{zaxisCreate} or @fref{vlistInqVarZaxis}.
    @Item  longname Longname of the Z-axis. The caller must allocate space for the
                    returned string. The maximum possible length, in characters, of
                    the string is given by the predefined constant @func{CDI_MAX_NAME}.

@Description
The function @func{zaxisInqLongname} returns the longname of a Z-axis.

@Result
@func{zaxisInqLongname} returns the longname of the Z-axis to the parameter longname.

@EndFunction
*/
void
zaxisInqLongname(int zaxisID, char *longname)
{
  int length = CDI_MAX_NAME;
  (void) cdiInqKeyString(zaxisID, CDI_GLOBAL, CDI_KEY_LONGNAME, longname, &length);
}

/*
@Function  zaxisInqUnits
@Title     Get the units of a Z-axis

@Prototype void zaxisInqUnits(int zaxisID, char *units)
@Parameter
    @Item  zaxisID  Z-axis ID, from a previous call to @fref{zaxisCreate} or @fref{vlistInqVarZaxis}.
    @Item  units    Units of the Z-axis. The caller must allocate space for the
                    returned string. The maximum possible length, in characters, of
                    the string is given by the predefined constant @func{CDI_MAX_NAME}.

@Description
The function @func{zaxisInqUnits} returns the units of a Z-axis.

@Result
@func{zaxisInqUnits} returns the units of the Z-axis to the parameter units.

@EndFunction
*/
void
zaxisInqUnits(int zaxisID, char *units)
{
  int length = CDI_MAX_NAME;
  (void) cdiInqKeyString(zaxisID, CDI_GLOBAL, CDI_KEY_UNITS, units, &length);
}

void
zaxisInqStdname(int zaxisID, char *stdname)
{
  int length = CDI_MAX_NAME;
  (void) cdiInqKeyString(zaxisID, CDI_GLOBAL, CDI_KEY_STDNAME, stdname, &length);
}

void
zaxisDefDatatype(int zaxisID, int datatype)
{
  cdiDefKeyInt(zaxisID, CDI_GLOBAL, CDI_KEY_DATATYPE, datatype);
}

int
zaxisInqDatatype(int zaxisID)
{
  int datatype = 0;
  cdiInqKeyInt(zaxisID, CDI_GLOBAL, CDI_KEY_DATATYPE, &datatype);
  return datatype;
}

void
zaxisDefPositive(int zaxisID, int positive)
{
  zaxis_t *zaxisptr = zaxis_to_pointer(zaxisID);

  if (zaxisptr->positive != (unsigned) positive)
    {
      zaxisptr->positive = (unsigned) positive;
      reshSetStatus(zaxisID, &zaxisOps, RESH_DESYNC_IN_USE);
    }
}

int
zaxisInqPositive(int zaxisID)
{
  zaxis_t *zaxisptr = zaxis_to_pointer(zaxisID);
  return (int) zaxisptr->positive;
}

void
zaxisDefScalar(int zaxisID)
{
  zaxis_t *zaxisptr = zaxis_to_pointer(zaxisID);

  zaxisptr->scalar = 1;
  reshSetStatus(zaxisID, &zaxisOps, RESH_DESYNC_IN_USE);
}

int
zaxisInqScalar(int zaxisID)
{
  zaxis_t *zaxisptr = zaxis_to_pointer(zaxisID);
  return zaxisptr->scalar;
}

/*
@Function  zaxisDefLevels
@Title     Define the levels of a Z-axis

@Prototype void zaxisDefLevels(int zaxisID, const double *levels)
@Parameter
    @Item  zaxisID  Z-axis ID, from a previous call to @fref{zaxisCreate}.
    @Item  levels   All levels of the Z-axis.

@Description
The function @func{zaxisDefLevels} defines the levels of a Z-axis.

@EndFunction
*/
void
zaxisDefLevels(int zaxisID, const double *levels)
{
  if (levels)
    {
      zaxis_t *zaxisptr = zaxis_to_pointer(zaxisID);
      const size_t size = (size_t) zaxisptr->size;
      xassert(size);

      if (zaxisptr->vals == NULL && size) zaxisptr->vals = (double *) Malloc(size * sizeof(double));

      double *vals = zaxisptr->vals;

      for (size_t ilev = 0; ilev < size; ++ilev) vals[ilev] = levels[ilev];

      reshSetStatus(zaxisID, &zaxisOps, RESH_DESYNC_IN_USE);
    }
}

void
zaxisDefCvals(int zaxisID, const char **cvals, int clen)
{
#ifndef USE_MPI
  if (cvals && clen)
    {
      zaxis_t *zaxisptr = zaxis_to_pointer(zaxisID);
      const size_t size = zaxisptr->size;
      xassert(size);

      zaxisptr->clength = clen;
      if (size) zaxisptr->cvals = (char **) Malloc(size * sizeof(char *));

      for (size_t ilev = 0; ilev < size; ++ilev)
        {
          zaxisptr->cvals[ilev] = (char *) Malloc(clen * sizeof(char));
          memcpy(zaxisptr->cvals[ilev], cvals[ilev], clen * sizeof(char));
        }
      reshSetStatus(zaxisID, &zaxisOps, RESH_DESYNC_IN_USE);
    }
#else
  Error("This function was disabled!");
#endif
}

/*
@Function  zaxisDefLevel
@Title     Define one level of a Z-axis

@Prototype void zaxisDefLevel(int zaxisID, int levelID, double level)
@Parameter
    @Item  zaxisID  Z-axis ID, from a previous call to @fref{zaxisCreate}.
    @Item  levelID  Level identifier.
    @Item  level    Level.

@Description
The function @func{zaxisDefLevel} defines one level of a Z-axis.

@EndFunction
*/
void
zaxisDefLevel(int zaxisID, int levelID, double level)
{
  zaxis_t *zaxisptr = zaxis_to_pointer(zaxisID);
  int size = zaxisptr->size;
  xassert(size);
  xassert(levelID >= 0 && levelID < size);

  if (zaxisptr->vals == NULL && size) zaxisptr->vals = (double *) Malloc((size_t) size * sizeof(double));

  if (levelID >= 0 && levelID < size) zaxisptr->vals[levelID] = level;

  reshSetStatus(zaxisID, &zaxisOps, RESH_DESYNC_IN_USE);
}

void
zaxisDefNlevRef(int zaxisID, int nlev)
{
  cdiDefKeyInt(zaxisID, CDI_GLOBAL, CDI_KEY_NLEV, nlev);
}

int
zaxisInqNlevRef(int zaxisID)
{
  int nlev = 0;
  cdiInqKeyInt(zaxisID, CDI_GLOBAL, CDI_KEY_NLEV, &nlev);
  return nlev;
}

/*
@Function  zaxisDefNumber
@Title     Define the reference number for a generalized Z-axis

@Prototype void zaxisDefNumber(int zaxisID, int number)
@Parameter
    @Item  zaxisID     Z-axis ID, from a previous call to @fref{zaxisCreate}.
    @Item  number      Reference number for a generalized Z-axis.

@Description
The function @func{zaxisDefNumber} defines the reference number for a generalized Z-axis.

@EndFunction
*/
void
zaxisDefNumber(int zaxisID, int number)
{
  cdiDefKeyInt(zaxisID, CDI_GLOBAL, CDI_KEY_NUMBEROFVGRIDUSED, number);
}

/*
@Function  zaxisInqNumber
@Title     Get the reference number to a generalized Z-axis

@Prototype int zaxisInqNumber(int zaxisID)
@Parameter
    @Item  zaxisID  Z-axis ID, from a previous call to @fref{zaxisCreate} or @fref{vlistInqVarZaxis}.

@Description
The function @func{zaxisInqNumber} returns the reference number to a generalized Z-axis.

@Result
@func{zaxisInqNumber} returns the reference number to a generalized Z-axis.
@EndFunction
*/
int
zaxisInqNumber(int zaxisID)
{
  int referenceNumber = 0;
  cdiInqKeyInt(zaxisID, CDI_GLOBAL, CDI_KEY_NUMBEROFVGRIDUSED, &referenceNumber);
  return referenceNumber;
}

/*
@Function  zaxisDefUUID
@Title     Define the UUID for a genralized Z-axis

@Prototype void zaxisDefUUID(int zaxisID, const char *uuid)
@Parameter
    @Item  zaxisID     Z-axis ID, from a previous call to @fref{zaxisCreate}.
    @Item  uuid        UUID for a generalized Z-axis.

@Description
The function @func{zaxisDefUUID} defines the UUID for a generalized  Z-axis.

@EndFunction
*/
void
zaxisDefUUID(int zaxisID, const unsigned char uuid[CDI_UUID_SIZE])
{
  cdiDefKeyBytes(zaxisID, CDI_GLOBAL, CDI_KEY_UUID, uuid, CDI_UUID_SIZE);

  reshSetStatus(zaxisID, &zaxisOps, RESH_DESYNC_IN_USE);
}

/*
@Function  zaxisInqUUID
@Title     Get the uuid to a generalized Z-axis

@Prototype void zaxisInqUUID(int zaxisID, char *uuid)
@Parameter
    @Item  zaxisID  Z-axis ID, from a previous call to @fref{zaxisCreate} or @fref{vlistInqVarZaxis}.
    @Item uuid A user supplied buffer of at least 16 bytes.

@Description
The function @func{zaxisInqUUID} returns the UUID to a generalized Z-axis.

@Result
@func{zaxisInqUUID} returns the UUID to a generalized Z-axis to the parameter uuid.
@EndFunction
*/
void
zaxisInqUUID(int zaxisID, unsigned char uuid[CDI_UUID_SIZE])
{
  memset(uuid, 0, CDI_UUID_SIZE);
  int length = CDI_UUID_SIZE;
  cdiInqKeyBytes(zaxisID, CDI_GLOBAL, CDI_KEY_UUID, uuid, &length);
}

/*
@Function  zaxisInqLevel
@Title     Get one level of a Z-axis

@Prototype double zaxisInqLevel(int zaxisID, int levelID)
@Parameter
    @Item  zaxisID  Z-axis ID, from a previous call to @fref{zaxisCreate} or @fref{vlistInqVarZaxis}.
    @Item  levelID  Level index (range: 0 to nlevel-1).

@Description
The function @func{zaxisInqLevel} returns one level of a Z-axis.

@Result
@func{zaxisInqLevel} returns the level of a Z-axis.
@EndFunction
*/
double
zaxisInqLevel(int zaxisID, int levelID)
{
  double level = 0;
  zaxis_t *zaxisptr = zaxis_to_pointer(zaxisID);
  if (zaxisptr->vals && levelID >= 0 && levelID < zaxisptr->size) level = zaxisptr->vals[levelID];

  return level;
}

double
zaxisInqLbound(int zaxisID, int levelID)
{
  double level = 0;
  zaxis_t *zaxisptr = zaxis_to_pointer(zaxisID);
  if (zaxisptr->lbounds && levelID >= 0 && levelID < zaxisptr->size) level = zaxisptr->lbounds[levelID];

  return level;
}

double
zaxisInqUbound(int zaxisID, int levelID)
{
  double level = 0;
  zaxis_t *zaxisptr = zaxis_to_pointer(zaxisID);
  if (zaxisptr->ubounds && levelID >= 0 && levelID < zaxisptr->size) level = zaxisptr->ubounds[levelID];

  return level;
}

const double *
zaxisInqLevelsPtr(int zaxisID)
{
  zaxis_t *zaxisptr = zaxis_to_pointer(zaxisID);
  return zaxisptr->vals;
}

#ifndef USE_MPI
char **
zaxisInqCValsPtr(int zaxisID)
{
  zaxis_t *zaxisptr = zaxis_to_pointer(zaxisID);
  return zaxisptr->cvals;
}
#endif

/*
@Function  zaxisInqLevels
@Title     Get all levels of a Z-axis

@Prototype void zaxisInqLevels(int zaxisID, double *levels)
@Parameter
    @Item  zaxisID  Z-axis ID, from a previous call to @fref{zaxisCreate} or @fref{vlistInqVarZaxis}.
    @Item  levels   Pointer to the location into which the levels are read.
                    The caller must allocate space for the returned values.

@Description
The function @func{zaxisInqLevels} returns all levels of a Z-axis.

@Result
@func{zaxisInqLevels} saves all levels to the parameter @func{levels}.
@EndFunction
*/
int
zaxisInqLevels(int zaxisID, double *levels)
{
  zaxis_t *zaxisptr = zaxis_to_pointer(zaxisID);

  int size = 0;
  if (zaxisptr->vals)
    {
      size = zaxisptr->size;

      if (levels)
        for (int i = 0; i < size; i++) levels[i] = zaxisptr->vals[i];
    }

  return size;
}

int
zaxisInqCLen(int zaxisID)
{
  int clen = 0;
#ifndef USE_MPI
  zaxis_t *zaxisptr = zaxis_to_pointer(zaxisID);
  if (zaxisptr->cvals && zaxisptr->clength) clen = zaxisptr->clength;
#endif

  return clen;
}

int
zaxisInqCVals(int zaxisID, char ***clevels)
{
  int size = 0;
#ifndef USE_MPI
  zaxis_t *zaxisptr = zaxis_to_pointer(zaxisID);
  if (zaxisptr->cvals)
    {
      size = zaxisptr->size;
      const size_t clen = zaxisptr->clength;
      if (size && clen)
        {
          (*clevels) = (char **) Malloc(size * sizeof(char *));
          for (int i = 0; i < size; i++)
            {
              (*clevels)[i] = (char *) Malloc(clen * sizeof(char));
              memcpy((*clevels)[i], zaxisptr->cvals[i], clen * sizeof(char));
            }
        }
    }
#endif

  return size;
}

int
zaxisInqLbounds(int zaxisID, double *lbounds)
{
  int size = 0;
  zaxis_t *zaxisptr = zaxis_to_pointer(zaxisID);
  if (zaxisptr->lbounds)
    {
      size = zaxisptr->size;

      if (lbounds)
        for (int i = 0; i < size; i++) lbounds[i] = zaxisptr->lbounds[i];
    }

  return size;
}

int
zaxisInqUbounds(int zaxisID, double *ubounds)
{
  int size = 0;
  zaxis_t *zaxisptr = zaxis_to_pointer(zaxisID);
  if (zaxisptr->ubounds)
    {
      size = zaxisptr->size;

      if (ubounds)
        for (int i = 0; i < size; i++) ubounds[i] = zaxisptr->ubounds[i];
    }

  return size;
}

int
zaxisInqWeights(int zaxisID, double *weights)
{
  int size = 0;
  zaxis_t *zaxisptr = zaxis_to_pointer(zaxisID);
  if (zaxisptr->weights)
    {
      size = zaxisptr->size;

      if (weights)
        for (int i = 0; i < size; i++) weights[i] = zaxisptr->weights[i];
    }

  return size;
}

int
zaxisInqLevelID(int zaxisID, double level)
{
  int levelID = CDI_UNDEFID;
  zaxis_t *zaxisptr = zaxis_to_pointer(zaxisID);
  if (zaxisptr->vals)
    {
      int size = zaxisptr->size;
      for (int i = 0; i < size; i++)
        if (fabs(level - zaxisptr->vals[i]) < DBL_EPSILON)
          {
            levelID = i;
            break;
          }
    }

  return levelID;
}

/*
@Function  zaxisInqType
@Title     Get the type of a Z-axis

@Prototype int zaxisInqType(int zaxisID)
@Parameter
    @Item  zaxisID  Z-axis ID, from a previous call to @fref{zaxisCreate} or @fref{vlistInqVarZaxis}.

@Description
The function @func{zaxisInqType} returns the type of a Z-axis.

@Result
@func{zaxisInqType} returns the type of the Z-axis,
one of the set of predefined CDI Z-axis types.
The valid CDI Z-axis types are @func{ZAXIS_GENERIC}, @func{ZAXIS_SURFACE},
@func{ZAXIS_HYBRID}, @func{ZAXIS_SIGMA}, @func{ZAXIS_PRESSURE}, @func{ZAXIS_HEIGHT},
@func{ZAXIS_ISENTROPIC}, @func{ZAXIS_ALTITUDE}, @func{ZAXIS_MEANSEA}, @func{ZAXIS_TOA},
@func{ZAXIS_SEA_BOTTOM}, @func{ZAXIS_ATMOSPHERE}, @func{ZAXIS_CLOUD_BASE},
@func{ZAXIS_CLOUD_TOP}, @func{ZAXIS_ISOTHERM_ZERO}, @func{ZAXIS_SNOW},
@func{ZAXIS_LAKE_BOTTOM}, @func{ZAXIS_SEDIMENT_BOTTOM}, @func{ZAXIS_SEDIMENT_BOTTOM_TA},
@func{ZAXIS_SEDIMENT_BOTTOM_TW}, @func{ZAXIS_MIX_LAYER},
@func{ZAXIS_DEPTH_BELOW_SEA} and @func{ZAXIS_DEPTH_BELOW_LAND}.

@EndFunction
*/
int
zaxisInqType(int zaxisID)
{
  zaxis_t *zaxisptr = zaxis_to_pointer(zaxisID);
  return zaxisptr->type;
}

/*
@Function  zaxisInqSize
@Title     Get the size of a Z-axis

@Prototype int zaxisInqSize(int zaxisID)
@Parameter
    @Item  zaxisID  Z-axis ID, from a previous call to @fref{zaxisCreate} or @fref{vlistInqVarZaxis}.

@Description
The function @func{zaxisInqSize} returns the size of a Z-axis.

@Result
@func{zaxisInqSize} returns the number of levels of a Z-axis.

@EndFunction
*/
int
zaxisInqSize(int zaxisID)
{
  zaxis_t *zaxisptr = zaxis_to_pointer(zaxisID);
  return zaxisptr->size;
}

void
cdiCheckZaxis(int zaxisID)
{
  zaxis_t *zaxisptr = zaxis_to_pointer(zaxisID);

  if (zaxisInqType(zaxisID) == ZAXIS_GENERIC && zaxisptr->vals)
    {
      int size = zaxisptr->size;
      if (size > 1)
        {
          /* check direction */
          if (!zaxisptr->direction)
            {
              int ups = 0, downs = 0;
              for (int i = 1; i < size; i++)
                {
                  ups += (zaxisptr->vals[i] > zaxisptr->vals[i - 1]);
                  downs += (zaxisptr->vals[i] < zaxisptr->vals[i - 1]);
                }
              if (ups == size - 1)
                {
                  zaxisptr->direction = LevelUp;
                }
              else if (downs == size - 1)
                {
                  zaxisptr->direction = LevelDown;
                }
              else /* !zaxisptr->direction */
                {
                  Warning("Direction undefined for zaxisID %d", zaxisID);
                }
            }
        }
    }
}

void
zaxisDefVct(int zaxisID, int size, const double *vct)
{
  zaxis_t *zaxisptr = zaxis_to_pointer(zaxisID);

  if (zaxisptr->vct == 0 || zaxisptr->vctsize != size)
    {
      zaxisptr->vctsize = size;
      zaxisptr->vct = (double *) Realloc(zaxisptr->vct, (size_t) size * sizeof(double));
    }

  if (vct) memcpy(zaxisptr->vct, vct, (size_t) size * sizeof(double));
  reshSetStatus(zaxisID, &zaxisOps, RESH_DESYNC_IN_USE);
}

void
zaxisInqVct(int zaxisID, double *vct)
{
  zaxis_t *zaxisptr = zaxis_to_pointer(zaxisID);
  memcpy(vct, zaxisptr->vct, (size_t) zaxisptr->vctsize * sizeof(double));
}

int
zaxisInqVctSize(int zaxisID)
{
  zaxis_t *zaxisptr = zaxis_to_pointer(zaxisID);
  return zaxisptr->vctsize;
}

const double *
zaxisInqVctPtr(int zaxisID)
{
  zaxis_t *zaxisptr = zaxis_to_pointer(zaxisID);
  return zaxisptr->vct;
}

void
zaxisDefLbounds(int zaxisID, const double *lbounds)
{
  zaxis_t *zaxisptr = zaxis_to_pointer(zaxisID);

  const size_t size = (size_t) zaxisptr->size;

  if (CDI_Debug)
    if (zaxisptr->lbounds) Warning("Lower bounds already defined for zaxisID = %d", zaxisID);

  if (zaxisptr->lbounds == NULL) zaxisptr->lbounds = (double *) Malloc(size * sizeof(double));

  if (lbounds) memcpy(zaxisptr->lbounds, lbounds, size * sizeof(double));
  reshSetStatus(zaxisID, &zaxisOps, RESH_DESYNC_IN_USE);
}

void
zaxisDefUbounds(int zaxisID, const double *ubounds)
{
  zaxis_t *zaxisptr = zaxis_to_pointer(zaxisID);

  const size_t size = (size_t) zaxisptr->size;

  if (CDI_Debug)
    if (zaxisptr->ubounds) Warning("Upper bounds already defined for zaxisID = %d", zaxisID);

  if (zaxisptr->ubounds == NULL) zaxisptr->ubounds = (double *) Malloc(size * sizeof(double));

  if (ubounds) memcpy(zaxisptr->ubounds, ubounds, size * sizeof(double));
  reshSetStatus(zaxisID, &zaxisOps, RESH_DESYNC_IN_USE);
}

void
zaxisDefWeights(int zaxisID, const double *weights)
{
  zaxis_t *zaxisptr = zaxis_to_pointer(zaxisID);

  const size_t size = (size_t) zaxisptr->size;

  if (CDI_Debug)
    if (zaxisptr->weights != NULL) Warning("Weights already defined for zaxisID = %d", zaxisID);

  if (zaxisptr->weights == NULL) zaxisptr->weights = (double *) Malloc(size * sizeof(double));

  memcpy(zaxisptr->weights, weights, size * sizeof(double));
  reshSetStatus(zaxisID, &zaxisOps, RESH_DESYNC_IN_USE);
}

void
zaxisChangeType(int zaxisID, int zaxistype)
{
  zaxis_t *zaxisptr = zaxis_to_pointer(zaxisID);
  zaxisptr->type = zaxistype;
}

void
zaxisResize(int zaxisID, int size)
{
  zaxis_t *zaxisptr = zaxis_to_pointer(zaxisID);

  xassert(size >= 0);

  zaxisptr->size = size;

  if (zaxisptr->vals) zaxisptr->vals = (double *) Realloc(zaxisptr->vals, (size_t) size * sizeof(double));
}

static inline void
zaxisCopyKeyStr(zaxis_t *zaxisptr1, zaxis_t *zaxisptr2, int key)
{
  cdi_key_t *keyp = find_key(&zaxisptr1->keys, key);
  if (keyp && keyp->type == KEY_BYTES)
    cdiDefVarKeyBytes(&zaxisptr2->keys, key, (const unsigned char *) keyp->v.s, (int) keyp->length);
}

int
zaxisDuplicate(int zaxisID)
{
  zaxis_t *zaxisptr = zaxis_to_pointer(zaxisID);

  int zaxistype = zaxisInqType(zaxisID);
  int zaxissize = zaxisInqSize(zaxisID);

  int zaxisIDnew = zaxisCreate(zaxistype, zaxissize);
  zaxis_t *zaxisptrnew = zaxis_to_pointer(zaxisIDnew);

  zaxis_copy(zaxisptrnew, zaxisptr);

  zaxisCopyKeyStr(zaxisptr, zaxisptrnew, CDI_KEY_NAME);
  zaxisCopyKeyStr(zaxisptr, zaxisptrnew, CDI_KEY_LONGNAME);
  zaxisCopyKeyStr(zaxisptr, zaxisptrnew, CDI_KEY_UNITS);

  if (zaxisptr->vals)
    {
      const size_t size = (size_t) zaxissize;
      zaxisptrnew->vals = (double *) Malloc(size * sizeof(double));
      memcpy(zaxisptrnew->vals, zaxisptr->vals, size * sizeof(double));
    }

  if (zaxisptr->lbounds)
    {
      const size_t size = (size_t) zaxissize;
      zaxisptrnew->lbounds = (double *) Malloc(size * sizeof(double));
      memcpy(zaxisptrnew->lbounds, zaxisptr->lbounds, size * sizeof(double));
    }

  if (zaxisptr->ubounds)
    {
      const size_t size = (size_t) zaxissize;
      zaxisptrnew->ubounds = (double *) Malloc(size * sizeof(double));
      memcpy(zaxisptrnew->ubounds, zaxisptr->ubounds, size * sizeof(double));
    }

  if (zaxisptr->vct)
    {
      const size_t size = (size_t) zaxisptr->vctsize;
      if (size)
        {
          zaxisptrnew->vctsize = (int) size;
          zaxisptrnew->vct = (double *) Malloc(size * sizeof(double));
          memcpy(zaxisptrnew->vct, zaxisptr->vct, size * sizeof(double));
        }
    }

  zaxisptrnew->atts.nelems = 0;
  cdiCopyAtts(zaxisID, CDI_GLOBAL, zaxisIDnew, CDI_GLOBAL);

  return zaxisIDnew;
}

static void
zaxisPrintKernel(zaxis_t *zaxisptr, FILE *fp)
{
  xassert(zaxisptr);

  int zaxisID = zaxisptr->self;
  int datatype = CDI_UNDEFID;
  cdiInqKeyInt(zaxisID, CDI_GLOBAL, CDI_KEY_DATATYPE, &datatype);

  int type = zaxisptr->type;
  int nlevels = zaxisptr->size;

  int dig = (datatype == CDI_DATATYPE_FLT64) ? 15 : 7;

  fprintf(fp, "zaxistype = %s\n", zaxisNamePtr(type));
  fprintf(fp, "size      = %d\n", nlevels);
  if (nlevels == 1)
    {
      const bool zscalar = (bool) zaxisptr->scalar;
      if (zscalar) fprintf(fp, "scalar    = true\n");
    }

  const char *string = cdiInqVarKeyString(&zaxisptr->keys, CDI_KEY_NAME);
  if (string[0]) fprintf(fp, "name      = %s\n", string);
  string = cdiInqVarKeyString(&zaxisptr->keys, CDI_KEY_LONGNAME);
  if (string[0]) fprintf(fp, "longname  = %s\n", string);
  string = cdiInqVarKeyString(&zaxisptr->keys, CDI_KEY_UNITS);
  if (string[0]) fprintf(fp, "units     = %s\n", string);

  if (zaxisptr->vals)
    {
      int nbyte0 = fprintf(fp, "levels    = ");
      int nbyte = nbyte0;
      for (int levelID = 0; levelID < nlevels; levelID++)
        {
          if (nbyte > 80)
            {
              fprintf(fp, "\n");
              fprintf(fp, "%*s", nbyte0, "");
              nbyte = nbyte0;
            }
          nbyte += fprintf(fp, "%.*g ", dig, zaxisptr->vals[levelID]);
        }
      fprintf(fp, "\n");
    }

  if (zaxisptr->lbounds && zaxisptr->ubounds)
    {
      int nbyte0 = fprintf(fp, "lbounds   = ");
      int nbyte = nbyte0;
      for (int levelID = 0; levelID < nlevels; levelID++)
        {
          if (nbyte > 80)
            {
              fprintf(fp, "\n");
              fprintf(fp, "%*s", nbyte0, "");
              nbyte = nbyte0;
            }
          nbyte += fprintf(fp, "%.*g ", dig, zaxisptr->lbounds[levelID]);
        }
      fprintf(fp, "\n");

      nbyte0 = fprintf(fp, "ubounds   = ");
      nbyte = nbyte0;
      for (int levelID = 0; levelID < nlevels; levelID++)
        {
          if (nbyte > 80)
            {
              fprintf(fp, "\n");
              fprintf(fp, "%*s", nbyte0, "");
              nbyte = nbyte0;
            }
          nbyte += fprintf(fp, "%.*g ", dig, zaxisptr->ubounds[levelID]);
        }
      fprintf(fp, "\n");
    }

  if (type == ZAXIS_HYBRID || type == ZAXIS_HYBRID_HALF)
    {
      int vctsize = zaxisptr->vctsize;
      const double *vct = zaxisptr->vct;
      fprintf(fp, "vctsize   = %d\n", vctsize);
      if (vctsize)
        {
          int nbyte0 = fprintf(fp, "vct       = ");
          int nbyte = nbyte0;
          for (int i = 0; i < vctsize; i++)
            {
              if (nbyte > 70 || i == vctsize / 2)
                {
                  fprintf(fp, "\n%*s", nbyte0, "");
                  nbyte = nbyte0;
                }
              nbyte += fprintf(fp, "%.15g ", vct[i]);
            }
          fprintf(fp, "\n");
        }
    }
}

static void
zaxisPrintP(void *voidptr, FILE *fp)
{
  zaxis_t *zaxisptr = (zaxis_t *) voidptr;

  xassert(zaxisptr);

  zaxisPrintKernel(zaxisptr, fp);
}

static int
zaxisCompareP(zaxis_t *z1, zaxis_t *z2)
{
  enum
  {
    differ = 1
  };
  int diff = 0;
  xassert(z1 && z2);

  diff |= (z1->type != z2->type)
          | (cdiInqVarKeyInt(&z1->keys, CDI_KEY_TYPEOFFIRSTFIXEDSURFACE)
             != cdiInqVarKeyInt(&z2->keys, CDI_KEY_TYPEOFFIRSTFIXEDSURFACE))
          | (cdiInqVarKeyInt(&z1->keys, CDI_KEY_DATATYPE) != cdiInqVarKeyInt(&z2->keys, CDI_KEY_DATATYPE))
          | (z1->direction != z2->direction) | (z1->size != z2->size) | (z1->vctsize != z2->vctsize)
          | (z1->positive != z2->positive);

  if (diff) return differ;

  int size = z1->size;
  int anyPresent = 0;
  int present = (z1->vals != NULL);
  diff |= (present ^ (z2->vals != NULL));
  anyPresent |= present;
  if (!diff && present)
    {
      const double *p = z1->vals, *q = z2->vals;
      for (int i = 0; i < size; i++) diff |= IS_NOT_EQUAL(p[i], q[i]);
    }

  present = (z1->lbounds != NULL);
  diff |= (present ^ (z2->lbounds != NULL));
  anyPresent |= present;
  if (!diff && present)
    {
      const double *p = z1->lbounds, *q = z2->lbounds;
      for (int i = 0; i < size; i++) diff |= IS_NOT_EQUAL(p[i], q[i]);
    }

  present = (z1->ubounds != NULL);
  diff |= (present ^ (z2->ubounds != NULL));
  anyPresent |= present;
  if (!diff && present)
    {
      const double *p = z1->ubounds, *q = z2->ubounds;
      for (int i = 0; i < size; ++i) diff |= IS_NOT_EQUAL(p[i], q[i]);
    }

  present = (z1->weights != NULL);
  diff |= (present ^ (z2->weights != NULL));
  anyPresent |= present;
  if (!diff && present)
    {
      const double *p = z1->weights, *q = z2->weights;
      for (int i = 0; i < size; ++i) diff |= IS_NOT_EQUAL(p[i], q[i]);
    }

  present = (z1->vct != NULL);
  diff |= (present ^ (z2->vct != NULL));
  if (!diff && present)
    {
      int vctsize = z1->vctsize;
      xassert(vctsize);
      const double *p = z1->vct, *q = z2->vct;
      for (int i = 0; i < vctsize; ++i) diff |= IS_NOT_EQUAL(p[i], q[i]);
    }

  if (anyPresent) xassert(size);

  diff |= strcmp(cdiInqVarKeyString(&z1->keys, CDI_KEY_NAME), cdiInqVarKeyString(&z2->keys, CDI_KEY_NAME))
          | strcmp(cdiInqVarKeyString(&z1->keys, CDI_KEY_LONGNAME), cdiInqVarKeyString(&z2->keys, CDI_KEY_LONGNAME))
          | strcmp(cdiInqVarKeyString(&z1->keys, CDI_KEY_STDNAME), cdiInqVarKeyString(&z2->keys, CDI_KEY_STDNAME))
          | strcmp(cdiInqVarKeyString(&z1->keys, CDI_KEY_UNITS), cdiInqVarKeyString(&z2->keys, CDI_KEY_UNITS));

  return diff != 0;
}

static int
zaxisTxCode(void *zaxisPtr)
{
  (void) zaxisPtr;
  return ZAXIS;
}

enum
{
  ZAXIS_PACK_INT_SELF,
  ZAXIS_PACK_INT_TYPE,
  ZAXIS_PACK_INT_SIZE,
  ZAXIS_PACK_INT_DIRECTION,
  ZAXIS_PACK_INT_VCTSIZE,
  ZAXIS_PACK_INT_MEMBERMASK,
  zaxisNint
};

enum
{
  vals = 1 << 0,
  lbounds = 1 << 1,
  ubounds = 1 << 2,
  weights = 1 << 3,
  vct = 1 << 4,
};

static int
zaxisGetMemberMask(zaxis_t *zaxisP)
{
  int memberMask = 0;

  if (zaxisP->vals) memberMask |= vals;
  if (zaxisP->lbounds) memberMask |= lbounds;
  if (zaxisP->ubounds) memberMask |= ubounds;
  if (zaxisP->weights) memberMask |= weights;
  if (zaxisP->vct) memberMask |= vct;

  return memberMask;
}

static int
zaxisGetPackSize(void *voidP, void *context)
{
  zaxis_t *zaxisP = (zaxis_t *) voidP;
  int packBufferSize = serializeGetSize(zaxisNint, CDI_DATATYPE_INT, context) + serializeGetSize(1, CDI_DATATYPE_UINT32, context);

  if (zaxisP->vals || zaxisP->lbounds || zaxisP->ubounds || zaxisP->weights) xassert(zaxisP->size);

  if (zaxisP->vals)
    packBufferSize
        += serializeGetSize(zaxisP->size, CDI_DATATYPE_FLT64, context) + serializeGetSize(1, CDI_DATATYPE_UINT32, context);

  if (zaxisP->lbounds)
    packBufferSize
        += serializeGetSize(zaxisP->size, CDI_DATATYPE_FLT64, context) + serializeGetSize(1, CDI_DATATYPE_UINT32, context);

  if (zaxisP->ubounds)
    packBufferSize
        += serializeGetSize(zaxisP->size, CDI_DATATYPE_FLT64, context) + serializeGetSize(1, CDI_DATATYPE_UINT32, context);

  if (zaxisP->weights)
    packBufferSize
        += serializeGetSize(zaxisP->size, CDI_DATATYPE_FLT64, context) + serializeGetSize(1, CDI_DATATYPE_UINT32, context);

  if (zaxisP->vct)
    {
      xassert(zaxisP->vctsize);
      packBufferSize
          += serializeGetSize(zaxisP->vctsize, CDI_DATATYPE_FLT64, context) + serializeGetSize(1, CDI_DATATYPE_UINT32, context);
    }

  packBufferSize += serializeKeysGetPackSize(&zaxisP->keys, context);

  packBufferSize += serializeGetSize(1, CDI_DATATYPE_UINT, context);

  return packBufferSize;
}

int
zaxisUnpack(char *unpackBuffer, int unpackBufferSize, int *unpackBufferPos, int originNamespace, void *context, int force_id)
{
  int intBuffer[zaxisNint], memberMask;
  uint32_t d;

  serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos, intBuffer, zaxisNint, CDI_DATATYPE_INT, context);
  serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos, &d, 1, CDI_DATATYPE_UINT32, context);

  xassert(cdiCheckSum(CDI_DATATYPE_INT, zaxisNint, intBuffer) == d);

  zaxisInit();

  zaxis_t *zaxisP = zaxisNewEntry(force_id ? namespaceAdaptKey(intBuffer[ZAXIS_PACK_INT_SELF], originNamespace) : CDI_UNDEFID);

  zaxisP->type = intBuffer[ZAXIS_PACK_INT_TYPE];
  zaxisP->size = intBuffer[ZAXIS_PACK_INT_SIZE];
  zaxisP->direction = intBuffer[ZAXIS_PACK_INT_DIRECTION];
  zaxisP->vctsize = intBuffer[ZAXIS_PACK_INT_VCTSIZE];
  memberMask = intBuffer[ZAXIS_PACK_INT_MEMBERMASK];

  if (memberMask & vals)
    {
      int size = zaxisP->size;
      xassert(size >= 0);

      zaxisP->vals = (double *) Malloc((size_t) size * sizeof(double));
      serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos, zaxisP->vals, size, CDI_DATATYPE_FLT64, context);
      serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos, &d, 1, CDI_DATATYPE_UINT32, context);
      xassert(cdiCheckSum(CDI_DATATYPE_FLT, size, zaxisP->vals) == d);
    }

  if (memberMask & lbounds)
    {
      int size = zaxisP->size;
      xassert(size >= 0);

      zaxisP->lbounds = (double *) Malloc((size_t) size * sizeof(double));
      serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos, zaxisP->lbounds, size, CDI_DATATYPE_FLT64, context);
      serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos, &d, 1, CDI_DATATYPE_UINT32, context);
      xassert(cdiCheckSum(CDI_DATATYPE_FLT, size, zaxisP->lbounds) == d);
    }

  if (memberMask & ubounds)
    {
      int size = zaxisP->size;
      xassert(size >= 0);

      zaxisP->ubounds = (double *) Malloc((size_t) size * sizeof(double));
      serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos, zaxisP->ubounds, size, CDI_DATATYPE_FLT64, context);
      serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos, &d, 1, CDI_DATATYPE_UINT32, context);
      xassert(cdiCheckSum(CDI_DATATYPE_FLT, size, zaxisP->ubounds) == d);
    }

  if (memberMask & weights)
    {
      int size = zaxisP->size;
      xassert(size >= 0);

      zaxisP->weights = (double *) Malloc((size_t) size * sizeof(double));
      serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos, zaxisP->weights, size, CDI_DATATYPE_FLT64, context);
      serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos, &d, 1, CDI_DATATYPE_UINT32, context);
      xassert(cdiCheckSum(CDI_DATATYPE_FLT, size, zaxisP->weights) == d);
    }

  if (memberMask & vct)
    {
      int size = zaxisP->vctsize;
      xassert(size >= 0);

      zaxisP->vct = (double *) Malloc((size_t) size * sizeof(double));
      serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos, zaxisP->vct, size, CDI_DATATYPE_FLT64, context);
      serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos, &d, 1, CDI_DATATYPE_UINT32, context);
      xassert(cdiCheckSum(CDI_DATATYPE_FLT64, size, zaxisP->vct) == d);
    }

  serializeKeysUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos, &zaxisP->keys, context);

  serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos, &zaxisP->positive, 1, CDI_DATATYPE_UINT, context);

  reshSetStatus(zaxisP->self, &zaxisOps, reshGetStatus(zaxisP->self, &zaxisOps) & ~RESH_SYNC_BIT);
  return zaxisP->self;
}

static void
zaxisPack(void *voidP, void *packBuffer, int packBufferSize, int *packBufferPos, void *context)
{
  zaxis_t *zaxisP = (zaxis_t *) voidP;
  int intBuffer[zaxisNint];
  int memberMask;
  uint32_t d;

  intBuffer[ZAXIS_PACK_INT_SELF] = zaxisP->self;
  intBuffer[ZAXIS_PACK_INT_TYPE] = zaxisP->type;
  intBuffer[ZAXIS_PACK_INT_SIZE] = zaxisP->size;
  intBuffer[ZAXIS_PACK_INT_DIRECTION] = zaxisP->direction;
  intBuffer[ZAXIS_PACK_INT_VCTSIZE] = zaxisP->vctsize;
  intBuffer[ZAXIS_PACK_INT_MEMBERMASK] = memberMask = zaxisGetMemberMask(zaxisP);

  serializePack(intBuffer, zaxisNint, CDI_DATATYPE_INT, packBuffer, packBufferSize, packBufferPos, context);
  d = cdiCheckSum(CDI_DATATYPE_INT, zaxisNint, intBuffer);
  serializePack(&d, 1, CDI_DATATYPE_UINT32, packBuffer, packBufferSize, packBufferPos, context);

  if (memberMask & vals)
    {
      xassert(zaxisP->size);
      serializePack(zaxisP->vals, zaxisP->size, CDI_DATATYPE_FLT64, packBuffer, packBufferSize, packBufferPos, context);
      d = cdiCheckSum(CDI_DATATYPE_FLT, zaxisP->size, zaxisP->vals);
      serializePack(&d, 1, CDI_DATATYPE_UINT32, packBuffer, packBufferSize, packBufferPos, context);
    }

  if (memberMask & lbounds)
    {
      xassert(zaxisP->size);
      serializePack(zaxisP->lbounds, zaxisP->size, CDI_DATATYPE_FLT64, packBuffer, packBufferSize, packBufferPos, context);
      d = cdiCheckSum(CDI_DATATYPE_FLT, zaxisP->size, zaxisP->lbounds);
      serializePack(&d, 1, CDI_DATATYPE_UINT32, packBuffer, packBufferSize, packBufferPos, context);
    }

  if (memberMask & ubounds)
    {
      xassert(zaxisP->size);

      serializePack(zaxisP->ubounds, zaxisP->size, CDI_DATATYPE_FLT64, packBuffer, packBufferSize, packBufferPos, context);
      d = cdiCheckSum(CDI_DATATYPE_FLT, zaxisP->size, zaxisP->ubounds);
      serializePack(&d, 1, CDI_DATATYPE_UINT32, packBuffer, packBufferSize, packBufferPos, context);
    }

  if (memberMask & weights)
    {
      xassert(zaxisP->size);

      serializePack(zaxisP->weights, zaxisP->size, CDI_DATATYPE_FLT64, packBuffer, packBufferSize, packBufferPos, context);
      d = cdiCheckSum(CDI_DATATYPE_FLT, zaxisP->size, zaxisP->weights);
      serializePack(&d, 1, CDI_DATATYPE_UINT32, packBuffer, packBufferSize, packBufferPos, context);
    }

  if (memberMask & vct)
    {
      xassert(zaxisP->vctsize);

      serializePack(zaxisP->vct, zaxisP->vctsize, CDI_DATATYPE_FLT64, packBuffer, packBufferSize, packBufferPos, context);
      d = cdiCheckSum(CDI_DATATYPE_FLT64, zaxisP->vctsize, zaxisP->vct);
      serializePack(&d, 1, CDI_DATATYPE_UINT32, packBuffer, packBufferSize, packBufferPos, context);
    }

  serializeKeysPack(&zaxisP->keys, packBuffer, packBufferSize, packBufferPos, context);

  serializePack(&zaxisP->positive, 1, CDI_DATATYPE_UINT, packBuffer, packBufferSize, packBufferPos, context);
}

void
cdiZaxisGetIndexList(unsigned nzaxis, int *zaxisResHs)
{
  reshGetResHListOfType(nzaxis, zaxisResHs, &zaxisOps);
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
