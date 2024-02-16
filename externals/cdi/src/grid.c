#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <assert.h>
#include <string.h>

#include "dmemory.h"
#include "cdi.h"
#include "cdi_cksum.h"
#include "cdi_int.h"
#include "cdi_uuid.h"
#include "grid.h"
#include "resource_handle.h"
#include "resource_unpack.h"
#include "namespace.h"
#include "serialize.h"
#include "vlist.h"

int (*proj_lonlat_to_lcc_func)(struct CDI_GridProjParams gpp, size_t, double *, double *) = NULL;
int (*proj_lcc_to_lonlat_func)(struct CDI_GridProjParams gpp, double, double, size_t, double *, double *) = NULL;
int (*proj_lonlat_to_stere_func)(struct CDI_GridProjParams gpp, size_t, double *, double *) = NULL;
int (*proj_stere_to_lonlat_func)(struct CDI_GridProjParams gpp, double, double, size_t, double *, double *) = NULL;

// the value in the second pair of brackets must match the length of the longest string (including terminating NUL)
static const char Grids[][17] = {
  /*  0 */ "undefined",
  /*  1 */ "generic",
  /*  2 */ "gaussian",
  /*  3 */ "gaussian_reduced",
  /*  4 */ "lonlat",
  /*  5 */ "spectral",
  /*  6 */ "fourier",
  /*  7 */ "gme",
  /*  8 */ "trajectory",
  /*  9 */ "unstructured",
  /* 10 */ "curvilinear",
  /* 11 */ "lcc",
  /* 12 */ "projection",
  /* 13 */ "characterXY",
};

// must match table below
enum xystdname_idx
{
  grid_xystdname_grid_latlon,
  grid_xystdname_latlon,
  grid_xystdname_projection,
  grid_xystdname_char,
};
static const char xystdname_tab[][2][24] = {
  [grid_xystdname_grid_latlon] = { "grid_longitude", "grid_latitude" },
  [grid_xystdname_latlon] = { "longitude", "latitude" },
  [grid_xystdname_projection] = { "projection_x_coordinate", "projection_y_coordinate" },
  [grid_xystdname_char] = { "region", "region" },
};

static int gridCompareP(void *gridptr1, void *gridptr2);
static void gridDestroyP(void *gridptr);
static void gridPrintP(void *gridptr, FILE *fp);
static int gridGetPackSize(void *gridptr, void *context);
static void gridPack(void *gridptr, void *buff, int size, int *position, void *context);
static int gridTxCode(void *gridptr);

static const resOps gridOps = { gridCompareP, gridDestroyP, gridPrintP, gridGetPackSize, gridPack, gridTxCode };

grid_t *
grid_to_pointer(int gridID)
{
  return (grid_t *) reshGetVal(gridID, &gridOps);
}

#define gridMark4Update(gridID) reshSetStatus(gridID, &gridOps, RESH_DESYNC_IN_USE)

static inline bool
grid_is_irregular(int gridType)
{
  return (gridType == GRID_UNSTRUCTURED || gridType == GRID_CURVILINEAR);
}

static bool
cdiInqAttConvertedToFloat(int gridID, int atttype, const char *attname, int attlen, double *attflt)
{
  bool status = true;

  if (atttype == CDI_DATATYPE_INT32)
    {
      int attint = 0;
      int *pattint = (attlen > 1) ? (int *) malloc(attlen * sizeof(int)) : &attint;
      cdiInqAttInt(gridID, CDI_GLOBAL, attname, attlen, pattint);
      for (int i = 0; i < attlen; ++i) attflt[i] = (double) pattint[i];
      if (attlen > 1) free(pattint);
    }
  else if (atttype == CDI_DATATYPE_FLT32 || atttype == CDI_DATATYPE_FLT64)
    {
      cdiInqAttFlt(gridID, CDI_GLOBAL, attname, attlen, attflt);
    }
  else
    {
      status = false;
    }

  return status;
}

static void
grid_axis_init(struct gridaxis_t *axisptr)
{
  axisptr->size = 0;
  axisptr->vals = NULL;
  axisptr->bounds = NULL;
  axisptr->flag = 0;
  axisptr->first = 0.0;
  axisptr->last = 0.0;
  axisptr->inc = 0.0;
#ifndef USE_MPI
  axisptr->clength = 0;
  axisptr->cvals = NULL;
#endif
  cdiInitKeys(&axisptr->keys);
}

enum cdiApplyRet
cdiGridApply(enum cdiApplyRet (*func)(int id, void *res, void *data), void *data)
{
  return cdiResHFilterApply(&gridOps, func, data);
}

void
grid_init(grid_t *gridptr)
{
  gridptr->self = CDI_UNDEFID;
  gridptr->type = CDI_UNDEFID;
  gridptr->datatype = CDI_UNDEFID;
  gridptr->proj = CDI_UNDEFID;
  gridptr->projtype = CDI_UNDEFID;
  gridptr->mask = NULL;
  gridptr->mask_gme = NULL;
  gridptr->size = 0;

  grid_axis_init(&gridptr->x);
  grid_axis_init(&gridptr->y);

  gridptr->area = NULL;
  gridptr->reducedPoints = NULL;
  gridptr->reducedPointsSize = 0;

  gridptr->gme.nd = 0;
  gridptr->gme.ni = 0;
  gridptr->gme.ni2 = 0;
  gridptr->gme.ni3 = 0;

  gridptr->trunc = 0;
  gridptr->nvertex = 0;
  gridptr->np = 0;
  gridptr->isCyclic = CDI_UNDEFID;

  gridptr->lcomplex = false;
  gridptr->hasdims = true;
  gridptr->name = NULL;
  gridptr->vtable = &cdiGridVtable;

  cdiInitKeys(&gridptr->keys);
  gridptr->atts.nalloc = MAX_ATTRIBUTES;
  gridptr->atts.nelems = 0;

  cdiDefVarKeyInt(&gridptr->keys, CDI_KEY_DATATYPE, CDI_DATATYPE_FLT64);
#ifdef HIRLAM_EXTENSIONS
  cdiDefVarKeyInt(&gridptr->keys, CDI_KEY_SCANNINGMODE, 64);
#endif

  gridptr->extraData = NULL;
}

static void
grid_free_components(grid_t *gridptr)
{
  void *p2free[] = { gridptr->mask,     gridptr->mask_gme, gridptr->x.vals,        gridptr->y.vals,
#ifndef USE_MPI
                     gridptr->x.cvals,  gridptr->y.cvals,
#endif
                     gridptr->x.bounds, gridptr->y.bounds, gridptr->reducedPoints, gridptr->area,   gridptr->name };

  for (size_t i = 0; i < sizeof(p2free) / sizeof(p2free[0]); ++i)
    if (p2free[i]) Free(p2free[i]);

  cdiDeleteVarKeys(&(gridptr->x.keys));
  cdiDeleteVarKeys(&(gridptr->y.keys));
  cdiDeleteVarKeys(&(gridptr->keys));
  /* 12 pio tests fail
  int gridID = gridptr->self;
  if (gridID != CDI_UNDEFID) cdiDeleteAtts(gridID, CDI_GLOBAL);
  */
}

void
grid_free(grid_t *gridptr)
{
  if (gridptr)
    {
      grid_free_components(gridptr);
      grid_init(gridptr);
    }
}

static grid_t *
gridNewEntry(cdiResH resH)
{
  grid_t *gridptr = (grid_t *) Malloc(sizeof(grid_t));
  grid_init(gridptr);

  if (resH == CDI_UNDEFID)
    gridptr->self = reshPut(gridptr, &gridOps);
  else
    {
      gridptr->self = resH;
      reshReplace(resH, gridptr, &gridOps);
    }

  return gridptr;
}

static void
gridInit(void)
{
  static bool gridInitialized = false;
  if (gridInitialized) return;
  gridInitialized = true;
}

static void
grid_copy_base_scalar_fields(grid_t *gridptrOrig, grid_t *gridptrDup)
{
  memcpy(gridptrDup, gridptrOrig, sizeof(grid_t));
  gridptrDup->self = CDI_UNDEFID;
  cdiInitKeys(&gridptrDup->keys);
  cdiCopyVarKeys(&gridptrOrig->keys, &gridptrDup->keys);
  cdiInitKeys(&gridptrDup->x.keys);
  cdiCopyVarKeys(&gridptrOrig->x.keys, &gridptrDup->x.keys);
  cdiInitKeys(&gridptrDup->y.keys);
  cdiCopyVarKeys(&gridptrOrig->y.keys, &gridptrDup->y.keys);
}

static grid_t *
grid_copy_base(grid_t *gridptrOrig)
{
  grid_t *gridptrDup = (grid_t *) Malloc(sizeof(*gridptrDup));
  gridptrOrig->vtable->copyScalarFields(gridptrOrig, gridptrDup);
  gridptrOrig->vtable->copyArrayFields(gridptrOrig, gridptrDup);
  return gridptrDup;
}

unsigned
cdiGridCount(void)
{
  return reshCountType(&gridOps);
}

static inline void
gridaxisSetKey(struct gridaxis_t *axisptr, int key, const char *name)
{
  if (find_key(&axisptr->keys, key) == NULL)
    cdiDefVarKeyBytes(&axisptr->keys, key, (const unsigned char *) name, (int) strlen(name) + 1);
}

void
cdiGridTypeInit(grid_t *gridptr, int gridtype, size_t size)
{
  gridptr->type = gridtype;
  gridptr->size = size;

  // clang-format off
  if      (gridtype == GRID_LONLAT)           gridptr->nvertex = 2;
  else if (gridtype == GRID_GAUSSIAN)         gridptr->nvertex = 2;
  else if (gridtype == GRID_GAUSSIAN_REDUCED) gridptr->nvertex = 2;
  else if (gridtype == GRID_CURVILINEAR)      gridptr->nvertex = 4;
  else if (gridtype == GRID_UNSTRUCTURED)     gridptr->x.size = size;
  // clang-format on

  switch (gridtype)
    {
    case GRID_LONLAT:
    case GRID_GAUSSIAN:
    case GRID_GAUSSIAN_REDUCED:
    case GRID_TRAJECTORY:
    case GRID_CURVILINEAR:
    case GRID_UNSTRUCTURED:
    case GRID_GME:
      {
        if (gridtype == GRID_TRAJECTORY)
          {
            gridaxisSetKey(&gridptr->x, CDI_KEY_NAME, "tlon");
            gridaxisSetKey(&gridptr->y, CDI_KEY_NAME, "tlat");
          }
        else
          {
            gridaxisSetKey(&gridptr->x, CDI_KEY_NAME, "lon");
            gridaxisSetKey(&gridptr->y, CDI_KEY_NAME, "lat");
          }

        gridaxisSetKey(&gridptr->x, CDI_KEY_LONGNAME, "longitude");
        gridaxisSetKey(&gridptr->y, CDI_KEY_LONGNAME, "latitude");

        gridaxisSetKey(&gridptr->x, CDI_KEY_UNITS, "degrees_east");
        gridaxisSetKey(&gridptr->y, CDI_KEY_UNITS, "degrees_north");

        gridaxisSetKey(&gridptr->x, CDI_KEY_STDNAME, xystdname_tab[grid_xystdname_latlon][0]);
        gridaxisSetKey(&gridptr->y, CDI_KEY_STDNAME, xystdname_tab[grid_xystdname_latlon][1]);

        break;
      }
#ifndef USE_MPI
    case GRID_CHARXY:
      {
        if (gridptr->x.cvals) gridaxisSetKey(&gridptr->x, CDI_KEY_STDNAME, xystdname_tab[grid_xystdname_char][0]);
        if (gridptr->y.cvals) gridaxisSetKey(&gridptr->y, CDI_KEY_STDNAME, xystdname_tab[grid_xystdname_char][1]);

        break;
      }
#endif
    case GRID_GENERIC:
    case GRID_PROJECTION:
      {
        gridaxisSetKey(&gridptr->x, CDI_KEY_NAME, "x");
        gridaxisSetKey(&gridptr->y, CDI_KEY_NAME, "y");
        if (gridtype == GRID_PROJECTION)
          {
            gridaxisSetKey(&gridptr->x, CDI_KEY_STDNAME, xystdname_tab[grid_xystdname_projection][0]);
            gridaxisSetKey(&gridptr->y, CDI_KEY_STDNAME, xystdname_tab[grid_xystdname_projection][1]);
            gridaxisSetKey(&gridptr->x, CDI_KEY_UNITS, "m");
            gridaxisSetKey(&gridptr->y, CDI_KEY_UNITS, "m");
          }
        break;
      }
    }
}

// used also in CDO
void
gridGenXvals(int xsize, double xfirst, double xlast, double xinc, double *restrict xvals)
{
  if (fabs(xinc) <= 0 && xsize > 1)
    {
      if (xfirst >= xlast)
        {
          while (xfirst >= xlast) xlast += 360;
          xinc = (xlast - xfirst) / (xsize);
        }
      else
        {
          xinc = (xlast - xfirst) / (xsize - 1);
        }
    }

  for (int i = 0; i < xsize; ++i) xvals[i] = xfirst + i * xinc;
}

static void
calc_gaussgrid(double *restrict yvals, int ysize, double yfirst, double ylast)
{
  double *restrict yw = (double *) malloc((size_t) ysize * sizeof(double));
  gaussianLatitudes((size_t) ysize, yvals, yw);
  free(yw);
  for (int i = 0; i < ysize; i++) yvals[i] = asin(yvals[i]) / M_PI * 180.0;

  if (yfirst < ylast && yfirst > -90.0 && ylast < 90.0)
    {
      int yhsize = ysize / 2;
      for (int i = 0; i < yhsize; i++)
        {
          const double ytmp = yvals[i];
          yvals[i] = yvals[ysize - i - 1];
          yvals[ysize - i - 1] = ytmp;
        }
    }
}

static void
gridGenYvalsGaussian(int ysize, double yfirst, double ylast, double *restrict yvals)
{
  const double deleps = 0.002;

  calc_gaussgrid(yvals, ysize, yfirst, ylast);

  if (!(IS_EQUAL(yfirst, 0) && IS_EQUAL(ylast, 0)))
    if (fabs(yvals[0] - yfirst) > deleps || fabs(yvals[ysize - 1] - ylast) > deleps)
      {
        bool lfound = false;
        int ny = (int) (180. / (fabs(ylast - yfirst) / (ysize - 1)) + 0.5);
        ny -= ny % 2;
        if (ny > ysize && ny < 4096)
          {
            double *ytmp = (double *) Malloc((size_t) ny * sizeof(double));
            calc_gaussgrid(ytmp, ny, yfirst, ylast);

            int i;
            for (i = 0; i < (ny - ysize); i++)
              if (fabs(ytmp[i] - yfirst) < deleps) break;
            int nstart = i;

            lfound = (nstart + ysize - 1) < ny && fabs(ytmp[nstart + ysize - 1] - ylast) < deleps;
            if (lfound)
              {
                for (i = 0; i < ysize; i++) yvals[i] = ytmp[i + nstart];
              }

            if (ytmp) Free(ytmp);
          }

        if (!lfound)
          {
            Warning("Cannot calculate gaussian latitudes for lat1 = %g latn = %g!", yfirst, ylast);
            for (int i = 0; i < ysize; i++) yvals[i] = 0;
            yvals[0] = yfirst;
            yvals[ysize - 1] = ylast;
          }
      }
}

static void
gridGenYvalsRegular(int ysize, double yfirst, double ylast, double yinc, double *restrict yvals)
{
  if (fabs(yinc) <= 0 && ysize > 1)
    {
      if (IS_EQUAL(yfirst, ylast) && IS_NOT_EQUAL(yfirst, 0)) ylast *= -1;

      if (yfirst > ylast)
        yinc = (yfirst - ylast) / (ysize - 1);
      else if (yfirst < ylast)
        yinc = (ylast - yfirst) / (ysize - 1);
      else
        {
          if (ysize % 2 != 0)
            {
              yinc = 180.0 / (ysize - 1);
              yfirst = -90;
            }
          else
            {
              yinc = 180.0 / ysize;
              yfirst = -90 + yinc / 2;
            }
        }
    }

  if (yfirst > ylast && yinc > 0) yinc = -yinc;

  for (int i = 0; i < ysize; i++) yvals[i] = yfirst + i * yinc;
}

// used also in CDO
void
gridGenYvals(int gridtype, int ysize, double yfirst, double ylast, double yinc, double *restrict yvals)
{
  if (gridtype == GRID_GAUSSIAN || gridtype == GRID_GAUSSIAN_REDUCED)
    {
      if (ysize > 2)
        {
          gridGenYvalsGaussian(ysize, yfirst, ylast, yvals);
        }
      else
        {
          yvals[0] = yfirst;
          yvals[ysize - 1] = ylast;
        }
    }
  // else if (gridtype == GRID_LONLAT || gridtype == GRID_GENERIC)
  else
    {
      gridGenYvalsRegular(ysize, yfirst, ylast, yinc, yvals);
    }
  /*
    else
    Error("unable to calculate values for %s grid!", gridNamePtr(gridtype));
  */
}

/*
@Function  gridCreate
@Title     Create a horizontal Grid

@Prototype int gridCreate(int gridtype, SizeType size)
@Parameter
    @Item  gridtype  The type of the grid, one of the set of predefined CDI grid types.
                     The valid CDI grid types are @func{GRID_GENERIC}, @func{GRID_LONLAT},
                     @func{GRID_GAUSSIAN}, @func{GRID_PROJECTION}, @func{GRID_SPECTRAL},
                     @func{GRID_GME}, @func{GRID_CURVILINEAR} and @func{GRID_UNSTRUCTURED}.
    @Item  size      Number of gridpoints.

@Description
The function @func{gridCreate} creates a horizontal Grid.

@Result
@func{gridCreate} returns an identifier to the Grid.

@Example
Here is an example using @func{gridCreate} to create a regular lon/lat Grid:

@Source
#include "cdi.h"
   ...
#define  nlon  12
#define  nlat   6
   ...
double lons[nlon] = {0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330};
double lats[nlat] = {-75, -45, -15, 15, 45, 75};
int gridID;
   ...
gridID = gridCreate(GRID_LONLAT, nlon*nlat);
gridDefXsize(gridID, nlon);
gridDefYsize(gridID, nlat);
gridDefXvals(gridID, lons);
gridDefYvals(gridID, lats);
   ...
@EndSource
@EndFunction
*/
int
gridCreate(int gridtype, SizeType size)
{
  if (CDI_Debug) Message("gridtype=%s  size=%zu", gridNamePtr(gridtype), size);

  xassert(size);
  gridInit();

  grid_t *gridptr = gridNewEntry(CDI_UNDEFID);
  if (!gridptr) Error("No memory");

  int gridID = gridptr->self;

  if (CDI_Debug) Message("gridID: %d", gridID);

  cdiGridTypeInit(gridptr, gridtype, (size_t) size);

  return gridID;
}

static void
gridDestroyKernel(grid_t *gridptr)
{
  xassert(gridptr);

  grid_free_components(gridptr);
  Free(gridptr);
}

/*
@Function  gridDestroy
@Title     Destroy a horizontal Grid

@Prototype void gridDestroy(int gridID)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate}.

@EndFunction
*/
void
gridDestroy(int gridID)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  gridptr->vtable->destroy(gridptr);
  reshRemove(gridID, &gridOps);
}

static void
gridDestroyP(void *gridptr)
{
  ((grid_t *) gridptr)->vtable->destroy((grid_t *) gridptr);
}

const char *
gridNamePtr(int gridtype)
{
  int size = (int) (sizeof(Grids) / sizeof(Grids[0]));

  const char *name = (gridtype >= 0 && gridtype < size) ? Grids[gridtype] : Grids[GRID_GENERIC];

  return name;
}

void
gridName(int gridtype, char *gridname)
{
  strcpy(gridname, gridNamePtr(gridtype));
}

/*
@Function  gridDefXname
@Title     Define the name of a X-axis

@Prototype void gridDefXname(int gridID, const char *name)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate}.
    @Item  name     Name of the X-axis.

@Description
The function @func{gridDefXname} defines the name of a X-axis.

@EndFunction
*/
void
gridDefXname(int gridID, const char *name)
{
  (void) cdiDefKeyString(gridID, CDI_XAXIS, CDI_KEY_NAME, name);
}

/*
@Function  gridDefXlongname
@Title     Define the longname of a X-axis

@Prototype void gridDefXlongname(int gridID, const char *longname)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate}.
    @Item  longname Longname of the X-axis.

@Description
The function @func{gridDefXlongname} defines the longname of a X-axis.

@EndFunction
*/
void
gridDefXlongname(int gridID, const char *longname)
{
  (void) cdiDefKeyString(gridID, CDI_XAXIS, CDI_KEY_LONGNAME, longname);
}

/*
@Function  gridDefXunits
@Title     Define the units of a X-axis

@Prototype void gridDefXunits(int gridID, const char *units)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate}.
    @Item  units    Units of the X-axis.

@Description
The function @func{gridDefXunits} defines the units of a X-axis.

@EndFunction
*/
void
gridDefXunits(int gridID, const char *units)
{
  (void) cdiDefKeyString(gridID, CDI_XAXIS, CDI_KEY_UNITS, units);
}

/*
@Function  gridDefYname
@Title     Define the name of a Y-axis

@Prototype void gridDefYname(int gridID, const char *name)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate}.
    @Item  name     Name of the Y-axis.

@Description
The function @func{gridDefYname} defines the name of a Y-axis.

@EndFunction
*/
void
gridDefYname(int gridID, const char *name)
{
  (void) cdiDefKeyString(gridID, CDI_YAXIS, CDI_KEY_NAME, name);
}

/*
@Function  gridDefYlongname
@Title     Define the longname of a Y-axis

@Prototype void gridDefYlongname(int gridID, const char *longname)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate}.
    @Item  longname Longname of the Y-axis.

@Description
The function @func{gridDefYlongname} defines the longname of a Y-axis.

@EndFunction
*/
void
gridDefYlongname(int gridID, const char *longname)
{
  (void) cdiDefKeyString(gridID, CDI_YAXIS, CDI_KEY_LONGNAME, longname);
}

/*
@Function  gridDefYunits
@Title     Define the units of a Y-axis

@Prototype void gridDefYunits(int gridID, const char *units)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate}.
    @Item  units    Units of the Y-axis.

@Description
The function @func{gridDefYunits} defines the units of a Y-axis.

@EndFunction
*/
void
gridDefYunits(int gridID, const char *units)
{
  (void) cdiDefKeyString(gridID, CDI_YAXIS, CDI_KEY_UNITS, units);
}

/*
@Function  gridInqXname
@Title     Get the name of a X-axis

@Prototype void gridInqXname(int gridID, char *name)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate} or @fref{vlistInqVarGrid}.
    @Item  name     Name of the X-axis. The caller must allocate space for the
                    returned string. The maximum possible length, in characters, of
                    the string is given by the predefined constant @func{CDI_MAX_NAME}.

@Description
The function @func{gridInqXname} returns the name of a X-axis.

@Result
@func{gridInqXname} returns the name of the X-axis to the parameter name.

@EndFunction
*/
void
gridInqXname(int gridID, char *name)
{
  int length = CDI_MAX_NAME;
  (void) cdiInqKeyString(gridID, CDI_XAXIS, CDI_KEY_NAME, name, &length);
}

/*
@Function  gridInqXlongname
@Title     Get the longname of a X-axis

@Prototype void gridInqXlongname(int gridID, char *longname)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate} or @fref{vlistInqVarGrid}.
    @Item  longname Longname of the X-axis. The caller must allocate space for the
                    returned string. The maximum possible length, in characters, of
                    the string is given by the predefined constant @func{CDI_MAX_NAME}.

@Description
The function @func{gridInqXlongname} returns the longname of a X-axis.

@Result
@func{gridInqXlongname} returns the longname of the X-axis to the parameter longname.

@EndFunction
*/
void
gridInqXlongname(int gridID, char *longname)
{
  int length = CDI_MAX_NAME;
  (void) cdiInqKeyString(gridID, CDI_XAXIS, CDI_KEY_LONGNAME, longname, &length);
}

/*
@Function  gridInqXunits
@Title     Get the units of a X-axis

@Prototype void gridInqXunits(int gridID, char *units)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate} or @fref{vlistInqVarGrid}.
    @Item  units    Units of the X-axis. The caller must allocate space for the
                    returned string. The maximum possible length, in characters, of
                    the string is given by the predefined constant @func{CDI_MAX_NAME}.

@Description
The function @func{gridInqXunits} returns the units of a X-axis.

@Result
@func{gridInqXunits} returns the units of the X-axis to the parameter units.

@EndFunction
*/
void
gridInqXunits(int gridID, char *units)
{
  int length = CDI_MAX_NAME;
  (void) cdiInqKeyString(gridID, CDI_XAXIS, CDI_KEY_UNITS, units, &length);
}

/*
@Function  gridInqYname
@Title     Get the name of a Y-axis

@Prototype void gridInqYname(int gridID, char *name)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate} or @fref{vlistInqVarGrid}.
    @Item  name     Name of the Y-axis. The caller must allocate space for the
                    returned string. The maximum possible length, in characters, of
                    the string is given by the predefined constant @func{CDI_MAX_NAME}.

@Description
The function @func{gridInqYname} returns the name of a Y-axis.

@Result
@func{gridInqYname} returns the name of the Y-axis to the parameter name.

@EndFunction
*/
void
gridInqYname(int gridID, char *name)
{
  int length = CDI_MAX_NAME;
  (void) cdiInqKeyString(gridID, CDI_YAXIS, CDI_KEY_NAME, name, &length);
}

/*
@Function  gridInqYlongname
@Title     Get the longname of a Y-axis

@Prototype void gridInqYlongname(int gridID, char *longname)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate} or @fref{vlistInqVarGrid}.
    @Item  longname Longname of the Y-axis. The caller must allocate space for the
                    returned string. The maximum possible length, in characters, of
                    the string is given by the predefined constant @func{CDI_MAX_NAME}.

@Description
The function @func{gridInqYlongname} returns the longname of a Y-axis.

@Result
@func{gridInqYlongname} returns the longname of the Y-axis to the parameter longname.

@EndFunction
*/
void
gridInqYlongname(int gridID, char *longname)
{
  int length = CDI_MAX_NAME;
  (void) cdiInqKeyString(gridID, CDI_YAXIS, CDI_KEY_LONGNAME, longname, &length);
}

/*
@Function  gridInqYunits
@Title     Get the units of a Y-axis

@Prototype void gridInqYunits(int gridID, char *units)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate} or @fref{vlistInqVarGrid}.
    @Item  units    Units of the Y-axis. The caller must allocate space for the
                    returned string. The maximum possible length, in characters, of
                    the string is given by the predefined constant @func{CDI_MAX_NAME}.

@Description
The function @func{gridInqYunits} returns the units of a Y-axis.

@Result
@func{gridInqYunits} returns the units of the Y-axis to the parameter units.

@EndFunction
*/
void
gridInqYunits(int gridID, char *units)
{
  int length = CDI_MAX_NAME;
  (void) cdiInqKeyString(gridID, CDI_YAXIS, CDI_KEY_UNITS, units, &length);
}

void
gridDefProj(int gridID, int projID)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  gridptr->proj = projID;

  if (gridptr->type == GRID_CURVILINEAR)
    {
      grid_t *projptr = grid_to_pointer(projID);
      const char *xdimname = cdiInqVarKeyStringPtr(&gridptr->x.keys, CDI_KEY_DIMNAME);
      const char *ydimname = cdiInqVarKeyStringPtr(&gridptr->y.keys, CDI_KEY_DIMNAME);
      if (xdimname && find_key(&projptr->x.keys, CDI_KEY_NAME)) cdiDefKeyString(projID, CDI_XAXIS, CDI_KEY_NAME, xdimname);
      if (ydimname && find_key(&projptr->y.keys, CDI_KEY_NAME)) cdiDefKeyString(projID, CDI_YAXIS, CDI_KEY_NAME, ydimname);
    }
}

int
gridInqProj(int gridID)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  return gridptr->proj;
}

int
gridInqProjType(int gridID)
{
  grid_t *gridptr = grid_to_pointer(gridID);

  int projtype = gridptr->projtype;
  if (projtype == -1)
    {
      char gmapname[CDI_MAX_NAME];
      int length = CDI_MAX_NAME;
      cdiInqKeyString(gridID, CDI_GLOBAL, CDI_KEY_GRIDMAP_NAME, gmapname, &length);
      if (gmapname[0])
        {
          // clang-format off
          if      (str_is_equal(gmapname, "rotated_latitude_longitude"))   projtype = CDI_PROJ_RLL;
          else if (str_is_equal(gmapname, "lambert_azimuthal_equal_area")) projtype = CDI_PROJ_LAEA;
          else if (str_is_equal(gmapname, "lambert_conformal_conic"))      projtype = CDI_PROJ_LCC;
          else if (str_is_equal(gmapname, "sinusoidal"))                   projtype = CDI_PROJ_SINU;
          else if (str_is_equal(gmapname, "polar_stereographic"))          projtype = CDI_PROJ_STERE;
          else if (str_is_equal(gmapname, "healpix"))                      projtype = CDI_PROJ_HEALPIX;
          // clang-format on
          gridptr->projtype = projtype;
        }
    }

  return projtype;
}

void
gridVerifyProj(int gridID)
{
  grid_t *gridptr = grid_to_pointer(gridID);

  int projtype = gridInqProjType(gridID);
  if (projtype == CDI_PROJ_RLL)
    {
      gridaxisSetKey(&gridptr->x, CDI_KEY_STDNAME, xystdname_tab[grid_xystdname_grid_latlon][0]);
      gridaxisSetKey(&gridptr->y, CDI_KEY_STDNAME, xystdname_tab[grid_xystdname_grid_latlon][1]);
      gridaxisSetKey(&gridptr->x, CDI_KEY_UNITS, "degrees");
      gridaxisSetKey(&gridptr->y, CDI_KEY_UNITS, "degrees");
    }
  else if (projtype == CDI_PROJ_LCC)
    {
      gridaxisSetKey(&gridptr->x, CDI_KEY_STDNAME, xystdname_tab[grid_xystdname_projection][0]);
      gridaxisSetKey(&gridptr->y, CDI_KEY_STDNAME, xystdname_tab[grid_xystdname_projection][1]);
      gridaxisSetKey(&gridptr->x, CDI_KEY_UNITS, "m");
      gridaxisSetKey(&gridptr->y, CDI_KEY_UNITS, "m");
    }
}

/*
@Function  gridInqType
@Title     Get the type of a Grid

@Prototype int gridInqType(int gridID)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate} or @fref{vlistInqVarGrid}.

@Description
The function @func{gridInqType} returns the type of a Grid.

@Result
@func{gridInqType} returns the type of the grid,
one of the set of predefined CDI grid types.
The valid CDI grid types are @func{GRID_GENERIC}, @func{GRID_LONLAT},
@func{GRID_GAUSSIAN}, @func{GRID_PROJECTION}, @func{GRID_SPECTRAL}, @func{GRID_GME},
@func{GRID_CURVILINEAR} and @func{GRID_UNSTRUCTURED}.

@EndFunction
*/
int
gridInqType(int gridID)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  return gridptr->type;
}

/*
@Function  gridInqSize
@Title     Get the size of a Grid

@Prototype SizeType gridInqSize(int gridID)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate} or @fref{vlistInqVarGrid}.

@Description
The function @func{gridInqSize} returns the size of a Grid.

@Result
@func{gridInqSize} returns the number of grid points of a Grid.

@EndFunction
*/
SizeType
gridInqSize(int gridID)
{
  grid_t *gridptr = grid_to_pointer(gridID);

  size_t size = gridptr->size;
  if (size == 0)
    {
      size_t xsize = gridptr->x.size;
      size_t ysize = gridptr->y.size;

      size = ysize ? xsize * ysize : xsize;

      gridptr->size = size;
    }

  return (SizeType) size;
}

static int
nsp2trunc(int nsp)
{
  /*  nsp = (trunc+1)*(trunc+1)              */
  /*      => trunc^2 + 3*trunc - (x-2) = 0   */
  /*                                         */
  /*  with:  y^2 + p*y + q = 0               */
  /*         y = -p/2 +- sqrt((p/2)^2 - q)   */
  /*         p = 3 and q = - (x-2)           */
  int trunc = (int) (sqrt(nsp * 4 + 1.) - 3) / 2;
  return trunc;
}

int
gridInqTrunc(int gridID)
{
  grid_t *gridptr = grid_to_pointer(gridID);

  if (gridptr->trunc == 0)
    {
      if (gridptr->type == GRID_SPECTRAL) gridptr->trunc = nsp2trunc(gridptr->size);
      /*
      else if      ( gridptr->type == GRID_GAUSSIAN )
        gridptr->trunc = nlat2trunc(gridptr->y.size);
      */
    }

  return gridptr->trunc;
}

void
gridDefTrunc(int gridID, int trunc)
{
  grid_t *gridptr = grid_to_pointer(gridID);

  if (gridptr->trunc != trunc)
    {
      gridMark4Update(gridID);
      gridptr->trunc = trunc;
    }
}

/*
@Function  gridDefXsize
@Title     Define the number of values of a X-axis

@Prototype void gridDefXsize(int gridID, SizeType xsize)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate}.
    @Item  xsize    Number of values of a X-axis.

@Description
The function @func{gridDefXsize} defines the number of values of a X-axis.

@EndFunction
*/
void
gridDefXsize(int gridID, SizeType xsize)
{
  grid_t *gridptr = grid_to_pointer(gridID);

  size_t gridSize = gridInqSize(gridID);
  if ((size_t) xsize > gridSize) Error("xsize %zu is greater then gridsize %zu", (size_t) xsize, gridSize);

  int gridType = gridInqType(gridID);
  if (gridType == GRID_UNSTRUCTURED && (size_t) xsize != gridSize)
    Error("xsize %zu must be equal to gridsize %zu for gridtype: %s", (size_t) xsize, gridSize, gridNamePtr(gridType));
  if (gridType == GRID_GAUSSIAN_REDUCED && xsize != 2 && (size_t) xsize != gridSize)
    Error("xsize %zu must be equal to gridsize %zu for gridtype: %s", (size_t) xsize, gridSize, gridNamePtr(gridType));

  if (gridptr->x.size != (size_t) xsize)
    {
      gridMark4Update(gridID);
      gridptr->x.size = (size_t) xsize;
    }

  if (gridType != GRID_UNSTRUCTURED && gridType != GRID_GAUSSIAN_REDUCED && gridType != GRID_PROJECTION)
    {
      size_t axisproduct = gridptr->x.size * gridptr->y.size;
      if (axisproduct > 0 && axisproduct != gridSize)
        Error("Inconsistent grid declaration! (xsize=%zu ysize=%zu gridsize=%zu)", gridptr->x.size, gridptr->y.size, gridSize);
    }
}

void
gridDefDatatype(int gridID, int datatype)
{
  cdiDefKeyInt(gridID, CDI_GLOBAL, CDI_KEY_DATATYPE, datatype);
}

int
gridInqDatatype(int gridID)
{
  int datatype = CDI_UNDEFID;
  cdiInqKeyInt(gridID, CDI_GLOBAL, CDI_KEY_DATATYPE, &datatype);
  return datatype;
}

/*
@Function  gridInqXsize
@Title     Get the number of values of a X-axis

@Prototype SizeType gridInqXsize(int gridID)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate} or @fref{vlistInqVarGrid}.

@Description
The function @func{gridInqXsize} returns the number of values of a X-axis.

@Result
@func{gridInqXsize} returns the number of values of a X-axis.

@EndFunction
*/
SizeType
gridInqXsize(int gridID)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  return (SizeType) gridptr->x.size;
}

/*
@Function  gridDefYsize
@Title     Define the number of values of a Y-axis

@Prototype void gridDefYsize(int gridID, SizeType ysize)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate}.
    @Item  ysize    Number of values of a Y-axis.

@Description
The function @func{gridDefYsize} defines the number of values of a Y-axis.

@EndFunction
*/
void
gridDefYsize(int gridID, SizeType ysize)
{
  grid_t *gridptr = grid_to_pointer(gridID);

  size_t gridSize = gridInqSize(gridID);

  if ((size_t) ysize > gridSize) Error("ysize %zu is greater then gridsize %zu", (size_t) ysize, gridSize);

  int gridType = gridInqType(gridID);
  if (gridType == GRID_UNSTRUCTURED && (size_t) ysize != gridSize)
    Error("ysize %zu must be equal gridsize %zu for gridtype: %s", gridNamePtr(gridType), (size_t) ysize, gridSize);

  if (gridptr->y.size != (size_t) ysize)
    {
      gridMark4Update(gridID);
      gridptr->y.size = (size_t) ysize;
    }

  if (gridType != GRID_UNSTRUCTURED && gridType != GRID_GAUSSIAN_REDUCED && gridType != GRID_PROJECTION)
    {
      size_t axisproduct = gridptr->x.size * gridptr->y.size;
      if (axisproduct > 0 && axisproduct != gridSize)
        Error("Inconsistent grid declaration! (xsize=%zu ysize=%zu gridsize=%zu)", gridptr->x.size, gridptr->y.size, gridSize);
    }
}

/*
@Function  gridInqYsize
@Title     Get the number of values of a Y-axis

@Prototype SizeType gridInqYsize(int gridID)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate} or @fref{vlistInqVarGrid}.

@Description
The function @func{gridInqYsize} returns the number of values of a Y-axis.

@Result
@func{gridInqYsize} returns the number of values of a Y-axis.

@EndFunction
*/
SizeType
gridInqYsize(int gridID)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  return (SizeType) gridptr->y.size;
}

/*
@Function  gridDefNP
@Title     Define the number of parallels between a pole and the equator

@Prototype void gridDefNP(int gridID, int np)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate}.
    @Item  np       Number of parallels between a pole and the equator.

@Description
The function @func{gridDefNP} defines the number of parallels between a pole and the equator
of a Gaussian grid.

@EndFunction
*/
void
gridDefNP(int gridID, int np)
{
  grid_t *gridptr = grid_to_pointer(gridID);

  if (gridptr->np != np)
    {
      gridMark4Update(gridID);
      gridptr->np = np;
    }
}

/*
@Function  gridInqNP
@Title     Get the number of parallels between a pole and the equator

@Prototype int gridInqNP(int gridID)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate} or @fref{vlistInqVarGrid}.

@Description
The function @func{gridInqNP} returns the number of parallels between a pole and the equator
of a Gaussian grid.

@Result
@func{gridInqNP} returns the number of parallels between a pole and the equator.

@EndFunction
*/
int
gridInqNP(int gridID)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  return gridptr->np;
}

/*
@Function
@Title

@Prototype
@Parameter
    @Item  Grid identifier

@EndFunction
*/
void
gridDefReducedPoints(int gridID, int reducedPointsSize, const int reducedPoints[])
{
  grid_t *gridptr = grid_to_pointer(gridID);

  gridptr->reducedPoints = (int *) Malloc((size_t) reducedPointsSize * sizeof(int));
  gridptr->reducedPointsSize = reducedPointsSize;
  memcpy(gridptr->reducedPoints, reducedPoints, (size_t) reducedPointsSize * sizeof(int));
  gridMark4Update(gridID);
}

/*
@Function
@Title

@Prototype
@Parameter
    @Item  Grid identifier

@EndFunction
*/
void
gridInqReducedPoints(int gridID, int *reducedPoints)
{
  grid_t *gridptr = grid_to_pointer(gridID);

  if (gridptr->reducedPoints == 0) Error("undefined pointer!");

  memcpy(reducedPoints, gridptr->reducedPoints, (size_t) gridptr->reducedPointsSize * sizeof(int));
}

static size_t
gridInqMaskSerialGeneric(grid_t *gridptr, mask_t **internalMask, int *restrict mask)
{
  size_t size = gridptr->size;

  if (CDI_Debug && size == 0) Warning("Size undefined for gridID = %d", gridptr->self);

  const mask_t *restrict mask_src = *internalMask;
  if (mask_src)
    {
      if (mask && size > 0)
        for (size_t i = 0; i < size; ++i) mask[i] = (int) mask_src[i];
    }
  else
    size = 0;

  return size;
}

static SizeType
gridInqMaskSerial(grid_t *gridptr, int *mask)
{
  return (SizeType) gridInqMaskSerialGeneric(gridptr, &gridptr->mask, mask);
}

int
gridInqMask(int gridID, int *mask)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  return gridptr->vtable->inqMask(gridptr, mask);
}

static void
gridDefMaskSerial(grid_t *gridptr, const int *mask)
{
  size_t size = gridptr->size;
  if (size == 0) Error("Size undefined for gridID = %d", gridptr->self);

  if (mask == NULL)
    {
      if (gridptr->mask)
        {
          Free(gridptr->mask);
          gridptr->mask = NULL;
        }
    }
  else
    {
      if (gridptr->mask == NULL)
        gridptr->mask = (mask_t *) Malloc(size * sizeof(mask_t));
      else if (CDI_Debug)
        Warning("grid mask already defined!");

      for (size_t i = 0; i < size; ++i) gridptr->mask[i] = (mask_t) (mask[i] != 0);
    }
}

void
gridDefMask(int gridID, const int *mask)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  gridptr->vtable->defMask(gridptr, mask);
  gridMark4Update(gridID);
}

static int
gridInqMaskGMESerial(grid_t *gridptr, int *mask_gme)
{
  return gridInqMaskSerialGeneric(gridptr, &gridptr->mask_gme, mask_gme);
}

int
gridInqMaskGME(int gridID, int *mask)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  return gridptr->vtable->inqMaskGME(gridptr, mask);
}

static void
gridDefMaskGMESerial(grid_t *gridptr, const int *mask)
{
  size_t size = gridptr->size;
  if (size == 0) Error("Size undefined for gridID = %d", gridptr->self);

  if (gridptr->mask_gme == NULL)
    gridptr->mask_gme = (mask_t *) Malloc(size * sizeof(mask_t));
  else if (CDI_Debug)
    Warning("mask already defined!");

  for (size_t i = 0; i < size; ++i) gridptr->mask_gme[i] = (mask_t) (mask[i] != 0);
}

void
gridDefMaskGME(int gridID, const int *mask)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  gridptr->vtable->defMaskGME(gridptr, mask);
  gridMark4Update(gridID);
}

static void
copy_darray(size_t n, const double *restrict in, double *restrict out)
{
#ifdef _OPENMP
#pragma omp parallel for if (n > 99999) default(shared) schedule(static)
#endif
  for (size_t i = 0; i < n; ++i) out[i] = in[i];
}

static SizeType
gridInqXValsSerial(grid_t *gridptr, double *xvals)
{
  int gridtype = gridptr->type;
  size_t size = grid_is_irregular(gridtype) ? gridptr->size : gridptr->x.size;

  if (CDI_Debug && size == 0) Warning("size undefined for gridID = %d", gridptr->self);

  if (gridptr->x.vals)
    {
      if (size && xvals)
        {
          const double *gridptr_xvals = gridptr->vtable->inqXValsPtr(gridptr);
          copy_darray(size, gridptr_xvals, xvals);
        }
    }
  else
    size = 0;

  return (SizeType) size;
}

static SizeType
gridInqXValsPartSerial(grid_t *gridptr, int start, SizeType length, double *xvals)
{
  int gridtype = gridptr->type;
  size_t size = grid_is_irregular(gridtype) ? gridptr->size : gridptr->x.size;

  if (CDI_Debug && size == 0) Warning("size undefined for gridID = %d", gridptr->self);

  if (gridptr->x.vals)
    {
      if (size && xvals && (size_t) length <= size)
        {
          const double *gridptr_xvals = gridptr->vtable->inqXValsPtr(gridptr);
          memcpy(xvals, gridptr_xvals + start, (size_t) length * sizeof(double));
        }
    }
  else
    size = 0;

  return (SizeType) size;
}

#ifndef USE_MPI
static SizeType
gridInqXCvalsSerial(grid_t *gridptr, char **xcvals)
{
  if (gridptr->type != GRID_CHARXY) Error("Function only valid for grid type 'GRID_CHARXY'.");

  size_t size = gridptr->x.size;
  size_t maxclength = 0;

  const char **gridptr_xcvals = gridptr->vtable->inqXCvalsPtr(gridptr);
  if (gridptr_xcvals && size && xcvals)
    {
      maxclength = gridptr->x.clength;
      for (size_t i = 0; i < size; i++) memcpy(xcvals[i], gridptr_xcvals[i], maxclength * sizeof(char));
    }

  return (SizeType) maxclength;
}

static int
gridInqXIscSerial(grid_t *gridptr)
{
  /*
  if ( gridptr->type != GRID_CHARXY )
    Error("Axis type is 'char' but grid is not type 'GRID_CHARXY'.");
  */
  return gridptr->x.clength;
}
#endif

/*
@Function  gridInqXvals
@Title     Get all values of a X-axis

@Prototype SizeType gridInqXvals(int gridID, double *xvals)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate} or @fref{vlistInqVarGrid}.
    @Item  xvals    Pointer to the location into which the X-values are read.
                    The caller must allocate space for the returned values.

@Description
The function @func{gridInqXvals} returns all values of the X-axis.

@Result
Upon successful completion @func{gridInqXvals} returns the number of values and
the values are stored in @func{xvals}.
Otherwise, 0 is returned and @func{xvals} is empty.

@EndFunction
*/
SizeType
gridInqXvals(int gridID, double *xvals)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  return gridptr->vtable->inqXVals(gridptr, xvals);
}

SizeType
gridInqXvalsPart(int gridID, int start, SizeType length, double *xvals)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  return gridptr->vtable->inqXValsPart(gridptr, start, length, xvals);
}

SizeType
gridInqXCvals(int gridID, char **xcvals)
{
  grid_t *gridptr = grid_to_pointer(gridID);
#ifndef USE_MPI
  return gridptr->vtable->inqXCvals(gridptr, xcvals);
#else
  return 0;
#endif
}

int
gridInqXIsc(int gridID)
{
  grid_t *gridptr = grid_to_pointer(gridID);
#ifndef USE_MPI
  return gridptr->vtable->inqXIsc(gridptr);
#else
  return 0;
#endif
}

static void
gridDefXValsSerial(grid_t *gridptr, const double *xvals)
{
  int gridtype = gridptr->type;
  size_t size = grid_is_irregular(gridtype) ? gridptr->size : gridptr->x.size;

  if (size == 0) Error("Size undefined for gridID = %d", gridptr->self);

  if (gridptr->x.vals && CDI_Debug) Warning("values already defined!");
  gridptr->x.vals = (double *) Realloc(gridptr->x.vals, size * sizeof(double));
  copy_darray(size, xvals, gridptr->x.vals);
}

#ifndef USE_MPI
static SizeType
gridInqYCvalsSerial(grid_t *gridptr, char **ycvals)
{
  if (gridptr->type != GRID_CHARXY) Error("Function only valid for grid type 'GRID_CHARXY'.");

  size_t size = gridptr->y.size;
  size_t maxclength = 0;

  const char **gridptr_ycvals = gridptr->vtable->inqYCvalsPtr(gridptr);
  if (gridptr_ycvals && size && ycvals)
    {
      maxclength = gridptr->y.clength;
      for (size_t i = 0; i < size; i++) memcpy(ycvals[i], gridptr_ycvals[i], maxclength * sizeof(char));
    }

  return (SizeType) maxclength;
}

static int
gridInqYIscSerial(grid_t *gridptr)
{
  // if ( gridptr->type != GRID_CHARXY ) Error("Axis type is 'char' but grid is not type 'GRID_CHARXY'.");
  return gridptr->y.clength;
}
#endif

/*
@Function  gridDefXvals
@Title     Define the values of a X-axis

@Prototype void gridDefXvals(int gridID, const double *xvals)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate}.
    @Item  xvals    X-values of the grid.

@Description
The function @func{gridDefXvals} defines all values of the X-axis.

@EndFunction
*/
void
gridDefXvals(int gridID, const double *xvals)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  gridptr->vtable->defXVals(gridptr, xvals);
  gridMark4Update(gridID);
}

static SizeType
gridInqYValsSerial(grid_t *gridptr, double *yvals)
{
  int gridtype = gridptr->type;
  size_t size = grid_is_irregular(gridtype) ? gridptr->size : gridptr->y.size;

  if (CDI_Debug && size == 0) Warning("size undefined for gridID = %d!", gridptr->self);

  if (gridptr->y.vals)
    {
      if (size && yvals)
        {
          const double *gridptr_yvals = gridptr->vtable->inqYValsPtr(gridptr);
          copy_darray(size, gridptr_yvals, yvals);
        }
    }
  else
    size = 0;

  return (SizeType) size;
}

static SizeType
gridInqYValsPartSerial(grid_t *gridptr, int start, SizeType length, double *yvals)
{
  int gridtype = gridptr->type;
  size_t size = grid_is_irregular(gridtype) ? gridptr->size : gridptr->y.size;

  if (CDI_Debug && size == 0) Warning("size undefined for gridID = %d!", gridptr->self);

  if (gridptr->y.vals)
    {
      if (size && yvals && (size_t) length <= size)
        {
          const double *gridptr_yvals = gridptr->vtable->inqYValsPtr(gridptr);
          memcpy(yvals, gridptr_yvals + start, (size_t) length * sizeof(double));
        }
    }
  else
    size = 0;

  return (SizeType) size;
}

/*
@Function  gridInqYvals
@Title     Get all values of a Y-axis

@Prototype SizeType gridInqYvals(int gridID, double *yvals)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate} or @fref{vlistInqVarGrid}.
    @Item  yvals    Pointer to the location into which the Y-values are read.
                    The caller must allocate space for the returned values.

@Description
The function @func{gridInqYvals} returns all values of the Y-axis.

@Result
Upon successful completion @func{gridInqYvals} returns the number of values and
the values are stored in @func{yvals}.
Otherwise, 0 is returned and @func{yvals} is empty.

@EndFunction
*/
SizeType
gridInqYvals(int gridID, double *yvals)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  return gridptr->vtable->inqYVals(gridptr, yvals);
}

SizeType
gridInqYvalsPart(int gridID, int start, SizeType size, double *yvals)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  return gridptr->vtable->inqYValsPart(gridptr, start, size, yvals);
}

SizeType
gridInqYCvals(int gridID, char **ycvals)
{
  grid_t *gridptr = grid_to_pointer(gridID);
#ifndef USE_MPI
  return gridptr->vtable->inqYCvals(gridptr, ycvals);
#else
  return 0;
#endif
}

int
gridInqYIsc(int gridID)
{
  grid_t *gridptr = grid_to_pointer(gridID);
#ifndef USE_MPI
  return gridptr->vtable->inqYIsc(gridptr);
#else
  return 0;
#endif
}

static void
gridDefYValsSerial(grid_t *gridptr, const double *yvals)
{
  int gridtype = gridptr->type;
  size_t size = grid_is_irregular(gridtype) ? gridptr->size : gridptr->y.size;

  if (size == 0) Error("Size undefined for gridID = %d!", gridptr->self);

  if (gridptr->y.vals && CDI_Debug) Warning("Values already defined!");

  gridptr->y.vals = (double *) Realloc(gridptr->y.vals, size * sizeof(double));
  copy_darray(size, yvals, gridptr->y.vals);
}

/*
@Function  gridDefYvals
@Title     Define the values of a Y-axis

@Prototype void gridDefYvals(int gridID, const double *yvals)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate}.
    @Item  yvals    Y-values of the grid.

@Description
The function @func{gridDefYvals} defines all values of the Y-axis.

@EndFunction
*/
void
gridDefYvals(int gridID, const double *yvals)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  gridptr->vtable->defYVals(gridptr, yvals);
  gridMark4Update(gridID);
}

static double
gridInqXValSerial(grid_t *gridptr, SizeType index)
{
  const double xval = gridptr->x.vals ? gridptr->x.vals[index] : 0;
  return xval;
}

double
gridInqXval(int gridID, SizeType index)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  return gridptr->vtable->inqXVal(gridptr, index);
}

static double
gridInqYValSerial(grid_t *gridptr, SizeType index)
{
  const double yval = gridptr->y.vals ? gridptr->y.vals[index] : 0;
  return yval;
}

/*
@Function
@Title

@Prototype
@Parameter
    @Item  Grid identifier

@EndFunction
*/
double
gridInqYval(int gridID, SizeType index)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  return gridptr->vtable->inqYVal(gridptr, index);
}

static double
grid_calc_increment(size_t size, const double *vals)
{
  if (size > 1)
    {
      double inc = (vals[size - 1] - vals[0]) / (size - 1);
      const double abs_inc = fabs(inc);
      for (size_t i = 1; i < size; ++i)
        if (fabs(fabs(vals[i - 1] - vals[i]) - abs_inc) > 0.01 * abs_inc)
          {
            inc = 0.0;
            break;
          }

      return inc;
    }

  return 0.0;
}

static double
grid_calc_increment_in_meter(SizeType size, const double *vals)
{
  if (size > 1)
    {
      const double inc = (vals[size - 1] - vals[0]) / (size - 1);
      return round(fabs(inc));
    }

  return 0.0;
}

static double
gridInqXIncBase(grid_t *gridptr)
{
  if (fabs(gridptr->x.inc) <= 0 && gridptr->x.vals)
    {
      size_t xsize = gridptr->x.size;
      if (xsize > 1)
        {
          const double *xvals = gridptr->vtable->inqXValsPtr(gridptr);
          gridptr->x.inc = grid_calc_increment(xsize, xvals);
        }
    }

  return gridptr->x.inc;
}

double
gridInqXincInMeter(int gridID)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  const double *xvals = gridptr->vtable->inqXValsPtr(gridptr);

  if (fabs(gridptr->x.inc) <= 0 && xvals)
    {
      size_t xsize = gridptr->x.size;
      if (xsize > 1) gridptr->x.inc = grid_calc_increment_in_meter(xsize, xvals);
    }

  return gridptr->x.inc;
}

/*
@Function
@Title

@Prototype
@Parameter
    @Item  Grid identifier

@EndFunction
*/
double
gridInqXinc(int gridID)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  return gridptr->vtable->inqXInc(gridptr);
}

static double
gridInqYIncBase(grid_t *gridptr)
{
  if (fabs(gridptr->y.inc) <= 0 && gridptr->y.vals)
    {
      size_t ysize = gridptr->y.size;
      if (ysize > 1)
        {
          const double *yvals = gridptr->vtable->inqYValsPtr(gridptr);
          gridptr->y.inc = grid_calc_increment(ysize, yvals);
        }
    }

  return gridptr->y.inc;
}

double
gridInqYincInMeter(int gridID)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  const double *yvals = gridptr->vtable->inqYValsPtr(gridptr);

  if (fabs(gridptr->y.inc) <= 0 && yvals)
    {
      size_t ysize = gridptr->y.size;
      if (ysize > 1) gridptr->y.inc = grid_calc_increment_in_meter(ysize, yvals);
    }

  return gridptr->y.inc;
}

/*
@Function
@Title

@Prototype
@Parameter
    @Item  Grid identifier

@EndFunction
*/
double
gridInqYinc(int gridID)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  return gridptr->vtable->inqYInc(gridptr);
}

/*
@Function
@Title

@Prototype
@Parameter
    @Item  Grid identifier

@EndFunction
*/
void
gridInqParamRLL(int gridID, double *xpole, double *ypole, double *angle)
{
  *xpole = 0;
  *ypole = 0;
  *angle = 0;

  static const char projection[] = "rotated_latitude_longitude";
  char name[CDI_MAX_NAME + 1];
  int length = CDI_MAX_NAME;
  cdiInqKeyString(gridID, CDI_GLOBAL, CDI_KEY_GRIDMAP_NAME, name, &length);
  if (name[0] && str_is_equal(name, projection))
    {
      int atttype, attlen;

      int natts, nfound = 0;
      cdiInqNatts(gridID, CDI_GLOBAL, &natts);

      for (int iatt = 0; iatt < natts; ++iatt)
        {
          cdiInqAtt(gridID, CDI_GLOBAL, iatt, name, &atttype, &attlen);
          if (attlen == 1)
            {
              double *attflt;
              // clang-format off
              if      (str_is_equal(name, "grid_north_pole_longitude")) attflt = xpole;
              else if (str_is_equal(name, "grid_north_pole_latitude") ) attflt = ypole;
              else if (str_is_equal(name, "north_pole_grid_longitude")) attflt = angle;
              else continue;
              // clang-format on
              bool valid = cdiInqAttConvertedToFloat(gridID, atttype, name, attlen, attflt);
              if ((nfound += valid) == 3) return;
            }
        }
    }
  else
    Warning("%s mapping parameter missing!", projection);
}

/*
@Function
@Title

@Prototype
@Parameter
    @Item  Grid identifier

@EndFunction
*/
void
gridDefParamRLL(int gridID, double xpole, double ypole, double angle)
{
  cdiDefKeyString(gridID, CDI_XAXIS, CDI_KEY_UNITS, "degrees");
  cdiDefKeyString(gridID, CDI_YAXIS, CDI_KEY_UNITS, "degrees");

  cdiDefKeyString(gridID, CDI_GLOBAL, CDI_KEY_GRIDMAP_VARNAME, "rotated_pole");

  const char *gmapname = "rotated_latitude_longitude";
  cdiDefKeyString(gridID, CDI_GLOBAL, CDI_KEY_GRIDMAP_NAME, gmapname);
  cdiDefAttTxt(gridID, CDI_GLOBAL, "grid_mapping_name", (int) (strlen(gmapname)), gmapname);
  cdiDefAttFlt(gridID, CDI_GLOBAL, "grid_north_pole_longitude", CDI_DATATYPE_FLT64, 1, &xpole);
  cdiDefAttFlt(gridID, CDI_GLOBAL, "grid_north_pole_latitude", CDI_DATATYPE_FLT64, 1, &ypole);
  if (IS_NOT_EQUAL(angle, 0)) cdiDefAttFlt(gridID, CDI_GLOBAL, "north_pole_grid_longitude", CDI_DATATYPE_FLT64, 1, &angle);

  grid_t *gridptr = grid_to_pointer(gridID);
  gridptr->projtype = CDI_PROJ_RLL;

  gridVerifyProj(gridID);
}

/*
@Function
@Title

@Prototype
@Parameter
    @Item  Grid identifier

@EndFunction
*/
void
gridInqParamGME(int gridID, int *nd, int *ni, int *ni2, int *ni3)
{
  grid_t *gridptr = grid_to_pointer(gridID);

  *nd = gridptr->gme.nd;
  *ni = gridptr->gme.ni;
  *ni2 = gridptr->gme.ni2;
  *ni3 = gridptr->gme.ni3;
}

/*
@Function
@Title

@Prototype
@Parameter
    @Item  Grid identifier

@EndFunction
*/
void
gridDefParamGME(int gridID, int nd, int ni, int ni2, int ni3)
{
  grid_t *gridptr = grid_to_pointer(gridID);

  if (gridptr->gme.nd != nd)
    {
      gridptr->gme.nd = nd;
      gridptr->gme.ni = ni;
      gridptr->gme.ni2 = ni2;
      gridptr->gme.ni3 = ni3;
      gridMark4Update(gridID);
    }
}

/*
@Function
@Title

@Prototype
@Parameter
    @Item  Grid identifier

@EndFunction
*/
void
gridChangeType(int gridID, int gridtype)
{
  grid_t *gridptr = grid_to_pointer(gridID);

  if (CDI_Debug) Message("Changed grid type from %s to %s", gridNamePtr(gridptr->type), gridNamePtr(gridtype));

  if (gridptr->type != gridtype)
    {
      gridptr->type = gridtype;
      gridMark4Update(gridID);
    }
}

static void
grid_check_cyclic(grid_t *gridptr)
{
  gridptr->isCyclic = 0;
  enum
  {
    numVertices = 4
  };
  size_t xsize = gridptr->x.size, ysize = gridptr->y.size;
  const double *xvals = gridptr->vtable->inqXValsPtr(gridptr), *yvals = gridptr->vtable->inqYValsPtr(gridptr),
               (*xbounds)[numVertices] = (const double(*)[numVertices]) gridptr->vtable->inqXBoundsPtr(gridptr);

  if (gridptr->type == GRID_GAUSSIAN || gridptr->type == GRID_LONLAT)
    {
      if (xvals && xsize > 1)
        {
          double xval1 = xvals[0];
          double xval2 = xvals[1];
          double xvaln = xvals[xsize - 1];
          if (xval2 < xval1) xval2 += 360;
          if (xvaln < xval1) xvaln += 360;

          if (IS_NOT_EQUAL(xval1, xvaln))
            {
              double xinc = xval2 - xval1;
              if (IS_EQUAL(xinc, 0)) xinc = (xvaln - xval1) / (xsize - 1);

              const double x0 = xvaln + xinc - 360;

              if (fabs(x0 - xval1) < 0.01 * xinc) gridptr->isCyclic = 1;
            }
        }
    }
  else if (gridptr->type == GRID_CURVILINEAR)
    {
      bool lcheck = true;
      if (yvals && xvals)
        {
          if ((fabs(yvals[0] - yvals[xsize - 1]) > fabs(yvals[0] - yvals[xsize * ysize - xsize]))
              && (fabs(yvals[xsize * ysize - xsize] - yvals[xsize * ysize - 1])
                  > fabs(yvals[xsize - 1] - yvals[xsize * ysize - 1])))
            lcheck = false;
        }
      else
        lcheck = false;

      if (lcheck && xvals && xsize > 1)
        {
          size_t nc = 0;
          for (size_t j = 0; j < ysize; ++j)
            {
              size_t i1 = j * xsize, i2 = j * xsize + 1, in = j * xsize + (xsize - 1);
              double val1 = xvals[i1], val2 = xvals[i2], valn = xvals[in];
              double xinc = fabs(val2 - val1);

              if (val1 < 1 && valn > 300) val1 += 360;
              if (valn < 1 && val1 > 300) valn += 360;
              if (val1 < -179 && valn > 120) val1 += 360;
              if (valn < -179 && val1 > 120) valn += 360;
              if (fabs(valn - val1) > 180) val1 += 360;

              double x0 = valn + copysign(xinc, val1 - valn);

              nc += fabs(x0 - val1) < 0.5 * xinc;
            }
          gridptr->isCyclic = nc > ysize / 2;
        }

      if (lcheck && xbounds && xsize > 1)
        {
          bool isCyclic = true;
          for (size_t j = 0; j < ysize; ++j)
            {
              size_t i1 = j * xsize, i2 = j * xsize + (xsize - 1);
              for (size_t k1 = 0; k1 < numVertices; ++k1)
                {
                  double val1 = xbounds[i1][k1];
                  for (size_t k2 = 0; k2 < numVertices; ++k2)
                    {
                      double val2 = xbounds[i2][k2];

                      if (val1 < 1 && val2 > 300) val1 += 360;
                      if (val2 < 1 && val1 > 300) val2 += 360;
                      if (val1 < -179 && val2 > 120) val1 += 360;
                      if (val2 < -179 && val1 > 120) val2 += 360;
                      if (fabs(val2 - val1) > 180) val1 += 360;

                      if (fabs(val1 - val2) < 0.001) goto foundCloseVertices;
                    }
                }
              // all vertices more than 0.001 degrees apart
              isCyclic = false;
              break;
            foundCloseVertices:;
            }
          gridptr->isCyclic = isCyclic;
        }
    }
}

int
gridIsCircular(int gridID)
{
  grid_t *gridptr = grid_to_pointer(gridID);

  if (gridptr->isCyclic == CDI_UNDEFID) grid_check_cyclic(gridptr);

  return gridptr->isCyclic;
}

static bool
compareXYvals(grid_t *gridRef, grid_t *gridTest)
{
  bool differ = false;
  int gridtype = gridTest->type;

  size_t xsizeTest = grid_is_irregular(gridtype) ? gridTest->size : gridTest->x.size;
  size_t xsizeRef = (size_t) gridRef->vtable->inqXVals(gridRef, NULL);
  if (xsizeTest != xsizeRef) return true;

  if (xsizeTest > 0)
    {
      const double *xvalsRef = gridRef->vtable->inqXValsPtr(gridRef);
      const double *xvalsTest = gridTest->vtable->inqXValsPtr(gridTest);
      if (!xvalsTest) return true;

      for (size_t i = 0; i < xsizeTest; ++i)
        if (fabs(xvalsTest[i] - xvalsRef[i]) > 1.e-10) return true;
    }

  size_t ysizeTest = grid_is_irregular(gridtype) ? gridTest->size : gridTest->y.size;
  size_t ysizeRef = (size_t) gridRef->vtable->inqYVals(gridRef, NULL);
  if (ysizeTest != ysizeRef) return true;

  if (ysizeTest > 0)
    {
      const double *yvalsRef = gridRef->vtable->inqYValsPtr(gridRef);
      const double *yvalsTest = gridTest->vtable->inqYValsPtr(gridTest);
      if (!yvalsTest) return true;

      for (size_t i = 0; i < ysizeTest; ++i)
        if (fabs(yvalsTest[i] - yvalsRef[i]) > 1.e-10) return true;
    }

  return differ;
}

static bool
compareXYvals2(grid_t *gridRef, grid_t *gridTest)
{
  size_t gridsize = gridTest->size;
  bool differ = ((gridTest->x.vals == NULL) ^ (gridRef->x.vals == NULL)) || ((gridTest->y.vals == NULL) ^ (gridRef->y.vals == NULL))
                || ((gridTest->x.bounds == NULL) ^ (gridRef->x.bounds == NULL))
                || ((gridTest->y.bounds == NULL) ^ (gridRef->y.bounds == NULL));

  typedef double (*inqVal)(grid_t * grid, SizeType index);
  inqVal inqXValRef = gridRef->vtable->inqXVal, inqYValRef = gridRef->vtable->inqYVal, inqXValTest = gridTest->vtable->inqXVal,
         inqYValTest = gridTest->vtable->inqYVal;

  if (!differ && gridTest->x.vals)
    differ = fabs(inqXValTest(gridTest, 0) - inqXValRef(gridRef, 0)) > 1.e-9
             || fabs(inqXValTest(gridTest, gridsize - 1) - inqXValRef(gridRef, gridsize - 1)) > 1.e-9;

  if (!differ && gridTest->y.vals)
    differ = fabs(inqYValTest(gridTest, 0) - inqYValRef(gridRef, 0)) > 1.e-9
             || fabs(inqYValTest(gridTest, gridsize - 1) - inqYValRef(gridRef, gridsize - 1)) > 1.e-9;

  return differ;
}

static bool
compare_bounds(const grid_t *grid, const grid_t *gridRef)
{
  bool differ = false;

  if ((grid->x.bounds && !gridRef->x.bounds) || (!grid->x.bounds && gridRef->x.bounds) || (grid->y.bounds && !gridRef->y.bounds)
      || (!grid->y.bounds && gridRef->y.bounds))
    differ = true;

  return differ;
}

static bool
compare_lonlat(int gridID, const grid_t *grid, const grid_t *gridRef)
{
  bool differ = false;
  /*
    printf("gridID      %d\n", gridID);
    printf("grid.xdef   %d\n", grid->x.flag);
    printf("grid.ydef   %d\n", grid->y.flag);
    printf("grid.xsize  %zu\n", grid->x.size);
    printf("grid.ysize  %zu\n", grid->y.size);
    printf("grid.xfirst %f\n", grid->x.first);
    printf("grid.yfirst %f\n", grid->y.first);
    printf("grid.xfirst %f\n", gridInqXval(gridID, 0));
    printf("grid.yfirst %f\n", gridInqYval(gridID, 0));
    printf("grid.xinc   %f\n", grid->x.inc);
    printf("grid.yinc   %f\n", grid->y.inc);
    printf("grid.xinc   %f\n", gridInqXinc(gridID));
    printf("grid.yinc   %f\n", gridInqYinc(gridID));
  */
  if (grid->x.size == gridRef->x.size && grid->y.size == gridRef->y.size)
    {
      if (grid->x.flag == 2 && grid->y.flag == 2)
        {
          if (!(IS_EQUAL(grid->x.first, 0) && IS_EQUAL(grid->x.last, 0) && IS_EQUAL(grid->x.inc, 0))
              && !(IS_EQUAL(grid->y.first, 0) && IS_EQUAL(grid->y.last, 0) && IS_EQUAL(grid->y.inc, 0))
              && IS_NOT_EQUAL(grid->x.first, grid->x.last) && IS_NOT_EQUAL(grid->y.first, grid->y.last))
            {
              if (IS_NOT_EQUAL(grid->x.first, gridInqXval(gridID, 0)) || IS_NOT_EQUAL(grid->y.first, gridInqYval(gridID, 0)))
                {
                  differ = true;
                }
              if (!differ && fabs(grid->x.inc) > 0 && fabs(fabs(grid->x.inc) - fabs(gridRef->x.inc)) > fabs(grid->x.inc / 1000))
                {
                  differ = true;
                }
              if (!differ && fabs(grid->y.inc) > 0 && fabs(fabs(grid->y.inc) - fabs(gridRef->y.inc)) > fabs(grid->y.inc / 1000))
                {
                  differ = true;
                }
            }
        }
      else if (grid->x.vals && grid->y.vals)
        differ = gridRef->vtable->compareXYFull((grid_t *) gridRef, (grid_t *) grid);

      if (!differ) differ = compare_bounds(grid, gridRef);
    }
  else
    differ = true;

  return differ;
}

static bool
compare_projection(int gridID, const grid_t *grid, const grid_t *gridRef)
{
  bool differ = compare_lonlat(gridID, grid, gridRef);

  if (!differ)
    {
      // printf(">%s< >%s<\n", cdiInqVarKeyString(&grid->keys, CDI_KEY_GRIDMAP_VARNAME), cdiInqVarKeyString(&gridRef->keys,
      // CDI_KEY_GRIDMAP_VARNAME)); printf(">%s< >%s<\n", cdiInqVarKeyString(&grid->keys, CDI_KEY_GRIDMAP_NAME),
      // cdiInqVarKeyString(&gridRef->keys, CDI_KEY_GRIDMAP_NAME));
      // if (!str_is_equal(cdiInqVarKeyString(&grid->keys, CDI_KEY_GRIDMAP_VARNAME), cdiInqVarKeyString(&gridRef->keys,
      // CDI_KEY_GRIDMAP_VARNAME))) return true; if (!str_is_equal(cdiInqVarKeyString(&grid->keys, CDI_KEY_GRIDMAP_NAME),
      // cdiInqVarKeyString(&gridRef->keys, CDI_KEY_GRIDMAP_NAME))) return true;
    }

  return differ;
}

static bool
compare_generic(const grid_t *grid, const grid_t *gridRef)
{
  bool differ = false;

  if (grid->x.size == gridRef->x.size && grid->y.size == gridRef->y.size)
    {
      if (grid->x.flag == 1 && grid->y.flag == 1 && grid->x.vals && grid->y.vals)
        differ = gridRef->vtable->compareXYFull((grid_t *) gridRef, (grid_t *) grid);
    }
  else if ((grid->y.size == 0 || grid->y.size == 1) && grid->x.size == gridRef->x.size * gridRef->y.size)
    {
    }
  else
    differ = true;

  return differ;
}

static bool
compare_gaussian(int gridID, const grid_t *grid, const grid_t *gridRef)
{
  const double cmp_eps = 0.0015;
  bool differ = false;

  if (grid->x.size == gridRef->x.size && grid->y.size == gridRef->y.size)
    {
      if (grid->x.flag == 2 && grid->y.flag == 2)
        {
          if (!(IS_EQUAL(grid->x.first, 0) && IS_EQUAL(grid->x.last, 0) && IS_EQUAL(grid->x.inc, 0))
              && !(IS_EQUAL(grid->y.first, 0) && IS_EQUAL(grid->y.last, 0)))
            if (fabs(grid->x.first - gridInqXval(gridID, 0)) > cmp_eps || fabs(grid->y.first - gridInqYval(gridID, 0)) > cmp_eps
                || (fabs(grid->x.inc) > 0 && fabs(fabs(grid->x.inc) - fabs(gridRef->x.inc)) > fabs(grid->x.inc / 1000)))
              {
                differ = true;
              }
        }
      else if (grid->x.vals && grid->y.vals)
        differ = gridRef->vtable->compareXYFull((grid_t *) gridRef, (grid_t *) grid);

      if (!differ) differ = compare_bounds(grid, gridRef);
    }
  else
    differ = true;

  return differ;
}

static bool
compare_curvilinear(const grid_t *grid, const grid_t *gridRef)
{
  bool differ = false;

  /*
    printf("gridID      %d\n", gridID);
    printf("grid.xsize  %d\n", grid->x.size);
    printf("grid.ysize  %d\n", grid->y.size);
    printf("grid.xfirst %f\n", grid->x.vals[0]);
    printf("grid.yfirst %f\n", grid->y.vals[0]);
    printf("grid xfirst %f\n", gridInqXval(gridID, 0));
    printf("grid yfirst %f\n", gridInqYval(gridID, 0));
    printf("grid.xlast  %f\n", grid->x.vals[grid->size-1]);
    printf("grid.ylast  %f\n", grid->y.vals[grid->size-1]);
    printf("grid xlast  %f\n", gridInqXval(gridID, grid->size-1));
    printf("grid ylast  %f\n", gridInqYval(gridID, grid->size-1));
    printf("grid.nv     %d\n", grid->nvertex);
    printf("grid nv     %d\n", gridInqNvertex(gridID));
  */
  if (grid->x.size == gridRef->x.size && grid->y.size == gridRef->y.size)
    differ = gridRef->vtable->compareXYAO((grid_t *) gridRef, (grid_t *) grid);

  return differ;
}

static bool
compare_unstructured(const grid_t *grid, const grid_t *gridRef, bool compareCoord)
{
  bool differ = false;

  unsigned char uuid1[CDI_UUID_SIZE] = { 0 };
  unsigned char uuid2[CDI_UUID_SIZE] = { 0 };
  int length = CDI_UUID_SIZE;
  cdiInqVarKeyBytes(&gridRef->keys, CDI_KEY_UUID, uuid1, &length);
  length = CDI_UUID_SIZE;
  cdiInqVarKeyBytes(&grid->keys, CDI_KEY_UUID, uuid2, &length);
  differ = ((!cdiUUIDIsNull(uuid1) || !cdiUUIDIsNull(uuid2)) && memcmp(uuid1, uuid2, CDI_UUID_SIZE));
  if (!differ)
    {
      int numberA = cdiInqVarKeyInt(&grid->keys, CDI_KEY_NUMBEROFGRIDUSED);
      int numberB = cdiInqVarKeyInt(&gridRef->keys, CDI_KEY_NUMBEROFGRIDUSED);
      int positionA = cdiInqVarKeyInt(&grid->keys, CDI_KEY_NUMBEROFGRIDINREFERENCE);
      int positionB = cdiInqVarKeyInt(&gridRef->keys, CDI_KEY_NUMBEROFGRIDINREFERENCE);
      if (compareCoord)
        {
          differ = (grid->nvertex != gridRef->nvertex || (numberA > 0 && positionA != positionB)
                    || gridRef->vtable->compareXYFull((grid_t *) gridRef, (grid_t *) grid));
        }
      else
        {
          if (((grid->x.vals == NULL) ^ (gridRef->x.vals == NULL)) && ((grid->y.vals == NULL) ^ (gridRef->y.vals == NULL)))
            {
              int nvertexA = grid->nvertex, nvertexB = gridRef->nvertex;
              differ = (nvertexA && nvertexB && (nvertexA != nvertexB))
                       || ((numberA && numberB && (numberA != numberB)) || (numberA && numberB && positionA != positionB));
            }
          else
            {
              differ = (grid->nvertex != gridRef->nvertex || numberA != numberB || (numberA > 0 && positionA != positionB)
                        || gridRef->vtable->compareXYAO((grid_t *) gridRef, (grid_t *) grid));
            }
        }
    }

  return differ;
}

static bool
gridCompare(int gridID, const grid_t *grid, bool compareCoord)
{
  bool differ = true;
  const grid_t *gridRef = grid_to_pointer(gridID);

  if (grid->type == gridRef->type || grid->type == GRID_GENERIC)
    {
      if (grid->size == gridRef->size)
        {
          differ = false;
          if (grid->type == GRID_LONLAT)
            {
              differ = compare_lonlat(gridID, grid, gridRef);
            }
          else if (grid->type == GRID_PROJECTION)
            {
              differ = compare_projection(gridID, grid, gridRef);
            }
          else if (grid->type == GRID_GENERIC)
            {
              differ = compare_generic(grid, gridRef);
            }
          else if (grid->type == GRID_GAUSSIAN)
            {
              differ = compare_gaussian(gridID, grid, gridRef);
            }
          else if (grid->type == GRID_CURVILINEAR)
            {
              differ = compare_curvilinear(grid, gridRef);
            }
          else if (grid->type == GRID_UNSTRUCTURED)
            {
              differ = compare_unstructured(grid, gridRef, compareCoord);
            }
        }
    }

  int scanningModeA = cdiInqVarKeyInt(&grid->keys, CDI_KEY_SCANNINGMODE);
  int scanningModeB = cdiInqVarKeyInt(&gridRef->keys, CDI_KEY_SCANNINGMODE);
  if (scanningModeA != scanningModeB)
    {
      // often grid definition may differ in UV-relativeToGrid
      differ = true;
#ifdef HIRLAM_EXTENSIONS
      if (cdiDebugExt >= 200)
        printf("gridCompare(gridID=%d): Differs: scanningModeA [%d] != scanningModeB(gridID) [%d]\n", gridID, scanningModeA,
               scanningModeB);
#endif  // HIRLAM_EXTENSIONS
    }

  return differ;
}

int
cmp_key_int(const cdi_keys_t *keysp1, const cdi_keys_t *keysp2, int key)
{
  int v1 = cdiInqVarKeyInt(keysp1, key);
  int v2 = cdiInqVarKeyInt(keysp2, key);
  return (v1 != v2);
}

int
gridCompareP(void *gridptr1, void *gridptr2)
{
  grid_t *g1 = (grid_t *) gridptr1;
  grid_t *g2 = (grid_t *) gridptr2;
  enum
  {
    equal = 0,
    differ = -1
  };

  xassert(g1);
  xassert(g2);

  if (cdiInqVarKeyInt(&g1->keys, CDI_KEY_DATATYPE) != cdiInqVarKeyInt(&g2->keys, CDI_KEY_DATATYPE)) return differ;
  if (g1->type != g2->type) return differ;
  if (g1->isCyclic != g2->isCyclic) return differ;
  if (g1->x.flag != g2->x.flag) return differ;
  if (g1->y.flag != g2->y.flag) return differ;
  if (g1->gme.nd != g2->gme.nd) return differ;
  if (g1->gme.ni != g2->gme.ni) return differ;
  if (g1->gme.ni2 != g2->gme.ni2) return differ;
  if (g1->gme.ni3 != g2->gme.ni3) return differ;
  if (cmp_key_int(&g1->keys, &g2->keys, CDI_KEY_NUMBEROFGRIDUSED)) return differ;
  if (cmp_key_int(&g1->keys, &g2->keys, CDI_KEY_NUMBEROFGRIDINREFERENCE)) return differ;
  if (g1->trunc != g2->trunc) return differ;
  if (g1->nvertex != g2->nvertex) return differ;
  if (g1->reducedPointsSize != g2->reducedPointsSize) return differ;
  if (g1->size != g2->size) return differ;
  if (g1->x.size != g2->x.size) return differ;
  if (g1->y.size != g2->y.size) return differ;
  if (g1->lcomplex != g2->lcomplex) return differ;

  if (IS_NOT_EQUAL(g1->x.first, g2->x.first)) return differ;
  if (IS_NOT_EQUAL(g1->y.first, g2->y.first)) return differ;
  if (IS_NOT_EQUAL(g1->x.last, g2->x.last)) return differ;
  if (IS_NOT_EQUAL(g1->y.last, g2->y.last)) return differ;
  if (IS_NOT_EQUAL(g1->x.inc, g2->x.inc)) return differ;
  if (IS_NOT_EQUAL(g1->y.inc, g2->y.inc)) return differ;
  if (cmp_key_int(&g1->keys, &g2->keys, CDI_KEY_SCANNINGMODE)) return differ;

  bool isIrregular = grid_is_irregular(g1->type);
  {
    const double *restrict g1_xvals = g1->vtable->inqXValsPtr(g1), *restrict g2_xvals = g2->vtable->inqXValsPtr(g2);
    if ((g1_xvals != NULL) ^ (g2_xvals != NULL)) return differ;
    if (g1_xvals)
      {
        size_t size = isIrregular ? g1->size : g1->x.size;
        xassert(size);
        for (size_t i = 0; i < size; i++)
          if (IS_NOT_EQUAL(g1_xvals[i], g2_xvals[i])) return differ;
      }
  }

  {
    const double *restrict g1_yvals = g1->vtable->inqYValsPtr(g1), *restrict g2_yvals = g2->vtable->inqYValsPtr(g2);
    if ((g1_yvals != NULL) ^ (g2_yvals != NULL)) return differ;
    if (g1_yvals)
      {
        size_t size = isIrregular ? g1->size : g1->y.size;
        xassert(size);
        for (size_t i = 0; i < size; i++)
          if (IS_NOT_EQUAL(g1_yvals[i], g2_yvals[i])) return differ;
      }
  }

  {
    const double *restrict g1_area = g1->vtable->inqAreaPtr(g1), *restrict g2_area = g2->vtable->inqAreaPtr(g2);
    if ((g1_area != NULL) ^ (g2_area != NULL)) return differ;
    if (g1_area)
      {
        size_t size = g1->size;
        xassert(size);

        for (size_t i = 0; i < size; i++)
          if (IS_NOT_EQUAL(g1_area[i], g2_area[i])) return differ;
      }
  }

  {
    const double *restrict g1_xbounds = g1->vtable->inqXBoundsPtr(g1), *restrict g2_xbounds = g2->vtable->inqXBoundsPtr(g2);
    if ((g1_xbounds != NULL) ^ (g2_xbounds != NULL)) return differ;
    if (g1_xbounds)
      {
        xassert(g1->nvertex);
        size_t size = g1->nvertex * (isIrregular ? g1->size : g1->x.size);
        xassert(size);

        for (size_t i = 0; i < size; i++)
          if (IS_NOT_EQUAL(g1_xbounds[i], g2_xbounds[i])) return differ;
      }
  }

  {
    const double *restrict g1_ybounds = g1->vtable->inqYBoundsPtr(g1), *restrict g2_ybounds = g2->vtable->inqYBoundsPtr(g2);
    if ((g1_ybounds != NULL) ^ (g2_ybounds != NULL)) return differ;
    if (g1_ybounds)
      {
        xassert(g1->nvertex);
        size_t size = g1->nvertex * (isIrregular ? g1->size : g1->y.size);
        xassert(size);

        for (size_t i = 0; i < size; i++)
          if (IS_NOT_EQUAL(g1_ybounds[i], g2_ybounds[i])) return differ;
      }
  }

  if (!str_is_equal(cdiInqVarKeyString(&g1->x.keys, CDI_KEY_NAME), cdiInqVarKeyString(&g2->x.keys, CDI_KEY_NAME))) return differ;
  if (!str_is_equal(cdiInqVarKeyString(&g1->y.keys, CDI_KEY_NAME), cdiInqVarKeyString(&g2->y.keys, CDI_KEY_NAME))) return differ;
  if (!str_is_equal(cdiInqVarKeyString(&g1->x.keys, CDI_KEY_LONGNAME), cdiInqVarKeyString(&g2->x.keys, CDI_KEY_LONGNAME)))
    return differ;
  if (!str_is_equal(cdiInqVarKeyString(&g1->y.keys, CDI_KEY_LONGNAME), cdiInqVarKeyString(&g2->y.keys, CDI_KEY_LONGNAME)))
    return differ;
  if (!str_is_equal(cdiInqVarKeyString(&g1->x.keys, CDI_KEY_UNITS), cdiInqVarKeyString(&g2->x.keys, CDI_KEY_UNITS))) return differ;
  if (!str_is_equal(cdiInqVarKeyString(&g1->y.keys, CDI_KEY_UNITS), cdiInqVarKeyString(&g2->y.keys, CDI_KEY_UNITS))) return differ;
  if (!str_is_equal(cdiInqVarKeyString(&g1->x.keys, CDI_KEY_STDNAME), cdiInqVarKeyString(&g2->x.keys, CDI_KEY_STDNAME)))
    return differ;
  if (!str_is_equal(cdiInqVarKeyString(&g1->y.keys, CDI_KEY_STDNAME), cdiInqVarKeyString(&g2->y.keys, CDI_KEY_STDNAME)))
    return differ;

  if (!str_is_equal(cdiInqVarKeyString(&g1->y.keys, CDI_KEY_REFERENCEURI), cdiInqVarKeyString(&g2->y.keys, CDI_KEY_REFERENCEURI)))
    return differ;

  if (g1->mask)
    {
      xassert(g1->size);
      if (!g2->mask) return differ;
      if (memcmp(g1->mask, g2->mask, g1->size * sizeof(mask_t))) return differ;
    }
  else if (g2->mask)
    return differ;

  if (g1->mask_gme)
    {
      xassert(g1->size);
      if (!g2->mask_gme) return differ;
      if (memcmp(g1->mask_gme, g2->mask_gme, g1->size * sizeof(mask_t))) return differ;
    }
  else if (g2->mask_gme)
    return differ;

  unsigned char uuid1[CDI_UUID_SIZE] = { 0 };
  unsigned char uuid2[CDI_UUID_SIZE] = { 0 };
  int length = CDI_UUID_SIZE;
  cdiInqVarKeyBytes(&g1->keys, CDI_KEY_UUID, uuid1, &length);
  length = CDI_UUID_SIZE;
  cdiInqVarKeyBytes(&g2->keys, CDI_KEY_UUID, uuid2, &length);
  if (memcmp(uuid1, uuid2, CDI_UUID_SIZE)) return differ;

  return equal;
}

static void
grid_complete(grid_t *grid)
{
  int gridID = grid->self;

  if (grid->datatype != CDI_UNDEFID) cdiDefKeyInt(gridID, CDI_GLOBAL, CDI_KEY_DATATYPE, grid->datatype);

  int gridtype = grid->type;
  switch (gridtype)
    {
    case GRID_LONLAT:
    case GRID_GAUSSIAN:
    case GRID_UNSTRUCTURED:
    case GRID_CURVILINEAR:
    case GRID_GENERIC:
    case GRID_PROJECTION:
    case GRID_CHARXY:
      {
        if (grid->x.size > 0) gridDefXsize(gridID, grid->x.size);
        if (grid->y.size > 0) gridDefYsize(gridID, grid->y.size);

        if (gridtype == GRID_GAUSSIAN) gridDefNP(gridID, grid->np);

        if (grid->nvertex > 0) gridDefNvertex(gridID, grid->nvertex);

        if (grid->x.flag == 2)
          {
            assert(gridtype != GRID_UNSTRUCTURED && gridtype != GRID_CURVILINEAR);
            double *xvals = (double *) Malloc(grid->x.size * sizeof(double));
            gridGenXvals(grid->x.size, grid->x.first, grid->x.last, grid->x.inc, xvals);
            grid->x.vals = xvals;
            // gridDefXinc(gridID, grid->x.inc);
          }

        if (grid->y.flag == 2)
          {
            assert(gridtype != GRID_UNSTRUCTURED && gridtype != GRID_CURVILINEAR);
            double *yvals = (double *) Malloc(grid->y.size * sizeof(double));
            gridGenYvals(gridtype, grid->y.size, grid->y.first, grid->y.last, grid->y.inc, yvals);
            grid->y.vals = yvals;
            // gridDefYinc(gridID, grid->y.inc);
          }

        if (grid->projtype == CDI_PROJ_RLL)
          {
            const char *name = cdiInqVarKeyString(&grid->x.keys, CDI_KEY_NAME);
            if (name[0] == 0 || name[0] == 'x') cdiDefKeyString(gridID, CDI_XAXIS, CDI_KEY_NAME, "rlon");
            name = cdiInqVarKeyString(&grid->y.keys, CDI_KEY_NAME);
            if (name[0] == 0 || name[0] == 'y') cdiDefKeyString(gridID, CDI_YAXIS, CDI_KEY_NAME, "rlat");
            name = cdiInqVarKeyString(&grid->x.keys, CDI_KEY_LONGNAME);
            if (name[0] == 0) cdiDefKeyString(gridID, CDI_XAXIS, CDI_KEY_LONGNAME, "longitude in rotated pole grid");
            name = cdiInqVarKeyString(&grid->y.keys, CDI_KEY_LONGNAME);
            if (name[0] == 0) cdiDefKeyString(gridID, CDI_YAXIS, CDI_KEY_LONGNAME, "latitude in rotated pole grid");
            name = cdiInqVarKeyString(&grid->x.keys, CDI_KEY_UNITS);
            if (name[0] == 0) cdiDefKeyString(gridID, CDI_XAXIS, CDI_KEY_UNITS, "degrees");
            name = cdiInqVarKeyString(&grid->y.keys, CDI_KEY_UNITS);
            if (name[0] == 0) cdiDefKeyString(gridID, CDI_YAXIS, CDI_KEY_UNITS, "degrees");
            cdiDefKeyString(gridID, CDI_XAXIS, CDI_KEY_STDNAME, xystdname_tab[grid_xystdname_grid_latlon][0]);
            cdiDefKeyString(gridID, CDI_YAXIS, CDI_KEY_STDNAME, xystdname_tab[grid_xystdname_grid_latlon][1]);
          }

        if (gridtype == GRID_UNSTRUCTURED)
          {
            int number = cdiInqVarKeyInt(&grid->keys, CDI_KEY_NUMBEROFGRIDUSED);
            if (number > 0)
              {
                cdiDefKeyInt(gridID, CDI_GLOBAL, CDI_KEY_NUMBEROFGRIDUSED, number);
                int position = cdiInqVarKeyInt(&grid->keys, CDI_KEY_NUMBEROFGRIDINREFERENCE);
                if (position > 0) cdiDefKeyInt(gridID, CDI_GLOBAL, CDI_KEY_NUMBEROFGRIDINREFERENCE, position);
              }
          }

        break;
      }
    case GRID_GAUSSIAN_REDUCED:
      {
        gridDefNP(gridID, grid->np);
        gridDefYsize(gridID, grid->y.size);
        if (grid->x.flag == 2)
          {
            double xvals[2] = { grid->x.first, grid->x.last };
            gridDefXsize(gridID, 2);
            gridDefXvals(gridID, xvals);
          }

        if (grid->y.flag == 2)
          {
            double *yvals = (double *) Malloc(grid->y.size * sizeof(double));
            gridGenYvals(gridtype, grid->y.size, grid->y.first, grid->y.last, grid->y.inc, yvals);
            grid->y.vals = yvals;
            // gridDefYinc(gridID, grid->y.inc);
          }
        break;
      }
    case GRID_SPECTRAL:
      {
        gridDefTrunc(gridID, grid->trunc);
        if (grid->lcomplex) gridDefComplexPacking(gridID, 1);
        break;
      }
    case GRID_FOURIER:
      {
        gridDefTrunc(gridID, grid->trunc);
        break;
      }
    case GRID_GME:
      {
        gridDefParamGME(gridID, grid->gme.nd, grid->gme.ni, grid->gme.ni2, grid->gme.ni3);
        break;
      }
      /*
    case GRID_GENERIC:
      {
        if ( grid->x.size > 0 && grid->y.size > 0 )
          {
            gridDefXsize(gridID, grid->x.size);
            gridDefYsize(gridID, grid->y.size);
            if ( grid->x.vals ) gridDefXvals(gridID, grid->x.vals);
            if ( grid->y.vals ) gridDefYvals(gridID, grid->y.vals);
          }
        break;
      }
      */
    case GRID_TRAJECTORY:
      {
        gridDefXsize(gridID, 1);
        gridDefYsize(gridID, 1);
        break;
      }
    default:
      {
        Error("Gridtype %s unsupported!", gridNamePtr(gridtype));
        break;
      }
    }
}

// Used only in iterator_grib.c
int
gridGenerate(const grid_t *grid)
{
  int gridtype = grid->type;
  int gridID = gridCreate(gridtype, grid->size);
  grid_t *restrict gridptr = grid_to_pointer(gridID);
  cdiCopyVarKey(&grid->keys, CDI_KEY_DATATYPE, &gridptr->keys);
  gridptr->x.size = grid->x.size;
  gridptr->y.size = grid->y.size;
  gridptr->np = grid->np;
  gridptr->nvertex = grid->nvertex;
  gridptr->x.flag = grid->x.flag;
  int valdef_group1 = 0;
  static const int valdef_group1_tab[]
      = { GRID_LONLAT, GRID_GAUSSIAN, GRID_UNSTRUCTURED, GRID_CURVILINEAR, GRID_GENERIC, GRID_PROJECTION };
  for (size_t i = 0; i < sizeof(valdef_group1_tab) / sizeof(valdef_group1_tab[0]); ++i)
    valdef_group1 |= (gridtype == valdef_group1_tab[i]);
  if (valdef_group1 && grid->x.flag == 1)
    {
      gridDefXvals(gridID, grid->x.vals);
      if (grid->x.bounds) gridDefXbounds(gridID, grid->x.bounds);
    }
  gridptr->x.first = grid->x.first;
  gridptr->x.last = grid->x.last;
  gridptr->x.inc = grid->x.inc;
  gridptr->y.flag = grid->y.flag;
  if ((valdef_group1 || gridtype == GRID_GAUSSIAN_REDUCED) && grid->y.flag == 1)
    {
      gridDefYvals(gridID, grid->y.vals);
      if (grid->y.bounds) gridDefYbounds(gridID, grid->y.bounds);
    }
  gridptr->y.first = grid->y.first;
  gridptr->y.last = grid->y.last;
  gridptr->y.inc = grid->y.inc;
  if (valdef_group1 && grid->area) gridDefArea(gridID, grid->area);

  cdiCopyVarKey(&grid->keys, CDI_KEY_NUMBEROFGRIDUSED, &gridptr->keys);
  cdiCopyVarKey(&grid->keys, CDI_KEY_NUMBEROFGRIDINREFERENCE, &gridptr->keys);
  cdiCopyVarKey(&grid->keys, CDI_KEY_REFERENCEURI, &gridptr->keys);

  cdiCopyVarKey(&grid->keys, CDI_KEY_SCANNINGMODE, &gridptr->keys);

  if (gridtype == GRID_PROJECTION) gridptr->name = strdup(grid->name);
  if (gridtype == GRID_GAUSSIAN_REDUCED) gridDefReducedPoints(gridID, grid->y.size, grid->reducedPoints);
  gridptr->trunc = grid->trunc;
  gridptr->lcomplex = grid->lcomplex;
  gridptr->gme.nd = grid->gme.nd;
  gridptr->gme.ni = grid->gme.ni;
  gridptr->gme.ni2 = grid->gme.ni2;
  gridptr->gme.ni3 = grid->gme.ni3;

  grid_complete(gridptr);

  cdiCopyVarKey(&grid->keys, CDI_KEY_UUID, &gridptr->keys);

  return gridID;
}

static void
grid_copy_base_array_fields(grid_t *gridptrOrig, grid_t *gridptrDup)
{
  size_t reducedPointsSize = (SizeType) gridptrOrig->reducedPointsSize;
  size_t gridsize = gridptrOrig->size;
  int gridtype = gridptrOrig->type;
  bool isIrregular = grid_is_irregular(gridtype);
  if (reducedPointsSize)
    {
      gridptrDup->reducedPoints = (int *) Malloc(reducedPointsSize * sizeof(int));
      memcpy(gridptrDup->reducedPoints, gridptrOrig->reducedPoints, reducedPointsSize * sizeof(int));
    }

  if (gridptrOrig->x.vals != NULL)
    {
      size_t size = isIrregular ? gridsize : gridptrOrig->x.size;
      gridptrDup->x.vals = (double *) Malloc(size * sizeof(double));
      memcpy(gridptrDup->x.vals, gridptrOrig->x.vals, size * sizeof(double));
    }

  if (gridptrOrig->y.vals != NULL)
    {
      size_t size = isIrregular ? gridsize : gridptrOrig->y.size;
      gridptrDup->y.vals = (double *) Malloc(size * sizeof(double));
      memcpy(gridptrDup->y.vals, gridptrOrig->y.vals, size * sizeof(double));
    }

  if (gridptrOrig->x.bounds != NULL)
    {
      size_t size = (isIrregular ? gridsize : gridptrOrig->x.size) * gridptrOrig->nvertex;
      gridptrDup->x.bounds = (double *) Malloc(size * sizeof(double));
      memcpy(gridptrDup->x.bounds, gridptrOrig->x.bounds, size * sizeof(double));
    }

  if (gridptrOrig->y.bounds != NULL)
    {
      size_t size = (isIrregular ? gridsize : gridptrOrig->y.size) * gridptrOrig->nvertex;
      gridptrDup->y.bounds = (double *) Malloc(size * sizeof(double));
      memcpy(gridptrDup->y.bounds, gridptrOrig->y.bounds, size * sizeof(double));
    }

  {
    const double *gridptrOrig_area = gridptrOrig->vtable->inqAreaPtr(gridptrOrig);
    if (gridptrOrig_area != NULL)
      {
        size_t size = gridsize;
        gridptrDup->area = (double *) Malloc(size * sizeof(double));
        memcpy(gridptrDup->area, gridptrOrig_area, size * sizeof(double));
      }
  }

  if (gridptrOrig->mask != NULL)
    {
      size_t size = gridsize;
      gridptrDup->mask = (mask_t *) Malloc(size * sizeof(mask_t));
      memcpy(gridptrDup->mask, gridptrOrig->mask, size * sizeof(mask_t));
    }

  if (gridptrOrig->mask_gme != NULL)
    {
      size_t size = gridsize;
      gridptrDup->mask_gme = (mask_t *) Malloc(size * sizeof(mask_t));
      memcpy(gridptrDup->mask_gme, gridptrOrig->mask_gme, size * sizeof(mask_t));
    }
}

/*
@Function  gridDuplicate
@Title     Duplicate a horizontal Grid

@Prototype int gridDuplicate(int gridID)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate} or @fref{vlistInqVarGrid}.

@Description
The function @func{gridDuplicate} duplicates a horizontal Grid.

@Result
@func{gridDuplicate} returns an identifier to the duplicated Grid.

@EndFunction
*/
int
gridDuplicate(int gridID)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  grid_t *gridptrnew = gridptr->vtable->copy(gridptr);
  int gridIDnew = reshPut(gridptrnew, &gridOps);
  gridptrnew->self = gridIDnew;
  return gridIDnew;
}

void
gridCompress(int gridID)
{
  grid_t *gridptr = grid_to_pointer(gridID);

  int gridtype = gridInqType(gridID);
  if (gridtype == GRID_UNSTRUCTURED)
    {
      if (gridptr->mask_gme != NULL)
        {
          size_t gridsize = gridInqSize(gridID);
          size_t nv = (size_t) gridptr->nvertex;
          double *restrict area = (double *) gridptr->vtable->inqAreaPtr(gridptr),
                           *restrict xvals = (double *) gridptr->vtable->inqXValsPtr(gridptr),
                           *restrict yvals = (double *) gridptr->vtable->inqYValsPtr(gridptr),
                           *restrict xbounds = (double *) gridptr->vtable->inqXBoundsPtr(gridptr),
                           *restrict ybounds = (double *) gridptr->vtable->inqYBoundsPtr(gridptr);
          mask_t *restrict mask_gme = gridptr->mask_gme;
          size_t *restrict selection = (size_t *) Malloc(gridsize * sizeof(selection[0]));
          size_t nselect;
          {
            size_t j = 0;
            for (size_t i = 0; i < gridsize; i++) selection[j] = i, j += (mask_gme[i] != 0);
            nselect = j;
          }
          selection = (size_t *) Realloc(selection, nselect * sizeof(selection[0]));
          if (xvals)
            for (size_t i = 0; i < nselect; i++) xvals[i] = xvals[selection[i]];
          if (yvals)
            for (size_t i = 0; i < nselect; i++) yvals[i] = yvals[selection[i]];
          if (area)
            for (size_t i = 0; i < nselect; i++) area[i] = area[selection[i]];
          if (xbounds)
            for (size_t i = 0; i < nselect; i++)
              for (size_t iv = 0; iv < nv; iv++) xbounds[i * nv + iv] = xbounds[selection[i] * nv + iv];
          if (ybounds)
            for (size_t i = 0; i < nselect; i++)
              for (size_t iv = 0; iv < nv; iv++) ybounds[i * nv + iv] = ybounds[selection[i] * nv + iv];
          Free(selection);

          /* fprintf(stderr, "grid compress %d %d %d\n", i, j, gridsize); */
          gridsize = nselect;
          gridptr->size = (int) gridsize;
          gridptr->x.size = (int) gridsize;
          gridptr->y.size = (int) gridsize;

          double **resizeP[] = { &gridptr->x.vals, &gridptr->y.vals, &gridptr->area, &gridptr->x.bounds, &gridptr->y.bounds };
          size_t newSize[] = { gridsize, gridsize, gridsize, nv * gridsize, nv * gridsize };
          for (size_t i = 0; i < sizeof(resizeP) / sizeof(resizeP[0]); ++i)
            if (*(resizeP[i])) *(resizeP[i]) = (double *) Realloc(*(resizeP[i]), newSize[i] * sizeof(double));

          Free(gridptr->mask_gme);
          gridptr->mask_gme = NULL;
          gridMark4Update(gridID);
        }
    }
  else
    Warning("Unsupported grid type: %s", gridNamePtr(gridtype));
}

static void
gridDefAreaSerial(grid_t *gridptr, const double *area)
{
  size_t size = gridptr->size;

  if (size == 0) Error("size undefined for gridID = %d", gridptr->self);

  if (gridptr->area == NULL)
    gridptr->area = (double *) Malloc(size * sizeof(double));
  else if (CDI_Debug)
    Warning("values already defined!");

  memcpy(gridptr->area, area, size * sizeof(double));
}

void
gridDefArea(int gridID, const double *area)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  gridptr->vtable->defArea(gridptr, area);
  gridMark4Update(gridID);
}

static void
gridInqAreaSerial(grid_t *gridptr, double *area)
{
  if (gridptr->area) memcpy(area, gridptr->area, gridptr->size * sizeof(double));
}

void
gridInqArea(int gridID, double *area)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  gridptr->vtable->inqArea(gridptr, area);
}

static int
gridInqPropPresenceBase(grid_t *gridptr, enum gridPropInq inq)
{
  bool present = false;
  switch (inq)
    {
    case GRID_PROP_MASK: present = gridptr->mask != NULL; break;
    case GRID_PROP_MASK_GME: present = gridptr->mask != NULL; break;
    case GRID_PROP_AREA: present = gridptr->area != NULL; break;
    case GRID_PROP_XVALS: present = gridptr->x.vals != NULL; break;
    case GRID_PROP_YVALS: present = gridptr->y.vals != NULL; break;
    case GRID_PROP_XBOUNDS: present = gridptr->x.bounds != NULL; break;
    case GRID_PROP_YBOUNDS: present = gridptr->y.bounds != NULL; break;
    }
  return present;
}

int
gridInqPropPresence(int gridID, enum gridPropInq inq)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  return gridptr->vtable->inqPropPresence(gridptr, inq);
}

int
gridHasArea(int gridID)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  return gridptr->vtable->inqPropPresence(gridptr, GRID_PROP_AREA);
}

static const double *
gridInqAreaPtrBase(grid_t *gridptr)
{
  return gridptr->area;
}

const double *
gridInqAreaPtr(int gridID)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  return gridptr->vtable->inqAreaPtr(gridptr);
}

void
gridDefNvertex(int gridID, int nvertex)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  if (gridptr->nvertex != nvertex)
    {
      gridptr->nvertex = nvertex;
      gridMark4Update(gridID);
    }
}

int
gridInqNvertex(int gridID)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  return gridptr->nvertex;
}

static void
gridDefBoundsGeneric(grid_t *gridptr, const double *bounds, size_t regularSize, double **field)
{
  bool isIrregular = grid_is_irregular(gridptr->type);
  size_t nvertex = (size_t) gridptr->nvertex;
  if (nvertex == 0)
    {
      Warning("nvertex undefined for gridID = %d. Cannot define bounds!", gridptr->self);
      return;
    }

  size_t size = nvertex * (isIrregular ? gridptr->size : regularSize);
  if (size == 0) Error("size undefined for gridID = %d", gridptr->self);

  if (*field == NULL && size)
    *field = (double *) Malloc(size * sizeof(double));
  else if (CDI_Debug)
    Warning("values already defined!");

  copy_darray(size, bounds, *field);
}

static void
gridDefXBoundsSerial(grid_t *gridptr, const double *xbounds)
{
  gridDefBoundsGeneric(gridptr, xbounds, gridptr->x.size, &gridptr->x.bounds);
}

/*
@Function  gridDefXbounds
@Title     Define the bounds of a X-axis

@Prototype void gridDefXbounds(int gridID, const double *xbounds)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate}.
    @Item  xbounds  X-bounds of the grid.

@Description
The function @func{gridDefXbounds} defines all bounds of the X-axis.

@EndFunction
*/
void
gridDefXbounds(int gridID, const double *xbounds)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  gridptr->vtable->defXBounds(gridptr, xbounds);
  gridMark4Update(gridID);
}

static SizeType
gridInqXBoundsSerial(grid_t *gridptr, double *xbounds)
{
  size_t nvertex = (size_t) gridptr->nvertex;

  bool isIrregular = grid_is_irregular(gridptr->type);
  size_t size = nvertex * (isIrregular ? gridptr->size : gridptr->x.size);

  if (gridptr->x.bounds)
    {
      if (size && xbounds)
        {
          const double *gridptr_xbounds = gridptr->vtable->inqXBoundsPtr(gridptr);
          copy_darray(size, gridptr_xbounds, xbounds);
        }
    }
  else
    size = 0;

  return (SizeType) size;
}

/*
@Function  gridInqXbounds
@Title     Get the bounds of a X-axis

@Prototype SizeType gridInqXbounds(int gridID, double *xbounds)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate} or @fref{vlistInqVarGrid}.
    @Item  xbounds  Pointer to the location into which the X-bounds are read.
                    The caller must allocate space for the returned values.

@Description
The function @func{gridInqXbounds} returns the bounds of the X-axis.

@Result
Upon successful completion @func{gridInqXbounds} returns the number of bounds and
the bounds are stored in @func{xbounds}.
Otherwise, 0 is returned and @func{xbounds} is empty.

@EndFunction
*/
SizeType
gridInqXbounds(int gridID, double *xbounds)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  return gridptr->vtable->inqXBounds(gridptr, xbounds);
}

static const double *
gridInqXBoundsPtrSerial(grid_t *gridptr)
{
  return gridptr->x.bounds;
}

const double *
gridInqXboundsPtr(int gridID)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  return gridptr->vtable->inqXBoundsPtr(gridptr);
}

static void
gridDefYBoundsSerial(grid_t *gridptr, const double *ybounds)
{
  gridDefBoundsGeneric(gridptr, ybounds, gridptr->y.size, &gridptr->y.bounds);
}

//----------------------------------------------------------------------------
// Parallel Version
//----------------------------------------------------------------------------

SizeType
gridInqXboundsPart(int gridID, int start, SizeType size, double *xbounds)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  const double *gridptr_xbounds = gridptr->vtable->inqXBoundsPtr(gridptr);
  if (gridptr_xbounds && size && xbounds) memcpy(xbounds, gridptr_xbounds + start, size * sizeof(double));

  return size;
}

SizeType
gridInqYboundsPart(int gridID, int start, SizeType size, double *ybounds)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  const double *gridptr_ybounds = gridptr->vtable->inqYBoundsPtr(gridptr);
  if (gridptr_ybounds && size && ybounds) memcpy(ybounds, gridptr_ybounds + start, size * sizeof(double));

  return size;
}

/*
@Function  gridDefYbounds
@Title     Define the bounds of a Y-axis

@Prototype void gridDefYbounds(int gridID, const double *ybounds)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate}.
    @Item  ybounds  Y-bounds of the grid.

@Description
The function @func{gridDefYbounds} defines all bounds of the Y-axis.

@EndFunction
*/
void
gridDefYbounds(int gridID, const double *ybounds)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  gridptr->vtable->defYBounds(gridptr, ybounds);
  gridMark4Update(gridID);
}

static SizeType
gridInqYBoundsSerial(grid_t *gridptr, double *ybounds)
{
  size_t nvertex = (size_t) gridptr->nvertex;

  bool isIrregular = grid_is_irregular(gridptr->type);
  size_t size = nvertex * (isIrregular ? gridptr->size : gridptr->y.size);

  if (gridptr->y.bounds)
    {
      if (size && ybounds)
        {
          const double *gridptr_ybounds = gridptr->vtable->inqYBoundsPtr(gridptr);
          copy_darray(size, gridptr_ybounds, ybounds);
        }
    }
  else
    size = 0;

  return (SizeType) size;
}

/*
@Function  gridInqYbounds
@Title     Get the bounds of a Y-axis

@Prototype SizeType gridInqYbounds(int gridID, double *ybounds)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate} or @fref{vlistInqVarGrid}.
    @Item  ybounds  Pointer to the location into which the Y-bounds are read.
                    The caller must allocate space for the returned values.

@Description
The function @func{gridInqYbounds} returns the bounds of the Y-axis.

@Result
Upon successful completion @func{gridInqYbounds} returns the number of bounds and
the bounds are stored in @func{ybounds}.
Otherwise, 0 is returned and @func{ybounds} is empty.

@EndFunction
*/
SizeType
gridInqYbounds(int gridID, double *ybounds)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  return gridptr->vtable->inqYBounds(gridptr, ybounds);
}

static const double *
gridInqYBoundsPtrSerial(grid_t *gridptr)
{
  return gridptr->y.bounds;
}

const double *
gridInqYboundsPtr(int gridID)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  return gridptr->vtable->inqYBoundsPtr(gridptr);
}

static void
printDblsPrefixAutoBrk(FILE *fp, int dig, const char prefix[], size_t nbyte0, size_t n, const double vals[])
{
  fputs(prefix, fp);
  size_t nbyte = nbyte0;
  for (size_t i = 0; i < n; i++)
    {
      if (nbyte > 80)
        {
          fprintf(fp, "\n%*s", (int) nbyte0, "");
          nbyte = nbyte0;
        }
      nbyte += (size_t) fprintf(fp, "%.*g ", dig, vals[i]);
    }
  fputs("\n", fp);
}

static inline void *
resizeBuffer(void **buf, size_t *bufSize, size_t reqSize)
{
  if (reqSize > *bufSize)
    {
      *buf = Realloc(*buf, reqSize);
      *bufSize = reqSize;
    }
  return *buf;
}

static void
gridPrintAttributes(FILE *fp, int gridID)
{
  int cdiID = gridID;
  int varID = CDI_GLOBAL;
  int atttype, attlen;
  char attname[CDI_MAX_NAME + 1];
  void *attBuf = NULL;
  size_t attBufSize = 0;

  int natts;
  cdiInqNatts(cdiID, varID, &natts);

  for (int iatt = 0; iatt < natts; ++iatt)
    {
      cdiInqAtt(cdiID, varID, iatt, attname, &atttype, &attlen);

      if (attlen == 0) continue;

      if (atttype == CDI_DATATYPE_TXT)
        {
          size_t attSize = (size_t) (attlen + 1) * sizeof(char);
          char *atttxt = (char *) resizeBuffer(&attBuf, &attBufSize, attSize);
          cdiInqAttTxt(cdiID, varID, attname, attlen, atttxt);
          atttxt[attlen] = 0;
          fprintf(fp, "ATTR_TXT: %s = \"%s\"\n", attname, atttxt);
        }
      else if (atttype == CDI_DATATYPE_INT8 || atttype == CDI_DATATYPE_UINT8 || atttype == CDI_DATATYPE_INT16
               || atttype == CDI_DATATYPE_UINT16 || atttype == CDI_DATATYPE_INT32 || atttype == CDI_DATATYPE_UINT32)
        {
          size_t attSize = (size_t) attlen * sizeof(int);
          int *attint = (int *) resizeBuffer(&attBuf, &attBufSize, attSize);
          cdiInqAttInt(cdiID, varID, attname, attlen, &attint[0]);
          if (attlen == 1)
            fprintf(fp, "ATTR_INT: %s =", attname);
          else
            fprintf(fp, "ATTR_INT_%d: %s =", attlen, attname);
          for (int i = 0; i < attlen; ++i) fprintf(fp, " %d", attint[i]);
          fprintf(fp, "\n");
        }
      else if (atttype == CDI_DATATYPE_FLT32 || atttype == CDI_DATATYPE_FLT64)
        {
          size_t attSize = (size_t) attlen * sizeof(double);
          double *attflt = (double *) resizeBuffer(&attBuf, &attBufSize, attSize);
          int dig = (atttype == CDI_DATATYPE_FLT64) ? 15 : 7;
          cdiInqAttFlt(cdiID, varID, attname, attlen, attflt);
          if (attlen == 1)
            fprintf(fp, "ATTR_FLT: %s =", attname);
          else
            fprintf(fp, "ATTR_FLT_%d: %s =", attlen, attname);
          for (int i = 0; i < attlen; ++i) fprintf(fp, " %.*g", dig, attflt[i]);
          fprintf(fp, "\n");
        }
    }

  Free(attBuf);
}

static void
gridPrintKernel(int gridID, int opt, FILE *fp)
{
  char attstr[CDI_MAX_NAME];
  char attstr2[CDI_MAX_NAME];
  size_t nxvals = gridInqXvals(gridID, NULL);
  size_t nyvals = gridInqYvals(gridID, NULL);

  int type = gridInqType(gridID);
  size_t gridsize = gridInqSize(gridID);
  size_t xsize = gridInqXsize(gridID);
  size_t ysize = gridInqYsize(gridID);
  int nvertex = gridInqNvertex(gridID);
  int datatype;
  cdiInqKeyInt(gridID, CDI_GLOBAL, CDI_KEY_DATATYPE, &datatype);

  int dig = (datatype == CDI_DATATYPE_FLT64) ? 15 : 7;

  fprintf(fp,
          "gridtype  = %s\n"
          "gridsize  = %zu\n",
          gridNamePtr(type), gridsize);

  if (type != GRID_GME)
    {
      if (type != GRID_UNSTRUCTURED && type != GRID_SPECTRAL && type != GRID_FOURIER)
        {
          if (xsize > 0) fprintf(fp, "xsize     = %zu\n", xsize);
          if (ysize > 0) fprintf(fp, "ysize     = %zu\n", ysize);
        }

      if (nxvals > 0)
        {
          int length = CDI_MAX_NAME;
          cdiInqKeyString(gridID, CDI_XAXIS, CDI_KEY_NAME, attstr, &length);
          if (attstr[0]) fprintf(fp, "xname     = %s\n", attstr);
          length = CDI_MAX_NAME;
          cdiInqKeyString(gridID, CDI_XAXIS, CDI_KEY_DIMNAME, attstr2, &length);
          if (attstr2[0] && !str_is_equal(attstr, attstr2)) fprintf(fp, "xdimname  = %s\n", attstr2);
          length = CDI_MAX_NAME;
          cdiInqKeyString(gridID, CDI_XAXIS, CDI_KEY_LONGNAME, attstr, &length);
          if (attstr[0]) fprintf(fp, "xlongname = %s\n", attstr);
          length = CDI_MAX_NAME;
          cdiInqKeyString(gridID, CDI_XAXIS, CDI_KEY_UNITS, attstr, &length);
          if (attstr[0]) fprintf(fp, "xunits    = %s\n", attstr);
        }

      if (nyvals > 0)
        {
          int length = CDI_MAX_NAME;
          cdiInqKeyString(gridID, CDI_YAXIS, CDI_KEY_NAME, attstr, &length);
          if (attstr[0]) fprintf(fp, "yname     = %s\n", attstr);
          length = CDI_MAX_NAME;
          cdiInqKeyString(gridID, CDI_YAXIS, CDI_KEY_DIMNAME, attstr2, &length);
          if (attstr2[0] && !str_is_equal(attstr, attstr2)) fprintf(fp, "ydimname  = %s\n", attstr2);
          length = CDI_MAX_NAME;
          cdiInqKeyString(gridID, CDI_YAXIS, CDI_KEY_LONGNAME, attstr, &length);
          if (attstr[0]) fprintf(fp, "ylongname = %s\n", attstr);
          length = CDI_MAX_NAME;
          cdiInqKeyString(gridID, CDI_YAXIS, CDI_KEY_UNITS, attstr, &length);
          if (attstr[0]) fprintf(fp, "yunits    = %s\n", attstr);
        }

      if (type == GRID_UNSTRUCTURED && nvertex > 0) fprintf(fp, "nvertex   = %d\n", nvertex);
    }

  switch (type)
    {
    case GRID_LONLAT:
    case GRID_GAUSSIAN:
    case GRID_GAUSSIAN_REDUCED:
    case GRID_GENERIC:
    case GRID_PROJECTION:
    case GRID_CURVILINEAR:
    case GRID_UNSTRUCTURED:
    case GRID_CHARXY:
      {
        if (type == GRID_GAUSSIAN || type == GRID_GAUSSIAN_REDUCED) fprintf(fp, "np        = %d\n", gridInqNP(gridID));

        if (type == GRID_UNSTRUCTURED)
          {
            int number = 0;
            cdiInqKeyInt(gridID, CDI_GLOBAL, CDI_KEY_NUMBEROFGRIDUSED, &number);
            if (number > 0)
              {
                fprintf(fp, "number    = %d\n", number);
                int position = 0;
                cdiInqKeyInt(gridID, CDI_GLOBAL, CDI_KEY_NUMBEROFGRIDINREFERENCE, &position);
                if (position >= 0) fprintf(fp, "position  = %d\n", position);
              }

            int length;
            if (CDI_NOERR == cdiInqKeyLen(gridID, CDI_GLOBAL, CDI_KEY_REFERENCEURI, &length))
              {
                char reference_link[8192];
                length = sizeof(reference_link);
                cdiInqKeyString(gridID, CDI_GLOBAL, CDI_KEY_REFERENCEURI, reference_link, &length);
                fprintf(fp, "uri       = %s\n", reference_link);
              }
          }

        if (nxvals > 0)
          {
            double xfirst = 0.0, xinc = 0.0;

            if (type == GRID_LONLAT || type == GRID_GAUSSIAN || type == GRID_PROJECTION || type == GRID_GENERIC)
              {
                xfirst = gridInqXval(gridID, 0);
                xinc = gridInqXinc(gridID);
              }

            if (IS_NOT_EQUAL(xinc, 0) && opt)
              {
                fprintf(fp,
                        "xfirst    = %.*g\n"
                        "xinc      = %.*g\n",
                        dig, xfirst, dig, xinc);
              }
            else
              {
                double *xvals = (double *) Malloc(nxvals * sizeof(double));
                gridInqXvals(gridID, xvals);
                static const char prefix[] = "xvals     = ";
                printDblsPrefixAutoBrk(fp, dig, prefix, sizeof(prefix) - 1, nxvals, xvals);
                Free(xvals);
              }
          }

        if (nyvals > 0)
          {
            double yfirst = 0.0, yinc = 0.0;

            if (type == GRID_LONLAT || type == GRID_GENERIC || type == GRID_PROJECTION || type == GRID_GENERIC)
              {
                yfirst = gridInqYval(gridID, 0);
                yinc = gridInqYinc(gridID);
              }

            if (IS_NOT_EQUAL(yinc, 0) && opt)
              {
                fprintf(fp,
                        "yfirst    = %.*g\n"
                        "yinc      = %.*g\n",
                        dig, yfirst, dig, yinc);
              }
            else
              {
                double *yvals = (double *) Malloc(nyvals * sizeof(double));
                gridInqYvals(gridID, yvals);
                static const char prefix[] = "yvals     = ";
                printDblsPrefixAutoBrk(fp, dig, prefix, sizeof(prefix) - 1, nyvals, yvals);
                Free(yvals);
              }
          }

        if (type == GRID_PROJECTION) gridPrintAttributes(fp, gridID);

        break;
      }
    case GRID_SPECTRAL:
      {
        fprintf(fp,
                "truncation = %d\n"
                "complexpacking = %d\n",
                gridInqTrunc(gridID), gridInqComplexPacking(gridID));
        break;
      }
    case GRID_FOURIER:
      {
        fprintf(fp, "truncation = %d\n", gridInqTrunc(gridID));
        break;
      }
    case GRID_GME:
      {
        int nd, ni, ni2, ni3;
        gridInqParamGME(gridID, &nd, &ni, &ni2, &ni3);
        fprintf(fp, "ni        = %d\n", ni);
        break;
      }
    default:
      {
        fprintf(stderr, "Unsupported grid type: %s\n", gridNamePtr(type));
        break;
      }
    }
}

void
gridPrintP(void *voidptr, FILE *fp)
{
  grid_t *gridptr = (grid_t *) voidptr;
  int gridID = gridptr->self;

  xassert(gridptr);

  gridPrintKernel(gridID, 0, fp);

  int datatype = CDI_UNDEFID;
  cdiInqKeyInt(gridID, CDI_GLOBAL, CDI_KEY_DATATYPE, &datatype);

  fprintf(fp,
          "datatype  = %d\n"
          "nd        = %d\n"
          "ni        = %d\n"
          "ni2       = %d\n"
          "ni3       = %d\n"
          "trunc     = %d\n"
          "lcomplex  = %d\n"
          "reducedPointsSize   = %d\n",
          datatype, gridptr->gme.nd, gridptr->gme.ni, gridptr->gme.ni2, gridptr->gme.ni3, gridptr->trunc, gridptr->lcomplex,
          gridptr->reducedPointsSize);
}

static const double *
gridInqXValsPtrSerial(grid_t *gridptr)
{
  return gridptr->x.vals;
}

#ifndef USE_MPI
static const char **
gridInqXCvalsPtrSerial(grid_t *gridptr)
{
  return (const char **) gridptr->x.cvals;
}
#endif

const double *
gridInqXvalsPtr(int gridID)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  return gridptr->vtable->inqXValsPtr(gridptr);
}

#ifndef USE_MPI
const char **
gridInqXCvalsPtr(int gridID)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  return gridptr->vtable->inqXCvalsPtr(gridptr);
}
#endif

static const double *
gridInqYValsPtrSerial(grid_t *gridptr)
{
  return gridptr->y.vals;
}

#ifndef USE_MPI
static const char **
gridInqYCvalsPtrSerial(grid_t *gridptr)
{
  return (const char **) gridptr->y.cvals;
}
#endif

const double *
gridInqYvalsPtr(int gridID)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  return gridptr->vtable->inqYValsPtr(gridptr);
}

#ifndef USE_MPI
const char **
gridInqYCvalsPtr(int gridID)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  return gridptr->vtable->inqYCvalsPtr(gridptr);
}
#endif

void
gridProjParamsInit(struct CDI_GridProjParams *gpp)
{
  // clang-format off
  gpp->mv      = CDI_Grid_Missval;   // Missing value
  gpp->lon_0   = CDI_Grid_Missval;   // The East longitude of the meridian which is parallel to the Y-axis
  gpp->lat_0   = CDI_Grid_Missval;   // Latitude of the projection origin
  gpp->lat_1   = CDI_Grid_Missval;   // First latitude from the pole at which the secant cone cuts the sphere
  gpp->lat_2   = CDI_Grid_Missval;   // Second latitude at which the secant cone cuts the sphere
  gpp->a       = CDI_Grid_Missval;   // Semi-major axis or earth radius in metres (optional)
  gpp->b       = CDI_Grid_Missval;   // Semi-minor axis in metres (optional)
  gpp->rf      = CDI_Grid_Missval;   // Inverse flattening (1/f) (optional)
  gpp->xval_0  = CDI_Grid_Missval;   // Longitude of the first grid point in degree (optional)
  gpp->yval_0  = CDI_Grid_Missval;   // Latitude of the first grid point in degree (optional)
  gpp->x_0     = CDI_Grid_Missval;   // False easting (optional)
  gpp->y_0     = CDI_Grid_Missval;   // False northing (optional)
  gpp->x_SP    = CDI_Grid_Missval;   // Longitude of southern pole
  gpp->y_SP    = CDI_Grid_Missval;   // Latitude of southern pole
  // clang-format on
}

static void
gridDefParamsCommon(int gridID, struct CDI_GridProjParams gpp)
{
  if (IS_NOT_EQUAL(gpp.a, gpp.mv))
    {
      if (IS_NOT_EQUAL(gpp.b, gpp.mv))
        {
          cdiDefAttFlt(gridID, CDI_GLOBAL, "semi_major_axis", CDI_DATATYPE_FLT64, 1, &gpp.a);
          cdiDefAttFlt(gridID, CDI_GLOBAL, "semi_minor_axis", CDI_DATATYPE_FLT64, 1, &gpp.b);
        }
      else
        {
          cdiDefAttFlt(gridID, CDI_GLOBAL, "earth_radius", CDI_DATATYPE_FLT64, 1, &gpp.a);
        }
    }
  if (IS_NOT_EQUAL(gpp.rf, gpp.mv)) cdiDefAttFlt(gridID, CDI_GLOBAL, "inverse_flattening", CDI_DATATYPE_FLT64, 1, &gpp.rf);
  if (IS_NOT_EQUAL(gpp.x_0, gpp.mv)) cdiDefAttFlt(gridID, CDI_GLOBAL, "false_easting", CDI_DATATYPE_FLT64, 1, &gpp.x_0);
  if (IS_NOT_EQUAL(gpp.y_0, gpp.mv)) cdiDefAttFlt(gridID, CDI_GLOBAL, "false_northing", CDI_DATATYPE_FLT64, 1, &gpp.y_0);
  if (IS_NOT_EQUAL(gpp.xval_0, gpp.mv))
    cdiDefAttFlt(gridID, CDI_GLOBAL, "longitudeOfFirstGridPointInDegrees", CDI_DATATYPE_FLT64, 1, &gpp.xval_0);
  if (IS_NOT_EQUAL(gpp.yval_0, gpp.mv))
    cdiDefAttFlt(gridID, CDI_GLOBAL, "latitudeOfFirstGridPointInDegrees", CDI_DATATYPE_FLT64, 1, &gpp.yval_0);
  if (IS_NOT_EQUAL(gpp.x_SP, gpp.mv))
    cdiDefAttFlt(gridID, CDI_GLOBAL, "longitudeOfSouthernPoleInDegrees", CDI_DATATYPE_FLT64, 1, &gpp.x_SP);
  if (IS_NOT_EQUAL(gpp.y_SP, gpp.mv))
    cdiDefAttFlt(gridID, CDI_GLOBAL, "latitudeOfSouthernPoleInDegrees", CDI_DATATYPE_FLT64, 1, &gpp.y_SP);
}

/*
@Function  gridDefParamsLCC
@Title     Define the parameters of a Lambert Conformal Conic grid

@Prototype void gridDefParamsLCC(int gridID, struct CDI_GridProjParams gridProjParams)
@Parameter
    @Item  gridID          Grid ID, from a previous call to @fref{gridCreate}.
    @Item  gridProjParams  Grid projection parameters.

@Description
The function @func{gridDefParamsLCC} defines the parameters of a Lambert Conformal Conic grid.

@EndFunction
*/
void
gridDefParamsLCC(int gridID, struct CDI_GridProjParams gpp)
{
  cdiDefKeyString(gridID, CDI_GLOBAL, CDI_KEY_GRIDMAP_VARNAME, "Lambert_Conformal");

  const char *gmapname = "lambert_conformal_conic";
  cdiDefKeyString(gridID, CDI_GLOBAL, CDI_KEY_GRIDMAP_NAME, gmapname);
  cdiDefAttTxt(gridID, CDI_GLOBAL, "grid_mapping_name", (int) (strlen(gmapname)), gmapname);
  int nlats = 0;
  double lats[2];
  lats[nlats++] = gpp.lat_1;
  if (IS_NOT_EQUAL(gpp.lat_1, gpp.lat_2)) lats[nlats++] = gpp.lat_2;
  cdiDefAttFlt(gridID, CDI_GLOBAL, "standard_parallel", CDI_DATATYPE_FLT64, nlats, lats);
  cdiDefAttFlt(gridID, CDI_GLOBAL, "longitude_of_central_meridian", CDI_DATATYPE_FLT64, 1, &gpp.lon_0);
  cdiDefAttFlt(gridID, CDI_GLOBAL, "latitude_of_projection_origin", CDI_DATATYPE_FLT64, 1, &gpp.lat_0);

  gridDefParamsCommon(gridID, gpp);

  grid_t *gridptr = grid_to_pointer(gridID);
  gridptr->projtype = CDI_PROJ_LCC;

  if (gridptr->type != GRID_PROJECTION) gridptr->type = GRID_PROJECTION;

  gridVerifyProj(gridID);
}

/*
@Function  gridInqParamsLCC
@Title     Get the parameter of a Lambert Conformal Conic grid

@Prototype void gridInqParamsLCC(int gridID, struct CDI_GridProjParams *gpp)
@Parameter
    @Item  gridID          Grid ID, from a previous call to @fref{gridCreate} or @fref{vlistInqVarGrid}.
    @Item  gridProjParams  Grid projection parameters.

@Description
The function @func{gridInqParamsLCC} returns the parameter of a Lambert Conformal Conic grid.

@EndFunction
*/
int
gridInqParamsLCC(int gridID, struct CDI_GridProjParams *gpp)
{
  int status = -1;
  if (gridInqType(gridID) != GRID_PROJECTION) return status;

  gridProjParamsInit(gpp);

  status = -2;
  const char *projection = "lambert_conformal_conic";
  char gmapname[CDI_MAX_NAME];
  int length = CDI_MAX_NAME;
  cdiInqKeyString(gridID, CDI_GLOBAL, CDI_KEY_GRIDMAP_NAME, gmapname, &length);
  if (gmapname[0] && str_is_equal(gmapname, projection))
    {
      char attname[CDI_MAX_NAME + 1];

      int natts;
      cdiInqNatts(gridID, CDI_GLOBAL, &natts);

      if (natts) status = 0;

      for (int iatt = 0; iatt < natts; ++iatt)
        {
          int atttype, attlen;
          cdiInqAtt(gridID, CDI_GLOBAL, iatt, attname, &atttype, &attlen);
          if (attlen > 2) continue;

          double attflt[2];
          if (cdiInqAttConvertedToFloat(gridID, atttype, attname, attlen, attflt))
            {
              // clang-format off
              if      (str_is_equal(attname, "earth_radius"))                       gpp->a      = attflt[0];
              else if (str_is_equal(attname, "semi_major_axis"))                    gpp->a      = attflt[0];
              else if (str_is_equal(attname, "semi_minor_axis"))                    gpp->b      = attflt[0];
              else if (str_is_equal(attname, "inverse_flattening"))                 gpp->rf     = attflt[0];
              else if (str_is_equal(attname, "longitude_of_central_meridian"))      gpp->lon_0  = attflt[0];
              else if (str_is_equal(attname, "latitude_of_projection_origin"))      gpp->lat_0  = attflt[0];
              else if (str_is_equal(attname, "false_easting"))                      gpp->x_0    = attflt[0];
              else if (str_is_equal(attname, "false_northing"))                     gpp->y_0    = attflt[0];
              else if (str_is_equal(attname, "longitudeOfFirstGridPointInDegrees")) gpp->xval_0 = attflt[0];
              else if (str_is_equal(attname, "latitudeOfFirstGridPointInDegrees"))  gpp->yval_0 = attflt[0];
              else if (str_is_equal(attname, "longitudeOfSouthernPoleInDegrees"))   gpp->x_SP   = attflt[0];
              else if (str_is_equal(attname, "latitudeOfSouthernPoleInDegrees"))    gpp->y_SP   = attflt[0];
              else if (str_is_equal(attname, "standard_parallel"))
                {
                  gpp->lat_1 = attflt[0];
                  gpp->lat_2 = (attlen == 2) ? attflt[1] : attflt[0];
                }
              // clang-format on
            }
        }
    }

  return status;
}

int
gridVerifyProjParamsLCC(struct CDI_GridProjParams *gpp)
{
  static bool lwarn = true;

  if (lwarn)
    {
      // lwarn = false;
      const char *projection = "lambert_conformal_conic";
      if (IS_EQUAL(gpp->lon_0, gpp->mv)) Warning("%s mapping parameter %s missing!", projection, "longitude_of_central_meridian");
      if (IS_EQUAL(gpp->lat_0, gpp->mv)) Warning("%s mapping parameter %s missing!", projection, "latitude_of_central_meridian");
      if (IS_EQUAL(gpp->lat_1, gpp->mv)) Warning("%s mapping parameter %s missing!", projection, "standard_parallel");
      if (IS_NOT_EQUAL(gpp->x_0, gpp->mv) && IS_NOT_EQUAL(gpp->y_0, gpp->mv)
          && (IS_EQUAL(gpp->xval_0, gpp->mv) || IS_EQUAL(gpp->yval_0, gpp->mv)))
        {
          if (proj_lcc_to_lonlat_func)
            {
              gpp->xval_0 = -gpp->x_0;
              gpp->yval_0 = -gpp->y_0;
              proj_lcc_to_lonlat_func(*gpp, 0.0, 0.0, (SizeType) 1, &gpp->xval_0, &gpp->yval_0);
            }
          if (IS_EQUAL(gpp->xval_0, gpp->mv) || IS_EQUAL(gpp->yval_0, gpp->mv))
            Warning("%s mapping parameter %s missing!", projection,
                    "longitudeOfFirstGridPointInDegrees and latitudeOfFirstGridPointInDegrees");
        }
    }

  return 0;
}

int
gridVerifyProjParamsSTERE(struct CDI_GridProjParams *gpp)
{
  static bool lwarn = true;

  if (lwarn)
    {
      // lwarn = false;
      const char *projection = "polar_stereographic";
      if (IS_EQUAL(gpp->lon_0, gpp->mv))
        Warning("%s mapping parameter %s missing!", projection, "straight_vertical_longitude_from_pole");
      if (IS_EQUAL(gpp->lat_0, gpp->mv)) Warning("%s mapping parameter %s missing!", projection, "latitude_of_projection_origin");
      if (IS_EQUAL(gpp->lat_1, gpp->mv)) Warning("%s mapping parameter %s missing!", projection, "standard_parallel");
      if (IS_NOT_EQUAL(gpp->x_0, gpp->mv) && IS_NOT_EQUAL(gpp->y_0, gpp->mv)
          && (IS_EQUAL(gpp->xval_0, gpp->mv) || IS_EQUAL(gpp->yval_0, gpp->mv)))
        {
          if (proj_stere_to_lonlat_func)
            {
              gpp->xval_0 = -gpp->x_0;
              gpp->xval_0 = -gpp->y_0;
              proj_stere_to_lonlat_func(*gpp, 0.0, 0.0, (SizeType) 1, &gpp->xval_0, &gpp->yval_0);
            }
          if (IS_EQUAL(gpp->xval_0, gpp->mv) || IS_EQUAL(gpp->yval_0, gpp->mv))
            Warning("%s mapping parameter %s missing!", projection,
                    "longitudeOfFirstGridPointInDegrees and latitudeOfFirstGridPointInDegrees");
        }
    }

  return 0;
}

/*
@Function  gridDefParamsSTERE
@Title     Define the parameter of a Polar stereographic grid

@Prototype void gridDefParamsSTERE(int gridID, struct CDI_GridProjParams gridProjParams)
@Parameter
    @Item  gridID          Grid ID, from a previous call to @fref{gridCreate}.
    @Item  gridProjParams  Grid projection parameters.

@Description
The function @func{gridDefParamsSTERE} defines the parameter of a Polar stereographic grid.

@EndFunction
*/
void
gridDefParamsSTERE(int gridID, struct CDI_GridProjParams gpp)
{
  cdiDefKeyString(gridID, CDI_GLOBAL, CDI_KEY_GRIDMAP_VARNAME, "Polar_Stereographic");

  const char *gmapname = "polar_stereographic";
  cdiDefKeyString(gridID, CDI_GLOBAL, CDI_KEY_GRIDMAP_NAME, gmapname);
  cdiDefAttTxt(gridID, CDI_GLOBAL, "grid_mapping_name", (int) (strlen(gmapname)), gmapname);
  cdiDefAttFlt(gridID, CDI_GLOBAL, "standard_parallel", CDI_DATATYPE_FLT64, 1, &gpp.lat_1);
  cdiDefAttFlt(gridID, CDI_GLOBAL, "straight_vertical_longitude_from_pole", CDI_DATATYPE_FLT64, 1, &gpp.lon_0);
  cdiDefAttFlt(gridID, CDI_GLOBAL, "latitude_of_projection_origin", CDI_DATATYPE_FLT64, 1, &gpp.lat_0);

  gridDefParamsCommon(gridID, gpp);

  grid_t *gridptr = grid_to_pointer(gridID);
  gridptr->projtype = CDI_PROJ_STERE;

  gridVerifyProj(gridID);
}

/*
@Function  gridInqParamsSTERE
@Title     Get the parameter of a Polar stereographic grid

@Prototype void gridInqParamsSTERE(int gridID, struct CDI_GridProjParams *gpp)
@Parameter
    @Item  gridID    Grid ID, from a previous call to @fref{gridCreate} or @fref{vlistInqVarGrid}.
    @Item  gridProjParams  Grid projection parameters.

@Description
The function @func{gridInqParamsSTERE} returns the parameter of a Polar stereographic grid.

@EndFunction
*/
int
gridInqParamsSTERE(int gridID, struct CDI_GridProjParams *gpp)
{
  int status = -1;
  if (gridInqType(gridID) != GRID_PROJECTION) return status;

  gridProjParamsInit(gpp);

  status = -2;
  const char *projection = "polar_stereographic";
  char gmapname[CDI_MAX_NAME];
  int length = CDI_MAX_NAME;
  cdiInqKeyString(gridID, CDI_GLOBAL, CDI_KEY_GRIDMAP_NAME, gmapname, &length);
  if (gmapname[0] && str_is_equal(gmapname, projection))
    {
      int atttype, attlen;
      char attname[CDI_MAX_NAME + 1];

      int natts;
      cdiInqNatts(gridID, CDI_GLOBAL, &natts);

      if (natts) status = 0;

      for (int iatt = 0; iatt < natts; ++iatt)
        {
          cdiInqAtt(gridID, CDI_GLOBAL, iatt, attname, &atttype, &attlen);
          if (attlen > 2) continue;

          double attflt[2];
          if (cdiInqAttConvertedToFloat(gridID, atttype, attname, attlen, attflt))
            {
              // clang-format off
              if      (str_is_equal(attname, "earth_radius"))                          gpp->a      = attflt[0];
              else if (str_is_equal(attname, "semi_major_axis"))                       gpp->a      = attflt[0];
              else if (str_is_equal(attname, "semi_minor_axis"))                       gpp->b      = attflt[0];
              else if (str_is_equal(attname, "inverse_flattening"))                    gpp->rf     = attflt[0];
              else if (str_is_equal(attname, "standard_parallel"))                     gpp->lat_1  = attflt[0];
              else if (str_is_equal(attname, "straight_vertical_longitude_from_pole")) gpp->lon_0  = attflt[0];
              else if (str_is_equal(attname, "latitude_of_projection_origin"))         gpp->lat_0  = attflt[0];
              else if (str_is_equal(attname, "false_easting"))                         gpp->x_0    = attflt[0];
              else if (str_is_equal(attname, "false_northing"))                        gpp->y_0    = attflt[0];
              else if (str_is_equal(attname, "longitudeOfFirstGridPointInDegrees"))    gpp->xval_0 = attflt[0];
              else if (str_is_equal(attname, "latitudeOfFirstGridPointInDegrees"))     gpp->yval_0 = attflt[0];
              // clang-format on
            }
        }
    }

  return status;
}

void
gridDefComplexPacking(int gridID, int lcomplex)
{
  grid_t *gridptr = grid_to_pointer(gridID);

  if (gridptr->lcomplex != lcomplex)
    {
      gridptr->lcomplex = lcomplex != 0;
      gridMark4Update(gridID);
    }
}

int
gridInqComplexPacking(int gridID)
{
  grid_t *gridptr = grid_to_pointer(gridID);

  return (int) gridptr->lcomplex;
}

void
gridDefHasDims(int gridID, int hasdims)
{
  grid_t *gridptr = grid_to_pointer(gridID);

  if (gridptr->hasdims != (hasdims != 0))
    {
      gridptr->hasdims = hasdims != 0;
      gridMark4Update(gridID);
    }
}

int
gridInqHasDims(int gridID)
{
  grid_t *gridptr = grid_to_pointer(gridID);

  return (int) gridptr->hasdims;
}

/*
@Function  gridDefNumber
@Title     Define the reference number for an unstructured grid

@Prototype void gridDefNumber(int gridID, int number)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate}.
    @Item  number   Reference number for an unstructured grid.

@Description
The function @func{gridDefNumber} defines the reference number for an unstructured grid.

@EndFunction
*/
void
gridDefNumber(int gridID, int number)
{
  cdiDefKeyInt(gridID, CDI_GLOBAL, CDI_KEY_NUMBEROFGRIDUSED, number);
}

/*
@Function  gridInqNumber
@Title     Get the reference number to an unstructured grid

@Prototype int gridInqNumber(int gridID)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate} or @fref{vlistInqVarGrid}.

@Description
The function @func{gridInqNumber} returns the reference number to an unstructured grid.

@Result
@func{gridInqNumber} returns the reference number to an unstructured grid.
@EndFunction
*/
int
gridInqNumber(int gridID)
{
  int number = 0;
  cdiInqKeyInt(gridID, CDI_GLOBAL, CDI_KEY_NUMBEROFGRIDUSED, &number);
  return number;
}

/*
@Function  gridDefPosition
@Title     Define the position of grid in the reference file

@Prototype void gridDefPosition(int gridID, int position)
@Parameter
    @Item  gridID     Grid ID, from a previous call to @fref{gridCreate}.
    @Item  position   Position of grid in the reference file.

@Description
The function @func{gridDefPosition} defines the position of grid in the reference file.

@EndFunction
*/
void
gridDefPosition(int gridID, int position)
{
  cdiDefKeyInt(gridID, CDI_GLOBAL, CDI_KEY_NUMBEROFGRIDINREFERENCE, position);
}

/*
@Function  gridInqPosition
@Title     Get the position of grid in the reference file

@Prototype int gridInqPosition(int gridID)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate} or @fref{vlistInqVarGrid}.

@Description
The function @func{gridInqPosition} returns the position of grid in the reference file.

@Result
@func{gridInqPosition} returns the position of grid in the reference file.
@EndFunction
*/
int
gridInqPosition(int gridID)
{
  int position = 0;
  cdiInqKeyInt(gridID, CDI_GLOBAL, CDI_KEY_NUMBEROFGRIDINREFERENCE, &position);
  return position;
}

/*
@Function  gridDefReference
@Title     Define the reference URI for an unstructured grid

@Prototype void gridDefReference(int gridID, const char *reference)
@Parameter
    @Item  gridID      Grid ID, from a previous call to @fref{gridCreate}.
    @Item  reference   Reference URI for an unstructured grid.

@Description
The function @func{gridDefReference} defines the reference URI for an unstructured grid.

@EndFunction
*/
void
gridDefReference(int gridID, const char *reference)
{
  if (reference)
    {
      cdiDefKeyString(gridID, CDI_GLOBAL, CDI_KEY_REFERENCEURI, reference);
      gridMark4Update(gridID);
    }
}

/*
@Function  gridInqReference
@Title     Get the reference URI to an unstructured grid

@Prototype char *gridInqReference(int gridID, char *reference)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate} or @fref{vlistInqVarGrid}.

@Description
The function @func{gridInqReference} returns the reference URI to an unstructured grid.

@Result
@func{gridInqReference} returns the reference URI to an unstructured grid.
@EndFunction
*/
int
gridInqReference(int gridID, char *reference)
{
  int length = 0;
  if (CDI_NOERR == cdiInqKeyLen(gridID, CDI_GLOBAL, CDI_KEY_REFERENCEURI, &length))
    {
      if (reference) cdiInqKeyString(gridID, CDI_GLOBAL, CDI_KEY_REFERENCEURI, reference, &length);
    }

  return length;
}

/*
@Function  gridDefUUID
@Title     Define the UUID for an unstructured grid

@Prototype void gridDefUUID(int gridID, const char *uuid)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate}.
    @Item  uuid     UUID for an unstructured grid.

@Description
The function @func{gridDefUUID} defines the UUID for an unstructured grid.

@EndFunction
*/
void
gridDefUUID(int gridID, const unsigned char uuid[CDI_UUID_SIZE])
{
  cdiDefKeyBytes(gridID, CDI_GLOBAL, CDI_KEY_UUID, uuid, CDI_UUID_SIZE);

  gridMark4Update(gridID);
}

/*
@Function  gridInqUUID
@Title     Get the UUID to an unstructured grid

@Prototype void gridInqUUID(int gridID, char *uuid)
@Parameter
    @Item  gridID   Grid ID, from a previous call to @fref{gridCreate} or @fref{vlistInqVarGrid}.

@Description
The function @func{gridInqUUID} returns the UUID to an unstructured grid.

@Result
@func{gridInqUUID} returns the UUID to an unstructured grid to the parameter uuid.
@EndFunction
*/
void
gridInqUUID(int gridID, unsigned char uuid[CDI_UUID_SIZE])
{
  memset(uuid, 0, CDI_UUID_SIZE);
  int length = CDI_UUID_SIZE;
  cdiInqKeyBytes(gridID, CDI_GLOBAL, CDI_KEY_UUID, uuid, &length);
}

void
cdiGridGetIndexList(unsigned ngrids, int *gridIndexList)
{
  reshGetResHListOfType(ngrids, gridIndexList, &gridOps);
}

static int
gridTxCode(void *voidP)
{
  grid_t *gridptr = (grid_t *) voidP;
  return gridptr->vtable->txCode;
}

enum
{
  GRID_PACK_INT_IDX_SELF,
  GRID_PACK_INT_IDX_TYPE,
  GRID_PACK_INT_IDX_IS_CYCLIC,
  GRID_PACK_INT_IDX_X_FLAG,
  GRID_PACK_INT_IDX_Y_FLAG,
  GRID_PACK_INT_IDX_GME_ND,
  GRID_PACK_INT_IDX_GME_NI,
  GRID_PACK_INT_IDX_GME_NI2,
  GRID_PACK_INT_IDX_GME_NI3,
  GRID_PACK_INT_IDX_TRUNC,
  GRID_PACK_INT_IDX_NVERTEX,
  GRID_PACK_INT_IDX_REDUCED_POINTS_SIZE,
  GRID_PACK_INT_IDX_SIZE,
  GRID_PACK_INT_IDX_X_SIZE,
  GRID_PACK_INT_IDX_Y_SIZE,
  GRID_PACK_INT_IDX_LCOMPLEX,
  GRID_PACK_INT_IDX_MEMBERMASK,
  /*
  GRID_PACK_INT_IDX_XTSTDNNAME,
  GRID_PACK_INT_IDX_YTSTDNNAME,
  GRID_PACK_INT_IDX_ISCANSNEGATIVELY,
  GRID_PACK_INT_IDX_JSCANSPOSITIVELY,
  GRID_PACK_INT_IDX_JPOINTSARECONSECUTIVE,
  */
  gridNint
};

enum
{
  GRID_PACK_DBL_IDX_X_FIRST,
  GRID_PACK_DBL_IDX_Y_FIRST,
  GRID_PACK_DBL_IDX_X_LAST,
  GRID_PACK_DBL_IDX_Y_LAST,
  GRID_PACK_DBL_IDX_X_INC,
  GRID_PACK_DBL_IDX_Y_INC,
  gridNdouble
};

enum
{
  gridHasMaskFlag = 1 << 0,
  gridHasGMEMaskFlag = 1 << 1,
  gridHasXValsFlag = 1 << 2,
  gridHasYValsFlag = 1 << 3,
  gridHasAreaFlag = 1 << 4,
  gridHasXBoundsFlag = 1 << 5,
  gridHasYBoundsFlag = 1 << 6,
  gridHasReducedPointsFlag = 1 << 7,
};

static int
gridGetComponentFlags(const grid_t *gridP)
{
  int flags = 0;
  for (int prop = 0; prop < GRID_PROP_YBOUNDS + 1; ++prop)
    flags |= (gridP->vtable->inqPropPresence((grid_t *) gridP, (enum gridPropInq) prop) << prop);
  flags |= (gridHasReducedPointsFlag & (int) ((unsigned) (gridP->reducedPoints == NULL) - 1U));
  return flags;
}

static int
gridGetPackSize(void *voidP, void *context)
{
  grid_t *gridP = (grid_t *) voidP;
  return gridP->vtable->getPackSize(gridP, context);
}

static int gridGetPackSizeScalars(grid_t *gridP, void *context);

static int gridGetPackSizeArrays(grid_t *gridP, void *context);

static int
gridGetPackSizeBase(grid_t *gridP, void *context)
{
  return gridP->vtable->getPackSizeScalars(gridP, context) + gridP->vtable->getPackSizeArrays(gridP, context);
}

static int
gridGetPackSizeScalars(grid_t *gridP, void *context)
{
  int packBuffSize = 0, ui32PackSize = serializeGetSize(1, CDI_DATATYPE_UINT32, context);

  packBuffSize += serializeGetSize(gridNint, CDI_DATATYPE_INT, context) + ui32PackSize;

  packBuffSize += serializeGetSize(gridNdouble, CDI_DATATYPE_FLT64, context) + ui32PackSize;

  packBuffSize += serializeKeysGetPackSize(&gridP->keys, context);
  packBuffSize += serializeKeysGetPackSize(&gridP->x.keys, context);
  packBuffSize += serializeKeysGetPackSize(&gridP->y.keys, context);

  return packBuffSize;
}

static int
gridGetPackSizeArrays(grid_t *gridP, void *context)
{
  int packBuffSize = 0, count, ui32PackSize = serializeGetSize(1, CDI_DATATYPE_UINT32, context);

  if (gridP->reducedPoints)
    {
      xassert(gridP->reducedPointsSize);
      packBuffSize += serializeGetSize(gridP->reducedPointsSize, CDI_DATATYPE_INT, context) + ui32PackSize;
    }

  if (gridP->vtable->inqXValsPtr(gridP))
    {
      if (gridP->type == GRID_UNSTRUCTURED || gridP->type == GRID_CURVILINEAR)
        count = gridP->size;
      else
        count = gridP->x.size;
      xassert(count);
      packBuffSize += serializeGetSize(count, CDI_DATATYPE_FLT64, context) + ui32PackSize;
    }

  if (gridP->vtable->inqYValsPtr(gridP))
    {
      if (gridP->type == GRID_UNSTRUCTURED || gridP->type == GRID_CURVILINEAR)
        count = gridP->size;
      else
        count = gridP->y.size;
      xassert(count);
      packBuffSize += serializeGetSize(count, CDI_DATATYPE_FLT64, context) + ui32PackSize;
    }

  if (gridP->vtable->inqAreaPtr(gridP))
    {
      xassert(gridP->size);
      packBuffSize += serializeGetSize(gridP->size, CDI_DATATYPE_FLT64, context) + ui32PackSize;
    }

  if (gridP->x.bounds)
    {
      xassert(gridP->nvertex);
      count = grid_is_irregular(gridP->type) ? gridP->size : gridP->x.size;
      xassert(count);
      packBuffSize += (serializeGetSize(gridP->nvertex * count, CDI_DATATYPE_FLT64, context) + ui32PackSize);
    }

  if (gridP->y.bounds)
    {
      xassert(gridP->nvertex);
      count = grid_is_irregular(gridP->type) ? gridP->size : gridP->y.size;
      xassert(count);
      packBuffSize += (serializeGetSize(gridP->nvertex * count, CDI_DATATYPE_FLT64, context) + ui32PackSize);
    }

  if (gridP->mask)
    {
      xassert(gridP->size);
      packBuffSize += serializeGetSize(gridP->size, CDI_DATATYPE_UCHAR, context) + ui32PackSize;
    }

  if (gridP->mask_gme)
    {
      xassert(gridP->size);
      packBuffSize += serializeGetSize(gridP->size, CDI_DATATYPE_UCHAR, context) + ui32PackSize;
    }

  return packBuffSize;
}

static grid_t *gridUnpackScalars(char *unpackBuffer, int unpackBufferSize, int *unpackBufferPos, int originNamespace, void *context,
                                 int force_id, int *memberMaskP);

static void gridUnpackArrays(grid_t *gridP, int memberMask, char *unpackBuffer, int unpackBufferSize, int *unpackBufferPos,
                             int originNamespace, void *context);

int
gridUnpack(char *unpackBuffer, int unpackBufferSize, int *unpackBufferPos, int originNamespace, void *context, int force_id)
{
  gridInit();
  int memberMask;
  grid_t *gridP
      = gridUnpackScalars(unpackBuffer, unpackBufferSize, unpackBufferPos, originNamespace, context, force_id, &memberMask);
  gridP->vtable->unpackArrays(gridP, memberMask, unpackBuffer, unpackBufferSize, unpackBufferPos, originNamespace, context);
  reshSetStatus(gridP->self, &gridOps, reshGetStatus(gridP->self, &gridOps) & ~RESH_SYNC_BIT);
  return gridP->self;
}

static grid_t *
gridUnpackScalars(char *unpackBuffer, int unpackBufferSize, int *unpackBufferPos, int originNamespace, void *context, int force_id,
                  int *memberMaskP)
{
  grid_t *gridP;
  uint32_t d;
  int memberMask;
  {
    int intBuffer[gridNint];
    serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos, intBuffer, gridNint, CDI_DATATYPE_INT, context);
    serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos, &d, 1, CDI_DATATYPE_UINT32, context);

    xassert(cdiCheckSum(CDI_DATATYPE_INT, gridNint, intBuffer) == d);
    int targetID = namespaceAdaptKey(intBuffer[0], originNamespace);
    gridP = gridNewEntry(force_id ? targetID : CDI_UNDEFID);

    xassert(!force_id || targetID == gridP->self);

    gridP->type = intBuffer[GRID_PACK_INT_IDX_TYPE];
    gridP->isCyclic = (signed char) intBuffer[GRID_PACK_INT_IDX_IS_CYCLIC];
    gridP->x.flag = (short) intBuffer[GRID_PACK_INT_IDX_X_FLAG];
    gridP->y.flag = (short) intBuffer[GRID_PACK_INT_IDX_Y_FLAG];
    gridP->gme.nd = intBuffer[GRID_PACK_INT_IDX_GME_ND];
    gridP->gme.ni = intBuffer[GRID_PACK_INT_IDX_GME_NI];
    gridP->gme.ni2 = intBuffer[GRID_PACK_INT_IDX_GME_NI2];
    gridP->gme.ni3 = intBuffer[GRID_PACK_INT_IDX_GME_NI3];
    gridP->trunc = intBuffer[GRID_PACK_INT_IDX_TRUNC];
    gridP->nvertex = intBuffer[GRID_PACK_INT_IDX_NVERTEX];
    gridP->reducedPointsSize = intBuffer[GRID_PACK_INT_IDX_REDUCED_POINTS_SIZE];
    gridP->size = intBuffer[GRID_PACK_INT_IDX_SIZE];
    gridP->x.size = intBuffer[GRID_PACK_INT_IDX_X_SIZE];
    gridP->y.size = intBuffer[GRID_PACK_INT_IDX_Y_SIZE];
    gridP->lcomplex = (bool) intBuffer[GRID_PACK_INT_IDX_LCOMPLEX];
    memberMask = intBuffer[GRID_PACK_INT_IDX_MEMBERMASK];
  }

  {
    double doubleBuffer[gridNdouble];
    serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos, doubleBuffer, gridNdouble, CDI_DATATYPE_FLT64, context);
    serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos, &d, 1, CDI_DATATYPE_UINT32, context);
    xassert(d == cdiCheckSum(CDI_DATATYPE_FLT, gridNdouble, doubleBuffer));

    gridP->x.first = doubleBuffer[GRID_PACK_DBL_IDX_X_FIRST];
    gridP->y.first = doubleBuffer[GRID_PACK_DBL_IDX_Y_FIRST];
    gridP->x.last = doubleBuffer[GRID_PACK_DBL_IDX_X_LAST];
    gridP->y.last = doubleBuffer[GRID_PACK_DBL_IDX_Y_LAST];
    gridP->x.inc = doubleBuffer[GRID_PACK_DBL_IDX_X_INC];
    gridP->y.inc = doubleBuffer[GRID_PACK_DBL_IDX_Y_INC];
  }

  serializeKeysUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos, &gridP->keys, context);
  serializeKeysUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos, &gridP->x.keys, context);
  serializeKeysUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos, &gridP->y.keys, context);

  *memberMaskP = memberMask;
  return gridP;
}

static void
gridUnpackArrays(grid_t *gridP, int memberMask, char *unpackBuffer, int unpackBufferSize, int *unpackBufferPos, int originNamespace,
                 void *context)
{
  UNUSED(originNamespace);
  uint32_t d;

  if (memberMask & gridHasReducedPointsFlag)
    {
      xassert(gridP->reducedPointsSize);
      gridP->reducedPoints = (int *) Malloc((size_t) gridP->reducedPointsSize * sizeof(int));
      serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos, gridP->reducedPoints, gridP->reducedPointsSize,
                      CDI_DATATYPE_INT, context);
      serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos, &d, 1, CDI_DATATYPE_UINT32, context);
      xassert(cdiCheckSum(CDI_DATATYPE_INT, gridP->reducedPointsSize, gridP->reducedPoints) == d);
    }

  bool isIrregular = grid_is_irregular(gridP->type);
  if (memberMask & gridHasXValsFlag)
    {
      int size = isIrregular ? gridP->size : gridP->x.size;

      gridP->x.vals = (double *) Malloc(size * sizeof(double));
      serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos, gridP->x.vals, size, CDI_DATATYPE_FLT64, context);
      serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos, &d, 1, CDI_DATATYPE_UINT32, context);
      xassert(cdiCheckSum(CDI_DATATYPE_FLT, size, gridP->x.vals) == d);
    }

  if (memberMask & gridHasYValsFlag)
    {
      int size = isIrregular ? gridP->size : gridP->y.size;

      gridP->y.vals = (double *) Malloc(size * sizeof(double));
      serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos, gridP->y.vals, size, CDI_DATATYPE_FLT64, context);
      serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos, &d, 1, CDI_DATATYPE_UINT32, context);
      xassert(cdiCheckSum(CDI_DATATYPE_FLT, size, gridP->y.vals) == d);
    }

  if (memberMask & gridHasAreaFlag)
    {
      int size = gridP->size;
      xassert(size);
      gridP->area = (double *) Malloc(size * sizeof(double));
      serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos, gridP->area, size, CDI_DATATYPE_FLT64, context);
      serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos, &d, 1, CDI_DATATYPE_UINT32, context);
      xassert(cdiCheckSum(CDI_DATATYPE_FLT, size, gridP->area) == d);
    }

  if (memberMask & gridHasXBoundsFlag)
    {
      int size = gridP->nvertex * (isIrregular ? gridP->size : gridP->x.size);
      xassert(size);

      gridP->x.bounds = (double *) Malloc(size * sizeof(double));
      serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos, gridP->x.bounds, size, CDI_DATATYPE_FLT64, context);
      serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos, &d, 1, CDI_DATATYPE_UINT32, context);
      xassert(cdiCheckSum(CDI_DATATYPE_FLT, size, gridP->x.bounds) == d);
    }

  if (memberMask & gridHasYBoundsFlag)
    {
      int size = gridP->nvertex * (isIrregular ? gridP->size : gridP->y.size);
      xassert(size);

      gridP->y.bounds = (double *) Malloc(size * sizeof(double));
      serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos, gridP->y.bounds, size, CDI_DATATYPE_FLT64, context);
      serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos, &d, 1, CDI_DATATYPE_UINT32, context);
      xassert(cdiCheckSum(CDI_DATATYPE_FLT, size, gridP->y.bounds) == d);
    }

  if (memberMask & gridHasMaskFlag)
    {
      int size = gridP->size;
      xassert(size);
      gridP->mask = (mask_t *) Malloc(size * sizeof(mask_t));
      serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos, gridP->mask, gridP->size, CDI_DATATYPE_UCHAR, context);
      serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos, &d, 1, CDI_DATATYPE_UINT32, context);
      xassert(cdiCheckSum(CDI_DATATYPE_UCHAR, gridP->size, gridP->mask) == d);
    }

  if (memberMask & gridHasGMEMaskFlag)
    {
      int size = gridP->size;
      xassert(size);
      gridP->mask_gme = (mask_t *) Malloc(size * sizeof(mask_t));
      serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos, gridP->mask_gme, gridP->size, CDI_DATATYPE_UCHAR, context);
      serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos, &d, 1, CDI_DATATYPE_UINT32, context);
      xassert(cdiCheckSum(CDI_DATATYPE_UCHAR, gridP->size, gridP->mask_gme) == d);
    }
}

void
gridPack(void *voidP, void *packBuffer, int packBufferSize, int *packBufferPos, void *context)
{
  grid_t *gridP = (grid_t *) voidP;
  gridP->vtable->pack(gridP, packBuffer, packBufferSize, packBufferPos, context);
}

static void
gridPackBase(grid_t *gridP, void *packBuffer, int packBufferSize, int *packBufferPos, void *context)
{
  int memberMask = gridP->vtable->packScalars(gridP, packBuffer, packBufferSize, packBufferPos, context);
  gridP->vtable->packArrays(gridP, memberMask, packBuffer, packBufferSize, packBufferPos, context);
}

static int
gridPackScalars(grid_t *gridP, void *packBuffer, int packBufferSize, int *packBufferPos, void *context)
{
  uint32_t d;
  int memberMask;

  {
    int intBuffer[gridNint];

    intBuffer[GRID_PACK_INT_IDX_SELF] = gridP->self;
    intBuffer[GRID_PACK_INT_IDX_TYPE] = gridP->type;
    intBuffer[GRID_PACK_INT_IDX_IS_CYCLIC] = gridP->isCyclic;
    intBuffer[GRID_PACK_INT_IDX_X_FLAG] = gridP->x.flag;
    intBuffer[GRID_PACK_INT_IDX_Y_FLAG] = gridP->y.flag;
    intBuffer[GRID_PACK_INT_IDX_GME_ND] = gridP->gme.nd;
    intBuffer[GRID_PACK_INT_IDX_GME_NI] = gridP->gme.ni;
    intBuffer[GRID_PACK_INT_IDX_GME_NI2] = gridP->gme.ni2;
    intBuffer[GRID_PACK_INT_IDX_GME_NI3] = gridP->gme.ni3;
    intBuffer[GRID_PACK_INT_IDX_TRUNC] = gridP->trunc;
    intBuffer[GRID_PACK_INT_IDX_NVERTEX] = gridP->nvertex;
    intBuffer[GRID_PACK_INT_IDX_REDUCED_POINTS_SIZE] = gridP->reducedPointsSize;
    intBuffer[GRID_PACK_INT_IDX_SIZE] = gridP->size;
    intBuffer[GRID_PACK_INT_IDX_X_SIZE] = gridP->x.size;
    intBuffer[GRID_PACK_INT_IDX_Y_SIZE] = gridP->y.size;
    intBuffer[GRID_PACK_INT_IDX_LCOMPLEX] = gridP->lcomplex;
    intBuffer[GRID_PACK_INT_IDX_MEMBERMASK] = memberMask = gridGetComponentFlags(gridP);

    serializePack(intBuffer, gridNint, CDI_DATATYPE_INT, packBuffer, packBufferSize, packBufferPos, context);
    d = cdiCheckSum(CDI_DATATYPE_INT, gridNint, intBuffer);
    serializePack(&d, 1, CDI_DATATYPE_UINT32, packBuffer, packBufferSize, packBufferPos, context);
  }

  {
    double doubleBuffer[gridNdouble];

    doubleBuffer[GRID_PACK_DBL_IDX_X_FIRST] = gridP->x.first;
    doubleBuffer[GRID_PACK_DBL_IDX_Y_FIRST] = gridP->y.first;
    doubleBuffer[GRID_PACK_DBL_IDX_X_LAST] = gridP->x.last;
    doubleBuffer[GRID_PACK_DBL_IDX_Y_LAST] = gridP->y.last;
    doubleBuffer[GRID_PACK_DBL_IDX_X_INC] = gridP->x.inc;
    doubleBuffer[GRID_PACK_DBL_IDX_Y_INC] = gridP->y.inc;

    serializePack(doubleBuffer, gridNdouble, CDI_DATATYPE_FLT64, packBuffer, packBufferSize, packBufferPos, context);
    d = cdiCheckSum(CDI_DATATYPE_FLT, gridNdouble, doubleBuffer);
    serializePack(&d, 1, CDI_DATATYPE_UINT32, packBuffer, packBufferSize, packBufferPos, context);
  }

  serializeKeysPack(&gridP->keys, packBuffer, packBufferSize, packBufferPos, context);
  serializeKeysPack(&gridP->x.keys, packBuffer, packBufferSize, packBufferPos, context);
  serializeKeysPack(&gridP->y.keys, packBuffer, packBufferSize, packBufferPos, context);

  return memberMask;
}

static void
gridPackArrays(grid_t *gridP, int memberMask, void *packBuffer, int packBufferSize, int *packBufferPos, void *context)
{
  uint32_t d;
  bool isIrregular = grid_is_irregular(gridP->type);

  if (memberMask & gridHasReducedPointsFlag)
    {
      int size = gridP->reducedPointsSize;
      xassert(size > 0);
      serializePack(gridP->reducedPoints, size, CDI_DATATYPE_INT, packBuffer, packBufferSize, packBufferPos, context);
      d = cdiCheckSum(CDI_DATATYPE_INT, size, gridP->reducedPoints);
      serializePack(&d, 1, CDI_DATATYPE_UINT32, packBuffer, packBufferSize, packBufferPos, context);
    }

  if (memberMask & gridHasXValsFlag)
    {
      int size = isIrregular ? gridP->size : gridP->x.size;
      xassert(size);

      const double *gridP_xvals = gridP->vtable->inqXValsPtr(gridP);
      serializePack(gridP_xvals, size, CDI_DATATYPE_FLT64, packBuffer, packBufferSize, packBufferPos, context);
      d = cdiCheckSum(CDI_DATATYPE_FLT, size, gridP_xvals);
      serializePack(&d, 1, CDI_DATATYPE_UINT32, packBuffer, packBufferSize, packBufferPos, context);
    }

  if (memberMask & gridHasYValsFlag)
    {
      int size = isIrregular ? gridP->size : gridP->y.size;
      xassert(size);
      const double *gridP_yvals = gridP->vtable->inqYValsPtr(gridP);
      serializePack(gridP_yvals, size, CDI_DATATYPE_FLT64, packBuffer, packBufferSize, packBufferPos, context);
      d = cdiCheckSum(CDI_DATATYPE_FLT, size, gridP_yvals);
      serializePack(&d, 1, CDI_DATATYPE_UINT32, packBuffer, packBufferSize, packBufferPos, context);
    }

  if (memberMask & gridHasAreaFlag)
    {
      int size = gridP->size;
      xassert(size);

      serializePack(gridP->area, size, CDI_DATATYPE_FLT64, packBuffer, packBufferSize, packBufferPos, context);
      d = cdiCheckSum(CDI_DATATYPE_FLT, size, gridP->area);
      serializePack(&d, 1, CDI_DATATYPE_UINT32, packBuffer, packBufferSize, packBufferPos, context);
    }

  if (memberMask & gridHasXBoundsFlag)
    {
      xassert(gridP->nvertex);
      int size = isIrregular ? gridP->nvertex * gridP->size : gridP->nvertex * gridP->x.size;
      xassert(size);

      serializePack(gridP->x.bounds, size, CDI_DATATYPE_FLT64, packBuffer, packBufferSize, packBufferPos, context);
      d = cdiCheckSum(CDI_DATATYPE_FLT, size, gridP->x.bounds);
      serializePack(&d, 1, CDI_DATATYPE_UINT32, packBuffer, packBufferSize, packBufferPos, context);
    }

  if (memberMask & gridHasYBoundsFlag)
    {
      xassert(gridP->nvertex);
      int size = isIrregular ? gridP->nvertex * gridP->size : gridP->nvertex * gridP->y.size;
      xassert(size);

      serializePack(gridP->y.bounds, size, CDI_DATATYPE_FLT64, packBuffer, packBufferSize, packBufferPos, context);
      d = cdiCheckSum(CDI_DATATYPE_FLT, size, gridP->y.bounds);
      serializePack(&d, 1, CDI_DATATYPE_UINT32, packBuffer, packBufferSize, packBufferPos, context);
    }

  if (memberMask & gridHasMaskFlag)
    {
      int size = gridP->size;
      xassert(size);
      serializePack(gridP->mask, size, CDI_DATATYPE_UCHAR, packBuffer, packBufferSize, packBufferPos, context);
      d = cdiCheckSum(CDI_DATATYPE_UCHAR, size, gridP->mask);
      serializePack(&d, 1, CDI_DATATYPE_UINT32, packBuffer, packBufferSize, packBufferPos, context);
    }

  if (memberMask & gridHasGMEMaskFlag)
    {
      int size = gridP->size;
      xassert(size);

      serializePack(gridP->mask_gme, size, CDI_DATATYPE_UCHAR, packBuffer, packBufferSize, packBufferPos, context);
      d = cdiCheckSum(CDI_DATATYPE_UCHAR, size, gridP->mask_gme);
      serializePack(&d, 1, CDI_DATATYPE_UINT32, packBuffer, packBufferSize, packBufferPos, context);
    }
}

struct gridCompareSearchState
{
  int resIDValue;
  const grid_t *queryKey;
};

static enum cdiApplyRet
gridCompareSearch(int id, void *res, void *data)
{
  struct gridCompareSearchState *state = (struct gridCompareSearchState *) data;
  (void) res;
  if (gridCompare(id, state->queryKey, true) == false)
    {
      state->resIDValue = id;
      return CDI_APPLY_STOP;
    }
  else
    return CDI_APPLY_GO_ON;
}

// Add grid (which must be Malloc'ed to vlist if not already found)
struct addIfNewRes
cdiVlistAddGridIfNew(int vlistID, grid_t *grid, int mode)
{
  /*
    mode: 0 search in vlist and grid table
          1 search in grid table only
          2 search in grid table only and don't store the grid in vlist
   */
  bool gridIsDefinedGlobal = false;
  bool gridIsDefined = false;
  int gridID = CDI_UNDEFID;
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  int ngrids = vlistptr->ngrids;

  if (mode == 0)
    for (int index = 0; index < ngrids; index++)
      {
        if ((gridID = vlistptr->gridIDs[index]) != CDI_UNDEFID)
          {
            if (gridCompare(gridID, grid, false) == false)
              {
                gridIsDefined = true;
                break;
              }
          }
        else
          Error("Internal problem: undefined gridID in vlist %d, position %u!", vlistID, index);
      }

  if (!gridIsDefined)
    {
      struct gridCompareSearchState query;
      query.queryKey = grid;  // = { .queryKey = grid };
      if ((gridIsDefinedGlobal = (cdiGridApply(gridCompareSearch, &query) == CDI_APPLY_STOP))) gridID = query.resIDValue;

      if (mode == 1 && gridIsDefinedGlobal)
        for (int index = 0; index < ngrids; index++)
          if (vlistptr->gridIDs[index] == gridID)
            {
              gridIsDefinedGlobal = false;
              break;
            }
    }

  if (!gridIsDefined)
    {
      if (!gridIsDefinedGlobal)
        {
          grid->self = gridID = reshPut(grid, &gridOps);
          grid_complete(grid);
        }
      if (mode < 2)
        {
          if (ngrids >= MAX_GRIDS_PS) Error("Internal limit exceeded, MAX_GRIDS_PS=%d needs to be increased!", MAX_GRIDS_PS);
          vlistptr->gridIDs[ngrids] = gridID;
          vlistptr->ngrids++;
        }
    }

  return (struct addIfNewRes){ .Id = gridID, .isNew = (!gridIsDefined && !gridIsDefinedGlobal) };
}

const struct gridVirtTable cdiGridVtable = {
  .destroy = gridDestroyKernel,
  .copy = grid_copy_base,
  .copyScalarFields = grid_copy_base_scalar_fields,
  .copyArrayFields = grid_copy_base_array_fields,
  .defXVals = gridDefXValsSerial,
  .defYVals = gridDefYValsSerial,
  .defMask = gridDefMaskSerial,
  .defMaskGME = gridDefMaskGMESerial,
  .defXBounds = gridDefXBoundsSerial,
  .defYBounds = gridDefYBoundsSerial,
  .defArea = gridDefAreaSerial,
  .inqXVal = gridInqXValSerial,
  .inqYVal = gridInqYValSerial,
  .inqXVals = gridInqXValsSerial,
  .inqXValsPart = gridInqXValsPartSerial,
  .inqYVals = gridInqYValsSerial,
  .inqYValsPart = gridInqYValsPartSerial,
  .inqXValsPtr = gridInqXValsPtrSerial,
  .inqYValsPtr = gridInqYValsPtrSerial,
#ifndef USE_MPI
  .inqXIsc = gridInqXIscSerial,
  .inqYIsc = gridInqYIscSerial,
  .inqXCvals = gridInqXCvalsSerial,
  .inqYCvals = gridInqYCvalsSerial,
  .inqXCvalsPtr = gridInqXCvalsPtrSerial,
  .inqYCvalsPtr = gridInqYCvalsPtrSerial,
#endif
  .inqXInc = gridInqXIncBase,
  .inqYInc = gridInqYIncBase,
  .compareXYFull = compareXYvals,
  .compareXYAO = compareXYvals2,
  .inqArea = gridInqAreaSerial,
  .inqAreaPtr = gridInqAreaPtrBase,
  .inqPropPresence = gridInqPropPresenceBase,
  .inqMask = gridInqMaskSerial,
  .inqMaskGME = gridInqMaskGMESerial,
  .inqXBounds = gridInqXBoundsSerial,
  .inqYBounds = gridInqYBoundsSerial,
  .inqXBoundsPtr = gridInqXBoundsPtrSerial,
  .inqYBoundsPtr = gridInqYBoundsPtrSerial,
  .txCode = GRID,
  .getPackSize = gridGetPackSizeBase,
  .getPackSizeScalars = gridGetPackSizeScalars,
  .getPackSizeArrays = gridGetPackSizeArrays,
  .unpackScalars = gridUnpackScalars,
  .unpackArrays = gridUnpackArrays,
  .pack = gridPackBase,
  .packScalars = gridPackScalars,
  .packArrays = gridPackArrays,
};

/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
