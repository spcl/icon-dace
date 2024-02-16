#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_LIBNETCDF
#include "stream_cdf.h"
#include "cdf_lazy_grid.h"

static struct gridVirtTable cdfLazyGridVtable;
double *cdfPendingLoad;
#ifdef HAVE_LIBPTHREAD
static pthread_once_t cdfLazyInitialized = PTHREAD_ONCE_INIT;
#else
static bool cdfLazyInitialized;
#endif

#ifdef HAVE_LIBPTHREAD
#define lock_lazy_load(plGrid) pthread_mutex_lock(&((plGrid)->loadSerialize))
#define unlock_lazy_load(plGrid) pthread_mutex_unlock(&((plGrid)->loadSerialize))
#define destroy_lazy_load_lock(plGrid) pthread_mutex_destroy(&((plGrid)->loadSerialize))
#define init_lazy_load_lock(plGrid) pthread_mutex_init(&((plGrid)->loadSerialize), NULL)
#else
#define lock_lazy_load(plGrid)
#define unlock_lazy_load(plGrid)
#define destroy_lazy_load_lock(plGrid)
#define init_lazy_load_lock(plGrid)
#endif

void
cdfLazyGridDestroy(struct cdfLazyGrid *lazyGrid)
{
  lazyGrid->base.extraData = NULL;
  if (lazyGrid->base.area == cdfPendingLoad) lazyGrid->base.area = NULL;
  if (lazyGrid->base.x.vals == cdfPendingLoad) lazyGrid->base.x.vals = NULL;
  if (lazyGrid->base.y.vals == cdfPendingLoad) lazyGrid->base.y.vals = NULL;
  if (lazyGrid->base.x.bounds == cdfPendingLoad) lazyGrid->base.x.bounds = NULL;
  if (lazyGrid->base.y.bounds == cdfPendingLoad) lazyGrid->base.y.bounds = NULL;
  destroy_lazy_load_lock(lazyGrid);
}

static void
cdfLazyGridDelete(grid_t *grid)
{
  struct cdfLazyGrid *cdfGrid = (struct cdfLazyGrid *) grid;
  void (*baseDestroy)(grid_t * grid) = cdfGrid->baseVtable->destroy;
  cdfLazyGridDestroy(cdfGrid);
  baseDestroy(grid);
}

static void
cdfLazyGridDestroyOnce(void)
{
  /*
#ifdef HAVE_MMAP
  size_t pgSize = cdiGetPageSize(false);
  munmap(cdfPendingLoad, pgSize);
#endif
  */
}

static void
cdfLazyGridDefArea(grid_t *grid, const double *area)
{
  struct cdfLazyGrid *cdfGrid = (struct cdfLazyGrid *) grid;
  lock_lazy_load(cdfGrid);
  if (grid->area == cdfPendingLoad) grid->area = NULL;
  cdfGrid->cellAreaGet.datasetNCId = -1;
  cdfGrid->cellAreaGet.varNCId = -1;
  cdfGrid->baseVtable->defArea(grid, area);
  unlock_lazy_load(cdfGrid);
}

static const double *
cdfLazyGridInqAreaPtr(grid_t *grid)
{
  struct cdfLazyGrid *lazyGrid = (struct cdfLazyGrid *) grid;
  lock_lazy_load(lazyGrid);
  if (grid->area == cdfPendingLoad)
    {
      grid->area = (double *) Malloc(grid->size * sizeof(double));
      cdf_get_var_double(lazyGrid->cellAreaGet.datasetNCId, lazyGrid->cellAreaGet.varNCId, grid->area);
    }
  unlock_lazy_load(lazyGrid);
  return lazyGrid->baseVtable->inqAreaPtr(grid);
}

static void
cdfLazyGridInqArea(grid_t *grid, double *area)
{
  grid->vtable->inqAreaPtr(grid);
  struct cdfLazyGrid *lazyGrid = (struct cdfLazyGrid *) grid;
  lazyGrid->baseVtable->inqArea(grid, area);
}

static void
cdfLazyLoadXYVals(struct xyValGet *valsGet, double **valsp)
{
  double *grid_vals = (double *) Malloc(valsGet->size * sizeof(double));
  *valsp = grid_vals;
  if (valsGet->ndims == 3)
    cdf_get_vara_double(valsGet->datasetNCId, valsGet->varNCId, valsGet->start, valsGet->count, grid_vals);
  else
    cdf_get_var_double(valsGet->datasetNCId, valsGet->varNCId, grid_vals);
  cdf_scale_add(valsGet->size, grid_vals, valsGet->addoffset, valsGet->scalefactor);
}

static const double *
cdfLazyGridInqXValsPtr(grid_t *grid)
{
  struct cdfLazyGrid *lazyGrid = (struct cdfLazyGrid *) grid;
  lock_lazy_load(lazyGrid);
  if (grid->x.vals == cdfPendingLoad) cdfLazyLoadXYVals(&lazyGrid->xValsGet, &grid->x.vals);
  unlock_lazy_load(lazyGrid);
  return lazyGrid->baseVtable->inqXValsPtr(grid);
}

static const double *
cdfLazyGridInqYValsPtr(grid_t *grid)
{
  struct cdfLazyGrid *lazyGrid = (struct cdfLazyGrid *) grid;
  lock_lazy_load(lazyGrid);
  if (grid->y.vals == cdfPendingLoad) cdfLazyLoadXYVals(&lazyGrid->yValsGet, &grid->y.vals);
  unlock_lazy_load(lazyGrid);
  return lazyGrid->baseVtable->inqYValsPtr(grid);
}

static double
cdfLazyGridInqXYVal(grid_t *grid, size_t index, const struct xyValGet *valsGet, double *vals,
                    const double *(*inqValsPtr)(grid_t *gridptr))
{
  size_t size = valsGet->size;
  double v;
  if (vals == cdfPendingLoad)
    {
      // prevent full load if only first/last values get inspected
      if (index == 0 || index == size - 1)
        {
          size_t indexND[3];
          if (valsGet->ndims == 3)
            {
              indexND[0] = 0;
              indexND[1] = index / valsGet->count[2];
              indexND[2] = index % valsGet->count[2];
            }
          else if (valsGet->ndims == 2)
            {
              indexND[0] = index / grid->x.size;
              indexND[1] = index % grid->x.size;
            }
          else
            indexND[0] = index;
          cdf_get_var1_double(valsGet->datasetNCId, valsGet->varNCId, indexND, &v);
        }
      else
        {
          const double *grid_vals = inqValsPtr(grid);
          v = grid_vals[index];
        }
    }
  else if (vals)
    v = vals[index];
  else
    v = 0.0;

  return v;
}

static void
cdfLazyGridDefXVals(grid_t *grid, const double *vals)
{
  struct cdfLazyGrid *cdfGrid = (struct cdfLazyGrid *) grid;
  lock_lazy_load(cdfGrid);
  if (grid->x.vals == cdfPendingLoad) grid->x.vals = NULL;
  cdfGrid->xValsGet.datasetNCId = -1;
  cdfGrid->xValsGet.varNCId = -1;
  cdfGrid->baseVtable->defXVals(grid, vals);
  unlock_lazy_load(cdfGrid);
}

static void
cdfLazyGridDefYVals(grid_t *grid, const double *vals)
{
  struct cdfLazyGrid *cdfGrid = (struct cdfLazyGrid *) grid;
  lock_lazy_load(cdfGrid);
  if (grid->y.vals == cdfPendingLoad) grid->y.vals = NULL;
  cdfGrid->yValsGet.datasetNCId = -1;
  cdfGrid->yValsGet.varNCId = -1;
  cdfGrid->baseVtable->defYVals(grid, vals);
  unlock_lazy_load(cdfGrid);
}

static double
cdfLazyGridInqXVal(grid_t *grid, SizeType index)
{
  struct cdfLazyGrid *lazyGrid = (struct cdfLazyGrid *) grid;
  lock_lazy_load(lazyGrid);
  const double rv = cdfLazyGridInqXYVal(grid, index, &lazyGrid->xValsGet, grid->x.vals, grid->vtable->inqXValsPtr);
  unlock_lazy_load(lazyGrid);
  return rv;
}

static double
cdfLazyGridInqYVal(grid_t *grid, SizeType index)
{
  struct cdfLazyGrid *lazyGrid = (struct cdfLazyGrid *) grid;
  lock_lazy_load(lazyGrid);
  const double rv = cdfLazyGridInqXYVal(grid, index, &lazyGrid->yValsGet, grid->y.vals, grid->vtable->inqYValsPtr);
  unlock_lazy_load(lazyGrid);
  return rv;
}

static bool
cdfLazyXYValGetCompare(struct cdfLazyGrid *lazyGridRef, struct cdfLazyGrid *lazyGridTest)
{
  struct xyValGet *valsGetXRef = &lazyGridRef->xValsGet, *valsGetYRef = &lazyGridRef->yValsGet,
                  *valsGetXTest = &lazyGridTest->xValsGet, *valsGetYTest = &lazyGridTest->yValsGet;
  if (valsGetXRef->datasetNCId == -1 || valsGetXTest->datasetNCId == -1 || valsGetYRef->datasetNCId == -1
      || valsGetYTest->datasetNCId == -1)
    return lazyGridRef->baseVtable->compareXYFull(&lazyGridRef->base, &lazyGridTest->base);

  return valsGetXRef->datasetNCId != valsGetXTest->datasetNCId || valsGetXRef->varNCId != valsGetXTest->varNCId
         || valsGetYRef->datasetNCId != valsGetYTest->datasetNCId || valsGetYRef->varNCId != valsGetYTest->varNCId;
}

static bool
cdfLazyCompareXYFull(grid_t *gridRef, grid_t *gridTest)
{
  bool diff;
  struct cdfLazyGrid *lazyGridRef = (struct cdfLazyGrid *) gridRef;
  if (gridTest->vtable == &cdfLazyGridVtable)
    diff = cdfLazyXYValGetCompare(lazyGridRef, (struct cdfLazyGrid *) gridTest);
  else
    diff = lazyGridRef->baseVtable->compareXYFull(gridRef, gridTest);
  return diff;
}

static bool
cdfLazyCompareXYAO(grid_t *gridRef, grid_t *gridTest)
{
  bool diff;
  struct cdfLazyGrid *lazyGridRef = (struct cdfLazyGrid *) gridRef;
  if (gridTest->vtable == &cdfLazyGridVtable)
    diff = cdfLazyXYValGetCompare(lazyGridRef, (struct cdfLazyGrid *) gridTest);
  else
    diff = lazyGridRef->baseVtable->compareXYAO(gridRef, gridTest);
  return diff;
}

static const double *
cdfLazyGridInqXBoundsPtr(grid_t *grid)
{
  struct cdfLazyGrid *lazyGrid = (struct cdfLazyGrid *) grid;
  lock_lazy_load(lazyGrid);
  if (grid->x.bounds == cdfPendingLoad)
    {
      grid->x.bounds = (double *) Malloc((size_t) grid->nvertex * grid->size * sizeof(double));
      cdf_get_var_double(lazyGrid->xBoundsGet.datasetNCId, lazyGrid->xBoundsGet.varNCId, grid->x.bounds);
    }
  unlock_lazy_load(lazyGrid);
  return lazyGrid->baseVtable->inqXBoundsPtr(grid);
}

static void
cdfLazyGridDefXBounds(grid_t *grid, const double *xbounds)
{
  struct cdfLazyGrid *cdfGrid = (struct cdfLazyGrid *) grid;
  lock_lazy_load(cdfGrid);
  if (grid->x.bounds == cdfPendingLoad) grid->x.bounds = NULL;
  cdfGrid->xBoundsGet.datasetNCId = -1;
  cdfGrid->xBoundsGet.varNCId = -1;
  cdfGrid->baseVtable->defXBounds(grid, xbounds);
  unlock_lazy_load(cdfGrid);
}

static void
cdfLazyGridDefYBounds(grid_t *grid, const double *ybounds)
{
  struct cdfLazyGrid *cdfGrid = (struct cdfLazyGrid *) grid;
  lock_lazy_load(cdfGrid);
  if (grid->y.bounds == cdfPendingLoad) grid->y.bounds = NULL;
  cdfGrid->yBoundsGet.datasetNCId = -1;
  cdfGrid->yBoundsGet.varNCId = -1;
  cdfGrid->baseVtable->defYBounds(grid, ybounds);
  unlock_lazy_load(cdfGrid);
}

static const double *
cdfLazyGridInqYBoundsPtr(grid_t *grid)
{
  struct cdfLazyGrid *lazyGrid = (struct cdfLazyGrid *) grid;
  lock_lazy_load(lazyGrid);
  if (grid->y.bounds == cdfPendingLoad)
    {
      grid->y.bounds = (double *) Malloc((size_t) grid->nvertex * grid->size * sizeof(double));
      cdf_get_var_double(lazyGrid->yBoundsGet.datasetNCId, lazyGrid->yBoundsGet.varNCId, grid->y.bounds);
    }
  unlock_lazy_load(lazyGrid);
  return lazyGrid->baseVtable->inqYBoundsPtr(grid);
}

static void
cdfLazyGridCopyScalarFields(grid_t *gridptrOrig, grid_t *gridptrDup)
{
  struct cdfLazyGrid *lazyGridDup = (struct cdfLazyGrid *) gridptrDup, *lazyGridOrig = (struct cdfLazyGrid *) gridptrOrig;
  lazyGridOrig->baseVtable->copyScalarFields(gridptrOrig, &lazyGridDup->base);
  lazyGridDup->baseVtable = lazyGridOrig->baseVtable;
  lazyGridDup->cellAreaGet = lazyGridOrig->cellAreaGet;
  lazyGridDup->xBoundsGet = lazyGridOrig->xBoundsGet;
  lazyGridDup->yBoundsGet = lazyGridOrig->yBoundsGet;
  lazyGridDup->xValsGet = lazyGridOrig->xValsGet;
  lazyGridDup->yValsGet = lazyGridOrig->yValsGet;
  init_lazy_load_lock(lazyGridDup);
}

static void
cdfLazyGridCopyArrayFields(grid_t *gridptrOrig, grid_t *gridptrDup)
{
  const size_t reducedPointsSize = (size_t) gridptrOrig->reducedPointsSize;
  const size_t gridsize = gridptrOrig->size;
  const int gridtype = gridptrOrig->type;
  const int irregular = (gridtype == GRID_CURVILINEAR || gridtype == GRID_UNSTRUCTURED);

  if (reducedPointsSize)
    {
      gridptrDup->reducedPoints = (int *) Malloc(reducedPointsSize * sizeof(int));
      memcpy(gridptrDup->reducedPoints, gridptrOrig->reducedPoints, reducedPointsSize * sizeof(int));
    }

  if (gridptrOrig->x.vals != NULL && gridptrOrig->x.vals != cdfPendingLoad)
    {
      const size_t size = irregular ? gridsize : gridptrOrig->x.size;
      gridptrDup->x.vals = (double *) Malloc(size * sizeof(double));
      memcpy(gridptrDup->x.vals, gridptrOrig->x.vals, size * sizeof(double));
    }

  if (gridptrOrig->y.vals != NULL && gridptrOrig->y.vals != cdfPendingLoad)
    {
      const size_t size = irregular ? gridsize : gridptrOrig->y.size;
      gridptrDup->y.vals = (double *) Malloc(size * sizeof(double));
      memcpy(gridptrDup->y.vals, gridptrOrig->y.vals, size * sizeof(double));
    }

  if (gridptrOrig->x.bounds != NULL && gridptrOrig->x.bounds != cdfPendingLoad)
    {
      const size_t size = (irregular ? gridsize : gridptrOrig->x.size) * (size_t) gridptrOrig->nvertex;
      gridptrDup->x.bounds = (double *) Malloc(size * sizeof(double));
      memcpy(gridptrDup->x.bounds, gridptrOrig->x.bounds, size * sizeof(double));
    }

  if (gridptrOrig->y.bounds != NULL && gridptrOrig->y.bounds != cdfPendingLoad)
    {
      const size_t size = (irregular ? gridsize : gridptrOrig->y.size) * (size_t) gridptrOrig->nvertex;
      gridptrDup->y.bounds = (double *) Malloc(size * sizeof(double));
      memcpy(gridptrDup->y.bounds, gridptrOrig->y.bounds, size * sizeof(double));
    }

  {
    if (gridptrOrig->area != NULL && gridptrOrig->area != cdfPendingLoad)
      {
        const size_t size = gridsize;
        gridptrDup->area = (double *) Malloc(size * sizeof(double));
        memcpy(gridptrDup->area, gridptrOrig->area, size * sizeof(double));
      }
  }

  if (gridptrOrig->mask != NULL)
    {
      const size_t size = gridsize;
      gridptrDup->mask = (mask_t *) Malloc(size * sizeof(mask_t));
      memcpy(gridptrDup->mask, gridptrOrig->mask, size * sizeof(mask_t));
    }

  if (gridptrOrig->mask_gme != NULL)
    {
      const size_t size = gridsize;
      gridptrDup->mask_gme = (mask_t *) Malloc(size * sizeof(mask_t));
      memcpy(gridptrDup->mask_gme, gridptrOrig->mask_gme, size * sizeof(mask_t));
    }
}

static grid_t *
cdfLazyGridCopy(grid_t *gridptrOrig)
{
  struct cdfLazyGrid *lazyGridDup = (struct cdfLazyGrid *) Malloc(sizeof(*lazyGridDup));
  gridptrOrig->vtable->copyScalarFields(gridptrOrig, &lazyGridDup->base);
  gridptrOrig->vtable->copyArrayFields(gridptrOrig, &lazyGridDup->base);
  return &lazyGridDup->base;
}

static void
cdfLazyGridInitOnce(void)
{
  cdfLazyGridVtable = cdiGridVtable;
  cdfLazyGridVtable.destroy = cdfLazyGridDelete;
  cdfLazyGridVtable.copy = cdfLazyGridCopy;
  cdfLazyGridVtable.copyScalarFields = cdfLazyGridCopyScalarFields;
  cdfLazyGridVtable.copyArrayFields = cdfLazyGridCopyArrayFields;
  cdfLazyGridVtable.defArea = cdfLazyGridDefArea;
  cdfLazyGridVtable.inqAreaPtr = cdfLazyGridInqAreaPtr;
  cdfLazyGridVtable.inqArea = cdfLazyGridInqArea;
  cdfLazyGridVtable.inqXValsPtr = cdfLazyGridInqXValsPtr;
  cdfLazyGridVtable.inqYValsPtr = cdfLazyGridInqYValsPtr;
  cdfLazyGridVtable.inqXVal = cdfLazyGridInqXVal;
  cdfLazyGridVtable.inqYVal = cdfLazyGridInqYVal;
  cdfLazyGridVtable.defXVals = cdfLazyGridDefXVals;
  cdfLazyGridVtable.defYVals = cdfLazyGridDefYVals;
  cdfLazyGridVtable.compareXYFull = cdfLazyCompareXYFull;
  cdfLazyGridVtable.compareXYAO = cdfLazyCompareXYAO;
  cdfLazyGridVtable.defXBounds = cdfLazyGridDefXBounds;
  cdfLazyGridVtable.defYBounds = cdfLazyGridDefYBounds;
  cdfLazyGridVtable.inqXBoundsPtr = cdfLazyGridInqXBoundsPtr;
  cdfLazyGridVtable.inqYBoundsPtr = cdfLazyGridInqYBoundsPtr;
  /* create inaccessible memory area, if possible, this serves as
   * dummy value for pointers to data not yet loaded */
  /*
#ifdef HAVE_MMAP
  {
    size_t pgSize = cdiGetPageSize(false);
    static const char devZero[] = "/dev/zero";
    int fd = open(devZero, O_RDWR);
    if (fd == -1)
      SysError("Could not open %s to map anonymous memory", devZero);
    void *cdfInvalid = mmap(NULL, pgSize, PROT_NONE, MAP_PRIVATE, fd, 0);
    if (cdfInvalid == MAP_FAILED)
      SysError("Could not mmap anonymous memory");
    cdfPendingLoad = cdfInvalid;
    int rc = close(fd);
    if (rc == -1)
      SysError("Could not close %s file handle %d after mapping anonymous"
               " memory", devZero, fd);
  }
#else
  */
  cdfPendingLoad = (double *) &cdfPendingLoad;
  //#endif
  atexit(cdfLazyGridDestroyOnce);
#ifndef HAVE_LIBPTHREAD
  cdfLazyInitialized = true;
#endif
}

static void
cdfBaseGridInit(grid_t *grid, int gridtype)
{
  grid_init(grid);
  cdiGridTypeInit(grid, gridtype, 0);
}

static void
cdfLazyGridInit(struct cdfLazyGrid *grid, int gridtype)
{
#ifdef HAVE_LIBPTHREAD
  pthread_once(&cdfLazyInitialized, cdfLazyGridInitOnce);
#else
  if (!cdfLazyInitialized) cdfLazyGridInitOnce();
#endif
  cdfBaseGridInit(&grid->base, gridtype);
  grid->baseVtable = grid->base.vtable;
  grid->cellAreaGet.datasetNCId = -1;
  grid->cellAreaGet.varNCId = -1;
  grid->xValsGet.datasetNCId = -1;
  grid->xValsGet.varNCId = -1;
  grid->yValsGet.datasetNCId = -1;
  grid->yValsGet.varNCId = -1;
  grid->xBoundsGet.datasetNCId = -1;
  grid->xBoundsGet.varNCId = -1;
  grid->yBoundsGet.datasetNCId = -1;
  grid->yBoundsGet.varNCId = -1;
  grid->base.vtable = &cdfLazyGridVtable;
  init_lazy_load_lock(grid);
}

void
cdfLazyGridRenew(struct cdfLazyGrid *restrict *restrict gridpptr, int gridtype)
{
  struct cdfLazyGrid *restrict grid = *gridpptr;
  if (!grid) *gridpptr = grid = (struct cdfLazyGrid *) Malloc(sizeof(*grid));
  cdfLazyGridInit(grid, gridtype);
}

void
cdfBaseGridRenew(struct cdfLazyGrid *restrict *restrict gridpptr, int gridtype)
{
  struct cdfLazyGrid *restrict grid = *gridpptr;
  if (!grid) *gridpptr = grid = (struct cdfLazyGrid *) Malloc(sizeof(grid_t));
  cdfBaseGridInit((grid_t *) grid, gridtype);
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
