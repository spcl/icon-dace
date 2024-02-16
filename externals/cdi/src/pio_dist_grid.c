/*
 * since curvilinear and irregular (icosahedral, finite elements etc.)
 * grids can themselves represent large chunks of data, this file
 * provides the necessary code to construct decomposed grids.
 */
#if defined(HAVE_CONFIG_H)
#include "config.h"
#endif

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>
#include <yaxt.h>
#ifdef HAVE_PPM_CORE
#include <ppm/ppm_uniform_partition.h>
#endif
#include <cdi.h>
#include <cdipio.h>

#include "cdi_int.h"
#include "dmemory.h"
#include "error.h"
#include "grid.h"
#include "namespace.h"
#include "resource_handle.h"
#include "resource_unpack.h"
#include "serialize.h"

#include "pio.h"
#include "pio_cdf_int.h"
#include "pio_comm.h"
#include "pio_conf.h"
#include "pio_dist_grid.h"
#include "pio_rpc.h"
#include "pio_server.h"
#include "pio_util.h"

#ifdef HAVE_PPM_DIST_ARRAY_H
#include <ppm/dist_array.h>
#ifdef HAVE_LIBNETCDF
#include "stream_cdf_postdef.h"
#include "cdf_int.h"
#endif

static struct gridVirtTable cdiPioDistGridVtable;

enum cdiPioGDsa
{
  cdiPioGDsaXVals,
  cdiPioGDsaYVals,
  cdiPioGDsaXBounds,
  cdiPioGDsaYBounds,
  cdiPioGDsaArea,
  cdiPioGDsaMask,
  cdiPioGDsaMaskGME,
  cdiPioGDsaNum
};

enum cdiPioGDdistType
{
  cdiPioGDdtX,
  cdiPioGDdtY,
  cdiPioGDdt2D,
  cdiPioGDdtNum
};

/* table of distribution used for each sub-array: for irregular
 * grids, everything is full size, for regular grids, only area and
 * mask are on xsize*ysize arrays but x/y values and bounds are of the
 * size of their respective axis. Note that the decomposition of
 * bounds has an undecomposed nvertex dimension.
 */
static const enum cdiPioGDdistType cdiPioGridDist[2][cdiPioGDsaNum] = { {
                                                                            [cdiPioGDsaXVals] = cdiPioGDdt2D,
                                                                            [cdiPioGDsaYVals] = cdiPioGDdt2D,
                                                                            [cdiPioGDsaXBounds] = cdiPioGDdt2D,
                                                                            [cdiPioGDsaYBounds] = cdiPioGDdt2D,
                                                                            [cdiPioGDsaArea] = cdiPioGDdt2D,
                                                                            [cdiPioGDsaMask] = cdiPioGDdt2D,
                                                                            [cdiPioGDsaMaskGME] = cdiPioGDdt2D,
                                                                        },
                                                                        {
                                                                            [cdiPioGDsaXVals] = cdiPioGDdtX,
                                                                            [cdiPioGDsaYVals] = cdiPioGDdtY,
                                                                            [cdiPioGDsaXBounds] = cdiPioGDdtX,
                                                                            [cdiPioGDsaYBounds] = cdiPioGDdtY,
                                                                            [cdiPioGDsaArea] = cdiPioGDdt2D,
                                                                            [cdiPioGDsaMask] = cdiPioGDdt2D,
                                                                            [cdiPioGDsaMaskGME] = cdiPioGDdt2D,
                                                                        } };

static const size_t cdiPioDistArrayElemSize[cdiPioGDsaNum]
    = { [cdiPioGDsaXVals] = sizeof(double),         [cdiPioGDsaYVals] = sizeof(double), [cdiPioGDsaXBounds] = sizeof(double),
        [cdiPioGDsaYBounds] = sizeof(double),       [cdiPioGDsaArea] = sizeof(double),  [cdiPioGDsaMask] = sizeof(unsigned char),
        [cdiPioGDsaMaskGME] = sizeof(unsigned char) };

static const int cdiPioDistArrayCDIDt[cdiPioGDsaNum] = {
  [cdiPioGDsaXVals] = CDI_DATATYPE_FLT,     [cdiPioGDsaYVals] = CDI_DATATYPE_FLT, [cdiPioGDsaXBounds] = CDI_DATATYPE_FLT,
  [cdiPioGDsaYBounds] = CDI_DATATYPE_FLT,   [cdiPioGDsaArea] = CDI_DATATYPE_FLT,  [cdiPioGDsaMask] = CDI_DATATYPE_UCHAR,
  [cdiPioGDsaMaskGME] = CDI_DATATYPE_UCHAR,
};

/* useful bits to compute often-used grid properties for X- and Y-axis,
 * specifically needed to reconstruct increment and */
struct axisReduction
{
  double first, last;
  int ownerRankFirst, ownerRankLast;
  Xt_redist dist2scanRedist;
};

enum
{
  cdiPioGDsaMaxRank = 2
};
struct cdiPioDistGridExtraData
{
  Xt_idxlist partDesc[cdiPioGDdtNum], distList[cdiPioGDdtNum];
  Xt_xmap defXmaps[cdiPioGDdtNum];
  Xt_redist defRedists[cdiPioGDsaNum];
  const enum cdiPioGDdistType *distTypes;
  const struct gridVirtTable *baseVtable;
  struct PPM_dist_mult_array *distData;
  struct axisReduction aReduce[2];
  struct PPM_extent local_chunks[cdiPioGDsaNum * cdiPioGDsaMaxRank];
  struct PPM_global_array_desc sub_arrays[cdiPioGDsaNum];
  bool rmaEnabled;
};

#endif

#ifdef HAVE_PPM_DIST_ARRAY_H
static void cdiPioDistGridInit(grid_t *gridptr, int gridtype, int size, int xsize, int ysize, int nvertex,
                               const int (*xy_decomposition)[2], Xt_idxlist partDesc2D, Xt_idxlist partDescX, Xt_idxlist partDescY);
#endif

cdiResH
cdiPioDistGridCreate(int gridtype, int size, int xsize, int ysize, int nvertex, const int xy_decomposition[][2],
                     Xt_idxlist partDesc2D, Xt_idxlist partDescX, Xt_idxlist partDescY)
{
#ifdef HAVE_PPM_DIST_ARRAY_H
  int gridID = gridCreate(gridtype, size);
  if (gridtype != GRID_UNSTRUCTURED)
    {
      gridDefXsize(gridID, xsize);
      gridDefYsize(gridID, ysize);
    }
  gridDefNvertex(gridID, nvertex);
  grid_t *gridptr = grid_to_pointer(gridID);
  cdiPioDistGridInit(gridptr, gridtype, size, xsize, ysize, nvertex, xy_decomposition, partDesc2D, partDescX, partDescY);
  return gridID;
#else
  (void) gridtype;
  (void) size;
  (void) xsize;
  (void) ysize;
  (void) nvertex;
  (void) xy_decomposition;
  (void) partDesc2D;
  (void) partDescX;
  (void) partDescY;
  Error("PPM distributed array is needed for distributed grids");
  return CDI_UNDEFID;
#endif
}

#ifdef HAVE_PPM_DIST_ARRAY_H
static bool
cdiPioDistGridSwitchSyncMode(int gridID, int mode)
{
  bool switched = false;
  grid_t *gridptr = grid_to_pointer(gridID);
  if (gridptr->vtable == &cdiPioDistGridVtable)
    {
      struct cdiPioDistGridExtraData *extraData = (struct cdiPioDistGridExtraData *) gridptr->extraData;
      PPM_dist_mult_array_set_sync_mode(extraData->distData, (enum PPM_dma_sync_mode) mode, 0);
      bool dataDefined = gridptr->x.vals || gridptr->y.vals || gridptr->x.bounds || gridptr->y.bounds || gridptr->area
                         || gridptr->mask || gridptr->mask_gme;
      if (dataDefined) PPM_dist_mult_array_expose(extraData->distData);
      extraData->rmaEnabled = mode == PPM_dma_sync_mode_passive_target;
      switched = true;
    }
  return switched;
}
#endif

void
cdiPioDistGridEnableIndividualQueries(int gridID)
{
#ifdef HAVE_PPM_DIST_ARRAY_H
  bool switched = cdiPioDistGridSwitchSyncMode(gridID, PPM_dma_sync_mode_passive_target);
#else
  (void) gridID;
  bool switched = false;
#endif
  if (!switched) Error("called for non-distributed grid gridID=%d", gridID);
}

void
cdiPioDistGridDisableIndividualQueries(int gridID)
{
#ifdef HAVE_PPM_DIST_ARRAY_H
  bool switched = cdiPioDistGridSwitchSyncMode(gridID, PPM_dma_sync_mode_local_only);
#else
  bool switched = false;
#endif
  if (!switched) Error("called for non-distributed grid gridID=%d", gridID);
}

bool
cdiPioDistGridIndividualQueriesEnabled(int gridID)
{
  bool enabled = true;
#ifdef HAVE_PPM_DIST_ARRAY_H
  grid_t *gridptr = grid_to_pointer(gridID);
  if (gridptr->vtable == &cdiPioDistGridVtable) enabled = ((struct cdiPioDistGridExtraData *) gridptr->extraData)->rmaEnabled;
#endif
  return enabled;
}

#ifdef HAVE_PPM_DIST_ARRAY_H

#ifdef HAVE_LIBNETCDF
#ifdef HAVE_PARALLEL_NC4
/* maximal spatial rank of variable */
enum
{
  maxCdfVarRank = 3
};
struct cdfPostDefPutVarADouble
{
  int fileID, ncvarid;
  size_t start[maxCdfVarRank], count[maxCdfVarRank];
  const double *values;
};

static void
cdiPioDistGridCdfDelayedPutVarDouble(void *data)
{
  struct cdfPostDefPutVarADouble *put = (struct cdfPostDefPutVarADouble *) data;
  cdf_put_vara_double(put->fileID, put->ncvarid, put->start, put->count, put->values);
}
#endif

static void
cdiPioParCdfPostDefActionGridProp(stream_t *streamptr, int gridID, int ncvarid, enum gridPropInq gridProp,
                                  struct cdfPostDefActionList **delayed)
{
  grid_t *gridptr = grid_to_pointer(gridID);
  int fileID = streamptr->fileID;
#ifdef HAVE_PARALLEL_NC4
  int ownerRank = cdiPioStream2Owner(streamptr->self);
  if (gridptr->vtable == &cdiPioDistGridVtable && ownerRank == CDI_PIO_COLLECTIVE_OPEN)
    {
      struct cdfPostDefPutVarADouble *put = (struct cdfPostDefPutVarADouble *) Malloc(sizeof(*put));
      put->fileID = fileID;
      put->ncvarid = ncvarid;
      const struct cdiPioDistGridExtraData *extraData = (struct cdiPioDistGridExtraData *) gridptr->extraData;
      size_t saMaxRank = gridptr->type == GRID_UNSTRUCTURED ? 1 : 2;
      const struct PPM_extent(*local_chunks)[saMaxRank] = (const struct PPM_extent(*)[saMaxRank]) extraData->local_chunks;
      size_t saIdx = 0;
      switch (gridProp)
        {
        case GRID_PROP_XVALS:
          put->values = gridptr->x.vals;
          saIdx = cdiPioGDsaXVals;
          break;
        case GRID_PROP_YVALS:
          put->values = gridptr->y.vals;
          saIdx = cdiPioGDsaYVals;
          break;
        case GRID_PROP_XBOUNDS:
          put->values = gridptr->x.bounds;
          saIdx = cdiPioGDsaXBounds;
          break;
        case GRID_PROP_YBOUNDS:
          put->values = gridptr->y.bounds;
          saIdx = cdiPioGDsaYBounds;
          break;
        case GRID_PROP_AREA:
          put->values = gridptr->area;
          saIdx = cdiPioGDsaArea;
          break;
        case GRID_PROP_MASK:
        case GRID_PROP_MASK_GME: Error("unsupported key: %d", (int) gridProp); break;
        }
      size_t aRank = extraData->sub_arrays[saIdx].a_rank;
      for (size_t dim = 0; dim < aRank; ++dim)
        {
          put->start[dim] = (size_t) local_chunks[saIdx][dim].first;
          put->count[dim] = (size_t) local_chunks[saIdx][dim].size;
        }
      if (gridProp == GRID_PROP_XBOUNDS || gridProp == GRID_PROP_YBOUNDS)
        {
          /* on NetCDF-level, bounds have an extra dimension that's
           * indirectly represent in the distributed multi-array, where the
           * nvertex dimension is replaced by a contiguous MPI datatype */
          put->start[aRank] = 0;
          put->count[aRank] = (size_t) gridptr->nvertex;
        }
      struct cdfPostDefAction addendum;
      addendum.data = put;
      addendum.action = cdiPioDistGridCdfDelayedPutVarDouble;
      addendum.cleanup = (void (*)(void *))(void (*)(void)) memFree;
      *delayed = cdfPostDefActionAdd(*delayed, addendum);
    }
  else
#endif
      if (gridptr->vtable == &cdiPioDistGridVtable)
    {
      const struct cdiPioDistGridExtraData *extraData = (struct cdiPioDistGridExtraData *) gridptr->extraData;
      /*
       * size_t saMaxRank = gridptr->type == GRID_UNSTRUCTURED ? 1 : 2;
       * const struct PPM_extent (*local_chunks)[saMaxRank]
       *   = (const struct PPM_extent (*)[saMaxRank])extraData->local_chunks;
       */
      const struct PPM_global_array_desc *sub_arrays = extraData->sub_arrays;
      struct PPM_dist_mult_array *distData = extraData->distData;
      size_t saIdx = 0;
      size_t stride = (gridProp == GRID_PROP_XBOUNDS || gridProp == GRID_PROP_YBOUNDS) ? (size_t) gridptr->nvertex : 1;
      switch (gridProp)
        {
        case GRID_PROP_XVALS: saIdx = cdiPioGDsaXVals; break;
        case GRID_PROP_YVALS: saIdx = cdiPioGDsaYVals; break;
        case GRID_PROP_AREA: saIdx = cdiPioGDsaArea; break;
        case GRID_PROP_XBOUNDS: saIdx = cdiPioGDsaXBounds; break;
        case GRID_PROP_YBOUNDS: saIdx = cdiPioGDsaYBounds; break;
        case GRID_PROP_MASK:
        case GRID_PROP_MASK_GME: Error("unsupported key: %d", (int) gridProp); break;
        }
      size_t aRank = sub_arrays[saIdx].a_rank;
      double *values = (double *) Malloc(sizeof(*values) * stride * (size_t) (PPM_extents_size(aRank, sub_arrays[saIdx].rect)));
      int32_t coord[cdiPioGDsaMaxRank];
      if (aRank == 2)
        {
          size_t xsize = (size_t) gridptr->x.size, ysize = (size_t) gridptr->y.size;
          for (size_t j = 0; j < ysize; ++j)
            {
              coord[0] = (int32_t) j;
              for (size_t i = 0; i < xsize; ++i)
                {
                  coord[1] = (int32_t) i;
                  PPM_dist_mult_array_get(distData, saIdx, coord, values + (j * xsize + i) * stride);
                }
            }
        }
      else /* aRank == 1 */
        {
          size_t dimsize = (size_t) sub_arrays[saIdx].rect[0].size;
          for (size_t i = 0; i < dimsize; ++i)
            {
              coord[0] = (int32_t) i;
              PPM_dist_mult_array_get(distData, saIdx, coord, values + i * stride);
            }
        }
      cdfPostDefActionAddPutVal(delayed, fileID, ncvarid, values, cdfDelayedPutVarDeepCleanup);
    }
  else
    cdfPostDefActionGridProp(streamptr, gridID, ncvarid, gridProp, delayed);
}

#endif

/* TODO: virtualize with namespaces */
static MPI_Comm cdiPioDistGridComm = MPI_COMM_NULL;

static void
cdiPioDistGridInitOnce(void)
{
  cdiPioDistGridVtable.copyScalarFields = cdiGridVtable.copyScalarFields;
  cdiPioDistGridVtable.compareXYAO = cdiGridVtable.compareXYAO;
  cdiPioDistGridVtable.getPackSize = cdiGridVtable.getPackSize;
  cdiPioDistGridVtable.getPackSizeScalars = cdiGridVtable.getPackSizeScalars;
  cdiPioDistGridVtable.unpackScalars = cdiGridVtable.unpackScalars;
  cdiPioDistGridVtable.pack = cdiGridVtable.pack;
  cdiPioDistGridVtable.packScalars = cdiGridVtable.packScalars;
  cdiPioDistGridVtable.inqPropPresence = cdiGridVtable.inqPropPresence;
  bool inClientGroup = !commInqIsProcIO();
  MPI_Comm comm = inClientGroup ? commInqCommModel() : commInqCommColl();
  if (inClientGroup)
    {
      xmpi(MPI_Comm_dup(comm, &cdiPioDistGridComm));
      /* our copy of the client communicator is fine to be used by
       * yaxt directly, since no conflicting use can happen */
      xt_mpi_comm_mark_exclusive(cdiPioDistGridComm);
    }
  else
    cdiPioDistGridComm = comm;
#ifdef HAVE_LIBNETCDF
  namespaceSwitchSet(NSSWITCH_CDF_POSTDEFACTION_GRID_PROP, NSSW_FUNC(cdiPioParCdfPostDefActionGridProp));
#endif
}

void
cdiPioDistGridFinalizeOnce(int namespace)
{
  (void) namespace;
  bool inClientGroup = !commInqIsProcIO();
  if (inClientGroup && cdiPioDistGridComm != MPI_COMM_NULL) xmpi(MPI_Comm_free(&cdiPioDistGridComm));
}

#ifdef HAVE_LIBPTHREAD
#include <pthread.h>
static pthread_once_t cdiPioDistGridInitOnceState = PTHREAD_ONCE_INIT;
#define CDI_PIO_DIST_GRID_INIT_ONCE() pthread_once(&cdiPioDistGridInitOnceState, cdiPioDistGridInitOnce)
#else
static bool cdiPioDistGridInitOnceState = false;
#define CDI_PIO_DIST_GRID_INIT_ONCE()         \
  do                                          \
    {                                         \
      if (!cdiPioDistGridInitOnceState)       \
        {                                     \
          cdiPioDistGridInitOnce();           \
          cdiPioDistGridInitOnceState = true; \
        }                                     \
    }                                         \
  while (0)
#endif

static void
cdiPioDistGridInit(grid_t *gridptr, int gridtype, int size, int xsize, int ysize, int nvertex, const int (*xy_decomposition)[2],
                   Xt_idxlist partDesc2D, Xt_idxlist partDescX, Xt_idxlist partDescY)
{
  CDI_PIO_DIST_GRID_INIT_ONCE();
  struct cdiPioDistGridExtraData *extraData = (struct cdiPioDistGridExtraData *) (gridptr->extraData = Malloc(sizeof(*extraData)));
  extraData->rmaEnabled = false;
  struct PPM_global_array_desc *sub_arrays = extraData->sub_arrays;
  sub_arrays[cdiPioGDsaXVals].element_dt = MPI_DOUBLE;
  sub_arrays[cdiPioGDsaYVals].element_dt = MPI_DOUBLE;
  xmpi(MPI_Type_contiguous(nvertex, MPI_DOUBLE, &sub_arrays[cdiPioGDsaXBounds].element_dt));
  xmpi(MPI_Type_commit(&sub_arrays[cdiPioGDsaXBounds].element_dt));
  sub_arrays[cdiPioGDsaYBounds].element_dt = sub_arrays[cdiPioGDsaXBounds].element_dt;
  sub_arrays[cdiPioGDsaArea].element_dt = MPI_DOUBLE;
  sub_arrays[cdiPioGDsaMask].element_dt = MPI_UNSIGNED_CHAR;
  sub_arrays[cdiPioGDsaMaskGME].element_dt = MPI_UNSIGNED_CHAR;

  bool irregular = false;
  size_t saMaxRank = 0;
  switch (gridtype)
    {
    case GRID_CURVILINEAR:
      saMaxRank = 2;
      for (size_t saIdx = 0; saIdx < cdiPioGDsaNum; ++saIdx)
        {
          sub_arrays[saIdx].a_rank = 2;
          sub_arrays[saIdx].rect[0].first = 0;
          sub_arrays[saIdx].rect[0].size = ysize;
          sub_arrays[saIdx].rect[1].first = 0;
          sub_arrays[saIdx].rect[1].size = xsize;
          extraData->defRedists[saIdx] = NULL;
        }
      irregular = true;
      break;
    case GRID_UNSTRUCTURED:
      for (size_t saIdx = 0; saIdx < cdiPioGDsaNum; ++saIdx) sub_arrays[saIdx].rect[0].size = size;
      irregular = true;
      /*-fallthrough*/
    case GRID_GENERIC:
    case GRID_GAUSSIAN:
    case GRID_LONLAT:
    case GRID_SPECTRAL:
    case GRID_FOURIER:
    case GRID_GME:
    case GRID_TRAJECTORY:
    case CDI_PROJ_LCC:
    case GRID_PROJECTION:
      saMaxRank = 1;
      for (size_t saIdx = 0; saIdx < cdiPioGDsaNum; ++saIdx)
        {
          sub_arrays[saIdx].rect[0].first = 0;
          sub_arrays[saIdx].a_rank = 1;
          extraData->defRedists[saIdx] = NULL;
        }
      if (irregular) break;
      /* FIXME: make the following three 2D */
      sub_arrays[cdiPioGDsaArea].rect[0].size = size;
      sub_arrays[cdiPioGDsaMask].rect[0].size = size;
      sub_arrays[cdiPioGDsaMaskGME].rect[0].size = size;
      sub_arrays[cdiPioGDsaXVals].rect[0].size = xsize;
      sub_arrays[cdiPioGDsaYVals].rect[0].size = ysize;
      sub_arrays[cdiPioGDsaXBounds].rect[0].size = xsize;
      sub_arrays[cdiPioGDsaYBounds].rect[0].size = ysize;
      break;
    case GRID_GAUSSIAN_REDUCED:
    default: Error("Unexpected or unsupported grid type %d", gridtype);
    }

  extraData->distTypes = cdiPioGridDist[!irregular];
  for (size_t i = 0; i < cdiPioGDdtNum; ++i)
    {
      extraData->distList[i] = NULL;
      extraData->defXmaps[i] = NULL;
    }

  bool inClientGroup = !commInqIsProcIO();
  int (*rankQuery)(void) = inClientGroup ? commInqRankModel : commInqRankColl;
  int commRank = rankQuery();
  {
    struct PPM_extent(*local_chunks)[saMaxRank] = (struct PPM_extent(*)[saMaxRank]) extraData->local_chunks;
    /* improve automatic finding of partitioning if saMaxRank == 2 */
    if (!xy_decomposition)
      {
        int (*commSizeQuery)(void) = inClientGroup ? cdiPioCommInqSizeClients : commInqSizeColl;
        int numClients = commSizeQuery();
        for (size_t j = 0; j < cdiPioGDsaNum; ++j)
          {
            size_t aRank = (size_t) sub_arrays[j].a_rank;
            local_chunks[j][0] = PPM_uniform_partition(sub_arrays[j].rect[0], numClients, commRank);
            for (size_t i = 1; i < aRank; ++i) local_chunks[j][i] = sub_arrays[j].rect[i];
          }
      }
    else
      {
        /* check if hints are in range */
        static const char dimName[2] = { 'X', 'Y' }, caller[] = "cdiPioDistGridCreate";
        for (size_t aRank = 0; aRank < saMaxRank; ++aRank)
          {
            int32_t dimUB = sub_arrays[cdiPioGDsaArea].rect[aRank].first + sub_arrays[cdiPioGDsaArea].rect[aRank].size;
            if (xy_decomposition[aRank][0] < 0 || xy_decomposition[aRank][0] > dimUB)
              Errorc("decomposition hint for %c start out "
                     "of range [0,%d]\n",
                     dimName[aRank], (int) dimUB - 1);
            if (xy_decomposition[aRank][1] < 0 || (xy_decomposition[aRank][0] + xy_decomposition[aRank][1] > dimUB))
              Errorc("decomposition hint value %d for %c size out "
                     "of range [0,%d]\n",
                     (int) xy_decomposition[aRank][1], dimName[aRank], (int) dimUB - xy_decomposition[aRank][0]);
          }
        if (irregular) /* for irregular grids, all arrays are decomposed identically */
          for (size_t j = 0; j < cdiPioGDsaNum; ++j)
            for (size_t i = 0; i < saMaxRank; ++i)
              {
                local_chunks[j][i].first = (int32_t) xy_decomposition[i][0];
                local_chunks[j][i].size = (int32_t) xy_decomposition[i][1];
              }
        else
          {
            /* handle 1d-decomposed arrays first */
            for (size_t saIdx = 0; saIdx < cdiPioGDsaArea; ++saIdx)
              {
                /* for the first 4 arrays, disttype is alternatively
                 * cdiPioGDdtX==0 and cdiPioGDdtY==1 */
                int distType = extraData->distTypes[saIdx];
                local_chunks[saIdx][0].first = (int32_t) xy_decomposition[distType][0];
                local_chunks[saIdx][0].size = (int32_t) xy_decomposition[distType][1];
              }
            /* then rest of 2D-decomposed parts */
            for (size_t saIdx = cdiPioGDsaArea; saIdx < cdiPioGDsaNum; ++saIdx)
              for (size_t dim = 0; dim < 2; ++dim)
                {
                  local_chunks[saIdx][dim].first = (int32_t) xy_decomposition[dim][0];
                  local_chunks[saIdx][dim].size = (int32_t) xy_decomposition[dim][1];
                }
          }
      }
    extraData->distData = PPM_dist_mult_array_new(cdiPioGDsaNum, sub_arrays, (struct PPM_extent *) local_chunks, cdiPioDistGridComm,
                                                  0, PPM_dma_sync_mode_local_only);
  }

  extraData->partDesc[cdiPioGDdt2D] = xt_idxlist_copy(partDesc2D);
  /* TODO: force x/y part desc presence for regular grids */
  extraData->partDesc[cdiPioGDdtX] = partDescX ? xt_idxlist_copy(partDescX) : NULL;
  extraData->partDesc[cdiPioGDdtY] = partDescY ? xt_idxlist_copy(partDescY) : NULL;
  extraData->baseVtable = gridptr->vtable;
  gridptr->vtable = &cdiPioDistGridVtable;
  /* determine from where to cache xfirst,xlast,yfirst,ylast from on
   * calls to gridDef[XY]Vals */
  {
    const struct PPM_extent(*local_chunks)[saMaxRank] = (const struct PPM_extent(*)[saMaxRank]) extraData->local_chunks;
    int32_t coord[cdiPioGDsaMaxRank];
    enum
    {
      IDX_FIRST,
      IDX_LAST,
      NUM_IDX
    };
    enum
    {
      NUM_COORD = cdiPioGDsaYVals - cdiPioGDsaXVals + 1
    };
    bool haveCoord[NUM_COORD][NUM_IDX];
    for (size_t saIdx = cdiPioGDsaXVals; saIdx <= cdiPioGDsaYVals; ++saIdx)
      {
        size_t a_rank = sub_arrays[saIdx].a_rank;
        for (size_t dim = 0; dim < a_rank; ++dim) coord[dim] = sub_arrays[saIdx].rect[dim].first;
        haveCoord[saIdx][IDX_FIRST] = PPM_coord_is_contained_in_extents(a_rank, coord, local_chunks[saIdx]);
        for (size_t dim = 0; dim < a_rank; ++dim) coord[dim] += sub_arrays[cdiPioGDsaXVals].rect[dim].size - 1;
        haveCoord[saIdx][IDX_LAST] = PPM_coord_is_contained_in_extents(a_rank, coord, local_chunks[saIdx]);
      }
    {
      int ownerRanks[NUM_COORD][NUM_IDX];
      for (size_t saIdx = cdiPioGDsaXVals; saIdx <= cdiPioGDsaYVals; ++saIdx)
        for (size_t i = 0; i < NUM_IDX; ++i) ownerRanks[saIdx][i] = haveCoord[saIdx][i] ? commRank : -1;
      enum
      {
        REDUCE_COUNT = NUM_COORD * NUM_IDX
      };
      xmpi(MPI_Allreduce(MPI_IN_PLACE, ownerRanks, REDUCE_COUNT, MPI_INT, MPI_MAX, cdiPioDistGridComm));
      for (size_t saIdx = cdiPioGDsaXVals; saIdx <= cdiPioGDsaYVals; ++saIdx)
        {
          extraData->aReduce[saIdx].ownerRankFirst = ownerRanks[saIdx][IDX_FIRST];
          extraData->aReduce[saIdx].ownerRankLast = ownerRanks[saIdx][IDX_LAST];
          extraData->aReduce[saIdx].first = (double) NAN;
          extraData->aReduce[saIdx].last = (double) NAN;
          extraData->aReduce[saIdx].dist2scanRedist = NULL;
        }
    }
  }
}

#define MAX(a, b) ((a) > (b) ? (a) : (b))

static Xt_redist
cdiPioDistGridLazyRedistCreate(grid_t *gridptr, enum cdiPioGDsa saIdx)
{
  struct cdiPioDistGridExtraData *restrict extraData = (struct cdiPioDistGridExtraData *) gridptr->extraData;
  enum cdiPioGDdistType distType = extraData->distTypes[saIdx];
  Xt_idxlist distList = extraData->distList[distType];
  Xt_xmap defXmap = extraData->defXmaps[distType];
  Xt_redist defRedist = extraData->defRedists[saIdx];
  const struct PPM_global_array_desc *sub_arrays = extraData->sub_arrays;
  if (!distList)
    {
      size_t saMaxRank = gridptr->type == GRID_UNSTRUCTURED ? 1 : 2;
      const struct PPM_extent(*local_chunks)[saMaxRank] = (const struct PPM_extent(*)[saMaxRank]) extraData->local_chunks;
      if (sub_arrays[saIdx].a_rank == 1)
        {
          struct Xt_stripe stripe = { .start = local_chunks[saIdx][0].first, .stride = 1, .nstrides = local_chunks[saIdx][0].size };
          distList = xt_idxstripes_new(&stripe, 1);
        }
      else /* sub_arrays[saIdx].a_rank > 1 */
        {
          size_t aRank = sub_arrays[saIdx].a_rank;
          Xt_int global_size[cdiPioGDsaMaxRank], local_start[cdiPioGDsaMaxRank];
          int local_size[cdiPioGDsaMaxRank];
          for (size_t i = 0; i < aRank; ++i)
            {
              global_size[i] = sub_arrays[saIdx].rect[i].size;
              local_start[i] = (Xt_int) local_chunks[saIdx][i].first;
              local_size[i] = (int) local_chunks[saIdx][i].size;
            }
          distList = xt_idxsection_new(0, (int) aRank, global_size, local_size, local_start);
        }
      extraData->distList[distType] = distList;
    }
  if (!defXmap)
    {
      MPI_Comm comm = cdiPioDistGridComm;
      xmap_new_func_ptr xmapNew = cdiPioGetConf()->xmap_new;
      extraData->defXmaps[distType] = defXmap = xmapNew(extraData->partDesc[distType], distList, comm);
    }

  if (!defRedist)
    {
      switch (saIdx)
        {
        case cdiPioGDsaXVals:
        case cdiPioGDsaMask:
        case cdiPioGDsaXBounds: defRedist = xt_redist_p2p_new(defXmap, sub_arrays[saIdx].element_dt); break;
        case cdiPioGDsaYVals:
        case cdiPioGDsaArea:
          if (distType == extraData->distTypes[cdiPioGDsaXVals])
            defRedist = cdiPioDistGridLazyRedistCreate(gridptr, cdiPioGDsaXVals);
          else
            defRedist = xt_redist_p2p_new(defXmap, MPI_DOUBLE);
          break;
        case cdiPioGDsaYBounds: defRedist = cdiPioDistGridLazyRedistCreate(gridptr, cdiPioGDsaXBounds); break;
        case cdiPioGDsaMaskGME: defRedist = cdiPioDistGridLazyRedistCreate(gridptr, cdiPioGDsaMask); break;
        default: abort();
        }
      extraData->defRedists[saIdx] = defRedist;
    }
  return defRedist;
}

static void
cdiPioDistGridDestroy(grid_t *gridptr)
{
  struct cdiPioDistGridExtraData *restrict extraData = (struct cdiPioDistGridExtraData *) gridptr->extraData;
  PPM_dist_mult_array_delete(extraData->distData);
  gridptr->x.vals = NULL;
  gridptr->y.vals = NULL;
  gridptr->x.bounds = NULL;
  gridptr->y.bounds = NULL;
  gridptr->area = NULL;
  gridptr->mask = NULL;
  gridptr->mask_gme = NULL;
  if (extraData->partDesc[cdiPioGDdtX]) xt_idxlist_delete(extraData->partDesc[cdiPioGDdtX]);
  if (extraData->partDesc[cdiPioGDdtY]) xt_idxlist_delete(extraData->partDesc[cdiPioGDdtY]);
  xt_idxlist_delete(extraData->partDesc[cdiPioGDdt2D]);
  xmpi(MPI_Type_free(&extraData->sub_arrays[cdiPioGDsaXBounds].element_dt));
  for (size_t i = 0; i < cdiPioGDdtNum; ++i)
    {
      if (extraData->defXmaps[i])
        {
          Xt_xmap xmap = extraData->defXmaps[i];
          xt_xmap_delete(xmap);
        }
      if (extraData->distList[i])
        {
          xt_idxlist_delete(extraData->distList[i]);
        }
    }
  for (size_t i = 0; i < cdiPioGDsaNum; ++i)
    if (extraData->defRedists[i])
      {
        Xt_redist redist = extraData->defRedists[i];
        extraData->defRedists[i] = NULL;
        xt_redist_delete(redist);
        for (size_t j = i + 1; j < cdiPioGDsaNum; ++j)
          if (extraData->defRedists[j] == redist) extraData->defRedists[j] = NULL;
      }
  for (size_t saIdx = cdiPioGDsaXVals; saIdx <= cdiPioGDsaYVals; ++saIdx)
    if (extraData->aReduce[saIdx].dist2scanRedist) xt_redist_delete(extraData->aReduce[saIdx].dist2scanRedist);
  gridptr->extraData = NULL;
  void (*baseDestroy)(grid_t * gridptr) = extraData->baseVtable->destroy;
  Free(extraData);
  baseDestroy(gridptr);
}

static grid_t *
cdiPioDistGridCopy(grid_t *gridptr)
{
  grid_t *gridptrDup = cdiGridVtable.copy(gridptr);
  struct cdiPioDistGridExtraData *restrict extraData = (struct cdiPioDistGridExtraData *) Malloc(sizeof(*extraData));
  memcpy(extraData, gridptr->extraData, sizeof(*extraData));
  for (size_t i = 0; i < cdiPioGDdtNum; ++i)
    {
      if (extraData->partDesc[i]) extraData->partDesc[i] = xt_idxlist_copy(extraData->partDesc[i]);
      if (extraData->distList[i]) extraData->distList[i] = xt_idxlist_copy(extraData->distList[i]);
      if (extraData->defXmaps[i]) extraData->defXmaps[i] = xt_xmap_copy(extraData->defXmaps[i]);
    }
  {
    bool nodup[cdiPioGDsaNum];
    for (size_t i = 0; i < cdiPioGDsaNum; ++i) nodup[i] = true;
    for (size_t i = 0; i < cdiPioGDsaNum; ++i)
      if (extraData->defRedists[i] && nodup[i])
        {
          Xt_redist redist = extraData->defRedists[i];
          extraData->defRedists[i] = xt_redist_copy(redist);
          for (size_t j = i + 1; j < cdiPioGDsaNum; ++j)
            if (extraData->defRedists[j] == redist)
              {
                extraData->defRedists[j] = extraData->defRedists[i];
                nodup[j] = false;
              }
        }
  }
  gridptrDup->extraData = extraData;
  gridptr->vtable->copyArrayFields(gridptr, gridptrDup);
  return gridptrDup;
}

static void
cdiPioDistGridCopyArrays(grid_t *gridptrOrig, grid_t *gridptrDup)
{
  /* TODO: develop scheme for distributed reducedPoints data */
  size_t reducedPointsSize = (size_t) gridptrOrig->reducedPointsSize;
  struct cdiPioDistGridExtraData *extraDataOrig = (struct cdiPioDistGridExtraData *) gridptrOrig->extraData,
                                 *extraDataDup
                                 = (struct cdiPioDistGridExtraData *) (gridptrDup->extraData = Malloc(sizeof(*extraDataDup)));
  for (size_t i = 0; i < cdiPioGDdtNum; ++i) extraDataDup->partDesc[i] = xt_idxlist_copy(extraDataOrig->partDesc[i]);
  extraDataDup->baseVtable = extraDataOrig->baseVtable;
  struct PPM_global_array_desc *sub_arrays = extraDataDup->sub_arrays;
  size_t saMaxRank = gridptrOrig->type == GRID_UNSTRUCTURED ? 1 : 2;
  struct PPM_extent(*local_chunks)[saMaxRank] = (struct PPM_extent(*)[saMaxRank]) extraDataDup->local_chunks;
  memcpy(sub_arrays, &extraDataOrig->sub_arrays, sizeof(extraDataDup->sub_arrays));
  memcpy(local_chunks, &extraDataOrig->local_chunks, sizeof(extraDataDup->local_chunks));
  extraDataDup->distData = PPM_dist_mult_array_new(cdiPioGDsaNum, sub_arrays, (struct PPM_extent *) local_chunks,
                                                   commInqCommModel(), 0, PPM_dma_sync_mode_passive_target);

  if (reducedPointsSize)
    {
      gridptrDup->reducedPoints = (int *) Malloc(reducedPointsSize * sizeof(int));
      memcpy(gridptrDup->reducedPoints, gridptrOrig->reducedPoints, reducedPointsSize * sizeof(int));
    }

  void *saPtr[cdiPioGDsaNum] = { [cdiPioGDsaXVals] = gridptrOrig->x.vals,     [cdiPioGDsaYVals] = gridptrOrig->y.vals,
                                 [cdiPioGDsaXBounds] = gridptrOrig->x.bounds, [cdiPioGDsaYBounds] = gridptrOrig->y.bounds,
                                 [cdiPioGDsaArea] = gridptrOrig->area,        [cdiPioGDsaMask] = gridptrOrig->mask,
                                 [cdiPioGDsaMaskGME] = gridptrOrig->mask_gme };

  size_t nvertex = (size_t) gridptrOrig->nvertex;
  for (size_t saIdx = 0; saIdx < cdiPioGDsaNum; ++saIdx)
    if (saPtr[saIdx])
      {
        size_t size = (size_t) (PPM_extents_size(sub_arrays[saIdx].a_rank, local_chunks[saIdx]));
        void *src = PPM_dist_mult_array_local_ptr(extraDataOrig->distData, saIdx),
             *dest = PPM_dist_mult_array_local_ptr(extraDataDup->distData, saIdx);
        if (saIdx == cdiPioGDsaXBounds || saIdx == cdiPioGDsaYBounds) size *= nvertex;
        memcpy(dest, src, size * cdiPioDistArrayElemSize[saIdx]);
        saPtr[saIdx] = dest;
      }

  gridptrDup->x.vals = saPtr[cdiPioGDsaXVals];
  gridptrDup->y.vals = saPtr[cdiPioGDsaYVals];
  gridptrDup->x.bounds = saPtr[cdiPioGDsaXBounds];
  gridptrDup->y.bounds = saPtr[cdiPioGDsaYBounds];
  gridptrDup->area = saPtr[cdiPioGDsaArea];
  gridptrDup->mask = saPtr[cdiPioGDsaMask];
  gridptrDup->mask_gme = saPtr[cdiPioGDsaMaskGME];
}

static void *
defDistData(grid_t *gridptr, const void *indata, enum cdiPioGDsa saIdx)
{
  struct cdiPioDistGridExtraData *restrict extraData = (struct cdiPioDistGridExtraData *) gridptr->extraData;
  Xt_redist redist = cdiPioDistGridLazyRedistCreate(gridptr, saIdx);
  bool dataDefined = gridptr->x.vals || gridptr->y.vals || gridptr->x.bounds || gridptr->y.bounds || gridptr->area || gridptr->mask
                     || gridptr->mask_gme;
  if (dataDefined) PPM_dist_mult_array_unexpose(extraData->distData);
  void *redistData = PPM_dist_mult_array_local_ptr(extraData->distData, saIdx);
  xt_redist_s_exchange1(redist, indata, redistData);
  PPM_dist_mult_array_expose(extraData->distData);
  return redistData;
}

#if 0
/* TODO: develop scheme for distributed reducedPoints data */
static void
cdiPioDistGridDefReducedPoints(int gridID, int reducedPointsSize, const int reducedPoints[])
{
}


static void
distGridInqReducedPoints(int gridID, int *reducedPoints)
{
}
#endif

static int
cdiPioDistGridGetPart(grid_t *gridptr, enum cdiPioGDsa saIdx, void *out, size_t elem_size)
{
  struct cdiPioDistGridExtraData *extraData = (struct cdiPioDistGridExtraData *) gridptr->extraData;
  enum cdiPioGDdistType distType = extraData->distTypes[saIdx];
  Xt_idxlist partDesc = extraData->partDesc[distType];
  size_t num_vals = (size_t) xt_idxlist_get_num_indices(partDesc);
  if (out)
    {
      Xt_int *indices = (Xt_int *) Malloc(num_vals * sizeof(indices[0]));
      xt_idxlist_get_indices(partDesc, indices);
      struct PPM_dist_mult_array *distData = extraData->distData;
      int32_t coord[cdiPioGDsaMaxRank];
      const struct PPM_global_array_desc *sub_array_desc = extraData->sub_arrays + saIdx;
      if (sub_array_desc->a_rank == 1)
        for (size_t i = 0; i < num_vals; ++i)
          {
            coord[0] = (int32_t) indices[i];
            PPM_dist_mult_array_get(distData, saIdx, coord, (unsigned char *) out + elem_size * i);
          }
      else
        {
          int xsize = (int) sub_array_desc->rect[1].size;
          for (size_t i = 0; i < num_vals; ++i)
            {
              coord[0] = (int32_t) (indices[i] / xsize);
              coord[1] = (int32_t) (indices[i] % xsize);
              PPM_dist_mult_array_get(distData, saIdx, coord, (unsigned char *) out + elem_size * i);
            }
        }
      Free(indices);
    }
  return (int) num_vals;
}

static int
inqMaskData(grid_t *gridptr, enum cdiPioGDsa saIdx, mask_t **localMaskData, int *mask)
{
  (void) saIdx;
  int size = 0;
  if (*localMaskData)
    {
      struct cdiPioDistGridExtraData *extraData = (struct cdiPioDistGridExtraData *) gridptr->extraData;
      size = xt_idxlist_get_num_indices(extraData->partDesc[cdiPioGDdt2D]);
      if (mask)
        {
          size_t numLocalMaskVals = (size_t) size;
          mask_t *convMask = Malloc(numLocalMaskVals * sizeof(*convMask));
#ifndef NDEBUG
          int getSize = cdiPioDistGridGetPart(gridptr, saIdx, convMask, sizeof(convMask[0]));
#endif
          assert(getSize == size);
          for (size_t i = 0; i < numLocalMaskVals; ++i) mask[i] = (int) convMask[i];
          Free(convMask);
        }
    }
  return size;
}

static mask_t *
defMaskData(grid_t *gridptr, const int *mask, enum cdiPioGDsa saIdx)
{
  mask_t *maskPtr;
  if (mask)
    {
      struct cdiPioDistGridExtraData *extraData = (struct cdiPioDistGridExtraData *) gridptr->extraData;
      size_t numLocalMaskVals = (size_t) xt_idxlist_get_num_indices(extraData->partDesc[cdiPioGDdt2D]);
      mask_t *convMask = Malloc(numLocalMaskVals * sizeof(*convMask));
      for (size_t i = 0; i < numLocalMaskVals; ++i) convMask[i] = (mask_t) (mask[i] != 0);
      maskPtr = defDistData(gridptr, convMask, saIdx);
      Free(convMask);
    }
  else
    maskPtr = NULL;
  return maskPtr;
}

static void
cdiPioDistGridDefMask(grid_t *gridptr, const int *mask)
{
  gridptr->mask = defMaskData(gridptr, mask, cdiPioGDsaMask);
}

static int
cdiPioDistGridInqMask(grid_t *gridptr, int *mask)
{
  return inqMaskData(gridptr, cdiPioGDsaMask, &gridptr->mask, mask);
}

static void
cdiPioDistGridDefMaskGME(grid_t *gridptr, const int *mask)
{
  gridptr->mask_gme = defMaskData(gridptr, mask, cdiPioGDsaMaskGME);
}

static int
cdiPioDistGridInqMaskGME(grid_t *gridptr, int *mask)
{
  return inqMaskData(gridptr, cdiPioGDsaMaskGME, &gridptr->mask_gme, mask);
}

static int
cdiPioDistGridInqXVals(grid_t *gridptr, double *xvals)
{
  return gridptr->x.vals ? cdiPioDistGridGetPart(gridptr, cdiPioGDsaXVals, xvals, sizeof(xvals[0])) : 0;
}

static double
cdiPioDistGridDefAReductions(grid_t *gridptr, size_t saIdx)
{
  struct cdiPioDistGridExtraData *extraData = gridptr->extraData;
  struct PPM_dist_mult_array *distData = extraData->distData;
  const struct PPM_global_array_desc *sub_arrays = extraData->sub_arrays;
  size_t aRank = sub_arrays[saIdx].a_rank;
  struct axisReduction *aReduce = extraData->aReduce + saIdx;
  bool inClientGroup = !commInqIsProcIO();
  int (*rankQuery)(void) = inClientGroup ? commInqRankModel : commInqRankColl;
  int commRank = rankQuery();
  if (commRank == aReduce->ownerRankFirst)
    {
      int coord[cdiPioGDsaMaxRank];
      for (size_t dim = 0; dim < aRank; ++dim) coord[dim] = sub_arrays[saIdx].rect[dim].first;
      PPM_dist_mult_array_get(distData, saIdx, coord, &aReduce->first);
    }
  else
    aReduce->first = (double) NAN;
  if (commRank == aReduce->ownerRankLast)
    {
      int coord[cdiPioGDsaMaxRank];
      for (size_t dim = 0; dim < aRank; ++dim)
        coord[dim] = sub_arrays[saIdx].rect[dim].first + sub_arrays[saIdx].rect[dim].size - 1;
      PPM_dist_mult_array_get(distData, saIdx, coord, &aReduce->last);
    }
  else
    aReduce->last = (double) NAN;
  xmpi(MPI_Bcast(&aReduce->first, 1, MPI_DOUBLE, aReduce->ownerRankFirst, cdiPioDistGridComm));
  xmpi(MPI_Bcast(&aReduce->last, 1, MPI_DOUBLE, aReduce->ownerRankLast, cdiPioDistGridComm));
  Xt_int gsize = (Xt_int) PPM_extents_size(aRank, sub_arrays[saIdx].rect);
  double refInc = 0.0;
  if (gsize > 1)
    {
      int (*commSizeQuery)(void) = inClientGroup ? cdiPioCommInqSizeClients : commInqSizeColl;
      int commSize = commSizeQuery();
      struct PPM_extent gridScanChunk
          = PPM_uniform_partition((struct PPM_extent){ .first = 0, .size = (int32_t) gsize }, commSize, commRank);
      Xt_redist dist2scanRedist;
      if (!(dist2scanRedist = aReduce->dist2scanRedist))
        {
          struct Xt_stripe stripes[2];
          int numStripes = 1;
          stripes[0] = (struct Xt_stripe){ .start = gridScanChunk.first, .stride = 1, .nstrides = gridScanChunk.size };
          if (commRank == 0)
            {
              stripes[numStripes++] = (struct Xt_stripe){ .start = gsize - (Xt_int) 1, .stride = 1, .nstrides = 1 };
              int nstrides = stripes[0].nstrides;
              gridScanChunk.size = stripes[0].nstrides = nstrides > 0 ? nstrides : 1;
            }
          else
            {
              stripes[0].start = (Xt_int) (stripes[0].start - 1);
              /* add predecessor reference if there is something to
               * compare it to */
              stripes[0].nstrides += stripes[0].nstrides > 0;
              gridScanChunk.size = stripes[0].nstrides;
            }
          Xt_idxlist scanList = xt_idxstripes_new(stripes, numStripes);
          enum cdiPioGDdistType distType = extraData->distTypes[saIdx];
          Xt_idxlist distList = extraData->distList[distType];
          xmap_new_func_ptr xmapNew = cdiPioGetConf()->xmap_new;
          Xt_xmap scanXmap = xmapNew(distList, scanList, cdiPioDistGridComm);
          xt_idxlist_delete(scanList);
          dist2scanRedist = aReduce->dist2scanRedist = xt_redist_p2p_new(scanXmap, MPI_DOUBLE);
          xt_xmap_delete(scanXmap);
        }
      size_t lsize = (size_t) gridScanChunk.size;
      double *restrict scanVals = Malloc((lsize + (size_t) (commRank == 0)) * sizeof(*scanVals));
      {
        double *distVals = PPM_dist_mult_array_local_ptr(distData, saIdx);
        xt_redist_s_exchange1(dist2scanRedist, distVals, scanVals);
      }
      if (commRank == 0)
        {
          if (saIdx == cdiPioGDsaXVals)
            {
              double xval0 = scanVals[0], xvalW = scanVals[lsize];
              refInc = (fabs(xvalW) - fabs(xval0)) / (double) (gsize - 1);
            }
          else if (saIdx == cdiPioGDsaYVals)
            {
              double yval0 = scanVals[0], yval1 = scanVals[1];
              refInc = yval1 - yval0;
            }
          else
            xabort("Internal error");
        }
      else
        refInc = (double) NAN;
      xmpi(MPI_Bcast(&refInc, 1, MPI_DOUBLE, 0, cdiPioDistGridComm));
      if (lsize > 1)
        {
          if (saIdx == cdiPioGDsaXVals)
            {
              double eps = 0.01 * refInc;
              for (size_t i = 1; i < lsize; i++)
                if (fabs(fabs(scanVals[i] - scanVals[i - 1]) - refInc) > eps)
                  {
                    refInc = 0;
                    break;
                  }
            }
          else if (saIdx == cdiPioGDsaYVals)
            {
              double absYinc = fabs(refInc), eps = 0.01 * absYinc;
              for (size_t i = 2; i < lsize; ++i)
                if (fabs(fabs(scanVals[i] - scanVals[i - 1]) - absYinc) > eps)
                  {
                    refInc = 0;
                    break;
                  }
            }
        }
      Free(scanVals);
      xmpi(MPI_Allreduce(MPI_IN_PLACE, &refInc, 1, MPI_DOUBLE, MPI_MIN, cdiPioDistGridComm));
    }
  return refInc;
}

static void
cdiPioDistGridDefXVals(grid_t *gridptr, const double *xvals)
{
  gridptr->x.vals = defDistData(gridptr, xvals, cdiPioGDsaXVals);
  gridptr->y.inc = cdiPioDistGridDefAReductions(gridptr, cdiPioGDsaXVals);
}

static int
cdiPioDistGridInqYVals(grid_t *gridptr, double *yvals)
{
  return gridptr->y.vals ? cdiPioDistGridGetPart(gridptr, cdiPioGDsaYVals, yvals, sizeof(yvals[0])) : 0;
}

static void
cdiPioDistGridDefYVals(grid_t *gridptr, const double *yvals)
{
  gridptr->y.vals = defDistData(gridptr, yvals, cdiPioGDsaYVals);
  gridptr->y.inc = cdiPioDistGridDefAReductions(gridptr, cdiPioGDsaYVals);
}

static double
cdiPioDistGridInqXorYVal(grid_t *gridptr, enum cdiPioGDsa saIdx, const double *vals, int index)
{
  double val = 0;
  if (vals)
    {
      struct cdiPioDistGridExtraData *extraData = (struct cdiPioDistGridExtraData *) gridptr->extraData;
      if (index == 0)
        {
          struct axisReduction *aReduce = extraData->aReduce + saIdx;
          return aReduce->first;
        }
      else if (index == gridptr->size - 1)
        {
          struct axisReduction *aReduce = extraData->aReduce + saIdx;
          return aReduce->last;
        }
      const struct PPM_global_array_desc *sub_array_desc = extraData->sub_arrays + saIdx;
      int32_t coord[cdiPioGDsaMaxRank];
      if (sub_array_desc->a_rank == 1)
        coord[0] = index;
      else
        {
          /* 2D decomposition */
          int32_t xsize = sub_array_desc->rect[1].size;
          coord[0] = index / xsize;
          coord[1] = index % xsize;
        }
      PPM_dist_mult_array_get(extraData->distData, saIdx, coord, &val);
    }
  return val;
}

static double
cdiPioDistGridInqXVal(grid_t *gridptr, int index)
{
  return cdiPioDistGridInqXorYVal(gridptr, cdiPioGDsaXVals, gridptr->x.vals, index);
}

static double
cdiPioDistGridInqYVal(grid_t *gridptr, int index)
{
  return cdiPioDistGridInqXorYVal(gridptr, cdiPioGDsaYVals, gridptr->y.vals, index);
}

static double
cdiPioDistGridInqXInc(grid_t *gridptr)
{
  /* since for distributed grids, x.inc is precomputed on gridDefXvals,
   * the value in x.inc is directly usable */
  return gridptr->x.inc;
}

static double
cdiPioDistGridInqYInc(grid_t *gridptr)
{
  /* since for distributed grids, y.inc is precomputed on gridDefYvals,
   * the value in y.inc is directly usable */
  return gridptr->y.inc;
}

#if 0
static void
dist_grid_check_cyclic(grid_t *gridptr)
{
}

static int
grid_cmp_xvals(const grid_t *gridID1, const grid_t *gridID2);

static int
grid_cmp_yvals(const grid_t *gridID1, const grid_t *gridID2);

int distGridGenerate(const grid_t *grid)
{
}


static void
distGridCompress(int gridID)
{
}

#endif

static void
cdiPioDistGridInqArea(grid_t *gridptr, double *area)
{
  if (gridptr->area) cdiPioDistGridGetPart(gridptr, cdiPioGDsaArea, area, sizeof(area[0]));
}

static void
cdiPioDistGridDefArea(grid_t *gridptr, const double *area)
{
  gridptr->area = defDistData(gridptr, area, cdiPioGDsaArea);
}

static const double *
cdiPioDistGridInqUnavailableDblPtr(grid_t *gridptr)
{
  (void) gridptr;
  Error("gridInq{Area|[XY]vals|[XY]bounds}Ptr unavailable for"
        " distributed grids");
  return NULL;
}

static void
cdiPioDistGridDefXBounds(grid_t *gridptr, const double *xbounds)
{
  gridptr->x.bounds = defDistData(gridptr, xbounds, cdiPioGDsaXBounds);
}

/*
 * static int
 * inqBounds(grid_t *gridptr, enum cdiPioGDsa saIdx, double *bounds)
 * {
 *   struct cdiPioDistGridExtraData *extraData
 *     = (struct cdiPioDistGridExtraData *)gridptr->extraData;
 *   enum cdiPioGDdistType distType = extraData->distTypes[saIdx];
 *   Xt_idxlist partDesc = extraData->partDesc[distType];
 *   size_t num_vals = (size_t)xt_idxlist_get_num_indices(partDesc),
 *     num_vertices = (size_t)gridptr->nvertex;
 *   int *indices = (int *)Malloc(num_vals * sizeof (indices[0]));
 *   xt_idxlist_get_indices(partDesc, indices);
 *   struct PPM_dist_mult_array *distData = extraData->distData;
 *   int32_t coord[2];
 *   const struct PPM_global_array_desc *sub_array_desc
 *     = extraData->sub_arrays + saIdx;
 *   if (sub_array_desc->a_rank == 2)
 *     {
 *       int xsize = (int)sub_array_desc->rect[1].size;
 *       for (size_t j = 0; j < num_vals; ++j)
 *       {
 *         coord[0] = (int32_t)(indices[j] / xsize);
 *         coord[1] = (int32_t)(indices[j] % xsize);
 *         PPM_dist_mult_array_get(distData, saIdx,
 *                                 coord, bounds + num_vertices * j);
 *       }
 *     }
 *   else
 *     for (size_t j = 0; j < num_vals; ++j)
 *       {
 *         coord[0] = (int32_t)indices[j];
 *         PPM_dist_mult_array_get(distData, saIdx,
 *                                 coord, bounds + num_vertices * j);
 *       }
 *   free(indices);
 *   return (int)(num_vals * num_vertices);
 * }
 */

static int
cdiPioDistGridInqXBounds(grid_t *gridptr, double *xbounds)
{
  size_t vertStride = (size_t) gridptr->nvertex * sizeof(*xbounds);
  return gridptr->x.bounds ? cdiPioDistGridGetPart(gridptr, cdiPioGDsaXBounds, xbounds, vertStride) : 0;
}

static void
cdiPioDistGridDefYBounds(grid_t *gridptr, const double *ybounds)
{
  gridptr->y.bounds = defDistData(gridptr, ybounds, cdiPioGDsaYBounds);
}

static int
cdiPioDistGridInqYBounds(grid_t *gridptr, double *ybounds)
{
  size_t vertStride = (size_t) gridptr->nvertex * sizeof(*ybounds);
  return gridptr->y.bounds ? cdiPioDistGridGetPart(gridptr, cdiPioGDsaYBounds, ybounds, vertStride) : 0;
}

#if 0
void distGridPrintMask(grid_t * gridptr, int index, int opt, FILE *fp)
{
}
#endif

static int
cdiPioDistGridPackSizeArrays(grid_t *gridptr, void *context)
{
  int clientRank = commInqRankModel(), numClients = cdiPioCommInqSizeClients(), numColl = commInqSizeColl(),
      collRank = cdiPioCollRank(clientRank, numClients, numColl);
  int packBufferSize = 0;
  if (clientRank == cdiPioClientRangeStart(collRank, numClients, numColl))
    {
      /* encode present flags in 1 int and number of aggregated ranks
       * in another */
      packBufferSize += serializeGetSize(2, CDI_DATATYPE_INT, context);
      struct cdiPioDistGridExtraData *extraData = (struct cdiPioDistGridExtraData *) gridptr->extraData;
      bool presence[cdiPioGDsaNum]
          = { [cdiPioGDsaXVals] = gridptr->x.vals != NULL,     [cdiPioGDsaYVals] = gridptr->y.vals != NULL,
              [cdiPioGDsaXBounds] = gridptr->x.bounds != NULL, [cdiPioGDsaYBounds] = gridptr->y.bounds != NULL,
              [cdiPioGDsaArea] = gridptr->area != NULL,        [cdiPioGDsaMask] = gridptr->mask != NULL,
              [cdiPioGDsaMaskGME] = gridptr->mask_gme != NULL };
      int lastClientInGroup = cdiPioClientRangeStart(collRank + 1, numClients, numColl) - 1,
          numClientsInGroup = lastClientInGroup - clientRank + 1;
      struct PPM_dist_mult_array *distData = extraData->distData;
      const struct PPM_global_array_desc *sub_arrays = extraData->sub_arrays;
      size_t saMaxRank = gridptr->type == GRID_UNSTRUCTURED ? 1 : 2;
      {
        packBufferSize +=
            /* description of outermost dimension chunks for all sending
             * ranks, 2D decomposition for both, regular and irregular grids
             * plus x-/y- axis decomposition for regular grids
             */
            serializeGetSize(numClientsInGroup * 2 * (int) saMaxRank, CDI_DATATYPE_INT32, context);
      }
      int nvertex = gridptr->nvertex;
      int extraDimLen[cdiPioGDsaNum] = {
        [cdiPioGDsaXVals] = 1, [cdiPioGDsaYVals] = 1, [cdiPioGDsaXBounds] = nvertex, [cdiPioGDsaYBounds] = nvertex,
        [cdiPioGDsaArea] = 1,  [cdiPioGDsaMask] = 1,  [cdiPioGDsaMaskGME] = 1,
      };

      struct PPM_extent(*local_chunks)[saMaxRank] = (struct PPM_extent(*)[saMaxRank]) extraData->local_chunks;
      for (size_t saIdx = 0; saIdx < cdiPioGDsaNum; ++saIdx)
        if (presence[saIdx])
          {
            size_t aRank = sub_arrays[saIdx].a_rank;
            int numElem = (int) (PPM_extents_size(aRank, local_chunks[saIdx]));
            for (int rank = clientRank + 1; rank <= lastClientInGroup; ++rank)
              {
                struct PPM_extent rankChunk[cdiPioGDsaMaxRank];
                PPM_dist_mult_array_rank_rect(distData, saIdx, rank, rankChunk);
                numElem += (int) (PPM_extents_size(aRank, rankChunk));
              }
            numElem *= extraDimLen[saIdx];
            packBufferSize +=
                /* data of chunks in client group for collector */
                +serializeGetSize(numElem, cdiPioDistArrayCDIDt[saIdx], context)
                /* check-sum of chunks */
                + serializeGetSize(1, CDI_DATATYPE_UINT32, context);
          }
    }
  return packBufferSize;
}

static void
cdiPioDistGridPackArrays(grid_t *gridptr, int memberMask, void *packBuffer, int packBufferSize, int *packBufferPos, void *context)
{
  (void) memberMask;
  int clientRank = commInqRankModel(), numClients = cdiPioCommInqSizeClients(), numColl = commInqSizeColl(),
      collRank = cdiPioCollRank(clientRank, numClients, numColl);

  const void *saPtr[cdiPioGDsaNum]
      = { [cdiPioGDsaXVals] = gridptr->x.vals,     [cdiPioGDsaYVals] = gridptr->y.vals, [cdiPioGDsaXBounds] = gridptr->x.bounds,
          [cdiPioGDsaYBounds] = gridptr->y.bounds, [cdiPioGDsaArea] = gridptr->area,    [cdiPioGDsaMask] = gridptr->mask,
          [cdiPioGDsaMaskGME] = gridptr->mask_gme };
  struct cdiPioDistGridExtraData *extraData = (struct cdiPioDistGridExtraData *) gridptr->extraData;
  struct PPM_global_array_desc *sub_arrays = extraData->sub_arrays;
  int groupLeader = cdiPioClientRangeStart(collRank, numClients, numColl);
  MPI_Comm aggComm = cdiPioInqCollClientIntraComm();
  int nvertex = gridptr->nvertex;
  int extraDimLen[cdiPioGDsaNum] = {
    [cdiPioGDsaXVals] = 1, [cdiPioGDsaYVals] = 1, [cdiPioGDsaXBounds] = nvertex, [cdiPioGDsaYBounds] = nvertex,
    [cdiPioGDsaArea] = 1,  [cdiPioGDsaMask] = 1,  [cdiPioGDsaMaskGME] = 1,
  };
  int type = gridptr->type;
  size_t saMaxRank = type == GRID_UNSTRUCTURED ? 1 : 2;
  const struct PPM_extent(*local_chunks)[saMaxRank] = (const struct PPM_extent(*)[saMaxRank]) extraData->local_chunks;
  if (clientRank == groupLeader)
    {
      {
        int presentFlags = 0;
        for (size_t saIdx = 0; saIdx < cdiPioGDsaNum; ++saIdx) presentFlags |= ((saPtr[saIdx] != NULL) << saIdx);
        serializePack(&presentFlags, 1, CDI_DATATYPE_INT, packBuffer, packBufferSize, packBufferPos, context);
      }
      int lastClientInGroup = cdiPioClientRangeStart(collRank + 1, numClients, numColl) - 1,
          numClientsInGroup = lastClientInGroup - clientRank + 1;
      serializePack(&numClientsInGroup, 1, CDI_DATATYPE_INT, packBuffer, packBufferSize, packBufferPos, context);
      struct PPM_dist_mult_array *distData = extraData->distData;
      struct PPM_extent *chunks;
      {
        size_t numChunksForGroup = saMaxRank * (size_t) numClientsInGroup;
        chunks = (struct PPM_extent *) (Malloc(numChunksForGroup * sizeof(*chunks)));
        /* area is always fully decomposed (only x-axis if grid is
         * unstructured, but then saMaxRank = 1) */
        memcpy(chunks, local_chunks[cdiPioGDsaArea], saMaxRank * sizeof(*chunks));
        for (int rank = clientRank + 1, rankOfs = 1; rank <= lastClientInGroup; ++rank, ++rankOfs)
          /* area is always 2D-decomposed */
          PPM_dist_mult_array_rank_rect(distData, cdiPioGDsaArea, rank, chunks + saMaxRank * (size_t) rankOfs);
        serializePack(chunks, 2 * (int) numChunksForGroup, CDI_DATATYPE_INT32, packBuffer, packBufferSize, packBufferPos, context);
      }
      int *shardSize = (int *) (Malloc((size_t) numClientsInGroup * sizeof(*shardSize)));
      const enum cdiPioGDdistType *distTypes = extraData->distTypes;
      void *shardBuf = NULL;
      for (size_t saIdx = 0; saIdx < cdiPioGDsaNum; ++saIdx)
        if (saPtr[saIdx] != NULL)
          {
            size_t elemSize = cdiPioDistArrayElemSize[saIdx];
            size_t aRank = sub_arrays[saIdx].a_rank;
            int numElem = PPM_extents_size(aRank, local_chunks[saIdx]);
            shardSize[0] = numElem;
            int maxShardSize = 0;
            size_t distOfs = (size_t) distTypes[saIdx] & (size_t) 1;
            for (int rank = clientRank + 1; rank <= lastClientInGroup; ++rank)
              {
                size_t rankOfs = (size_t) (rank - clientRank);
                /* 0 for x and 2d decomposed data, 1 for y-decompositions */
                struct PPM_extent *rankChunk = chunks + saMaxRank * rankOfs + distOfs;
                numElem += (shardSize[rankOfs] = (int) (PPM_extents_size(aRank, rankChunk)));
                maxShardSize = MAX(shardSize[rankOfs], maxShardSize);
              }
            if (maxShardSize > 0) shardBuf = Realloc(shardBuf, elemSize * (size_t) maxShardSize * (size_t) extraDimLen[saIdx]);
            struct cdiCheckSumState crcState;
            cdiCheckSumRStart(&crcState);
            /* pack own shard */
            serializePack(saPtr[saIdx], shardSize[0] * extraDimLen[saIdx], cdiPioDistArrayCDIDt[saIdx], packBuffer, packBufferSize,
                          packBufferPos, context);
            cdiCheckSumRAdd(&crcState, cdiPioDistArrayCDIDt[saIdx], shardSize[0] * extraDimLen[saIdx], saPtr[saIdx]);
            /* recv and pack other shards in client group */
            for (int rank = clientRank + 1; rank <= lastClientInGroup; ++rank)
              {
                size_t rankOfs = (size_t) (rank - clientRank);
                xmpi(MPI_Recv(shardBuf, shardSize[rankOfs], sub_arrays[saIdx].element_dt, rank, DIST_DATA_AGG, aggComm,
                              MPI_STATUS_IGNORE));
                serializePack(shardBuf, shardSize[rankOfs] * extraDimLen[saIdx], cdiPioDistArrayCDIDt[saIdx], packBuffer,
                              packBufferSize, packBufferPos, context);
                cdiCheckSumRAdd(&crcState, cdiPioDistArrayCDIDt[saIdx], shardSize[rankOfs] * extraDimLen[saIdx], shardBuf);
              }
            uint32_t d = cdiCheckSumRValue(crcState);
            /* append check-sum of shards to packed data */
            serializePack(&d, 1, CDI_DATATYPE_UINT32, packBuffer, packBufferSize, packBufferPos, context);
          }
      Free(shardSize);
      Free(chunks);
      Free(shardBuf);
    }
  else
    for (size_t saIdx = 0; saIdx < cdiPioGDsaNum; ++saIdx)
      if (saPtr[saIdx] != NULL)
        {
          size_t aRank = sub_arrays[saIdx].a_rank;
          int numElem = PPM_extents_size(aRank, local_chunks[saIdx]);
          xmpi(MPI_Send((void *) saPtr[saIdx], numElem, sub_arrays[saIdx].element_dt, groupLeader, DIST_DATA_AGG, aggComm));
        }
}

static enum cdiApplyRet
cdiPioDistGridPackAssistFunc(int id, void *res, void *data)
{
  UNUSED(id);
  UNUSED(data);
  grid_t *gridptr = (grid_t *) res;
  int resStatus = reshGetStatus(id, NULL);
  if (gridptr->vtable == &cdiPioDistGridVtable && (resStatus & RESH_SYNC_BIT))
    {
      cdiPioDistGridPackArrays(gridptr, 0, NULL, 0, NULL, NULL);
      reshSetStatus(id, NULL, resStatus & ~RESH_SYNC_BIT);
    }
  return CDI_APPLY_GO_ON;
}

void
cdiPioDistGridPackAssist(void)
{
  cdiGridApply(cdiPioDistGridPackAssistFunc, NULL);
}

int
cdiPioDistGridUnpack(char *unpackBuffer, int unpackBufferSize, int *unpackBufferPos, int originNamespace, void *context,
                     int force_id)
{
  CDI_PIO_DIST_GRID_INIT_ONCE();
  int memberMask;
  grid_t *gridP = cdiPioDistGridVtable.unpackScalars(unpackBuffer, unpackBufferSize, unpackBufferPos, originNamespace, context,
                                                     force_id, &memberMask);
  cdiPioDistGridVtable.unpackArrays(gridP, memberMask, unpackBuffer, unpackBufferSize, unpackBufferPos, originNamespace, context);
  reshSetStatus(gridP->self, NULL, reshGetStatus(gridP->self, NULL) & ~RESH_SYNC_BIT);
  return gridP->self;
}

static void
cdiPioDistGridUnpackArrays(grid_t *gridptr, int memberMask, char *unpackBuffer, int unpackBufferSize, int *unpackBufferPos,
                           int originNamespace, void *context)
{
  UNUSED(originNamespace);
  UNUSED(memberMask);
  int presentFlags, numClientsInGroup;
  serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos, &presentFlags, 1, CDI_DATATYPE_INT, context);
  serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos, &numClientsInGroup, 1, CDI_DATATYPE_INT, context);
  int type = gridptr->type;
  bool irregular = (type == GRID_CURVILINEAR || type == GRID_UNSTRUCTURED);
  size_t saMaxRank = type == GRID_UNSTRUCTURED ? 1 : 2;
  size_t chunksPerRank = saMaxRank, numChunksForGroup = chunksPerRank * (size_t) numClientsInGroup;
  struct PPM_extent *chunks = (struct PPM_extent *) (Malloc(numChunksForGroup * sizeof(*chunks)));
  serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos, chunks, 2 * (int) numChunksForGroup, CDI_DATATYPE_INT32,
                  context);
  int nvertex = gridptr->nvertex;
  {
    Xt_idxlist partDesc2D, partDescX, partDescY;
    if (type == GRID_UNSTRUCTURED)
      {
        struct Xt_stripe *stripes = (struct Xt_stripe *) (Malloc((size_t) numClientsInGroup * sizeof(*stripes)));
        for (size_t rankOfs = 0; rankOfs < (size_t) numClientsInGroup; ++rankOfs)
          stripes[rankOfs]
              = (struct Xt_stripe){ .start = (Xt_int) chunks[rankOfs].first, .stride = 1, .nstrides = (int) chunks[rankOfs].size };
        partDesc2D = xt_idxstripes_new(stripes, numClientsInGroup);
        partDescY = partDescX = NULL;
        Free(stripes);
      }
    else if (type == GRID_CURVILINEAR)
      {
        Xt_int global_size[cdiPioGDsaMaxRank] = { gridptr->y.size, gridptr->x.size }, local_start[cdiPioGDsaMaxRank];
        int local_size[cdiPioGDsaMaxRank];
        Xt_idxlist *partDescClient2D = (Xt_idxlist *) (Malloc((size_t) numClientsInGroup * sizeof(*partDescClient2D)));
        for (size_t rankOfs = 0; rankOfs < (size_t) numClientsInGroup; ++rankOfs)
          {
            for (size_t dim = 0; dim < cdiPioGDsaMaxRank; ++dim)
              {
                local_size[dim] = chunks[cdiPioGDsaMaxRank * rankOfs + dim].size;
                local_start[dim] = (Xt_int) chunks[cdiPioGDsaMaxRank * rankOfs + dim].first;
              }
            partDescClient2D[rankOfs] = xt_idxsection_new(0, cdiPioGDsaMaxRank, global_size, local_size, local_start);
          }
        partDesc2D = xt_idxlist_collection_new(partDescClient2D, numClientsInGroup);
        partDescY = partDescX = NULL;
        for (size_t rankOfs = 0; rankOfs < (size_t) numClientsInGroup; ++rankOfs) xt_idxlist_delete(partDescClient2D[rankOfs]);
        Free(partDescClient2D);
      }
    else
      {
        Xt_int global_size[cdiPioGDsaMaxRank] = { gridptr->y.size, gridptr->x.size }, local_start[cdiPioGDsaMaxRank];
        int local_size[cdiPioGDsaMaxRank];
        /* describe 2D decomposition via section for later collection
         * construction */
        Xt_idxlist *partDescClient2D = (Xt_idxlist *) (Malloc((size_t) numClientsInGroup * sizeof(*partDescClient2D)));
        /* describe x- and y-decomposition via stripes */
        struct Xt_stripe(*stripes)[numClientsInGroup]
            = (struct Xt_stripe(*)[numClientsInGroup])(Malloc((size_t) numClientsInGroup * 2 * sizeof(*stripes)));
        for (size_t rankOfs = 0; rankOfs < (size_t) numClientsInGroup; ++rankOfs)
          {
            for (size_t dim = 0; dim < cdiPioGDsaMaxRank; ++dim)
              {
                size_t chunkDimRankOfs = cdiPioGDsaMaxRank * rankOfs + dim;
                local_size[dim] = chunks[chunkDimRankOfs].size;
                local_start[dim] = (Xt_int) chunks[chunkDimRankOfs].first;
              }
            partDescClient2D[rankOfs] = xt_idxsection_new(0, cdiPioGDsaMaxRank, global_size, local_size, local_start);
            for (size_t i = cdiPioGDdtX; i <= cdiPioGDdtY; ++i)
              stripes[i][rankOfs] = (struct Xt_stripe){ .start = (Xt_int) chunks[cdiPioGDsaMaxRank * rankOfs + i].first,
                                                        .stride = 1,
                                                        .nstrides = (int) chunks[cdiPioGDsaMaxRank * rankOfs + i].size };
          }
        partDesc2D = xt_idxlist_collection_new(partDescClient2D, numClientsInGroup);
        partDescX = xt_idxstripes_new(stripes[cdiPioGDdtX], numClientsInGroup);
        partDescY = xt_idxstripes_new(stripes[cdiPioGDdtY], numClientsInGroup);
        Free(stripes);
        Free(partDescClient2D);
      }
    cdiPioDistGridInit(gridptr, type, gridptr->size, gridptr->x.size, gridptr->y.size, nvertex, NULL, partDesc2D, partDescX,
                       partDescY);
    if (!irregular)
      {
        xt_idxlist_delete(partDescY);
        xt_idxlist_delete(partDescX);
      }
    xt_idxlist_delete(partDesc2D);
  }
  void *buf = NULL;
  int extraDimLen[cdiPioGDsaNum] = {
    [cdiPioGDsaXVals] = 1, [cdiPioGDsaYVals] = 1, [cdiPioGDsaXBounds] = nvertex, [cdiPioGDsaYBounds] = nvertex,
    [cdiPioGDsaArea] = 1,  [cdiPioGDsaMask] = 1,  [cdiPioGDsaMaskGME] = 1,
  };
  /* since we count on Fortran interoperability in other places, void *
   * will have a compatible representation to any type */
  typedef void (*saDefFuncType)(grid_t *, void *);
  static const saDefFuncType saDefFunc[cdiPioGDsaNum] = {
    [cdiPioGDsaXVals] = (saDefFuncType) cdiPioDistGridDefXVals,     [cdiPioGDsaYVals] = (saDefFuncType) cdiPioDistGridDefYVals,
    [cdiPioGDsaXBounds] = (saDefFuncType) cdiPioDistGridDefXBounds, [cdiPioGDsaYBounds] = (saDefFuncType) cdiPioDistGridDefYBounds,
    [cdiPioGDsaArea] = (saDefFuncType) cdiPioDistGridDefArea,       [cdiPioGDsaMask] = (saDefFuncType) cdiPioDistGridDefMask,
    [cdiPioGDsaMaskGME] = (saDefFuncType) cdiPioDistGridDefMaskGME,
  };
  const enum cdiPioGDdistType *distTypes = cdiPioGridDist[!irregular];
  struct cdiPioDistGridExtraData *extraData = (struct cdiPioDistGridExtraData *) gridptr->extraData;
  struct PPM_global_array_desc *sub_arrays = extraData->sub_arrays;
  for (size_t saIdx = 0; saIdx < cdiPioGDsaNum; ++saIdx)
    if (presentFlags & (1 << saIdx))
      {
        /* choose chunks to use by decomposition for saIdx */
        size_t decoOfs = (size_t) distTypes[saIdx] % (size_t) 2;
        /* sum chunk sizes aggregated on this collector */
        size_t aRank = sub_arrays[saIdx].a_rank;
        int numElem = 0;
        for (size_t rankOfs = 0; rankOfs < (size_t) numClientsInGroup; ++rankOfs)
          numElem += PPM_extents_size(aRank, chunks + decoOfs + chunksPerRank * rankOfs);
        numElem *= extraDimLen[saIdx];
        size_t elemSize = cdiPioDistArrayElemSize[saIdx];
        buf = Realloc(buf, (size_t) numElem * elemSize);
        serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos, buf, numElem, cdiPioDistArrayCDIDt[saIdx], context);
        uint32_t d = cdiCheckSum(cdiPioDistArrayCDIDt[saIdx], numElem, buf), dRef;
        serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos, &dRef, 1, CDI_DATATYPE_UINT32, context);
        xassert(d == dRef);
        saDefFunc[saIdx](gridptr, buf);
      }
  Free(chunks);
  Free(buf);
}

static inline bool
compareGridVals(size_t numVals, const double *restrict valsRef, const double *restrict valsTest)
{
  bool differ = false;
  for (size_t i = 0; i < numVals; ++i)
    if (fabs(valsTest[i] - valsRef[i]) > 1.e-10)
      {
        differ = true;
        break;
      }
  return differ;
}

static bool
cdiPioDistGridCompareXYvals(grid_t *gridRef, grid_t *gridTest)
{
  bool differ = false;
  int xsizeTest = gridTest->x.size, ysizeTest = gridTest->y.size;
  if (gridRef->vtable == &cdiPioDistGridVtable && gridTest->vtable == &cdiPioDistGridVtable)
    {
      if (!differ && xsizeTest > 0 && xsizeTest == gridRef->vtable->inqXVals(gridRef, NULL))
        {
          const double *restrict xvalsRef = gridRef->x.vals, *restrict xvalsTest = gridTest->x.vals;
          differ = compareGridVals((size_t) xsizeTest, xvalsRef, xvalsTest);
        }

      if (!differ && ysizeTest > 0 && ysizeTest == gridRef->vtable->inqYVals(gridRef, NULL))
        {
          const double *restrict yvalsRef = gridRef->y.vals, *restrict yvalsTest = gridTest->y.vals;
          differ = compareGridVals((size_t) ysizeTest, yvalsRef, yvalsTest);
        }
      int differReduce = differ;
      xmpi(MPI_Allreduce(MPI_IN_PLACE, &differReduce, 1, MPI_INT, MPI_LOR, commInqCommCalc()));
      differ = differReduce;
    }
  else
    {
      if (xsizeTest != gridRef->vtable->inqXVals(gridRef, NULL) || ysizeTest != gridRef->vtable->inqYVals(gridRef, NULL))
        {
          differ = true;
        }
      else
        {
          /* warn about grid compare of distributed vs. undistributed grid */
          Warning("efficiency: comparing data of distributed "
                  "vs. non-distributed grid!");
          double *valsRefCopy = Malloc((size_t) MAX(xsizeTest, ysizeTest) * sizeof(double));
          /* fetch X data */
          gridRef->vtable->inqXVals(gridRef, valsRefCopy);
          differ = compareGridVals((size_t) xsizeTest, valsRefCopy, gridTest->vtable->inqXValsPtr(gridTest));
          if (!differ)
            {
              gridRef->vtable->inqYVals(gridRef, valsRefCopy);
              differ = compareGridVals((size_t) ysizeTest, valsRefCopy, gridTest->vtable->inqYValsPtr(gridTest));
            }
          /* todo: handle this collectively, i.e. only compare
           * indices from chunk and allreduce result */
          Free(valsRefCopy);
        }
    }
  return differ;
}

static struct gridVirtTable cdiPioDistGridVtable = {
  .destroy = cdiPioDistGridDestroy,
  .copy = cdiPioDistGridCopy,
  .copyArrayFields = cdiPioDistGridCopyArrays,
  .defXVals = cdiPioDistGridDefXVals,
  .defYVals = cdiPioDistGridDefYVals,
  .defMask = cdiPioDistGridDefMask,
  .defMaskGME = cdiPioDistGridDefMaskGME,
  .defXBounds = cdiPioDistGridDefXBounds,
  .defYBounds = cdiPioDistGridDefYBounds,
  .defArea = cdiPioDistGridDefArea,
  .inqXVal = cdiPioDistGridInqXVal,
  .inqYVal = cdiPioDistGridInqYVal,
  .inqXVals = cdiPioDistGridInqXVals,
  .inqYVals = cdiPioDistGridInqYVals,
  .inqXValsPtr = cdiPioDistGridInqUnavailableDblPtr,
  .inqYValsPtr = cdiPioDistGridInqUnavailableDblPtr,
  .inqXInc = cdiPioDistGridInqXInc,
  .inqYInc = cdiPioDistGridInqYInc,
  .inqArea = cdiPioDistGridInqArea,
  .inqAreaPtr = cdiPioDistGridInqUnavailableDblPtr,
  .inqMask = cdiPioDistGridInqMask,
  .inqMaskGME = cdiPioDistGridInqMaskGME,
  .inqXBounds = cdiPioDistGridInqXBounds,
  .inqYBounds = cdiPioDistGridInqYBounds,
  .inqXBoundsPtr = cdiPioDistGridInqUnavailableDblPtr,
  .inqYBoundsPtr = cdiPioDistGridInqUnavailableDblPtr,
  .compareXYFull = cdiPioDistGridCompareXYvals,
  .txCode = DIST_GRID,
  .getPackSizeArrays = cdiPioDistGridPackSizeArrays,
  .packArrays = cdiPioDistGridPackArrays,
  .unpack = cdiPioDistGridUnpack,
  .unpackArrays = cdiPioDistGridUnpackArrays,
};
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
