#if defined(HAVE_CONFIG_H)
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>

#ifdef USE_MPI
#include <mpi.h>
#include <yaxt.h>
#include "error.h"
#include "pio_util.h"
#endif
#include "dmemory.h"

#if defined USE_MPI && defined HAVE_PPM_CORE
#include <ppm/ppm_uniform_partition.h>
#include <core/ppm_combinatorics.h>
#endif

#include "cdi.h"
#ifdef USE_MPI
#include "cdipio.h"
#endif
#include "pio_write_setup_grid.h"
#include "cdi_uuid.h"
#include "simple_model_helper.h"

#if USE_MPI
void
findPartition2D(int npart[2], int num_parts)
{
  const uint64_t rscale = 256;
  uint32_t *factors = NULL;
  xassert(num_parts > 0);
  int numFactors = PPM_prime_factorization_32((uint32_t) num_parts, &factors);
  /* try to distribute prime factors on dimensions to get
   * approx. 2 times as many parts in x dim than y dim */
  const uint64_t optimumRatio = rscale * 2;
  npart[0] = num_parts, npart[1] = 1;
  uint_fast32_t npart_attempt[2];
  uint64_t bestRatio = (uint64_t) num_parts * rscale, bestDiff = (uint64_t) llabs((long long) (bestRatio - optimumRatio));
  /* test all assignments of factors to dimensions, starting with
   * only one assigned to x dim (omitting 0 because that would
   * always give npart[1] > npart[0] */
  for (int assign2X = 1; assign2X <= numFactors; ++assign2X)
    {
      uint_fast32_t pattern = (UINT32_C(1) << assign2X) - 1, lastPattern = pattern << (numFactors - assign2X);
      do
        {
          npart_attempt[0] = 1;
          npart_attempt[1] = 1;
          /* loop over all factors */
          for (uint_fast32_t i = 0; i < (uint_fast32_t) numFactors; ++i)
            {
              uint_fast32_t dim_idx = (pattern >> i) & 1;
              npart_attempt[dim_idx] *= factors[i];
            }
          uint64_t ratio = ((uint64_t) npart_attempt[0] * rscale) / (uint64_t) npart_attempt[1];
          uint64_t diff = (uint64_t) llabs((long long) (ratio - optimumRatio));
          if (diff < bestDiff)
            {
              npart[0] = (int) npart_attempt[0];
              npart[1] = (int) npart_attempt[1];
              bestDiff = diff;
              bestRatio = ratio;
            }
          {
            uint_fast32_t t;
#if HAVE_DECL___BUILTIN_CTZ
            t = pattern | (pattern - 1);
            pattern = (t + 1) | (((~t & -~t) - 1) >> (__builtin_ctz((unsigned) pattern) + 1));
#else
            t = (pattern | (pattern - 1)) + 1;
            pattern = t | ((((t & -t) / (pattern & -pattern)) >> 1) - 1);
#endif
          }
        }
      while (pattern <= lastPattern);
    }
  Free(factors);
}
#endif

int
setupGrid(const struct model_config *setup, MPI_Comm comm)
{
  int nlon = setup->nlon, nlat = setup->nlat;
  int gridID;
  int rank = 0;
#if USE_MPI
  xmpi(MPI_Comm_rank(comm, &rank));
#else
  (void) comm;
#endif
  if (setup->flags & PIO_WRITE_CONFIG_USE_DIST_GRID_FLAG)
    {
#if USE_MPI
      if (setup->flags & PIO_WRITE_CONFIG_CREATE_CURVILINEAR_GRID_FLAG)
        {
          size_t gridsize = (size_t) nlon * (size_t) nlat;
          Xt_int grid_global_shape[2] = { (Xt_int) nlat, (Xt_int) nlon }, grid_local_chunk_start[2];
          int nprocs;
          xmpi(MPI_Comm_size(comm, &nprocs));
          int nprocs2D[2], grid_local_chunk_shape[2];
          findPartition2D(nprocs2D, nprocs);
          {
            struct PPM_extent set_interval_y = { 0, nlat },
                              part_interval_y = PPM_uniform_partition(set_interval_y, nprocs2D[1], rank / nprocs2D[0]);
            grid_local_chunk_shape[0] = part_interval_y.size;
            grid_local_chunk_start[0] = part_interval_y.first;
          }
          {
            struct PPM_extent set_interval_x = { 0, nlon },
                              part_interval_x = PPM_uniform_partition(set_interval_x, nprocs2D[0], rank % nprocs2D[0]);
            grid_local_chunk_shape[1] = part_interval_x.size;
            grid_local_chunk_start[1] = part_interval_x.first;
          }
          {
            Xt_idxlist partDesc2D
                = xt_idxsection_new((Xt_int) 0, 2, grid_global_shape, grid_local_chunk_shape, grid_local_chunk_start);
            const int gridDecoChunks[2][2] = { [0] = {
                [0] = (int)grid_local_chunk_start[0],
                [1] = (int)grid_local_chunk_shape[0],
              }, [1] = {
                [0] = (int)grid_local_chunk_start[1],
                [1] = (int)grid_local_chunk_shape[1],
              } };
            gridID = cdiPioDistGridCreate(GRID_CURVILINEAR, (int) gridsize, nlon, nlat, 0, gridDecoChunks, partDesc2D, NULL, NULL);
            xt_idxlist_delete(partDesc2D);
          }
          size_t grid_chunk_size = (size_t) grid_local_chunk_shape[0] * (size_t) grid_local_chunk_shape[1];
          {
            double *grid_chunk = (double *) Malloc(2 * grid_chunk_size * sizeof(*grid_chunk));
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
            const size_t xyRange[2][2] = { { (size_t) grid_local_chunk_start[1], (size_t) grid_local_chunk_start[0] },
                                           { (size_t) grid_local_chunk_shape[1], (size_t) grid_local_chunk_shape[0] } };
            computeCurvilinearChunk(grid_chunk, region, (size_t) nlon, (size_t) nlat, xyRange);
            gridDefXvals(gridID, grid_chunk + grid_chunk_size);
            gridDefYvals(gridID, grid_chunk);
            Free(grid_chunk);
          }
        }
      else
        {
          fputs("Error: distributed lat-lon grids are unsupported"
                " at this moment.\n",
                stderr);
          abort();
        }
#else
      fputs("Error: distributed grids are unsupported for"
            " non-MPI configurations.\n",
            stderr);
      abort();
#endif
    }
  else
    {
      if (setup->flags & PIO_WRITE_CONFIG_CREATE_CURVILINEAR_GRID_FLAG)
        gridID = createLocalCurvilinearGrid(nlon, nlat);
      else
        gridID = createGlobalLatLonGrid(nlon, nlat);
    }
  if (setup->flags & PIO_WRITE_CONFIG_CREATE_UUID_FLAG)
    {
      unsigned char uuid[CDI_UUID_SIZE];
      if (rank == 0) cdiCreateUUID(uuid);
#if USE_MPI
      MPI_Bcast(uuid, CDI_UUID_SIZE, MPI_UNSIGNED_CHAR, 0, comm);
#endif
      gridDefUUID(gridID, uuid);
    }
  return gridID;
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
