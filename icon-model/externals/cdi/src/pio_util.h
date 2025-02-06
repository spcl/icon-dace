#ifndef PIO_UTIL_
#define PIO_UTIL_

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <mpi.h>

#ifndef ERROR_H
#include "error.h"
#endif

#include "cdi.h"

#define MAXDEBUG 3

#define ddebug 0

#define debugString "#####"

void cdiAbortC_MPI(const char *caller, const char *filename, const char *functionname, int line, const char *errorString,
                   va_list ap) __attribute__((noreturn));

void cdiPioWarning(const char *caller, const char *fmt, va_list ap);

static inline int
callsToMPIAreAllowed()
{
  int init_flag = 0, finished_flag = 0;
  return MPI_Initialized(&init_flag) == MPI_SUCCESS && init_flag && MPI_Finalized(&finished_flag) == MPI_SUCCESS && !finished_flag;
}

static inline int
getMPICommWorldRank()
{
  int rank = -1;
  if (callsToMPIAreAllowed()) MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  return rank;
}

#define xdebug(fmt, ...)                                                                                                      \
  if (ddebug)                                                                                                                 \
    {                                                                                                                         \
      int rank = getMPICommWorldRank();                                                                                       \
      fprintf(stderr, "%s pe%d in %s, %s, line %d: " fmt "\n", debugString, rank, __func__, __FILE__, __LINE__, __VA_ARGS__); \
    }

#define xdebug3(fmt, ...)                                                                                     \
  if (ddebug == MAXDEBUG)                                                                                     \
    {                                                                                                         \
      int rank = getMPICommWorldRank();                                                                       \
      fprintf(stderr, "pe%d in %s, %s, line %d: " fmt "\n", rank, __func__, __FILE__, __LINE__, __VA_ARGS__); \
    }

void pcdiXMPI(int iret, const char *, int);
#define xmpi(ret)                                                        \
  do                                                                     \
    {                                                                    \
      int tmpIRet = (ret);                                               \
      if (tmpIRet != MPI_SUCCESS) pcdiXMPI(tmpIRet, __FILE__, __LINE__); \
    }                                                                    \
  while (0)

void pcdiXMPIStat(int, const char *, int, MPI_Status *);
#define xmpiStat(ret, stat)                                                    \
  do                                                                           \
    {                                                                          \
      int tmpIRet = (ret);                                                     \
      if (tmpIRet != MPI_SUCCESS) pcdiXMPIStat(ret, __FILE__, __LINE__, stat); \
    }                                                                          \
  while (0)

void pcdiXMPIStats(int, const char *, int, int, MPI_Status *);
#define xmpiStats(ret, n, stats)                                                    \
  do                                                                                \
    {                                                                               \
      int tmpIRet = (ret);                                                          \
      if (tmpIRet != MPI_SUCCESS) pcdiXMPIStats(ret, __FILE__, __LINE__, n, stats); \
    }                                                                               \
  while (0)

static inline int
sum_int(size_t n, int *a)
{
  int sum = 0;
  for (size_t i = 0; i < n; ++i) sum += a[i];
  return sum;
}

#ifndef NDEBUG
void cdiPioAssertConsistentIntVec(size_t n, const int *restrict a, MPI_Comm comm, int commRank, int rankRangeStart,
                                  int rankRangeEnd);
#else
#define cdiPioAssertConsistentIntVec(n, a, comm, r, rs, re) \
  do                                                        \
    {                                                       \
    }                                                       \
  while (0)
#endif

void printArray(const char *, const char *, const void *, int, int, const char *, const char *, int);
#define xprintArray(ps, array, n, datatype) \
  if (ddebug) printArray(debugString, ps, array, n, datatype, __func__, __FILE__, __LINE__)

#define xprintArray3(ps, array, n, datatype) \
  if (ddebug == MAXDEBUG) printArray(debugString, ps, array, n, datatype, __func__, __FILE__, __LINE__)

/**
 * @return number of dimensions
 */
int cdiPioQueryVarDims(int varShape[3], int vlistID, int varID);

/**
 * Computes simple, balanced 1D decomposition of weighted elements e.
 *
 * @param                 n number of elements to balance
 * @param[in] weightPfxSums points to array of the \a n+1 partial
 *                          weight sums, where
 *                          weightPfxSums[i]-weightPfxSums[i-1] =
 *                          weight of \$e_{i-1}\$
 * @param            nparts number of parts to generate for partition
 * @param[out]   separators pointer to array of size nparts+1,
 *                          initialized by this function such that
 *                          part i consists of elements indexed as
 *                          separators[i] to separators[i+1]-1
 */
void cdiPioDeco1D_CCP(size_t n, const size_t *restrict weightPfxSums, size_t nparts, size_t *restrict separators);

static inline size_t
cdiPioElemSizeInference(size_t varID, const int *conversion)
{
  int conv = conversion ? conversion[varID] : CDI_DATATYPE_FLT64;
  size_t elemSize;
  switch (conv)
    {
    case CDI_DATATYPE_FLT32: elemSize = sizeof(float); break;
    case CDI_DATATYPE_FLT64: elemSize = sizeof(double); break;
    default: Error("Invalid conversion specification: %d\n", conv); elemSize = (size_t) -1;
    }
  return elemSize;
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
