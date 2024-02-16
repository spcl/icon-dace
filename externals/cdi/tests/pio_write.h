#ifndef PIO_WRITE_H
#define PIO_WRITE_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef USE_MPI
#include <mpi.h>
#else
typedef int MPI_Comm;
#endif

#include <stdbool.h>

enum
{
  PIO_WRITE_CONFIG_CHECKSUM_BIT = 0,
  PIO_WRITE_CONFIG_CREATE_UUID_BIT,
  PIO_WRITE_CONFIG_CREATE_CURVILINEAR_GRID_BIT,
  PIO_WRITE_CONFIG_USE_DIST_GRID_BIT,
  PIO_WRITE_CONFIG_PRESET_DECOMPOSITION_BIT,
  PIO_WRITE_CONFIG_CHECKSUM_FLAG = 1 << PIO_WRITE_CONFIG_CHECKSUM_BIT,
  PIO_WRITE_CONFIG_CREATE_UUID_FLAG = 1 << PIO_WRITE_CONFIG_CREATE_UUID_BIT,
  PIO_WRITE_CONFIG_CREATE_CURVILINEAR_GRID_FLAG = 1 << PIO_WRITE_CONFIG_CREATE_CURVILINEAR_GRID_BIT,
  PIO_WRITE_CONFIG_USE_DIST_GRID_FLAG = 1 << PIO_WRITE_CONFIG_USE_DIST_GRID_BIT,
  PIO_WRITE_CONFIG_PRESET_DECOMPOSITION_FLAG = 1 << PIO_WRITE_CONFIG_PRESET_DECOMPOSITION_BIT,
};

struct model_config
{
  int nlon, nlat, nts, max_nlev, nvars;
  int filetype, datatype;
  int flags;
  int taxistype, taxisunit;
  const char *suffix, *prefix;
};

void modelRun(const struct model_config *setup, MPI_Comm comm);

#endif
