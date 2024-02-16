#ifndef PIO_WRITE_SETUP_GRID_H
#define PIO_WRITE_SETUP_GRID_H

#if defined(HAVE_CONFIG_H)
#include "config.h"
#endif

#ifdef USE_MPI
#include <mpi.h>
#endif

#include "pio_write.h"

#ifdef USE_MPI
void findPartition2D(int npart[2], int num_parts);
#endif

int setupGrid(const struct model_config *setup, MPI_Comm comm);

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
