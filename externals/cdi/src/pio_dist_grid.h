#ifndef PIO_DIST_GRID_H
#define PIO_DIST_GRID_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_PPM_DIST_ARRAY_H
void cdiPioDistGridPackAssist(void);

int cdiPioDistGridUnpack(char *unpackBuffer, int unpackBufferSize, int *unpackBufferPos, int originNamespace, void *context,
                         int force_id);

void cdiPioDistGridFinalizeOnce(int namespace);

#endif

#endif
