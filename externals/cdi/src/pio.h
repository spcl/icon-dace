#ifndef CDI_PIO_H
#define CDI_PIO_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include <mpi.h>

#include "namespace.h"
#include "pio_conf.h"

/* indices of extra keys put into cdiPioExtraNSKeys */

enum
{
  cdiPioEKFileWritingFinalize,
  cdiPioEKConf,
  cdiPioEKComms,
  cdiPioExtraNSKeysSize
};

extern int cdiPioExtraNSKeys[cdiPioExtraNSKeysSize];

static inline struct cdiPioConf *
cdiPioGetConf(void)
{
  return (struct cdiPioConf *) namespaceSwitchGet(cdiPioExtraNSKeys[cdiPioEKConf]).data;
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
