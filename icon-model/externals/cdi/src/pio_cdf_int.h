#ifndef PIO_CDF_INT_H
#define PIO_CDF_INT_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_PARALLEL_NC4
#include <setjmp.h>

#include "cdf_int.h"

enum
{
  CDI_PIO_COLLECTIVE_OPEN = -1,
};

struct cdiPioNcCreateLongJmpRetBuf
{
  sigjmp_buf jmpBuf;
  int openRank;
};

#if !defined TLS && defined HAVE_PTHREAD
extern pthread_key_t cdiPioCdfJmpKey;
#else
extern TLS struct cdiPioNcCreateLongJmpRetBuf *cdiPioCdfJmpBuf;
#endif

void cdiPioEnableNetCDFParAccess(void);

#endif /* HAVE_PARALLEL_NC4 */
#endif /* PIO_CDF_INT_H */

/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
