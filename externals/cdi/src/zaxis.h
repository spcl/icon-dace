#ifndef ZAXIS_H
#define ZAXIS_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "cdi_key.h"
#include "cdi_att.h"

// clang-format off
typedef struct
{
  double    *vals;
#ifndef USE_MPI
  char     **cvals;
  int        clength;
#endif
  double    *lbounds;
  double    *ubounds;
  double    *weights;
  int        self;
  int        scalar;
  int        type;
  int        size;
  int        direction;
  int        vctsize;
  unsigned   positive;
  double    *vct;
  cdi_keys_t keys;
  cdi_atts_t atts;
}
zaxis_t;
// clang-format on

void zaxisGetTypeDescription(int zaxisType, int *outPositive, const char **outName, const char **outLongName,
                             const char **outStdName,
                             const char **outUnit);  // The returned const char* point to static storage. Don't free or modify them.

unsigned cdiZaxisCount(void);

zaxis_t *zaxis_to_pointer(int zaxisID);

void cdiZaxisGetIndexList(unsigned numIDs, int *IDs);

int zaxisUnpack(char *unpackBuffer, int unpackBufferSize, int *unpackBufferPos, int originNamespace, void *context, int force_id);

const resOps *getZaxisOps(void);

const char *zaxisInqNamePtr(int zaxisID);

const double *zaxisInqLevelsPtr(int zaxisID);
#ifndef USE_MPI
char **zaxisInqCValsPtr(int zaxisID);
#endif
void zaxisResize(int zaxisID, int size);

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
