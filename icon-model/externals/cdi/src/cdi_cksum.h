#ifndef CDI_CKSUM_H_
#define CDI_CKSUM_H_

#include <inttypes.h>

/* single checksum computation over one array */
uint32_t cdiCheckSum(int type, int count, const void *data);

/* composable check-sum computation,
 * 0. datatype,
 * 1. init,
 * 2. partial, appendable computation, and
 * 3. final checksum-computation
 */
struct cdiCheckSumState
{
  uint32_t sum;
  off_t len;
};

void cdiCheckSumRStart(struct cdiCheckSumState *state);
void cdiCheckSumRAdd(struct cdiCheckSumState *state, int type, int count, const void *data);
uint32_t cdiCheckSumRValue(struct cdiCheckSumState state);

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
