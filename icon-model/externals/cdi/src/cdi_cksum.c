#include <inttypes.h>
#include <sys/types.h>
#include <stdlib.h>

#include "cdi_cksum.h"
#include "cksum.h"
#include "error.h"
#include "serialize.h"

uint32_t
cdiCheckSum(int type, int count, const void *buffer)
{
  uint32_t s = 0U;
  xassert(count >= 0);
  size_t elemSize = (size_t) serializeGetSizeInCore(1, type, NULL);
  memcrc_r_eswap(&s, (const unsigned char *) buffer, (size_t) count, elemSize);
  s = memcrc_finish(&s, (off_t) (elemSize * (size_t) count));
  return s;
}

void
cdiCheckSumRStart(struct cdiCheckSumState *state)
{
  state->sum = 0U;
  state->len = 0;
}

void
cdiCheckSumRAdd(struct cdiCheckSumState *state, int type, int count, const void *buffer)
{
  size_t elemSize = (size_t) serializeGetSizeInCore(1, type, NULL);
  memcrc_r_eswap(&state->sum, (const unsigned char *) buffer, (size_t) count, elemSize);
  state->len += (off_t) (elemSize * (size_t) count);
}

uint32_t
cdiCheckSumRValue(struct cdiCheckSumState state)
{
  return memcrc_finish(&state.sum, state.len);
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
