#ifndef CDI_UUID_H
#define CDI_UUID_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "cdi.h"

// clang-format off
#ifdef __cplusplus
extern "C" {
#endif

enum {
  uuidNumHexChars = 36,
};

static inline
int cdiUUIDIsNull(const unsigned char uuid[])
{
  int isNull = 1;
  for (size_t i = 0; i < CDI_UUID_SIZE; ++i) isNull &= (uuid[i] == 0);
  return isNull;
}

void cdiCreateUUID(unsigned char uuid[CDI_UUID_SIZE]);

int cdiUUID2Str(const unsigned char uuid[], char uuidstr[]);
int cdiStr2UUID(const char *uuidstr, unsigned char uuid[]);

#ifdef __cplusplus
}
#endif
// clang-format on

#endif
