#ifndef BINARY_H
#define BINARY_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <inttypes.h>

#include "dtypes.h"

#ifndef HOST_ENDIANNESS
#ifdef __cplusplus
static const uint32_t HOST_ENDIANNESS_temp[1] = { UINT32_C(0x00030201) };
#define HOST_ENDIANNESS (((const unsigned char *) HOST_ENDIANNESS_temp)[0])
#else
#define HOST_ENDIANNESS (((const unsigned char *) &(const uint32_t[1]){ UINT32_C(0x00030201) })[0])
#endif
#endif

UINT32 get_UINT32(unsigned char *x);
UINT32 get_SUINT32(unsigned char *x);
UINT64 get_UINT64(unsigned char *x);
UINT64 get_SUINT64(unsigned char *x);

size_t binReadF77Block(int fileID, int byteswap);
void binWriteF77Block(int fileID, int byteswap, size_t blocksize);

int binReadInt32(int fileID, int byteswap, size_t size, INT32 *ptr);
int binReadInt64(int fileID, int byteswap, size_t size, INT64 *ptr);

int binWriteInt32(int fileID, int byteswap, size_t size, INT32 *ptr);
int binWriteInt64(int fileID, int byteswap, size_t size, INT64 *ptr);

int binWriteFlt32(int fileID, int byteswap, size_t size, FLT32 *ptr);
int binWriteFlt64(int fileID, int byteswap, size_t size, FLT64 *ptr);

#endif /* BINARY_H */
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
