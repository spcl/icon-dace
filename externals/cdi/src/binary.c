#include "cdi.h"
#include "error.h"
#include "file.h"
#include "swap.h"
#include "binary.h"

UINT32
get_UINT32(unsigned char *x)
{
  // clang-format off
  switch (HOST_ENDIANNESS)
    {
    case CDI_BIGENDIAN:
      return (((UINT32)x[0])<<24) + (((UINT32)x[1])<<16) + (((UINT32)x[2])<< 8) + (UINT32)x[3];
    case CDI_LITTLEENDIAN:
      return (((UINT32)x[3])<<24) + (((UINT32)x[2])<<16) + (((UINT32)x[1])<< 8) + (UINT32)x[0];
    default:
      Error("Unhandled endianness %d", HOST_ENDIANNESS);
      return UINT32_C(0xFFFFFFFF);
    }
  // clang-format on
}

UINT32
get_SUINT32(unsigned char *x)
{
  // clang-format off
  switch (HOST_ENDIANNESS)
    {
    case CDI_BIGENDIAN:
      return (((UINT32)x[3])<<24) + (((UINT32)x[2])<<16) + (((UINT32)x[1])<< 8) + (UINT32)x[0];
    case CDI_LITTLEENDIAN:
      return (((UINT32)x[0])<<24) + (((UINT32)x[1])<<16) + (((UINT32)x[2])<< 8) + (UINT32)x[3];
    default:
      Error("Unhandled endianness %d", HOST_ENDIANNESS);
      return UINT32_C(0xFFFFFFFF);
    }
  // clang-format on
}

UINT64
get_UINT64(unsigned char *x)
{
  // clang-format off
  switch (HOST_ENDIANNESS)
    {
    case CDI_BIGENDIAN:
      return (((UINT64)x[0])<<56) + (((UINT64)x[1])<<48) + (((UINT64)x[2])<<40) + (((UINT64)x[3])<<32) +
             (((UINT64)x[4])<<24) + (((UINT64)x[5])<<16) + (((UINT64)x[6])<< 8) + (UINT64)x[7];
    case CDI_LITTLEENDIAN:
      return (((UINT64)x[7])<<56) + (((UINT64)x[6])<<48) + (((UINT64)x[5])<<40) + (((UINT64)x[4])<<32) +
             (((UINT64)x[3])<<24) + (((UINT64)x[2])<<16) + (((UINT64)x[1])<< 8) + (UINT64)x[0];
    default:
      Error("Unhandled endianness %d", HOST_ENDIANNESS);
      return UINT64_C(0xFFFFFFFFFFFFFFFF);
    }
  // clang-format on
}

UINT64
get_SUINT64(unsigned char *x)
{
  // clang-format off
  switch (HOST_ENDIANNESS)
    {
    case CDI_BIGENDIAN:
      return (((UINT64)x[7])<<56) + (((UINT64)x[6])<<48) + (((UINT64)x[5])<<40) + (((UINT64)x[4])<<32) +
             (((UINT64)x[3])<<24) + (((UINT64)x[2])<<16) + (((UINT64)x[1])<< 8) + (UINT64)x[0];
    case CDI_LITTLEENDIAN:
      return (((UINT64)x[0])<<56) + (((UINT64)x[1])<<48) + (((UINT64)x[2])<<40) + (((UINT64)x[3])<<32) +
             (((UINT64)x[4])<<24) + (((UINT64)x[5])<<16) + (((UINT64)x[6])<< 8) + (UINT64)x[7];
    default:
      Error("Unhandled endianness %d", HOST_ENDIANNESS);
      return UINT64_C(0xFFFFFFFFFFFFFFFF);
    }
  // clang-format on
}

size_t
binReadF77Block(int fileID, int byteswap)
{
  unsigned char f77block[4];
  size_t blocklen = 0;

  if (fileRead(fileID, f77block, 4) == 4)
    {
      blocklen = byteswap ? get_SUINT32(f77block) : get_UINT32(f77block);
    }

  return blocklen;
}

void
binWriteF77Block(int fileID, int byteswap, size_t blocksize)
{
  static const unsigned int s[4] = { 0, 8, 16, 24 };
  const unsigned long ublocksize = (unsigned long) blocksize;
  unsigned char f77block[4];

  switch (HOST_ENDIANNESS)
    {
    case CDI_BIGENDIAN:
      if (byteswap)
        {
          for (int i = 0; i <= 3; ++i) f77block[i] = (unsigned char) (ublocksize >> s[i]);
        }
      else
        {
          for (int i = 0; i <= 3; ++i) f77block[3 - i] = (unsigned char) (ublocksize >> s[i]);
        }
      break;
    case CDI_LITTLEENDIAN:
      if (byteswap)
        {
          for (int i = 0; i <= 3; ++i) f77block[3 - i] = (unsigned char) (ublocksize >> s[i]);
        }
      else
        {
          for (int i = 0; i <= 3; ++i) f77block[i] = (unsigned char) (ublocksize >> s[i]);
        }
      break;
    default: Error("Unhandled endianness %d", HOST_ENDIANNESS);
    }

  if (fileWrite(fileID, f77block, 4) != 4) Error("Write failed on %s", fileInqName(fileID));
}

int
binReadInt32(int fileID, int byteswap, size_t size, INT32 *ptr)
{
  if (sizeof(INT32) != 4) Error("Not implemented for %d byte integer!", sizeof(INT32));

  fileRead(fileID, (void *) ptr, 4 * size);
  if (byteswap) swap4byte(ptr, size);

  return 0;
}

int
binReadInt64(int fileID, int byteswap, size_t size, INT64 *ptr)
{
  if (sizeof(INT64) != 8) Error("Not implemented for %d byte integer!", sizeof(INT64));

  fileRead(fileID, (void *) ptr, 8 * size);
  if (byteswap) swap8byte(ptr, size);

  return 0;
}

int
binWriteInt32(int fileID, int byteswap, size_t size, INT32 *ptr)
{
  if (sizeof(INT32) != 4) Error("Not implemented for %d byte integer!", sizeof(INT32));

  if (byteswap) swap4byte(ptr, size);
  fileWrite(fileID, (void *) ptr, 4 * size);

  return 0;
}

int
binWriteInt64(int fileID, int byteswap, size_t size, INT64 *ptr)
{
  if (sizeof(INT64) != 8) Error("Not implemented for %d byte integer!", sizeof(INT64));

  if (byteswap) swap8byte(ptr, size);
  fileWrite(fileID, (void *) ptr, 8 * size);

  return 0;
}

int
binWriteFlt32(int fileID, int byteswap, size_t size, FLT32 *ptr)
{
  if (sizeof(FLT32) != 4) Error("Not implemented for %d byte float!", sizeof(FLT32));

  if (byteswap) swap4byte(ptr, size);
  fileWrite(fileID, (void *) ptr, 4 * size);

  return 0;
}

int
binWriteFlt64(int fileID, int byteswap, size_t size, FLT64 *ptr)
{
  if (sizeof(FLT64) != 8) Error("Not implemented for %d byte float!", sizeof(FLT64));

  if (byteswap) swap8byte(ptr, size);
  fileWrite(fileID, (void *) ptr, 8 * size);

  return 0;
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
