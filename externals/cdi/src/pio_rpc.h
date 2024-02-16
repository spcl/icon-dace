#ifndef PIO_RPC_H
#define PIO_RPC_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <limits.h>
#include <stdbool.h>

#include <mpi.h>
#include <yaxt.h>

enum collectorCommandTags
{
  FINALIZE,
  RESOURCES,
  STREAMOPEN,
  STREAMFLUSHCLOSE,
  STREAMCLOSE,
  STREAMDEFVLIST,
  STREAM_DEF_DECOMPOSED_VLIST,
  WRITETS,
  BLOCK_XFER,
  DIST_DATA_AGG,
  ASSERT_CONSISTENCY,
  IDXLIST_PRESET,
  IDXLIST_PRESET_SIZE,
};

#define MAXWINBUFFERSIZE ((size_t) 2048 * 1024 * 1024)

enum
{
  DATA_HEADER_FLOAT = 4,
  DATA_HEADER_DOUBLE = 8,
  numRPCFuncs = 1,
  STREAMDEFTIMESTEP = -1,
  HEADERSIZEMARKER = -numRPCFuncs - 1,
  PARTDESCMARKER = -numRPCFuncs - 2,
};
enum
{
  MINFUNCID = -numRPCFuncs,
  MAXFUNCID = -1,
  defVlistNInts = 2,
};
extern const char *const funcMap[numRPCFuncs];

struct headerSize
{
  int numDataEntries, numRPCEntries;
};

struct dataRecord
{
  int varID, nmiss;
};

union funcArgs
{
  struct
  {
    int streamID, tsID;
  } streamNewTimestep;
};

enum
{
  /* number of unsigned vector elements to pack an Xt_uid into */
  uidInts = (sizeof(Xt_uid) + sizeof(unsigned) - 1) / sizeof(unsigned),
  packUIDShift = (sizeof(unsigned) % sizeof(Xt_uid)) * CHAR_BIT,
};
/* Describes offset and ID of serialized partition descriptor.
 * partDescMarker == PARTDESCMARKER, always. */
struct partDescRecord
{
  unsigned packedUID[uidInts];
};

static inline void
packXTUID(unsigned *restrict packedUID, Xt_uid uid)
{
  for (size_t i = 0; i < uidInts; ++i)
    {
      packedUID[i] = uid & ~0U;
      uid >>= packUIDShift;
    }
}

static inline Xt_uid
unpackXTUID(const unsigned *packedUID)
{
  Xt_uid uid = 0;
  for (size_t i = uidInts; i > 0; --i) uid = (uid << packUIDShift) | packedUID[i - 1];
  return uid;
}

struct winHeaderEntry
{
  int id;
  union
  {
    struct headerSize headerSize;
    struct dataRecord dataRecord;
    union funcArgs funcArgs;
    struct partDescRecord partDesc;
  } specific;
  int offset;
};

/* round size to next multiple of factor */
static inline size_t
roundUpToMultiple(size_t size, size_t factor)
{
  return (size + factor - 1) / factor * factor;
}

enum
{
  /* align window base addresses and sizes to this value */
  PIO_WIN_ALIGN = sizeof(double),
};

struct clientBufSize
{
  size_t bufSize;
  int numDataRecords, numRPCRecords;
};

struct collSpec
{
  const Xt_idxlist *partDesc;
  int *conversion;
  int numClients;
  int numServers;
  bool sendRPCData;
};

struct clientBufSize computeClientStreamBufSize(int streamID, const struct collSpec *collector);

struct partDescPreset
{
  Xt_idxlist *lists;
  Xt_uid *uids;
  int *conversion;
};

void cdiPioDestroyPartDescPreset(size_t nRanks, size_t nVars, struct partDescPreset *deco);

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
