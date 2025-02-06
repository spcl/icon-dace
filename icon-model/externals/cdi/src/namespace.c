#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE 600 /* PTHREAD_MUTEX_RECURSIVE */
#endif

#include <limits.h>
#include <stdio.h>

#include "cdi.h"
#include "cdi_int.h"
#include "dmemory.h"
#include "namespace.h"
#include "resource_handle.h"
#include "serialize.h"
#include "error.h"
#include "file.h"
#include "vlist.h"
#ifdef HAVE_LIBNETCDF
#include "cdf_int.h"
#include "stream_cdf.h"
#include "stream_cdf_postdef.h"
#endif

static unsigned nNamespaces = 1;
static int activeNamespace = 0;

#ifdef HAVE_LIBNETCDF
// clang-format off
#define CDI_NETCDF_SWITCHES                              \
  { .func = (void (*)(void)) nc__create },               \
  { .func = (void (*)(void)) cdf_def_var_serial },       \
  { .func = (void (*)(void)) cdi_nc_enddef_serial },     \
  { .func = (void (*)(void)) cdi_nc__enddef_serial },    \
  { .func = (void (*)(void)) cdfDefTimestep },           \
  { .func = (void (*)(void)) cdfDefCoordinateVars },     \
  { .func = (void (*)(void)) cdfPostDefActionGridProp }
// clang-format on
#else
#define CDI_NETCDF_SWITCHES
#endif

// clang-format off
#define defaultSwitches {                                       \
    { .func = (void (*)(void)) cdiAbortC_serial },              \
    { .func = (void (*)(void)) cdiWarning },                    \
    { .func = (void (*)(void)) serializeGetSizeInCore },        \
    { .func = (void (*)(void)) serializePackInCore },           \
    { .func = (void (*)(void)) serializeUnpackInCore },         \
    { .func = (void (*)(void)) fileOpen_serial },               \
    { .func = (void (*)(void)) fileWrite },                     \
    { .func = (void (*)(void)) fileClose_serial },              \
    { .func = (void (*)(void)) cdiStreamOpenDefaultDelegate },  \
    { .func = (void (*)(void)) cdiStreamDefVlist_ },            \
    { .func = (void (*)(void)) cdiStreamSetupVlist_ },          \
    { .func = (void (*)(void)) cdiStreamWriteVar_ },            \
    { .func = (void (*)(void)) cdiStreamWriteVarChunk_ },       \
    { .func = (void (*)(void)) 0 },                             \
    { .func = (void (*)(void)) 0 },                             \
    { .func = (void (*)(void)) cdiStreamCloseDefaultDelegate }, \
    { .func = (void (*)(void)) cdiStreamDefTimestep_ },         \
    { .func = (void (*)(void)) cdiStreamSync_ },                \
    { .func = (void (*)(void)) cdiVlistDestroy_ },              \
    CDI_NETCDF_SWITCHES                                         \
    }
// clang-format on

#if defined(SX) || defined(__cplusplus)
static const union namespaceSwitchValue defaultSwitches_[NUM_NAMESPACE_SWITCH] = defaultSwitches;
#endif

enum namespaceStatus
{
  NAMESPACE_STATUS_INUSE,
  NAMESPACE_STATUS_UNUSED,
};

static union namespaceSwitchValue initialSwitches[NUM_NAMESPACE_SWITCH] = defaultSwitches;

static struct Namespace
{
  enum namespaceStatus resStage;
  unsigned numSwitches;
  union namespaceSwitchValue *switches;
} initialNamespace = { .resStage = NAMESPACE_STATUS_INUSE, .numSwitches = NUM_NAMESPACE_SWITCH, .switches = initialSwitches };

static struct Namespace *namespaces = &initialNamespace;

static unsigned namespacesSize = 1;

#if defined(HAVE_LIBPTHREAD)
#include <pthread.h>

static pthread_once_t namespaceOnce = PTHREAD_ONCE_INIT;
static pthread_mutex_t namespaceMutex;

static void
namespaceInitialize(void)
{
  pthread_mutexattr_t ma;
  pthread_mutexattr_init(&ma);
  pthread_mutexattr_settype(&ma, PTHREAD_MUTEX_RECURSIVE);
  pthread_mutex_init(&namespaceMutex, &ma);
  pthread_mutexattr_destroy(&ma);
}

#define NAMESPACE_LOCK() pthread_mutex_lock(&namespaceMutex)
#define NAMESPACE_UNLOCK() pthread_mutex_unlock(&namespaceMutex)
#define NAMESPACE_INIT() pthread_once(&namespaceOnce, namespaceInitialize)

#else

#define NAMESPACE_INIT() \
  do                     \
    {                    \
    }                    \
  while (0)
#define NAMESPACE_LOCK()
#define NAMESPACE_UNLOCK()

#endif

enum
{
  intbits = sizeof(int) * CHAR_BIT,
  nspbits = 4,
  idxbits = intbits - nspbits,
  nspmask = (int) ((((unsigned) 1 << nspbits) - 1) << idxbits),
  idxmask = (1 << idxbits) - 1,
};

enum
{
  NUM_NAMESPACES = 1 << nspbits,
  NUM_IDX = 1 << idxbits,
};

int
namespaceIdxEncode(namespaceTuple_t tin)
{
  xassert(tin.nsp < NUM_NAMESPACES && tin.idx < NUM_IDX);
  return (tin.nsp << idxbits) + tin.idx;
}

int
namespaceIdxEncode2(int nsp, int idx)
{
  xassert(nsp < NUM_NAMESPACES && idx < NUM_IDX);
  return (nsp << idxbits) + idx;
}

namespaceTuple_t
namespaceResHDecode(int resH)
{
  namespaceTuple_t tin;

  tin.idx = resH & idxmask;
  tin.nsp = (int) (((unsigned) (resH & nspmask)) >> idxbits);

  return tin;
}

int
namespaceNew(void)
{
  int newNamespaceID = -1;
  NAMESPACE_INIT();
  NAMESPACE_LOCK();
  if (namespacesSize > nNamespaces)
    {
      /* namespace is already available and only needs reinitialization */
      for (unsigned i = 0; i < namespacesSize; ++i)
        if (namespaces[i].resStage == NAMESPACE_STATUS_UNUSED)
          {
            newNamespaceID = (int) i;
            break;
          }
    }
  else if (namespacesSize == 1)
    {
      /* make room for additional namespace */
      struct Namespace *newNameSpaces = (struct Namespace *) Malloc(((size_t) namespacesSize + 1) * sizeof(namespaces[0]));
      memcpy(newNameSpaces, namespaces, sizeof(namespaces[0]));
      namespaces = newNameSpaces;
      ++namespacesSize;
      newNamespaceID = 1;
    }
  else if (namespacesSize < NUM_NAMESPACES)
    {
      /* make room for additional namespace */
      newNamespaceID = (int) namespacesSize;
      namespaces = (struct Namespace *) Realloc(namespaces, ((size_t) namespacesSize + 1) * sizeof(namespaces[0]));
      ++namespacesSize;
    }
  else /* implicit: namespacesSize >= NUM_NAMESPACES */
    {
      NAMESPACE_UNLOCK();
      return -1;
    }
  xassert(newNamespaceID >= 0 && newNamespaceID < NUM_NAMESPACES);
  ++nNamespaces;
  namespaces[newNamespaceID].resStage = NAMESPACE_STATUS_INUSE;
  namespaces[newNamespaceID].numSwitches = NUM_NAMESPACE_SWITCH;
  enum
  {
    initialNSSWSize = sizeof(union namespaceSwitchValue) * NUM_NAMESPACE_SWITCH
  };
  namespaces[newNamespaceID].switches = (union namespaceSwitchValue *) Malloc(initialNSSWSize);
#if defined(SX) || defined(__cplusplus)
  memcpy(namespaces[newNamespaceID].switches, defaultSwitches_, initialNSSWSize);
#else
  memcpy(namespaces[newNamespaceID].switches, (union namespaceSwitchValue[NUM_NAMESPACE_SWITCH]) defaultSwitches, initialNSSWSize);
#endif
  reshListCreate(newNamespaceID);
  NAMESPACE_UNLOCK();
  return newNamespaceID;
}

void
namespaceDelete(int namespaceID)
{
  NAMESPACE_INIT();
  NAMESPACE_LOCK();
  xassert(namespaceID >= 0 && (unsigned) namespaceID < namespacesSize && nNamespaces);
  reshListDestruct(namespaceID);
  if (namespaces[namespaceID].switches != initialSwitches) Free(namespaces[namespaceID].switches);
  namespaces[namespaceID].resStage = NAMESPACE_STATUS_UNUSED;
  --nNamespaces;
  NAMESPACE_UNLOCK();
}

int
namespaceGetNumber(void)
{
  return (int) nNamespaces;
}

void
namespaceSetActive(int nId)
{
  xassert((unsigned) nId < namespacesSize && namespaces[nId].resStage != NAMESPACE_STATUS_UNUSED);
  activeNamespace = nId;
}

int
namespaceGetActive(void)
{
  return activeNamespace;
}

int
namespaceAdaptKey(int originResH, int originNamespace)
{
  if (originResH == CDI_UNDEFID) return CDI_UNDEFID;

  namespaceTuple_t tin = { .idx = originResH & idxmask, .nsp = (int) (((unsigned) (originResH & nspmask)) >> idxbits) };
  xassert(tin.nsp == originNamespace);

  int nsp = namespaceGetActive();

  return namespaceIdxEncode2(nsp, tin.idx);
}

int
namespaceAdaptKey2(int originResH)
{

  if (originResH == CDI_UNDEFID) return CDI_UNDEFID;

  namespaceTuple_t tin = { .idx = originResH & idxmask, .nsp = (int) (((unsigned) (originResH & nspmask)) >> idxbits) };
  int nsp = namespaceGetActive();
  return namespaceIdxEncode2(nsp, tin.idx);
}

void
namespaceSwitchSet(int sw, union namespaceSwitchValue value)
{
  xassert(sw > NSSWITCH_NO_SUCH_SWITCH);
  int nsp = namespaceGetActive();
  NAMESPACE_LOCK();
  if (namespaces[nsp].numSwitches <= (unsigned) sw)
    {
      if (namespaces[nsp].switches != initialSwitches)
        namespaces[nsp].switches
            = (union namespaceSwitchValue *) Realloc(namespaces[nsp].switches, ((unsigned) sw + 1) * sizeof value);
      else
        {
          void *temp = Malloc(((unsigned) sw + 1) * sizeof value);
          memcpy(temp, (void *) namespaces[nsp].switches, namespaces[nsp].numSwitches * sizeof value);
          namespaces[nsp].switches = (union namespaceSwitchValue *) temp;
        }
      namespaces[nsp].numSwitches = (unsigned) sw + 1;
    }
  namespaces[nsp].switches[sw] = value;
  NAMESPACE_UNLOCK();
}

union namespaceSwitchValue
namespaceSwitchGet(int sw)
{
  int nsp = namespaceGetActive();
  xassert(sw > NSSWITCH_NO_SUCH_SWITCH && (unsigned) sw < namespaces[nsp].numSwitches);
  NAMESPACE_LOCK();
  union namespaceSwitchValue sw_val = namespaces[nsp].switches[sw];
  NAMESPACE_UNLOCK();
  return sw_val;
}

int
cdiNamespaceSwitchNewKey(void)
{
  static unsigned reservedKeys = 0;
#if defined(HAVE_LIBPTHREAD)
  static pthread_mutex_t keyMutex = PTHREAD_MUTEX_INITIALIZER;
  pthread_mutex_lock(&keyMutex);
#endif
  if (reservedKeys >= INT_MAX - NUM_NAMESPACE_SWITCH - 1) Error("pool of available namespace switch keys exhausted!");
  int newKey = (int) (reservedKeys++) + NUM_NAMESPACE_SWITCH;
#if defined(HAVE_LIBPTHREAD)
  pthread_mutex_unlock(&keyMutex);
#endif
  return newKey;
}

void
cdiReset(void)
{
  NAMESPACE_INIT();
  NAMESPACE_LOCK();
  for (unsigned namespaceID = 0; namespaceID < namespacesSize; ++namespaceID)
    if (namespaces[namespaceID].resStage != NAMESPACE_STATUS_UNUSED) namespaceDelete((int) namespaceID);
  if (namespaces != &initialNamespace)
    {
      Free(namespaces);
      namespaces = &initialNamespace;
      namespaces[0].resStage = NAMESPACE_STATUS_UNUSED;
    }
  namespacesSize = 1;
  nNamespaces = 0;
  NAMESPACE_UNLOCK();
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
