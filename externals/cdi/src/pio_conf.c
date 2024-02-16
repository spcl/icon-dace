#include <errno.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>

#include <yaxt.h>

#include "dmemory.h"
#include "error.h"
#include "resource_handle.h"

#include "cdipio.h"
#include "pio_conf.h"

static int cdiPioConfCompareP(void *conf1, void *conf2);
static void cdiPioConfDestroyP(void *conf);
static void cdiPioConfPrintP(void *cdiPioConfPtr, FILE *fp);

struct intCodeStrMap
{
  const char text[28];
  int code;
};

static size_t
mapSearchStr(const struct intCodeStrMap map[], size_t mapSize, const char *str)
{
  size_t retval = SIZE_MAX;
  for (size_t i = 0; i < mapSize; ++i)
    if (!strcmp(str, map[i].text))
      {
        retval = i;
        break;
      }
  return retval;
}

static size_t
mapSearchCode(const struct intCodeStrMap map[], size_t mapSize, int code)
{
  size_t retval = SIZE_MAX;
  for (size_t i = 0; i < mapSize; ++i)
    if (code == map[i].code)
      {
        retval = i;
        break;
      }
  return retval;
}

static const struct intCodeStrMap modeMap[] = {
  { "PIO_NONE", PIO_NONE },
  { "PIO_MPI", PIO_MPI },
  { "PIO_FPGUARD", PIO_FPGUARD },
  { "PIO_ASYNCH", PIO_ASYNCH },
  { "PIO_WRITER", PIO_WRITER },
  { "PIO_MPI_FW_ORDERED", PIO_MPI_FW_ORDERED },
  { "PIO_MPI_FW_AT_ALL", PIO_MPI_FW_AT_ALL },
  { "PIO_MPI_FW_AT_REBLOCK", PIO_MPI_FW_AT_REBLOCK },
};

int
cdiPioStr2IOMode(const char *modeStr)
{
  size_t idx = mapSearchStr(modeMap, sizeof(modeMap) / sizeof(modeMap[0]), modeStr);
  int mode = (idx != SIZE_MAX) ? modeMap[idx].code : -1;
  return mode;
}

const char *
cdiPioIOMode2Str(int IOMode)
{
  if (IOMode < PIO_MINIOMODE || IOMode > PIO_MAXIOMODE)
    return "";
  else
    return modeMap[IOMode].text;
}

static const struct intCodeStrMap roleMap[] = {
  { "PIO_ROLE_CLIENT", PIO_ROLE_CLIENT },   { "PIO_ROLE_COLLECTOR", PIO_ROLE_COLLECTOR },
  { "PIO_ROLE_WRITER", PIO_ROLE_WRITER },   { "PIO_ROLE_WRITER_COLLECTOR", PIO_ROLE_WRITER_COLLECTOR },
  { "PIO_ROLE_FPGUARD", PIO_ROLE_FPGUARD },
};

#if 0
static int
cdiPioStr2CSRole(const char *roleStr)
{
  size_t idx = mapSearchStr(roleMap, sizeof (roleMap) / sizeof (roleMap[0]),
                            roleStr);
  int role = (idx != SIZE_MAX) ? roleMap[idx].code : -1;
  return role;
}
#endif

static const char *
cdiPioCSRole2Str(int role)
{
  size_t pos = mapSearchCode(roleMap, sizeof(roleMap) / sizeof(roleMap[0]), role);
  const char *roleStr = (pos != SIZE_MAX) ? roleMap[pos].text : "";
  return roleStr;
}

const resOps cdiPioConfOps = { cdiPioConfCompareP, cdiPioConfDestroyP, cdiPioConfPrintP,
                               /* serialization of configuration is not supported */
                               0, 0, 0 };

static int
cdiPioConfCompareP(void *p1, void *p2)
{
  struct cdiPioConf *a = p1, *b = p2;
  bool callBackDifference = false;
  for (size_t i = 0; i < CDIPIO_NUM_CALLBACKS; ++i) callBackDifference |= (a->callbacks[i] != b->callbacks[i]);
  return (a->IOMode != b->IOMode) | (a->clientServerRole != b->clientServerRole) | (a->partInflate != b->partInflate)
         | callBackDifference;
}

static void
cdiPioConfDestroyP(void *conf)
{
  Free(conf);
}

static void
cdiPioConfPrintP(void *cdiPioConfPtr, FILE *fp)
{
  struct cdiPioConf *conf = cdiPioConfPtr;
  const char *iomodeStr = cdiPioIOMode2Str(conf->IOMode), *CSRoleStr = cdiPioCSRole2Str(conf->clientServerRole);
  if (!iomodeStr[0]) iomodeStr = "(invalid!)";
  if (!CSRoleStr[0]) CSRoleStr = "(invalid!)";
  fprintf(fp,
          "configuration object %p\n"
          "IOMode = %s\n"
          "client/server = %s\n"
          "part data imbalance = %f\n"
          "aligning of block buffers to large pages is %sabled\n"
          "record aggregation buffer size %zu\n"
          "write aggregation buffer size %zu\n"
          "stripe conversion of index lists is %sabled\n"
          "caching of YAXT redists is %sabled\n"
          "caching of YAXT xmaps is %sabled\n"
          "callback after setup of communication = %p\n",
          cdiPioConfPtr, iomodeStr, CSRoleStr, conf->partInflate, conf->largePageAlign ? "en" : "dis",
          (size_t) conf->recordAggBufLimMB * 1024 * 1024, conf->writeAggBufLim, conf->stripify ? "en" : "dis",
          conf->cacheRedists ? "en" : "dis", conf->cacheXmaps ? "en" : "dis",
          (void *) conf->callbacks[CDIPIO_CALLBACK_POSTCOMMSETUP]);
}

static inline size_t
findWriteAccumBufsize()
{
  size_t buffersize = (size_t) 16 * 1024 * 1024;
  const char *p = getenv("BUFSIZE");
  if (p)
    {
      char *end;
      long temp = strtol(p, &end, 0);
      if ((errno == ERANGE && (temp == LONG_MAX || temp == LONG_MIN)) || (errno != 0 && temp == 0))
        {
          perror("failed to parse BUFSIZE environment variable");
        }
      else if (temp > 0 && (unsigned long) temp > buffersize)
        buffersize = (unsigned long) temp;
    }
  return buffersize;
}

int
cdiPioConfCreate(void)
{
  struct cdiPioConf *conf = Malloc(sizeof(*conf));
  conf->IOMode = PIO_NONE;
  conf->clientServerRole = PIO_ROLE_CLIENT;
  conf->partInflate = 1.1f;
  for (size_t i = 0; i < CDIPIO_NUM_CALLBACKS; ++i) conf->callbacks[i] = cdiPioNoPostCommSetup;
  conf->largePageAlign = false;
  conf->cacheRedists = conf->cacheXmaps = true;
  conf->recordAggBufLimMB = 128;
  conf->writeAggBufLim = findWriteAccumBufsize();
  conf->maxPathLen = 2 * (size_t) PATH_MAX;
  conf->aioQueueDepth = 4U;
  conf->xmap_new = xt_xmap_dist_dir_new;
  conf->stripify = true;
  conf->batchedRMA = true;
  int resH = reshPut(conf, &cdiPioConfOps);
  /* configuration objects are never forwarded */
  reshSetStatus(resH, &cdiPioConfOps, reshGetStatus(resH, &cdiPioConfOps) & ~RESH_SYNC_BIT);
  return resH;
}

void
cdiPioConfDestroy(int confResH)
{
  struct cdiPioConf *conf = reshGetVal(confResH, &cdiPioConfOps);
  cdiPioConfDestroyP(conf);
  reshRemove(confResH, &cdiPioConfOps);
}

void
cdiPioConfSetPartInflate(int confResH, float partInflate)
{
  struct cdiPioConf *conf = reshGetVal(confResH, &cdiPioConfOps);
  xassert(partInflate >= 1.0f);
  conf->partInflate = partInflate;
}

float
cdiPioConfGetPartInflate(int confResH)
{
  struct cdiPioConf *conf = reshGetVal(confResH, &cdiPioConfOps);
  return conf->partInflate;
}

void
cdiPioConfSetIOMode(int confResH, int IOMode)
{
  struct cdiPioConf *conf = reshGetVal(confResH, &cdiPioConfOps);
  conf->IOMode = IOMode;
}

int
cdiPioConfGetIOMode(int confResH)
{
  struct cdiPioConf *conf = reshGetVal(confResH, &cdiPioConfOps);
  return conf->IOMode;
}

void
cdiPioConfSetCSRole(int confResH, int CSRole)
{
  struct cdiPioConf *conf = reshGetVal(confResH, &cdiPioConfOps);
  conf->clientServerRole = CSRole;
}

int
cdiPioConfGetCSRole(int confResH)
{
  struct cdiPioConf *conf = reshGetVal(confResH, &cdiPioConfOps);
  return conf->clientServerRole;
}

void
cdiPioConfSetCallBackActions(int confResH, int trigger, void (*callback)(void))
{
  struct cdiPioConf *conf = reshGetVal(confResH, &cdiPioConfOps);
  if (trigger >= 0 && trigger < CDIPIO_NUM_CALLBACKS)
    {
      conf->callbacks[trigger] = callback;
    }
  else
    Error("invalid trigger callback query: %d", trigger);
}

void (*cdiPioConfGetCallBackActions(int confResH, int trigger))(void)
{
  void (*callback)(void) = (void (*)(void)) 0;
  struct cdiPioConf *conf = reshGetVal(confResH, &cdiPioConfOps);
  if (trigger >= 0 && trigger < CDIPIO_NUM_CALLBACKS)
    {
      callback = conf->callbacks[trigger];
    }
  else
    Error("invalid trigger callback query: %d", trigger);
  return callback;
}

void
cdiPioConfSetPostCommSetupActions(int confResH, void (*postCommSetupActions)(void))
{
  cdiPioConfSetCallBackActions(confResH, CDIPIO_CALLBACK_POSTCOMMSETUP, postCommSetupActions);
}

void (*cdiPioConfGetPostCommSetupActions(int confResH))(void)
{
  return cdiPioConfGetCallBackActions(confResH, CDIPIO_CALLBACK_POSTCOMMSETUP);
}

void
cdiPioConfSetLargePageAlign(int confResH, int largePageAlign)
{
  struct cdiPioConf *conf = reshGetVal(confResH, &cdiPioConfOps);
  conf->largePageAlign = largePageAlign != 0;
}

int
cdiPioConfGetLargePageAlign(int confResH)
{
  struct cdiPioConf *conf = reshGetVal(confResH, &cdiPioConfOps);
  return conf->largePageAlign;
}

void
cdiPioConfSetRedistCache(int confResH, int doCache)
{
  struct cdiPioConf *conf = reshGetVal(confResH, &cdiPioConfOps);
  conf->cacheRedists = doCache;
}

int
cdiPioConfGetRedistCache(int confResH)
{
  struct cdiPioConf *conf = reshGetVal(confResH, &cdiPioConfOps);
  return conf->cacheRedists;
}

void
cdiPioConfSetXmapCache(int confResH, int doCache)
{
  struct cdiPioConf *conf = reshGetVal(confResH, &cdiPioConfOps);
  conf->cacheXmaps = doCache;
}

int
cdiPioConfGetXmapCache(int confResH)
{
  struct cdiPioConf *conf = reshGetVal(confResH, &cdiPioConfOps);
  return conf->cacheXmaps;
}

void
cdiPioConfSetRecordAggBufLim(int confResH, int lim_mb)
{
  struct cdiPioConf *conf = reshGetVal(confResH, &cdiPioConfOps);
  if (lim_mb > 0)
    conf->recordAggBufLimMB = (unsigned) lim_mb;
  else
    Error("unexpected negative buffer size value %d requested", lim_mb);
}

int
cdiPioConfGetRecordAggBufLim(int confResH)
{
  struct cdiPioConf *conf = reshGetVal(confResH, &cdiPioConfOps);
  return (int) conf->recordAggBufLimMB;
}

void
cdiPioConfSetWriteAggBufLim(int confResH, int lim_mb)
{
  struct cdiPioConf *conf = reshGetVal(confResH, &cdiPioConfOps);
  if (lim_mb > 0)
    conf->writeAggBufLim = (size_t) lim_mb * 1024 * 1024;
  else
    Error("unexpected negative buffer size value %d requested", lim_mb);
}

int
cdiPioConfGetWriteAggBufLim(int confResH)
{
  struct cdiPioConf *conf = reshGetVal(confResH, &cdiPioConfOps);
  return (int) (conf->recordAggBufLimMB / (1024 * 1024));
}

void
cdiPioConfSetAioQueueDepth(int confResH, int queue_depth)
{
  struct cdiPioConf *conf = reshGetVal(confResH, &cdiPioConfOps);
  if (queue_depth >= 1 && queue_depth <= USHRT_MAX)
    conf->aioQueueDepth = (unsigned short) queue_depth;
  else
    Error("out of range AIO queue size value %d requested (MIN: 1, MAX: %d)", queue_depth, USHRT_MAX);
}

int
cdiPioConfGetAioQueueDepth(int confResH)
{
  struct cdiPioConf *conf = reshGetVal(confResH, &cdiPioConfOps);
  return conf->aioQueueDepth;
}

void
cdiPioConfSetMaxPathLen(int confResH, int max_path_len)
{
  struct cdiPioConf *conf = reshGetVal(confResH, &cdiPioConfOps);
  if (max_path_len >= 1)
    conf->maxPathLen = max_path_len;
  else
    Error("negative value %d for maximum path length requested", max_path_len);
}

int
cdiPioConfGetMaxPathLen(int confResH)
{
  struct cdiPioConf *conf = reshGetVal(confResH, &cdiPioConfOps);
  return conf->maxPathLen;
}

void
cdiPioConfSetXmapNew(int confResH, xmap_new_func_ptr xmap_new)
{
  struct cdiPioConf *conf = reshGetVal(confResH, &cdiPioConfOps);
  conf->xmap_new = xmap_new;
}

xmap_new_func_ptr
cdiPioConfGetXmapNew(int confResH)
{
  struct cdiPioConf *conf = reshGetVal(confResH, &cdiPioConfOps);
  return conf->xmap_new;
}

void
cdiPioConfSetStripeConversion(int confResH, int doConversion)
{
  struct cdiPioConf *conf = reshGetVal(confResH, &cdiPioConfOps);
  conf->stripify = doConversion;
}

int
cdiPioConfGetStripeConversion(int confResH)
{
  struct cdiPioConf *conf = reshGetVal(confResH, &cdiPioConfOps);
  return conf->stripify;
}

void
cdiPioConfSetBatchedRMA(int confResH, int doBatchedRMA)
{
  struct cdiPioConf *conf = reshGetVal(confResH, &cdiPioConfOps);
  conf->batchedRMA = doBatchedRMA;
}

int
cdiPioConfGetBatchedRMA(int confResH)
{
  struct cdiPioConf *conf = reshGetVal(confResH, &cdiPioConfOps);
  return conf->batchedRMA;
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
