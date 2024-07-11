#include <limits.h>

#include "dmemory.h"
#include "cdi.h"
#include "cdi_int.h"
#include "model.h"
#include "namespace.h"
#include "resource_handle.h"
#include "resource_unpack.h"
#include "serialize.h"

#undef CDI_UNDEFID
#define CDI_UNDEFID -1

typedef struct
{
  int self;
  int instID;
  int modelgribID;
  char *name;
} model_t;

static void modelInit(void);

static int modelCompareP(void *modelptr1, void *modelptr2);
static void modelDestroyP(void *modelptr);
static void modelPrintP(void *modelptr, FILE *fp);
static int modelGetSizeP(void *modelptr, void *context);
static void modelPackP(void *modelptr, void *buff, int size, int *position, void *context);
static int modelTxCode(void *modelptr);

static const resOps modelOps = { modelCompareP, modelDestroyP, modelPrintP, modelGetSizeP, modelPackP, modelTxCode };

static void
modelDefaultValue(model_t *modelptr)
{
  modelptr->self = CDI_UNDEFID;
  modelptr->instID = CDI_UNDEFID;
  modelptr->modelgribID = CDI_UNDEFID;
  modelptr->name = NULL;
}

static model_t *
modelNewEntry(cdiResH resH, int instID, int modelgribID, const char *name)
{
  model_t *modelptr = (model_t *) Malloc(sizeof(model_t));
  modelDefaultValue(modelptr);
  if (resH == CDI_UNDEFID)
    modelptr->self = reshPut(modelptr, &modelOps);
  else
    {
      modelptr->self = resH;
      reshReplace(resH, modelptr, &modelOps);
    }
  modelptr->instID = instID;
  modelptr->modelgribID = modelgribID;
  if (name && *name) modelptr->name = strdup(name);

  return (modelptr);
}

void
modelDefaultEntries(void)
{
  int instID;
  enum
  {
    nDefModels = 10
  };
  cdiResH resH[nDefModels];

  instID = institutInq(0, 0, "ECMWF", NULL);
  /* (void)    modelDef(instID, 131, "ERA15"); */
  /* (void)    modelDef(instID, 199, "ERA40"); */

  instID = institutInq(98, 232, "MPIMET", NULL);
  resH[0] = modelDef(instID, 64, "ECHAM5.4");
  resH[1] = modelDef(instID, 63, "ECHAM5.3");
  resH[2] = modelDef(instID, 62, "ECHAM5.2");
  resH[3] = modelDef(instID, 61, "ECHAM5.1");

  instID = institutInq(98, 255, "MPIMET", NULL);
  resH[4] = modelDef(instID, 60, "ECHAM5.0");
  resH[5] = modelDef(instID, 50, "ECHAM4");
  resH[6] = modelDef(instID, 110, "MPIOM1");

  instID = institutInq(0, 0, "DWD", NULL);
  resH[7] = modelDef(instID, 149, "GME");

  instID = institutInq(0, 0, "MCH", NULL);
  //(void)  = modelDef(instID, 137, "COSMO");
  resH[8] = modelDef(instID, 255, "COSMO");

  instID = institutInq(0, 1, "NCEP", NULL);
  resH[9] = modelDef(instID, 80, "T62L28MRF");

  /* pre-defined models are not synchronized */
  for (int i = 0; i < nDefModels; i++) reshSetStatus(resH[i], &modelOps, RESH_IN_USE);
}

static void
modelInit(void)
{
  static bool modelInitialized = false;
  if (modelInitialized) return;
}

struct modelLoc
{
  const char *name;
  int instID, modelgribID, resID;
};

static enum cdiApplyRet
findModelByID(int resID, void *res, void *data)
{
  model_t *modelptr = (model_t *) res;
  struct modelLoc *ret = (struct modelLoc *) data;
  int instID = ret->instID, modelgribID = ret->modelgribID;
  if (modelptr->instID == instID && modelptr->modelgribID == modelgribID)
    {
      ret->resID = resID;
      return CDI_APPLY_STOP;
    }
  else
    return CDI_APPLY_GO_ON;
}

static enum cdiApplyRet
findModelByName(int resID, void *res, void *data)
{
  model_t *modelptr = (model_t *) res;
  struct modelLoc *ret = (struct modelLoc *) data;
  int instID = ret->instID, modelgribID = ret->modelgribID;
  const char *name = ret->name;
  if ((instID == -1 || modelptr->instID == instID) && (modelgribID == 0 || modelptr->modelgribID == modelgribID) && modelptr->name)
    {
      const char *p = name, *q = modelptr->name;
      while (*p != '\0' && *p == *q) ++p, ++q;
      if (*p == '\0' || *q == '\0')
        {
          ret->resID = resID;
          return CDI_APPLY_STOP;
        }
    }
  return CDI_APPLY_GO_ON;
}

int
modelInq(int instID, int modelgribID, const char *name)
{
  modelInit();

  struct modelLoc searchState = { .name = name, .instID = instID, .modelgribID = modelgribID, .resID = CDI_UNDEFID };
  if (name && *name)
    cdiResHFilterApply(&modelOps, findModelByName, &searchState);
  else
    cdiResHFilterApply(&modelOps, findModelByID, &searchState);
  return searchState.resID;
}

int
modelDef(int instID, int modelgribID, const char *name)
{
  model_t *modelptr;

  modelInit();

  modelptr = modelNewEntry(CDI_UNDEFID, instID, modelgribID, name);

  return modelptr->self;
}

int
modelInqInstitut(int modelID)
{
  model_t *modelptr = NULL;

  modelInit();

  if (modelID != CDI_UNDEFID) modelptr = (model_t *) reshGetVal(modelID, &modelOps);

  return modelptr ? modelptr->instID : CDI_UNDEFID;
}

int
modelInqGribID(int modelID)
{
  model_t *modelptr = NULL;

  modelInit();

  if (modelID != CDI_UNDEFID) modelptr = (model_t *) reshGetVal(modelID, &modelOps);

  return modelptr ? modelptr->modelgribID : CDI_UNDEFID;
}

const char *
modelInqNamePtr(int modelID)
{
  model_t *modelptr = NULL;

  modelInit();

  if (modelID != CDI_UNDEFID) modelptr = (model_t *) reshGetVal(modelID, &modelOps);

  return modelptr ? modelptr->name : NULL;
}

static int
modelCompareP(void *modelptr1, void *modelptr2)
{
  model_t *model1 = (model_t *) modelptr1, *model2 = (model_t *) modelptr2;
  int diff = (namespaceResHDecode(model1->instID).idx != namespaceResHDecode(model2->instID).idx)
             | (model1->modelgribID != model2->modelgribID) | !str_is_equal(model1->name, model2->name);
  return diff;
}

void
modelDestroyP(void *modelptr)
{
  model_t *mp = (model_t *) modelptr;
  if (mp->name) Free(mp->name);
  Free(mp);
}

void
modelPrintP(void *modelptr, FILE *fp)
{
  model_t *mp = (model_t *) modelptr;
  fprintf(fp,
          "#\n"
          "# modelID %d\n"
          "#\n"
          "self          = %d\n"
          "instID        = %d\n"
          "modelgribID   = %d\n"
          "name          = %s\n",
          mp->self, mp->self, mp->instID, mp->modelgribID, mp->name ? mp->name : "NN");
}

static int
modelTxCode(void *modelptr)
{
  (void) modelptr;
  return MODEL;
}

enum
{
  MODEL_PACK_INT_SELF,
  MODEL_PACK_INT_INSTID,
  MODEL_PACK_INT_MODELGRIBID,
  MODEL_PACK_INT_NAMELEN,
  modelNints,
};

static int
modelGetSizeP(void *modelptr, void *context)
{
  model_t *p = (model_t *) modelptr;
  size_t txsize = (size_t) serializeGetSize(modelNints, CDI_DATATYPE_INT, context)
                  + (size_t) serializeGetSize(p->name ? (int) strlen(p->name) : 0, CDI_DATATYPE_TXT, context);
  xassert(txsize <= INT_MAX);
  return (int) txsize;
}

static void
modelPackP(void *modelptr, void *buf, int size, int *position, void *context)
{
  model_t *p = (model_t *) modelptr;
  int tempbuf[modelNints];
  tempbuf[MODEL_PACK_INT_SELF] = p->self;
  tempbuf[MODEL_PACK_INT_INSTID] = p->instID;
  tempbuf[MODEL_PACK_INT_MODELGRIBID] = p->modelgribID;
  tempbuf[MODEL_PACK_INT_NAMELEN] = p->name ? (int) strlen(p->name) : 0;
  serializePack(tempbuf, modelNints, CDI_DATATYPE_INT, buf, size, position, context);
  if (p->name) serializePack(p->name, tempbuf[MODEL_PACK_INT_NAMELEN], CDI_DATATYPE_TXT, buf, size, position, context);
}

int
modelUnpack(void *buf, int size, int *position, int originNamespace, void *context, int force_id)
{
#define adaptKey(key) (namespaceAdaptKey((key), originNamespace))
  int tempbuf[modelNints];
  char *name;
  serializeUnpack(buf, size, position, tempbuf, modelNints, CDI_DATATYPE_INT, context);
  if (tempbuf[MODEL_PACK_INT_NAMELEN] != 0)
    {
      size_t len = (size_t) tempbuf[MODEL_PACK_INT_NAMELEN];
      name = (char *) Malloc(len + 1);
      serializeUnpack(buf, size, position, name, tempbuf[MODEL_PACK_INT_NAMELEN], CDI_DATATYPE_TXT, context);
      name[len] = '\0';
    }
  else
    {
      name = (char *) "";
    }
  int targetID = adaptKey(tempbuf[MODEL_PACK_INT_SELF]);
  model_t *mp = modelNewEntry(force_id ? targetID : CDI_UNDEFID, adaptKey(tempbuf[MODEL_PACK_INT_INSTID]),
                              tempbuf[MODEL_PACK_INT_MODELGRIBID], name);
  if (tempbuf[MODEL_PACK_INT_NAMELEN] != 0) Free(name);
  xassert(!force_id || (mp->self == adaptKey(tempbuf[0])));
  reshSetStatus(mp->self, &modelOps, reshGetStatus(mp->self, &modelOps) & ~RESH_SYNC_BIT);
#undef adaptKey
  return mp->self;
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
