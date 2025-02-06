#include "dmemory.h"
#include "cdi.h"
#include "cdi_int.h"
#include "vlist.h"
#include "vlist_var.h"
#include "namespace.h"
#include "serialize.h"

enum
{
  VLISTVAR_PACK_INT_IDX_FLAG,
  VLISTVAR_PACK_INT_IDX_GRIDID,
  VLISTVAR_PACK_INT_IDX_ZAXISID,
  VLISTVAR_PACK_INT_IDX_TIMETYPE,
  VLISTVAR_PACK_INT_IDX_DATATYPE,
  VLISTVAR_PACK_INT_IDX_PARAM,
  VLISTVAR_PACK_INT_IDX_INSTID,
  VLISTVAR_PACK_INT_IDX_MODELID,
  VLISTVAR_PACK_INT_IDX_TABLEID,
  VLISTVAR_PACK_INT_IDX_MISSVALUSED,
  VLISTVAR_PACK_INT_IDX_COMPTYPE,
  VLISTVAR_PACK_INT_IDX_COMPLEVEL,
  VLISTVAR_PACK_INT_IDX_NLEVS,
  vlistvarNint
};

enum
{
  VLIST_VAR_PACK_DBL_MISSVAL,
  vlistvar_ndbls,
};

int
vlistVarGetPackSize(vlist_t *p, int varID, void *context)
{
  var_t *var = p->vars + varID;
  int varsize
      = serializeGetSize(vlistvarNint, CDI_DATATYPE_INT, context) + serializeGetSize(vlistvar_ndbls, CDI_DATATYPE_FLT64, context);

  if (var->levinfo) varsize += serializeGetSize(4 * zaxisInqSize(var->zaxisID), CDI_DATATYPE_INT, context);
  varsize += serializeKeysGetPackSize(&var->keys, context);
  varsize += cdiAttsGetSize(p, varID, context);

  return varsize;
}

void
vlistVarPack(vlist_t *p, int varID, char *buf, int size, int *position, void *context)
{
  var_t *var = p->vars + varID;
  int nlevs;
  {
    int tempbuf[vlistvarNint];
    tempbuf[VLISTVAR_PACK_INT_IDX_FLAG] = var->flag;
    tempbuf[VLISTVAR_PACK_INT_IDX_GRIDID] = var->gridID;
    tempbuf[VLISTVAR_PACK_INT_IDX_ZAXISID] = var->zaxisID;
    tempbuf[VLISTVAR_PACK_INT_IDX_TIMETYPE] = var->timetype;
    tempbuf[VLISTVAR_PACK_INT_IDX_DATATYPE] = var->datatype;
    tempbuf[VLISTVAR_PACK_INT_IDX_PARAM] = var->param;
    tempbuf[VLISTVAR_PACK_INT_IDX_INSTID] = var->instID;
    tempbuf[VLISTVAR_PACK_INT_IDX_MODELID] = var->modelID;
    tempbuf[VLISTVAR_PACK_INT_IDX_TABLEID] = var->tableID;
    tempbuf[VLISTVAR_PACK_INT_IDX_MISSVALUSED] = (int) var->missvalused;
    tempbuf[VLISTVAR_PACK_INT_IDX_COMPTYPE] = var->comptype;
    tempbuf[VLISTVAR_PACK_INT_IDX_COMPLEVEL] = var->complevel;
    nlevs = var->levinfo ? zaxisInqSize(var->zaxisID) : 0;
    tempbuf[VLISTVAR_PACK_INT_IDX_NLEVS] = nlevs;
    serializePack(tempbuf, vlistvarNint, CDI_DATATYPE_INT, buf, size, position, context);
  }
  {
    double dtempbuf[vlistvar_ndbls];
    dtempbuf[VLIST_VAR_PACK_DBL_MISSVAL] = var->missval;
    serializePack(dtempbuf, vlistvar_ndbls, CDI_DATATYPE_FLT64, buf, size, position, context);
  }
  if (nlevs)
    {
      int *levbuf = (int *) malloc(nlevs * sizeof(int));
      for (int levID = 0; levID < nlevs; ++levID) levbuf[levID] = var->levinfo[levID].flag;
      serializePack(levbuf, nlevs, CDI_DATATYPE_INT, buf, size, position, context);
      for (int levID = 0; levID < nlevs; ++levID) levbuf[levID] = var->levinfo[levID].index;
      serializePack(levbuf, nlevs, CDI_DATATYPE_INT, buf, size, position, context);
      for (int levID = 0; levID < nlevs; ++levID) levbuf[levID] = var->levinfo[levID].mlevelID;
      serializePack(levbuf, nlevs, CDI_DATATYPE_INT, buf, size, position, context);
      for (int levID = 0; levID < nlevs; ++levID) levbuf[levID] = var->levinfo[levID].flevelID;
      free(levbuf);
    }

  serializeKeysPack(&var->keys, buf, size, position, context);

  cdiAttsPack(p, varID, buf, size, position, context);
}

void
vlistVarUnpack(int vlistID, char *buf, int size, int *position, int originNamespace, void *context)
{
#define adaptKey(key) (namespaceAdaptKey((key), originNamespace))
  double dtempbuf[vlistvar_ndbls];
  int tempbuf[vlistvarNint];
  vlist_t *vlistptr = vlist_to_pointer(vlistID);
  serializeUnpack(buf, size, position, tempbuf, vlistvarNint, CDI_DATATYPE_INT, context);
  serializeUnpack(buf, size, position, dtempbuf, vlistvar_ndbls, CDI_DATATYPE_FLT64, context);

  /* ------------------------------------------- */
  /* NOTE: Tile sets  currently not supported!!! */
  /* ------------------------------------------- */

  int newvar = vlistDefVar(vlistID, adaptKey(tempbuf[VLISTVAR_PACK_INT_IDX_GRIDID]),
                           adaptKey(tempbuf[VLISTVAR_PACK_INT_IDX_ZAXISID]), tempbuf[VLISTVAR_PACK_INT_IDX_TIMETYPE]);
  vlistDefVarDatatype(vlistID, newvar, tempbuf[VLISTVAR_PACK_INT_IDX_DATATYPE]);
  vlistDefVarInstitut(vlistID, newvar, adaptKey(tempbuf[VLISTVAR_PACK_INT_IDX_INSTID]));
  vlistDefVarModel(vlistID, newvar, adaptKey(tempbuf[VLISTVAR_PACK_INT_IDX_MODELID]));
  vlistDefVarTable(vlistID, newvar, tempbuf[VLISTVAR_PACK_INT_IDX_TABLEID]);
  // FIXME: changing the table might change the param code
  vlistDefVarParam(vlistID, newvar, tempbuf[VLISTVAR_PACK_INT_IDX_PARAM]);
  if (tempbuf[VLISTVAR_PACK_INT_IDX_MISSVALUSED]) vlistDefVarMissval(vlistID, newvar, dtempbuf[VLIST_VAR_PACK_DBL_MISSVAL]);
  vlistDefVarCompType(vlistID, newvar, tempbuf[VLISTVAR_PACK_INT_IDX_COMPTYPE]);
  vlistDefVarCompLevel(vlistID, newvar, tempbuf[VLISTVAR_PACK_INT_IDX_COMPLEVEL]);
  const int nlevs = tempbuf[VLISTVAR_PACK_INT_IDX_NLEVS];
  var_t *var = vlistptr->vars + newvar;
  if (nlevs)
    {
      int i, flagSetLev = 0;
      cdiVlistCreateVarLevInfo(vlistptr, newvar);

      int *levbuf = (int *) malloc(nlevs * sizeof(int));
      serializeUnpack(buf, size, position, levbuf, nlevs, CDI_DATATYPE_INT, context);
      for (i = 0; i < nlevs; ++i) vlistDefFlag(vlistID, newvar, i, levbuf[i]);
      for (i = 0; i < nlevs; ++i)
        if (levbuf[i] == tempbuf[0]) flagSetLev = i;
      vlistDefFlag(vlistID, newvar, flagSetLev, levbuf[flagSetLev]);
      serializeUnpack(buf, size, position, levbuf, nlevs, CDI_DATATYPE_INT, context);
      for (i = 0; i < nlevs; ++i) vlistDefIndex(vlistID, newvar, i, levbuf[i]);
      serializeUnpack(buf, size, position, levbuf, nlevs, CDI_DATATYPE_INT, context);
      for (i = 0; i < nlevs; ++i) var->levinfo[i].mlevelID = levbuf[i];
      serializeUnpack(buf, size, position, levbuf, nlevs, CDI_DATATYPE_INT, context);
      for (i = 0; i < nlevs; ++i) var->levinfo[i].flevelID = levbuf[i];
      free(levbuf);
    }

  serializeKeysUnpack(buf, size, position, &var->keys, context);

  cdiAttsUnpack(vlistID, newvar, buf, size, position, context);
#undef adaptKey
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
