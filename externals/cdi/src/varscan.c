#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "cdi.h"
#include "cdi_int.h"
#include "cdi_uuid.h"
#include "cdi_key.h"
#include "dmemory.h"
#include "resource_handle.h"
#include "varscan.h"
#include "vlist.h"
#include "zaxis.h"
#include "subtype.h"

static size_t Vctsize = 0;
static double *Vct = NULL;

static int numberOfVerticalLevels = 0;
static int numberOfVerticalGrid = 0;
static unsigned char uuidVGrid[CDI_UUID_SIZE];

typedef struct
{
  int level1;
  int level2;
  int recID;
  int lindex;
} leveltable_t;

typedef struct
{
  int subtypeIndex;  //  corresponding tile in subtype_t structure (subtype->self)
  int nlevels;
  int levelTableSize;
  leveltable_t *levelTable;
} subtypetable_t;

typedef struct
{
  int varID;
  int param;
  int prec;
  int tsteptype;
  VarScanKeys scanKeys;
  int gridID;
  int zaxistype;
  int ltype1;  // GRIB first level type
  int ltype2;  // GRIB second level type
  int hasBounds;
  int level_sf;
  int level_unit;
  int zaxisID;

  int nsubtypes_alloc;
  int nsubtypes;
  subtypetable_t *recordTable;  // ~ two-dimensional record list [nsubtypes_alloc][levelTableSize]

  int instID;
  int modelID;
  int tableID;
  int comptype;   // compression type
  int complevel;  // compression level
  bool lmissval;
  double missval;
  char *name;

  // meta-data for specification of tiles (currently only GRIB-API:
  subtype_t *tiles;

  cdi_keys_t keys;

  int opt_grib_nentries;                // current no. key-value pairs
  int opt_grib_kvpair_size;             // current allocated size
  opt_key_val_pair_t *opt_grib_kvpair;  // (optional) list of keyword/value pairs
} vartable_t;

static vartable_t *vartable;
static int varTableSize = 0;
static int varTableUsed = 0;

static void
paramInitEntry(int varID, int param)
{
  vartable[varID].varID = varID;
  vartable[varID].param = param;
  vartable[varID].prec = 0;
  vartable[varID].tsteptype = TSTEP_INSTANT;
  varScanKeysInit(&vartable[varID].scanKeys);
  vartable[varID].gridID = CDI_UNDEFID;
  vartable[varID].zaxistype = 0;
  vartable[varID].ltype1 = 0;
  vartable[varID].ltype2 = -1;
  vartable[varID].hasBounds = 0;
  vartable[varID].level_sf = 0;
  vartable[varID].level_unit = 0;
  vartable[varID].recordTable = NULL;
  vartable[varID].nsubtypes_alloc = 0;
  vartable[varID].nsubtypes = 0;
  vartable[varID].instID = CDI_UNDEFID;
  vartable[varID].modelID = CDI_UNDEFID;
  vartable[varID].tableID = CDI_UNDEFID;
  cdiInitKeys(&vartable[varID].keys);
  vartable[varID].comptype = CDI_COMPRESS_NONE;
  vartable[varID].complevel = 1;
  vartable[varID].lmissval = false;
  vartable[varID].missval = 0;
  vartable[varID].name = NULL;
  vartable[varID].tiles = NULL;
}

// Test if a variable specified by the given meta-data has already been registered in "vartable".
static int
varGetEntry(int param, int gridID, int zaxistype, int ltype1, int tsteptype, const char *name, const VarScanKeys *scanKeys,
            const var_tile_t *tiles)
{
  for (int varID = 0; varID < varTableSize; ++varID)
    {
      // testing for "param" implicitly checks if we are beyond the current vartable size:
      if (vartable[varID].param == param)
        {
          const int no_of_tiles = tiles ? tiles->numberOfTiles : -1;
          const int vt_no_of_tiles
              = vartable[varID].tiles ? subtypeGetGlobalDataP(vartable[varID].tiles, SUBTYPE_ATT_NUMBER_OF_TILES) : -1;
          if ((vartable[varID].zaxistype == zaxistype) && (vartable[varID].ltype1 == ltype1)
              && (vartable[varID].tsteptype == tsteptype)
              && (scanKeys == NULL || varScanKeysIsEqual(&vartable[varID].scanKeys, scanKeys)) && (vartable[varID].gridID == gridID)
              && (vt_no_of_tiles == no_of_tiles))
            {
              if (name && name[0] && vartable[varID].name && vartable[varID].name[0])
                {
                  if (str_is_equal(name, vartable[varID].name)) return varID;
                }
              else
                {
                  return varID;
                }
            }
        }
    }

  return -1;
}

static void
varFree(void)
{
  if (CDI_Debug) Message("call to varFree");

  for (int varID = 0; varID < varTableUsed; ++varID)
    {
      if (vartable[varID].recordTable)
        {
          for (int isub = 0; isub < vartable[varID].nsubtypes_alloc; isub++) Free(vartable[varID].recordTable[isub].levelTable);
          Free(vartable[varID].recordTable);
        }

      if (vartable[varID].name) Free(vartable[varID].name);
      if (vartable[varID].tiles) subtypeDestroyPtr(vartable[varID].tiles);

      cdi_keys_t *keysp = &(vartable[varID].keys);
      cdiDeleteVarKeys(keysp);

      if (vartable[varID].opt_grib_kvpair)
        {
          for (int i = 0; i < vartable[varID].opt_grib_nentries; i++)
            {
              if (vartable[varID].opt_grib_kvpair[i].keyword) Free(vartable[varID].opt_grib_kvpair[i].keyword);
            }
          Free(vartable[varID].opt_grib_kvpair);
        }
      vartable[varID].opt_grib_nentries = 0;
      vartable[varID].opt_grib_kvpair_size = 0;
      vartable[varID].opt_grib_kvpair = NULL;
    }

  if (vartable) Free(vartable);
  vartable = NULL;
  varTableSize = 0;
  varTableUsed = 0;

  if (Vct) Free(Vct);
  Vct = NULL;
  Vctsize = 0;
}

// Search for a tile subtype with subtypeIndex == tile_index.
static int
tileGetEntry(int varID, int tile_index)
{
  for (int isub = 0; isub < vartable[varID].nsubtypes; isub++)
    if (vartable[varID].recordTable[isub].subtypeIndex == tile_index) return isub;
  return CDI_UNDEFID;
}

/* Resizes vartable:recordTable data structure, if necessary. */
static int
tileNewEntry(int varID)
{
  int tileID = 0;
  if (vartable[varID].nsubtypes_alloc == 0)
    {
      /* create table for the first time. */
      vartable[varID].nsubtypes_alloc = 2;
      vartable[varID].nsubtypes = 0;
      vartable[varID].recordTable = (subtypetable_t *) Malloc((size_t) vartable[varID].nsubtypes_alloc * sizeof(subtypetable_t));
      if (vartable[varID].recordTable == NULL) SysError("Allocation of leveltable failed!");

      for (int isub = 0; isub < vartable[varID].nsubtypes_alloc; isub++)
        {
          vartable[varID].recordTable[isub].levelTable = NULL;
          vartable[varID].recordTable[isub].levelTableSize = 0;
          vartable[varID].recordTable[isub].nlevels = 0;
          vartable[varID].recordTable[isub].subtypeIndex = CDI_UNDEFID;
        }
    }
  else
    {
      /* data structure large enough; find a free entry. */
      while (tileID < vartable[varID].nsubtypes_alloc)
        {
          if (vartable[varID].recordTable[tileID].levelTable == NULL) break;
          tileID++;
        }
    }

  /* If the table overflows, double its size. */
  if (tileID == vartable[varID].nsubtypes_alloc)
    {
      tileID = vartable[varID].nsubtypes_alloc;
      vartable[varID].nsubtypes_alloc *= 2;
      vartable[varID].recordTable = (subtypetable_t *) Realloc(vartable[varID].recordTable,
                                                               (size_t) vartable[varID].nsubtypes_alloc * sizeof(subtypetable_t));
      if (vartable[varID].recordTable == NULL) SysError("Reallocation of leveltable failed");
      for (int isub = tileID; isub < vartable[varID].nsubtypes_alloc; isub++)
        {
          vartable[varID].recordTable[isub].levelTable = NULL;
          vartable[varID].recordTable[isub].levelTableSize = 0;
          vartable[varID].recordTable[isub].nlevels = 0;
          vartable[varID].recordTable[isub].subtypeIndex = CDI_UNDEFID;
        }
    }

  return tileID;
}

static int
levelNewEntry(int varID, int level1, int level2, int tileID)
{
  int levelID = 0;
  int levelTableSize = vartable[varID].recordTable[tileID].levelTableSize;
  leveltable_t *levelTable = vartable[varID].recordTable[tileID].levelTable;

  // Look for a free slot in levelTable. (Create the table the first time through).
  if (!levelTableSize)
    {
      levelTableSize = 2;
      levelTable = (leveltable_t *) Malloc((size_t) levelTableSize * sizeof(leveltable_t));
      for (int i = 0; i < levelTableSize; i++) levelTable[i].recID = CDI_UNDEFID;
    }
  else
    {
      while (levelID < levelTableSize && levelTable[levelID].recID != CDI_UNDEFID) ++levelID;
    }

  // If the table overflows, double its size.
  if (levelID == levelTableSize)
    {
      levelTable = (leveltable_t *) Realloc(levelTable, (size_t) (levelTableSize *= 2) * sizeof(leveltable_t));
      for (int i = levelID; i < levelTableSize; i++) levelTable[i].recID = CDI_UNDEFID;
    }

  levelTable[levelID].level1 = level1;
  levelTable[levelID].level2 = level2;
  levelTable[levelID].lindex = levelID;

  vartable[varID].recordTable[tileID].nlevels = levelID + 1;
  vartable[varID].recordTable[tileID].levelTableSize = levelTableSize;
  vartable[varID].recordTable[tileID].levelTable = levelTable;

  return levelID;
}

#define UNDEF_PARAM -4711

static int
paramNewEntry(int param)
{
  int varID = 0;

  // Look for a free slot in vartable. (Create the table the first time through).
  if (!varTableSize)
    {
      varTableSize = 2;
      vartable = (vartable_t *) Malloc((size_t) varTableSize * sizeof(vartable_t));
      if (vartable == NULL)
        {
          Message("varTableSize = %d", varTableSize);
          SysError("Allocation of vartable failed");
        }

      for (int i = 0; i < varTableSize; i++)
        {
          vartable[i].param = UNDEF_PARAM;
          vartable[i].opt_grib_kvpair = NULL;
          vartable[i].opt_grib_kvpair_size = 0;
          vartable[i].opt_grib_nentries = 0;
        }
    }
  else
    {
      while (varID < varTableSize)
        {
          if (vartable[varID].param == UNDEF_PARAM) break;
          varID++;
        }
    }

  // If the table overflows, double its size.
  if (varID == varTableSize)
    {
      vartable = (vartable_t *) Realloc(vartable, (size_t) (varTableSize *= 2) * sizeof(vartable_t));
      for (int i = varID; i < varTableSize; i++)
        {
          vartable[i].param = UNDEF_PARAM;
          vartable[i].opt_grib_kvpair = NULL;
          vartable[i].opt_grib_kvpair_size = 0;
          vartable[i].opt_grib_nentries = 0;
        }
    }

  paramInitEntry(varID, param);

  return varID;
}

// Append tile set to a subtype. Return index of the new tile (i.e. the "entry->self" value).
static int
varInsertTileSubtype(vartable_t *vptr, const var_tile_t *tiles)
{
  if (tiles == NULL) return 0;

  // first, generate a subtype based on the info in "tiles".
  subtype_t *subtype_ptr;
  subtypeAllocate(&subtype_ptr, SUBTYPE_TILES);
  subtypeDefGlobalDataP(subtype_ptr, SUBTYPE_ATT_TOTALNO_OF_TILEATTR_PAIRS, tiles->totalno_of_tileattr_pairs);
  subtypeDefGlobalDataP(subtype_ptr, SUBTYPE_ATT_TILE_CLASSIFICATION, tiles->tileClassification);
  subtypeDefGlobalDataP(subtype_ptr, SUBTYPE_ATT_NUMBER_OF_TILES, tiles->numberOfTiles);

  // Here, we create a tile set for comparison that contains only one tile/attribute pair (based on "tiles").
  struct subtype_entry_t *entry = subtypeEntryInsert(subtype_ptr);
  subtypeDefEntryDataP(entry, SUBTYPE_ATT_NUMBER_OF_ATTR, tiles->numberOfAttributes);
  subtypeDefEntryDataP(entry, SUBTYPE_ATT_TILEINDEX, tiles->tileindex);
  subtypeDefEntryDataP(entry, SUBTYPE_ATT_TILEATTRIBUTE, tiles->attribute);

  if (vptr->tiles == NULL)
    {
      vptr->tiles = subtype_ptr;
      return 0;
    }
  else
    {
      tilesetInsertP(vptr->tiles, subtype_ptr);
      subtypeDestroyPtr(subtype_ptr);
      return vptr->tiles->nentries - 1;
    }
}

void
varAddRecord(int recID, int param, int gridID, int zaxistype, int hasBounds, int level1, int level2, int level_sf, int level_unit,
             int prec, int *pvarID, int *plevelID, int tsteptype, int ltype1, int ltype2, const char *name,
             const VarScanKeys *scanKeys, const var_tile_t *tiles, int *tile_index)
{
  int varID = (CDI_Split_Ltype105 != 1 || zaxistype != ZAXIS_HEIGHT)
                  ? varGetEntry(param, gridID, zaxistype, ltype1, tsteptype, name, scanKeys, tiles)
                  : CDI_UNDEFID;

  if (varID == CDI_UNDEFID)
    {
      varTableUsed++;
      varID = paramNewEntry(param);
      vartable[varID].gridID = gridID;
      vartable[varID].zaxistype = zaxistype;
      vartable[varID].ltype1 = ltype1;
      vartable[varID].ltype2 = ltype2;
      vartable[varID].hasBounds = hasBounds;
      vartable[varID].level_sf = level_sf;
      vartable[varID].level_unit = level_unit;
      vartable[varID].tsteptype = tsteptype;
      if (scanKeys) vartable[varID].scanKeys = *scanKeys;

      if (name && name[0]) vartable[varID].name = strdup(name);
    }
  else
    {
      char paramstr[32];
      cdiParamToString(param, paramstr, sizeof(paramstr));

      if (vartable[varID].gridID != gridID)
        {
          Message("param = %s gridID = %d", paramstr, gridID);
          Error("horizontal grid must not change for same parameter!");
        }
      if (vartable[varID].zaxistype != zaxistype)
        {
          Message("param = %s zaxistype = %d", paramstr, zaxistype);
          Error("zaxistype must not change for same parameter!");
        }
    }

  if (prec > vartable[varID].prec) vartable[varID].prec = prec;

  // append current tile to tile subtype info.
  const int this_tile = varInsertTileSubtype(&vartable[varID], tiles);
  int tileID = tileGetEntry(varID, this_tile);
  if (tile_index) (*tile_index) = this_tile;
  if (tileID == CDI_UNDEFID)
    {
      tileID = tileNewEntry((int) varID);
      vartable[varID].recordTable[tileID].subtypeIndex = this_tile;
      vartable[varID].nsubtypes++;
    }

  // append current level to level table info
  const int levelID = levelNewEntry(varID, level1, level2, tileID);
  if (CDI_Debug)
    Message("vartable[%d].recordTable[%d].levelTable[%d].recID = %d; level1,2=%d,%d", varID, tileID, levelID, recID, level1,
            level2);
  vartable[varID].recordTable[tileID].levelTable[levelID].recID = recID;

  *pvarID = (int) varID;
  *plevelID = levelID;
}

/*
static
int dblcmp(const void *s1, const void *s2)
{
  int cmp = 0;

  if      ( *((double *) s1) < *((double *) s2) ) cmp = -1;
  else if ( *((double *) s1) > *((double *) s2) ) cmp =  1;

  return cmp;
}
*/
static int
cmpLevelTable(const void *s1, const void *s2)
{
  int cmp = 0;
  const leveltable_t *x = (const leveltable_t *) s1;
  const leveltable_t *y = (const leveltable_t *) s2;
  // printf("%g %g  %d %d\n", x->leve11, y->level1, x, y);
  if (x->level1 < y->level1)
    cmp = -1;
  else if (x->level1 > y->level1)
    cmp = 1;

  return cmp;
}

static int
cmpLevelTableInv(const void *s1, const void *s2)
{
  int cmp = 0;
  const leveltable_t *x = (const leveltable_t *) s1;
  const leveltable_t *y = (const leveltable_t *) s2;
  // printf("%g %g  %d %d\n", x->leve11, y->level1, x, y);
  if (x->level1 < y->level1)
    cmp = 1;
  else if (x->level1 > y->level1)
    cmp = -1;

  return cmp;
}

void
varCopyKeys(int vlistID, int varID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);
  cdiInitKeys(&vlistptr->vars[varID].keys);
  cdiCopyVarKeys(&vartable[varID].keys, &vlistptr->vars[varID].keys);
}

struct cdi_generate_varinfo
{
  int varid;
  const char *name;
};

/*
static int
cdi_generate_cmp_varname(const void *s1, const void *s2)
{
  const struct cdi_generate_varinfo *x = (const struct cdi_generate_varinfo *) s1, *y = (const struct cdi_generate_varinfo *) s2;
  return strcmp(x->name, y->name);
}
*/

void
cdi_generate_vars(stream_t *streamptr)
{
  const int vlistID = streamptr->vlistID;

  int *varids = (int *) Malloc(varTableUsed * sizeof(int));
  for (int varID = 0; varID < varTableUsed; varID++) varids[varID] = (int) varID;
  /*
  if ( streamptr->sortname )
    {
      size_t varID;
      for (varID = 0; varID < varTableUsed; varID++)
        if (!vartable[varID].name) break;

      if ( varID == varTableUsed )
        {
          struct cdi_generate_varinfo *varInfo
            = (struct cdi_generate_varinfo *) Malloc((size_t)varTableUsed * sizeof(struct cdi_generate_varinfo));

          for (size_t varID = 0; varID < varTableUsed; varID++)
            {
              varInfo[varID].varid = varids[varID];
              varInfo[varID].name = vartable[varids[varID]].name;
            }
          qsort(varInfo, varTableUsed, sizeof(varInfo[0]), cdi_generate_cmp_varname);
          for (size_t varID = 0; varID < varTableUsed; varID++)
            {
              varids[varID] = varInfo[varID].varid;
            }
          Free(varInfo);
        }
    }
  */
  for (int index = 0; index < varTableUsed; index++)
    {
      const int varid = varids[index];

      const int gridID = vartable[varid].gridID;
      const int param = vartable[varid].param;
      const int ltype1 = vartable[varid].ltype1;
      const int ltype2 = vartable[varid].ltype2;
      int zaxistype = vartable[varid].zaxistype;
      if (ltype1 == 0 && zaxistype == ZAXIS_GENERIC && cdiDefaultLeveltype != -1) zaxistype = cdiDefaultLeveltype;
      const int hasBounds = vartable[varid].hasBounds;
      const int prec = vartable[varid].prec;
      int instID = vartable[varid].instID;
      int modelID = vartable[varid].modelID;
      int tableID = vartable[varid].tableID;
      const int tsteptype = vartable[varid].tsteptype;
      const int comptype = vartable[varid].comptype;

      const double level_sf = (vartable[varid].level_sf != 0) ? (1.0 / vartable[varid].level_sf) : 1;

      /* consistency check: test if all subtypes have the same levels: */
      const int nlevels = vartable[varid].recordTable[0].nlevels;
      for (int isub = 1; isub < vartable[varid].nsubtypes; isub++)
        {
          if (vartable[varid].recordTable[isub].nlevels != nlevels)
            {
              fprintf(stderr,
                      "var \"%s\": isub = %d / %d :: "
                      "nlevels = %d, vartable[varid].recordTable[isub].nlevels = %d\n",
                      vartable[varid].name, isub, vartable[varid].nsubtypes, nlevels, vartable[varid].recordTable[isub].nlevels);
              Error("zaxis size must not change for same parameter!");
            }

          const leveltable_t *t1 = vartable[varid].recordTable[isub - 1].levelTable;
          const leveltable_t *t2 = vartable[varid].recordTable[isub].levelTable;
          for (int ilev = 0; ilev < nlevels; ilev++)
            if ((t1[ilev].level1 != t2[ilev].level1) || (t1[ilev].level2 != t2[ilev].level2)
                || (t1[ilev].lindex != t2[ilev].lindex))
              {
                fprintf(stderr,
                        "var \"%s\", varID=%d: isub = %d / %d :: "
                        "nlevels = %d, vartable[varid].recordTable[isub].nlevels = %d\n",
                        vartable[varid].name, varid, isub, vartable[varid].nsubtypes, nlevels,
                        vartable[varid].recordTable[isub].nlevels);
                Message("t1[ilev].level1=%d / t2[ilev].level1=%d", t1[ilev].level1, t2[ilev].level1);
                Message("t1[ilev].level2=%d / t2[ilev].level2=%d", t1[ilev].level2, t2[ilev].level2);
                Message("t1[ilev].lindex=%d / t2[ilev].lindex=%d", t1[ilev].lindex, t2[ilev].lindex);
                Error("zaxis type must not change for same parameter!");
              }
        }
      leveltable_t *levelTable = vartable[varid].recordTable[0].levelTable;

      if (ltype1 == 0 && zaxistype == ZAXIS_GENERIC && nlevels == 1 && levelTable[0].level1 == 0) zaxistype = ZAXIS_SURFACE;

      double *dlevels = (double *) Malloc(nlevels * sizeof(double));

      /*
      if ( hasBounds && zaxistype != ZAXIS_HYBRID && zaxistype != ZAXIS_HYBRID_HALF )
        for (int levelID = 0; levelID < nlevels; levelID++)
          dlevels[levelID] = (level_sf*levelTable[levelID].level1 + level_sf*levelTable[levelID].level2) / 2.0;
      else
      */
      for (int levelID = 0; levelID < nlevels; levelID++) dlevels[levelID] = level_sf * levelTable[levelID].level1;

      if (nlevels > 1)
        {
          bool linc = true, ldec = true, lsort = false;
          for (int levelID = 1; levelID < nlevels; levelID++)
            {
              // check increasing of levels
              linc &= (dlevels[levelID] > dlevels[levelID - 1]);
              // check decreasing of levels
              ldec &= (dlevels[levelID] < dlevels[levelID - 1]);
            }
          /*
           * always sort pressure z-axis to ensure
           * levelTable[levelID1].level1 < levelTable[levelID2].level1 <=> levelID1 > levelID2
           * unless already sorted in decreasing order
           */
          if ((!linc && !ldec) && zaxistype == ZAXIS_PRESSURE)
            {
              qsort(levelTable, nlevels, sizeof(leveltable_t), cmpLevelTableInv);
              lsort = true;
            }
          /*
           * always sort hybrid and depth-below-land z-axis to ensure
           * levelTable[levelID1].level1 < levelTable[levelID2].level1 <=> levelID1 < levelID2
           * unless already sorted in increasing order
           */
          else if ((!linc && !ldec) || zaxistype == ZAXIS_HYBRID || zaxistype == ZAXIS_DEPTH_BELOW_LAND)
            {
              qsort(levelTable, nlevels, sizeof(leveltable_t), cmpLevelTable);
              lsort = true;
            }

          if (lsort)
            {
              /*
              if ( hasBounds && zaxistype != ZAXIS_HYBRID && zaxistype != ZAXIS_HYBRID_HALF )
                for (int levelID = 0; levelID < nlevels; levelID++)
                  dlevels[levelID] = (level_sf*levelTable[levelID].level1 + level_sf*levelTable[levelID].level2) / 2.0;
              else
              */
              for (int levelID = 0; levelID < nlevels; levelID++) dlevels[levelID] = level_sf * levelTable[levelID].level1;
            }
        }

      double *dlevels1 = NULL;
      double *dlevels2 = NULL;
      if (hasBounds)
        {
          dlevels1 = (double *) Malloc(nlevels * sizeof(double));
          for (int levelID = 0; levelID < nlevels; levelID++) dlevels1[levelID] = level_sf * levelTable[levelID].level1;
          dlevels2 = (double *) Malloc(nlevels * sizeof(double));
          for (int levelID = 0; levelID < nlevels; levelID++) dlevels2[levelID] = level_sf * levelTable[levelID].level2;
        }

      const char **cvals = NULL;
      const char *unitptr = cdiUnitNamePtr(vartable[varid].level_unit);
      int zaxisID = varDefZaxis(vlistID, zaxistype, (int) nlevels, dlevels, cvals, 0, hasBounds, dlevels1, dlevels2, (int) Vctsize,
                                Vct, NULL, NULL, unitptr, 0, 0, ltype1, ltype2);

      if (CDI_CMOR_Mode && nlevels == 1 && zaxistype != ZAXIS_HYBRID) zaxisDefScalar(zaxisID);

      if (zaxisInqType(zaxisID) == ZAXIS_REFERENCE)
        {
          if (numberOfVerticalLevels > 0) cdiDefKeyInt(zaxisID, CDI_GLOBAL, CDI_KEY_NLEV, numberOfVerticalLevels);
          if (numberOfVerticalGrid > 0) cdiDefKeyInt(zaxisID, CDI_GLOBAL, CDI_KEY_NUMBEROFVGRIDUSED, numberOfVerticalGrid);
          if (!cdiUUIDIsNull(uuidVGrid)) cdiDefKeyBytes(zaxisID, CDI_GLOBAL, CDI_KEY_UUID, uuidVGrid, CDI_UUID_SIZE);
        }

      if (hasBounds) Free(dlevels1);
      if (hasBounds) Free(dlevels2);
      Free(dlevels);

      // define new subtype for tile set
      int tilesetID = CDI_UNDEFID;
      if (vartable[varid].tiles) tilesetID = vlistDefTileSubtype(vlistID, vartable[varid].tiles);

      // generate new variable
      int varID = stream_new_var(streamptr, gridID, zaxisID, tilesetID);
      varID = vlistDefVarTiles(vlistID, gridID, zaxisID, TIME_VARYING, tilesetID);

      vlistDefVarTsteptype(vlistID, varID, tsteptype);
      vlistDefVarParam(vlistID, varID, param);
      vlistDefVarDatatype(vlistID, varID, prec);
      vlistDefVarCompType(vlistID, varID, comptype);

      varCopyKeys(vlistID, varID);

      if (vartable[varid].lmissval) vlistDefVarMissval(vlistID, varID, vartable[varid].missval);
      if (vartable[varid].name) cdiDefKeyString(vlistID, varID, CDI_KEY_NAME, vartable[varid].name);

      vlist_t *vlistptr = vlist_to_pointer(vlistID);
      for (int i = 0; i < vartable[varid].opt_grib_nentries; i++)
        {
          resize_opt_grib_entries(&vlistptr->vars[varID], vlistptr->vars[varID].opt_grib_nentries + 1);
          vlistptr->vars[varID].opt_grib_nentries += 1;
          const int idx = vlistptr->vars[varID].opt_grib_nentries - 1;

          vlistptr->vars[varID].opt_grib_kvpair[idx] = vartable[varid].opt_grib_kvpair[i];
          vlistptr->vars[varID].opt_grib_kvpair[idx].keyword = NULL;
          if (vartable[varid].opt_grib_kvpair[i].keyword)
            vlistptr->vars[varID].opt_grib_kvpair[idx].keyword = strdupx(vartable[varid].opt_grib_kvpair[i].keyword);
          vlistptr->vars[varID].opt_grib_kvpair[i].update = true;
        }
      // note: if the key is not defined, we do not throw an error!

      if (CDI_Default_TableID != CDI_UNDEFID)
        {
          int pdis, pcat, pnum;
          cdiDecodeParam(param, &pnum, &pcat, &pdis);
          char name[CDI_MAX_NAME];
          name[0] = 0;
          char longname[CDI_MAX_NAME];
          longname[0] = 0;
          char units[CDI_MAX_NAME];
          units[0] = 0;
          tableInqEntry(CDI_Default_TableID, pnum, -1, name, longname, units);
          if (name[0])
            {
              if (tableID != CDI_UNDEFID)
                {
                  cdiDefKeyString(vlistID, varID, CDI_KEY_NAME, name);
                  if (longname[0]) cdiDefKeyString(vlistID, varID, CDI_KEY_LONGNAME, longname);
                  if (units[0]) cdiDefKeyString(vlistID, varID, CDI_KEY_UNITS, units);
                }
              else
                tableID = CDI_Default_TableID;
            }
          if (CDI_Default_ModelID != CDI_UNDEFID) modelID = CDI_Default_ModelID;
          if (CDI_Default_InstID != CDI_UNDEFID) instID = CDI_Default_InstID;
        }

      if (instID != CDI_UNDEFID) vlistDefVarInstitut(vlistID, varID, instID);
      if (modelID != CDI_UNDEFID) vlistDefVarModel(vlistID, varID, modelID);
      if (tableID != CDI_UNDEFID) vlistDefVarTable(vlistID, varID, tableID);
    }

  for (int index = 0; index < varTableUsed; index++)
    {
      int varid = varids[index];
      int nlevels = vartable[varid].recordTable[0].nlevels;

      int nsub = (vartable[varid].nsubtypes >= 0) ? vartable[varid].nsubtypes : 0;
      for (int isub = 0; isub < nsub; isub++)
        {
          sleveltable_t *streamRecordTable = streamptr->vars[index].recordTable + isub;
          leveltable_t *vartableLevelTable = vartable[varid].recordTable[isub].levelTable;
          for (int levelID = 0; levelID < nlevels; levelID++)
            {
              streamRecordTable->recordID[levelID] = vartableLevelTable[levelID].recID;
              int lindex;
              for (lindex = 0; lindex < nlevels; lindex++)
                if (levelID == vartableLevelTable[lindex].lindex) break;
              if (lindex == nlevels) Error("Internal problem! lindex not found.");
              streamRecordTable->lindex[levelID] = (int) lindex;
            }
        }
    }

  Free(varids);

  varFree();
}

void
varDefVCT(size_t vctsize, double *vctptr)
{
  if (Vct == NULL && vctptr != NULL && vctsize > 0)
    {
      Vctsize = vctsize;
      Vct = (double *) Malloc(vctsize * sizeof(double));
      memcpy(Vct, vctptr, vctsize * sizeof(double));
    }
}

void
varDefZAxisReference(int nhlev, int nvgrid, unsigned char uuid[CDI_UUID_SIZE])
{
  numberOfVerticalLevels = nhlev;
  numberOfVerticalGrid = nvgrid;
  memcpy(uuidVGrid, uuid, CDI_UUID_SIZE);
}

bool
zaxis_compare(int zaxisID, int zaxistype, int nlevels, const double *levels, const double *lbounds, const double *ubounds,
              const char *longname, const char *units, int ltype1, int ltype2)
{
  bool differ = true;

  int ltype1_0 = 0, ltype2_0 = -1;
  cdiInqKeyInt(zaxisID, CDI_GLOBAL, CDI_KEY_TYPEOFFIRSTFIXEDSURFACE, &ltype1_0);
  cdiInqKeyInt(zaxisID, CDI_GLOBAL, CDI_KEY_TYPEOFSECONDFIXEDSURFACE, &ltype2_0);
  const bool ltype1IsEqual = (ltype1 == ltype1_0);
  const bool ltype2IsEqual = (ltype2 == ltype2_0);
  const bool hasBounds = (lbounds && ubounds);

  if (ltype1IsEqual && ltype2IsEqual && (zaxistype == zaxisInqType(zaxisID) || zaxistype == ZAXIS_GENERIC))
    {
      const bool hasBoundsZ = (zaxisInqLbounds(zaxisID, NULL) > 0 && zaxisInqUbounds(zaxisID, NULL) > 0);
      if (nlevels == zaxisInqSize(zaxisID) && hasBoundsZ == hasBounds)
        {
          const double *dlevels = zaxisInqLevelsPtr(zaxisID);
          if (dlevels && levels)
            {
              int levelID;
              for (levelID = 0; levelID < nlevels; levelID++)
                {
                  if (fabs(dlevels[levelID] - levels[levelID]) > 1.e-9) break;
                }
              if (levelID == nlevels) differ = false;
            }

          if (!differ && hasBounds)
            {
              double *bounds = (double *) malloc(2 * nlevels * sizeof(double));
              zaxisInqLbounds(zaxisID, bounds);
              zaxisInqUbounds(zaxisID, bounds + nlevels);
              for (int levelID = 0; levelID < nlevels; levelID++)
                {
                  if (fabs(lbounds[levelID] - bounds[levelID]) > 1.e-9
                      || fabs(ubounds[levelID] - bounds[levelID + nlevels]) > 1.e-9)
                    {
                      differ = true;
                      break;
                    }
                }
              free(bounds);
            }

          if (!differ)
            {
              if (longname && longname[0])
                {
                  char zlongname[CDI_MAX_NAME];
                  int length = CDI_MAX_NAME;
                  cdiInqKeyString(zaxisID, CDI_GLOBAL, CDI_KEY_LONGNAME, zlongname, &length);
                  if (zlongname[0] && !str_is_equal(longname, zlongname)) differ = true;
                }
              if (units && units[0])
                {
                  char zunits[CDI_MAX_NAME];
                  int length = CDI_MAX_NAME;
                  cdiInqKeyString(zaxisID, CDI_GLOBAL, CDI_KEY_UNITS, zunits, &length);
                  if (zunits[0] && !str_is_equal(units, zunits)) differ = true;
                }
            }
        }
    }

  return differ;
}

struct varDefZAxisSearchState
{
  int resIDValue;
  int zaxistype;
  int nlevels;
  const double *levels;
  const double *lbounds;
  const double *ubounds;
  const char *longname;
  const char *units;
  int ltype1;
  int ltype2;
};

static enum cdiApplyRet
varDefZAxisSearch(int id, void *res, void *data)
{
  struct varDefZAxisSearchState *state = (struct varDefZAxisSearchState *) data;
  (void) res;
  if (zaxis_compare(id, state->zaxistype, state->nlevels, state->levels, state->lbounds, state->ubounds, state->longname,
                    state->units, state->ltype1, state->ltype2)
      == false)
    {
      state->resIDValue = id;
      return CDI_APPLY_STOP;
    }
  else
    return CDI_APPLY_GO_ON;
}

int
varDefZaxis(int vlistID, int zaxistype, int nlevels, const double *levels, const char **cvals, size_t clength, bool hasBounds,
            const double *levels1, const double *levels2, int vctsize, const double *vct, char *name, const char *longname,
            const char *units, int prec, int mode, int ltype1, int ltype2)
{
  /*
    mode: 0 search in vlist and zaxis table
          1 search in zaxis table
   */
  int zaxisID = CDI_UNDEFID;
  bool zaxisdefined = false;
  bool zaxisglobdefined = false;
  vlist_t *vlistptr = vlist_to_pointer(vlistID);
  int nzaxis = vlistptr->nzaxis;

  if (ltype2 == 255) ltype2 = -1;

  if (mode == 0)
    for (int index = 0; index < nzaxis; index++)
      {
        zaxisID = vlistptr->zaxisIDs[index];

        if (!zaxis_compare(zaxisID, zaxistype, nlevels, levels, levels1, levels2, longname, units, ltype1, ltype2))
          {
            zaxisdefined = true;
            break;
          }
      }

  if (!zaxisdefined)
    {
      struct varDefZAxisSearchState query;
      query.zaxistype = zaxistype;
      query.nlevels = nlevels;
      query.levels = levels;
      query.lbounds = levels1;
      query.ubounds = levels2;
      query.longname = longname;
      query.units = units;
      query.ltype1 = ltype1;
      query.ltype2 = ltype2;

      if ((zaxisglobdefined = (cdiResHFilterApply(getZaxisOps(), varDefZAxisSearch, &query) == CDI_APPLY_STOP)))
        zaxisID = query.resIDValue;

      if (mode == 1 && zaxisglobdefined)
        for (int index = 0; index < nzaxis; index++)
          if (vlistptr->zaxisIDs[index] == zaxisID)
            {
              zaxisglobdefined = false;
              break;
            }
    }

  if (!zaxisdefined)
    {
      if (!zaxisglobdefined)
        {
          zaxisID = zaxisCreate(zaxistype, nlevels);
          if (levels) zaxisDefLevels(zaxisID, levels);
          if (hasBounds)
            {
              zaxisDefLbounds(zaxisID, levels1);
              zaxisDefUbounds(zaxisID, levels2);
            }

          if (cvals != NULL && nlevels != 0 && clength != 0) zaxisDefCvals(zaxisID, cvals, (int) clength);

          if ((zaxistype == ZAXIS_HYBRID || zaxistype == ZAXIS_HYBRID_HALF) && vctsize > 0) zaxisDefVct(zaxisID, vctsize, vct);

          if (name && name[0]) cdiDefKeyString(zaxisID, CDI_GLOBAL, CDI_KEY_NAME, name);
          if (longname && longname[0]) cdiDefKeyString(zaxisID, CDI_GLOBAL, CDI_KEY_LONGNAME, longname);
          if (units && units[0]) cdiDefKeyString(zaxisID, CDI_GLOBAL, CDI_KEY_UNITS, units);
          zaxisDefDatatype(zaxisID, prec);
          cdiDefKeyInt(zaxisID, CDI_GLOBAL, CDI_KEY_TYPEOFFIRSTFIXEDSURFACE, ltype1);
          if (ltype2 != -1) cdiDefKeyInt(zaxisID, CDI_GLOBAL, CDI_KEY_TYPEOFSECONDFIXEDSURFACE, ltype2);
        }

      vlistptr->zaxisIDs[nzaxis] = zaxisID;
      vlistptr->nzaxis++;
    }

  return zaxisID;
}

void
varDefMissval(int varID, double missval)
{
  vartable[varID].lmissval = true;
  vartable[varID].missval = missval;
}

void
varDefCompType(int varID, int comptype)
{
  if (vartable[varID].comptype == CDI_COMPRESS_NONE) vartable[varID].comptype = comptype;
}

void
varDefCompLevel(int varID, int complevel)
{
  vartable[varID].complevel = complevel;
}

int
varInqInst(int varID)
{
  return vartable[varID].instID;
}

void
varDefInst(int varID, int instID)
{
  vartable[varID].instID = instID;
}

int
varInqModel(int varID)
{
  return vartable[varID].modelID;
}

void
varDefModel(int varID, int modelID)
{
  vartable[varID].modelID = modelID;
}

int
varInqTable(int varID)
{
  return vartable[varID].tableID;
}

void
varDefTable(int varID, int tableID)
{
  vartable[varID].tableID = tableID;
}

void
varDefKeyInt(int varID, int key, int value)
{
  cdi_keys_t *keysp = &(vartable[varID].keys);
  cdiDefVarKeyInt(keysp, key, value);
}

void
varDefKeyBytes(int varID, int key, const unsigned char *bytes, int length)
{
  cdi_keys_t *keysp = &(vartable[varID].keys);
  cdiDefVarKeyBytes(keysp, key, bytes, length);
}

void
varDefKeyString(int varID, int key, const char *string)
{
  int length = strlen(string) + 1;
  cdi_keys_t *keysp = &(vartable[varID].keys);
  cdiDefVarKeyBytes(keysp, key, (const unsigned char *) string, length);
}

#ifdef HAVE_LIBGRIB_API
// Resizes and initializes opt_grib_kvpair data structure.
static void
resize_vartable_opt_grib_entries(vartable_t *var, int nentries)
{
  if (var->opt_grib_kvpair_size < nentries)
    {
      if (CDI_Debug) Message("resize data structure, %d -> %d", var->opt_grib_kvpair_size, nentries);

      int new_size = ((2 * var->opt_grib_kvpair_size) > nentries) ? (2 * var->opt_grib_kvpair_size) : nentries;
      if (CDI_Debug) Message("resize vartable opt_grib_entries array to size %d", new_size);
      opt_key_val_pair_t *tmp = (opt_key_val_pair_t *) Malloc((size_t) new_size * sizeof(opt_key_val_pair_t));
      for (int i = 0; i < var->opt_grib_kvpair_size; i++)
        {
          tmp[i] = var->opt_grib_kvpair[i];
        }
      for (int i = var->opt_grib_kvpair_size; i < new_size; i++)
        {
          tmp[i].int_val = 0;
          tmp[i].dbl_val = 0;
          tmp[i].update = false;
          tmp[i].keyword = NULL;
        }  // for
      var->opt_grib_kvpair_size = new_size;
      Free(var->opt_grib_kvpair);
      var->opt_grib_kvpair = tmp;
    }
}
#endif

#ifdef HAVE_LIBGRIB_API
void
varDefOptGribInt(int varID, int tile_index, long lval, const char *keyword)
{
  int idx = -1;
  for (int i = 0; i < vartable[varID].opt_grib_nentries; i++)
    {
      if (str_is_equal(keyword, vartable[varID].opt_grib_kvpair[i].keyword)
          && (vartable[varID].opt_grib_kvpair[i].data_type == t_int)
          && (vartable[varID].opt_grib_kvpair[i].subtype_index == tile_index))
        idx = i;
    }

  if (idx == -1)
    {
      resize_vartable_opt_grib_entries(&vartable[varID], vartable[varID].opt_grib_nentries + 1);
      vartable[varID].opt_grib_nentries += 1;
      idx = vartable[varID].opt_grib_nentries - 1;
    }
  else
    {
      if (vartable[varID].opt_grib_kvpair[idx].keyword) Free(vartable[varID].opt_grib_kvpair[idx].keyword);
    }
  vartable[varID].opt_grib_kvpair[idx].data_type = t_int;
  vartable[varID].opt_grib_kvpair[idx].int_val = (int) lval;
  vartable[varID].opt_grib_kvpair[idx].keyword = strdupx(keyword);
  vartable[varID].opt_grib_kvpair[idx].subtype_index = tile_index;
}
#endif

#ifdef HAVE_LIBGRIB_API
void
varDefOptGribDbl(int varID, int tile_index, double dval, const char *keyword)
{
  int idx = -1;
  for (int i = 0; i < vartable[varID].opt_grib_nentries; i++)
    {
      if (str_is_equal(keyword, vartable[varID].opt_grib_kvpair[i].keyword)
          && (vartable[varID].opt_grib_kvpair[i].data_type == t_double)
          && (vartable[varID].opt_grib_kvpair[i].subtype_index == tile_index))
        idx = i;
    }

  if (idx == -1)
    {
      resize_vartable_opt_grib_entries(&vartable[varID], vartable[varID].opt_grib_nentries + 1);
      vartable[varID].opt_grib_nentries += 1;
      idx = vartable[varID].opt_grib_nentries - 1;
    }
  else
    {
      if (vartable[varID].opt_grib_kvpair[idx].keyword) Free(vartable[varID].opt_grib_kvpair[idx].keyword);
    }
  vartable[varID].opt_grib_kvpair[idx].data_type = t_double;
  vartable[varID].opt_grib_kvpair[idx].dbl_val = dval;
  vartable[varID].opt_grib_kvpair[idx].keyword = strdupx(keyword);
  vartable[varID].opt_grib_kvpair[idx].subtype_index = tile_index;
}
#endif

#ifdef HAVE_LIBGRIB_API
int
varOptGribNentries(int varID)
{
  int nentries = vartable[varID].opt_grib_nentries;
  return nentries;
}
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
