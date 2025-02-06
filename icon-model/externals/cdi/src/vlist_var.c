#include <float.h> /* FLT_MAX */
#include "dmemory.h"
#include "cdi.h"
#include "cdi_int.h"
#include "vlist.h"
#include "vlist_var.h"
#include "resource_handle.h"
#include "tablepar.h"
#include "namespace.h"
#include "serialize.h"
#include "error.h"

static void
vlistvar_init_entry(var_t *vlistvars, int varID)
{
  var_t *varptr = &vlistvars[varID];
  varptr->flag = false;
  varptr->lvalidrange = false;
  varptr->xyz = 5;  // xyzStorVals[5] == 321
  varptr->missvalused = false;
  varptr->mvarID = varID;
  varptr->fvarID = varID;
  varptr->param = 0;
  varptr->gridID = CDI_UNDEFID;
  varptr->zaxisID = CDI_UNDEFID;
  varptr->timetype = CDI_UNDEFID;
  varptr->tsteptype = TSTEP_INSTANT;
  varptr->datatype = CDI_UNDEFID;
  varptr->instID = CDI_UNDEFID;
  varptr->modelID = CDI_UNDEFID;
  varptr->tableID = CDI_UNDEFID;
  varptr->timave = 0;
  varptr->nsb = 0;
  varptr->missval = CDI_Default_Missval;
  varptr->validrange[0] = VALIDMISS;
  varptr->validrange[1] = VALIDMISS;
  varptr->levinfo = NULL;
  varptr->comptype = CDI_COMPRESS_NONE;
  varptr->complevel = 1;
  varptr->keys.nalloc = MAX_KEYS;
  varptr->keys.nelems = 0;
  for (int i = 0; i < MAX_KEYS; ++i) varptr->keys.value[i].length = 0;
  varptr->atts.nalloc = MAX_ATTRIBUTES;
  varptr->atts.nelems = 0;
  varptr->subtypeID = CDI_UNDEFID;
  varptr->opt_grib_nentries = 0;
  varptr->opt_grib_kvpair_size = 0;
  varptr->opt_grib_kvpair = NULL;
}

static int
vlistvarNewEntry(int vlistID)
{
  int varID = 0;
  vlist_t *vlistptr = vlist_to_pointer(vlistID);
  int vlistvarSize = vlistptr->varsAllocated;
  var_t *vlistvars = vlistptr->vars;
  // Look for a free slot in vlistvar. (Create the table the first time through).
  if (!vlistvarSize)
    {
      vlistvarSize = 2;
      vlistvars = (var_t *) Malloc((size_t) vlistvarSize * sizeof(var_t));
      for (int i = 0; i < vlistvarSize; ++i) vlistvars[i].isUsed = false;
    }
  else
    {
      while (varID < vlistvarSize && vlistvars[varID].isUsed) ++varID;
    }
  // If the table overflows, double its size.
  if (varID == vlistvarSize)
    {
      vlistvars = (var_t *) Realloc(vlistvars, (size_t) (vlistvarSize *= 2) * sizeof(var_t));
      for (int i = varID; i < vlistvarSize; ++i) vlistvars[i].isUsed = false;
    }

  vlistptr->varsAllocated = vlistvarSize;
  vlistptr->vars = vlistvars;

  vlistvar_init_entry(vlistvars, varID);
  vlistvars[varID].isUsed = true;

  return varID;
}

static var_t *
vlistptr_get_varptr(const char *caller, vlist_t *vlistptr, int varID)
{
  if (varID < 0 || varID >= vlistptr->nvars || !vlistptr->vars[varID].isUsed) Errorc("varID %d undefined!", varID);
  return &vlistptr->vars[varID];
}

int
vlistDefVarTiles(int vlistID, int gridID, int zaxisID, int timetype, int tilesetID)
{
  if (CDI_Debug) Message("gridID = %d  zaxisID = %d  timetype = %d", gridID, zaxisID, timetype);

  int varID = vlistvarNewEntry(vlistID);

  vlist_t *vlistptr = vlist_to_pointer(vlistID);
  vlistptr->nvars++;
  var_t *varptr = &vlistptr->vars[varID];
  varptr->gridID = gridID;
  varptr->zaxisID = zaxisID;
  varptr->timetype = timetype;
  varptr->subtypeID = tilesetID;

  if (timetype < 0)
    {
      Warning("Unexpected time type %d, set to TIME_VARYING!", timetype);
      varptr->timetype = TIME_VARYING;
    }

  // Compatibility for release 1.8.3
  if (timetype > 1 && timetype < 15)
    {
      varptr->timetype = TIME_VARYING;
      varptr->tsteptype = timetype;
      static bool printInfo = true;
      if (printInfo)
        {
          printInfo = false;
          fprintf(stdout, "CDI info: The vlistDefVar() function was called with an invalid value for the timetype parameter.\n");
          fprintf(stdout, "CDI info:    This may be an obsolete form of using the vlistDefVar() function.\n");
          fprintf(stdout, "CDI info:    The correct form is:\n");
          fprintf(stdout, "CDI info:       varID = vlistDefVar(vlistID, gridID, zaxisID, timetype)\n");
          fprintf(stdout, "CDI info:       vlistDefVarTsteptype(vlistID, varID, tsteptype)\n");
          fprintf(stdout, "CDI info:          with: timetype = TIME_CONSTANT | TIME_VARYING\n");
          fprintf(stdout, "CDI info:                tsteptype = TSTEP_AVG .... TSTEP_SUM\n");
        }
    }

  vlistAdd2GridIDs(vlistptr, gridID);
  vlistAdd2ZaxisIDs(vlistptr, zaxisID);
  vlistAdd2SubtypeIDs(vlistptr, tilesetID);

  varptr->param = cdiEncodeParam(-(varID + 1), 255, 255);
  reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);

  return varID;
}

/*
@Function  vlistDefVar
@Title     Define a Variable

@Prototype int vlistDefVar(int vlistID, int gridID, int zaxisID, int timetype)
@Parameter
    @Item  vlistID   Variable list ID, from a previous call to @fref{vlistCreate}.
    @Item  gridID    Grid ID, from a previous call to @fref{gridCreate}.
    @Item  zaxisID   Z-axis ID, from a previous call to @fref{zaxisCreate}.
    @Item  timetype  One of the set of predefined CDI timestep types.
                     The valid CDI timestep types are @func{TIME_CONSTANT} and @func{TIME_VARYING}.

@Description
The function @func{vlistDefVar} adds a new variable to vlistID.

@Result
@func{vlistDefVar} returns an identifier to the new variable.

@Example
Here is an example using @func{vlistCreate} to create a variable list
and add a variable with @func{vlistDefVar}.

@Source
#include "cdi.h"
   ...
int vlistID, varID;
   ...
vlistID = vlistCreate();
varID = vlistDefVar(vlistID, gridID, zaxisID, TIME_VARYING);
   ...
streamDefVlist(streamID, vlistID);
   ...
vlistDestroy(vlistID);
   ...
@EndSource
@EndFunction
*/
int
vlistDefVar(int vlistID, int gridID, int zaxisID, int timetype)
{
  // call "vlistDefVarTiles" with a trivial tile index:
  return vlistDefVarTiles(vlistID, gridID, zaxisID, timetype, CDI_UNDEFID);
}

void
cdiVlistCreateVarLevInfo(vlist_t *vlistptr, int varID)
{
  var_t *varptr = vlistptr_get_varptr(__func__, vlistptr, varID);
  xassert(varID >= 0 && varID < vlistptr->nvars && varptr->levinfo == NULL);

  int zaxisID = varptr->zaxisID;
  size_t nlevs = (size_t) zaxisInqSize(zaxisID);

  varptr->levinfo = (levinfo_t *) Malloc(nlevs * sizeof(levinfo_t));

  for (size_t levID = 0; levID < nlevs; ++levID) varptr->levinfo[levID] = DEFAULT_LEVINFO((int) levID);
}

/*
@Function  vlistDefVarParam
@Title     Define the parameter number of a Variable

@Prototype void vlistDefVarParam(int vlistID, int varID, int param)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate}.
    @Item  varID    Variable identifier.
    @Item  param    Parameter number.

@Description
The function @func{vlistDefVarParam} defines the parameter number of a variable.

@EndFunction
*/
void
vlistDefVarParam(int vlistID, int varID, int param)
{
  var_t *varptr = vlistptr_get_varptr(__func__, vlist_to_pointer(vlistID), varID);

  if (varptr->param != param)
    {
      varptr->param = param;
      reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
    }
}

/*
@Function  vlistDefVarCode
@Title     Define the code number of a Variable

@Prototype void vlistDefVarCode(int vlistID, int varID, int code)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate}.
    @Item  varID    Variable identifier.
    @Item  code     Code number.

@Description
The function @func{vlistDefVarCode} defines the code number of a variable.

@EndFunction
*/
void
vlistDefVarCode(int vlistID, int varID, int code)
{
  var_t *varptr = vlistptr_get_varptr(__func__, vlist_to_pointer(vlistID), varID);

  int pnum, pcat, pdis;
  cdiDecodeParam(varptr->param, &pnum, &pcat, &pdis);
  int newParam = cdiEncodeParam(code, pcat, pdis);
  if (varptr->param != newParam)
    {
      varptr->param = newParam;
      reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
    }
}

void
vlistInqVar(int vlistID, int varID, int *gridID, int *zaxisID, int *timetype)
{
  var_t *varptr = vlistptr_get_varptr(__func__, vlist_to_pointer(vlistID), varID);

  *gridID = varptr->gridID;
  *zaxisID = varptr->zaxisID;
  *timetype = varptr->timetype;
}

/*
@Function  vlistInqVarGrid
@Title     Get the Grid ID of a Variable

@Prototype int vlistInqVarGrid(int vlistID, int varID)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate} or @fref{streamInqVlist}.
    @Item  varID    Variable identifier.

@Description
The function @func{vlistInqVarGrid} returns the grid ID of a Variable.

@Result
@func{vlistInqVarGrid} returns the grid ID of the Variable.

@EndFunction
*/
int
vlistInqVarGrid(int vlistID, int varID)
{
  var_t *varptr = vlistptr_get_varptr(__func__, vlist_to_pointer(vlistID), varID);
  return varptr->gridID;
}

/*
@Function  vlistInqVarZaxis
@Title     Get the Zaxis ID of a Variable

@Prototype int vlistInqVarZaxis(int vlistID, int varID)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate} or @fref{streamInqVlist}.
    @Item  varID    Variable identifier.

@Description
The function @func{vlistInqVarZaxis} returns the zaxis ID of a variable.

@Result
@func{vlistInqVarZaxis} returns the zaxis ID of the variable.

@EndFunction
*/
int
vlistInqVarZaxis(int vlistID, int varID)
{
  var_t *varptr = vlistptr_get_varptr(__func__, vlist_to_pointer(vlistID), varID);
  return varptr->zaxisID;
}

/*
@Function  vlistInqVarSubtype
@Title     Get the Subtype ID of a Variable

@Description
The function @func{vlistInqVarSubtype} returns the subtype ID of a variable.

@Result
@func{vlistInqVarSubtype} returns the subtype ID of the variable.

@EndFunction
*/
int
vlistInqVarSubtype(int vlistID, int varID)
{
  var_t *varptr = vlistptr_get_varptr(__func__, vlist_to_pointer(vlistID), varID);
  return varptr->subtypeID;
}

/*
@Function  vlistInqVarParam
@Title     Get the parameter number of a Variable

@Prototype int vlistInqVarParam(int vlistID, int varID)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate} or @fref{streamInqVlist}.
    @Item  varID    Variable identifier.

@Description
The function @func{vlistInqVarParam} returns the parameter number of a variable.

@Result
@func{vlistInqVarParam} returns the parameter number of the variable.

@EndFunction
*/
int
vlistInqVarParam(int vlistID, int varID)
{
  var_t *varptr = vlistptr_get_varptr(__func__, vlist_to_pointer(vlistID), varID);
  return varptr->param;
}

/*
@Function  vlistInqVarCode
@Title     Get the Code number of a Variable

@Prototype int vlistInqVarCode(int vlistID, int varID)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate} or @fref{streamInqVlist}.
    @Item  varID    Variable identifier.

@Description
The function @func{vlistInqVarCode} returns the code number of a variable.

@Result
@func{vlistInqVarCode} returns the code number of the variable.

@EndFunction
*/
int
vlistInqVarCode(int vlistID, int varID)
{
  var_t *varptr = vlistptr_get_varptr(__func__, vlist_to_pointer(vlistID), varID);

  int pdis, pcat, pnum;
  cdiDecodeParam(varptr->param, &pnum, &pcat, &pdis);
  int code = pnum;
  if (pdis != 255) code = -varID - 1;  // GRIB2 Parameter

  int tableID = varptr->tableID;
  if (code < 0 && tableID != -1)
    {
      char name[CDI_MAX_NAME];
      int length = CDI_MAX_NAME;
      (void) cdiInqKeyString(vlistID, varID, CDI_KEY_NAME, name, &length);

      if (name[0]) tableInqParCode(tableID, name, &code);
    }

  return code;
}

/*
@Function  vlistInqVarName
@Title     Get the name of a Variable

@Prototype void vlistInqVarName(int vlistID, int varID, char *name)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate} or @fref{streamInqVlist}.
    @Item  varID    Variable identifier.
    @Item  name     Returned variable name. The caller must allocate space for the
                    returned string. The maximum possible length, in characters, of
                    the string is given by the predefined constant @func{CDI_MAX_NAME}.

@Description
The function @func{vlistInqVarName} returns the name of a variable.

@Result
@func{vlistInqVarName} returns the name of the variable to the parameter name if available,
otherwise the result is an empty string.

@EndFunction
*/
void
vlistInqVarName(int vlistID, int varID, char *name)
{
  int length = CDI_MAX_NAME;
  if (CDI_NOERR != cdiInqKeyString(vlistID, varID, CDI_KEY_NAME, name, &length))
    {
      var_t *varptr = vlistptr_get_varptr(__func__, vlist_to_pointer(vlistID), varID);

      int pdis, pcat, pnum;
      cdiDecodeParam(varptr->param, &pnum, &pcat, &pdis);
      if (pdis == 255)
        {
          int code = pnum;
          int tableID = varptr->tableID;
          tableInqEntry(tableID, code, -1, name, NULL, NULL);
          if (!name[0]) snprintf(name, 8, "var%d", code);
        }
      else
        {
          snprintf(name, CDI_MAX_NAME, "param%d.%d.%d", pnum, pcat, pdis);
        }
    }
}

/*
@Function vlistCopyVarName
@Tatle    Get the name of a Variable in a safe way

@Prototype char* vlistCopyVarName(int vlistId, int varId)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate} or @fref{streamInqVlist}.
    @Item  varID    Variable identifier.

@Return A pointer to a malloc'ed string. Must be cleaned up with Free().

@Description
This is the buffer overflow immune version of vlistInqVarName().
The memory for the returned string is allocated to fit the string via Malloc().

@EndFunction
*/
char *
vlistCopyVarName(int vlistID, int varID)
{
  // If a name is set in the variable description, use that.
  char name[CDI_MAX_NAME];
  int length = CDI_MAX_NAME;
  (void) cdiInqKeyString(vlistID, varID, CDI_KEY_NAME, name, &length);
  if (name[0]) return strdup(name);

  var_t *varptr = vlistptr_get_varptr(__func__, vlist_to_pointer(vlistID), varID);

  // Otherwise we check if we should use the table of parameter descriptions.
  int discipline, category, number;
  cdiDecodeParam(varptr->param, &number, &category, &discipline);
  char *result = NULL;
  if (discipline == 255)
    {
      tableInqEntry(varptr->tableID, number, -1, name, NULL, NULL);
      if (name[0])
        result = strdup(name);
      else
        {
          // No luck, fall back to outputting a name of the format "var<num>".
          size_t len = 3 + 3 * sizeof(int) * CHAR_BIT / 8 + 2;
          result = (char *) Malloc(len);
          snprintf(result, len, "var%d", number);
        }
    }
  else
    {
      size_t len = 5 + 2 + 3 * (3 * sizeof(int) * CHAR_BIT + 1) + 1;
      result = (char *) Malloc(len);
      snprintf(result, len, "param%d.%d.%d", number, category, discipline);
    }
  // Finally, we fall back to outputting a name of the format "param<num>.<cat>.<dis>".
  return result;
}

/*
@Function  vlistInqVarLongname
@Title     Get the longname of a Variable

@Prototype void vlistInqVarLongname(int vlistID, int varID, char *longname)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate} or @fref{streamInqVlist}.
    @Item  varID    Variable identifier.
    @Item  longname Long name of the variable. The caller must allocate space for the
                    returned string. The maximum possible length, in characters, of
                    the string is given by the predefined constant @func{CDI_MAX_NAME}.

@Description
The function @func{vlistInqVarLongname} returns the longname of a variable if available,
otherwise the result is an empty string.

@Result
@func{vlistInqVarLongname} returns the longname of the variable to the parameter longname.

@EndFunction
*/
void
vlistInqVarLongname(int vlistID, int varID, char *longname)
{
  int length = CDI_MAX_NAME;
  if (CDI_NOERR != cdiInqKeyString(vlistID, varID, CDI_KEY_LONGNAME, longname, &length))
    {
      var_t *varptr = vlistptr_get_varptr(__func__, vlist_to_pointer(vlistID), varID);

      int pdis, pcat, pnum;
      cdiDecodeParam(varptr->param, &pnum, &pcat, &pdis);
      if (pdis == 255) tableInqEntry(varptr->tableID, pnum, -1, NULL, longname, NULL);
    }
}

/*
@Function  vlistInqVarStdname
@Title     Get the standard name of a Variable

@Prototype void vlistInqVarStdname(int vlistID, int varID, char *stdname)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate} or @fref{streamInqVlist}.
    @Item  varID    Variable identifier.
    @Item  stdname  Standard name of the variable. The caller must allocate space for the
                    returned string. The maximum possible length, in characters, of
                    the string is given by the predefined constant @func{CDI_MAX_NAME}.

@Description
The function @func{vlistInqVarStdname} returns the standard name of a variable if available,
otherwise the result is an empty string.

@Result
@func{vlistInqVarStdname} returns the standard name of the variable to the parameter stdname.

@EndFunction
*/
void
vlistInqVarStdname(int vlistID, int varID, char *stdname)
{
  int length = CDI_MAX_NAME;
  (void) cdiInqKeyString(vlistID, varID, CDI_KEY_STDNAME, stdname, &length);
}

// obsolete function
int
vlistInqVarTypeOfGeneratingProcess(int vlistID, int varID)
{
  static bool printInfo = true;
  if (printInfo) printInfo = cdiObsoleteInfo(__func__, "cdiInqKeyInt");
  int typeOfGeneratingProcess = 0;
  cdiInqKeyInt(vlistID, varID, CDI_KEY_TYPEOFGENERATINGPROCESS, &typeOfGeneratingProcess);
  return typeOfGeneratingProcess;
}

// obsolete function
void
vlistDefVarTypeOfGeneratingProcess(int vlistID, int varID, int typeOfGeneratingProcess)
{
  static bool printInfo = true;
  if (printInfo) printInfo = cdiObsoleteInfo(__func__, "cdiDefKeyInt");
  cdiDefKeyInt(vlistID, varID, CDI_KEY_TYPEOFGENERATINGPROCESS, typeOfGeneratingProcess);
}

// obsolete function
void
vlistDefVarProductDefinitionTemplate(int vlistID, int varID, int productDefinitionTemplate)
{
  static bool printInfo = true;
  if (printInfo) printInfo = cdiObsoleteInfo(__func__, "cdiDefKeyInt");
  cdiDefKeyInt(vlistID, varID, CDI_KEY_PRODUCTDEFINITIONTEMPLATE, productDefinitionTemplate);
}

/*
@Function  vlistInqVarUnits
@Title     Get the units of a Variable

@Prototype void vlistInqVarUnits(int vlistID, int varID, char *units)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate} or @fref{streamInqVlist}.
    @Item  varID    Variable identifier.
    @Item  units    Units of the variable. The caller must allocate space for the
                    returned string. The maximum possible length, in characters, of
                    the string is given by the predefined constant @func{CDI_MAX_NAME}.

@Description
The function @func{vlistInqVarUnits} returns the units of a variable if available,
otherwise the result is an empty string.

@Result
@func{vlistInqVarUnits} returns the units of the variable to the parameter units.

@EndFunction
*/
void
vlistInqVarUnits(int vlistID, int varID, char *units)
{
  int length = CDI_MAX_NAME;
  if (CDI_NOERR != cdiInqKeyString(vlistID, varID, CDI_KEY_UNITS, units, &length))
    {
      var_t *varptr = vlistptr_get_varptr(__func__, vlist_to_pointer(vlistID), varID);

      int pdis, pcat, pnum;
      cdiDecodeParam(varptr->param, &pnum, &pcat, &pdis);
      if (pdis == 255) tableInqEntry(varptr->tableID, pnum, -1, NULL, NULL, units);
    }
}

// used in MPIOM !
int
vlistInqVarID(int vlistID, int code)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  for (int varID = 0; varID < vlistptr->nvars; ++varID)
    {
      int param = vlistptr->vars[varID].param;
      int pdis, pcat, pnum;
      cdiDecodeParam(param, &pnum, &pcat, &pdis);
      if (pnum == code) return varID;
    }

  return CDI_UNDEFID;
}

SizeType
vlistInqVarSize(int vlistID, int varID)
{
  int zaxisID, gridID, timetype;
  vlistInqVar(vlistID, varID, &gridID, &zaxisID, &timetype);

  SizeType nlevs = (SizeType) zaxisInqSize(zaxisID);
  SizeType gridsize = gridInqSize(gridID);

  return gridsize * nlevs;
}

/*
@Function  vlistInqVarDatatype
@Title     Get the data type of a Variable

@Prototype int vlistInqVarDatatype(int vlistID, int varID)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate} or @fref{streamInqVlist}.
    @Item  varID    Variable identifier.

@Description
The function @func{vlistInqVarDatatype} returns the data type of a variable.

@Result
@func{vlistInqVarDatatype} returns an identifier to the data type of the variable.
The valid CDI data types are @func{CDI_DATATYPE_PACK8}, @func{CDI_DATATYPE_PACK16}, @func{CDI_DATATYPE_PACK24},
@func{CDI_DATATYPE_FLT32}, @func{CDI_DATATYPE_FLT64}, @func{CDI_DATATYPE_INT8}, @func{CDI_DATATYPE_INT16} and
@func{CDI_DATATYPE_INT32}.

@EndFunction
*/
int
vlistInqVarDatatype(int vlistID, int varID)
{
  var_t *varptr = vlistptr_get_varptr(__func__, vlist_to_pointer(vlistID), varID);
  return varptr->datatype;
}

int
vlistInqVarNumber(int vlistID, int varID)
{
  var_t *varptr = vlistptr_get_varptr(__func__, vlist_to_pointer(vlistID), varID);
  return (varptr->datatype == CDI_DATATYPE_CPX32 || varptr->datatype == CDI_DATATYPE_CPX64) ? CDI_COMP : CDI_REAL;
}

static bool
check_range(double value, double lowerBound, double upperBound)
{
  return (value >= lowerBound && value <= upperBound);
}

/*
@Function  vlistDefVarDatatype
@Title     Define the data type of a Variable

@Prototype void vlistDefVarDatatype(int vlistID, int varID, int datatype)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate}.
    @Item  varID    Variable identifier.
    @Item  datatype The data type identifier.
                    The valid CDI data types are @func{CDI_DATATYPE_PACK8}, @func{CDI_DATATYPE_PACK16},
                    @func{CDI_DATATYPE_PACK24}, @func{CDI_DATATYPE_FLT32}, @func{CDI_DATATYPE_FLT64},
                    @func{CDI_DATATYPE_INT8}, @func{CDI_DATATYPE_INT16} and @func{CDI_DATATYPE_INT32}.

@Description
The function @func{vlistDefVarDatatype} defines the data type of a variable.

@EndFunction
*/
void
vlistDefVarDatatype(int vlistID, int varID, int datatype)
{
  var_t *varptr = vlistptr_get_varptr(__func__, vlist_to_pointer(vlistID), varID);

  if (varptr->datatype != datatype)
    {
      varptr->datatype = datatype;

      if (!varptr->missvalused)
        {
          double missval = varptr->missval;
          bool missvalIsDefault = DBL_IS_EQUAL(missval, CDI_Default_Missval);
          if (missvalIsDefault)
            {
              // clang-format off
              switch (datatype)
                {
                case CDI_DATATYPE_INT8:   varptr->missval = -SCHAR_MAX; break;
                case CDI_DATATYPE_UINT8:  varptr->missval =  UCHAR_MAX; break;
                case CDI_DATATYPE_INT16:  varptr->missval = -SHRT_MAX;  break;
                case CDI_DATATYPE_UINT16: varptr->missval =  USHRT_MAX; break;
                case CDI_DATATYPE_INT32:  varptr->missval = -INT_MAX;   break;
                case CDI_DATATYPE_UINT32: varptr->missval =  UINT_MAX;  break;
                case CDI_DATATYPE_FLT32:  varptr->missval =  (float) CDI_Default_Missval;  break;
                }
              // clang-format on
            }
          else
            {
              // clang-format off
              switch (datatype)
                {
                case CDI_DATATYPE_INT8:   varptr->missval = check_range(missval, -SCHAR_MAX, SCHAR_MAX) ? missval : -SCHAR_MAX; break;
                case CDI_DATATYPE_UINT8:  varptr->missval = check_range(missval,          0, UCHAR_MAX) ? missval :  UCHAR_MAX; break;
                case CDI_DATATYPE_INT16:  varptr->missval = check_range(missval,  -SHRT_MAX,  SHRT_MAX) ? missval : -SHRT_MAX;  break;
                case CDI_DATATYPE_UINT16: varptr->missval = check_range(missval,          0, USHRT_MAX) ? missval :  USHRT_MAX; break;
                case CDI_DATATYPE_INT32:  varptr->missval = check_range(missval,   -INT_MAX,   INT_MAX) ? missval : -INT_MAX;   break;
                case CDI_DATATYPE_UINT32: varptr->missval = check_range(missval,          0,  UINT_MAX) ? missval :  UINT_MAX;  break;
                case CDI_DATATYPE_FLT32:  varptr->missval = check_range(missval,   -FLT_MAX,   FLT_MAX) ? missval :  CDI_Default_Missval;  break;
                }
              // clang-format on
            }
        }
      reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
    }
}

void
vlistDefVarInstitut(int vlistID, int varID, int instID)
{
  var_t *varptr = vlistptr_get_varptr(__func__, vlist_to_pointer(vlistID), varID);
  if (varptr->instID != instID)
    {
      varptr->instID = instID;
      reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
    }
}

int
vlistInqVarInstitut(int vlistID, int varID)
{
  var_t *varptr = vlistptr_get_varptr(__func__, vlist_to_pointer(vlistID), varID);
  return varptr->instID;
}

void
vlistDefVarModel(int vlistID, int varID, int modelID)
{
  var_t *varptr = vlistptr_get_varptr(__func__, vlist_to_pointer(vlistID), varID);
  if (varptr->modelID != modelID)
    {
      varptr->modelID = modelID;
      reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
    }
}

int
vlistInqVarModel(int vlistID, int varID)
{
  var_t *varptr = vlistptr_get_varptr(__func__, vlist_to_pointer(vlistID), varID);
  return varptr->modelID;
}

void
vlistDefVarTable(int vlistID, int varID, int tableID)
{
  var_t *varptr = vlistptr_get_varptr(__func__, vlist_to_pointer(vlistID), varID);

  if (varptr->tableID != tableID)
    {
      varptr->tableID = tableID;
      int tablenum = tableInqNum(tableID);
      int pnum, pcat, pdis;
      cdiDecodeParam(varptr->param, &pnum, &pcat, &pdis);
      varptr->param = cdiEncodeParam(pnum, tablenum, pdis);
      reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
    }
}

int
vlistInqVarTable(int vlistID, int varID)
{
  var_t *varptr = vlistptr_get_varptr(__func__, vlist_to_pointer(vlistID), varID);
  return varptr->tableID;
}

/*
@Function  vlistDefVarName
@Title     Define the name of a Variable

@Prototype void vlistDefVarName(int vlistID, int varID, const char *name)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate}.
    @Item  varID    Variable identifier.
    @Item  name     Name of the variable.

@Description
The function @func{vlistDefVarName} defines the name of a variable.

@EndFunction
*/
void
vlistDefVarName(int vlistID, int varID, const char *name)
{
  (void) cdiDefKeyString(vlistID, varID, CDI_KEY_NAME, name);
}

/*
@Function  vlistDefVarLongname
@Title     Define the long name of a Variable

@Prototype void vlistDefVarLongname(int vlistID, int varID, const char *longname)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate}.
    @Item  varID    Variable identifier.
    @Item  longname Long name of the variable.

@Description
The function @func{vlistDefVarLongname} defines the long name of a variable.

@EndFunction
*/
void
vlistDefVarLongname(int vlistID, int varID, const char *longname)
{
  if (longname) (void) cdiDefKeyString(vlistID, varID, CDI_KEY_LONGNAME, longname);
}

/*
@Function  vlistDefVarStdname
@Title     Define the standard name of a Variable

@Prototype void vlistDefVarStdname(int vlistID, int varID, const char *stdname)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate}.
    @Item  varID    Variable identifier.
    @Item  stdname  Standard name of the variable.

@Description
The function @func{vlistDefVarStdname} defines the standard name of a variable.

@EndFunction
*/
void
vlistDefVarStdname(int vlistID, int varID, const char *stdname)
{
  if (stdname) (void) cdiDefKeyString(vlistID, varID, CDI_KEY_STDNAME, stdname);
}

/*
@Function  vlistDefVarUnits
@Title     Define the units of a Variable

@Prototype void vlistDefVarUnits(int vlistID, int varID, const char *units)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate}.
    @Item  varID    Variable identifier.
    @Item  units    Units of the variable.

@Description
The function @func{vlistDefVarUnits} defines the units of a variable.

@EndFunction
*/
void
vlistDefVarUnits(int vlistID, int varID, const char *units)
{
  if (units) (void) cdiDefKeyString(vlistID, varID, CDI_KEY_UNITS, units);
}

/*
@Function  vlistInqVarMissval
@Title     Get the missing value of a Variable

@Prototype double vlistInqVarMissval(int vlistID, int varID)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate} or @fref{streamInqVlist}.
    @Item  varID    Variable identifier.

@Description
The function @func{vlistInqVarMissval} returns the missing value of a variable.

@Result
@func{vlistInqVarMissval} returns the missing value of the variable.

@EndFunction
*/
double
vlistInqVarMissval(int vlistID, int varID)
{
  var_t *varptr = vlistptr_get_varptr(__func__, vlist_to_pointer(vlistID), varID);
  return varptr->missval;
}

/*
@Function  vlistDefVarMissval
@Title     Define the missing value of a Variable

@Prototype void vlistDefVarMissval(int vlistID, int varID, double missval)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate}.
    @Item  varID    Variable identifier.
    @Item  missval  Missing value.

@Description
The function @func{vlistDefVarMissval} defines the missing value of a variable.

@EndFunction
*/
void
vlistDefVarMissval(int vlistID, int varID, double missval)
{
  var_t *varptr = vlistptr_get_varptr(__func__, vlist_to_pointer(vlistID), varID);
  varptr->missval = missval;
  varptr->missvalused = true;
}

int
vlistInqVarValidrange(int vlistID, int varID, double *validrange)
{
  var_t *varptr = vlistptr_get_varptr(__func__, vlist_to_pointer(vlistID), varID);

  if (validrange != NULL && varptr->lvalidrange)
    {
      validrange[0] = varptr->validrange[0];
      validrange[1] = varptr->validrange[1];
    }

  return (int) varptr->lvalidrange;
}

void
vlistDefVarValidrange(int vlistID, int varID, const double *validrange)
{
  var_t *varptr = vlistptr_get_varptr(__func__, vlist_to_pointer(vlistID), varID);
  varptr->validrange[0] = validrange[0];
  varptr->validrange[1] = validrange[1];
  varptr->lvalidrange = true;
  reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
}

void
vlistDefVarTimetype(int vlistID, int varID, int timetype)
{
  var_t *varptr = vlistptr_get_varptr(__func__, vlist_to_pointer(vlistID), varID);
  if (varptr->timetype != timetype)
    {
      varptr->timetype = timetype;
      reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
    }
}

int
vlistInqVarTimetype(int vlistID, int varID)
{
  var_t *varptr = vlistptr_get_varptr(__func__, vlist_to_pointer(vlistID), varID);
  return varptr->timetype;
}

void
vlistDefVarTsteptype(int vlistID, int varID, int tsteptype)
{
  var_t *varptr = vlistptr_get_varptr(__func__, vlist_to_pointer(vlistID), varID);
  if (varptr->tsteptype != tsteptype)
    {
      varptr->tsteptype = tsteptype;
      reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
    }
}

/*
@Function  vlistInqVarTsteptype
@Title     Get the timestep type of a Variable

@Prototype int vlistInqVarTsteptype(int vlistID, int varID)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate} or @fref{streamInqVlist}.
    @Item  varID    Variable identifier.

@Description
The function @func{vlistInqVarTsteptype} returns the timestep type of a Variable.

@Result
@func{vlistInqVarTsteptype} returns the timestep type of the Variable,
one of the set of predefined CDI timestep types.
The valid CDI timestep types are @func{TSTEP_INSTANT},
@func{TSTEP_ACCUM}, @func{TSTEP_AVG}, @func{TSTEP_MAX}, @func{TSTEP_MIN} and @func{TSTEP_SD}.

@EndFunction
*/
int
vlistInqVarTsteptype(int vlistID, int varID)
{
  var_t *varptr = vlistptr_get_varptr(__func__, vlist_to_pointer(vlistID), varID);
  return varptr->tsteptype;
}

int
vlistInqVarMissvalUsed(int vlistID, int varID)
{
  var_t *varptr = vlistptr_get_varptr(__func__, vlist_to_pointer(vlistID), varID);
  return (int) varptr->missvalused;
}

void
vlistDefFlag(int vlistID, int varID, int levID, int flag)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);
  var_t *varptr = vlistptr_get_varptr(__func__, vlistptr, varID);

  levinfo_t li = DEFAULT_LEVINFO(levID);
  if (varptr->levinfo)
    ;
  else if (flag != li.flag)
    cdiVlistCreateVarLevInfo(vlistptr, varID);
  else
    return;

  varptr->levinfo[levID].flag = flag;
  varptr->flag = 0;

  int nlevs = zaxisInqSize(varptr->zaxisID);
  for (int levelID = 0; levelID < nlevs; ++levelID)
    {
      if (varptr->levinfo[levelID].flag)
        {
          varptr->flag = 1;
          break;
        }
    }

  reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
}

int
vlistInqFlag(int vlistID, int varID, int levID)
{
  var_t *varptr = vlistptr_get_varptr(__func__, vlist_to_pointer(vlistID), varID);

  if (varptr->levinfo)
    return varptr->levinfo[levID].flag;
  else
    {
      levinfo_t li = DEFAULT_LEVINFO(levID);
      return li.flag;
    }
}

int
vlistFindVar(int vlistID, int fvarID)
{
  int varID;
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  for (varID = 0; varID < vlistptr->nvars; ++varID)
    {
      if (vlistptr->vars[varID].fvarID == fvarID) break;
    }

  if (varID == vlistptr->nvars)
    {
      varID = -1;
      Warning("varID not found for fvarID %d in vlistID %d!", fvarID, vlistID);
    }

  return varID;
}

int
vlistFindLevel(int vlistID, int fvarID, int flevelID)
{
  int levelID = -1;
  int varID = vlistFindVar(vlistID, fvarID);
  var_t *varptr = vlistptr_get_varptr(__func__, vlist_to_pointer(vlistID), varID);

  if (varID != -1)
    {
      int nlevs = zaxisInqSize(varptr->zaxisID);
      for (levelID = 0; levelID < nlevs; ++levelID)
        {
          if (varptr->levinfo[levelID].flevelID == flevelID) break;
        }

      if (levelID == nlevs)
        {
          levelID = -1;
          Warning("levelID not found for fvarID %d and levelID %d in vlistID %d!", fvarID, flevelID, vlistID);
        }
    }

  return levelID;
}

int
vlistMergedVar(int vlistID, int varID)
{
  var_t *varptr = vlistptr_get_varptr(__func__, vlist_to_pointer(vlistID), varID);
  return varptr->mvarID;
}

int
vlistMergedLevel(int vlistID, int varID, int levelID)
{
  var_t *varptr = vlistptr_get_varptr(__func__, vlist_to_pointer(vlistID), varID);

  if (varptr->levinfo)
    return varptr->levinfo[levelID].mlevelID;
  else
    {
      levinfo_t li = DEFAULT_LEVINFO(levelID);
      return li.mlevelID;
    }
}

void
vlistDefIndex(int vlistID, int varID, int levelID, int index)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);
  var_t *varptr = vlistptr_get_varptr(__func__, vlistptr, varID);

  levinfo_t li = DEFAULT_LEVINFO(levelID);
  if (varptr->levinfo)
    ;
  else if (index != li.index)
    cdiVlistCreateVarLevInfo(vlistptr, varID);
  else
    return;

  varptr->levinfo[levelID].index = index;
  reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
}

int
vlistInqIndex(int vlistID, int varID, int levelID)
{
  var_t *varptr = vlistptr_get_varptr(__func__, vlist_to_pointer(vlistID), varID);

  if (varptr->levinfo)
    return varptr->levinfo[levelID].index;
  else
    {
      levinfo_t li = DEFAULT_LEVINFO(levelID);
      return li.index;
    }
}

void
vlistChangeVarZaxis(int vlistID, int varID, int zaxisID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);
  var_t *varptr = vlistptr_get_varptr(__func__, vlistptr, varID);

  int nlevs1 = zaxisInqSize(varptr->zaxisID);
  int nlevs2 = zaxisInqSize(zaxisID);

  if (nlevs1 != nlevs2) Error("Number of levels must not change!");

  int nvars = vlistptr->nvars;
  int found = 0;
  int oldZaxisID = varptr->zaxisID;
  for (int i = 0; i < varID; ++i) found |= (vlistptr->vars[i].zaxisID == oldZaxisID);
  for (int i = varID + 1; i < nvars; ++i) found |= (vlistptr->vars[i].zaxisID == oldZaxisID);

  if (found)
    {
      int nzaxis = vlistptr->nzaxis;
      for (int i = 0; i < nzaxis; ++i)
        if (vlistptr->zaxisIDs[i] == oldZaxisID) vlistptr->zaxisIDs[i] = zaxisID;
    }
  else
    vlistAdd2ZaxisIDs(vlistptr, zaxisID);

  varptr->zaxisID = zaxisID;
  reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
}

void
vlistChangeVarGrid(int vlistID, int varID, int gridID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);
  var_t *varptr = vlistptr_get_varptr(__func__, vlistptr, varID);

  int nvars = vlistptr->nvars;
  int index;
  for (index = 0; index < nvars; index++)
    if (index != varID)
      if (vlistptr->vars[index].gridID == vlistptr->vars[varID].gridID) break;

  if (index == nvars)
    {
      for (index = 0; index < vlistptr->ngrids; index++)
        if (vlistptr->gridIDs[index] == vlistptr->vars[varID].gridID) vlistptr->gridIDs[index] = gridID;
    }
  else
    vlistAdd2GridIDs(vlistptr, gridID);

  varptr->gridID = gridID;
  reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
}

void
vlistDefVarCompType(int vlistID, int varID, int comptype)
{
  var_t *varptr = vlistptr_get_varptr(__func__, vlist_to_pointer(vlistID), varID);

  if (varptr->comptype != comptype)
    {
      varptr->comptype = comptype;
      reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
    }
}

int
vlistInqVarCompType(int vlistID, int varID)
{
  var_t *varptr = vlistptr_get_varptr(__func__, vlist_to_pointer(vlistID), varID);
  return varptr->comptype;
}

void
vlistDefVarCompLevel(int vlistID, int varID, int complevel)
{
  var_t *varptr = vlistptr_get_varptr(__func__, vlist_to_pointer(vlistID), varID);

  if (varptr->complevel != complevel)
    {
      varptr->complevel = complevel;
      reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
    }
}

int
vlistInqVarCompLevel(int vlistID, int varID)
{
  var_t *varptr = vlistptr_get_varptr(__func__, vlist_to_pointer(vlistID), varID);
  return varptr->complevel;
}

void
vlistDefVarNSB(int vlistID, int varID, int nsb)
{
  var_t *varptr = vlistptr_get_varptr(__func__, vlist_to_pointer(vlistID), varID);
  if (varptr->nsb != nsb)
    {
      varptr->nsb = nsb;
      reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
    }
}

int
vlistInqVarNSB(int vlistID, int varID)
{
  var_t *varptr = vlistptr_get_varptr(__func__, vlist_to_pointer(vlistID), varID);
  return varptr->nsb;
}

static int
vlistEncodeXyz(const int dimorder[3])
{
  return (short) (dimorder[0] * 100 + dimorder[1] * 10 + dimorder[2]);
}

static void
vlistDecodeXyz(int xyz, int outDimorder[3])
{
  outDimorder[0] = xyz / 100, xyz = xyz % 100;
  outDimorder[1] = xyz / 10;
  outDimorder[2] = xyz % 10;
}

static const short xyzStorVals[] = { 123, 132, 213, 231, 312, 321 };
enum
{
  numXYZStorVals = sizeof(xyzStorVals) / sizeof(xyzStorVals[0])
};

void
vlistDefVarXYZ(int vlistID, int varID, int xyz)
{
  var_t *varptr = vlistptr_get_varptr(__func__, vlist_to_pointer(vlistID), varID);

  if (xyz == 3) xyz = 321;

  // check xyz dimension order
  {
    int dimorder[3];
    vlistDecodeXyz(xyz, dimorder);
    int dimx = 0, dimy = 0, dimz = 0;
    for (int id = 0; id < 3; ++id)
      {
        switch (dimorder[id])
          {
          case 1: dimx++; break;
          case 2: dimy++; break;
          case 3: dimz++; break;
          default: dimorder[id] = 0; break;  // Ensure that we assign a valid dimension to this position.
          }
      }
    if (dimz > 1 || dimy > 1 || dimx > 1)
      xyz = 321;  // ZYX
    else
      {
        // clang-format off
        if (dimz == 0) for (int id = 0; id < 3; ++id) if (dimorder[id] == 0) { dimorder[id] = 3; break; }
        if (dimy == 0) for (int id = 0; id < 3; ++id) if (dimorder[id] == 0) { dimorder[id] = 2; break; }
        if (dimx == 0) for (int id = 0; id < 3; ++id) if (dimorder[id] == 0) { dimorder[id] = 1; break; }
        // clang-format on
        xyz = vlistEncodeXyz(dimorder);
      }
  }

  assert(xyz == 123 || xyz == 312 || xyz == 231 || xyz == 321 || xyz == 132 || xyz == 213);

  for (size_t i = 0; i < numXYZStorVals; ++i)
    if (xyz == xyzStorVals[i])
      {
        xyz = (int) i;
        break;
      }
  varptr->xyz = (signed char) xyz;
  reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
}

void
vlistInqVarDimorder(int vlistID, int varID, int outDimorder[3])
{
  var_t *varptr = vlistptr_get_varptr(__func__, vlist_to_pointer(vlistID), varID);
  vlistDecodeXyz(xyzStorVals[varptr->xyz], outDimorder);
}

int
vlistInqVarXYZ(int vlistID, int varID)
{
  var_t *varptr = vlistptr_get_varptr(__func__, vlist_to_pointer(vlistID), varID);
  return xyzStorVals[varptr->xyz];
}

int
vlistVarCompare(vlist_t *a, int varIDA, vlist_t *b, int varIDB)
{
  xassert(a && b && varIDA >= 0 && varIDA < a->nvars && varIDB >= 0 && varIDB < b->nvars);
  var_t *pva = a->vars + varIDA, *pvb = b->vars + varIDB;
#define FCMP(f) ((pva->f) != (pvb->f))
#define FCMPFLT(f) (IS_NOT_EQUAL((pva->f), (pvb->f)))
#define FCMPSTR(fs) ((pva->fs) != (pvb->fs) && strcmp((pva->fs), (pvb->fs)))
#define FCMP2(f) (namespaceResHDecode(pva->f).idx != namespaceResHDecode(pvb->f).idx)
  int diff = (int) FCMP(fvarID) | FCMP(mvarID) | FCMP(flag) | FCMP(param) | FCMP(datatype) | FCMP(timetype) | FCMP(tsteptype)
             | FCMP(xyz) | FCMP2(gridID) | FCMP2(zaxisID) | FCMP2(instID) | FCMP2(modelID) | FCMP2(tableID) | FCMP(missvalused)
             | FCMPFLT(missval) | FCMP(comptype) | FCMP(complevel) | FCMP(lvalidrange) | FCMPFLT(validrange[0])
             | FCMPFLT(validrange[1]);
#undef FCMP
#undef FCMPFLT
#undef FCMPSTR
#undef FCMP2
  if ((diff |= ((pva->levinfo == NULL) ^ (pvb->levinfo == NULL)))) return 1;

  if (pva->levinfo)
    {
      int zaxisID = pva->zaxisID;
      size_t nlevs = (size_t) zaxisInqSize(zaxisID);
      diff |= (memcmp(pva->levinfo, pvb->levinfo, sizeof(levinfo_t) * nlevs) != 0);
      if (diff) return 1;
    }

  size_t natts = a->vars[varIDA].atts.nelems;
  if (natts != b->vars[varIDB].atts.nelems) return 1;
  for (size_t attID = 0; attID < natts; ++attID) diff |= cdi_att_compare(&a->vars[varIDA].atts, &b->vars[varIDB].atts, (int) attID);

  size_t nkeys = a->vars[varIDA].keys.nelems;
  if (nkeys != b->vars[varIDB].keys.nelems) return 1;
  for (size_t keyID = 0; keyID < nkeys; ++keyID) diff |= cdi_key_compare(&a->vars[varIDA].keys, &b->vars[varIDB].keys, (int) keyID);

  return diff;
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
