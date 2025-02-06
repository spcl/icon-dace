#include "cdi.h"
#include "cdi_int.h"
#include "zaxis.h"
#include "grid.h"
#include "vlist.h"
#include "resource_unpack.h"

static cdi_keys_t *
vlist_get_keysp(vlist_t *vlistptr, int varID)
{
  if (varID == CDI_GLOBAL) return &vlistptr->keys;
  if (varID >= 0 && varID < vlistptr->nvars) return &(vlistptr->vars[varID].keys);

  return NULL;
}

static cdi_keys_t *
grid_get_keysp(grid_t *gridptr, int varID)
{
  if (varID == CDI_GLOBAL) return &gridptr->keys;
  if (varID == CDI_XAXIS) return &gridptr->x.keys;
  if (varID == CDI_YAXIS) return &gridptr->y.keys;

  return NULL;
}

static cdi_keys_t *
zaxis_get_keysp(zaxis_t *zaxisptr, int varID)
{
  return (varID == CDI_GLOBAL) ? &zaxisptr->keys : NULL;
}

static cdi_key_t *
new_key(cdi_keys_t *keysp, int key)
{
  xassert(keysp != NULL);

  if (keysp->nelems == keysp->nalloc) return NULL;

  cdi_key_t *keyp = &(keysp->value[keysp->nelems]);
  keysp->nelems++;

  keyp->key = key;
  keyp->type = 0;
  keyp->length = 0;
  keyp->v.s = NULL;

  return keyp;
}

cdi_key_t *
find_key(cdi_keys_t *keysp, int key)
{
  xassert(keysp != NULL);

  if (keysp->nelems == 0) return NULL;

  for (size_t keyid = 0; keyid < keysp->nelems; keyid++)
    {
      cdi_key_t *keyp = &(keysp->value[keyid]);
      if (keyp->key == key) return keyp;  // Normal return
    }

  return NULL;
}

static const cdi_key_t *
find_key_const(const cdi_keys_t *keysp, int key)
{
  xassert(keysp != NULL);

  if (keysp->nelems == 0) return NULL;

  for (size_t keyid = 0; keyid < keysp->nelems; keyid++)
    {
      const cdi_key_t *keyp = &(keysp->value[keyid]);
      if (keyp->key == key) return keyp;  // Normal return
    }

  return NULL;
}

static cdi_keys_t *
cdi_get_keysp(int objID, int varID)
{
  if (reshGetTxCode(objID) == GRID) return grid_get_keysp(grid_to_pointer(objID), varID);
  if (reshGetTxCode(objID) == DIST_GRID) return grid_get_keysp(grid_to_pointer(objID), varID);
  if (reshGetTxCode(objID) == ZAXIS) return zaxis_get_keysp(zaxis_to_pointer(objID), varID);
  if (reshGetTxCode(objID) == VLIST) return vlist_get_keysp(vlist_to_pointer(objID), varID);

  return NULL;
}

int
cdi_key_compare(cdi_keys_t *keyspa, cdi_keys_t *keyspb, int keynum)
{
  xassert(keynum >= 0 && keynum < (int) keyspa->nelems && keynum < (int) keyspb->nelems);
  cdi_key_t *keypa = keyspa->value + keynum, *keypb = keyspb->value + keynum;

  if (keypa->key != keypb->key) return 1;
  if (keypa->type != keypb->type) return 1;
  if (keypa->length != keypb->length) return 1;

  if (keypa->type == KEY_BYTES) return (memcmp(keypa->v.s, keypb->v.s, keypa->length) != 0);
  if (keypa->type == KEY_FLOAT) return (IS_NOT_EQUAL(keypa->v.d, keypb->v.d));
  if (keypa->type == KEY_INT) return (keypa->v.i != keypb->v.i);

  return 0;
}

static void
cdi_delete_key(cdi_key_t *keyp)
{
  if (keyp != NULL && keyp->length)  // key in use
    {
      keyp->length = 0;
      if (keyp->type == KEY_BYTES)
        {
          if (keyp->v.s) free(keyp->v.s);
          keyp->v.s = NULL;
        }
      else if (keyp->type == KEY_FLOAT)
        {
          keyp->v.d = 0.0;
        }
      else if (keyp->type == KEY_INT)
        {
          keyp->v.i = 0;
        }
    }
}

void
cdiDeleteVarKeys(cdi_keys_t *keysp)
{
  int nelems = keysp ? (int) keysp->nelems : 0;
  for (int keyid = 0; keyid < nelems; keyid++)
    {
      cdi_delete_key(&(keysp->value[keyid]));
    }

  keysp->nelems = 0;
}

void
cdiDeleteKeys(int cdiID, int varID)
{
  cdi_keys_t *keysp = cdi_get_keysp(cdiID, varID);
  xassert(keysp != NULL);

  cdiDeleteVarKeys(keysp);
}

void
cdiPrintVarKeys(cdi_keys_t *keysp)
{
  int nelems = keysp ? (int) keysp->nelems : 0;
  for (int keyid = 0; keyid < nelems; keyid++)
    {
      cdi_key_t *keyp = &(keysp->value[keyid]);
      if (keyp->length == 0) continue;
      if (keyp->type == KEY_BYTES)
        {
          fprintf(stdout, "%d key %d length %d value %s\n", keyid + 1, keyp->key, keyp->length, keyp->v.s);
        }
      else if (keyp->type == KEY_FLOAT)
        {
          fprintf(stdout, "%d key %d value %g\n", keyid + 1, keyp->key, keyp->v.d);
        }
      else if (keyp->type == KEY_INT)
        {
          fprintf(stdout, "%d key %d value %d\n", keyid + 1, keyp->key, keyp->v.i);
        }
    }
}

void
cdiPrintKeys(int cdiID, int varID)
{
  cdi_keys_t *keysp = cdi_get_keysp(cdiID, varID);
  xassert(keysp != NULL);

  cdiPrintVarKeys(keysp);
}

//  cdiInqKeyLen: Get the length of the string representation of the key
int
cdiInqKeyLen(int cdiID, int varID, int key, int *length)
{
  int status = -1;

  const cdi_keys_t *keysp = cdi_get_keysp(cdiID, varID);
  xassert(keysp != NULL);

  const cdi_key_t *keyp = find_key_const(keysp, key);
  if (keyp != NULL && keyp->length > 0)
    {
      *length = keyp->length;
      status = CDI_NOERR;
    }

  return status;
}

static void
cdi_define_key(const cdi_key_t *keyp, cdi_keys_t *keysp)
{
  // clang-format off
  if      (keyp->type == KEY_INT)   cdiDefVarKeyInt(keysp, keyp->key, keyp->v.i);
  else if (keyp->type == KEY_FLOAT) cdiDefVarKeyFloat(keysp, keyp->key, keyp->v.d);
  else if (keyp->type == KEY_BYTES) cdiDefVarKeyBytes(keysp, keyp->key, keyp->v.s, keyp->length);
  // clang-format on
}

int
cdiDeleteKey(int cdiID, int varID, int key)
{
  int status = CDI_NOERR;

  cdi_keys_t *keysp = cdi_get_keysp(cdiID, varID);
  xassert(keysp != NULL);

  cdi_delete_key(find_key(keysp, key));

  return status;
}

void
cdiCopyVarKeys(const cdi_keys_t *keysp1, cdi_keys_t *keysp2)
{
  for (size_t keyid = 0; keyid < keysp1->nelems; keyid++)
    {
      const cdi_key_t *keyp = &(keysp1->value[keyid]);
      if (keyp->length > 0) cdi_define_key(keyp, keysp2);
    }
}

int
cdiCopyKeys(int cdiID1, int varID1, int cdiID2, int varID2)
{
  int status = CDI_NOERR;

  cdi_keys_t *keysp1 = cdi_get_keysp(cdiID1, varID1);
  xassert(keysp1 != NULL);

  cdi_keys_t *keysp2 = cdi_get_keysp(cdiID2, varID2);
  xassert(keysp2 != NULL);

  cdiCopyVarKeys(keysp1, keysp2);

  return status;
}

int
cdiCopyVarKey(const cdi_keys_t *keysp1, int key, cdi_keys_t *keysp2)
{
  int status = CDI_NOERR;

  const cdi_key_t *keyp = find_key_const(keysp1, key);
  if (keyp == NULL) return -1;

  if (keyp->length > 0) cdi_define_key(keyp, keysp2);

  return status;
}

int
cdiCopyKey(int cdiID1, int varID1, int key, int cdiID2)
{
  cdi_keys_t *keysp1 = cdi_get_keysp(cdiID1, varID1);
  xassert(keysp1 != NULL);

  cdi_keys_t *keysp2 = cdi_get_keysp(cdiID2, varID1);
  xassert(keysp2 != NULL);

  return cdiCopyVarKey(keysp1, key, keysp2);
}

void
cdiDefVarKeyInt(cdi_keys_t *keysp, int key, int value)
{
  cdi_key_t *keyp = find_key(keysp, key);
  if (keyp == NULL) keyp = new_key(keysp, key);

  if (keyp != NULL)
    {
      // if ( keyp->v.i != value )
      {
        keyp->type = KEY_INT;
        keyp->v.i = value;
        keyp->length = 1;
      }
    }
}

/*
@Function  cdiDefKeyInt
@Title     Define an integer value from a key

@Prototype int cdiDefKeyInt(int cdiID, int varID, int key, int value)
@Parameter
    @Item  cdiID    CDI object ID (vlistID, gridID, zaxisID).
    @Item  varID    Variable identifier or CDI_GLOBAL.
    @Item  key      The key to be searched.
    @Item  value    An integer where the data will be read.

@Description
The function @func{cdiDefKeyInt} defines an integer value from a key.

@Result
@func{cdiDefKeyInt} returns CDI_NOERR if OK.

@EndFunction
*/
int
cdiDefKeyInt(int cdiID, int varID, int key, int value)
{
  int status = CDI_NOERR;

  cdi_keys_t *keysp = cdi_get_keysp(cdiID, varID);
  xassert(keysp != NULL);

  cdiDefVarKeyInt(keysp, key, value);

  return status;
}

/*
@Function  cdiInqKeyInt
@Title     Get an integer value from a key

@Prototype int cdiInqKeyInt(int cdiID, int varID, int key, int *value)
@Parameter
    @Item  cdiID    CDI object ID (vlistID, gridID, zaxisID).
    @Item  varID    Variable identifier or CDI_GLOBAL.
    @Item  key      The key to be searched..
    @Item  value    The address of an integer where the data will be retrieved.

@Description
The function @func{cdiInqKeyInt} gets an integer value from a key.

@Result
@func{cdiInqKeyInt} returns CDI_NOERR if key is available.

@EndFunction
*/
int
cdiInqKeyInt(int cdiID, int varID, int key, int *value)
{
  int status = -1;

  // if (varID != CDI_GLOBAL) status = cdiInqKeyInt(cdiID, CDI_GLOBAL, key, value);

  const cdi_keys_t *keysp = cdi_get_keysp(cdiID, varID);
  xassert(keysp != NULL);

  const cdi_key_t *keyp = find_key_const(keysp, key);
  if (keyp != NULL && keyp->length == 1)  // key in use
    {
      if (keyp->type == KEY_INT)
        {
          *value = keyp->v.i;
          status = CDI_NOERR;
        }
    }

  return status;
}

int
cdiInqVarKeyInt(const cdi_keys_t *keysp, int key)
{
  int value = 0;

  const cdi_key_t *keyp = find_key_const(keysp, key);
  if (keyp && keyp->type == KEY_INT) value = keyp->v.i;

  return value;
}

void
cdiDefVarKeyFloat(cdi_keys_t *keysp, int key, double value)
{
  cdi_key_t *keyp = find_key(keysp, key);
  if (keyp == NULL) keyp = new_key(keysp, key);

  if (keyp != NULL)
    {
      keyp->type = KEY_FLOAT;
      keyp->v.d = value;
      keyp->length = 1;
    }
}

/*
@Function  cdiDefKeyFloat
@Title     Define a floating point value from a key

@Prototype int cdiDefKeyFloat(int cdiID, int varID, int key, double value)
@Parameter
    @Item  cdiID    CDI object ID (vlistID, gridID, zaxisID).
    @Item  varID    Variable identifier or CDI_GLOBAL.
    @Item  key      The key to be searched
    @Item  value    A double where the data will be read

@Description
The function @func{cdiDefKeyFloat} defines a CDI floating point value from a key.

@Result
@func{cdiDefKeyFloat} returns CDI_NOERR if OK.

@EndFunction
*/
int
cdiDefKeyFloat(int cdiID, int varID, int key, double value)
{
  int status = CDI_NOERR;

  cdi_keys_t *keysp = cdi_get_keysp(cdiID, varID);
  xassert(keysp != NULL);

  cdiDefVarKeyFloat(keysp, key, value);

  return status;
}

/*
@Function  cdiInqKeyFloat
@Title     Get a floating point value from a key

@Prototype int cdiInqKeyFloat(int cdiID, int varID, int key, double *value)
@Parameter
    @Item  cdiID    CDI object ID (vlistID, gridID, zaxisID).
    @Item  varID    Variable identifier or CDI_GLOBAL.
    @Item  key      The key to be searched.
    @Item  value    The address of a double where the data will be retrieved.

@Description
The function @func{cdiInqKeyFloat} gets a floating point value from a key.

@Result
@func{cdiInqKeyFloat} returns CDI_NOERR if key is available.

@EndFunction
*/
int
cdiInqKeyFloat(int cdiID, int varID, int key, double *value)
{
  int status = -1;

  // if (varID != CDI_GLOBAL) status = cdiInqKeyFloat(cdiID, CDI_GLOBAL, key, value);

  const cdi_keys_t *keysp = cdi_get_keysp(cdiID, varID);
  xassert(keysp != NULL);

  const cdi_key_t *keyp = find_key_const(keysp, key);
  if (keyp != NULL && keyp->length == 1)  // key in use
    {
      if (keyp->type == KEY_FLOAT)
        {
          *value = keyp->v.d;
          status = CDI_NOERR;
        }
    }

  return status;
}

void
cdiDefVarKeyBytes(cdi_keys_t *keysp, int key, const unsigned char *bytes, int length)
{
  cdi_key_t *keyp = find_key(keysp, key);
  if (keyp == NULL) keyp = new_key(keysp, key);

  if (keyp != NULL)
    {
      if (keyp->length != 0 && keyp->length != length)
        {
          if (keyp->v.s) free(keyp->v.s);
          keyp->length = 0;
        }
      if (keyp->length == 0)
        {
          keyp->v.s = (unsigned char *) malloc((size_t) length);
          keyp->length = length;
        }

      memcpy(keyp->v.s, bytes, length);
      keyp->type = KEY_BYTES;
    }
}

/*
@Function  cdiDefKeyBytes
@Title     Define a byte array from a key

@Prototype int cdiDefKeyBytes(int cdiID, int varID, int key, const unsigned char *bytes, int length)
@Parameter
    @Item  cdiID    CDI object ID (vlistID, gridID, zaxisID).
    @Item  varID    Variable identifier or CDI_GLOBAL.
    @Item  key      The key to be searched.
    @Item  bytes    The address of a byte array where the data will be read.
    @Item  length   Length of the byte array

@Description
The function @func{cdiDefKeyBytes} defines a byte array from a key.

@Result
@func{cdiDefKeyBytes} returns CDI_NOERR if OK.

@EndFunction
*/
int
cdiDefKeyBytes(int cdiID, int varID, int key, const unsigned char *bytes, int length)
{
  int status = CDI_NOERR;

  cdi_keys_t *keysp = cdi_get_keysp(cdiID, varID);
  xassert(keysp != NULL);

  cdiDefVarKeyBytes(keysp, key, bytes, length);

  return status;
}

int
cdiInqVarKeyBytes(const cdi_keys_t *keysp, int key, unsigned char *bytes, int *length)
{
  int status = -1;

  const cdi_key_t *keyp = find_key_const(keysp, key);
  if (keyp != NULL && keyp->length > 0)  // key in use
    {
      if (keyp->type == KEY_BYTES)
        {
          if (keyp->length < *length) *length = keyp->length;
          memcpy(bytes, keyp->v.s, *length);
          status = CDI_NOERR;
        }
    }

  return status;
}

//  cdiInqKeyBytes: Get a byte array from a key
/*
@Function  cdiInqKeyBytes
@Title     Get a byte array from a key

@Prototype int cdiInqKeyBytes(int cdiID, int varID, int key, unsigned char *bytes, int *length)
@Parameter
    @Item  cdiID    CDI object ID (vlistID, gridID, zaxisID).
    @Item  varID    Variable identifier or CDI_GLOBAL.
    @Item  key      The key to be searched.
    @Item  bytes    The address of a byte array where the data will be retrieved.
                    The caller must allocate space for the returned byte array.
    @Item  length   The allocated length of the byte array on input.
@Description
The function @func{cdiInqKeyBytes} gets a byte array from a key.

@Result
@func{cdiInqKeyBytes} returns CDI_NOERR if key is available.

@EndFunction
*/
int
cdiInqKeyBytes(int cdiID, int varID, int key, unsigned char *bytes, int *length)
{
  xassert(bytes != NULL);
  xassert(length != NULL);

  // if (varID != CDI_GLOBAL) status = cdiInqKeyBytes(cdiID, CDI_GLOBAL, key, bytes, length);

  const cdi_keys_t *keysp = cdi_get_keysp(cdiID, varID);
  xassert(keysp != NULL);

  return cdiInqVarKeyBytes(keysp, key, bytes, length);
}

/*
@Function  cdiDefKeyString
@Title     Define a string from a key

@Prototype int cdiDefKeyString(int cdiID, int varID, int key, const char *string)
@Parameter
    @Item  cdiID    CDI object ID (vlistID, gridID, zaxisID).
    @Item  varID    Variable identifier or CDI_GLOBAL.
    @Item  key      The key to be searched.
    @Item  string   The address of a string where the data will be read.

@Description
The function @func{cdiDefKeyString} defines a text string from a key.

@Result
@func{cdiDefKeyString} returns CDI_NOERR if OK.

@Example
Here is an example using @func{cdiDefKeyString} to define the name of a variable:

@Source
#include "cdi.h"
   ...
int vlistID, varID, status;
   ...
vlistID = vlistCreate();
varID = vlistDefVar(vlistID, gridID, zaxisID, TIME_VARYING);
   ...
status = cdiDefKeyString(vlistID, varID, CDI_KEY_NAME, "temperature");
   ...
@EndSource
@EndFunction
*/
int
cdiDefKeyString(int cdiID, int varID, int key, const char *string)
{
  xassert(string != NULL);

  int length = strlen(string) + 1;
  int status = cdiDefKeyBytes(cdiID, varID, key, (const unsigned char *) string, length);

  return status;
}

/*
@Function  cdiInqKeyString
@Title     Get a string from a key

@Prototype int cdiInqKeyString(int cdiID, int varID, int key, char *string, int *length)
@Parameter
    @Item  cdiID    CDI object ID (vlistID, gridID, zaxisID).
    @Item  varID    Variable identifier or CDI_GLOBAL.
    @Item  key      The key to be searched.
    @Item  string   The address of a string where the data will be retrieved.
                    The caller must allocate space for the returned string.
    @Item  length   The allocated length of the string on input.
@Description
The function @func{cdiInqKeyString} gets a text string from a key.

@Result
@func{cdiInqKeyString} returns CDI_NOERR if key is available.

@Example
Here is an example using @func{cdiInqKeyString} to get the name of the first variable:

@Source
#include "cdi.h"
   ...
#define STRLEN 256
   ...
int streamID, vlistID, varID, status;
int length = STRLEN;
char name[STRLEN];
   ...
streamID = streamOpenRead(...);
vlistID = streamInqVlist(streamID);
   ...
varID = 0;
status = cdiInqKeyString(vlistID, varID, CDI_KEY_NAME, name, &length);
   ...
@EndSource
@EndFunction
*/
int
cdiInqKeyString(int cdiID, int varID, int key, char *string, int *length)
{
  xassert(string != NULL);
  xassert(length != NULL);

  int maxlength = *length;
  if (maxlength > 0) string[0] = '\0';

  int status = cdiInqKeyBytes(cdiID, varID, key, (unsigned char *) string, length);
  if (CDI_NOERR == status)
    string[maxlength - 1] = '\0';
  else
    *length = 0;

  return status;
}

const char *
cdiInqVarKeyStringPtr(const cdi_keys_t *keysp, int key)
{
  const cdi_key_t *keyp = find_key_const(keysp, key);
  if (keyp != NULL)  // key in use
    {
      if (keyp->type == KEY_BYTES) return (const char *) keyp->v.s;
    }

  return NULL;
}

void
cdiInitKeys(cdi_keys_t *keysp)
{
  keysp->nalloc = MAX_KEYS;
  keysp->nelems = 0;
  for (int i = 0; i < MAX_KEYS; ++i) keysp->value[i].length = 0;
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
