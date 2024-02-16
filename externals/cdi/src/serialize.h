#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifndef SERIALIZE_H
#define SERIALIZE_H

#include <string.h>

#include "cdi.h"
#ifndef CDI_CKSUM_H_
#include "cdi_cksum.h"
#endif
#ifndef CDI_KEY_H_
#include "cdi_key.h"
#endif
#ifndef ERROR_H
#include "error.h"
#endif

/*
 * Generic interfaces for (de-)marshalling
 */
int serializeGetSize(int count, int datatype, void *context);
void serializePack(const void *data, int count, int datatype, void *buf, int buf_size, int *position, void *context);
void serializeUnpack(const void *buf, int buf_size, int *position, void *data, int count, int datatype, void *context);

/*
 * (de-)marshalling function for key/value structures
 */
static inline int
serializeKeysGetPackSize(const cdi_keys_t *keysp, void *context)
{
  int packBuffSize = 0;

  int nelems = keysp->nelems;
  packBuffSize += serializeGetSize(1, CDI_DATATYPE_INT, context);
  for (int keyid = 0; keyid < nelems; keyid++)
    {
      const cdi_key_t *keyp = &(keysp->value[keyid]);
      int type = keyp->type;
      packBuffSize += serializeGetSize(1, CDI_DATATYPE_INT, context);  // key
      packBuffSize += serializeGetSize(1, CDI_DATATYPE_INT, context);  // type
      if (type == KEY_BYTES)
        {
          int length = keyp->length;
          packBuffSize += serializeGetSize(1, CDI_DATATYPE_INT, context) + serializeGetSize(length, CDI_DATATYPE_TXT, context);
        }
      else if (type == KEY_INT)
        {
          packBuffSize += serializeGetSize(1, CDI_DATATYPE_INT, context);
        }
      else if (type == KEY_FLOAT)
        {
          packBuffSize += serializeGetSize(1, CDI_DATATYPE_FLT64, context);
        }
    }
  packBuffSize += serializeGetSize(1, CDI_DATATYPE_UINT32, context);
  return packBuffSize;
}

static inline void
serializeKeysPack(const cdi_keys_t *keysp, void *buf, int buf_size, int *position, void *context)
{
  uint32_t d = 0;

  int nelems = keysp->nelems;
  serializePack(&nelems, 1, CDI_DATATYPE_INT, buf, buf_size, position, context);
  for (int keyid = 0; keyid < nelems; keyid++)
    {
      const cdi_key_t *keyp = &(keysp->value[keyid]);
      int key = keyp->key;
      int type = keyp->type;
      serializePack(&key, 1, CDI_DATATYPE_INT, buf, buf_size, position, context);
      serializePack(&type, 1, CDI_DATATYPE_INT, buf, buf_size, position, context);
      if (type == KEY_BYTES)
        {
          int length = keyp->length;
          serializePack(&length, 1, CDI_DATATYPE_INT, buf, buf_size, position, context);
          serializePack(keyp->v.s, length, CDI_DATATYPE_TXT, buf, buf_size, position, context);
          d ^= cdiCheckSum(CDI_DATATYPE_TXT, length, keyp->v.s);
        }
      else if (type == KEY_INT)
        {
          serializePack(&keyp->v.i, 1, CDI_DATATYPE_INT, buf, buf_size, position, context);
        }
      else if (type == KEY_FLOAT)
        {
          serializePack(&keyp->v.d, 1, CDI_DATATYPE_FLT64, buf, buf_size, position, context);
        }
    }

  serializePack(&d, 1, CDI_DATATYPE_UINT32, buf, buf_size, position, context);
}

static inline void
serializeKeysUnpack(const void *buf, int buf_size, int *position, cdi_keys_t *keysp, void *context)
{
  uint32_t d, d2 = 0;
  void *buffer = NULL;
  int buffersize = 0;

  int nelems;
  serializeUnpack(buf, buf_size, position, &nelems, 1, CDI_DATATYPE_INT, context);
  for (int i = 0; i < nelems; ++i)
    {
      int key, type;
      serializeUnpack(buf, buf_size, position, &key, 1, CDI_DATATYPE_INT, context);
      serializeUnpack(buf, buf_size, position, &type, 1, CDI_DATATYPE_INT, context);
      if (type == KEY_BYTES)
        {
          int length;
          serializeUnpack(buf, buf_size, position, &length, 1, CDI_DATATYPE_INT, context);
          if (length > buffersize)
            {
              buffersize = length;
              buffer = realloc(buffer, buffersize);
            }
          serializeUnpack(buf, buf_size, position, buffer, length, CDI_DATATYPE_TXT, context);
          cdiDefVarKeyBytes(keysp, key, (unsigned char *) buffer, length);
          d2 ^= cdiCheckSum(CDI_DATATYPE_TXT, length, buffer);
        }
      else if (type == KEY_INT)
        {
          int ival;
          serializeUnpack(buf, buf_size, position, &ival, 1, CDI_DATATYPE_INT, context);
          cdiDefVarKeyInt(keysp, key, ival);
        }
      else if (type == KEY_FLOAT)
        {
          double dval;
          serializeUnpack(buf, buf_size, position, &dval, 1, CDI_DATATYPE_FLT64, context);
          cdiDefVarKeyFloat(keysp, key, dval);
        }
    }
  serializeUnpack(buf, buf_size, position, &d, 1, CDI_DATATYPE_UINT32, context);
  xassert(d == d2);
  if (buffer) free(buffer);
}

/*
 * (de-)marshalling function for common data structures
 */
static inline int
serializeStrTabGetPackSize(const char **strTab, int numStr, void *context)
{
  xassert(numStr >= 0);
  int packBuffSize = 0;
  for (size_t i = 0; i < (size_t) numStr; ++i)
    {
      size_t len = strlen(strTab[i]);
      packBuffSize += serializeGetSize(1, CDI_DATATYPE_INT, context) + serializeGetSize((int) len, CDI_DATATYPE_TXT, context);
    }
  packBuffSize += serializeGetSize(1, CDI_DATATYPE_UINT32, context);
  return packBuffSize;
}

static inline void
serializeStrTabPack(const char **strTab, int numStr, void *buf, int buf_size, int *position, void *context)
{
  uint32_t d = 0;
  xassert(numStr >= 0);
  for (size_t i = 0; i < (size_t) numStr; ++i)
    {
      int len = (int) strlen(strTab[i]);
      serializePack(&len, 1, CDI_DATATYPE_INT, buf, buf_size, position, context);
      serializePack(strTab[i], len, CDI_DATATYPE_TXT, buf, buf_size, position, context);
      d ^= cdiCheckSum(CDI_DATATYPE_TXT, len, strTab[i]);
    }
  serializePack(&d, 1, CDI_DATATYPE_UINT32, buf, buf_size, position, context);
}

static inline void
serializeStrTabUnpack(const void *buf, int buf_size, int *position, char **strTab, int numStr, void *context)
{
  uint32_t d, d2 = 0;
  xassert(numStr >= 0);
  for (size_t i = 0; i < (size_t) numStr; ++i)
    {
      int len;
      serializeUnpack(buf, buf_size, position, &len, 1, CDI_DATATYPE_INT, context);
      serializeUnpack(buf, buf_size, position, strTab[i], len, CDI_DATATYPE_TXT, context);
      strTab[i][len] = '\0';
      d2 ^= cdiCheckSum(CDI_DATATYPE_TXT, len, strTab[i]);
    }
  serializeUnpack(buf, buf_size, position, &d, 1, CDI_DATATYPE_UINT32, context);
  xassert(d == d2);
}

/*
 * Interfaces for marshalling within a single memory domain
 */
int serializeGetSizeInCore(int count, int datatype, void *context);
void serializePackInCore(const void *data, int count, int datatype, void *buf, int buf_size, int *position, void *context);
void serializeUnpackInCore(const void *buf, int buf_size, int *position, void *data, int count, int datatype, void *context);

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
