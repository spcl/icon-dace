#ifndef CDI_KEY_H
#define CDI_KEY_H

#include "cdi_limits.h"

// CDI key
typedef struct
{
  int key;     // CDI key
  int type;    // KEY_INT, KEY_FLOAT, KEY_BYTES
  int length;  // number of bytes in v.s
  union
  {
    int i;
    double d;
    unsigned char *s;
  } v;
} cdi_key_t;

typedef struct
{
  size_t nalloc;  // number allocated >= nelems
  size_t nelems;  // length of the array
  cdi_key_t value[MAX_KEYS];
} cdi_keys_t;

enum
{
  KEY_INT = 1,
  KEY_FLOAT,
  KEY_BYTES
};

void cdiDefVarKeyInt(cdi_keys_t *keysp, int key, int value);
void cdiDefVarKeyFloat(cdi_keys_t *keysp, int key, double value);
void cdiDefVarKeyBytes(cdi_keys_t *keysp, int key, const unsigned char *bytes, int length);
int cdiInqVarKeyInt(const cdi_keys_t *keysp, int key);
int cdiInqVarKeyBytes(const cdi_keys_t *keysp, int key, unsigned char *bytes, int *length);

cdi_key_t *find_key(cdi_keys_t *keysp, int key);
const char *cdiInqVarKeyStringPtr(const cdi_keys_t *keysp, int key);

static inline const char *
cdiInqVarKeyString(const cdi_keys_t *keysp, int key)
{
  const char *string = cdiInqVarKeyStringPtr(keysp, key);
  if (string == NULL) string = "";
  return string;
}

int cdiCopyVarKey(const cdi_keys_t *keysp1, int key, cdi_keys_t *keysp2);
void cdiCopyVarKeys(const cdi_keys_t *keysp1, cdi_keys_t *keysp2);
void cdiDeleteVarKeys(cdi_keys_t *keysp);
void cdiDeleteKeys(int cdiID, int varID);
void cdiPrintKeys(int cdiID, int varID);

void cdiInitKeys(cdi_keys_t *keysp);

int cdi_key_compare(cdi_keys_t *keyspa, cdi_keys_t *keyspb, int keynum);

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
