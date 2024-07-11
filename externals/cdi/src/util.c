#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#define _XOPEN_SOURCE 600

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cdi.h"
#include "cdi_int.h"
#include "cdi_uuid.h"
#include "dmemory.h"

static const char uuidFmt[] = "%02hhx%02hhx%02hhx%02hhx-"
                              "%02hhx%02hhx-%02hhx%02hhx-%02hhx%02hhx-"
                              "%02hhx%02hhx%02hhx%02hhx%02hhx%02hhx";

int
cdiUUID2Str(const unsigned char *uuid, char *uuidstr)
{
  if (uuid == NULL || uuidstr == NULL) return 0;

  int iret = snprintf(uuidstr, uuidNumHexChars + 1, uuidFmt, uuid[0], uuid[1], uuid[2], uuid[3], uuid[4], uuid[5], uuid[6], uuid[7],
                      uuid[8], uuid[9], uuid[10], uuid[11], uuid[12], uuid[13], uuid[14], uuid[15]);

  if (iret != uuidNumHexChars)
    {
      uuidstr[0] = 0;
      iret = -1;
    }

  return iret;
}

int
cdiStr2UUID(const char *uuidstr, unsigned char *uuid)
{
  if (uuid == NULL || uuidstr == NULL || strlen(uuidstr) != uuidNumHexChars) return -1;

  int iret = sscanf(uuidstr, uuidFmt, &uuid[0], &uuid[1], &uuid[2], &uuid[3], &uuid[4], &uuid[5], &uuid[6], &uuid[7], &uuid[8],
                    &uuid[9], &uuid[10], &uuid[11], &uuid[12], &uuid[13], &uuid[14], &uuid[15]);
  if (iret != CDI_UUID_SIZE) return -1;

  return iret;
}

// Returns a malloc'ed string that escapes all spaces and backslashes with backslashes.
char *
cdiEscapeSpaces(const char *string)
{
  // How much memory do we need?
  size_t escapeCount = 0, length = 0;
  for (; string[length]; ++length) escapeCount += string[length] == ' ' || string[length] == '\\';

  char *result = (char *) Malloc(length + escapeCount + 1);
  if (!result) return NULL;

  // Do the escaping.
  for (size_t in = 0, out = 0; in < length; ++out, ++in)
    {
      if (string[in] == ' ' || string[in] == '\\') result[out++] = '\\';
      result[out] = string[in];
    }
  result[length + escapeCount] = 0;  // termination!
  return result;
}

// input: a space terminated string that may contain escaped characters
// output: a new zero terminated string with the escape characters removed
//*outStringEnd points to the terminating character upon return.
char *
cdiUnescapeSpaces(const char *string, const char **outStringEnd)
{
  // How much memory do we need?
  size_t escapeCount = 0, length = 0;
  for (const char *current = string; *current && *current != ' '; current++)
    {
      if (*current == '\\')
        {
          current++, escapeCount++;
          if (!current) return NULL;
        }
      length++;
    }

  char *result = (char *) Malloc(length + 1);
  if (!result) return NULL;

  // Do the unescaping.
  for (size_t in = 0, out = 0; out < length;)
    {
      if (string[in] == '\\') in++;
      result[out++] = string[in++];
    }
  result[length] = 0;  // termination!
  if (outStringEnd) *outStringEnd = &string[length + escapeCount];
  return result;
}

#if defined(HAVE_DECL_UUID_GENERATE) && defined(HAVE_UUID_UUID_H)
#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif
#include <uuid/uuid.h>
void
cdiCreateUUID(unsigned char uuid[CDI_UUID_SIZE])
{
  static int uuid_seeded = 0;
  static char uuid_rand_state[31 * sizeof(long)];
  char *caller_rand_state;
  if (uuid_seeded)
    caller_rand_state = setstate(uuid_rand_state);
  else
    {
      struct timeval tv;
      int status = gettimeofday(&tv, NULL);
      if (status != 0)
        {
          perror("uuid random seed generation failed!");
          exit(1);
        }
      unsigned seed = (unsigned) (tv.tv_sec ^ tv.tv_usec);
      caller_rand_state = initstate(seed, uuid_rand_state, sizeof(uuid_rand_state));
      uuid_seeded = 1;
    }
  uuid_generate(uuid);
  setstate(caller_rand_state);
}
#elif defined(HAVE_DECL_UUID_CREATE) && defined(HAVE_UUID_H)
#ifdef HAVE_DECL_UUID_MAKE_V5
#include <uuid.h>
void
cdiCreateUUID(unsigned char *uuid)
{
  static const char error_stage[][16]
      = { "uuid_create", "uuid_create", "uuid_load", "uuid_make", "uuid_export", "uuid_destroy1", "uuid_destroy2" };
  uuid_t *objuuid = NULL, *nsuuid = NULL;
  int stage = 0;
  uuid_rc_t status;
  if ((status = uuid_create(&objuuid)) == UUID_RC_OK)
    {
      ++stage;
      if ((status = uuid_create(&nsuuid)) == UUID_RC_OK)
        {
          ++stage;
          if ((status = uuid_load(nsuuid, "ns:OID")) == UUID_RC_OK)
            {
              ++stage;
              if ((status = uuid_make(objuuid, UUID_MAKE_V5, nsuuid, cdiLibraryVersion())) == UUID_RC_OK)
                {
                  ++stage;
                  size_t datalen = CDI_UUID_SIZE;
                  status = uuid_export(objuuid, UUID_FMT_BIN, &uuid, &datalen);
                }
            }
        }
    }
  if (status != UUID_RC_OK) Error("failed to generate UUID at stage %s\n", error_stage[stage]);
  stage = 5;
  if ((status = uuid_destroy(nsuuid)) != UUID_RC_OK) Error("failed to generate UUID at stage %s\n", error_stage[stage]);
  ++stage;
  if ((status = uuid_destroy(objuuid)) != UUID_RC_OK) Error("failed to generate UUID at stage %s\n", error_stage[stage]);
}
#else
#include <inttypes.h>
typedef uint8_t u_int8_t;
typedef uint16_t u_int16_t;
typedef uint32_t u_int32_t;
#include <uuid.h>
void
cdiCreateUUID(unsigned char *uuid)
{
  uint32_t status;
  uuid_create((uuid_t *) (void *) uuid, &status);
  if (status != uuid_s_ok)
    {
      perror("uuid generation failed!");
      exit(1);
    }
}
#endif
#else
#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif
void
cdiCreateUUID(unsigned char *uuid)
{
  static int uuid_seeded = 0;
#ifndef _SX
  static char uuid_rand_state[31 * sizeof(long)];
  char *caller_rand_state;
  if (uuid_seeded)
    caller_rand_state = setstate(uuid_rand_state);
  else
    {
#ifdef HAVE_SYS_TIME_H
      struct timeval tv;
      int status = gettimeofday(&tv, NULL);
      if (status != 0)
        {
          perror("failed seed generation!");
          exit(1);
        }
      unsigned seed = tv.tv_sec ^ tv.tv_usec;
#else
      unsigned seed = 0;
#endif
      caller_rand_state = initstate(seed, uuid_rand_state, sizeof(uuid_rand_state));
      uuid_seeded = 1;
    }
  for (size_t i = 0; i < CDI_UUID_SIZE; ++i) uuid[i] = (unsigned char) random();
#else
  unsigned short caller_rand_state[3];
  {
    static unsigned short our_rand_state[3];
    if (!uuid_seeded)
      {
#ifdef HAVE_SYS_TIME_H
        struct timeval tv;
        int status = gettimeofday(&tv, NULL);
        if (status != 0)
          {
            perror("failed seed generation!");
            exit(1);
          }
        unsigned seed = tv.tv_sec ^ tv.tv_usec;
#else
        unsigned seed = 0;
#endif
        our_rand_state[0] = 0x330E;
        our_rand_state[1] = (unsigned short) (seed & 0xFFFFU);
        our_rand_state[2] = (unsigned short) ((seed >> 16) & 0xFFFFU);
      }
    unsigned short *p = seed48(our_rand_state);
    uuid_seeded = 1;
    memcpy(caller_rand_state, p, sizeof(caller_rand_state));
  }
  for (size_t i = 0; i < CDI_UUID_SIZE; ++i) uuid[i] = (unsigned char) lrand48();
#endif
  /* encode variant into msb of octet 8 */
  uuid[8] = (unsigned char) ((uuid[8] & 0x3f) | (1 << 7));
  /* encode version 4 ((pseudo-)random uuid) into msb of octet 7 */
  uuid[7] = (unsigned char) ((uuid[7] & 0x0f) | (4 << 4));
#ifndef _SX
  setstate(caller_rand_state);
#else
  seed48(caller_rand_state);
#endif
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
