#ifndef CDI_ACROSS_H
#define CDI_ACROSS_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_ACROSS

#include "cdi_int.h"

#define ACROSS_DEFAULT_PORT "13859"

typedef struct
{
  char *expid;
  int expver;
} across_info_t;

int across_connect(const char *path, char filemode, stream_t *streamptr);
void across_disconnect(int sock);
int across_write_grib_message(stream_t *streamptr, const void *gribbuffer, size_t nbytes);

#endif

#endif /* CDI_ACROSS_H */
